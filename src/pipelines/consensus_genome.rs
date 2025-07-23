use std::fs::File;
use std::sync::Arc;
use std::collections::HashMap;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::Technology;
use std::process::Command;
use std::time::Instant;
use tokio::fs::File as TokioFile;
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::mpsc::Receiver;
use serde::Serialize;
use crate::utils::command::{generate_cli, check_versions};
use crate::utils::file::{extension_remover, file_path_manipulator, write_parse_output_to_temp, write_vecu8_to_file};
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, parse_and_filter_fastq_id, validate_sequence, validate_sequence_parallel};
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd, parse_fastq, parse_bytes, y_junction, bytes_to_lines};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand, KRAKEN2_TAG, BCFTOOLS_TAG, BcftoolsSubcommand, MAFFT_TAG, QUAST_TAG, SEQKIT_TAG, SeqkitSubcommand};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::command::bcftools::BcftoolsConfig;
use crate::utils::command::seqkit::SeqkitConfig;
use crate::utils::db::{get_index, retrieve_h5_seq};
use tokio::sync::mpsc;
use futures::future::try_join_all;
use crate::utils::streams::ToBytes;
use crate::config::defs::RunConfig;
use crate::utils::command::quast::QuastConfig;
use crate::utils::stats::{parse_samtools_stats, parse_samtools_depth, compute_depth_stats, parse_seqkit_stats, parse_ercc_stats, compute_allele_counts, compute_coverage_bins};
use crate::utils::vcf::parse_vcf_stream;
use crate::utils::plotting::plot_depths;

const ERCC_FASTA: &str = "ercc_sequences.fasta";

#[derive(Serialize)]
struct Stats {
    sample_name: String,
    depth_avg: f64,
    depth_q25: f64,
    depth_q50: f64,
    depth_q75: f64,
    depth_frac_above_10x: f64,
    depth_frac_above_25x: f64,
    depth_frac_above_50x: f64,
    depth_frac_above_100x: f64,
    allele_counts: HashMap<char, u64>,
    total_reads: u64,
    mapped_reads: u64,
    mapped_paired: Option<u64>,
    paired_inward: Option<u64>,
    paired_outward: Option<u64>,
    paired_other_orientation: Option<u64>,
    ercc_mapped_reads: Option<u64>,
    ercc_mapped_paired: Option<u64>,
    ref_snps: u64,
    ref_mnps: u64,
    ref_indels: u64,
    n_actg: u64,
    n_missing: u64,
    n_gap: u64,
    n_ambiguous: u64,
    coverage_breadth: f64,
    max_aligned_length: usize,
    total_length: usize,
    coverage_bin_size: f64,
    coverage: Vec<(usize, f64, f64, u8, u8)>,
}

pub async fn run(config: &RunConfig) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", config.args.module);

    let cwd = config.cwd.clone();
    let ram_temp_dir = config.ram_temp_dir.clone();
    let mut temp_files = Vec::new();

    // Initialize cleanup tasks
    let mut cleanup_tasks: Vec<tokio::task::JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<tokio::sync::oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();
    let mut quast_write_tasks: Vec<tokio::task::JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut stats_tasks: Vec<tokio::task::JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut ercc_stats_task: Option<tokio::task::JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>> = None;

    // External tools check
    check_versions(vec![SAMTOOLS_TAG, MINIMAP2_TAG, FASTP_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, SEQKIT_TAG, QUAST_TAG]).await?;

    // Arguments and files check
    let file1_path: PathBuf = match &config.args.file1 {
        Some(file) => {
            let file1_full_path = file_path_manipulator(&PathBuf::from(file), &cwd, None, None, "");
            if file1_full_path.exists() {
                file1_full_path
            } else {
                return Err(anyhow!("File1 path does not exist: {}", file1_full_path.display()));
            }
        }
        None => {
            return Err(anyhow!("File1 path required"));
        }
    };

    let sample_base: String;
    let file1_r1r2 = r1r2_base(&file1_path);
    match file1_r1r2.file_name {
        Some(prefix) => { sample_base = prefix; }
        None => {
            eprintln!("No R1 tag found. Using bare file 1 stem as sample_base.");
            sample_base = file1_path.to_string_lossy().into_owned();
        },
    }
    if !file1_path.exists() {
        return Err(anyhow!("Specified file1 {:?} does not exist.", file1_path));
    }

    let sample_base_buf: PathBuf = PathBuf::from(&sample_base);
    let (no_ext_sample_base_buf, _) = extension_remover(&sample_base_buf);
    let no_ext_sample_base = no_ext_sample_base_buf.to_string_lossy().into_owned();
    
    let file2_path: Option<PathBuf> = match &config.args.file2 {
        Some(file) => {
            let file2_full_path = file_path_manipulator(&PathBuf::from(file), &cwd.clone(), None, None, "");
            if file2_full_path.exists() {
                Some(file2_full_path)
            } else {
                eprintln!("File2 path does not exist: {}", file2_full_path.display());
                None
            }
        }
        None => {
            None
        }
    };

    let technology = config.args.technology.clone();
    
    let ref_db = config.args.ref_db.clone().ok_or_else(|| {
        anyhow!("HDF5 database file must be given (-d).")
    })?;
    let ref_db_path = PathBuf::from(&ref_db);

    let ercc_path = file_path_manipulator(&PathBuf::from(ERCC_FASTA), &cwd, None, None, "");
    if !ercc_path.exists() {
        return Err(anyhow!("Specified ercc {:?} does not exist.", ercc_path));
    }

    //*****************
    // Input Validation

    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(&sample_base), &cwd.clone(), None, Some("validated"), "_");
    let rx = read_and_interleave_sequences(file1_path, file2_path, Some(technology.clone()), config.args.max_reads, config.args.min_read_len, config.args.max_read_len)?;
    let val_rx_stream = ReceiverStream::new(rx);
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.args.buffer_size,
        config.args.stall_threshold.try_into().unwrap(),
        Some(config.args.stream_sleep_ms),
        50,
    )
        .await?;
    cleanup_receivers.push(val_done_rx);

    let mut streams_iter = val_streams.into_iter();
    let val_fastp_stream = streams_iter.next().unwrap();
    let val_pigz_stream = streams_iter.next().unwrap();

    //Pigz stream to intermediate file output
    let val_pigz_args = generate_cli(PIGZ_TAG, &config, None)?;
    let (mut val_pigz_child, val_pigz_stream_task, val_pigz_err_task) = stream_to_cmd(val_pigz_stream, PIGZ_TAG, val_pigz_args, StreamDataType::IlluminaFastq, config.args.verbose).await?;
    cleanup_tasks.push(val_pigz_stream_task);
    cleanup_tasks.push(val_pigz_err_task);

    let val_pigz_out_stream = parse_child_output(
        &mut val_pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.args.buffer_size,
    ).await?;
    let val_pigz_write_task = tokio::spawn(stream_to_file(
        val_pigz_out_stream,
        validated_interleaved_file_path,
    ));
    cleanup_tasks.push(val_pigz_write_task);

    // Fastp stream
    let val_fastp_args = generate_cli(FASTP_TAG, &config, None)?;
    let (mut val_fastp_child, val_fastp_stream_task, val_fastp_err_task) = stream_to_cmd(val_fastp_stream, FASTP_TAG, val_fastp_args, StreamDataType::IlluminaFastq, config.args.verbose).await?;
    cleanup_tasks.push(val_fastp_stream_task);
    cleanup_tasks.push(val_fastp_err_task);
    let val_fastp_out_stream = parse_child_output(
        &mut val_fastp_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.args.buffer_size / 4,
    ).await?;
    let val_fastp_out_stream = ReceiverStream::new(val_fastp_out_stream);

    //*****************
    // Fetch reference
    let index_start = Instant::now();
    let h5_index = get_index(&config.args).await?;
    println!("Index retrieve time: {} milliseconds.", index_start.elapsed().as_millis());

    //*****************
    // Host Removal
    let (_host_accession, host_seq) = retrieve_h5_seq(config.args.host_accession.clone(), config.args.host_sequence.clone(), Some(&ref_db_path), Some(&h5_index)).await?;
    let host_ref_temp = NamedTempFile::new_in(&ram_temp_dir)?;
    let host_ref_fasta_path = host_ref_temp.path().to_path_buf();
    temp_files.push(host_ref_temp);
    let host_ref_write_task = write_vecu8_to_file(host_seq.clone(), &host_ref_fasta_path, config.args.buffer_size).await?;
    host_ref_write_task.await?;  // This has to be done immediately to make sure the host query minimap2 can read it

    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(val_fastp_out_stream, None).await?;
    cleanup_tasks.push(host_query_write_task);
    let host_minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(host_ref_fasta_path.clone(), host_query_pipe_path.clone())))?;
    let (mut host_minimap2_child, host_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, host_minimap2_args, config.args.verbose).await?;
    cleanup_tasks.push(host_minimap2_err_task);

    let host_minimap2_out_stream = parse_child_output(
        &mut host_minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.args.buffer_size / 4,
    ).await?;

    let host_samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([("-f".to_string(), Some("4".to_string())),]),
    };
    let host_samtools_args_view = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&host_samtools_config_view),
    )?;

    let (mut host_samtools_child_view, host_samtools_task_view, host_samtools_err_task_view) = stream_to_cmd(host_minimap2_out_stream, SAMTOOLS_TAG, host_samtools_args_view, StreamDataType::JustBytes, config.args.verbose).await?;
    let host_samtools_out_stream_view = parse_child_output(
        &mut host_samtools_child_view,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.args.buffer_size / 4,
    ).await?;
    cleanup_tasks.push(host_samtools_task_view);
    cleanup_tasks.push(host_samtools_err_task_view);

    let host_samtools_config_fastq = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let host_samtools_args_fastq = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&host_samtools_config_fastq),
    )?;

    let (mut host_samtools_child_fastq, host_samtools_task_fastq, host_samtools_err_task_fastq) = stream_to_cmd(host_samtools_out_stream_view, SAMTOOLS_TAG, host_samtools_args_fastq, StreamDataType::JustBytes, config.args.verbose).await?;
    let host_samtools_out_stream_fastq = parse_child_output(
        &mut host_samtools_child_fastq,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.args.buffer_size / 4,
    ).await?;
    let host_samtools_out_stream_fastq = ReceiverStream::new(host_samtools_out_stream_fastq);
    cleanup_tasks.push(host_samtools_task_fastq);
    cleanup_tasks.push(host_samtools_err_task_fastq);

    let no_host_file_path = file_path_manipulator(&no_ext_sample_base_buf, &cwd.clone(), None, Some("no_host.fq"), "_");

    let stats_seqkit_config_stats = SeqkitConfig {
        subcommand: SeqkitSubcommand::Stats,
        subcommand_fields: HashMap::from([]),
    };
    let stats_seqkit_args_stats = generate_cli(
        SEQKIT_TAG,
        &config,
        Some(&stats_seqkit_config_stats),
    )?;

    let (host_streams, host_done_rx) = t_junction(
        host_samtools_out_stream_fastq,
        3,
        config.args.buffer_size,
        config.args.stall_threshold.try_into().unwrap(),
        Some(config.args.stream_sleep_ms),
        50,
    )
        .await?;
    cleanup_receivers.push(host_done_rx);

    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().unwrap();
    let no_host_file_stream = streams_iter.next().unwrap();
    let no_host_count_stream = streams_iter.next().unwrap();


    let (mut no_host_seqkit_child_stats, no_host_seqkit_task_stats, no_host_seqkit_err_task_stats) = stream_to_cmd(no_host_count_stream, SEQKIT_TAG, stats_seqkit_args_stats, StreamDataType::JustBytes, config.args.verbose).await?;
    let no_host_seqkit_out_stream_stats = parse_child_output(
        &mut no_host_seqkit_child_stats,
        ChildStream::Stdout,
        ParseMode::Lines,
        config.args.buffer_size / 4,
    ).await?;
    stats_tasks.push(no_host_seqkit_task_stats);
    cleanup_tasks.push(no_host_seqkit_err_task_stats);


    let host_samtools_write_task = tokio::spawn(stream_to_file(
        no_host_file_stream,
        PathBuf::from(no_host_file_path),
    ));
    stats_tasks.push(host_samtools_write_task);

    let no_host_output_stream = ReceiverStream::new(no_host_output_stream);

    // Declare file paths as Option<PathBuf> before the match
    let mut target_ref_fasta_path: Option<PathBuf> = None;
    let mut align_bam_path: Option<PathBuf> = None;
    let mut align_query_pipe_path: Option<PathBuf> = None;
    let mut consensus_file_path: Option<PathBuf> = None;
    let mut consensus_eval_stream : Option<Receiver<ParseOutput>>  = None;
    let mut consensus_bam_eval_stream : Option<Receiver<ParseOutput>>  = None;
    let mut ercc_stats_stream: Option<Receiver<ParseOutput>>  = None;
    let mut consensus_bam_stats_stream : Option<Receiver<ParseOutput>>  = None;
    let mut consensus_bam_depth_stream : Option<Receiver<ParseOutput>>  = None;
    let mut consensus_stats_stream : Option<Receiver<ParseOutput>>  = None;
    let mut call_bcftools_stats_stream : Option<Receiver<ParseOutput>>  = None;

    //*****************
    // Split by Technology
    match technology {
        Technology::Illumina => {
            eprintln!("Technology: Illumina");

            //*****************
            // Get Target reference sequence
            let (_filter_align_accession, filter_align_seq) = retrieve_h5_seq(config.args.ref_accession.clone(), config.args.ref_sequence.clone(), Some(&ref_db_path), Some(&h5_index)).await?;
            let target_ref_temp = NamedTempFile::with_suffix_in(".fasta", &config.ram_temp_dir)?;
            target_ref_fasta_path = Some(target_ref_temp.path().to_path_buf());
            temp_files.push(target_ref_temp);

            let seq_len = filter_align_seq.len();
            let seq_arc = Arc::new(filter_align_seq.clone());
            let validation_handle = if seq_len > 1_000_000 {
                tokio::spawn(validate_sequence_parallel(seq_arc.clone(), b"ACGTN", config.args.threads.min(16).max(1)))
            } else {
                tokio::spawn(async move { validate_sequence(&seq_arc, b"ACGTN") })
            };

            
            //*****************
            // ERCC
            let (ercc_streams, ercc_done_rx) = t_junction(
                no_host_output_stream,
                2,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            ).await?;
            cleanup_receivers.push(ercc_done_rx);

            let mut streams_iter = ercc_streams.into_iter();
            let ercc_stream = streams_iter.next().unwrap();
            let ercc_bypass_stream = streams_iter.next().unwrap();
            let ercc_stream = ReceiverStream::new(ercc_stream);

            let (ercc_query_write_task, ercc_query_pipe_path) = write_parse_output_to_temp(ercc_stream, None).await?;
            cleanup_tasks.push(ercc_query_write_task);
            let ercc_minimap2_args = generate_cli(
                MINIMAP2_TAG,
                &config,
                Some(&(ercc_path, ercc_query_pipe_path.clone())),
            )?;

            let (mut ercc_minimap2_child, ercc_minimap2_err_task) = spawn_cmd(
                MINIMAP2_TAG,
                ercc_minimap2_args,
                config.args.verbose,
            ).await?;
            let ercc_minimap2_out_stream = parse_child_output(
                &mut ercc_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_minimap2_err_task);

            let ercc_samtools_config_view = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::View,
                subcommand_fields: HashMap::from([]),
            };
            let ercc_samtools_args_view = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&ercc_samtools_config_view),
            )?;

            let (mut ercc_samtools_child_view, ercc_samtools_task_view, ercc_samtools_err_task_view) = stream_to_cmd(
                ercc_minimap2_out_stream,
                SAMTOOLS_TAG,
                ercc_samtools_args_view,
                StreamDataType::JustBytes,
                config.args.verbose,
            ).await?;
            let ercc_samtools_out_stream_view = parse_child_output(
                &mut ercc_samtools_child_view,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_samtools_task_view);
            cleanup_tasks.push(ercc_samtools_err_task_view);

            let ercc_samtools_config_stats = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Stats,
                subcommand_fields: HashMap::from([("-".to_string(), None)]),
            };
            let ercc_samtools_args_stats = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&ercc_samtools_config_stats),
            )?;

            let ercc_stats_file_path = no_ext_sample_base.clone() + "_stats.txt";
            let (mut ercc_samtools_child_stats, ercc_samtools_task_stats, ercc_samtools_err_task_stats) = stream_to_cmd(
                ercc_samtools_out_stream_view,
                SAMTOOLS_TAG,
                ercc_samtools_args_stats,
                StreamDataType::JustBytes,
                config.args.verbose,
            ).await?;

            let ercc_samtools_out_stream_stats = parse_child_output(
                &mut ercc_samtools_child_stats,
                ChildStream::Stdout,
                ParseMode::Lines,
                config.args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_samtools_task_stats);
            cleanup_tasks.push(ercc_samtools_err_task_stats);

            let ercc_samtools_out_stream_stats = ReceiverStream::new(ercc_samtools_out_stream_stats);
            let (ercc_streams, ercc_done_rx) = t_junction(
                ercc_samtools_out_stream_stats,
                2,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            ).await?;
            cleanup_receivers.push(ercc_done_rx);

            let mut streams_iter = ercc_streams.into_iter();
            let ercc_file_stream = streams_iter.next().unwrap();
            ercc_stats_stream = Some(streams_iter.next().unwrap());

            match ercc_stats_stream {
                Some(stream) => {
                    // Convert byte stream to line stream
                    let (line_stream, line_task) = bytes_to_lines(stream, config.args.buffer_size / 4).await?;
                    stats_tasks.push(line_task);

                    // Spawn ERCC stats parsing task
                    ercc_stats_task = Some(tokio::spawn(async move {
                        let result = parse_ercc_stats(line_stream).await;
                        if result.is_err() {
                            eprintln!("Warning: Failed to parse ERCC stats: {:?}", result);
                        }
                        result
                    }));

                    // Write ERCC stats to file
                    let ercc_stats_write_task = tokio::spawn(stream_to_file(
                        ercc_file_stream,
                        PathBuf::from(ercc_stats_file_path),
                    ));
                    cleanup_tasks.push(ercc_stats_write_task);
                }
                None => {
                    eprintln!("Warning: ercc_stats_stream is not available, skipping ERCC stats parsing and file write");
                }
            }

            // Check that the target ref is ACGTN, get handle here
            validation_handle.await?;
            let target_ref_write_task = write_vecu8_to_file(filter_align_seq.clone(), target_ref_fasta_path.as_ref().unwrap(), config.args.buffer_size).await?;
            println!("Temporary FASTA written to: {:?}", target_ref_fasta_path.as_ref().unwrap());
            target_ref_write_task.await?;  // Must make sure this file is written

            //*****************
            // Filter Reads
            let mut filter_reads_out_stream: ReceiverStream<ParseOutput>;
            if config.args.dont_filter_reads {
                filter_reads_out_stream = ReceiverStream::new(ercc_bypass_stream);
            } else {
                eprintln!("Filtering");
                let ercc_bypass_stream = ReceiverStream::new(ercc_bypass_stream);

                let(
                    filter_query_write_task,
                    filter_query_pipe_path
                ) = write_parse_output_to_temp(ercc_bypass_stream, None).await?;
                cleanup_tasks.push(filter_query_write_task);

                let filter_minimap2_args = generate_cli(
                    MINIMAP2_TAG,
                    &config,
                    Some(&(
                        target_ref_fasta_path.as_ref().unwrap().clone(),
                        filter_query_pipe_path
                    ))
                )?;
                let (mut filter_minimap2_child, filter_minimap2_err_task) = spawn_cmd(
                    MINIMAP2_TAG,
                    filter_minimap2_args,
                    config.args.verbose
                ).await?;
                let filter_minimap2_out_stream = parse_child_output(
                    &mut filter_minimap2_child,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    config.args.buffer_size / 4,
                ).await?;
                cleanup_tasks.push(filter_minimap2_err_task);

                // Sort output
                let filter_samtools_config_sort = SamtoolsConfig {
                    subcommand: SamtoolsSubcommand::Sort,
                    subcommand_fields: HashMap::from([
                        ("-n".to_string(), None),
                        ("-O".to_string(), Some("BAM".to_string())),
                        ("-".to_string(), None)
                    ]),
                };
                let filter_samtools_args_sort = generate_cli(
                    SAMTOOLS_TAG,
                    &config,
                    Some(&filter_samtools_config_sort),
                )?;
                let (
                    mut filter_samtools_child_sort,
                    filter_samtools_task_sort,
                    filter_samtools_err_task_sort
                ) = stream_to_cmd(
                    filter_minimap2_out_stream,
                    SAMTOOLS_TAG,
                    filter_samtools_args_sort,
                    StreamDataType::JustBytes,
                    config.args.verbose
                ).await?;
                cleanup_tasks.push(filter_samtools_task_sort);
                cleanup_tasks.push(filter_samtools_err_task_sort);
                let filter_samtools_out_stream_sort = parse_child_output(
                    &mut filter_samtools_child_sort,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    config.args.buffer_size / 4,
                ).await?;

                // Convert to FASTQ
                let filter_samtools_config_fastq = SamtoolsConfig {
                    subcommand: SamtoolsSubcommand::Fastq,
                    subcommand_fields: HashMap::from([("-".to_string(), None)]),
                };
                let filter_samtools_args_fastq = generate_cli(
                    SAMTOOLS_TAG,
                    &config,
                    Some(&filter_samtools_config_fastq),
                )?;
                let (
                    mut filter_samtools_child_fastq,
                    filter_samtools_task_fastq,
                    filter_samtools_err_task_fastq
                ) = stream_to_cmd(
                    filter_samtools_out_stream_sort,
                    SAMTOOLS_TAG,
                    filter_samtools_args_fastq,
                    StreamDataType::JustBytes,
                    config.args.verbose
                ).await?;
                cleanup_tasks.push(filter_samtools_task_fastq);
                cleanup_tasks.push(filter_samtools_err_task_fastq);
                let filter_samtools_out_stream_fastq = parse_child_output(
                    &mut filter_samtools_child_fastq,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    config.args.buffer_size / 4,
                ).await?;
                let filter_samtools_out_stream_fastq = ReceiverStream::new(filter_samtools_out_stream_fastq);

                // Kraken2
                let kraken2_report_path = file_path_manipulator(
                    &PathBuf::from(&no_ext_sample_base_buf),
                    &cwd.clone(),
                    None,
                    Some("kraken2_report.txt"),
                    "_"
                );
                let kraken2_classified_temp = NamedTempFile::new()?;
                let kraken2_classified_pipe_path = kraken2_classified_temp.path().to_path_buf();
                if kraken2_classified_pipe_path.exists() {
                    std::fs::remove_file(&kraken2_classified_pipe_path)?;
                }
                Command::new("mkfifo")
                    .arg(&kraken2_classified_pipe_path)
                    .status()?
                    .success()
                    .then(|| ())
                    .ok_or_else(|| anyhow!("Failed to create mkfifo for Kraken2 classified output"))?;

                let (
                    kraken2_query_write_task,
                    kraken2_query_pipe_path
                ) = write_parse_output_to_temp(filter_samtools_out_stream_fastq, None).await?;
                cleanup_tasks.push(kraken2_query_write_task);

                let filter_reads_kraken2_config = Kraken2Config {
                    report_path: kraken2_report_path,
                    classified_path: kraken2_classified_pipe_path.clone(),
                    fastq_path: kraken2_query_pipe_path.clone(),
                };

                let filter_reads_kraken2_args = generate_cli(
                    KRAKEN2_TAG,
                    &config,
                    Some(&filter_reads_kraken2_config)
                )?;
                let (_filter_kraken2_child, filter_kraken2_err_task) = spawn_cmd(
                    KRAKEN2_TAG,
                    filter_reads_kraken2_args,
                    config.args.verbose
                ).await?;
                cleanup_tasks.push(filter_kraken2_err_task);

                let kraken2_classified_stream = TokioFile::open(&kraken2_classified_pipe_path).await?;
                let parse_rx = parse_fastq(kraken2_classified_stream, config.args.buffer_size).await?;
                let filter_fn = |id: &str| id.contains("kraken:taxid|2697049");

                let (filtered_rx, filter_task) = parse_and_filter_fastq_id(parse_rx, config.args.buffer_size, filter_fn.clone());
                cleanup_tasks.push(filter_task);

                let (parse_output_tx, parse_output_rx) = mpsc::channel(config.args.buffer_size);
                tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(filtered_rx);
                    while let Some(record) = stream.next().await {
                        let bytes = record.to_bytes()?;
                        if parse_output_tx.send(ParseOutput::Bytes(bytes)).await.is_err() {
                            eprintln!("Failed to send ParseOutput::Bytes");
                            break;
                        }
                    }
                    Ok::<(), anyhow::Error>(())
                });

                filter_reads_out_stream = ReceiverStream::new(parse_output_rx);
            }


            let align_fastq_path = file_path_manipulator(
                &PathBuf::from(&no_ext_sample_base_buf),
                &cwd.clone(),
                None,
                Some("filtered.fq"),
                "_"
            );

            let (align_streams, align_done_rx) = t_junction(
                filter_reads_out_stream,
                2,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            )
                .await?;
            cleanup_receivers.push(align_done_rx);

            let mut streams_iter = align_streams.into_iter();
            let align_query_stream = streams_iter.next().unwrap();
            let align_file_stream = streams_iter.next().unwrap();
            let align_query_stream = ReceiverStream::new(align_query_stream);

            let align_file_write_task = tokio::spawn(stream_to_file(
                align_file_stream,
                align_fastq_path.clone(),
            ));
            cleanup_tasks.push(align_file_write_task);


            //*****************
            // Align Reads to Target
            let (align_query_write_task, align_query_pipe_path_temp) = write_parse_output_to_temp(
                align_query_stream,
                None
            ).await?;
            align_query_pipe_path = Some(align_query_pipe_path_temp.clone());
            quast_write_tasks.push(align_query_write_task);

            let align_minimap2_args = generate_cli(
                MINIMAP2_TAG,
                &config,
                Some(&(
                    target_ref_fasta_path.as_ref().unwrap().clone(),
                    align_query_pipe_path_temp
                ))
            )?;
            let (mut align_minimap2_child, align_minimap2_err_task) = spawn_cmd(
                MINIMAP2_TAG,
                align_minimap2_args,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(align_minimap2_err_task);
            let align_minimap2_out_stream = parse_child_output(
                &mut align_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;

            let align_samtools_config_sort = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Sort,
                subcommand_fields: HashMap::from([
                    ("-O".to_string(), Some("sam".to_string())),
                    ("-".to_string(), None)
                ]),
            };
            let align_samtools_args_sort = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&align_samtools_config_sort),
            )?;
            let (
                mut align_samtools_child_sort,
                align_samtools_task_sort,
                align_samtools_err_task_sort
            ) = stream_to_cmd(
                align_minimap2_out_stream,
                SAMTOOLS_TAG,
                align_samtools_args_sort,
                StreamDataType::JustBytes,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(align_samtools_task_sort);
            cleanup_tasks.push(align_samtools_err_task_sort);
            let align_samtools_out_stream_sort = parse_child_output(
                &mut align_samtools_child_sort,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;

            //*****************
            // Make Consensus
            align_bam_path = Some(file_path_manipulator(
                &no_ext_sample_base_buf,
                &cwd.clone(),
                None,
                Some("align.sam"),
                "_"
            ));
            let align_samtools_out_stream_sort = ReceiverStream::new(align_samtools_out_stream_sort);

            let (consensus_bam_streams, consensus_bam_done_rx) = t_junction(
                align_samtools_out_stream_sort,
                6,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            )
                .await?;
            cleanup_receivers.push(consensus_bam_done_rx);

            let mut streams_iter = consensus_bam_streams.into_iter();
            let consensus_bam_output_stream = streams_iter.next().unwrap();
            let consensus_bam_file_stream = streams_iter.next().unwrap();
            let consensus_bam_call_stream = streams_iter.next().unwrap();
            consensus_bam_eval_stream = Some(streams_iter.next().unwrap());
            consensus_bam_stats_stream = Some(streams_iter.next().unwrap());
            consensus_bam_depth_stream = Some(streams_iter.next().unwrap());

            let consensus_bam_write_task = tokio::spawn(stream_to_file(
                consensus_bam_file_stream,
                align_bam_path.as_ref().unwrap().clone(),
            ));
            quast_write_tasks.push(consensus_bam_write_task);

            let consensus_samtools_config = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Consensus,
                subcommand_fields: HashMap::from([("-".to_string(), None)]),
            };
            let consensus_samtools_args = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&consensus_samtools_config),
            )?;
            consensus_file_path = Some(file_path_manipulator(
                &PathBuf::from(&no_ext_sample_base_buf),
                &cwd.clone(),
                None,
                Some("consensus.fa"),
                "_"
            ));
            let (
                mut consensus_samtools_child,
                consensus_samtools_task_sort,
                consensus_samtools_err_task_sort
            ) = stream_to_cmd(
                consensus_bam_output_stream,
                SAMTOOLS_TAG,
                consensus_samtools_args,
                StreamDataType::JustBytes,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(consensus_samtools_task_sort);
            cleanup_tasks.push(consensus_samtools_err_task_sort);
            let consensus_samtools_out_stream = parse_child_output(
                &mut consensus_samtools_child,
                ChildStream::Stdout,
                ParseMode::Fasta,
                config.args.buffer_size / 4,
            ).await?;

            let consensus_samtools_out_stream = ReceiverStream::new(consensus_samtools_out_stream);

            let (consensus_streams, consensus_done_rx) = t_junction(
                consensus_samtools_out_stream,
                4,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            )
                .await?;
            cleanup_receivers.push(consensus_done_rx);
            let mut streams_iter = consensus_streams.into_iter();
            let consensus_file_stream = streams_iter.next().unwrap();
            let consensus_realign_stream = streams_iter.next().unwrap();
            consensus_eval_stream = Some(streams_iter.next().unwrap());
            consensus_stats_stream = Some(streams_iter.next().unwrap());


            let consensus_fa_write_task = tokio::spawn(stream_to_file(
                consensus_file_stream,
                consensus_file_path.as_ref().unwrap().clone()
            ));
            quast_write_tasks.push(consensus_fa_write_task);

            //*****************
            // Call Variants
            let call_bcftools_config_mpileup = BcftoolsConfig {
                subcommand: BcftoolsSubcommand::Mpileup,
                subcommand_fields: HashMap::from([
                    (
                        "-f".to_string(),
                        Some(target_ref_fasta_path.as_ref().unwrap().to_string_lossy().into_owned())
                    ),
                    ("-".to_string(), None),
                ])
            };
            let call_bcftools_args_mpileup = generate_cli(
                BCFTOOLS_TAG,
                &config,
                Some(&call_bcftools_config_mpileup),
            )?;
            let (
                mut call_bcftools_child_mpileup,
                call_bcftools_task_mpileup,
                call_bcftools_err_task_mpileup
            ) = stream_to_cmd(
                consensus_bam_call_stream,
                BCFTOOLS_TAG,
                call_bcftools_args_mpileup,
                StreamDataType::JustBytes,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(call_bcftools_task_mpileup);
            cleanup_tasks.push(call_bcftools_err_task_mpileup);
            let call_bcftools_out_stream_mpileup = parse_child_output(
                &mut call_bcftools_child_mpileup,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;

            let call_bcftools_config_call = BcftoolsConfig {
                subcommand: BcftoolsSubcommand::Call,
                subcommand_fields: HashMap::from([
                    ("--ploidy".to_string(), Some("1".to_string())),
                    ("-m".to_string(), None),
                    ("-v".to_string(), None),
                    ("-P".to_string(), Some(config.args.bcftools_call_theta.to_string())),
                    ("-".to_string(), None),
                ])
            };
            let call_bcftools_args_call = generate_cli(
                BCFTOOLS_TAG,
                &config,
                Some(&call_bcftools_config_call),
            )?;

            let (
                mut call_bcftools_child_call,
                call_bcftools_task_call,
                call_bcftools_err_task_call
            ) = stream_to_cmd(
                call_bcftools_out_stream_mpileup,
                BCFTOOLS_TAG,
                call_bcftools_args_call,
                StreamDataType::JustBytes,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(call_bcftools_task_call);
            cleanup_tasks.push(call_bcftools_err_task_call);
            let call_bcftools_out_stream_call = parse_child_output(
                &mut call_bcftools_child_call,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;

            let called_variants_path = file_path_manipulator(
                &no_ext_sample_base_buf,
                &cwd.clone(),
                None,
                Some("called.vcf"),
                "_"
            );
            let call_bcftools_config_view = BcftoolsConfig {
                subcommand: BcftoolsSubcommand::View,
                subcommand_fields: HashMap::from([
                    ("-i".to_string(), Some(format!("DP>={}", config.args.min_depth))),
                    ("-O".to_string(), Some("v".to_string())),
                    ("-".to_string(), None),
                ])
            };
            let call_bcftools_args_view = generate_cli(
                BCFTOOLS_TAG,
                &config,
                Some(&call_bcftools_config_view),
            )?;

            let (
                mut call_bcftools_child_view,
                call_bcftools_task_view,
                call_bcftools_err_task_view
            ) = stream_to_cmd(
                call_bcftools_out_stream_call,
                BCFTOOLS_TAG,
                call_bcftools_args_view,
                StreamDataType::JustBytes,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(call_bcftools_task_view);
            cleanup_tasks.push(call_bcftools_err_task_view);

            let call_bcftools_out_stream_view = parse_child_output(
                &mut call_bcftools_child_view,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;

            let call_bcftools_out_stream_view = ReceiverStream::new(call_bcftools_out_stream_view);
            let (call_bcftools_out_streams, call_bcftools_done_rx) = t_junction(
                call_bcftools_out_stream_view,
                2,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            )
                .await?;
            cleanup_receivers.push(call_bcftools_done_rx);

            let mut streams_iter = call_bcftools_out_streams.into_iter();
            let call_bcftools_file_stream = streams_iter.next().unwrap();
            call_bcftools_stats_stream  = Some(streams_iter.next().unwrap());
            let called_variants_write_task = tokio::spawn(stream_to_file(
                call_bcftools_file_stream,
                called_variants_path
            ));
            cleanup_tasks.push(called_variants_write_task);

            //*****************
            // Realign Consensus
            let reference_file = TokioFile::open(target_ref_fasta_path.as_ref().unwrap()).await?;
            let reference_rx = parse_bytes(reference_file, config.args.buffer_size).await?;

            let realign_streams = vec![consensus_realign_stream, reference_rx];
            let (combined_rx, combined_task) = y_junction(realign_streams, config.args.buffer_size).await?;
            cleanup_tasks.push(combined_task);

            let realign_consensus_path = file_path_manipulator(
                &no_ext_sample_base_buf,
                &cwd.clone(),
                None,
                Some("consensus_realigned.fa"),
                "_"
            );
            let realign_mafft_args = generate_cli(MAFFT_TAG, &config, None)?;
            let (
                mut realign_consensus_mafft_child,
                realign_consensus_mafft_task,
                realign_consensus_mafft_err_task
            ) = stream_to_cmd(
                combined_rx,
                MAFFT_TAG,
                realign_mafft_args,
                StreamDataType::JustBytes,
                config.args.verbose
            ).await?;
            cleanup_tasks.push(realign_consensus_mafft_task);
            cleanup_tasks.push(realign_consensus_mafft_err_task);
            let consensus_samtools_out_stream = parse_child_output(
                &mut realign_consensus_mafft_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;

            let realign_consensus_write_task = tokio::spawn(stream_to_file(
                consensus_samtools_out_stream,
                PathBuf::from(realign_consensus_path),
            ));
            cleanup_tasks.push(realign_consensus_write_task);


        } // end tech == illumina
        Technology::ONT => {
            target_ref_fasta_path = None;
            align_bam_path = None;
            align_query_pipe_path = None;
            consensus_file_path = None;
            return Err(anyhow!("Minion not ready"));
        }
    }


    //*****************
    // Calculate Statistics
    let results = try_join_all(stats_tasks).await?;
    for result in results {
        result?;
    }

    let stats_samtools_config_stats = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Stats,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let stats_samtools_args_stats = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&stats_samtools_config_stats),
    )?;

    let stats_samtools_out_stream_stats = match consensus_bam_stats_stream {
        Some(stream) => {
            let (mut stats_samtools_child_stats, stats_samtools_task_stats, stats_samtools_err_task_stats) = stream_to_cmd(
                stream,
                SAMTOOLS_TAG,
                stats_samtools_args_stats.clone(),
                StreamDataType::JustBytes,
                config.args.verbose,
            ).await?;
            cleanup_tasks.push(stats_samtools_task_stats);
            cleanup_tasks.push(stats_samtools_err_task_stats);
            parse_child_output(&mut stats_samtools_child_stats, ChildStream::Stdout, ParseMode::Lines, config.args.buffer_size / 4).await?
        }
        None => return Err(anyhow!("consensus_bam_stats_stream is not available")),
    };
    let samtools_stats_out = parse_samtools_stats(stats_samtools_out_stream_stats).await?;

    let depth_samtools_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Depth,
        subcommand_fields: HashMap::from([
            ("-aa".to_string(), None),
            ("-d".to_string(), Some("0".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let depth_samtools_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&depth_samtools_config),
    )?;

    let depth_samtools_out_stream = match consensus_bam_depth_stream {
        Some(stream) => {
            let (mut depth_samtools_child, depth_samtools_task, depth_samtools_err_task) = stream_to_cmd(
                stream,
                SAMTOOLS_TAG,
                depth_samtools_args,
                StreamDataType::JustBytes,
                config.args.verbose,
            ).await?;
            cleanup_tasks.push(depth_samtools_task);
            cleanup_tasks.push(depth_samtools_err_task);
            parse_child_output(
                &mut depth_samtools_child,
                ChildStream::Stdout,
                ParseMode::Lines,
                config.args.buffer_size / 4,
            ).await?
        }
        None => {
            return Err(anyhow!("consensus_bam_stats_stream is not available"));
        }
    };

    let depth_map = parse_samtools_depth(depth_samtools_out_stream).await?;
    if depth_map.is_empty() {
        return Err(anyhow!("No depth data found"));
    }
    // NB: for this pipeline, for now at least, there should only be one chrom, so takingt he first one
    let first_chr = depth_map.keys().next().ok_or_else(|| anyhow!("No chromosomes found"))?.clone();
    let first_chr_depth_map = depth_map.get(&first_chr).unwrap();
    let depths: Vec<u32> = first_chr_depth_map.values().copied().collect();
    let samtools_depth_stats = compute_depth_stats(&depths)?;

    let depth_plot_path = file_path_manipulator(
        &PathBuf::from(&no_ext_sample_base_buf),
        &cwd.clone(),
        None,
        Some("depth.png"),
        "_"
    );

    plot_depths(&first_chr_depth_map, &no_ext_sample_base, &depth_plot_path)?;

    let seqkit_stats = parse_seqkit_stats(no_host_seqkit_out_stream_stats).await?;

    let ercc_stats = if let Technology::Illumina = technology {
        match ercc_stats_task {
            Some(task) => match task.await {
                Ok(Ok(stats)) => stats,
                Ok(Err(e)) => return Err(anyhow!("Failed to parse ERCC stats: {}", e)),
                Err(e) => return Err(anyhow!("ERCC stats task failed: {}", e)),
            },
            None => return Err(anyhow!("ERCC stats task not initialized for Illumina technology")),
        }
    } else {
        HashMap::new() // Empty stats for ONT
    };

    let allele_counts = if let Some(stream) = consensus_stats_stream {
        compute_allele_counts(stream).await?
    } else {
        HashMap::new()
    };

    let (ref_snps, ref_mnps, ref_indels) = if let Some(stream) = call_bcftools_stats_stream {
        parse_vcf_stream(stream).await?
    } else {
        (0, 0, 0)
    };

    let n_actg = allele_counts.iter().filter(|&(k, _)| "ACTGU".contains(*k)).map(|(_, &v)| v).sum::<u64>();
    let n_missing = allele_counts.get(&'N').copied().unwrap_or(0);
    let n_gap = allele_counts.get(&'-').copied().unwrap_or(0);
    let n_ambiguous = allele_counts.iter().filter(|&(k, _)| !"ACTGUN-".contains(*k)).map(|(_, &v)| v).sum::<u64>();

    let coverage_breadth = if !depths.is_empty() { depths.iter().filter(|&&d| d > 0).count() as f64 / depths.len() as f64 } else { 0.0 };
    let max_aligned_length = depths.len();
    let total_length = depths.len();
    let (coverage_bin_size, coverage) = if !depths.is_empty() { compute_coverage_bins(&depths, 500) } else { (0.0, Vec::new()) };

    let stats = Stats {
        sample_name: no_ext_sample_base.clone(),
        depth_avg: samtools_depth_stats.get("depth_avg").copied().unwrap_or(0.0),
        depth_q25: samtools_depth_stats.get("depth_q.25").copied().unwrap_or(0.0),
        depth_q50: samtools_depth_stats.get("depth_q.5").copied().unwrap_or(0.0),
        depth_q75: samtools_depth_stats.get("depth_q.75").copied().unwrap_or(0.0),
        depth_frac_above_10x: samtools_depth_stats.get("depth_frac_above_10x").copied().unwrap_or(0.0),
        depth_frac_above_25x: samtools_depth_stats.get("depth_frac_above_25x").copied().unwrap_or(0.0),
        depth_frac_above_50x: samtools_depth_stats.get("depth_frac_above_50x").copied().unwrap_or(0.0),
        depth_frac_above_100x: samtools_depth_stats.get("depth_frac_above_100x").copied().unwrap_or(0.0),
        allele_counts,
        total_reads: seqkit_stats.get("num_seqs").and_then(|s| s.parse::<u64>().ok()).unwrap_or(0),
        mapped_reads: samtools_stats_out.get("reads mapped").and_then(|s| s.parse::<u64>().ok()).unwrap_or(0),
        mapped_paired: samtools_stats_out.get("reads mapped and paired").and_then(|s| s.parse::<u64>().ok()),
        paired_inward: samtools_stats_out.get("inward oriented pairs").and_then(|s| s.parse::<u64>().ok()).map(|v| v * 2),
        paired_outward: samtools_stats_out.get("outward oriented pairs").and_then(|s| s.parse::<u64>().ok()).map(|v| v * 2),
        paired_other_orientation: samtools_stats_out.get("pairs with other orientation").and_then(|s| s.parse::<u64>().ok()).map(|v| v * 2),
        ercc_mapped_reads: ercc_stats.get("ercc_mapped_reads").copied(),
        ercc_mapped_paired: ercc_stats.get("ercc_mapped_paired").copied(),
        ref_snps,
        ref_mnps,
        ref_indels,
        n_actg,
        n_missing,
        n_gap,
        n_ambiguous,
        coverage_breadth,
        max_aligned_length,
        total_length,
        coverage_bin_size,
        coverage,
    };

    let stats_file_path = cwd.join(format!("{}_stats.json", no_ext_sample_base));
    let mut stats_file = File::create(&stats_file_path)?;
    serde_json::to_writer_pretty(&mut stats_file, &stats)?;


    //*****************
    // Assembly Evaluation

    // Makes sure the needed files for Quast are written
    let results = try_join_all(quast_write_tasks).await?;
    for result in results {
        result?;
    }

    if [
        target_ref_fasta_path.as_ref(),
        align_bam_path.as_ref(),
        consensus_file_path.as_ref()
    ].iter().any(|opt| opt.is_none()) {
        return Err(anyhow!("One or more required paths are not set"));
    }

    let quast_config = QuastConfig {
        ref_fasta: target_ref_fasta_path.unwrap().to_string_lossy().into_owned(),
        ref_bam: align_bam_path.unwrap().to_string_lossy().into_owned(),
        assembly_fasta: consensus_file_path.unwrap().to_string_lossy().into_owned(),
    };

    let assembly_eval_quast_args = generate_cli(
        QUAST_TAG,
        &config,
        Some(&quast_config),
    )?;

    let (_assembly_eval_quast_child, assembly_eval_quast_err_task) = spawn_cmd(QUAST_TAG, assembly_eval_quast_args, config.args.verbose).await?;
    cleanup_tasks.push(assembly_eval_quast_err_task);


    //*****************
    // Cleanup, hanging tasks
    let results = try_join_all(cleanup_tasks).await?;
    for result in results {
        result?;
    }
    for receiver in cleanup_receivers {
        receiver.await??;
    }
    drop(temp_files);
    println!("Finished generating consensus genome");

    Ok(())
}