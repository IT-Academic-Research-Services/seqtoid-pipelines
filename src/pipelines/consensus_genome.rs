use std::fs;
use std::collections::HashMap;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::{Path, PathBuf};
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::Technology;
use std::process::Command;
use std::time::Instant;
use tokio::fs::File as TokioFile;
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::mpsc::Receiver;
use crate::utils::command::{generate_cli, check_versions};
use crate::utils::file::{extension_remover, file_path_manipulator, write_parse_output_to_temp, write_vecu8_to_file};
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, parse_and_filter_fastq_id, validate_sequence, CONTIG_THRESHOLDS, SequenceRecord};
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd, parse_fastq, parse_bytes, y_junction, convert_stream};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand, KRAKEN2_TAG, BCFTOOLS_TAG, BcftoolsSubcommand, MAFFT_TAG, NUCMER_TAG, NUCMER_DELTA, SHOW_COORDS_TAG};
use crate::utils::sequence::valid_bases::DNA_WITH_N;
use crate::utils::metrics::write_metrics_to_tsv;
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::command::bcftools::BcftoolsConfig;
use crate::utils::command::nucmer::NucmerConfig;
use crate::utils::db::{get_index, retrieve_h5_seq};
use tokio::sync::mpsc;
use futures::future::try_join_all;
use crate::utils::streams::ToBytes;
use crate::config::defs::RunConfig;

const ERCC_FASTA: &str = "ercc_sequences.fasta";

pub async fn run(config: &RunConfig) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", config.args.module);

    let cwd = config.cwd.clone();
    let ram_temp_dir = config.ram_temp_dir.clone();
    let mut temp_files = Vec::new();

    // Initialize cleanup tasks
    let mut cleanup_tasks: Vec<tokio::task::JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<tokio::sync::oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();

    // External tools check
    let _tool_versions = check_versions(vec![SAMTOOLS_TAG, MINIMAP2_TAG, FASTP_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, NUCMER_TAG]).await?;

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
    let val_pigz_args = generate_cli(PIGZ_TAG, &config.args, None)?;
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
    let val_fastp_args = generate_cli(FASTP_TAG, &config.args, None)?;
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
    cleanup_tasks.push(host_ref_write_task);

    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(val_fastp_out_stream, None).await?;
    cleanup_tasks.push(host_query_write_task);
    let host_minimap2_args = generate_cli(MINIMAP2_TAG, &config.args, Some(&(host_ref_fasta_path.clone(), host_query_pipe_path.clone())))?;
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
        &config.args,
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
        &config.args,
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

    let no_host_file_path = file_path_manipulator(&PathBuf::from(&sample_base), &cwd.clone(), None, Some("no_host"), "_");

    let (host_streams, host_done_rx) = t_junction(
        host_samtools_out_stream_fastq,
        2,
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

    let host_samtools_write_task = tokio::spawn(stream_to_file(
        no_host_file_stream,
        PathBuf::from(no_host_file_path),
    ));
    cleanup_tasks.push(host_samtools_write_task);

    let no_host_output_stream = ReceiverStream::new(no_host_output_stream);

    // Declare file paths as Option<PathBuf> before the match
    let mut target_ref_fasta_path: Option<PathBuf> = None;
    let mut align_bam_path: Option<PathBuf> = None;
    let mut align_query_pipe_path: Option<PathBuf> = None;
    let mut consensus_file_path: Option<PathBuf> = None;
    let mut consensus_eval_stream : Option<Receiver<ParseOutput>>  = None;
    let mut consensus_bam_eval_stream : Option<Receiver<ParseOutput>>  = None;

    //*****************
    // Split by Technology
    match technology {
        Technology::Illumina => {
            eprintln!("Illumina");

            //*****************
            // Get Target reference sequence
            let (_filter_align_accession, filter_align_seq) = retrieve_h5_seq(config.args.ref_accession.clone(), config.args.ref_sequence.clone(), Some(&ref_db_path), Some(&h5_index)).await?;
            let target_ref_temp = NamedTempFile::with_suffix_in(".fasta", &config.ram_temp_dir)?;
            target_ref_fasta_path = Some(target_ref_temp.path().to_path_buf());
            temp_files.push(target_ref_temp);
            let target_ref_write_task = write_vecu8_to_file(filter_align_seq.clone(), target_ref_fasta_path.as_ref().unwrap(), config.args.buffer_size).await?;
            cleanup_tasks.push(target_ref_write_task);

            //*****************
            // ERCC
            let (ercc_streams, ercc_done_rx) = t_junction(
                no_host_output_stream,
                2,
                config.args.buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
            )
                .await?;
            cleanup_receivers.push(ercc_done_rx);

            let mut streams_iter = ercc_streams.into_iter();
            let ercc_stream = streams_iter.next().unwrap();
            let ercc_bypass_stream = streams_iter.next().unwrap();
            let ercc_stream = ReceiverStream::new(ercc_stream);

            let (ercc_query_write_task, ercc_query_pipe_path) = write_parse_output_to_temp(ercc_stream, None).await?;
            cleanup_tasks.push(ercc_query_write_task);
            let ercc_minimap2_args = generate_cli(MINIMAP2_TAG, &config.args, Some(&(ercc_path, ercc_query_pipe_path.clone())))?;

            let (mut ercc_minimap2_child, ercc_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, ercc_minimap2_args, config.args.verbose).await?;
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
                &config.args,
                Some(&ercc_samtools_config_view),
            )?;

            let (mut ercc_samtools_child_view, ercc_samtools_task_view, ercc_samtools_err_task_view) = stream_to_cmd(ercc_minimap2_out_stream, SAMTOOLS_TAG, ercc_samtools_args_view, StreamDataType::JustBytes, config.args.verbose).await?;
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
                &config.args,
                Some(&ercc_samtools_config_stats),
            )?;

            let ercc_stats_file_path = no_ext_sample_base.clone() + "_stats.txt";
            let (mut ercc_samtools_child_stats, ercc_samtools_task_stats, ercc_samtools_err_task_stats) = stream_to_cmd(ercc_samtools_out_stream_view, SAMTOOLS_TAG, ercc_samtools_args_stats, StreamDataType::JustBytes, config.args.verbose).await?;

            let ercc_samtools_out_stream_stats = parse_child_output(
                &mut ercc_samtools_child_stats,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_samtools_task_stats);
            cleanup_tasks.push(ercc_samtools_err_task_stats);

            let ercc_stats_write_task = tokio::spawn(stream_to_file(
                ercc_samtools_out_stream_stats,
                PathBuf::from(ercc_stats_file_path),
            ));
            cleanup_tasks.push(ercc_stats_write_task);

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
                    &config.args,
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
                    &config.args,
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
                    &config.args,
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
                    &config.args,
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
            cleanup_tasks.push(align_query_write_task);

            let align_minimap2_args = generate_cli(
                MINIMAP2_TAG,
                &config.args,
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
                    ("-O".to_string(), Some("bam".to_string())),
                    ("-".to_string(), None)
                ]),
            };
            let align_samtools_args_sort = generate_cli(
                SAMTOOLS_TAG,
                &config.args,
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
                Some("align.bam"),
                "_"
            ));
            let align_samtools_out_stream_sort = ReceiverStream::new(align_samtools_out_stream_sort);

            let (consensus_bam_streams, consensus_bam_done_rx) = t_junction(
                align_samtools_out_stream_sort,
                4,
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


            let consensus_bam_write_task = tokio::spawn(stream_to_file(
                consensus_bam_file_stream,
                align_bam_path.as_ref().unwrap().clone(),
            ));
            cleanup_tasks.push(consensus_bam_write_task);

            let consensus_samtools_config = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Consensus,
                subcommand_fields: HashMap::from([("-".to_string(), None)]),
            };
            let consensus_samtools_args = generate_cli(
                SAMTOOLS_TAG,
                &config.args,
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
                3,
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


            let consensus_fa_write_task = tokio::spawn(stream_to_file(
                consensus_file_stream,
                consensus_file_path.as_ref().unwrap().clone()
            ));
            cleanup_tasks.push(consensus_fa_write_task);

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
                &config.args,
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
                &config.args,
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
                &config.args,
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

            let called_variants_write_task = tokio::spawn(stream_to_file(
                call_bcftools_out_stream_view,
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
            let realign_mafft_args = generate_cli(MAFFT_TAG, &config.args, None)?;
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
    // Assembly Evaluation

    let nucmer_delta_buf:  PathBuf = PathBuf::from(&NUCMER_DELTA);
    if nucmer_delta_buf.clone().exists() {
        fs::remove_file(nucmer_delta_buf.clone())?;
    }

    let assembly_eval_nucmer_config = NucmerConfig {
        ref_fasta: target_ref_fasta_path.as_ref().unwrap().clone(),
        assembly_fasta: consensus_file_path.as_ref().unwrap().clone(),
    };

    let assembly_eval_args = generate_cli(
        NUCMER_TAG,
        &config.args,
        Some(&assembly_eval_nucmer_config),
    )?;

    let (_assembly_eval_nucmer_child, assembly_eval_nucmer_err_task) = spawn_cmd(NUCMER_TAG, assembly_eval_args, config.args.verbose).await?;
    assembly_eval_nucmer_err_task.await?;

    if !nucmer_delta_buf.clone().exists() {
        return Err(anyhow!("File alignment.delta missing after nucmer run."));
    }

    let assembly_show_coords_args = generate_cli(
        SHOW_COORDS_TAG,
        &config.args,
        None,
    )?;

    let (mut assembly_show_coords_child, assembly_show_coords_err_task) = spawn_cmd(SHOW_COORDS_TAG, assembly_show_coords_args, config.args.verbose).await?;
    let assembly_show_coords_out_stream = parse_child_output(
        &mut assembly_show_coords_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.args.buffer_size / 4,
    ).await?;
    cleanup_tasks.push(assembly_show_coords_err_task);


    let assembly_eval_report_path = file_path_manipulator(
        &no_ext_sample_base_buf,
        &cwd.clone(),
        None,
        Some("assembly_eval_report.tsv"),
        "_"
    );


    let consensus_eval_stream = consensus_eval_stream.ok_or_else(|| anyhow!("Consensus eval stream is None"))?;
    let (fasta_rx, convert_task) = convert_stream(consensus_eval_stream, config.args.buffer_size).await?;
    cleanup_tasks.push(convert_task);

    let consensus_bam_eval_stream = consensus_bam_eval_stream.ok_or_else(|| anyhow!("Consensus BAM eval stream is None"))?;

    write_metrics_to_tsv(fasta_rx,
                         &target_ref_fasta_path.as_ref().unwrap().clone(),
                         consensus_bam_eval_stream,
                         assembly_show_coords_out_stream,
                         &nucmer_delta_buf.clone(),
                         &CONTIG_THRESHOLDS,
                         &assembly_eval_report_path).await?;



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