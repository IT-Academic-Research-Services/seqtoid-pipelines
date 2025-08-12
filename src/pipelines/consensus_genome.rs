use std::io;
use std::fs::File;
use std::sync::Arc;
use std::io::Write;
use std::collections::HashMap;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::Technology;
use std::process::Command;
use tokio::task::JoinHandle;
use std::time::Instant;
use tokio::fs::File as TokioFile;
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::mpsc::Receiver;
use tokio::sync::Notify;
use serde::Serialize;
use crate::utils::command::{generate_cli, check_versions};
use crate::utils::file::{extension_remover, file_path_manipulator, write_parse_output_to_temp, write_vecu8_to_file};
use crate::utils::sambam::stream_sam_alignment_counter;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, parse_and_filter_fastq_id, validate_sequence, validate_sequence_parallel, SequenceRecord, parse_header, stream_record_counter};
use crate::utils::streams::{t_junction, stream_to_cmd, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd, parse_fastq, parse_bytes, y_junction, bytes_to_lines};
use crate::config::defs::{PipelineError, StreamDataType, PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand, KRAKEN2_TAG, BCFTOOLS_TAG, BcftoolsSubcommand, MAFFT_TAG, QUAST_TAG, SEQKIT_TAG, SeqkitSubcommand};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::command::bcftools::BcftoolsConfig;
use crate::utils::command::seqkit::SeqkitConfig;
use crate::utils::db::{get_index, retrieve_h5_seq};
use tokio::sync::{mpsc, oneshot};
use futures::future::try_join_all;
use fxhash::FxHashMap as FxHashMap;
use crate::utils::streams::ToBytes;
use crate::config::defs::RunConfig;
use crate::utils::command::quast::QuastConfig;
use crate::utils::stats::{parse_samtools_stats, parse_samtools_depth, compute_depth_stats, parse_seqkit_stats, parse_ercc_stats, compute_allele_counts, compute_coverage_bins};
use crate::utils::vcf::parse_vcf_stream;
use crate::utils::plotting::plot_depths;


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


async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {

    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated"),
        "_",
    );

    let rx = read_and_interleave_sequences(
        file1_path,
        file2_path,
        Some(config.args.technology.clone()),
        config.args.max_reads,
        config.args.min_read_len,
        config.args.max_read_len,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    if val_streams.is_empty() {
        return Err(PipelineError::EmptyStream);
    }

    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = vec![val_done_rx];
    let mut streams_iter = val_streams.into_iter();
    let val_fastp_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let val_pigz_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let val_pigz_args = generate_cli(PIGZ_TAG, &config, None)
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut val_pigz_child, val_pigz_stream_task, val_pigz_err_task) = stream_to_cmd(
        config.clone(),
        val_pigz_stream,
        PIGZ_TAG,
        val_pigz_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(val_pigz_stream_task);
    cleanup_tasks.push(val_pigz_err_task);

    let val_pigz_out_stream = parse_child_output(
        &mut val_pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let val_pigz_write_task = tokio::spawn(stream_to_file(val_pigz_out_stream, validated_interleaved_file_path));
    cleanup_tasks.push(val_pigz_write_task);

    let val_fastp_args = generate_cli(FASTP_TAG, &config, None)
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut val_fastp_child, val_fastp_stream_task, val_fastp_err_task) = stream_to_cmd(
        config.clone(),
        val_fastp_stream,
        FASTP_TAG,
        val_fastp_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(val_fastp_stream_task);
    cleanup_tasks.push(val_fastp_err_task);

    let val_fastp_out_stream = parse_child_output(
        &mut val_fastp_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;
    // let val_fastp_out_stream = ReceiverStream::new(val_fastp_out_stream);

    Ok((ReceiverStream::new(val_fastp_out_stream), cleanup_tasks, cleanup_receivers))
}


async fn fetch_host_reference<'a>(
    config: &RunConfig,
    ref_db_path: Option<PathBuf>,
    ram_temp_dir: &PathBuf,
    h5_index: Option<&'a FxHashMap<[u8; 24], u64>>,
) -> Result<(PathBuf, NamedTempFile,JoinHandle<Result<(), anyhow::Error>>), PipelineError> {
    let (_host_accession, host_seq) = retrieve_h5_seq(
        config.args.host_accession.clone(),
        config.args.host_sequence.clone(),
        ref_db_path.as_ref(),
        h5_index,
    )
        .await
        .map_err(|e| PipelineError::ReferenceRetrievalFailed(e.to_string()))?;

    let host_ref_temp = NamedTempFile::new_in(ram_temp_dir)
        .map_err(|e| PipelineError::Other(e.into()))?;
    let host_ref_fasta_path = host_ref_temp.path().to_path_buf();

    let host_ref_write_task = write_vecu8_to_file(host_seq.clone(), &host_ref_fasta_path, config.base_buffer_size)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok((host_ref_fasta_path, host_ref_temp, host_ref_write_task))
}

async fn fetch_target_reference<'a>(
    config: &RunConfig,
    ref_db_path: Option<PathBuf>,
    ram_temp_dir: &PathBuf,
    h5_index: Option<&'a FxHashMap<[u8; 24], u64>>,
) -> Result<(PathBuf, NamedTempFile,JoinHandle<Result<(), anyhow::Error>>), PipelineError> {

    let (_target_accession, target_seq) = retrieve_h5_seq(
        config.args.ref_accession.clone(),
        config.args.ref_sequence.clone(),
        ref_db_path.as_ref(),
        h5_index,
    )
        .await
        .map_err(|e| PipelineError::ReferenceRetrievalFailed(e.to_string()))?;

    let target_ref_temp = NamedTempFile::with_suffix_in(".fasta", ram_temp_dir)
        .map_err(|e| PipelineError::Other(e.into()))?;
    let target_ref_fasta_path = target_ref_temp.path().to_path_buf();
    let target_ref_write_task = write_vecu8_to_file(target_seq.clone(), &target_ref_fasta_path, config.base_buffer_size)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok((target_ref_fasta_path, target_ref_temp,  target_ref_write_task))
}

async fn align_to_host(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    host_ref_path: PathBuf,
    no_host_file_path: PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, ReceiverStream<ParseOutput>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(input_stream, None)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    let mut cleanup_tasks = vec![host_query_write_task];

    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(host_ref_path, host_query_pipe_path)))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut minimap2_child, minimap2_err_task) = spawn_cmd(
        config.clone(),
        MINIMAP2_TAG,
        minimap2_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = parse_child_output(
        &mut minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-f".to_string(), Some("4".to_string())), // unmapped, since what does not map to host is what we want to filter
            ("--no-PG".to_string(), None),
            ("-h".to_string(), None),                 // Include header
        ]),
    };
    let samtools_args_view = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_child_view, samtools_task_view, samtools_err_task_view) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_view);
    cleanup_tasks.push(samtools_err_task_view);

    let samtools_out_stream_view = parse_child_output(
        &mut samtools_child_view,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let samtools_config_fastq = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let samtools_args_fastq = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_fastq))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_child_fastq, samtools_task_fastq, samtools_err_task_fastq) = stream_to_cmd(
        config.clone(),
        samtools_out_stream_view,
        SAMTOOLS_TAG,
        samtools_args_fastq,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_fastq);
    cleanup_tasks.push(samtools_err_task_fastq);

    let samtools_out_stream_fastq = parse_child_output(
        &mut samtools_child_fastq,
        ChildStream::Stdout,
        ParseMode::Fastq,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    let samtools_out_stream_fastq = ReceiverStream::new(samtools_out_stream_fastq);

    let (host_streams, host_done_rx) = t_junction(
        samtools_out_stream_fastq,
        4,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::IlluminaFastq,
        "align_to_host".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    let mut cleanup_receivers = vec![host_done_rx];
    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_count_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_check_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Check task for host-removed empty FASTQ
    let (check_tx, mut check_rx) = mpsc::channel(1);
    let check_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(no_host_check_stream);
        let mut count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Fastq(_) => {
                    count = 1; // Count only the first record
                    // Continue consuming to keep receiver alive
                }
                _ => return Err(anyhow!("Unexpected item type in no_host_check_stream")),
            }
        }
        check_tx.send(count).await
            .map_err(|e| anyhow!("Failed to send check count: {}", e))?;
        Ok(())
    });
    cleanup_tasks.push(check_task);

    let check_count = check_rx.recv().await.ok_or(PipelineError::EmptyStream)?;
    if check_count == 0 {
        let drain_rxs = vec![no_host_output_stream, no_host_file_stream, no_host_count_stream];
        let drain_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = drain_rxs.into_iter().map(|rx| {
            tokio::spawn(async move {
                let mut stream = ReceiverStream::new(rx);
                while stream.next().await.is_some() {}
                Ok(())
            })
        }).collect();
        cleanup_tasks.extend(drain_tasks);

        return Err(PipelineError::EmptyStream);
    }

    let pigz_args = generate_cli(PIGZ_TAG, &config, None)
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut pigz_child, pigz_stream_task, pigz_err_task) = stream_to_cmd(
        config.clone(),
        no_host_file_stream,
        PIGZ_TAG,
        pigz_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(pigz_stream_task);
    cleanup_tasks.push(pigz_err_task);

    let pigz_out_stream = parse_child_output(
        &mut pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let pigz_write_task = tokio::spawn(stream_to_file(pigz_out_stream, no_host_file_path));
    cleanup_tasks.push(pigz_write_task);

    Ok((ReceiverStream::new(no_host_output_stream), ReceiverStream::new(no_host_count_stream), cleanup_tasks, cleanup_receivers))
}


async fn process_ercc(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    ercc_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base: &str,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let ercc_stats_file_path = file_path_manipulator(
        &PathBuf::from(no_ext_sample_base),
        Some(out_dir),
        None,
        Some("ercc_stats.txt"),
        "_",
    );

    let (ercc_streams, ercc_done_rx) = t_junction(
        input_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::IlluminaFastq,
        "process_ercc_bypass".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    if ercc_streams.len() < 3 {
        return Err(PipelineError::EmptyStream);
    }

    let mut streams_iter = ercc_streams.into_iter();
    let bypass_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let alignment_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let stats_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    cleanup_receivers.push(ercc_done_rx);

    // Process alignment stream for ERCC mapping
    let (ercc_query_write_task, ercc_query_pipe_path) = write_parse_output_to_temp(
        ReceiverStream::new(alignment_stream),
        None,
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    cleanup_tasks.push(ercc_query_write_task);

    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(ercc_path, ercc_query_pipe_path)))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut minimap2_child, minimap2_err_task) = spawn_cmd(
        config.clone(),
        MINIMAP2_TAG,
        minimap2_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = parse_child_output(
        &mut minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    if !minimap2_child.wait().await.map_err(|e| PipelineError::Other(e.into()))?.success() {
        return Err(PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: "Non-zero exit".to_string(),
        });
    }

    // Split SAM stream for check and further processing
    let (sam_streams, sam_done_rx) = t_junction(
        ReceiverStream::new(minimap2_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::JustBytes,
        "process_ercc_empty_check".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    if sam_streams.len() < 2 {
        return Err(PipelineError::EmptyStream);
    }

    cleanup_receivers.push(sam_done_rx);
    let mut sam_streams_iter = sam_streams.into_iter();
    let sam_check_stream = sam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_view_stream = sam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Check task for empty ERCC alignments
    let (count_tx, mut count_rx) = mpsc::channel(1);
    let count_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(sam_check_stream);
        let mut count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(_) => {
                    count = 1; // Count only the first SAM alignment
                    // Continue consuming to keep receiver alive
                }
                _ => return Err(anyhow!("Unexpected item type in sam_check_stream")),
            }
        }
        count_tx
            .send(count)
            .await
            .map_err(|e| anyhow!("Failed to send SAM alignment count: {}", e))?;
        Ok(())
    });
    cleanup_tasks.push(count_task);

    let alignment_count = count_rx
        .recv()
        .await
        .ok_or(PipelineError::EmptyStream)?;

    if alignment_count == 0 {  // Case of no ERCC spike ins
        // Create a dummy stream for the return value
        let (dummy_tx, dummy_rx) = mpsc::channel(config.base_buffer_size);
        drop(dummy_tx); // Ensure stream is empty but valid
        // Drain the bypass and stats streams
        let drain_task = tokio::spawn(async move {
            let mut stream = ReceiverStream::new(bypass_output_stream);
            while stream.next().await.is_some() {}
            Ok(())
        });
        let stats_drain_task = tokio::spawn(async move {
            let mut stream = ReceiverStream::new(stats_stream);
            while stream.next().await.is_some() {}
            Ok(())
        });
        cleanup_tasks.push(drain_task);
        cleanup_tasks.push(stats_drain_task);
        let zero_stats = HashMap::from([
            ("ercc_mapped_reads".to_string(), 0),
            ("ercc_mapped_paired".to_string(), 0),
        ]);
        let stats_write_task = tokio::spawn(async move {
            let mut file = File::create(&ercc_stats_file_path)?;
            writeln!(file, "ercc_mapped_reads: 0\nercc_mapped_paired: 0")?;
            Ok(())
        });
        cleanup_tasks.push(stats_write_task);
        return Ok((
            ReceiverStream::new(dummy_rx),
            Some(tokio::spawn(async { Ok(zero_stats) })),
            cleanup_tasks,
            cleanup_receivers,
        ));
    }

    // Non-zero case: process stats
    let samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([("--no-PG".to_string(), None)]),
    };
    let samtools_args_view = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_child_view, samtools_task_view, samtools_err_task_view) = stream_to_cmd(
        config.clone(),
        sam_view_stream,
        SAMTOOLS_TAG,
        samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_view);
    cleanup_tasks.push(samtools_err_task_view);

    let samtools_out_stream_view = parse_child_output(
        &mut samtools_child_view,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    if !samtools_child_view
        .wait()
        .await
        .map_err(|e| PipelineError::Other(e.into()))?
        .success()
    {
        return Err(PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: "Non-zero exit".to_string(),
        });
    }

    let stats_samtools_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Stats,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let stats_samtools_args = generate_cli(SAMTOOLS_TAG, &config, Some(&stats_samtools_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut stats_samtools_child, stats_samtools_task, stats_samtools_err_task) = stream_to_cmd(
        config.clone(),
        samtools_out_stream_view,
        SAMTOOLS_TAG,
        stats_samtools_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(stats_samtools_task);
    cleanup_tasks.push(stats_samtools_err_task);

    let stats_samtools_out_stream = parse_child_output(
        &mut stats_samtools_child,
        ChildStream::Stdout,
        ParseMode::Lines,
        config.base_buffer_size / 2,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    if !stats_samtools_child
        .wait()
        .await
        .map_err(|e| PipelineError::Other(e.into()))?
        .success()
    {
        return Err(PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: "Non-zero exit on ERCC stats".to_string(),
        });
    }

    let (stats_streams, stats_done_rx) = t_junction(
        ReceiverStream::new(stats_samtools_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::JustBytes,
        "process_ercc_stats".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    if stats_streams.len() < 2 {
        return Err(PipelineError::EmptyStream);
    }

    cleanup_receivers.push(stats_done_rx);
    let mut stats_streams_iter = stats_streams.into_iter();
    let stats_file_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let stats_parse_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let stats_write_task = tokio::spawn(stream_to_file(stats_file_stream, ercc_stats_file_path));
    cleanup_tasks.push(stats_write_task);

    let ercc_stats_task = Some(tokio::spawn(parse_ercc_stats(stats_parse_stream)));

    Ok((
        ReceiverStream::new(bypass_output_stream),
        ercc_stats_task,
        cleanup_tasks,
        cleanup_receivers,
    ))
}


/// Filters input FASTQ stream through Kraken2
///
/// # Arguments
///
/// * `config` - reading stream
/// * `input_stream` - ReceiverStream<ParseOutput>,
/// * 'targer_ref_path' -
///
/// # Returns
/// samtools_sort_out_stream: Result<(ReceiverStream<ParseOutput>, <--- FASTQ uncompressed
/// cleanup_tasks
/// quast_write_tasks
///
async fn filter_with_kraken(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    target_ref_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
    target_taxid: &str,
) -> Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];

    let kraken2_report_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("kraken2_report.txt"),
        "_",
    );
    let final_compressed_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("classified_filtered.fq.gz"),
        "_",
    );
    let kraken2_pipe_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("kraken2_classified_pipe.fq"),
        "_",
    );

    // Remove existing pipe if exists to avoid errors
    let _ = tokio::fs::remove_file(&kraken2_pipe_path).await;

    // Create named pipe for Kraken2 classified output
    Command::new("mkfifo")
        .arg(&kraken2_pipe_path)
        .status()
        .map_err(|e| PipelineError::Other(anyhow::anyhow!("Failed to create named pipe: {}", e)))?;

    // Write input stream to temp file for Kraken2
    let (kraken_query_write_task, kraken_query_pipe_path) = write_parse_output_to_temp(input_stream, None)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    cleanup_tasks.push(kraken_query_write_task);

    // Run Kraken2 with named pipe for classified output
    let kraken2_config = Kraken2Config {
        report_path: kraken2_report_path.clone(),
        classified_path: kraken2_pipe_path.clone(),
        fastq_path: kraken_query_pipe_path.clone(),
    };
    let kraken2_args = generate_cli(KRAKEN2_TAG, &config, Some(&kraken2_config))?;

    let (mut _kraken2_child, kraken2_err_task) = spawn_cmd(
        config.clone(),
        KRAKEN2_TAG,
        kraken2_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KRAKEN2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(kraken2_err_task);

    // Stream classified FASTQ from named pipe using parse_fastq
    let kraken2_file = tokio::fs::File::open(&kraken2_pipe_path)
        .await
        .map_err(|e| PipelineError::Other(anyhow::anyhow!("Failed to open named pipe: {}", e)))?;
    let parse_rx = parse_fastq(kraken2_file, config.base_buffer_size)
        .await
        .map_err(|e| PipelineError::Other(anyhow::anyhow!("Failed to parse FASTQ from named pipe: {}", e)))?;

    // Filter using parse_and_filter_fastq_id
    let taxid = target_taxid;
    let pattern = format!("kraken:taxid|{}", taxid);
    let filter_fn = move |id: &str| id.contains(&pattern);
    let (filtered_rx, filter_task) = parse_and_filter_fastq_id(parse_rx, config.base_buffer_size, filter_fn.clone());
    cleanup_tasks.push(filter_task);

    // Convert SequenceRecord back to ParseOutput::Fastq for downstream compatibility
    let (parse_output_tx, parse_output_rx) = mpsc::channel(config.base_buffer_size);
    let conversion_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(filtered_rx);
        let mut count = 0;
        while let Some(record) = stream.next().await {
            validate_sequence(&Arc::new(record.seq().to_vec()), b"ACGTN")
                .map_err(|e| anyhow::anyhow!("Sequence validation failed at record {}: {}", count + 1, e))?;
            parse_output_tx.send(ParseOutput::Fastq(record)).await
                .map_err(|e| anyhow::anyhow!("Failed to send ParseOutput at record {}: {}", count + 1, e))?;
            count += 1;
        }
        Ok(())
    });
    cleanup_tasks.push(conversion_task);


    let filtered_stream = ReceiverStream::new(parse_output_rx);

    // split stream for output and compression
    let (kraken_streams, kraken_done_rx) = t_junction(
        filtered_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::IlluminaFastq,
        "filter_reads_output".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    if kraken_streams.len() < 2 {
        return Err(PipelineError::EmptyStream);
    }

    let mut cleanup_receivers = vec![kraken_done_rx];
    let mut streams_iter = kraken_streams.into_iter();
    let kraken_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let kraken_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Run pigz to compress filtered output
    let pigz_args = generate_cli(PIGZ_TAG, &config, None)?;
    let (mut pigz_child, pigz_stream_task, pigz_err_task) = stream_to_cmd(
        config.clone(),
        kraken_file_stream,
        PIGZ_TAG,
        pigz_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(pigz_stream_task);
    cleanup_tasks.push(pigz_err_task);

    let pigz_out_stream = parse_child_output(
        &mut pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let pigz_write_task = tokio::spawn(stream_to_file(pigz_out_stream, final_compressed_path));
    cleanup_tasks.push(pigz_write_task);

    let cleanup_pipe_task = tokio::spawn(async move {
        tokio::fs::remove_file(&kraken2_pipe_path)
            .await
            .map_err(|e| anyhow::anyhow!("Failed to remove named pipe: {}", e))?;
        Ok(())
    });
    cleanup_tasks.push(cleanup_pipe_task);

    Ok((ReceiverStream::new(kraken_output_stream), cleanup_tasks, cleanup_receivers))
}


/// Runs minimap2 to align to the target reference, then sorts.
///
/// # Arguments
///
/// * `config` - reading stream
/// * `input_stream` - stream buffer size
/// * 'targer_ref_path' -
///
/// # Returns
/// samtools_sort_out_stream: Result<(ReceiverStream<ParseOutput>, <--- SAM uncompressed
/// cleanup_tasks
/// quast_write_tasks
///
async fn align_to_target(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    target_ref_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf
) -> Result<(ReceiverStream<ParseOutput>,  Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, Vec<JoinHandle<Result<(), anyhow::Error>>>, PathBuf), PipelineError> {

    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];
    let mut quast_write_tasks = vec![];

    let (align_query_write_task, align_query_pipe_path) = write_parse_output_to_temp(input_stream, None)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    quast_write_tasks.push(align_query_write_task);

    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(target_ref_path, align_query_pipe_path)))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut minimap2_child, minimap2_err_task) = spawn_cmd(
        config.clone(),
        MINIMAP2_TAG,
        minimap2_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = parse_child_output(
        &mut minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-O".to_string(), Some("sam".to_string())),
            ("-".to_string(), None)
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = parse_child_output(
        &mut samtools_sort_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let samtools_sort_out_stream = ReceiverStream::new(samtools_sort_out_stream);

    let align_bam_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("target_aligned.bam"),
        "_",
    );

    let (sam_streams, sam_done_rx) = t_junction(
        samtools_sort_out_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::JustBytes,
        "align_to_target".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(sam_done_rx);
    let mut streams_iter = sam_streams.into_iter();
    let sam_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let align_samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-h".to_string(), None),
            ("-b".to_string(), None),
            ("-o".to_string(), Some(align_bam_path.display().to_string())),
            ("-".to_string(), None),
        ])
    };
    let align_samtools_args_view = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&align_samtools_config_view),
    )?;
    let (mut _align_samtools_child_view, align_samtools_view_stream_task, align_samtools_view_err_task) = stream_to_cmd(
        config.clone(),
        sam_file_stream,
        SAMTOOLS_TAG,
        align_samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool:  SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    quast_write_tasks.push(align_samtools_view_stream_task);
    cleanup_tasks.push(align_samtools_view_err_task);


    Ok((
        ReceiverStream::new(sam_output_stream),
        cleanup_tasks,
        cleanup_receivers,
        quast_write_tasks,
        align_bam_path
    ))
}


async fn generate_consensus(
    config: Arc<RunConfig>,
    bam_stream: ReceiverStream<ParseOutput>,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, ReceiverStream<ParseOutput>, PathBuf, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];
    let mut quast_write_tasks = vec![];

    let samtools_consensus_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Consensus,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let samtools_consensus_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_consensus_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_consensus_child, samtools_consensus_task, samtools_consensus_err_task) = stream_to_cmd(
        config.clone(),
        bam_stream.into_inner(),
        SAMTOOLS_TAG,
        samtools_consensus_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_consensus_task);
    cleanup_tasks.push(samtools_consensus_err_task);

    let samtools_consensus_out_stream = parse_child_output(  // FASTQ output
        &mut samtools_consensus_child,
        ChildStream::Stdout,
        ParseMode::Fasta,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    let samtools_consensus_out_stream = ReceiverStream::new(samtools_consensus_out_stream);

    let (consensus_streams, consensus_done_rx) = t_junction(
        samtools_consensus_out_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::JustBytes,
        "generate_consensus".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(consensus_done_rx);
    let mut streams_iter = consensus_streams.into_iter();
    let consensus_realign_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let consensus_stats_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let consensus_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let consensus_file_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("consensus.fa"),
        "_",
    );

    let consensus_write_task = tokio::spawn(stream_to_file(
        consensus_file_stream,
        consensus_file_path.clone(),
    ));
    quast_write_tasks.push(consensus_write_task);

    Ok((
        ReceiverStream::new(consensus_realign_stream),
        ReceiverStream::new(consensus_stats_stream),
        consensus_file_path,
        cleanup_tasks,
        cleanup_receivers, quast_write_tasks,
    ))
}


async fn call_variants(
    config: Arc<RunConfig>,
    bam_stream: ReceiverStream<ParseOutput>,
    target_ref_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, PathBuf, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];

    let bcftools_mpileup_config = BcftoolsConfig {
        subcommand: BcftoolsSubcommand::Mpileup,
        subcommand_fields: HashMap::from([("-f".to_string(), Some(target_ref_path.to_string_lossy().into_owned())), ("-".to_string(), None),]),
    };
    let bcftools_mpileup_args = generate_cli(BCFTOOLS_TAG, &config, Some(&bcftools_mpileup_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bcftools_mpileup_child, bcftools_mpileup_task, bcftools_mpileup_err_task) = stream_to_cmd(
        config.clone(),
        bam_stream.into_inner(),
        BCFTOOLS_TAG,
        bcftools_mpileup_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bcftools_mpileup_task);
    cleanup_tasks.push(bcftools_mpileup_err_task);

    let bcftools_mpileup_out_stream = parse_child_output(
        &mut bcftools_mpileup_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let bcftools_call_config = BcftoolsConfig {
        subcommand: BcftoolsSubcommand::Call,
        subcommand_fields: HashMap::from([
            ("--ploidy".to_string(), Some("1".to_string())),
            ("-m".to_string(), None),
            ("-v".to_string(), None),
            ("-P".to_string(), Some(config.args.bcftools_call_theta.to_string())),
            ("-Ov".to_string(), None),
            ("-".to_string(), None),
        ])
    };
    let bcftools_call_args = generate_cli(BCFTOOLS_TAG, &config, Some(&bcftools_call_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bcftools_call_child, bcftools_call_task, bcftools_call_err_task) = stream_to_cmd(
        config.clone(),
        bcftools_mpileup_out_stream,
        BCFTOOLS_TAG,
        bcftools_call_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bcftools_call_task);
    cleanup_tasks.push(bcftools_call_err_task);

    let bcftools_call_out_stream = parse_child_output(
        &mut bcftools_call_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    let bcftools_call_out_stream = ReceiverStream::new(bcftools_call_out_stream);

    let (vcf_streams, vcf_done_rx) = t_junction(
        bcftools_call_out_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::JustBytes,
        "call_variants".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(vcf_done_rx);

    let mut streams_iter = vcf_streams.into_iter();
    let vcf_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let vcf_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let vcf_file_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("variants.vcf"),
        "_",
    );

    let vcf_write_task = tokio::spawn(stream_to_file(
        vcf_file_stream,
        vcf_file_path.clone(),
    ));
    cleanup_tasks.push(vcf_write_task);

    Ok((ReceiverStream::new(vcf_output_stream), vcf_file_path, cleanup_tasks, cleanup_receivers))
}


async fn realign_consensus_to_ref(
    config: Arc<RunConfig>,
    consensus_realign_stream: Receiver<ParseOutput>,
    target_ref_fasta_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let reference_file = TokioFile::open(&target_ref_fasta_path).await
        .map_err(|e| PipelineError::Other(e.into()))?;
    let reference_rx = parse_bytes(reference_file, config.base_buffer_size).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let realign_streams = vec![consensus_realign_stream, reference_rx];
    let (combined_rx, combined_task) = y_junction(realign_streams, config.base_buffer_size).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let realign_consensus_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("consensus_realigned.fa"),
        "_"
    );
    let realign_mafft_args = generate_cli(MAFFT_TAG, &config, None)
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (
        mut realign_consensus_mafft_child,
        realign_consensus_mafft_task,
        realign_consensus_mafft_err_task
    ) = stream_to_cmd(config.clone(),
                      combined_rx,
                      MAFFT_TAG,
                      realign_mafft_args,
                      StreamDataType::JustBytes,
                      config.args.verbose
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;

    let mut cleanup_tasks = vec![combined_task, realign_consensus_mafft_task, realign_consensus_mafft_err_task];

    let consensus_samtools_out_stream = parse_child_output(
        &mut realign_consensus_mafft_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;

    // keep alive task
    let mafft_wait_task = tokio::spawn(async move {
        let status = realign_consensus_mafft_child.wait().await?;
        if !status.success() {
            return Err(anyhow!("MAFFT exited with non-zero status: {}", status));
        }
        Ok(())
    });
    cleanup_tasks.push(mafft_wait_task);

    let realign_consensus_write_task = tokio::spawn(stream_to_file(
        consensus_samtools_out_stream,
        realign_consensus_path,
    ));
    cleanup_tasks.push(realign_consensus_write_task);



    Ok((cleanup_tasks))
}


async fn calculate_statistics(
    config: Arc<RunConfig>,
    no_ext_sample_base: &str,
    consensus_bam_stats_stream: Option<ReceiverStream<ParseOutput>>,
    consensus_bam_depth_stream: Option<ReceiverStream<ParseOutput>>,
    no_host_seqkit_out_stream_stats: Receiver<ParseOutput>,
    ercc_stats_task: Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>>,
    consensus_stats_stream: Option<ReceiverStream<ParseOutput>>,
    call_bcftools_stats_stream: Option<ReceiverStream<ParseOutput>>,
    out_dir: &PathBuf,
    technology: Technology,
) -> Result<(), anyhow::Error> {
    let mut local_cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();


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
            let (mut stats_samtools_child_stats, stats_samtools_task_stats, stats_samtools_err_task_stats) = stream_to_cmd(config.clone(),
                                                                                                                           stream.into_inner(),
                                                                                                                           SAMTOOLS_TAG,
                                                                                                                           stats_samtools_args_stats.clone(),
                                                                                                                           StreamDataType::JustBytes,
                                                                                                                           config.args.verbose,
            ).await?;
            local_cleanup_tasks.push(stats_samtools_task_stats);
            local_cleanup_tasks.push(stats_samtools_err_task_stats);
            parse_child_output(&mut stats_samtools_child_stats, ChildStream::Stdout, ParseMode::Lines, config.base_buffer_size / 2).await?
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
            let (mut depth_samtools_child, depth_samtools_task, depth_samtools_err_task) = stream_to_cmd(config.clone(),
                                                                                                         stream.into_inner(),
                                                                                                         SAMTOOLS_TAG,
                                                                                                         depth_samtools_args,
                                                                                                         StreamDataType::JustBytes,
                                                                                                         config.args.verbose,
            ).await?;
            local_cleanup_tasks.push(depth_samtools_task);
            local_cleanup_tasks.push(depth_samtools_err_task);
            parse_child_output(
                &mut depth_samtools_child,
                ChildStream::Stdout,
                ParseMode::Lines,
                config.base_buffer_size / 2,
            ).await?
        }
        None => {
            return Err(anyhow!("consensus_bam_depth_stream is not available"));
        }
    };

    let depth_map = parse_samtools_depth(depth_samtools_out_stream).await?;
    if depth_map.is_empty() {
        return Err(anyhow!("No depth data found"));
    }
    let first_chr = depth_map.keys().next().ok_or_else(|| anyhow!("No chromosomes found"))?.clone();
    let first_chr_depth_map = depth_map.get(&first_chr).unwrap();
    let depths: Vec<u32> = first_chr_depth_map.values().copied().collect();
    let samtools_depth_stats = compute_depth_stats(&depths)?;

    let depth_plot_path = file_path_manipulator(
        &PathBuf::from(no_ext_sample_base),
        Some(out_dir),
        None,
        Some("depth.png"),
        "_"
    );

    plot_depths(&first_chr_depth_map, no_ext_sample_base, &depth_plot_path)?;

    let seqkit_stats = parse_seqkit_stats(no_host_seqkit_out_stream_stats).await?;

    let ercc_stats = if let Technology::Illumina = technology {
        match ercc_stats_task {
            Some(task) => match task.await? {
                Ok(stats) => stats,
                Err(e) => return Err(anyhow!("Failed to parse ERCC stats: {}", e)),
            },
            None => return Err(anyhow!("ERCC stats task not initialized for Illumina technology")),
        }
    } else {
        HashMap::new()
    };

    let allele_counts = if let Some(stream) = consensus_stats_stream {
        compute_allele_counts(stream.into_inner()).await?
    } else {
        HashMap::new()
    };

    let (ref_snps, ref_mnps, ref_indels) = if let Some(stream) = call_bcftools_stats_stream {
        parse_vcf_stream(stream.into_inner()).await?
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
        sample_name: no_ext_sample_base.to_string(),
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

    let stats_file_path = out_dir.join(format!("{}_stats.json", no_ext_sample_base));
    let mut stats_file = File::create(&stats_file_path)?;
    serde_json::to_writer_pretty(&mut stats_file, &stats)?;

    try_join_all(local_cleanup_tasks).await?.into_iter().collect::<Result<Vec<_>, _>>()?;

    Ok(())
}

async fn evaluate_assembly(
    config: Arc<RunConfig>,
    target_ref_fasta_path: PathBuf,
    align_bam_path: PathBuf,
    consensus_file_path: PathBuf,
) -> Result<(Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];
    let quast_config = QuastConfig {
        ref_fasta: target_ref_fasta_path.to_string_lossy().into_owned(),
        ref_bam: align_bam_path.to_string_lossy().into_owned(),
        assembly_fasta: consensus_file_path.to_string_lossy().into_owned(),
    };

    let assembly_eval_quast_args = generate_cli(
        QUAST_TAG,
        &*config,
        Some(&quast_config),
    ).map_err(|e| PipelineError::ToolExecution {
        tool: QUAST_TAG.to_string(),
        error: e.to_string(),
    })?;

    let (mut assembly_eval_quast_child, assembly_eval_quast_err_task) = spawn_cmd(
        config.clone(),
        QUAST_TAG,
        assembly_eval_quast_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: QUAST_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(assembly_eval_quast_err_task);

    // Await child exit in a task to keep it alive
    let quast_wait_task = tokio::spawn(async move {
        let status = assembly_eval_quast_child.wait().await?;
        if !status.success() {
            return Err(anyhow!("QUAST exited with non-zero status: {}", status));
        }
        Ok(())
    });
    cleanup_tasks.push(quast_wait_task);


    Ok((cleanup_tasks))
}





pub async fn run(config: Arc<RunConfig>) -> Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let ram_temp_dir = config.ram_temp_dir.clone();
    let out_dir = config.out_dir.clone();
    let mut temp_files: Vec<NamedTempFile>  = Vec::new();
    let mut cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();
    let mut quast_write_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut stats_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut ercc_stats_task: Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>> = None;

    let mut target_ref_fasta_path: Option<PathBuf> = None;
    let mut consensus_file_path: Option<PathBuf> = None;
    let mut align_bam_path: Option<PathBuf> = None;
    let mut align_sam_stats_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut align_sam_depth_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut align_sam_eval_stream: Option<Receiver<ParseOutput>> = None;
    let mut consensus_stats_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut call_bcftools_stats_stream: Option<ReceiverStream<ParseOutput>> = None;



    // External tools check
    check_versions(vec![SAMTOOLS_TAG, MINIMAP2_TAG, FASTP_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, SEQKIT_TAG, QUAST_TAG]).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Arguments and files check
    let file1_path: PathBuf = match &config.args.file1 {
        Some(file) => {
            let file1_full_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd), None, None, "");
            if file1_full_path.exists() {
                file1_full_path
            } else {
                return Err(PipelineError::FileNotFound(file1_full_path));
            }
        }
        None => return Err(PipelineError::InvalidConfig("File1 path required".to_string())),
    };

    let sample_base: String;
    let file1_r1r2 = r1r2_base(&file1_path);
    sample_base = match file1_r1r2.file_name {
        Some(prefix) => prefix,
        None => {
            eprintln!("No R1 tag found. Using bare file 1 stem as sample_base.");
            file1_path.to_string_lossy().into_owned()
        }
    };

    if !file1_path.exists() {
        return Err(PipelineError::FileNotFound(file1_path));
    }

    let sample_base_buf: PathBuf = PathBuf::from(&sample_base);
    let (no_ext_sample_base_buf, _) = extension_remover(&sample_base_buf);
    let no_ext_sample_base = no_ext_sample_base_buf.to_string_lossy().into_owned();

    let file2_path: Option<PathBuf> = match &config.args.file2 {
        Some(file) => {
            let file2_full_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd.clone()), None, None, "");
            if file2_full_path.exists() {
                Some(file2_full_path)
            } else {
                eprintln!("File2 path does not exist: {}", file2_full_path.display());
                None
            }
        }
        None => None,
    };

    let technology = config.args.technology.clone();

    let ref_db_path: Option<PathBuf> = config.args.ref_db.as_ref().map(PathBuf::from);

    let ercc_path = file_path_manipulator(&PathBuf::from(&config.args.ercc_sequences), Some(&cwd), None, None, "");
    if !ercc_path.exists() {
        return Err(PipelineError::FileNotFound(ercc_path));
    }

    // Input Validation
    let (val_fastp_out_stream, validate_cleanup_tasks, validate_cleanup_receivers) = validate_input(
        config.clone(),
        file1_path,
        file2_path,
        sample_base_buf.clone(),
        &out_dir,
    ).await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // Retrieve Index: if ref_db is None will return a None quickly
    let index_start = Instant::now();
    let h5_index = get_index(&config.args)
        .await
        .map_err(|e| PipelineError::ReferenceRetrievalFailed(e.to_string()))?;
    if h5_index.is_some() {
        println!("Index retrieve time: {} milliseconds.", index_start.elapsed().as_millis());
    }

    // Fetch host reference
    let (host_ref_fasta_path, host_ref_temp, host_ref_write_task) = fetch_host_reference(
        &*config,
        ref_db_path.clone(),
        &ram_temp_dir,
        h5_index.as_ref()
    ).await?;
    host_ref_write_task
        .await
        .map_err(|e| PipelineError::Other(e.into()))?
        .map_err(|e| PipelineError::Other(e))?;
    temp_files.push(host_ref_temp);


    // Host Removal
    let no_host_file_path = file_path_manipulator(
        &no_ext_sample_base_buf,
        Some(&out_dir),
        None,
        Some("no_host.fq.gz"),
        "_",
    );

    //
    let (no_host_output_stream, no_host_seqkit_out_stream_stats, no_host_cleanup_tasks, no_host_cleanup_receivers) = align_to_host(
        config.clone(),
        val_fastp_out_stream,
        host_ref_fasta_path,
        no_host_file_path,
    ).await?;
    cleanup_tasks.extend(no_host_cleanup_tasks);
    cleanup_receivers.extend(no_host_cleanup_receivers);

    // Counting stats for the host-removed reads
    let stats_seqkit_config_stats = SeqkitConfig {
        subcommand: SeqkitSubcommand::Stats,
        subcommand_fields: HashMap::from([]),
    };
    let stats_seqkit_args_stats = generate_cli(
        SEQKIT_TAG,
        &config,
        Some(&stats_seqkit_config_stats),
    )?;
    let no_host_seqkit_out_stream_stats = no_host_seqkit_out_stream_stats.into_inner();
    let (mut no_host_seqkit_child_stats, no_host_seqkit_task_stats, no_host_seqkit_err_task_stats) = stream_to_cmd(config.clone(), no_host_seqkit_out_stream_stats, SEQKIT_TAG, stats_seqkit_args_stats, StreamDataType::JustBytes, config.args.verbose).await?;
    let no_host_seqkit_out_stream_stats = parse_child_output(
        &mut no_host_seqkit_child_stats,
        ChildStream::Stdout,
        ParseMode::Lines,
        config.base_buffer_size / 2,
    ).await?;
    stats_tasks.push(no_host_seqkit_task_stats);
    cleanup_tasks.push(no_host_seqkit_err_task_stats);


    //*****************
    // Split by Technology
    match technology {
        Technology::Illumina => {
            eprintln!("Technology: Illumina");


            // Fetch target reference
            let (target_ref_fasta_path_inner, target_ref_temp, target_ref_write_task) = fetch_target_reference(
                &*config,
                ref_db_path.clone(),
                &ram_temp_dir,
                h5_index.as_ref()
            ).await?;
            target_ref_fasta_path = Some(target_ref_fasta_path_inner.clone());

            // ERCC
            let (no_host_ercc_stream, ercc_stats_out_task, mut ercc_cleanup_tasks, mut ercc_cleanup_receivers) = process_ercc(
                config.clone(),
                no_host_output_stream,
                ercc_path,
                &out_dir,
                &no_ext_sample_base,
            ).await?;
            ercc_stats_task = ercc_stats_out_task;
            cleanup_tasks.append(&mut ercc_cleanup_tasks);
            cleanup_receivers.append(&mut ercc_cleanup_receivers);


            target_ref_write_task
                .await
                .map_err(|e| PipelineError::Other(e.into()))?
                .map_err(|e| PipelineError::Other(e))?;
            temp_files.push(target_ref_temp);


            // Filter Reads
            let (filter_reads_out_stream, filter_reads_cleanup_tasks, filter_reads_cleanup_receivers) = filter_with_kraken(
                config.clone(),
                no_host_ercc_stream,
                target_ref_fasta_path_inner.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
                config.args.ref_taxid.as_ref().expect("ref_taxid must be set"),
            ).await?;
            cleanup_tasks.extend(filter_reads_cleanup_tasks);
            cleanup_receivers.extend(filter_reads_cleanup_receivers);


            // Align Reads to Target
            let (sam_output_stream, align_cleanup_tasks, align_cleanup_receivers, align_quast_tasks, align_bam_path_inner) = align_to_target(
                config.clone(),
                filter_reads_out_stream,
                target_ref_fasta_path_inner.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
            ).await?;
            cleanup_tasks.extend(align_cleanup_tasks);
            cleanup_receivers.extend(align_cleanup_receivers);
            quast_write_tasks.extend(align_quast_tasks);

            align_bam_path = Some(align_bam_path_inner);

            // Split SAM streams for bypass, stats, etc
            let (align_sam_streams, align_sam_done_rx) = t_junction(
                sam_output_stream,
                5,
                config.base_buffer_size,
                config.args.stall_threshold,
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::JustBytes,
                "pipeline_aligned_sam_split".to_string(),
                None
            )
                .await
                .map_err(|_| PipelineError::StreamDataDropped)?;
            cleanup_receivers.push(align_sam_done_rx);
            let mut align_streams_iter = align_sam_streams.into_iter();
            let align_sam_output_stream = align_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
            let align_sam_call_stream = align_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
            align_sam_stats_stream = Some(ReceiverStream::new(align_streams_iter.next().ok_or(PipelineError::EmptyStream)?));
            align_sam_depth_stream = Some(ReceiverStream::new(align_streams_iter.next().ok_or(PipelineError::EmptyStream)?));
            align_sam_eval_stream = Some(align_streams_iter.next().ok_or(PipelineError::EmptyStream).unwrap());

            let mut align_sam_output_stream = ReceiverStream::new(align_sam_output_stream);
            let mut align_sam_call_stream = ReceiverStream::new(align_sam_call_stream);


            // Make Consensus
            let (consensus_realign_stream, consensus_stats_stream_rx, consensus_file_path_x, consensus_cleanup_tasks, consensus_cleanup_receivers, consensus_quast_tasks) = generate_consensus(
                config.clone(),
                align_sam_output_stream,
                &out_dir,
                &no_ext_sample_base_buf,
            ).await?;

            consensus_file_path = Some(consensus_file_path_x);
            cleanup_tasks.extend(consensus_cleanup_tasks);
            cleanup_receivers.extend(consensus_cleanup_receivers);
            quast_write_tasks.extend(consensus_quast_tasks);
            consensus_stats_stream = Some(consensus_stats_stream_rx);

            // Call Variants
            let (call_bcftools_stats_stream_out, _, call_cleanup_tasks, call_cleanup_receivers) = call_variants(
                config.clone(),
                align_sam_call_stream,
                target_ref_fasta_path_inner.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
            ).await?;
            cleanup_tasks.extend(call_cleanup_tasks);
            cleanup_receivers.extend(call_cleanup_receivers);
            call_bcftools_stats_stream = Some(call_bcftools_stats_stream_out);

            // Realign Consensus to Ref
            let realign_cleanup_tasks = realign_consensus_to_ref(
                config.clone(),
                consensus_realign_stream.into_inner(),
                target_ref_fasta_path_inner.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
            ).await?;
            cleanup_tasks.extend(realign_cleanup_tasks);


        }

        Technology::ONT => {
            eprintln!("Technology: ONT not ready");
        }
    }


    let results = try_join_all(stats_tasks).await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }



    // // Calculate Statistics
    calculate_statistics(
        config.clone(),
        &no_ext_sample_base,
        align_sam_stats_stream,
        align_sam_depth_stream,
        no_host_seqkit_out_stream_stats,
        ercc_stats_task,
        consensus_stats_stream,
        call_bcftools_stats_stream,
        &out_dir,
        technology,
    ).await?;

    // Assembly Evaluation
    let results = try_join_all(quast_write_tasks).await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }

    if let (Some(target_ref_fasta_path), Some(align_bam_path), Some(consensus_file_path)) = (target_ref_fasta_path, align_bam_path, consensus_file_path) {
        evaluate_assembly(
            config,
            target_ref_fasta_path,
            align_bam_path,
            consensus_file_path,
        ).await?;
    }

    // Cleanup
    let results = try_join_all(cleanup_tasks).await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }
    for receiver in cleanup_receivers {
        receiver.await
            .map_err(|e| PipelineError::Other(e.into()))?
            .map_err(|e| PipelineError::Other(e))?;
    }
    drop(temp_files);

    println!("Finished generating consensus genome");
    Ok(())
}



pub async fn old_run(config: &RunConfig) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", config.args.module);

    if config.args.ref_taxid == None {
        return Err(anyhow!("The ref_taxid argument must be set for this pipeline. Consult https://www.ncbi.nlm.nih.gov/taxonomy for the taxid you need."));
    }

    let cwd = config.cwd.clone();
    let ram_temp_dir = config.ram_temp_dir.clone();
    let out_dir = config.out_dir.clone();
    let mut temp_files = Vec::new();
    let config_arc = Arc::new((*config).clone());

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
            let file1_full_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd), None, None, "");
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
            let file2_full_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd.clone()), None, None, "");
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

    let ref_db_path: Option<PathBuf> = config.args.ref_db.as_ref().map(|ref_db| PathBuf::from(ref_db));

    let ercc_path = file_path_manipulator(&PathBuf::from(&config.args.ercc_sequences), Some(&cwd), None, None, "");
    if !ercc_path.exists() {
        return Err(anyhow!("Specified ercc {:?} does not exist.", ercc_path));
    }

    //*****************
    // Input Validation

    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(&sample_base), Some(&out_dir), None, Some("validated"), "_");
    let rx = read_and_interleave_sequences(file1_path, file2_path, Some(technology.clone()), config.args.max_reads, config.args.min_read_len, config.args.max_read_len)?;
    let val_rx_stream = ReceiverStream::new(rx);
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold.try_into().unwrap(),
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None
    )
        .await?;
    cleanup_receivers.push(val_done_rx);

    let mut streams_iter = val_streams.into_iter();
    let val_fastp_stream = streams_iter.next().unwrap();
    let val_pigz_stream = streams_iter.next().unwrap();

    // Pigz stream to intermediate file output
    let val_pigz_args = generate_cli(PIGZ_TAG, &config, None)?;
    let (mut val_pigz_child, val_pigz_stream_task, val_pigz_err_task) = stream_to_cmd(config_arc.clone(), val_pigz_stream, PIGZ_TAG, val_pigz_args, StreamDataType::IlluminaFastq, config.args.verbose).await?;
    cleanup_tasks.push(val_pigz_stream_task);
    cleanup_tasks.push(val_pigz_err_task);

    let val_pigz_out_stream = parse_child_output(
        &mut val_pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await?;
    let val_pigz_write_task = tokio::spawn(stream_to_file(
        val_pigz_out_stream,
        validated_interleaved_file_path,
    ));
    cleanup_tasks.push(val_pigz_write_task);

    // Fastp stream
    let val_fastp_args = generate_cli(FASTP_TAG, &config, None)?;
    let (mut val_fastp_child, val_fastp_stream_task, val_fastp_err_task) = stream_to_cmd(config_arc.clone(), val_fastp_stream, FASTP_TAG, val_fastp_args, StreamDataType::IlluminaFastq, config.args.verbose).await?;
    cleanup_tasks.push(val_fastp_stream_task);
    cleanup_tasks.push(val_fastp_err_task);
    let val_fastp_out_stream = parse_child_output(
        &mut val_fastp_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await?;
    let val_fastp_out_stream = ReceiverStream::new(val_fastp_out_stream);

    //*****************
    // Fetch reference
    let index_start = Instant::now();
    let h5_index = get_index(&config.args).await?;
    println!("Index retrieve time: {} milliseconds.", index_start.elapsed().as_millis());

    //*****************
    // Host Removal
    let (_host_accession, host_seq) = retrieve_h5_seq(
        config.args.host_accession.clone(),
        config.args.host_sequence.clone(),
        ref_db_path.as_ref(),
        h5_index.as_ref(),
    ).await?;
    let host_ref_temp = NamedTempFile::new_in(&ram_temp_dir)?;
    let host_ref_fasta_path = host_ref_temp.path().to_path_buf();
    temp_files.push(host_ref_temp);
    let host_ref_write_task = write_vecu8_to_file(host_seq.clone(), &host_ref_fasta_path, config.base_buffer_size).await?;
    host_ref_write_task.await?; // This has to be done immediately to make sure the host query minimap2 can read

    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(val_fastp_out_stream, None).await?;
    cleanup_tasks.push(host_query_write_task);
    let host_minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(host_ref_fasta_path.clone(), host_query_pipe_path.clone())))?;
    let (mut host_minimap2_child, host_minimap2_err_task) = spawn_cmd(config_arc.clone(), MINIMAP2_TAG, host_minimap2_args, config.args.verbose).await?;
    cleanup_tasks.push(host_minimap2_err_task);

    let host_minimap2_out_stream = parse_child_output(
        &mut host_minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
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

    let (mut host_samtools_child_view, host_samtools_task_view, host_samtools_err_task_view) = stream_to_cmd(config_arc.clone(), host_minimap2_out_stream, SAMTOOLS_TAG, host_samtools_args_view, StreamDataType::JustBytes, config.args.verbose).await?;
    let host_samtools_out_stream_view = parse_child_output(
        &mut host_samtools_child_view,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
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

    let (mut host_samtools_child_fastq, host_samtools_task_fastq, host_samtools_err_task_fastq) = stream_to_cmd(config_arc.clone(), host_samtools_out_stream_view, SAMTOOLS_TAG, host_samtools_args_fastq, StreamDataType::JustBytes, config.args.verbose).await?;
    let host_samtools_out_stream_fastq = parse_child_output(
        &mut host_samtools_child_fastq,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await?;
    let host_samtools_out_stream_fastq = ReceiverStream::new(host_samtools_out_stream_fastq);
    cleanup_tasks.push(host_samtools_task_fastq);
    cleanup_tasks.push(host_samtools_err_task_fastq);

    let no_host_file_path = file_path_manipulator(&no_ext_sample_base_buf, Some(&out_dir), None, Some("no_host.fq.gz"), "_");

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
        config.base_buffer_size,
        config.args.stall_threshold.try_into().unwrap(),
        Some(config.args.stream_sleep_ms),
        50,
        StreamDataType::IlluminaFastq,
        "align_to_host".to_string(),
        None
    )
        .await?;
    cleanup_receivers.push(host_done_rx);

    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().unwrap();
    let no_host_file_stream = streams_iter.next().unwrap();
    let no_host_count_stream = streams_iter.next().unwrap();

    let (mut no_host_seqkit_child_stats, no_host_seqkit_task_stats, no_host_seqkit_err_task_stats) = stream_to_cmd(config_arc.clone(), no_host_count_stream, SEQKIT_TAG, stats_seqkit_args_stats, StreamDataType::JustBytes, config.args.verbose).await?;
    let no_host_seqkit_out_stream_stats = parse_child_output(
        &mut no_host_seqkit_child_stats,
        ChildStream::Stdout,
        ParseMode::Lines,
        config.base_buffer_size / 2,
    ).await?;
    stats_tasks.push(no_host_seqkit_task_stats);
    cleanup_tasks.push(no_host_seqkit_err_task_stats);

    let host_pigz_args = generate_cli(PIGZ_TAG, &config, None)?;
    let (mut host_pigz_child, host_pigz_stream_task, host_pigz_err_task) = stream_to_cmd(config_arc.clone(), no_host_file_stream, PIGZ_TAG, host_pigz_args, StreamDataType::IlluminaFastq, config.args.verbose).await?;
    cleanup_tasks.push(host_pigz_stream_task);
    cleanup_tasks.push(host_pigz_err_task);

    let host_pigz_out_stream = parse_child_output(
        &mut host_pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await?;
    let host_pigz_write_task = tokio::spawn(stream_to_file(
        host_pigz_out_stream,
        no_host_file_path,
    ));
    cleanup_tasks.push(host_pigz_write_task);

    let no_host_output_stream = ReceiverStream::new(no_host_output_stream);

    // Declare file paths as Option<PathBuf> before the match
    let mut target_ref_fasta_path: Option<PathBuf> = None;
    let mut align_bam_path: Option<PathBuf> = None;
    let mut align_query_pipe_path: Option<PathBuf> = None;
    let mut consensus_file_path: Option<PathBuf> = None;
    let mut consensus_eval_stream: Option<Receiver<ParseOutput>> = None;
    let mut consensus_bam_eval_stream: Option<Receiver<ParseOutput>> = None;
    let mut ercc_stats_stream: Option<Receiver<ParseOutput>> = None;
    let mut consensus_bam_stats_stream: Option<Receiver<ParseOutput>> = None;
    let mut consensus_bam_depth_stream: Option<Receiver<ParseOutput>> = None;
    let mut consensus_stats_stream: Option<Receiver<ParseOutput>> = None;
    let mut call_bcftools_stats_stream: Option<Receiver<ParseOutput>> = None;

    //*****************
    // Split by Technology
    match technology {
        Technology::Illumina => {
            eprintln!("Technology: Illumina");

            //*****************
            // Get Target reference sequence
            let (_filter_align_accession, filter_align_seq) = retrieve_h5_seq(
                config.args.ref_accession.clone(),
                config.args.ref_sequence.clone(),
                ref_db_path.as_ref(),
                h5_index.as_ref(),
            ).await?;
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
                config.base_buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::IlluminaFastq,
                "process_ercc".to_string(),
                None
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

            let (mut ercc_minimap2_child, ercc_minimap2_err_task) = spawn_cmd(config_arc.clone(),
                                                                              MINIMAP2_TAG,
                                                                              ercc_minimap2_args,
                                                                              config.args.verbose,
            ).await?;
            let ercc_minimap2_out_stream = parse_child_output(
                &mut ercc_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.base_buffer_size,
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

            let (mut ercc_samtools_child_view, ercc_samtools_task_view, ercc_samtools_err_task_view) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
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

            let ercc_stats_file_path = out_dir.join(format!("{}_ercc_stats.txt", no_ext_sample_base));
            let (mut ercc_samtools_child_stats, ercc_samtools_task_stats, ercc_samtools_err_task_stats) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size / 2,
            ).await?;
            cleanup_tasks.push(ercc_samtools_task_stats);
            cleanup_tasks.push(ercc_samtools_err_task_stats);

            let ercc_samtools_out_stream_stats = ReceiverStream::new(ercc_samtools_out_stream_stats);
            let (ercc_streams, ercc_done_rx) = t_junction(
                ercc_samtools_out_stream_stats,
                2,
                config.base_buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::JustBytes,
                "process_ercc_bypass".to_string(),
                None
            ).await?;
            cleanup_receivers.push(ercc_done_rx);

            let mut streams_iter = ercc_streams.into_iter();
            let ercc_file_stream = streams_iter.next().unwrap();
            ercc_stats_stream = Some(streams_iter.next().unwrap());

            match ercc_stats_stream {
                Some(stream) => {
                    // Convert byte stream to line stream
                    let (line_stream, line_task) = bytes_to_lines(stream, config.base_buffer_size / 2).await?;
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
            let target_ref_write_task = write_vecu8_to_file(filter_align_seq.clone(), target_ref_fasta_path.as_ref().unwrap(), config.base_buffer_size).await?;
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

                let (
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
                let (mut filter_minimap2_child, filter_minimap2_err_task) = spawn_cmd(config_arc.clone(),
                                                                                      MINIMAP2_TAG,
                                                                                      filter_minimap2_args,
                                                                                      config.args.verbose
                ).await?;
                let filter_minimap2_out_stream = parse_child_output(
                    &mut filter_minimap2_child,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    config.base_buffer_size,
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
                ) = stream_to_cmd(config_arc.clone(),
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
                    config.base_buffer_size,
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
                ) = stream_to_cmd(config_arc.clone(),
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
                    config.base_buffer_size,
                ).await?;
                let filter_samtools_out_stream_fastq = ReceiverStream::new(filter_samtools_out_stream_fastq);

                // Kraken2
                let kraken2_report_path = file_path_manipulator(
                    &PathBuf::from(&no_ext_sample_base_buf),
                    Some(&out_dir),
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
                let (_filter_kraken2_child, filter_kraken2_err_task) = spawn_cmd(config_arc.clone(),
                                                                                 KRAKEN2_TAG,
                                                                                 filter_reads_kraken2_args,
                                                                                 config.args.verbose
                ).await?;
                cleanup_tasks.push(filter_kraken2_err_task);

                let kraken2_classified_stream = TokioFile::open(&kraken2_classified_pipe_path).await?;
                let parse_rx = parse_fastq(kraken2_classified_stream, config.base_buffer_size).await?;

                let taxid = config.args.ref_taxid.as_ref().expect("ref_taxid should be Some");
                let pattern = format!("kraken:taxid|{}", taxid);
                let filter_fn = move |id: &str| id.contains(&pattern);

                let (filtered_rx, filter_task) = parse_and_filter_fastq_id(parse_rx, config.base_buffer_size, filter_fn.clone());
                cleanup_tasks.push(filter_task);

                let (parse_output_tx, parse_output_rx) = mpsc::channel(config.base_buffer_size);
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
                Some(&out_dir),
                None,
                Some("filtered.fq.gz"),
                "_"
            );

            let (align_streams, align_done_rx) = t_junction(
                filter_reads_out_stream,
                2,
                config.base_buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::IlluminaFastq,
                "filter_reads_output".to_string(),
                None
            )
                .await?;
            cleanup_receivers.push(align_done_rx);

            let mut streams_iter = align_streams.into_iter();
            let align_query_stream = streams_iter.next().unwrap();
            let align_file_stream = streams_iter.next().unwrap();
            let align_query_stream = ReceiverStream::new(align_query_stream);

            let align_pigz_args = generate_cli(PIGZ_TAG, &config, None)?;
            let (mut align_pigz_child, align_pigz_stream_task, align_pigz_err_task) = stream_to_cmd(config_arc.clone(), align_file_stream, PIGZ_TAG, align_pigz_args, StreamDataType::IlluminaFastq, config.args.verbose).await?;
            cleanup_tasks.push(align_pigz_stream_task);
            cleanup_tasks.push(align_pigz_err_task);

            let align_pigz_out_stream = parse_child_output(
                &mut align_pigz_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.base_buffer_size,
            ).await?;
            let align_pigz_write_task = tokio::spawn(stream_to_file(
                align_pigz_out_stream,
                align_fastq_path,
            ));
            cleanup_tasks.push(align_pigz_write_task);

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
            let (mut align_minimap2_child, align_minimap2_err_task) = spawn_cmd(config_arc.clone(),
                                                                                MINIMAP2_TAG,
                                                                                align_minimap2_args,
                                                                                config.args.verbose
            ).await?;
            cleanup_tasks.push(align_minimap2_err_task);
            let align_minimap2_out_stream = parse_child_output(
                &mut align_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                config.base_buffer_size,
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
            ) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
            ).await?;

            //*****************
            // Make Consensus
            align_bam_path = Some(file_path_manipulator(
                &no_ext_sample_base_buf,
                Some(&out_dir),
                None,
                Some("align.sam"),
                "_"
            ));
            let align_samtools_out_stream_sort = ReceiverStream::new(align_samtools_out_stream_sort);

            let (consensus_bam_streams, consensus_bam_done_rx) = t_junction(
                align_samtools_out_stream_sort,
                6,
                config.base_buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::JustBytes,
                "generate_consensus".to_string(),
                None
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
                Some(&out_dir),
                None,
                Some("consensus.fa"),
                "_"
            ));
            let (
                mut consensus_samtools_child,
                consensus_samtools_task_sort,
                consensus_samtools_err_task_sort
            ) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
            ).await?;

            let consensus_samtools_out_stream = ReceiverStream::new(consensus_samtools_out_stream);

            let (consensus_streams, consensus_done_rx) = t_junction(
                consensus_samtools_out_stream,
                4,
                config.base_buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::JustBytes,
                "generate_consensus".to_string(),
                None
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
            ) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
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
            ) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
            ).await?;

            let called_variants_path = file_path_manipulator(
                &no_ext_sample_base_buf,
                Some(&out_dir),
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
            ) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
            ).await?;

            let call_bcftools_out_stream_view = ReceiverStream::new(call_bcftools_out_stream_view);
            let (call_bcftools_out_streams, call_bcftools_done_rx) = t_junction(
                call_bcftools_out_stream_view,
                2,
                config.base_buffer_size,
                config.args.stall_threshold.try_into().unwrap(),
                Some(config.args.stream_sleep_ms),
                50,
                StreamDataType::JustBytes,
                "call_variants".to_string(),
                None
            )
                .await?;
            cleanup_receivers.push(call_bcftools_done_rx);

            let mut streams_iter = call_bcftools_out_streams.into_iter();
            let call_bcftools_file_stream = streams_iter.next().unwrap();
            call_bcftools_stats_stream = Some(streams_iter.next().unwrap());
            let called_variants_write_task = tokio::spawn(stream_to_file(
                call_bcftools_file_stream,
                called_variants_path
            ));
            cleanup_tasks.push(called_variants_write_task);

            //*****************
            // Realign Consensus
            let reference_file = TokioFile::open(target_ref_fasta_path.as_ref().unwrap()).await?;
            let reference_rx = parse_bytes(reference_file, config.base_buffer_size).await?;

            let realign_streams = vec![consensus_realign_stream, reference_rx];
            let (combined_rx, combined_task) = y_junction(realign_streams, config.base_buffer_size).await?;
            cleanup_tasks.push(combined_task);

            let realign_consensus_path = file_path_manipulator(
                &no_ext_sample_base_buf,
                Some(&out_dir),
                None,
                Some("consensus_realigned.fa"),
                "_"
            );
            let realign_mafft_args = generate_cli(MAFFT_TAG, &config, None)?;
            let (
                mut realign_consensus_mafft_child,
                realign_consensus_mafft_task,
                realign_consensus_mafft_err_task
            ) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size,
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
            let (mut stats_samtools_child_stats, stats_samtools_task_stats, stats_samtools_err_task_stats) = stream_to_cmd(config_arc.clone(),
                                                                                                                           stream,
                                                                                                                           SAMTOOLS_TAG,
                                                                                                                           stats_samtools_args_stats.clone(),
                                                                                                                           StreamDataType::JustBytes,
                                                                                                                           config.args.verbose,
            ).await?;
            cleanup_tasks.push(stats_samtools_task_stats);
            cleanup_tasks.push(stats_samtools_err_task_stats);
            parse_child_output(&mut stats_samtools_child_stats, ChildStream::Stdout, ParseMode::Lines, config.base_buffer_size / 2).await?
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
            let (mut depth_samtools_child, depth_samtools_task, depth_samtools_err_task) = stream_to_cmd(config_arc.clone(),
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
                config.base_buffer_size / 2,
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
    // NB: for this pipeline, for now at least, there should only be one chrom, so taking the first one
    let first_chr = depth_map.keys().next().ok_or_else(|| anyhow!("No chromosomes found"))?.clone();
    let first_chr_depth_map = depth_map.get(&first_chr).unwrap();
    let depths: Vec<u32> = first_chr_depth_map.values().copied().collect();
    let samtools_depth_stats = compute_depth_stats(&depths)?;

    let depth_plot_path = file_path_manipulator(
        &PathBuf::from(&no_ext_sample_base_buf),
        Some(&out_dir),
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

    let stats_file_path = out_dir.join(format!("{}_stats.json", no_ext_sample_base));
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

    let (_assembly_eval_quast_child, assembly_eval_quast_err_task) = spawn_cmd(config_arc.clone(), QUAST_TAG, assembly_eval_quast_args, config.args.verbose).await?;
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
