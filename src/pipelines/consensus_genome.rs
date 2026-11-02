//! Consensus genome pipeline

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;
use tokio_stream::StreamExt;

use anyhow::{anyhow, Result};
use futures::future::try_join_all;
use log::{self, debug, error, info};
use serde::Serialize;
use std::path::PathBuf;
use tempfile::{NamedTempFile, TempDir};
use tokio::fs::File as TokioFile;
use tokio::sync::mpsc::Receiver;
use tokio::sync::{mpsc, oneshot};
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;

use crate::cli::Technology;
use crate::config::defs::{
    BcftoolsSubcommand, PairingMode, PipelineError, SamtoolsSubcommand, SeqkitSubcommand,
    StreamDataType, BCFTOOLS_TAG, KRAKEN2_TAG, MAFFT_TAG, MINIMAP2_TAG, PIGZ_TAG, QUAST_TAG,
    SAMTOOLS_TAG, SEQKIT_TAG,
};
use crate::config::defs::{ReadStats, RunConfig};
use crate::utils::command::bcftools::BcftoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::command::minimap2::{minimap2_index_prep, Minimap2Config};
use crate::utils::command::quast::QuastConfig;
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::seqkit::SeqkitConfig;
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::fastx::{
    concatenate_paired_reads, parse_and_filter_fastq_id, parse_byte_stream_to_fastq, read_fastq,
    SequenceRecord
};

use crate::utils::file::{
    choose_temp_dir, file_path_manipulator, validate_file_inputs, write_parse_output_to_temp_fifo,
};
use crate::utils::plotting::plot_depths;
use crate::utils::stats::{
    compute_allele_counts, compute_coverage_bins, compute_depth_stats, parse_ercc_stats,
    parse_samtools_depth, parse_samtools_stats, parse_seqkit_stats,
};
use crate::utils::streams::ParseOutput;
use crate::utils::streams::{
    fanout_to_channels, join_with_error_handling, parse_bytes, parse_child_output, parse_fastq,
    spawn_cmd, stream_to_cmd, stream_to_file, y_junction, ChildStream, ParseMode, create_fifo
};
use crate::utils::vcf::count_variants_from_bcftools_stats;


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

/// Validates the input FASTQ files and returns an interleaved stream of sequence records.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `file1_path`: path to the R1 (or single-end) FASTQ file
/// * `file2_path`: optional path to the R2 FASTQ file
/// * `sample_base_buf`: base path for the sample
/// * `out_dir`: directory where output files should be saved
///
/// # Returns
///
/// Result containing:
/// - interleaved FASTQ stream
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - handle for the read stats task
async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        JoinHandle<anyhow::Result<ReadStats, anyhow::Error>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();
    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated.fq.gz"),
        "_",
    );

    let (rx, val_count_task) = read_fastq(
        file1_path,
        file2_path,
        PairingMode::Strict,
        Some(config.args.technology.clone()),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        "validate_input",
        &config,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);

    // Split the byte stream for output and quast write
    let (val_streams, val_router_handle) = fanout_to_channels(
        val_rx_stream,
        2,
        "validate_input",
        &config,
        StreamDataType::IlluminaFastq,
    )
        .await?;

    cleanup_receivers.push(val_router_handle);

    if val_streams.len() < 2 {
        return Err(PipelineError::EmptyStream);
    }

    // Consume Vec to get Receivers
    let mut val_streams_iter = val_streams.into_iter();
    let val_out_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let val_pigz_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let val_pigz_args =
        generate_cli(PIGZ_TAG, &config, None).map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (val_pigz_child, val_pigz_stream_task, val_pigz_err_task) = stream_to_cmd(
        config.clone(),
        val_pigz_stream,
        PIGZ_TAG,
        val_pigz_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(val_pigz_stream_task);
    cleanup_tasks.push(val_pigz_err_task);

    let val_pigz_out_stream = {
        let mut guard = val_pigz_child.lock().await;
        parse_child_output(
            &mut guard, // Pass locked Child
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: PIGZ_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let val_pigz_write_task = tokio::spawn(stream_to_file(
        val_pigz_out_stream,
        validated_interleaved_file_path,
    ));
    cleanup_tasks.push(val_pigz_write_task);

    Ok((
        ReceiverStream::new(val_out_stream),
        cleanup_tasks,
        cleanup_receivers,
        val_count_task,
    ))
}

/// Aligns reads against a host index to filter out host sequences.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `host_index_path`: path to the minimap2 host index
/// * `no_host_file_path`: path where non-host reads will be written
///
/// # Returns
///
/// Result containing:
/// - stream of non-host FASTQ records
/// - stream of host-aligned BAM records
/// - vector of cleanup tasks
/// - vector of cleanup receivers
async fn align_to_host(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>, // FASTQ raw byte stream
    host_index_path: PathBuf,                  // minimap2 .mmi index file
    no_host_file_path: PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let minimap_preset = match config.args.technology {
        Technology::Illumina => "sr",
        Technology::ONT => "map-ont",
    };

    let minimap2_config = Minimap2Config::new()
        .option(
            "-t",
            config
                .thread_allocation(MINIMAP2_TAG, None)
                .to_string(),
        )
        .option("-ax", minimap_preset)
        .positional(host_index_path.to_string_lossy().to_string())
        .positional("-");

    let minimap2_args =
        generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (minimap2_child, minimap2_stream_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(), // Convert to Receiver<ParseOutput> for streaming
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::JustBytes, // FASTQ bytes
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_stream_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Samtools view to filter unmapped reads (-f 4), output as uncompressed BAM (-b -u)
    let samtools_config_view = SamtoolsConfig::new(SamtoolsSubcommand::View)
        .option("-@", config.thread_allocation(SAMTOOLS_TAG, Some("view")).to_string())
        .flag("--no-PG")
        .option("-f", "4") // Unmapped reads
        .flag("-h") // Include header
        .flag("-b") // BAM output
        .flag("-u"); // Uncompressed BAM;
    let samtools_args_view = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (samtools_child_view, samtools_task_view, samtools_err_task_view) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_view);
    cleanup_tasks.push(samtools_err_task_view);

    let samtools_out_stream_view = {
        let mut guard = samtools_child_view.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Samtools fastq to convert uncompressed BAM to FASTQ bytes
    let samtools_config_fastq = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
        .option("-@", config.thread_allocation(SAMTOOLS_TAG, Some("fastq")).to_string())
        .option("-c", "6")
        .positional("-");
    let samtools_args_fastq = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_fastq))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (samtools_child_fastq, samtools_task_fastq, samtools_err_task_fastq) = stream_to_cmd(
        config.clone(),
        samtools_out_stream_view,
        SAMTOOLS_TAG,
        samtools_args_fastq,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_fastq);
    cleanup_tasks.push(samtools_err_task_fastq);

    let samtools_out_stream_fastq = {
        let mut guard = samtools_child_fastq.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    let samtools_out_stream_fastq = ReceiverStream::new(samtools_out_stream_fastq);

    // Split FASTQ byte stream for output, file write, count, and check
    let (host_streams, host_router_handle) = fanout_to_channels(
        samtools_out_stream_fastq,
        4,
        "align_to_host",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(host_router_handle);

    // Consume Vec to get Receivers
    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_count_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_check_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Check task for host-removed empty FASTQ (adapted for bytes)
    let (check_tx, mut check_rx) = mpsc::channel(1);
    let check_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(no_host_check_stream);
        let mut count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(bytes) => {
                    if !bytes.is_empty() {
                        count = 1;
                        // Spawn background drain for remaining items to avoid blocking
                        tokio::spawn(async move {
                            while stream.next().await.is_some() {}
                            Ok::<(), anyhow::Error>(())
                        });
                        break; // Exit loop early to unblock sender
                    }
                }
                _ => return Err(anyhow!("Unexpected item type in no_host_check_stream")),
            }
        }
        check_tx
            .send(count)
            .await
            .map_err(|e| anyhow!("Failed to send check count: {}", e))?;
        Ok(())
    });
    cleanup_tasks.push(check_task);

    let check_count = check_rx.recv().await.ok_or(PipelineError::EmptyStream)?;
    if check_count == 0 {
        let drain_rxs = vec![
            no_host_output_stream,
            no_host_file_stream,
            no_host_count_stream,
        ];
        let drain_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = drain_rxs
            .into_iter()
            .map(|rx| {
                tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(rx);
                    while stream.next().await.is_some() {}
                    Ok(())
                })
            })
            .collect();
        cleanup_tasks.extend(drain_tasks);
        return Err(PipelineError::EmptyStream);
    }

    // Pigz for compression
    let pigz_args =
        generate_cli(PIGZ_TAG, &config, None).map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (pigz_child, pigz_stream_task, pigz_err_task) = stream_to_cmd(
        config.clone(),
        no_host_file_stream,
        PIGZ_TAG,
        pigz_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(pigz_stream_task);
    cleanup_tasks.push(pigz_err_task);

    let pigz_out_stream = {
        let mut guard = pigz_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: PIGZ_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    let pigz_write_task = tokio::spawn(stream_to_file(pigz_out_stream, no_host_file_path));
    cleanup_tasks.push(pigz_write_task);

    Ok((
        ReceiverStream::new(no_host_output_stream),
        ReceiverStream::new(no_host_count_stream),
        cleanup_tasks,
        cleanup_receivers,
    ))
}

/// Quantifies ERCC spike-in sequences from a side-branch copy of the
/// host-removed FASTQ stream.
///
/// This mirrors the WDL QuantifyERCCs task:
///
///   host-removed FASTQs
///       -> minimap2 against ERCC reference
///       -> samtools view (BAM)
///       -> samtools stats
///
/// ERCC processing does NOT filter or transform the main read stream.
/// The caller is responsible for fanout before calling this function.
///
/// Zero ERCC mappings are valid and are handled naturally by samtools stats.
///
/// # Returns
///
/// - handle for parsed ERCC statistics
/// - cleanup tasks
/// - cleanup receivers
async fn process_ercc(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    ercc_index_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base: &str,
) -> Result<
    (
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

    // ------------------------------------------------------------
    // Align host-removed reads to the ERCC reference.
    //
    // This stream is already a dedicated ERCC branch created by
    // the caller. There is deliberately no bypass stream here.
    // ------------------------------------------------------------
    let minimap_preset = match config.args.technology {
        Technology::Illumina => "sr",
        Technology::ONT => "map-ont",
    };

    let minimap2_config = Minimap2Config::new()
        .option(
            "-t",
            config
                .thread_allocation(MINIMAP2_TAG, None)
                .to_string(),
        )
        .option("-ax", minimap_preset)
        .positional(ercc_index_path.to_string_lossy().to_string())
        .positional("-");

    let minimap2_args =
        generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (minimap2_child, minimap2_stream_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(minimap2_stream_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // ------------------------------------------------------------
    // Convert minimap2 SAM output to BAM.
    //
    // Match the WDL:
    //
    //   samtools view --no-PG -bo ercc_mapped.bam
    //
    // Do NOT filter for mapped reads here. samtools stats determines
    // the mapped-read counts from the BAM flags.
    // ------------------------------------------------------------
    let samtools_config_view =
        SamtoolsConfig::new(SamtoolsSubcommand::View)
            .option(
                "-@",
                config
                    .thread_allocation(SAMTOOLS_TAG, Some("view"))
                    .to_string(),
            )
            .flag("--no-PG")
            .flag("-b")
            .flag("-u");

    let samtools_args_view =
        generate_cli(
            SAMTOOLS_TAG,
            &config,
            Some(&samtools_config_view),
        )
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

    let (
        samtools_child_view,
        samtools_task_view,
        samtools_err_task_view,
    ) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(samtools_task_view);
    cleanup_tasks.push(samtools_err_task_view);

    let samtools_out_stream_view = {
        let mut guard = samtools_child_view.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // ------------------------------------------------------------
    // Generate ERCC statistics.
    //
    // This runs regardless of how many reads mapped to ERCC.
    // Zero mapped reads are a normal samtools-stats result, not an
    // empty-stream condition.
    // ------------------------------------------------------------
    let stats_samtools_config =
        SamtoolsConfig::new(SamtoolsSubcommand::Stats)
            .positional("-");

    let stats_samtools_args =
        generate_cli(
            SAMTOOLS_TAG,
            &config,
            Some(&stats_samtools_config),
        )
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

    let (
        stats_samtools_child,
        stats_samtools_task,
        stats_samtools_err_task,
    ) = stream_to_cmd(
        config.clone(),
        samtools_out_stream_view,
        SAMTOOLS_TAG,
        stats_samtools_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(stats_samtools_task);
    cleanup_tasks.push(stats_samtools_err_task);

    let stats_samtools_out_stream = {
        let mut guard = stats_samtools_child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Lines,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // ------------------------------------------------------------
    // samtools stats output has two consumers:
    //
    //   1. write <sample>_ercc_stats.txt
    //   2. parse the values needed by calculate_statistics()
    // ------------------------------------------------------------
    let (stats_streams, stats_router_handle) = fanout_to_channels(
        ReceiverStream::new(stats_samtools_out_stream),
        2,
        "process_ercc_stats",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(stats_router_handle);

    let mut stats_streams_iter = stats_streams.into_iter();

    let stats_file_stream = stats_streams_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    let stats_parse_stream = stats_streams_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    let stats_write_task = tokio::spawn(stream_to_file(
        stats_file_stream,
        ercc_stats_file_path,
    ));

    cleanup_tasks.push(stats_write_task);

    let ercc_stats_task =
        Some(tokio::spawn(parse_ercc_stats(stats_parse_stream)));

    Ok((
        ercc_stats_task,
        cleanup_tasks,
        cleanup_receivers,
    ))
}

/// Filters reads using Kraken2 to identify those belonging to a target taxon.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `out_dir`: directory where Kraken2 outputs will be written
/// * `no_ext_sample_base_buf`: sample base name without extension
/// * `target_taxid`: taxid of the target organism
///
/// # Returns
///
/// Result containing:
/// - stream of Kraken2-assigned FASTQ records
/// - vector of cleanup tasks
/// - vector of cleanup receivers
async fn filter_with_kraken(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    target_index_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
    target_taxid: &str,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<NamedTempFile>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = vec![];
    let cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();
    let mut temp_dirs: Vec<TempDir> = Vec::new();
    let named_temp_files: Vec<NamedTempFile> = Vec::new();

    // -------------------------------------------------------------------------
    // WDL FilterReads increment 2:
    //
    // host-removed reads
    //   -> minimap2 against target
    //   -> samtools sort -n
    //   -> samtools fastq -G 12 -0 /dev/null -s /dev/null -N
    //   -> route explicit /1 and /2 mate identities
    //   -> R1 FIFO + R2 FIFO
    //   -> kraken2 --paired
    //
    // samtools retains only complete category-1/category-2 pairs:
    //
    //   -G 12         exclude records where read + mate are both unmapped
    //   -0 /dev/null  discard category-0 records
    //   -s /dev/null  discard singleton category-1/category-2 records
    //
    // -N is used instead of the WDL's -n only so Rust can identify the
    // category-1/category-2 output records as /1 and /2. Those temporary
    // suffixes are removed before Kraken.
    //
    // No mate identity is inferred from record position or adjacency.
    //
    // Increment 1 validates paired Kraken input accounting.
    // Increment 2 additionally validates Kraken paired classified output.
    // This test still stops intentionally before downstream handoff.
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // minimap2 -ax sr TARGET -
    // -------------------------------------------------------------------------

    let minimap2_config = Minimap2Config::new()
        .option(
            "-t",
            config
                .thread_allocation(MINIMAP2_TAG, None)
                .to_string(),
        )
        .option("-ax", "sr")
        .positional(target_index_path.to_string_lossy().to_string())
        .positional("-");

    let minimap2_args =
        generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config))
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?;

    let (
        minimap2_child,
        minimap2_stream_task,
        minimap2_err_task,
    ) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(minimap2_stream_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // -------------------------------------------------------------------------
    // samtools sort -n
    // -------------------------------------------------------------------------

    let sort_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    let sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("sort"))
                .to_string(),
        )
        .option("-m", "2G")
        .flag("-n")
        .flag("-u")
        .option("-O", "bam")
        .option(
            "-T",
            sort_temp_dir.path().to_string_lossy().to_string(),
        )
        .positional("-");

    let sort_args =
        generate_cli(SAMTOOLS_TAG, &config, Some(&sort_config))
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

    let (
        sort_child,
        sort_stream_task,
        sort_err_task,
    ) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(sort_stream_task);
    cleanup_tasks.push(sort_err_task);

    let sorted_bam_stream = {
        let mut guard = sort_child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    temp_dirs.push(sort_temp_dir);

    // -------------------------------------------------------------------------
    // samtools fastq
    //
    // WDL paired-end selection semantics:
    //
    //   -G 12
    //   -0 /dev/null
    //   -s /dev/null
    //
    // WDL uses -1/-2 to separate the surviving complete pairs.
    //
    // For this FIFO test we instead allow those surviving category-1 and
    // category-2 records to stdout and use -N to expose /1 and /2 so Rust
    // can route them to independently-owned FIFO writers.
    // -------------------------------------------------------------------------

    let fastq_config = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("fastq"))
                .to_string(),
        )
        .option("-G", "12")
        .option("-0", "/dev/null")
        .option("-s", "/dev/null")
        .flag("-N")
        .option("-c", "6")
        .positional("-");

    let fastq_args =
        generate_cli(SAMTOOLS_TAG, &config, Some(&fastq_config))
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

    debug!(
        "FilterReads starting samtools fastq -G 12 -0 /dev/null -s /dev/null -N -> explicit mate router"
    );

    let (
        fastq_child,
        fastq_stream_task,
        fastq_err_task,
    ) = stream_to_cmd(
        config.clone(),
        sorted_bam_stream,
        SAMTOOLS_TAG,
        fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let fastq_out_stream = {
        let mut guard = fastq_child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // -------------------------------------------------------------------------
    // Create R1/R2 FIFOs.
    // -------------------------------------------------------------------------

    let fifo_temp_dir = tempfile::Builder::new()
        .prefix("consensus_filter_reads_")
        .tempdir_in(&config.ram_temp_dir)
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to create FilterReads FIFO directory in {}: {}",
                config.ram_temp_dir.display(),
                e
            ))
        })?;

    let fq1_fifo = fifo_temp_dir.path().join("paired1.fq");
    let fq2_fifo = fifo_temp_dir.path().join("paired2.fq");

    create_fifo(&fq1_fifo)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    create_fifo(&fq2_fifo)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    debug!(
        "FilterReads paired query FIFOs created: R1={} R2={}",
        fq1_fifo.display(),
        fq2_fifo.display()
    );

    temp_dirs.push(fifo_temp_dir);

    // -------------------------------------------------------------------------
    // Explicit mate router.
    //
    // Every surviving record is routed according to the /1 or /2 marker
    // produced by samtools fastq -N.
    // -------------------------------------------------------------------------

    let fifo_channel_buffer = crate::utils::system::compute_buffer_size(
        &config,
        "filter_with_kraken_paired_router",
        StreamDataType::IlluminaFastq,
        1.6,
    );

    let (r1_tx, r1_rx) = mpsc::channel::<ParseOutput>(fifo_channel_buffer);
    let (r2_tx, r2_rx) = mpsc::channel::<ParseOutput>(fifo_channel_buffer);

    let mate_router_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(fastq_out_stream);

        let mut total_count: u64 = 0;
        let mut r1_count: u64 = 0;
        let mut r2_count: u64 = 0;

        while let Some(item) = stream.next().await {
            let record = match item {
                ParseOutput::Fastq(record) => record,
                _ => {
                    return Err(anyhow!(
                        "FilterReads paired router received non-FASTQ item"
                    ));
                }
            };

            let (
                id,
                desc,
                seq,
                qual,
            ) = match record {
                SequenceRecord::Fastq {
                    id,
                    desc,
                    seq,
                    qual,
                } => (id, desc, seq, qual),
                _ => {
                    return Err(anyhow!(
                        "FilterReads paired router received non-FASTQ SequenceRecord"
                    ));
                }
            };

            total_count += 1;

            if let Some(base_id) = id.strip_suffix("/1") {
                let output = ParseOutput::Fastq(
                    SequenceRecord::Fastq {
                        id: base_id.to_string(),
                        desc,
                        seq,
                        qual,
                    }
                );

                r1_tx
                    .send(output)
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "FilterReads R1 FIFO channel closed after {} R1 records",
                            r1_count
                        )
                    })?;

                r1_count += 1;
            } else if let Some(base_id) = id.strip_suffix("/2") {
                let output = ParseOutput::Fastq(
                    SequenceRecord::Fastq {
                        id: base_id.to_string(),
                        desc,
                        seq,
                        qual,
                    }
                );

                r2_tx
                    .send(output)
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "FilterReads R2 FIFO channel closed after {} R2 records",
                            r2_count
                        )
                    })?;

                r2_count += 1;
            } else {
                return Err(anyhow!(
                    "FilterReads FASTQ record lacks explicit /1 or /2 mate marker: '{}'",
                    id
                ));
            }

            if total_count % 100_000 == 0 {
                debug!(
                    "FilterReads mate router: {} records routed (R1={}, R2={})",
                    total_count,
                    r1_count,
                    r2_count
                );
            }
        }

        drop(r1_tx);
        drop(r2_tx);

        if r1_count != r2_count {
            return Err(anyhow!(
                "FilterReads mate imbalance after samtools fastq: R1={} R2={} total={}",
                r1_count,
                r2_count,
                total_count
            ));
        }

        if total_count != r1_count + r2_count {
            return Err(anyhow!(
                "FilterReads mate accounting mismatch: total={} R1={} R2={}",
                total_count,
                r1_count,
                r2_count
            ));
        }

        debug!(
            "FilterReads mate router completed cleanly: total={} R1={} R2={}",
            total_count,
            r1_count,
            r2_count
        );

        Ok::<(), anyhow::Error>(())
    });

    // -------------------------------------------------------------------------
    // Independent FIFO writer tasks.
    // -------------------------------------------------------------------------

    let r1_fifo_for_writer = fq1_fifo.clone();
    let r2_fifo_for_writer = fq2_fifo.clone();

    let r1_write_task = tokio::spawn(async move {
        crate::utils::streams::write_to_fifo(
            r1_rx,
            r1_fifo_for_writer,
        )
            .await
    });

    let r2_write_task = tokio::spawn(async move {
        crate::utils::streams::write_to_fifo(
            r2_rx,
            r2_fifo_for_writer,
        )
            .await
    });

    // -------------------------------------------------------------------------
    // Kraken paired classified output goes to temporary files.
    //
    // Kraken requires '#' in --classified-out when using --paired.
    // -------------------------------------------------------------------------

    let classified_temp_dir = tempfile::Builder::new()
        .prefix("consensus_kraken_classified_")
        .tempdir_in(&config.ram_temp_dir)
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to create Kraken classified temp dir in {}: {}",
                config.ram_temp_dir.display(),
                e
            ))
        })?;

    let classified_path =
        classified_temp_dir.path().join("classified#.fq");
    let classified_r1_path =
        classified_temp_dir.path().join("classified_1.fq");
    let classified_r2_path =
        classified_temp_dir.path().join("classified_2.fq");
    temp_dirs.push(classified_temp_dir);

    // -------------------------------------------------------------------------
    // Kraken2 --paired R1_FIFO R2_FIFO
    // -------------------------------------------------------------------------

    let kraken_report_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("kraken2_report.txt"),
        "_",
    );

    let kraken2_config = Kraken2Config {
        report_path: kraken_report_path.clone(),
        classified_path: classified_path.clone(),
        fastq_path: fq1_fifo.clone(),
        fastq_path_2: Some(fq2_fifo.clone()),
    };

    let kraken2_args =
        generate_cli(
            KRAKEN2_TAG,
            &config,
            Some(&kraken2_config),
        )
            .map_err(|e| PipelineError::ToolExecution {
                tool: KRAKEN2_TAG.to_string(),
                error: e.to_string(),
            })?;

    debug!(
        "FilterReads starting paired Kraken2: R1={} R2={}",
        fq1_fifo.display(),
        fq2_fifo.display()
    );

    let (
        mut kraken2_child,
        kraken2_err_task,
    ) = spawn_cmd(
        config.clone(),
        KRAKEN2_TAG,
        kraken2_args,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KRAKEN2_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(kraken2_err_task);

    // -------------------------------------------------------------------------
    // Drain Kraken normal stdout.
    // -------------------------------------------------------------------------

    let kraken_stdout = kraken2_child
        .stdout
        .take()
        .ok_or_else(|| {
            PipelineError::Other(anyhow!(
                "Failed to capture Kraken2 stdout"
            ))
        })?;

    let stdout_drain_task = tokio::spawn(async move {
        use tokio::io::{AsyncBufReadExt, BufReader};

        let mut lines = BufReader::new(kraken_stdout).lines();
        let mut count = 0u64;

        let target = "\t2697049\t";
        let mut target_count = 0u64;

        while let Some(line) = lines.next_line().await? {
            if count < 10 {
                debug!("KRAKEN CLASSIFICATION OUTPUT: {}", line);
            }

            if line.contains(target) {
                target_count += 1;

                if target_count <= 5 {
                    debug!(
                "KRAKEN TARGET CLASSIFICATION: {}",
                line
            );
                }
            }

            count += 1;
        }

        debug!(
    "Kraken classification output: {} lines, {} assigned directly to target taxid",
    count,
    target_count
);

        debug!("Kraken classification output: {} lines", count);

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(stdout_drain_task);

    // -------------------------------------------------------------------------
    // Wait for paired Kraken.
    // -------------------------------------------------------------------------

    let kraken_status = kraken2_child
        .wait()
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KRAKEN2_TAG.to_string(),
            error: format!(
                "Failed waiting for paired Kraken2: {}",
                e
            ),
        })?;

    if !kraken_status.success() {
        return Err(PipelineError::ToolExecution {
            tool: KRAKEN2_TAG.to_string(),
            error: format!(
                "Paired Kraken2 exited with status {:?}",
                kraken_status.code()
            ),
        });
    }

    // -------------------------------------------------------------------------
    // TEST ACCOUNTING CHECK
    // -------------------------------------------------------------------------

    debug!(
        "FilterReads accounting test: Kraken exited successfully; \
         waiting for mate router and both FIFO writers"
    );

    match mate_router_task.await {
        Ok(Ok(())) => {
            debug!(
                "FilterReads accounting test: mate router completed successfully"
            );
        }
        Ok(Err(e)) => {
            return Err(PipelineError::Other(anyhow!(
                "FilterReads mate router failed: {}",
                e
            )));
        }
        Err(e) => {
            return Err(PipelineError::Other(anyhow!(
                "FilterReads mate router task join failed: {}",
                e
            )));
        }
    }

    match r1_write_task.await {
        Ok(Ok(())) => {
            debug!(
                "FilterReads accounting test: R1 FIFO writer completed successfully"
            );
        }
        Ok(Err(e)) => {
            return Err(PipelineError::Other(anyhow!(
                "FilterReads R1 FIFO writer failed: {}",
                e
            )));
        }
        Err(e) => {
            return Err(PipelineError::Other(anyhow!(
                "FilterReads R1 FIFO writer task join failed: {}",
                e
            )));
        }
    }

    match r2_write_task.await {
        Ok(Ok(())) => {
            debug!(
                "FilterReads accounting test: R2 FIFO writer completed successfully"
            );
        }
        Ok(Err(e)) => {
            return Err(PipelineError::Other(anyhow!(
                "FilterReads R2 FIFO writer failed: {}",
                e
            )));
        }
        Err(e) => {
            return Err(PipelineError::Other(anyhow!(
                "FilterReads R2 FIFO writer task join failed: {}",
                e
            )));
        }
    }

    info!(
        "FilterReads increment-1 ACCOUNTING SUCCESS: \
         mate router + R1 writer + R2 writer + paired Kraken2 all completed; \
         report={}",
        kraken_report_path.display()
    );

    // -------------------------------------------------------------------------
    // Increment 2: validate paired Kraken classified output.
    //
    // With --paired, Kraken expands classified#.fq into:
    //
    //   classified_1.fq
    //   classified_2.fq
    //
    // Kraken has already exited, so both files are complete and can be
    // reopened without introducing another FIFO lifecycle into this test.
    // Consume the two FASTQ streams in lockstep and require:
    //
    //   * equal numbers of R1 and R2 classified records
    //   * matching mate IDs for every classified pair
    //   * matching Kraken taxid annotations for both mates
    //
    // Also count the exact target-taxid pairs that the legacy implementation
    // would retain via parse_and_filter_fastq_id().
    // -------------------------------------------------------------------------

    let classified_r1_file = tokio::fs::File::open(&classified_r1_path)
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to open Kraken classified R1 FASTQ {}: {}",
                classified_r1_path.display(),
                e
            ))
        })?;

    let classified_r2_file = tokio::fs::File::open(&classified_r2_path)
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to open Kraken classified R2 FASTQ {}: {}",
                classified_r2_path.display(),
                e
            ))
        })?;

    let classified_r1_rx = parse_fastq(
        classified_r1_file,
        &config,
        StreamDataType::IlluminaFastq,
    )
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to parse Kraken classified R1 FASTQ: {}",
                e
            ))
        })?;

    let classified_r2_rx = parse_fastq(
        classified_r2_file,
        &config,
        StreamDataType::IlluminaFastq,
    )
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to parse Kraken classified R2 FASTQ: {}",
                e
            ))
        })?;

    let mut classified_r1_stream = ReceiverStream::new(classified_r1_rx);
    let mut classified_r2_stream = ReceiverStream::new(classified_r2_rx);

    let target_pattern = format!("kraken:taxid|{}", target_taxid);

    let mut classified_pairs: u64 = 0;
    let mut target_pairs: u64 = 0;

    loop {
        let r1_item = classified_r1_stream.next().await;
        let r2_item = classified_r2_stream.next().await;

        match (r1_item, r2_item) {
            (None, None) => break,

            (Some(_), None) => {
                return Err(PipelineError::Other(anyhow!(
                    "Kraken classified mate imbalance: R1 contains more records than R2 after {} complete pairs",
                    classified_pairs
                )));
            }

            (None, Some(_)) => {
                return Err(PipelineError::Other(anyhow!(
                    "Kraken classified mate imbalance: R2 contains more records than R1 after {} complete pairs",
                    classified_pairs
                )));
            }

            (Some(r1_item), Some(r2_item)) => {
                let r1_record = match r1_item {
                    ParseOutput::Fastq(record) => record,
                    _ => {
                        return Err(PipelineError::Other(anyhow!(
                            "Kraken classified R1 parser emitted non-FASTQ item"
                        )));
                    }
                };

                let r2_record = match r2_item {
                    ParseOutput::Fastq(record) => record,
                    _ => {
                        return Err(PipelineError::Other(anyhow!(
                            "Kraken classified R2 parser emitted non-FASTQ item"
                        )));
                    }
                };

                let r1_id = r1_record.id();
                let r2_id = r2_record.id();




                if classified_pairs < 100 {
                    debug!(
            "FilterReads classified pair {}: R1 id={:?}",
            classified_pairs + 1,
            r1_record.id(),
        );

                    debug!(
            "FilterReads classified pair {}: R2 id={:?}",
            classified_pairs + 1,
            r2_record.id(),
        );
                }





                if r1_id != r2_id {
                    return Err(PipelineError::Other(anyhow!(
                        "Kraken classified mate ID mismatch at pair {}: R1='{}' R2='{}'",
                        classified_pairs + 1,
                        r1_id,
                        r2_id
                    )));
                }

                let r1_is_target = r1_id.contains(&target_pattern);
                let r2_is_target = r2_id.contains(&target_pattern);

                if r1_is_target != r2_is_target {
                    return Err(PipelineError::Other(anyhow!(
                        "Kraken classified target-taxid disagreement at pair {}: R1='{}' R2='{}' target_taxid={}",
                        classified_pairs + 1,
                        r1_id,
                        r2_id,
                        target_taxid
                    )));
                }

                classified_pairs += 1;

                if r1_is_target {
                    target_pairs += 1;
                }

                if classified_pairs % 100_000 == 0 {
                    debug!(
                        "FilterReads classified-output validation: {} classified pairs, {} exact target-taxid pairs",
                        classified_pairs,
                        target_pairs
                    );
                }
            }
        }
    }

    info!(
        "FilterReads increment-2 CLASSIFIED OUTPUT SUCCESS: \
         classified_pairs={} classified_records={} exact_target_taxid_pairs={} \
         exact_target_taxid_records={} target_taxid={} R1={} R2={}",
        classified_pairs,
        classified_pairs * 2,
        target_pairs,
        target_pairs * 2,
        target_taxid,
        classified_r1_path.display(),
        classified_r2_path.display()
    );

    // -------------------------------------------------------------------------
    // TEST STOP
    // -------------------------------------------------------------------------

    Err(PipelineError::Other(anyhow!(
        "FilterReads increment-2 paired classified-output validation completed successfully; \
         intentional stop before target-pair downstream handoff"
    )))
}


/// Aligns reads against a target organism's reference index and sorts the resulting BAM.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `target_index_path`: path to the minimap2 target index
/// * `out_dir`: directory where alignment outputs will be saved
/// * `no_ext_sample_base_buf`: sample base name without extension
///
/// # Returns
///
/// Result containing:
/// - stream of target-aligned BAM records
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - vector of QUAST write tasks
/// - path to the generated BAM file
async fn align_to_target(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>, // FASTQ SequenceRecord stream
    target_index_path: PathBuf,                // Minimap2 .mmi index
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        PathBuf,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];
    let mut quast_write_tasks = vec![];

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    let minimap_preset = match config.args.technology {
        Technology::Illumina => "sr",
        Technology::ONT => "map-ont",
    };

    let minimap2_config = Minimap2Config::new()
        .option(
            "-t",
            config
                .thread_allocation(MINIMAP2_TAG, None)
                .to_string(),
        )
        .option("-ax", minimap_preset)
        .positional(target_index_path.to_string_lossy().to_string())
        .positional("-");

    let minimap2_args =
        generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (minimap2_child, _minimap2_stream_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(), // Convert to Receiver<ParseOutput> for streaming
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::JustBytes, // FASTQ bytes, not structured records
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option("-@", config.thread_allocation(SAMTOOLS_TAG, Some("sort")).to_string())
        .option("-m", "2G")
        .flag("-u") // Uncompressed
        .option("-O", "bam") // BAM output
        .option("-T", temp_dir.path().to_string_lossy().to_string())
        .positional("-");
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_out_stream = ReceiverStream::new(samtools_sort_out_stream);

    let align_bam_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("target_aligned.bam"),
        "_",
    );

    let (sam_streams, sam_router_handle) = fanout_to_channels(
        samtools_sort_out_stream,
        3,
        "align_to_target",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(sam_router_handle);
    let mut streams_iter = sam_streams.into_iter();
    let sam_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_check_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    //  check stream convert to SAM for non-empty chk
    let sam_view_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&SamtoolsConfig::new(SamtoolsSubcommand::View)
            .option("-@", config.thread_allocation(SAMTOOLS_TAG, Some("view")).to_string())
            .flag("--no-PG")
            .flag("-h")
            .positional("-")),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (sam_view_child, sam_view_task, sam_view_err_task) = stream_to_cmd(
        config.clone(),
        sam_check_stream,
        SAMTOOLS_TAG,
        sam_view_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(sam_view_task);
    cleanup_tasks.push(sam_view_err_task);

    let sam_check_stream = {
        let mut guard = sam_view_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    //check for alignment records (not just ehader)
    let (check_tx, mut check_rx) = mpsc::channel(1);
    let check_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(sam_check_stream);
        let mut alignment_count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(line) => {
                    let line_str = String::from_utf8_lossy(&line);
                    if !line_str.starts_with('@') && !line_str.trim().is_empty() {
                        alignment_count = 1;
                        // Spawn background drain to avoid blocking
                        tokio::spawn(async move {
                            while stream.next().await.is_some() {}
                            Ok::<(), anyhow::Error>(())
                        });
                        break; // Exit early to unblock sender
                    }
                }
                _ => return Err(anyhow!("Unexpected item type in sam_check_stream")),
            }
        }
        check_tx
            .send(alignment_count)
            .await
            .map_err(|e| anyhow!("Failed to send alignment count: {}", e))?;
        Ok::<(), anyhow::Error>(())
    });
    cleanup_tasks.push(check_task);

    let alignment_count = check_rx.recv().await.ok_or(PipelineError::EmptyStream)?;
    if alignment_count == 0 {
        let drain_rxs = vec![sam_output_stream, sam_file_stream];
        let drain_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = drain_rxs
            .into_iter()
            .map(|rx| {
                tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(rx);
                    while stream.next().await.is_some() {}
                    Ok(())
                })
            })
            .collect();
        cleanup_tasks.extend(drain_tasks);
        return Err(PipelineError::Other(anyhow!(
            "No alignment records in BAM stream"
        )));
    }

    let bam_write_task = tokio::spawn(stream_to_file(sam_file_stream, align_bam_path.clone()));
    quast_write_tasks.push(bam_write_task);

    Ok((
        ReceiverStream::new(sam_output_stream),
        cleanup_tasks,
        cleanup_receivers,
        quast_write_tasks,
        align_bam_path,
    ))
}

/// Generates a consensus sequence from a BAM stream of aligned reads.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `bam_stream`: stream of uncompressed BAM records
/// * `out_dir`: directory where consensus outputs will be written
/// * `no_ext_sample_base_buf`: sample base name without extension
///
/// # Returns
///
/// Result containing:
/// - stream for consensus realignment
/// - stream for consensus statistics
/// - path to the generated consensus FASTA file
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - vector of additional consensus tasks
async fn generate_consensus(
    config: Arc<RunConfig>,
    bam_stream: ReceiverStream<ParseOutput>, // Uncompressed BAM byte stream
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        ReceiverStream<ParseOutput>,
        PathBuf,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];
    let mut quast_write_tasks = vec![];

    // Add primer trimming with samtools ampliconclip
    let primer_bed_path =
        config.args.primer_bed_path.clone().ok_or_else(|| {
            PipelineError::InvalidConfig("Primer BED file not provided".to_string())
        })?;
    let trim_config = SamtoolsConfig::new(SamtoolsSubcommand::Ampliconclip)
        .option("-@", config.thread_allocation(SAMTOOLS_TAG, Some("ampliconclip")).to_string())
        .option("-b", primer_bed_path)
        .flag("--both-ends")
        .option("-o", "-")
        .positional("-");
    let trim_args = generate_cli(SAMTOOLS_TAG, &config, Some(&trim_config)).map_err(|e| {
        PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        }
    })?;

    let (trim_child, trim_task, trim_err_task) = stream_to_cmd(
        config.clone(),
        bam_stream.into_inner(),
        SAMTOOLS_TAG,
        trim_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(trim_task);
    cleanup_tasks.push(trim_err_task);

    let trimmed_bam_stream = {
        let mut guard = trim_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_consensus_config = SamtoolsConfig::new(SamtoolsSubcommand::Consensus)
        .option("-f", "fasta")
        .option("-m", "simple")
        .option("--min-BQ", config.args.quality.to_string())
        .option("-d", config.args.min_depth.to_string())
        .positional("-");
    let samtools_consensus_args =
        generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_consensus_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (samtools_consensus_child, samtools_consensus_task, samtools_consensus_err_task) =
        stream_to_cmd(
            config.clone(),
            ReceiverStream::new(trimmed_bam_stream).into_inner(),
            SAMTOOLS_TAG,
            samtools_consensus_args,
            StreamDataType::JustBytes,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;
    cleanup_tasks.push(samtools_consensus_task);
    cleanup_tasks.push(samtools_consensus_err_task);

    let samtools_consensus_out_stream = {
        let mut guard = samtools_consensus_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Fasta, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_consensus_out_stream = ReceiverStream::new(samtools_consensus_out_stream);

    let (consensus_streams, consensus_router_handle) = fanout_to_channels(
        samtools_consensus_out_stream,
        3,
        "generate_consensus",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(consensus_router_handle);
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
        cleanup_receivers,
        quast_write_tasks,
    ))
}

/// Calls variants from a BAM stream and a target reference FASTA.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `bam_stream`: stream of uncompressed BAM records
/// * `target_ref_path`: path to the target reference FASTA file
/// * `out_dir`: directory where variant calling outputs will be written
/// * `no_ext_sample_base_buf`: sample base name without extension
///
/// # Returns
///
/// Result containing:
/// - stream of variant calling statistics
/// - path to the generated VCF file
/// - vector of cleanup tasks
/// - vector of cleanup receivers
async fn call_variants(
    config: Arc<RunConfig>,
    bam_stream: ReceiverStream<ParseOutput>, //Uncompressed BAM byte stream
    target_ref_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        PathBuf,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];

    let bcftools_mpileup_config =
        BcftoolsConfig::new(BcftoolsSubcommand::Mpileup)
            .option("-a", "AD")
            .option("-d", "10000")
            .option("-L", "100000000")
            .option("-Q", config.args.quality.to_string())
            .option(
                "--threads",
                RunConfig::thread_allocation(
                    &config,
                    BCFTOOLS_TAG,
                    Some("mpileup"),
                )
                    .to_string(),
            )
            .option(
                "-f",
                target_ref_path.to_string_lossy().into_owned(),
            )
            .flag("-Ou")
            .positional("-");
    let bcftools_mpileup_args = generate_cli(BCFTOOLS_TAG, &config, Some(&bcftools_mpileup_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (bcftools_mpileup_child, bcftools_mpileup_task, bcftools_mpileup_err_task) = stream_to_cmd(
        config.clone(),
        bam_stream.into_inner(),
        BCFTOOLS_TAG,
        bcftools_mpileup_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bcftools_mpileup_task);
    cleanup_tasks.push(bcftools_mpileup_err_task);

    let bcftools_mpileup_out_stream = {
        let mut guard = bcftools_mpileup_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BCFTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let bcftools_call_config =
        BcftoolsConfig::new(BcftoolsSubcommand::Call)
            .option("--ploidy", "1")
            .flag("-m")
            .flag("-v")
            .option(
                "-P",
                config.args.bcftools_call_theta.to_string(),
            )
            .flag("-Ou")
            .positional("-");

    let bcftools_call_args = generate_cli(BCFTOOLS_TAG, &config, Some(&bcftools_call_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (bcftools_call_child, bcftools_call_task, bcftools_call_err_task) = stream_to_cmd(
        config.clone(),
        bcftools_mpileup_out_stream,
        BCFTOOLS_TAG,
        bcftools_call_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bcftools_call_task);
    cleanup_tasks.push(bcftools_call_err_task);

    // Parse the bcftools call output
    let bcftools_call_out_stream = {
        let mut guard = bcftools_call_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BCFTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    let bcftools_call_out_stream = ReceiverStream::new(bcftools_call_out_stream);

    let (vcf_streams, vcf_router_handle) = fanout_to_channels(
        bcftools_call_out_stream,
        2,
        "call_variants",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(vcf_router_handle);

    let mut streams_iter = vcf_streams.into_iter();
    let vcf_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let vcf_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let vcf_file_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("variants.bcf"),
        "_",
    );

    let vcf_write_task = tokio::spawn(stream_to_file(vcf_file_stream, vcf_file_path.clone()));
    cleanup_tasks.push(vcf_write_task);

    Ok((
        ReceiverStream::new(vcf_output_stream),
        vcf_file_path,
        cleanup_tasks,
        cleanup_receivers,
    ))
}

/// Realigns the generated consensus sequence back to the target reference.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `consensus_realign_stream`: receiver for consensus sequence bytes
/// * `target_ref_fasta_path`: path to the target reference FASTA file
/// * `out_dir`: directory where realignment outputs will be written
/// * `no_ext_sample_base_buf`: sample base name without extension
///
/// # Returns
///
/// Result containing a vector of cleanup tasks
async fn realign_consensus_to_ref(
    config: Arc<RunConfig>,
    consensus_realign_stream: Receiver<ParseOutput>, // FASTA
    target_ref_fasta_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<Vec<JoinHandle<Result<(), anyhow::Error>>>, PipelineError> {
    let reference_file = TokioFile::open(&target_ref_fasta_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    let reference_rx = parse_bytes(reference_file, &config, StreamDataType::JustBytes)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let realign_streams = vec![consensus_realign_stream, reference_rx];
    let (combined_rx, combined_task) = y_junction(realign_streams, config.base_buffer_size)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let realign_consensus_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("consensus_realigned.fa"),
        "_",
    );
    let realign_mafft_args =
        generate_cli(MAFFT_TAG, &config, None).map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (
        realign_consensus_mafft_child,
        realign_consensus_mafft_task,
        realign_consensus_mafft_err_task,
    ) = stream_to_cmd(
        config.clone(),
        combined_rx,
        MAFFT_TAG,
        realign_mafft_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;

    let mut cleanup_tasks = vec![
        combined_task,
        realign_consensus_mafft_task,
        realign_consensus_mafft_err_task,
    ];

    let consensus_samtools_out_stream = {
        let mut guard = realign_consensus_mafft_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MAFFT_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let realign_consensus_write_task = tokio::spawn(stream_to_file(
        consensus_samtools_out_stream,
        realign_consensus_path,
    ));
    cleanup_tasks.push(realign_consensus_write_task);

    Ok(cleanup_tasks)
}

/// Calculates pipeline statistics and writes them to a JSON file.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `no_ext_sample_base`: sample base name without extension
/// * `consensus_bam_stats_stream`: optional stream of BAM statistics
/// * `consensus_bam_depth_stream`: optional stream of BAM depth information
/// * `no_host_seqkit_out_stream_stats`: stream of Seqkit stats for non-host reads
/// * `ercc_stats_task`: optional handle for the ERCC count task
/// * `consensus_stats_stream`: optional stream of consensus sequence statistics
/// * `call_bcftools_stats_stream`: optional stream of variant calling statistics
/// * `out_dir`: directory where statistics JSON will be written
/// * `technology`: the sequencing technology used
///
/// # Returns
///
/// Result<()>: success or error
async fn calculate_statistics(
    config: Arc<RunConfig>,
    no_ext_sample_base: &str,
    consensus_bam_stats_stream: Option<ReceiverStream<ParseOutput>>, //Uncompressed BAM byte stream
    consensus_bam_depth_stream: Option<ReceiverStream<ParseOutput>>, //Uncompressed BAM byte stream
    no_host_seqkit_out_stream_stats: Receiver<ParseOutput>,
    ercc_stats_task: Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>>,
    consensus_stats_stream: Option<ReceiverStream<ParseOutput>>,
    call_bcftools_stats_stream: Option<ReceiverStream<ParseOutput>>, // Uncompressed BCF stream (-Ou)
    out_dir: &PathBuf,
    technology: Technology,
) -> Result<(), anyhow::Error> {
    let mut local_cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();

    let stats_samtools_config_stats = SamtoolsConfig::new(SamtoolsSubcommand::Stats)
        .positional("-");
    let stats_samtools_args_stats =
        generate_cli(SAMTOOLS_TAG, &config, Some(&stats_samtools_config_stats))?;

    let stats_samtools_out_stream_stats = match consensus_bam_stats_stream {
        Some(stream) => {
            let (
                stats_samtools_child_stats,
                stats_samtools_task_stats,
                stats_samtools_err_task_stats,
            ) = stream_to_cmd(
                config.clone(),
                stream.into_inner(),
                SAMTOOLS_TAG,
                stats_samtools_args_stats.clone(),
                StreamDataType::JustBytes,
                config.args.verbose,
                None,
            )
                .await?;
            local_cleanup_tasks.push(stats_samtools_task_stats);
            local_cleanup_tasks.push(stats_samtools_err_task_stats);
            {
                let mut guard = stats_samtools_child_stats.lock().await;
                parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Lines, &config)
                    .await?
            }
        }
        None => return Err(anyhow!("consensus_bam_stats_stream is not available")),
    };
    let samtools_stats_out = parse_samtools_stats(stats_samtools_out_stream_stats).await?;
    let samtools_summary = samtools_stats_out.summary;

    let depth_samtools_config = SamtoolsConfig::new(SamtoolsSubcommand::Depth)
        .flag("-aa")
        .option("-d", "0")
        .positional("-");
    let depth_samtools_args = generate_cli(SAMTOOLS_TAG, &config, Some(&depth_samtools_config))?;

    let depth_samtools_out_stream = match consensus_bam_depth_stream {
        Some(stream) => {
            let (depth_samtools_child, depth_samtools_task, depth_samtools_err_task) =
                stream_to_cmd(
                    config.clone(),
                    stream.into_inner(),
                    SAMTOOLS_TAG,
                    depth_samtools_args,
                    StreamDataType::JustBytes,
                    config.args.verbose,
                    None,
                )
                    .await?;
            local_cleanup_tasks.push(depth_samtools_task);
            local_cleanup_tasks.push(depth_samtools_err_task);
            {
                let mut guard = depth_samtools_child.lock().await;
                parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Lines, &config)
                    .await?
            }
        }
        None => {
            return Err(anyhow!("consensus_bam_depth_stream is not available"));
        }
    };

    let depth_map = parse_samtools_depth(depth_samtools_out_stream).await?;
    if depth_map.is_empty() {
        return Err(anyhow!("No depth data found"));
    }
    let first_chr = depth_map
        .keys()
        .next()
        .ok_or_else(|| anyhow!("No chromosomes found"))?
        .clone();
    let first_chr_depth_map = depth_map.get(&first_chr).unwrap();
    let depths: Vec<u32> = first_chr_depth_map.values().copied().collect();
    let samtools_depth_stats = compute_depth_stats(&depths)?;

    let depth_plot_path = file_path_manipulator(
        &PathBuf::from(no_ext_sample_base),
        Some(out_dir),
        None,
        Some("depth.png"),
        "_",
    );

    plot_depths(&first_chr_depth_map, no_ext_sample_base, &depth_plot_path)?;

    let seqkit_stats = parse_seqkit_stats(no_host_seqkit_out_stream_stats).await?;

    let ercc_stats = if let Technology::Illumina = technology {
        match ercc_stats_task {
            Some(task) => match task.await? {
                Ok(stats) => stats,
                Err(e) => return Err(anyhow!("Failed to parse ERCC stats: {}", e)),
            },
            None => {
                return Err(anyhow!(
                    "ERCC stats task not initialized for Illumina technology"
                ))
            }
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
        let bcftools_stats_args = vec!["stats".to_string(), "-".to_string()];
        let (bcftools_stats_child, bcftools_stats_stream_task, bcftools_stats_err_task) =
            stream_to_cmd(
                config.clone(),
                stream.into_inner(),
                BCFTOOLS_TAG,
                bcftools_stats_args,
                StreamDataType::JustBytes,
                config.args.verbose,
                None,
            )
                .await?;

        local_cleanup_tasks.push(bcftools_stats_stream_task);
        local_cleanup_tasks.push(bcftools_stats_err_task);

        // Parse the stdout of bcftools stats as lines
        let bcftools_stats_out_rx = {
            let mut guard = bcftools_stats_child.lock().await;
            parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Lines, &config).await?
        };

        count_variants_from_bcftools_stats(bcftools_stats_out_rx).await?
    } else {
        (0, 0, 0)
    };

    let n_actg = allele_counts
        .iter()
        .filter(|&(k, _)| "ACTGU".contains(*k))
        .map(|(_, &v)| v)
        .sum::<u64>();
    let n_missing = allele_counts.get(&'N').copied().unwrap_or(0);
    let n_gap = allele_counts.get(&'-').copied().unwrap_or(0);
    let n_ambiguous = allele_counts
        .iter()
        .filter(|&(k, _)| !"ACTGUN-".contains(*k))
        .map(|(_, &v)| v)
        .sum::<u64>();

    let coverage_breadth = if !depths.is_empty() {
        depths.iter().filter(|&&d| d > 0).count() as f64 / depths.len() as f64
    } else {
        0.0
    };
    let max_aligned_length = depths.len();
    let total_length = depths.len();
    let (coverage_bin_size, coverage) = if !depths.is_empty() {
        compute_coverage_bins(&depths, 500)
    } else {
        (0.0, Vec::new())
    };

    let stats = Stats {
        sample_name: no_ext_sample_base.to_string(),
        depth_avg: samtools_depth_stats
            .get("depth_avg")
            .copied()
            .unwrap_or(0.0),
        depth_q25: samtools_depth_stats
            .get("depth_q.25")
            .copied()
            .unwrap_or(0.0),
        depth_q50: samtools_depth_stats
            .get("depth_q.5")
            .copied()
            .unwrap_or(0.0),
        depth_q75: samtools_depth_stats
            .get("depth_q.75")
            .copied()
            .unwrap_or(0.0),
        depth_frac_above_10x: samtools_depth_stats
            .get("depth_frac_above_10x")
            .copied()
            .unwrap_or(0.0),
        depth_frac_above_25x: samtools_depth_stats
            .get("depth_frac_above_25x")
            .copied()
            .unwrap_or(0.0),
        depth_frac_above_50x: samtools_depth_stats
            .get("depth_frac_above_50x")
            .copied()
            .unwrap_or(0.0),
        depth_frac_above_100x: samtools_depth_stats
            .get("depth_frac_above_100x")
            .copied()
            .unwrap_or(0.0),
        allele_counts,
        total_reads: seqkit_stats
            .get("num_seqs")
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0),
        mapped_reads: samtools_summary
            .get("reads mapped")
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0),
        mapped_paired: samtools_summary
            .get("reads mapped and paired")
            .and_then(|s| s.parse::<u64>().ok()),
        paired_inward: samtools_summary
            .get("inward oriented pairs")
            .and_then(|s| s.parse::<u64>().ok())
            .map(|v| v * 2),
        paired_outward: samtools_summary
            .get("outward oriented pairs")
            .and_then(|s| s.parse::<u64>().ok())
            .map(|v| v * 2),
        paired_other_orientation: samtools_summary
            .get("pairs with other orientation")
            .and_then(|s| s.parse::<u64>().ok())
            .map(|v| v * 2),
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

    try_join_all(local_cleanup_tasks)
        .await?
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;

    Ok(())
}

/// Evaluates the generated assembly using QUAST.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `target_ref_fasta_path`: path to the target reference FASTA file
/// * `align_bam_path`: path to the BAM file of reads aligned to the target
/// * `consensus_file_path`: path to the generated consensus FASTA file
///
/// # Returns
///
/// Result containing a vector of cleanup tasks
async fn evaluate_assembly(
    config: Arc<RunConfig>,
    target_ref_fasta_path: PathBuf,
    align_bam_path: PathBuf,
    consensus_file_path: PathBuf,
) -> Result<Vec<JoinHandle<Result<(), anyhow::Error>>>, PipelineError> {
    let mut cleanup_tasks = vec![];
    let quast_config = QuastConfig {
        ref_fasta: target_ref_fasta_path.to_string_lossy().into_owned(),
        ref_bam: align_bam_path.to_string_lossy().into_owned(),
        assembly_fasta: consensus_file_path.to_string_lossy().into_owned(),
    };

    let assembly_eval_quast_args =
        generate_cli(QUAST_TAG, &*config, Some(&quast_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: QUAST_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (mut assembly_eval_quast_child, assembly_eval_quast_err_task) = spawn_cmd(
        config.clone(),
        QUAST_TAG,
        assembly_eval_quast_args,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
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

    Ok(cleanup_tasks)
}

/// The main entry point for the consensus genome pipeline.
///
/// # Arguments
///
/// * `config`: the run configuration
///
/// # Returns
///
/// Result<()>: success or error
pub async fn run(config: Arc<RunConfig>) -> Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let ram_temp_dir = config.ram_temp_dir.clone();
    let out_dir = config.out_dir.clone();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let mut temp_dirs: Vec<TempDir> = Vec::new();
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
    let mut consensus_stats_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut call_bcftools_stats_stream: Option<ReceiverStream<ParseOutput>> = None;

    info!("Starting consensus genome pipeline.");

    // External tools check
    check_versions(
        vec![
            SAMTOOLS_TAG,
            MINIMAP2_TAG,
            SAMTOOLS_TAG,
            KRAKEN2_TAG,
            BCFTOOLS_TAG,
            MAFFT_TAG,
            SEQKIT_TAG,
            QUAST_TAG,
        ],
        &out_dir.clone(),
        &config,
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let (file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base) =
        validate_file_inputs(&config, &cwd, true).await?;

    let technology = config.args.technology.clone();

    let (
        _host_ref_fasta_path,
        host_ref_index_path,
        host_ref_temp,
        host_index_temp,
        host_temp_dir,
        host_ref_tasks,
    ) = minimap2_index_prep(
        &config,
        &ram_temp_dir,
        config.args.host_sequence.clone(),
        config.args.host_mm2_index.clone(),
        "host",
    )
        .await?;
    try_join_all(host_ref_tasks).await?;
    if let Some(temp) = host_ref_temp {
        temp_files.push(temp);
    }
    if let Some(index_temp) = host_index_temp {
        temp_files.push(index_temp);
    }
    if let Some(index_temp) = host_temp_dir {
        temp_dirs.push(index_temp);
    }

    let (
        target_ref_fasta_path_inner,
        target_ref_index_path,
        target_ref_temp,
        target_index_temp,
        target_temp_dir,
        target_ref_tasks,
    ) = minimap2_index_prep(
        &config,
        &ram_temp_dir,
        config.args.target_sequence.clone(),
        config.args.target_index.clone(),
        "target",
    )
        .await?;
    try_join_all(target_ref_tasks).await?;
    if let Some(temp) = target_ref_temp {
        temp_files.push(temp);
    }
    if let Some(index_temp) = target_index_temp {
        temp_files.push(index_temp);
    }
    if let Some(index_temp) = target_temp_dir {
        temp_dirs.push(index_temp);
    }

    let target_fasta = target_ref_fasta_path_inner.ok_or_else(|| {
        PipelineError::InvalidConfig(
            "Target FASTA required for filter_with_kraken, call_variants, and evaluate_assembly"
                .to_string(),
        )
    })?;
    target_ref_fasta_path = Some(target_fasta.clone());

    // Input Validation
    let (val_out_stream, validate_cleanup_tasks, validate_cleanup_receivers, val_count_task) =
        validate_input(
            config.clone(),
            file1_path,
            file2_path,
            no_ext_sample_base_buf.clone(),
            &out_dir,
        )
            .await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // Host Removal
    let no_host_file_path = file_path_manipulator(
        &no_ext_sample_base_buf,
        Some(&out_dir),
        None,
        Some("no_host.fq.gz"),
        "_",
    );

    let (
        no_host_output_stream,
        no_host_seqkit_out_stream_stats,
        no_host_cleanup_tasks,
        no_host_cleanup_receivers,
    ) = align_to_host(
        config.clone(),
        val_out_stream,
        host_ref_index_path,
        no_host_file_path,
    )
        .await?;
    cleanup_tasks.extend(no_host_cleanup_tasks);
    cleanup_receivers.extend(no_host_cleanup_receivers);

    // Counting stats for the host-removed reads
    let stats_seqkit_config_stats =
        SeqkitConfig::new(SeqkitSubcommand::Stats)
            .positional("-");

    let stats_seqkit_args_stats =
        generate_cli(
            SEQKIT_TAG,
            &config,
            Some(&stats_seqkit_config_stats),
        )?;

    let no_host_seqkit_out_stream_stats = no_host_seqkit_out_stream_stats.into_inner();
    let (no_host_seqkit_child_stats, no_host_seqkit_task_stats, no_host_seqkit_err_task_stats) =
        stream_to_cmd(
            config.clone(),
            no_host_seqkit_out_stream_stats,
            SEQKIT_TAG,
            stats_seqkit_args_stats,
            StreamDataType::JustBytes,
            config.args.verbose,
            None,
        )
            .await?;
    let no_host_seqkit_out_stream_stats = {
        let mut guard = no_host_seqkit_child_stats.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Lines, &config).await?
    };
    stats_tasks.push(no_host_seqkit_task_stats);
    cleanup_tasks.push(no_host_seqkit_err_task_stats);

    // Split by Technology
    match technology {
        Technology::Illumina => {
            info!("Technology: Illumina");
            let (
                _ercc_fasta_path,
                ercc_index_path,
                ercc_ref_temp,
                ercc_index_temp,
                ercc_temp_dir,
                ercc_ref_tasks,
            ) = minimap2_index_prep(
                &config,
                &ram_temp_dir,
                config.args.ercc_sequence.clone(), // Use ercc_sequences as sequence path
                config.args.ercc_index.clone(),
                "ercc",
            )
                .await?;
            try_join_all(ercc_ref_tasks).await?;
            if let Some(temp) = ercc_ref_temp {
                temp_files.push(temp);
            }
            if let Some(index_temp) = ercc_index_temp {
                temp_files.push(index_temp);
            }
            if let Some(index_temp) = ercc_temp_dir {
                temp_dirs.push(index_temp);
            }

            let (no_host_branches, no_host_router_done) = fanout_to_channels(
                no_host_output_stream,
                2,
                "consensus_no_host_split",
                &config,
                StreamDataType::JustBytes,
            ).await?;

            cleanup_receivers.push(no_host_router_done);

            let mut no_host_branches = no_host_branches.into_iter();

            let filter_reads_input_stream =
                ReceiverStream::new(no_host_branches.next().ok_or(PipelineError::EmptyStream)?);

            let ercc_input_stream =
                ReceiverStream::new(no_host_branches.next().ok_or(PipelineError::EmptyStream)?);


            // ERCC
            let (
                ercc_stats_out_task,
                mut ercc_cleanup_tasks,
                mut ercc_cleanup_receivers,
            ) = process_ercc(
                config.clone(),
                ercc_input_stream,
                ercc_index_path,
                &out_dir,
                &no_ext_sample_base,
            )
                .await?;

            ercc_stats_task = ercc_stats_out_task;
            cleanup_tasks.append(&mut ercc_cleanup_tasks);
            cleanup_receivers.append(&mut ercc_cleanup_receivers);

            // Filter Reads
            let (
                filter_reads_out_stream,
                filter_reads_cleanup_tasks,
                filter_reads_cleanup_receivers,
                filter_reads_temp_files,
                filter_reads_temp_dirs,
            ) = filter_with_kraken(
                config.clone(),
                filter_reads_input_stream,
                target_ref_index_path.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
                config
                    .args
                    .target_taxid
                    .as_ref()
                    .expect("target_taxid must be set"),
            )
                .await?;

            cleanup_tasks.extend(filter_reads_cleanup_tasks);
            cleanup_receivers.extend(filter_reads_cleanup_receivers);
            temp_files.extend(filter_reads_temp_files);
            temp_dirs.extend(filter_reads_temp_dirs);

            // Align Reads to Target
            let (
                sam_output_stream,
                align_cleanup_tasks,
                align_cleanup_receivers,
                align_quast_tasks,
                align_bam_path_inner,
            ) = align_to_target(
                config.clone(),
                filter_reads_out_stream,
                target_ref_index_path,
                &out_dir,
                &no_ext_sample_base_buf,
            )
                .await?;
            cleanup_tasks.extend(align_cleanup_tasks);
            cleanup_receivers.extend(align_cleanup_receivers);
            quast_write_tasks.extend(align_quast_tasks);
            align_bam_path = Some(align_bam_path_inner);

            // Split SAM streams for bypass, stats, etc
            let (align_sam_streams, align_sam_router_handle) = fanout_to_channels(
                sam_output_stream,
                4,
                "pipeline_aligned_sam_split",
                &config,
                StreamDataType::JustBytes,
            )
                .await?;
            cleanup_receivers.push(align_sam_router_handle);
            let mut align_streams_iter = align_sam_streams.into_iter();
            let align_sam_output_stream = align_streams_iter
                .next()
                .ok_or(PipelineError::EmptyStream)?;
            let align_sam_call_stream = align_streams_iter
                .next()
                .ok_or(PipelineError::EmptyStream)?;
            align_sam_stats_stream = Some(ReceiverStream::new(
                align_streams_iter
                    .next()
                    .ok_or(PipelineError::EmptyStream)?,
            ));
            align_sam_depth_stream = Some(ReceiverStream::new(
                align_streams_iter
                    .next()
                    .ok_or(PipelineError::EmptyStream)?,
            ));

            let align_sam_output_stream = ReceiverStream::new(align_sam_output_stream);
            let align_sam_call_stream = ReceiverStream::new(align_sam_call_stream);

            // Make Consensus
            let (
                consensus_realign_stream,
                consensus_stats_stream_rx,
                consensus_file_path_x,
                consensus_cleanup_tasks,
                consensus_cleanup_receivers,
                consensus_quast_tasks,
            ) = generate_consensus(
                config.clone(),
                align_sam_output_stream,
                &out_dir,
                &no_ext_sample_base_buf,
            )
                .await?;
            consensus_file_path = Some(consensus_file_path_x);
            cleanup_tasks.extend(consensus_cleanup_tasks);
            cleanup_receivers.extend(consensus_cleanup_receivers);
            quast_write_tasks.extend(consensus_quast_tasks);
            consensus_stats_stream = Some(consensus_stats_stream_rx);

            // Call Variants
            let (call_bcftools_stats_stream_out, _, call_cleanup_tasks, call_cleanup_receivers) =
                call_variants(
                    config.clone(),
                    align_sam_call_stream,
                    target_fasta.clone(),
                    &out_dir,
                    &no_ext_sample_base_buf,
                )
                    .await?;
            cleanup_tasks.extend(call_cleanup_tasks);
            cleanup_receivers.extend(call_cleanup_receivers);
            call_bcftools_stats_stream = Some(call_bcftools_stats_stream_out);

            // Realign Consensus to Ref
            let realign_cleanup_tasks = realign_consensus_to_ref(
                config.clone(),
                consensus_realign_stream.into_inner(),
                target_fasta.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
            )
                .await?;
            cleanup_tasks.extend(realign_cleanup_tasks);
        }
        Technology::ONT => {
            error!("Technology: ONT not ready");
        }
    }

    // Join stats tasks
    let results = try_join_all(stats_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }

    // Calculate Statistics
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
    )
        .await?;

    // Assembly Evaluation
    let results = try_join_all(quast_write_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }

    if let (Some(target_ref_fasta_path), Some(align_bam_path), Some(consensus_file_path)) =
        (target_ref_fasta_path, align_bam_path, consensus_file_path)
    {
        evaluate_assembly(
            config,
            target_ref_fasta_path,
            align_bam_path,
            consensus_file_path,
        )
            .await?;
    }

    // Cleanup
    let results = try_join_all(cleanup_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }
    for receiver in cleanup_receivers {
        receiver
            .await
            .map_err(|e| PipelineError::Other(e.into()))?
            .map_err(|e| PipelineError::Other(e))?;
    }
    let _input_read_stats = join_with_error_handling(val_count_task).await?;
    drop(temp_files);
    drop(temp_dirs);

    info!("Finished consensus genome pipeline.");
    Ok(())
}
