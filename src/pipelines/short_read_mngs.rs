
use std::path::PathBuf;
use std::sync::Arc;
use std::io::{BufReader, BufRead};
use std::hash::Hasher;
use std::cmp::{Ord, Ordering, PartialOrd, Eq, PartialEq};
use std::collections::{HashMap, BinaryHeap};
use std::cmp::Reverse;

use anyhow::{anyhow, Result};
use log::{self, LevelFilter, debug, info, error, warn};
use futures::future::try_join_all;
use rand::prelude::*;
use rand_core::{RngCore, OsRng};
use tokio::fs;
use tokio::time::{sleep, Duration, Instant};
use tokio::sync::oneshot;
use tokio::sync::mpsc;
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;
use tokio::sync::Notify;
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufWriter, BufReader as TokioBufReader};
use tokio::fs::{OpenOptions as TokioOpenOptions};
use tempfile::NamedTempFile;
use tempfile::TempDir;
use serde_json::Value;
use twox_hash::XxHash64;

use crate::config::defs::{PipelineError, RunConfig, StreamDataType, ReadStats, MINIMAP2_TAG, BOWTIE2_TAG, SAMTOOLS_TAG, FASTP_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, SEQKIT_TAG, QUAST_TAG, HISAT2_TAG, SamtoolsSubcommand, KALLISTO_TAG, KallistoSubcommand, STAR_TAG, SamtoolsStats, CZID_DEDUP_TAG};
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file, available_space_for_path};
use crate::utils::fastx::{raw_read_count, read_fastq, stream_record_counter, SequenceRecord, compare_read_ids, parse_header};
use crate::utils::streams::{t_junction, ParseOutput, join_with_error_handling, stream_to_cmd, parse_child_output, ChildStream, ParseMode, stream_to_file, read_child_output_to_vec, spawn_cmd, parse_fastq, ChannelReader, write_to_fifo, deinterleave_fastq_stream, create_fifo, interleave_fastq_streams, ToBytes};
use crate::utils::command::bowtie2::{Bowtie2Config, bowtie2_index_prep};
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::kallisto::KallistoConfig;
use crate::utils::command::hisat2::{Hisat2Config, hisat2_index_prep};
use crate::utils::streams::{deinterleave_fastq_stream_to_fifos};
use crate::utils::command::star::{StarConfig, star_index_prep};
use crate::utils::command::minimap2::{Minimap2ArgGenerator, Minimap2Config, minimap2_index_prep};
use crate::utils::stats::parse_samtools_stats;
use crate::utils::plotting::plot_insert_sizes;
use crate::utils::command::czid_dedup::CzidDedupConfig;
use crate::utils::paf::PafRecord;

#[derive(Debug)]
pub struct KallistoResults {
    pub ercc_counts: Vec<(String, f64)>, // (target_id, est_counts)
    pub transcript_to_gene: Vec<(String, String)>, // (transcript_id, gene_id)
}

// SampleItem for reservoir sampling
#[derive(Debug, Clone)]
struct SampleItem {
    key: f64, // For A-Res sampling
    records: Vec<SequenceRecord>, // 1 for SE, 2 for PE
    weight: u64, // Duplicate count
    hash_key: Option<u64>, // For disk mode
}

impl PartialEq for SampleItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}
impl Eq for SampleItem {}
impl Ord for SampleItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.partial_cmp(&other.key).unwrap_or(Ordering::Equal)
    }
}
impl PartialOrd for SampleItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

//Memory estimates for dedup_and_subsample, for switching between RAM and disk based
const EST_BYTES_PER_UNIQUE: u64 = 1024;  // ~0.5–1KB for SequenceRecord pair + overhead
const RAM_BUFFER_GB: u64 = 10;  // Min available RAM to stay in-memory

/// Called read_fastq in single or paired FASTQ's and streams interleaved output
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `file1_path` - Path to R1 or single ended FASTQ.
/// * `file2_path` - Optional path to R2.
///  * `sample_base_buf` - PathBuf of sample basename.
/// * `out_dir` - Base dir for output files.
///
/// # Returns
/// ParseOutput validated interleaved FASTQ stream, vecs of cleaup takss and recoevers, handle for raw and validated read counts.
async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> anyhow::Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>, JoinHandle<anyhow::Result<u64, anyhow::Error>>, JoinHandle<anyhow::Result<ReadStats, anyhow::Error>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();


    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated.fq.gz"),
        "_",
    );

    let raw_count_task = raw_read_count(file1_path.clone(), file2_path.clone());


    let (rx, val_count_task) = read_fastq(
        file1_path,
        file2_path,
        Some(config.args.technology.clone()),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        config.base_buffer_size,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);

    // Split the byte stream for fastp and write
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(val_done_rx);

    let mut val_streams_iter = val_streams.into_iter();
    let val_ercc_bowtie2_filter_out_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let val_file_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let val_file_write_task = write_byte_stream_to_file(
        &validated_interleaved_file_path,
        ReceiverStream::new(val_file_stream),
        Some(config.base_buffer_size),
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(val_file_write_task);

    Ok((ReceiverStream::new(val_ercc_bowtie2_filter_out_stream),  cleanup_tasks, cleanup_receivers, raw_count_task, val_count_task))
}

/// bowtie2 filter function where the passing stream contains the unmapped reads
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
/// * `bt2_index_path` - Path to Bowtie2 index.
/// * `paired` - Whether the input is paired-end.
/// * `bowtie2_options` - Additional Bowtie2 options as a HashMap (e.g., HashMap::from([("--very-sensitive-local".to_string(), None)])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - unmapped FASTQ stream.
/// - Optional receiver for the total mapped count (u64) if `count_mapped` is true.
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.

async fn bowtie2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    bowtie2_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
    compute_insert_stats: bool,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, Option<oneshot::Receiver<SamtoolsStats>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let paired_compute_insert_stats = paired && compute_insert_stats;

    // BT2
    let bt2_config_view = Bowtie2Config {
        bt2_index_path: bt2_index_path.clone(),
        option_fields: bowtie2_options,
    };

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bt2_child, bt2_stream_task, bt2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        BOWTIE2_TAG,
        bt2_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bt2_stream_task);
    cleanup_tasks.push(bt2_err_task);

    let bt2_out_stream = {
        let mut guard = bt2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BOWTIE2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted for fastq extraction
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out_stream,
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

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 } + if paired_compute_insert_stats { 1 } else { 0 };
    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = t_junction(
        bam_rx_stream,
        num_tees,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "bowtie2_bam_split".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let bam_file_stream = if output_bam_path.is_some() {
        Some(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?)
    } else {
        None
    };
    let insert_stats_stream = if paired_compute_insert_stats {
        Some(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?)
    } else {
        None
    };

    // Handle BAM output if requested
    if let (Some(bam_path), Some(stream)) = (output_bam_path, bam_file_stream) {
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count mapped reads
    let (count_tx, count_rx) = oneshot::channel();
    let count_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(count_stream);
        let mut line_count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(line) if !line.is_empty() => line_count += 1,
                _ => continue,
            }
        }
        let count = line_count as u64 / if paired { 8 } else { 4 };
        count_tx.send(count).map_err(|_| anyhow!("Failed to send count"))?;
        Ok(())
    });
    cleanup_tasks.push(count_task);

    // Samtools fastq to extract unmapped reads
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Compute insert size metrics if paired and compute_insert_stats is true
    let insert_stats_rx = if paired_compute_insert_stats {
        let insert_stats_stream = insert_stats_stream.unwrap();

        // Split insert_stats_stream for counting and processing
        let (mut stats_streams, stats_done_rx) = t_junction(
            ReceiverStream::new(insert_stats_stream),
            2,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes,
            "insert_stats_split".to_string(),
            None,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?;
        cleanup_receivers.push(stats_done_rx);

        let mut stats_streams_iter = stats_streams.into_iter();
        let count_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let process_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

        // Count mapped reads in BAM
        let (count_tx, count_rx) = oneshot::channel();
        let count_task = tokio::spawn(async move {
            let mut stream = ReceiverStream::new(count_stream);
            let mut line_count = 0;
            while let Some(item) = stream.next().await {
                match item {
                    ParseOutput::Bytes(line) if !line.is_empty() => line_count += 1,
                    _ => continue,
                }
            }
            count_tx.send(line_count).map_err(|_| anyhow!("Failed to send BAM line count"))?;
            Ok(())
        });
        cleanup_tasks.push(count_task);

        let bam_line_count = count_rx
            .await
            .map_err(|e| PipelineError::Other(anyhow!("BAM count receiver failed: {}", e)))?;

        if bam_line_count == 0 {
            warn!("No mapped reads in BAM for insert size stats; skipping samtools stats");
            let (stats_tx, stats_rx) = oneshot::channel();
            stats_tx
                .send(SamtoolsStats {
                    summary: HashMap::new(),
                    insert_sizes: Vec::new(),
                })
                .map_err(|_| PipelineError::Other(anyhow!("Failed to send empty stats")))?;
            Some(stats_rx)
        } else {
            // Coordinate-sort BAM for samtools stats
            let samtools_coord_sort_config = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Sort,
                subcommand_fields: HashMap::from([
                    ("-u".to_string(), None),
                    ("-O".to_string(), Some("bam".to_string())),
                    ("-".to_string(), None),
                ]),
            };
            let samtools_coord_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_coord_sort_config))
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;

            let (mut coord_sort_child, coord_sort_task, coord_sort_err_task) = stream_to_cmd(
                config.clone(),
                process_stream,
                SAMTOOLS_TAG,
                samtools_coord_sort_args,
                StreamDataType::JustBytes,
                config.args.verbose,
            )
                .await
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;
            cleanup_tasks.push(coord_sort_task);
            cleanup_tasks.push(coord_sort_err_task);

            let coord_sort_out_stream = {
                let mut guard = coord_sort_child.lock().await;
                parse_child_output(
                    &mut guard,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    config.base_buffer_size,
                )
                    .await
                    .map_err(|e| PipelineError::ToolExecution {
                        tool: SAMTOOLS_TAG.to_string(),
                        error: e.to_string(),
                    })?
            };

            // Run samtools stats
            let samtools_stats_config = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Stats,
                subcommand_fields: HashMap::new(),
            };
            let samtools_stats_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_stats_config))
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;

            let (mut stats_child, stats_stream_task, stats_err_task) = stream_to_cmd(
                config.clone(),
                coord_sort_out_stream,
                SAMTOOLS_TAG,
                samtools_stats_args,
                StreamDataType::JustBytes,
                config.args.verbose,
            )
                .await
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;
            cleanup_tasks.push(stats_stream_task);
            cleanup_tasks.push(stats_err_task);

            let stats_out_stream = {
                let mut guard = stats_child.lock().await;
                parse_child_output(
                    &mut guard,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    config.base_buffer_size,
                )
                    .await
                    .map_err(|e| PipelineError::ToolExecution {
                        tool: SAMTOOLS_TAG.to_string(),
                        error: e.to_string(),
                    })?
            };

            // Parse samtools stats
            let (stats_tx, stats_rx) = oneshot::channel();
            let stats_task = tokio::spawn(async move {
                let stats = parse_samtools_stats(stats_out_stream)
                    .await
                    .map_err(|e| anyhow!("Failed to parse samtools stats: {}", e))?;
                stats_tx.send(stats).map_err(|_| anyhow!("Failed to send samtools stats"))?;
                Ok(())
            });
            cleanup_tasks.push(stats_task);

            Some(stats_rx)
        }
    } else {
        None
    };

    Ok((
        ReceiverStream::new(unmapped_fastq_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        insert_stats_rx,
    ))
}


/// QC's input stream using FASTP
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns
///
async fn fastp_qc(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<(ReceiverStream<ParseOutput>,  Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, oneshot::Receiver<Result<u64, anyhow::Error>>), PipelineError>{
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let qc_fastp_config_view = FastpConfig {

        //These default QC thresholds are loosely based on the pre-2022 CZI pipeline using PriceSeq & LZW
        command_fields: HashMap::from([
            ("--dont_eval_duplication".to_string(), None),
            ("--length_required".to_string(), Some("35".to_string())),
            ("--qualified_quality_phred".to_string(), Some("17".to_string())),
            ("--unqualified_percent_limit".to_string(), Some("15".to_string())),
            ("--n_base_limit".to_string(), Some("15".to_string())),
            ("--complexity_threshold".to_string(), Some("60".to_string())),
        ]),
    };

    let qc_fastp_args = generate_cli(FASTP_TAG, &config, Some(&qc_fastp_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut qc_fastp_child, qc_fastp_stream_task, qc_fastp_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        FASTP_TAG,
        qc_fastp_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(qc_fastp_stream_task);
    cleanup_tasks.push(qc_fastp_err_task);

    let qc_fastp_out_stream = {
        let mut guard = qc_fastp_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Tee the byte stream for counting
    let tee_count_input = ReceiverStream::new(qc_fastp_out_stream);
    let (mut tee_count_streams, tee_count_done_rx) = t_junction(
        tee_count_input,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "fastp_count".to_string(),
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "t_junction".to_string(),
            error: e.to_string(),
        })?;
    cleanup_receivers.push(tee_count_done_rx);

    let mut tee_count_streams_iter = tee_count_streams.into_iter();
    let qc_out_stream = tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let qc_count_stream = tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (count_tx, count_result_rx) = oneshot::channel::<Result<u64, anyhow::Error>>();

    let count_task = tokio::spawn(async move {
        match stream_record_counter(qc_count_stream, false).await {
            Ok(count) => {
                debug!("Fastp output reads: {}", count);
                let _ = count_tx.send(Ok(count)); // Send the count result
                Ok(())
            }
            Err(e) => {
                let _ = count_tx.send(Err(e));
                Err(anyhow!("Failed to count fastp output reads"))
            },
        }
    });
    cleanup_tasks.push(count_task);

    Ok((ReceiverStream::new(qc_out_stream), cleanup_tasks, cleanup_receivers, count_result_rx))
}

/// Runs Kallisto quant
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns
async fn kallisto_quant(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    out_dir: PathBuf,
    paired: bool,
    sample_base_buf: PathBuf,
) -> Result<(
    oneshot::Sender<KallistoResults>,
    oneshot::Receiver<KallistoResults>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    JoinHandle<Result<(), anyhow::Error>>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Convert sample_base_buf to &str
    let sample_base = sample_base_buf
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| PipelineError::IOError(format!("Invalid sample base path: {}", sample_base_buf.display())))?
        .to_string();

    // Create Kallisto output dir (idempotent, fast on NVMe)
    let kallisto_out_dir = out_dir.join("kallisto");
    fs::create_dir_all(&kallisto_out_dir).await
        .map_err(|e| PipelineError::IOError(format!("Failed to create Kallisto output directory: {}", e)))?;
    debug!("Kallisto output dir created/confirmed: {}", kallisto_out_dir.display());

    let kallisto_index = config.args.kallisto_index.clone()
        .ok_or(PipelineError::MissingArgument("kallisto_index required".to_string()))?;

    // Convert raw byte stream to FASTQ
    let byte_rx = input_stream.into_inner();
    let byte_reader = ChannelReader::new(byte_rx);
    let fastq_rx = parse_fastq(byte_reader, config.base_buffer_size).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "parse_fastq".to_string(),
            error: e.to_string(),
        })?;
    let fastq_stream = ReceiverStream::new(fastq_rx);

    debug!("Deinterleaving FASTQ stream for sample: {}", sample_base);
    let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        fastq_stream,
        &sample_base,
        paired,
    ).await?;

    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle {
        cleanup_tasks.push(r2_handle);
    }

    let kallisto_config = if paired {
        KallistoConfig {
            subcommand: KallistoSubcommand::Quant,
            subcommand_fields: HashMap::from([
                ("R1".to_string(), Some(r1_fifo.to_string_lossy().to_string())),
                ("R2".to_string(), Some(r2_fifo.to_string_lossy().to_string())),
            ]),
            output_dir: kallisto_out_dir.clone(),
            reproducible: false,
        }
    } else {
        KallistoConfig {
            subcommand: KallistoSubcommand::Quant,
            subcommand_fields: HashMap::from([
                ("--single".to_string(), None),
                ("-l".to_string(), Some("200".to_string())),
                ("-s".to_string(), Some("20".to_string())),
                ("R1".to_string(), Some(r1_fifo.to_string_lossy().to_string())),
            ]),
            output_dir: kallisto_out_dir.clone(),
            reproducible: false,
        }
    };

    let kallisto_args = generate_cli(KALLISTO_TAG, &config, Some(&kallisto_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    debug!("Spawning Kallisto with args: {:?}", kallisto_args);

    let (mut kallisto_child, kallisto_err_task) = spawn_cmd(
        config.clone(),
        KALLISTO_TAG,
        kallisto_args,
        config.args.verbose,
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(kallisto_err_task);

    // Create Kallisto exit task
    let r1_fifo_clone = r1_fifo.clone();
    let r2_fifo_clone = r2_fifo.clone();
    let kallisto_exit_task = tokio::spawn(async move {
        let status = kallisto_child.wait().await
            .map_err(|e| anyhow!("Kallisto wait failed: {}", e))?;
        if !status.success() {
            return Err(anyhow!("Kallisto exited with code: {:?}", status.code()));
        }
        // Clean up FIFOs
        tokio::fs::remove_file(&r1_fifo_clone).await
            .map_err(|e| anyhow!("Failed to remove R1 FIFO {}: {}", r1_fifo_clone.display(), e))?;
        if paired {
            tokio::fs::remove_file(&r2_fifo_clone).await
                .map_err(|e| anyhow!("Failed to remove R2 FIFO {}: {}", r2_fifo_clone.display(), e))?;
        }
        Ok(())
    });

    // Create oneshot channel for results
    let (ercc_tx, ercc_rx) = oneshot::channel::<KallistoResults>();

    Ok((ercc_tx, ercc_rx, cleanup_tasks, cleanup_receivers, kallisto_exit_task))
}

/// Parses results of kallisto
///
/// # Arguments
///
/// * 'kallisto_out_dir': PathBuf to output dir of kallisto
/// *  'ercc_tx': oneshot::Sender<KallistoResults>,
///
/// # Returns
async fn kallisto_results(
    kallisto_out_dir: PathBuf,
    ercc_tx: oneshot::Sender<KallistoResults>,
) -> Result<JoinHandle<Result<(), anyhow::Error>>, PipelineError> {
    let abundance_path = kallisto_out_dir.join("abundance.tsv");
    let parse_task = tokio::spawn(async move {
        // Retry reading abundance.tsv (handles NVMe flush delays)
        let mut retries = 5;
        loop {
            match fs::metadata(&abundance_path).await {
                Ok(meta) if meta.len() > 0 => {
                    debug!("Found abundance.tsv: {} bytes", meta.len());
                    break;
                }
                _ => {
                    if retries == 0 {
                        return Err(anyhow!("No abundance.tsv found at {}", abundance_path.display()));
                    }
                    warn!("abundance.tsv not ready (retry {}); sleeping 2s", retries);
                    sleep(Duration::from_secs(2)).await;
                    retries -= 1;
                }
            }
        }

        let abundance_file = tokio::fs::File::open(&abundance_path).await
            .map_err(|e| anyhow!("Failed to open abundance.tsv: {}", e))?;
        let reader = tokio::io::BufReader::new(abundance_file);
        let mut lines = reader.lines();
        let mut ercc_counts = Vec::new();
        let mut transcript_to_gene = Vec::new();

        debug!("Starting to read abundance.tsv from: {}", abundance_path.display());

        if lines.next_line().await?.is_none() {
            return Err(anyhow!("Empty abundance.tsv"));
        }

        // Parse each line: target_id, length, eff_length, est_counts, tpm
        while let Some(line) = lines.next_line().await? {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                return Err(anyhow!("Invalid abundance.tsv line: {}", line));
            }
            let target_id = fields[0].to_string();
            let est_counts = fields[3].parse::<f64>()
                .map_err(|e| anyhow!("Failed to parse est_counts in line {}: {}", line, e))?;
            ercc_counts.push((target_id.clone(), est_counts));
            transcript_to_gene.push((target_id, format!("gene_{}", fields[0])));
        }

        debug!("Finished reading kallisto's abundance.tsv: {} transcripts", ercc_counts.len());
        let results = KallistoResults {
            ercc_counts,
            transcript_to_gene,
        };
        ercc_tx.send(results)
            .map_err(|_| anyhow!("Failed to send Kallisto results"))?;
        debug!("Sent Kallisto results");
        Ok(())
    });

    Ok(parse_task)
}

/// HISAT2 filter
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream.
/// * `hisat2_index_path` - Path to HISAT2 index.
/// * `paired` - Whether the input is paired-end.
/// * `hisat2_options` - Additional HISAT2 options as a HashMap (e.g., HashMap::from([("--no-spliced-alignment".to_string(), None)])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - Unmapped FASTQ stream.
/// - Optional receiver for the total mapped count (u64) if count is needed.
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.

async fn hisat2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    hisat2_index_path: PathBuf,
    paired: bool,
    hisat2_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Deinterleave and write to files in CWD
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let (r1_rx, r2_rx_opt, deint_handle) = deinterleave_fastq_stream(
        input_stream,
        paired,
        config.base_buffer_size,
    ).await.map_err(|e| PipelineError::Other(e))?;

    let r1_path = cwd.join("unmapped_r1.fq");
    let r1_write_task = write_byte_stream_to_file(
        &r1_path,
        ReceiverStream::new(r1_rx),
        Some(config.base_buffer_size),
    ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
    cleanup_tasks.push(r1_write_task);

    let r2_path_opt = if paired {
        let r2_path = cwd.join("unmapped_r2.fq");
        let r2_write_task = write_byte_stream_to_file(
            &r2_path,
            ReceiverStream::new(r2_rx_opt.unwrap()),
            Some(config.base_buffer_size),
        ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(r2_write_task);
        Some(r2_path)
    } else {
        None
    };

    deint_handle.await.map_err(|e| PipelineError::Other(e.into()))??;

    let hisat2_config = Hisat2Config {
        hisat2_index_path: hisat2_index_path.clone(),
        option_fields: hisat2_options,
        r1_path: r1_path.to_string_lossy().to_string(),
        r2_path: r2_path_opt.as_ref().map(|p| p.to_string_lossy().to_string()),
    };

    let hisat2_args = generate_cli(HISAT2_TAG, &config, Some(&hisat2_config))
        .map_err(|e| PipelineError::ToolExecution { tool: HISAT2_TAG.to_string(), error: e.to_string() })?;

    let (mut hisat2_child, hisat2_err_task) = spawn_cmd(
        config.clone(),
        HISAT2_TAG,
        hisat2_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: HISAT2_TAG.to_string(),
        error: e.to_string(),
    })?;
    
    let hisat2_out_stream = parse_child_output(
        &mut hisat2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: HISAT2_TAG.to_string(),
        error: e.to_string(),
    })?;

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted (required for paired-end fastq extraction)
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        hisat2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?
    };

    // number of streams for t_junction (unmapped + count + optional BAM)
    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = if num_tees > 1 {
        t_junction(
            bam_rx_stream,
            num_tees,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes,
            "hisat2_bam_split".to_string(),
            None,
        ).await.map_err(|_| PipelineError::StreamDataDropped)?
    } else {
        (vec![bam_rx_stream.into_inner()], oneshot::channel::<Result<(), anyhow::Error>>().1)
    };
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let mapped_flag = if paired { "-F13".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None), // Count
            (mapped_flag, None),
            ("-".to_string(), None),
        ]),
    };
    let samtools_count_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_count_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        bam_count_stream,
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (count_tx, count_rx) = oneshot::channel::<u64>();

    let count_future = tokio::spawn(async move {
        let mut count_child = count_child_arc.lock().await;
        let count_lines = read_child_output_to_vec(&mut count_child, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        let _ = count_tx.send(mapped_count);
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Unmapped stream
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Extract unmapped FASTQ
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?
    };


    Ok((ReceiverStream::new(unmapped_fastq_stream), count_rx, cleanup_tasks, cleanup_receivers))
}


/// STAR filter
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream (interleaved).
/// * `star_index_dir` - Path to STAR index directory.
/// * `paired` - Whether the input is paired-end.
/// * `star_options` - Additional STAR options as a HashMap (e.g., HashMap::from([("--outFilterMultimapNmax".to_string(), Some("5".to_string()))])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - unmapped FASTQ stream.
/// - Receiver for the total mapped count (u64).
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.
async fn star_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    star_index_path: PathBuf,
    paired: bool,
    star_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Deinterleave to FIFOs
    let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        input_stream,
        "star_filter",
        paired,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: "deinterleave".to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle {
        cleanup_tasks.push(r2_handle);
    }

    // STAR config
    let star_config = StarConfig {
        star_index_dir: star_index_prep(&star_index_path, &std::env::current_dir().unwrap())?,
        option_fields: star_options,
        r1_fifo: r1_fifo.clone(),
        r2_fifo: if paired { Some(r2_fifo.clone()) } else { None },
    };

    let star_args = generate_cli(STAR_TAG, &config, Some(&star_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: STAR_TAG.to_string(),
            error: e.to_string(),
        })?;

    // Spawn STAR (no stdin stream, uses FIFOs)
    let (mut star_child, star_err_task) = spawn_cmd(
        config.clone(),
        STAR_TAG,
        star_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: STAR_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(star_err_task);

    let star_out_stream = parse_child_output(
        &mut star_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: STAR_TAG.to_string(),
        error: e.to_string(),
    })?;

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted (required for paired-end fastq extraction)
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        star_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?
    };

    // number of streams for t_junction (unmapped + count + optional BAM)
    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = if num_tees > 1 {
        t_junction(
            bam_rx_stream,
            num_tees,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes,
            "star_bam_split".to_string(),
            None,
        ).await.map_err(|_| PipelineError::StreamDataDropped)?
    } else {
        (vec![bam_rx_stream.into_inner()], oneshot::channel::<Result<(), anyhow::Error>>().1)
    };
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let mapped_flag = if paired { "-F13".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None), // Count
            (mapped_flag, None),
            ("-".to_string(), None),
        ]),
    };
    let samtools_count_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_count_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        bam_count_stream,
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (count_tx, count_rx) = oneshot::channel::<u64>();

    let count_future = tokio::spawn(async move {
        let mut count_child = count_child_arc.lock().await;
        let count_lines = read_child_output_to_vec(&mut count_child, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        let _ = count_tx.send(mapped_count);
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Unmapped stream
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Extract unmapped FASTQ
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?
    };

    // Cleanup FIFOs
    let r1_fifo_cleanup = r1_fifo.clone();
    let r2_fifo_cleanup = r2_fifo.clone();
    let fifo_cleanup_task = tokio::spawn(async move {
        tokio::fs::remove_file(&r1_fifo_cleanup).await.ok();
        if paired {
            tokio::fs::remove_file(&r2_fifo_cleanup).await.ok();
        }
        Ok(())
    });
    cleanup_tasks.push(fifo_cleanup_task);

    Ok((ReceiverStream::new(unmapped_fastq_stream), count_rx, cleanup_tasks, cleanup_receivers))
}


/// Minimap2 filter
///
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream (paired or single-end).
/// * `paired` - Whether the input is paired-end.
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - Interleaved FASTQ stream : unmapped reads.
/// - Receiver for the total mapped read count
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
/// - Optional temporary FASTA file for host reference.
/// - Optional temporary minimap2 index file.
async fn minimap2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    output_bam_path: Option<PathBuf>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        oneshot::Receiver<u64>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Option<NamedTempFile>,
        Option<NamedTempFile>,
        Option<TempDir>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let (_host_ref_fasta_path, host_index_path, host_ref_temp, host_index_temp, host_temp_dir,  mut host_ref_tasks) =
        minimap2_index_prep(
            &config,
            &config.ram_temp_dir,
            config.args.host_sequence.clone(),
            config.args.host_index.clone(),
            "host",
        )
            .await?;
    try_join_all(host_ref_tasks).await?; // Have to wait on this to make sure index is in place ahead of running minimap2
    debug!("Host index : {}", host_index_path.display());

    let minimap2_config = Minimap2Config {
        minimap2_index_path: host_index_path,
        option_fields: HashMap::from([
            ("-ax".to_string(), Some("sr".to_string())), // Short-read preset
            ("-k".to_string(), Some("15".to_string())), // Default k-mer
            ("-w".to_string(), Some("10".to_string())), // Default window
            ("-s".to_string(), Some("30".to_string())), // Match HISAT2's score threshold
            ("-A".to_string(), Some("2".to_string())), // Match score
            ("-B".to_string(), Some("4".to_string())), // Mismatch penalty
            ("-O".to_string(), Some("4,24".to_string())), // Gap open
            ("-E".to_string(), Some("2,1".to_string())), // Gap extension
        ]),
    };
    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;


    let (mut minimap2_child, minimap2_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::IlluminaFastq, // Assuming interleaved FASTQ
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
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

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };
    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);
    let (bam_streams, bam_done_rx) = t_junction(
        bam_rx_stream,
        num_tees,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "minimap2_bam_split".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let mapped_flag = if paired { "-F12".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None),
            (mapped_flag, None),
            ("-".to_string(), None),
        ]),
    };
    let samtools_count_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_count_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        bam_count_stream,
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (count_tx, count_rx) = oneshot::channel::<u64>();
    let count_future = tokio::spawn(async move {
        let mut guard = count_child_arc.lock().await;
        let count_lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        count_tx.send(mapped_count).map_err(|_| anyhow!("Failed to send mapped count"))?;
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Extract unmapped FASTQ
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None),
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    Ok((
        ReceiverStream::new(unmapped_fastq_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        host_ref_temp,
        host_index_temp,
        host_temp_dir
    ))
}



// Deduplicate using czid-dedup with FIFOs
async fn czid_dedup_dedup(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    out_dir: PathBuf,
    sample_base_buf: PathBuf,
    prefix_length: u32,
) -> Result<(
    ReceiverStream<ParseOutput>,
    oneshot::Receiver<u64>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Check if input stream is empty or near-empty
    let (pre_tx, pre_rx) = mpsc::channel(config.base_buffer_size);
    let (t_junc_streams, t_junc_done_rx) = t_junction(
        input_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,
        "czid_dedup_pre_check".to_string(),
        None,
    ).await?;
    cleanup_receivers.push(t_junc_done_rx);
    let mut t_junc_iter = t_junc_streams.into_iter();
    let input_stream = ReceiverStream::new(t_junc_iter.next().ok_or(PipelineError::EmptyStream)?);
    let pre_count_stream = ReceiverStream::new(t_junc_iter.next().ok_or(PipelineError::EmptyStream)?);
    let pre_count_task = tokio::spawn(async move {
        let mut records = Vec::new();
        let mut pre_count_stream = pre_count_stream;
        while let Some(record) = pre_count_stream.next().await {
            records.push(format!("{:?}", record));
        }
        let count = records.len() as u64;
        Ok::<u64, anyhow::Error>(count)
    });
    let pre_count = pre_count_task.await??; // Unwrap Result<u64, anyhow::Error>
    if pre_count <= 4 {
        let (count_tx, count_rx) = oneshot::channel();
        count_tx.send(0).map_err(|_| anyhow!("Failed to send empty count"))?;
        drop(pre_tx); // Close sender to signal empty stream
        return Ok((ReceiverStream::new(pre_rx), count_rx, cleanup_tasks, cleanup_receivers));
    }
    // Create output FIFOs
    let num_files = if paired { 2 } else { 1 };
    let mut output_fifos = Vec::with_capacity(num_files);
    for i in 0..num_files {
        let suffix = if paired { if i == 0 { "r1" } else { "r2" } } else { "se" };
        let output_fifo = config.ram_temp_dir.join(format!("czid_output_{}_{}.fq", sample_base_buf.file_stem().unwrap().to_string_lossy(), suffix));
        create_fifo(&output_fifo).await?;
        output_fifos.push(output_fifo);
    }

    // Synchronization for FIFO writers
    let notify = Arc::new(Notify::new());

    // Deinterleave input stream to input FIFOs
    let sample_base = sample_base_buf.file_stem().unwrap().to_string_lossy();
    let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        input_stream,
        &sample_base,
        paired,
    ).await?;
    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(handle) = r2_write_handle {
        cleanup_tasks.push(handle);
    }

    // Input FIFOs for czid-dedup
    let input_fifos = if paired {
        vec![r1_fifo.clone(), r2_fifo.clone()]
    } else {
        vec![r1_fifo.clone()]
    };

    // Debug FIFO contents
    let r1_fifo_debug = r1_fifo.clone();
    let r2_fifo_debug = r2_fifo.clone();
    let debug_fifo_task = tokio::spawn(async move {
        for fifo in [&r1_fifo_debug, &r2_fifo_debug] {
            if let Ok(data) = fs::read(fifo).await {
                debug!("FIFO {} ok", fifo.display());
            } else {
                warn!("FIFO {} is empty or inaccessible", fifo.display());
            }
        }
        Ok(())
    });
    cleanup_tasks.push(debug_fifo_task);

    // Notify that writers are ready
    let notify_clone = notify.clone();
    let writer_ready_task = tokio::spawn(async move {
        tokio::time::sleep(Duration::from_millis(10)).await;
        notify_clone.notify_one();
        Ok(())
    });
    cleanup_tasks.push(writer_ready_task);

    // Optional cluster outputs
    let cluster_out = out_dir.join("czid_clusters.csv");
    let cluster_size_out = out_dir.join("czid_cluster_sizes.csv");

    // Config and args
    let czid_config = CzidDedupConfig {
        input_paths: input_fifos,
        output_paths: output_fifos.clone(),
        prefix_length: Some(prefix_length),
        cluster_output: Some(cluster_out),
        cluster_size_output: Some(cluster_size_out),
    };

    let czid_args = generate_cli(CZID_DEDUP_TAG, &config, Some(&czid_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: CZID_DEDUP_TAG.to_string(),
            error: e.to_string(),
        })?;

    // Wait for writers to be ready
    notify.notified().await;

    // Spawn czid-dedup
    let (mut child, err_task) = spawn_cmd(
        config.clone(),
        CZID_DEDUP_TAG,
        czid_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: CZID_DEDUP_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(err_task);

    let exit_task = tokio::spawn(async move {
        let status = child.wait().await?;
        if !status.success() {
            // Handle empty input case gracefully
            if status.code() == Some(2) {
                warn!("czid-dedup exited with code 2 (likely empty input); treating as empty output");
                return Ok(());
            }
            return Err(anyhow!("czid-dedup failed with exit code: {:?}", status.code()));
        }
        Ok(())
    });
    cleanup_tasks.push(exit_task);

    // Read from output FIFOs
    let (dedup_rx, count_rx) = if paired {
        let (r1_rx, r1_count_task) = read_fastq(
            output_fifos[0].clone(),
            None,
            None,
            0, // max_reads
            None,
            None,
            config.base_buffer_size,
        ).map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;
        let (r2_rx, r2_count_task) = read_fastq(
            output_fifos[1].clone(),
            None,
            None,
            0, // max_reads
            None,
            None,
            config.base_buffer_size,
        ).map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;
        let (inter_rx, inter_task) = interleave_fastq_streams(r1_rx, r2_rx, config.base_buffer_size).await?;
        cleanup_tasks.push(inter_task);
        // Validate pair counts and send count
        let (count_tx, count_rx) = oneshot::channel();
        let count_validation_task = tokio::spawn(async move {
            let r1_result = r1_count_task.await.map_err(|e| anyhow!("R1 count task failed: {}", e))?;
            let r2_result = r2_count_task.await.map_err(|e| anyhow!("R2 count task failed: {}", e))?;
            let (r1_stats, r2_stats) = match (r1_result, r2_result) {
                (Ok(r1_stats), Ok(r2_stats)) => (r1_stats, r2_stats),
                (Err(e), _) | (_, Err(e)) => {
                    error!("read_fastq failed for R1 or R2: {}", e);
                    count_tx.send(0).map_err(|_| anyhow!("Send failed"))?;
                    return Ok(());
                }
            };
            if r1_stats.validated != r2_stats.validated {
                return Err(anyhow!("Mismatched pair counts: R1={}, R2={}", r1_stats.validated, r2_stats.validated));
            }
            count_tx.send(r1_stats.validated / 2)
                .map_err(|_| anyhow!("Send failed"))?;
            Ok(())
        });
        cleanup_tasks.push(count_validation_task);
        (inter_rx, count_rx)
    } else {
        let (rx, count_task) = read_fastq(
            output_fifos[0].clone(),
            None,
            None,
            0, // max_reads
            None,
            None,
            config.base_buffer_size,
        ).map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;
        let (count_tx, count_rx) = oneshot::channel();
        let count_wrap = tokio::spawn(async move {
            let stats = count_task.await.map_err(|e| anyhow!("Count task failed: {}", e))?;
            match stats {
                Ok(stats) => {
                    count_tx.send(stats.validated)
                        .map_err(|_| anyhow!("Send failed"))?;
                    Ok(())
                }
                Err(e) => {
                    count_tx.send(0).map_err(|_| anyhow!("Send failed"))?;
                    Ok(())
                }
            }
        });
        cleanup_tasks.push(count_wrap);
        (rx, count_rx)
    };

    // FIFO cleanup with logging
    let fifo_cleanup = tokio::spawn(async move {
        for path in output_fifos.iter() {
            if let Err(e) = fs::remove_file(path).await {
                error!("Failed to remove FIFO {}: {}", path.display(), e);
            }
        }
        Ok(())
    });
    cleanup_tasks.push(fifo_cleanup);

    Ok((ReceiverStream::new(dedup_rx), count_rx, cleanup_tasks, cleanup_receivers))
}

async fn fastp_dedup(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    out_dir: PathBuf,
    sample_base_buf: PathBuf,
    prefix_length: Option<u32>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        oneshot::Receiver<u64>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Configure fastp for deduplication only
    let mut command_fields = HashMap::from([
        ("-q".to_string(), Some("20".to_string())),
        ("--qualified_quality_phred".to_string(), Some("15".to_string())),
        ("--complexity_threshold".to_string(), Some("30".to_string())),
        ("--stdin".to_string(), None), // Input from stdin
        ("--stdout".to_string(), None), // Output to stdout
        ("--dedup".to_string(), None), // Enable deduplication
        ("--dup_calc_accuracy".to_string(), Some("3".to_string())), // Low memory (~1-2GB)
        ("-w".to_string(), Some(RunConfig::thread_allocation(&config, FASTP_TAG, None).to_string())),
        ("-j".to_string(), Some(out_dir.join(format!("{}_dedup.json", sample_base_buf.file_name().unwrap().to_string_lossy())).to_string_lossy().to_string())), // JSON report
    ]);

    // Handle paired-end
    if paired {
        command_fields.insert("--interleaved_in".to_string(), None);
    }

    // Handle prefix-length by trimming (if specified)
    if let Some(len) = prefix_length {
        command_fields.insert("--cut_tail".to_string(), Some(len.to_string()));
        command_fields.insert("--cut_front".to_string(), Some("0".to_string()));
        debug!("Trimming reads to {}bp for prefix-based deduplication", len);
    }

    let fastp_config = FastpConfig { command_fields };
    let fastp_args = generate_cli(FASTP_TAG, &config, Some(&fastp_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;

    // Stream to fastp
    let (mut fastp_child, fastp_stream_task, fastp_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        FASTP_TAG,
        fastp_args,
        StreamDataType::JustBytes, // Or IlluminaFastq if parsing FASTQ records
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(fastp_stream_task);
    cleanup_tasks.push(fastp_err_task);

    // Parse output stream
    let dedup_stream = {
        let mut guard = fastp_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Count unique reads from JSON report
    let json_path = out_dir.join(format!("{}_dedup.json", sample_base_buf.file_name().unwrap().to_string_lossy()));
    let (count_tx, count_rx) = oneshot::channel();
    let count_task = tokio::spawn(async move {
        let json_content = match fs::read_to_string(&json_path).await {
            Ok(content) => content,
            Err(e) => {
                return Err(anyhow!("JSON read error"));
            }
        };
        let parsed: Value = match serde_json::from_str(&json_content) {
            Ok(v) => v,
            Err(e) => {
                error!("JSON parse error: {}", e);
                count_tx.send(0).map_err(|_| anyhow!("Send failed"))?; // Fallback
                return Ok(());
            }
        };

        let uniques = parsed["duplication"]["unique_reads"].as_u64().unwrap_or_else(|| {
            warn!("Missing unique_reads; falling back to 0");
            0
        });
        let unique_pairs = uniques / if paired { 2 } else { 1 };
        count_tx.send(unique_pairs).map_err(|_| anyhow!("Send failed"))?;
        Ok(())
    });
    cleanup_tasks.push(count_task);

    Ok((
        ReceiverStream::new(dedup_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
    ))
}


/// Combined de-duplication and subsambler.
/// Uses and A-Res algorithm and hashing to creat a min heap BinaryHeap
/// Runs in memory when possible, but has a disk option when
/// file size and memory constraints demand.
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Strucutred FASTQ stream
/// * 'paired` - Bool controlling single or paired end data
/// * `seed` - Seed value for rng allowing possibility of deterministic output
/// * `max_subsample` - Limit of reads for subsampling
/// * `prefix_len` - Length of prefix subsample of read
/// * `out_dir` - Output directory for intermediate files
///
/// # Returns
/// Result of:
/// Tuple of:
/// Output stream
/// Count stream
/// Vec of tkio tasks
/// Vec of tokio receivers
/// with Pipeline Error
async fn dedup_and_subsample(
    config: Arc<RunConfig>,
    input_stream: mpsc::Receiver<ParseOutput>,
    paired: bool,
    seed: u64, // Used for hashing, rng from config
    max_subsample: u64,
    prefix_len: Option<usize>,
    out_dir: PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let mut count_map: HashMap<u64, u64> = HashMap::new();
    let mut dedup_map: HashMap<u64, (u64, Vec<SequenceRecord>)> = HashMap::new();
    let mut total_count = 0u64;

    // RAM-based deduplication
    let mut stream = ReceiverStream::new(input_stream);
    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Fastq(record) => {
                let mut records = vec![record.clone()];
                let mut seq_bytes = (*record.seq()).to_vec(); // Clone Vec<u8> to owned Vec<u8>
                if paired {
                    if let Some(next_item) = stream.next().await {
                        match next_item {
                            ParseOutput::Fastq(next_record) => {
                                records.push(next_record.clone());
                                seq_bytes.extend_from_slice(&next_record.seq()); // Extend with &[u8]
                            }
                            _ => {
                                return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in paired stream".to_string()));
                            }
                        }
                    } else {
                        return Err(PipelineError::InvalidFastqFormat("Unpaired read in paired stream".to_string()));
                    }
                }
                let mut hasher = XxHash64::with_seed(seed);
                let write_len = prefix_len.map(|len| len.min(seq_bytes.len())).unwrap_or(seq_bytes.len());
                hasher.write(&seq_bytes[0..write_len]);
                let hash_key = hasher.finish();
                let entry = dedup_map.entry(hash_key).or_insert((0, records.clone()));
                entry.0 += 1;
                entry.1 = records;
                total_count += 1;
            }
            _ => {
                return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in dedup stream".to_string()));
            }
        }
    }
    count_map = dedup_map.iter().map(|(&k, &(w, _))| (k, w)).collect();
    debug!("RAM-based count complete: {} unique hashes from {} total {}", count_map.len(), total_count, if paired { "pairs" } else { "reads" });

    let subsample_size = max_subsample.min(count_map.len() as u64);

    let mut heap: BinaryHeap<Reverse<SampleItem>> = BinaryHeap::new();
    let mut rng = config.rng.clone(); // Clone the seeded RNG from config

    if subsample_size < count_map.len() as u64 {
        for (&hash_key, &weight) in &count_map {
            let key = rng.random::<f64>().powf(1.0 / weight as f64);
            if heap.len() < subsample_size as usize {
                heap.push(Reverse(SampleItem { key, records: vec![], weight, hash_key: Some(hash_key) }));
            } else if key > heap.peek().unwrap().0.key {
                heap.pop();
                heap.push(Reverse(SampleItem { key, records: vec![], weight, hash_key: Some(hash_key) }));
            }
        }
    } else {
        // All uniques
        for (&hash_key, &weight) in &count_map {
            heap.push(Reverse(SampleItem { key: 0.0, records: vec![], weight, hash_key: Some(hash_key) }));
        }
    }

    let mut sampled: Vec<SampleItem> = heap.into_iter().map(|r| r.0).collect();
    sampled.sort_by(|a, b| a.key.partial_cmp(&b.key).unwrap_or(Ordering::Equal));

    for item in sampled.iter_mut() {
        if let Some(hash_key) = item.hash_key {
            if let Some((_, records)) = dedup_map.get(&hash_key) {
                item.records = records.clone();
            }
        }
    }

    let tsv_path = out_dir.join("duplicate_cluster_sizes.tsv");
    let dedup_map_clone = dedup_map.clone();
    let tsv_write_task = tokio::spawn(async move {
        let mut file = TokioOpenOptions::new().write(true).create(true).open(&tsv_path).await.map_err(|e| anyhow!("TSV open error: {}", e))?;
        for (_, &(weight, ref records)) in &dedup_map_clone {
            if !records.is_empty() {
                file.write_all(format!("{}\t{}\n", records[0].id(), weight).as_bytes()).await?;
            }
        }
        Ok(())
    });
    cleanup_tasks.push(tsv_write_task);

    let unique_count = sampled.len() as u64;
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(unique_count).map_err(|_| PipelineError::Other(anyhow!("Failed to send subsample count")))?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let send_task = tokio::spawn(async move {
        for item in sampled {
            for record in item.records {
                tx.send(ParseOutput::Fastq(record)).await.map_err(|_| anyhow!("Send failed"))?;
            }
        }
        Ok(())
    });
    cleanup_tasks.push(send_task);

    Ok((ReceiverStream::new(rx), count_rx, cleanup_tasks, cleanup_receivers))
}

/// Alignment of deduped/subsampled stream to a non-host DB (canoncially NT)
/// for output in PAF format, allowing downstream metagenomic calling.
///
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream (paired or single-end).
///
/// # Returns
/// Tuple:
/// - Interleaved FASTQ stream : unmapped reads.
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
/// - Optional temporary FASTA file for non-host reference.
/// - Optional temporary minimap2 index file.
async fn minimap2_non_host_align(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Option<NamedTempFile>,
        Option<NamedTempFile>,
        Option<TempDir>
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let (_non_host_ref_fasta_path, non_host_index_path, non_host_ref_temp, non_host_index_temp, non_host_temp_dir,mut non_host_ref_tasks) =
        minimap2_index_prep(
            &config,
            &config.ram_temp_dir,
            config.args.target_sequence.clone(),
            config.args.target_index.clone(),
            "non-host",
        )
            .await?;
    try_join_all(non_host_ref_tasks).await?;

    debug!("Non-host index: {}", non_host_index_path.display());

    let minimap2_config = Minimap2Config {
        minimap2_index_path: non_host_index_path,
        option_fields: HashMap::from([
            ("-cx".to_string(), Some("sr".to_string())), // Short-read preset, PAF out
            ("--secondary".to_string(), Some("yes".to_string())), // allow secondary alignments
        ]),
    };

    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut minimap2_child, minimap2_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Lines,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    Ok((
        ReceiverStream::new(minimap2_out_stream),
        cleanup_tasks,
        cleanup_receivers,
        non_host_ref_temp,
        non_host_index_temp,
        non_host_temp_dir
    ))
}

/// Converts a PAF stream to an m8 stream, splits it with t_junction for file writing,
/// and returns the main m8 stream along with cleanup handles.
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream (paired or single-end).
///
/// # Returns
/// Tuple:
/// - m8 stream
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
pub async fn paf_to_m8(
    config: Arc<RunConfig>,
    mut input_stream: ReceiverStream<ParseOutput>,
    m8_path: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
)> {
    // Channel for converted m8 lines
    let (m8_tx, m8_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);

    // Clone config for closure — fixes E0382
    let config_clone = config.clone();

    // Spawn PAF → m8 conversion task
    let conversion_handle = tokio::spawn(async move {
        while let Some(item) = input_stream.next().await {
            let line_bytes = match item {
                ParseOutput::Bytes(b) => b,
                _ => {
                    debug!("Skipping non-bytes item in PAF stream");
                    continue;
                }
            };

            let line = match String::from_utf8((*line_bytes).clone()) {
                Ok(s) => s,
                Err(e) => {
                    debug!("Invalid UTF-8 in PAF line: {}", e);
                    continue;
                }
            };

            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            let record = match PafRecord::parse_line(&line) {
                Ok(r) => r,
                Err(e) => {
                    debug!("Failed to parse PAF line: {} — {}", e, line);
                    continue;
                }
            };

            // Use user-provided nt_db_size — no fallback
            let genome_size = config_clone.args.nt_db_size as f64;
            let m8_line = record.to_m8_line(genome_size);
            let m8_bytes = (m8_line + "\n").into_bytes();

            if m8_tx.send(ParseOutput::Bytes(Arc::new(m8_bytes))).await.is_err() {
                return Err(anyhow!("m8 channel send failed"));
            }
        }
        Ok(())
    });

    let m8_stream = ReceiverStream::new(m8_rx);

    // Split stream: one for file, one for downstream
    let (mut streams, done_rx) = t_junction(
        m8_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "paf_to_m8_stream".to_string(),
        None,
    )
        .await?;

    let main_stream = streams.remove(0);
    let file_stream = streams.remove(0);

    // Write to file
    let write_task = write_byte_stream_to_file(
        &m8_path,
        ReceiverStream::new(file_stream),
        Some(config.base_buffer_size),
    )
        .await?;

    Ok((
        ReceiverStream::new(main_stream),
        vec![conversion_handle, write_task],
        vec![done_rx],
    ))
}


/// Run function for Short Read mNGS pipelines
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
///
/// # Returns
/// Result<(), PipelineError>
/// Run function for Short Read mNGS pipelines
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
///
/// # Returns
/// Result<(), PipelineError>
pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();
    let mut cleanup_tasks: Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>> = Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let mut temp_dirs: Vec<TempDir> = Vec::new();

    info!("Starting short read mNGS pipeline.");

    // *******************
    // Setup and Validation
    // *******************

    // External tools check
    check_versions(vec![BOWTIE2_TAG, MINIMAP2_TAG, KALLISTO_TAG, CZID_DEDUP_TAG])
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Check required files
    let host_bowtie2_index: String = config.args.host_bowtie2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_bowtie2_index is required".to_string()))?;

    let (file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;
    let paired = file2_path.is_some();

    let seed = config.args.seed.unwrap_or_else(|| {
        let mut bytes = [0u8; 8];
        OsRng.fill_bytes(&mut bytes);
        let random_seed = u64::from_le_bytes(bytes);
        random_seed
    });

    // Input Validation
    let (val_out_stream, validate_cleanup_tasks, validate_cleanup_receivers, raw_count_task, val_count_task) = validate_input(
        config.clone(),
        file1_path,
        file2_path,
        no_ext_sample_base_buf.clone(),
        &out_dir,
    ).await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // *******************
    // Host Filtering Stage
    // *******************

    // ERCC bt2 filtering and count
    let ercc_bt2_index_path = bowtie2_index_prep(&config.args.ercc_bowtie2_index, &cwd)?;
    let ercc_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
    let (ercc_bt2_out_stream, ercc_count_rx, ercc_bt2_cleanup_tasks, ercc_bt2_cleanup_receivers, _ercc_bt2_insert_stats_rx) = bowtie2_filter(
        config.clone(),
        val_out_stream,
        ercc_bt2_index_path,
        paired,
        ercc_bt2_options,
        None,
        false
    ).await?;
    cleanup_tasks.extend(ercc_bt2_cleanup_tasks);
    cleanup_receivers.extend(ercc_bt2_cleanup_receivers);

    // QC with fastp
    let (qc_fastp_out_stream, qc_cleanup_tasks, qc_cleanup_receivers, qc_count_result_rx) = fastp_qc(
        config.clone(),
        ercc_bt2_out_stream,
    ).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);
    cleanup_receivers.extend(qc_cleanup_receivers);

    // Split for Kallisto and Bowtie2
    let (mut kallisto_streams, kallisto_split_done_rx) = t_junction(
        qc_fastp_out_stream,
        2,
        config.base_buffer_size * 10,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "kallisto_split".to_string(),
        None,
    ).await?;
    cleanup_receivers.push(kallisto_split_done_rx);

    let mut kallisto_streams_iter = kallisto_streams.into_iter();
    let kallisto_bypass_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let kallisto_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let kallisto_stream = ReceiverStream::new(kallisto_stream);

    let (kallisto_ercc_tx, kallisto_ercc_rx, kallisto_cleanup_tasks, kallisto_cleanup_receivers, kallisto_exit_task) = kallisto_quant(
        config.clone(),
        kallisto_stream,
        out_dir.clone(),
        paired,
        no_ext_sample_base_buf.clone(),
    ).await?;
    cleanup_tasks.extend(kallisto_cleanup_tasks);
    cleanup_receivers.extend(kallisto_cleanup_receivers);

    // Host filtering: bt2
    let host_bt2_index_path = bowtie2_index_prep(host_bowtie2_index, &cwd)?;
    let host_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
    let (host_bt2_out_stream, host_bt2_count_rx, host_bt2_cleanup_tasks, host_bt2_cleanup_receivers, host_bt2_insert_stats_rx) = bowtie2_filter(
        config.clone(),
        ReceiverStream::new(kallisto_bypass_stream),
        host_bt2_index_path,
        paired,
        host_bt2_options,
        None,
        true
    ).await?;
    cleanup_tasks.extend(host_bt2_cleanup_tasks);
    cleanup_receivers.extend(host_bt2_cleanup_receivers);


    // Host filtering: mm2
    let (host_mm2_out_stream, host_mm2_count_rx, mut host_mm2_cleanup_tasks, mut host_mm2_cleanup_receivers, host_ref_temp, host_index_temp, host_temp_dir) = minimap2_filter(
        config.clone(),
        host_bt2_out_stream,
        paired,
        None, // No BAM output; set to Some(PathBuf::from("host_mm2.bam")) if needed
    )
        .await?;
    cleanup_tasks.append(&mut host_mm2_cleanup_tasks);
    cleanup_receivers.append(&mut host_mm2_cleanup_receivers);

    if let Some(temp) = host_ref_temp {
        temp_files.push(temp);
    }
    if let Some(temp) = host_index_temp {
        temp_files.push(temp);
    }
    if let Some(temp) = host_temp_dir {
        temp_dirs.push(temp);
    }

    // If host is no huma, run an additional filter stage using a human reference
    let post_filter_stream  = if config.args.human_host {
        host_mm2_out_stream
    }
    else {
        let human_bowtie2_index: String = config.args.human_bowtie2_index.clone();

        // human filtering: bt2
        let human_bt2_index_path = bowtie2_index_prep(human_bowtie2_index, &cwd)?;
        let human_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
        let (human_bt2_out_stream, human_bt2_count_rx, human_bt2_cleanup_tasks, human_bt2_cleanup_receivers, _human_bt2_insert_stats_rx) = bowtie2_filter(
            config.clone(),
            host_mm2_out_stream,
            human_bt2_index_path,
            paired,
            human_bt2_options,
            None,
            false
        ).await?;
        cleanup_receivers.extend(human_bt2_cleanup_receivers);

        // human filtering: mm2
        let (human_mm2_out_stream, human_mm2_count_rx, mut human_mm2_cleanup_tasks, mut human_mm2_cleanup_receivers, human_ref_temp, human_index_temp, human_index_dir) = minimap2_filter(
            config.clone(),
            human_bt2_out_stream,
            paired,
            None, // No BAM output; set to Some(PathBuf::from("human_mm2.bam")) if needed
        )
            .await?;
        cleanup_tasks.append(&mut human_mm2_cleanup_tasks);
        cleanup_receivers.append(&mut human_mm2_cleanup_receivers);

        if let Some(temp) = human_ref_temp {
            temp_files.push(temp);
        }
        if let Some(temp) = human_index_temp {
            temp_files.push(temp);
        }
        if let Some(temp) = human_index_dir {
            temp_dirs.push(temp);
        }

        human_mm2_out_stream

    };



    let (dedup_stream, dedup_count_rx, dedup_cleanup_tasks, dedup_cleanup_receivers) = dedup_and_subsample(
        config.clone(),
        post_filter_stream.into_inner(),
        paired,
        seed,
        config.args.max_subsample as u64,
        Some(75), // Prefix length for deduplication. Hardcoded fior now
        out_dir.clone(),
    ).await?;
    cleanup_tasks.extend(dedup_cleanup_tasks);
    cleanup_receivers.extend(dedup_cleanup_receivers);




    // *******************
    // Non host Alignment
    // *******************

    // Minimap2 non_host alignment

    let (non_host_mm2_out_stream, mut non_host_mm2_cleanup_tasks, mut non_host_mm2_cleanup_receivers, non_host_ref_temp, non_host_index_temp, non_host_index_dir) = minimap2_non_host_align(
        config.clone(),
        dedup_stream,
    )
        .await?;
    cleanup_tasks.append(&mut non_host_mm2_cleanup_tasks);
    cleanup_receivers.append(&mut non_host_mm2_cleanup_receivers);

    if let Some(temp) = non_host_ref_temp {
        temp_files.push(temp);
    }
    if let Some(index_temp) = non_host_index_temp {
        temp_files.push(index_temp);
    }
    if let Some(index_temp) = non_host_index_dir {
        temp_dirs.push(index_temp);
    }

    let m8_file_path = file_path_manipulator(
        &PathBuf::from(&no_ext_sample_base_buf),
        Some(&out_dir),
        None,
        Some(".m8"),
        "",
    );

    let (m8_stream, mut m8_cleanup_tasks, mut m8_cleanup_receivers) = paf_to_m8(
        config.clone(),
        non_host_mm2_out_stream,
        m8_file_path,
    )
        .await?;
    cleanup_tasks.append(&mut m8_cleanup_tasks);
    cleanup_receivers.append(&mut m8_cleanup_receivers);


    // Write test stream
    let test_write_task = tokio::spawn(stream_to_file(
        m8_stream.into_inner(),
        out_dir.join("non_host.m8"),
    ));
    cleanup_tasks.push(test_write_task);

    // *******************
    // Results retrieval
    // *******************

    let raw_count = join_with_error_handling(raw_count_task).await?;
    info!("Processed {} raw reads (additive from R1 and R2 if paired)", raw_count);

    let stats = join_with_error_handling(val_count_task).await?;
    info!("Processed {} validated, {} undersized, {} oversized reads",
             stats.validated, stats.undersized, stats.oversized);

    let _qc_fastp_read_count = match qc_count_result_rx.await {
        Ok(Ok(count)) => {
            info!("Received fastp read count: {}", count);
            count
        }
        Ok(Err(e)) => {
            error!("Error getting fastp read count: {}", e);
            0
        }
        Err(e) => {
            error!("Count receiver dropped: {}", e);
            0
        }
    };


    let ercc_mapped_count = ercc_count_rx.await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC count receiver failed: {}", e)))?;

    let host_bt2_counts = host_bt2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host bt2 counts receiver failed: {}", e)))?;
    info!("Host bt2: mapped counts: {:?}", host_bt2_counts);

    if let Some(insert_stats_rx) = host_bt2_insert_stats_rx {
        match insert_stats_rx.await {
            Ok(insert_stats) if !insert_stats.insert_sizes.is_empty() => {
                info!("Host bt2 insert size stats: {:?}", insert_stats.insert_sizes);
                let output_path = out_dir.join("insert_size_histogram.png");
                let sample_name = no_ext_sample_base_buf
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .to_string();
                let plot_task = tokio::spawn(async move {
                    plot_insert_sizes(&insert_stats.insert_sizes, &sample_name, &output_path)
                        .map_err(|e| anyhow!("Failed to plot insert sizes: {}", e))?;
                    Ok(())
                });
                cleanup_tasks.push(plot_task);
            }
            Ok(_) => warn!("No insert size stats available; skipping plotting"),
            Err(e) => error!("Host bt2 insert stats receiver failed: {}", e),
        }
    }

    let host_mm2_counts = host_mm2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host minimap2 counts receiver failed: {}", e)))?;
    info!("Host minimap2: mapped counts: {:?}", host_mm2_counts);

    let dedup_count = dedup_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("dedup count receiver failed: {}", e)))?;
    info!("Dedup: unique reads/pairs: {}", dedup_count);

    // Await Kallisto exit and process results
    join_with_error_handling(kallisto_exit_task).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let kallisto_results_task = kallisto_results(out_dir.join("kallisto"), kallisto_ercc_tx)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    join_with_error_handling(kallisto_results_task).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Collect Kallisto results
    let kallisto_ercc_counts = kallisto_ercc_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC counts receiver failed: {}", e)))?;
    // info!("Kallisto ERCC counts: {:?}", kallisto_ercc_counts);

    // *******************
    // Cleanup
    // *******************
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

    drop(temp_files);
    drop(temp_dirs);

    info!("Finished short read mNGS pipeline.");
    Ok(())
}
