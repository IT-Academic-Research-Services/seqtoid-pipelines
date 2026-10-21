//! Long Read Metagenomics Pipeline

use tokio_stream::StreamExt;
use std::path::PathBuf;
use std::sync::Arc;
use std::collections::HashMap;

use anyhow::anyhow;
use futures::future::try_join_all;
use log::{debug, info, warn};
use tempfile::{NamedTempFile, TempDir};
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;

use crate::config::defs::{PipelineError, RunConfig, PairingMode, ReadStats, StreamDataType,
                           MINIMAP2_TAG, KALLISTO_TAG, SPADES_TAG, MAKEBLASTDB_TAG, BLASTN_TAG,
                          BLASTX_TAG, SORT_TAG, SAMTOOLS_TAG, FASTP_TAG, HISAT2_TAG,
                          SamtoolsSubcommand, FLYE_TAG};
use crate::cli::args::{Technology, LibraryType};
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::minimap2::{Minimap2Config};
use crate::utils::command::samtools::{SamtoolsConfig};
use crate::utils::command::minimap2::{minimap2_index_prep};
use crate::utils::command::flye::FlyeConfig;
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file,
                         choose_temp_dir};
use crate::utils::fastx::{read_fastq, StreamRecordStats, stream_record_counter, fastq_to_fasta_stream};
use crate::utils::sambam::filter_clipped_alignments;
use crate::utils::streams::{fanout_to_channels, join_with_error_handling, ParseOutput,
                            stream_to_cmd, parse_child_output, ChildStream, ParseMode,
                            ChannelReader, parse_fastq, spawn_cmd};
use crate::utils::subsample::subsample_uniform_reservoir;


struct ReadsToContigsOutputs {
    reads_to_contigs_sam: PathBuf,
    reads_to_contigs_bam: PathBuf,
    reads_to_contigs_bai: PathBuf,
    reads_to_contigs_tsv: PathBuf,
    non_contigs_fastq: PathBuf,
    non_contigs_fasta: PathBuf,
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
/// - handle for the validated read count task
async fn validate_input(
    config: Arc<RunConfig>,
    input_fastq: PathBuf,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> anyhow::Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
        Vec<tokio::sync::oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
        JoinHandle<anyhow::Result<ReadStats, anyhow::Error>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let validated_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("sample_validated.fastq.gz"),
        "_",
    );


    let (validated_rx, val_count_task) = read_fastq(
        input_fastq.clone(),
        None,
        PairingMode::Strict,
        Some(Technology::ONT),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        "run_validate_input",
        &config,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(validated_rx);

    let (val_streams, val_router_done_rx) = fanout_to_channels(
        val_rx_stream,
        2,
        "validate_input",
        &config,
        StreamDataType::IlluminaFastq,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(val_router_done_rx);

    let mut val_streams_iter = val_streams.into_iter();

    let val_out_stream =
        ReceiverStream::new(val_streams_iter.next().ok_or(PipelineError::EmptyStream)?);

    let val_file_stream =
        ReceiverStream::new(val_streams_iter.next().ok_or(PipelineError::EmptyStream)?);

    let val_file_write_task = write_byte_stream_to_file(
        &validated_file_path,
        val_file_stream,
        config.clone(),
        StreamDataType::OntFastq,
        "validate_input",
        false
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(val_file_write_task);


    Ok((
        val_out_stream,
        cleanup_tasks,
        cleanup_receivers,
        val_count_task,
    ))
}


async fn fastp_qc(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
        Vec<tokio::sync::oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
        tokio::sync::oneshot::Receiver<anyhow::Result<StreamRecordStats, anyhow::Error>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let fastp_config = FastpConfig {
        command_fields: HashMap::from([
            ("--disable_adapter_trimming".to_string(), None),
            (
                "--qualified_quality_phred".to_string(),
                Some("9".to_string()),
            ),
            (
                "--length_required".to_string(),
                Some("100".to_string()),
            ),
            ("--low_complexity_filter".to_string(), None),
            (
                "--complexity_threshold".to_string(),
                Some("30".to_string()),
            ),
            ("--dont_eval_duplication".to_string(), None),
        ]),
        paired: false,
    };

    let fastp_args =
        generate_cli(FASTP_TAG, &config, Some(&fastp_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (fastp_child, fastp_stream_task, fastp_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        FASTP_TAG,
        fastp_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(fastp_stream_task);
    cleanup_tasks.push(fastp_err_task);

    let fastp_out_stream = {
        let mut guard = fastp_child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Tee the fastp output:
    //   branch 1 -> next pipeline stage
    //   branch 2 -> streamed read/base counter
    let tee_count_input = ReceiverStream::new(fastp_out_stream);

    let (tee_streams, tee_router_done_rx) = fanout_to_channels(
        tee_count_input,
        2,
        "fastp_count",
        &config,
        StreamDataType::OntFastq,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "fanout_to_channels".to_string(),
            error: e.to_string(),
        })?;

    cleanup_receivers.push(tee_router_done_rx);

    let mut tee_streams_iter = tee_streams.into_iter();

    let qc_out_stream = ReceiverStream::new(
        tee_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let qc_count_stream = ReceiverStream::new(
        tee_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let (count_tx, count_result_rx) =
        tokio::sync::oneshot::channel::<Result<StreamRecordStats, anyhow::Error>>();

    let count_task = tokio::spawn(async move {
        match stream_record_counter(qc_count_stream.into_inner(), false).await {
            Ok(stats) => {
                debug!(
                    "Fastp output reads: {}, bases: {}",
                    stats.reads,
                    stats.bases
                );

                let _ = count_tx.send(Ok(stats));
                Ok(())
            }

            Err(e) => {
                let _ = count_tx.send(Err(e));
                Err(anyhow!("Failed to count fastp output reads and bases"))
            }
        }
    });

    cleanup_tasks.push(count_task);

    Ok((
        qc_out_stream,
        cleanup_tasks,
        cleanup_receivers,
        count_result_rx,
    ))
}

async fn minimap2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    host_index_path: PathBuf,
    output_bam_path: PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        tokio::sync::oneshot::Receiver<StreamRecordStats>,
        Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
        Vec<tokio::sync::oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // ------------------------------------------------------------------
    // minimap2
    //
    // The WDL selects splice-aware alignment for RNA and map-ont for DNA.
    // Input comes directly from the previous FASTQ stream.
    // ------------------------------------------------------------------

    let library_type = config.args.library_type;
    let minimap_preset = match library_type {
        LibraryType::DNA => "map-ont",
        LibraryType::RNA => "splice",
    };

    let minimap2_config = Minimap2Config {
        minimap2_index_path: host_index_path,
        r1_path: None,
        r2_path: None,
        option_fields: HashMap::from([(
            "-ax".to_string(),
            Some(minimap_preset.to_string()),
        )]),
        num_threads: None,
    };

    let minimap2_args =
        generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (minimap2_child, minimap2_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::OntFastq,
        config.args.verbose,
        None,
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
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // ------------------------------------------------------------------
    // The minimap2 SAM stream feeds two independent samtools processes:
    //
    //   1. fastq -n -f 4 -> downstream unmapped FASTQ stream
    //   2. view -b        -> host-filter BAM file
    // ------------------------------------------------------------------
    let (sam_streams, sam_router_done_rx) = fanout_to_channels(
        ReceiverStream::new(minimap2_out_stream),
        2,
        "host_filter",
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(sam_router_done_rx);

    let mut sam_streams_iter = sam_streams.into_iter();

    // ------------------------------------------------------------------
    // Branch 1: unmapped FASTQ
    // ------------------------------------------------------------------
    let unmapped_sam_stream = ReceiverStream::new(
        sam_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-f".to_string(), Some("4".to_string())),
            ("-".to_string(), None),
        ]),
    };

    let samtools_fastq_args =
        generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (
        samtools_fastq_child,
        samtools_fastq_task,
        samtools_fastq_err_task,
    ) = stream_to_cmd(
        config.clone(),
        unmapped_sam_stream.into_inner(),
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(samtools_fastq_task);
    cleanup_tasks.push(samtools_fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = samtools_fastq_child.lock().await;

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

    // ------------------------------------------------------------------
    // Tee unmapped FASTQ:
    //
    //   branch 1 -> next pipeline stage
    //   branch 2 -> read/base counter
    // ------------------------------------------------------------------
    let (count_streams, count_router_done_rx) = fanout_to_channels(
        ReceiverStream::new(unmapped_fastq_stream),
        2,
        "host_filter_count",
        &config,
        StreamDataType::OntFastq,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(count_router_done_rx);

    let mut count_streams_iter = count_streams.into_iter();

    let host_filtered_stream = ReceiverStream::new(
        count_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let host_filtered_count_stream = ReceiverStream::new(
        count_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let (count_tx, count_rx) =
        tokio::sync::oneshot::channel::<StreamRecordStats>();

    let count_task = tokio::spawn(async move {
        let stats = stream_record_counter(
            host_filtered_count_stream.into_inner(),
            false,
        )
            .await?;

        debug!(
            "Host filter output: {} reads, {} bases",
            stats.reads,
            stats.bases
        );

        let _ = count_tx.send(stats);

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(count_task);

    // ------------------------------------------------------------------
    // Branch 2: BAM
    // ------------------------------------------------------------------
    let mapped_sam_stream = ReceiverStream::new(
        sam_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let samtools_view_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-S".to_string(), None),
            ("-b".to_string(), None),
            ("-".to_string(), None),
        ]),
    };

    let samtools_view_args =
        generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_view_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (
        samtools_view_child,
        samtools_view_task,
        samtools_view_err_task,
    ) = stream_to_cmd(
        config.clone(),
        mapped_sam_stream.into_inner(),
        SAMTOOLS_TAG,
        samtools_view_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(samtools_view_task);
    cleanup_tasks.push(samtools_view_err_task);

    let bam_stream = {
        let mut guard = samtools_view_child.lock().await;

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

    let bam_write_task = write_byte_stream_to_file(
        &output_bam_path,
        ReceiverStream::new(bam_stream),
        config.clone(),
        StreamDataType::JustBytes,
        "host_filter_bam",
        false,
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(bam_write_task);

    Ok((
        host_filtered_stream,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
    ))
}


async fn run_reads_to_contigs(
    config: Arc<RunConfig>,
    input_fastq: PathBuf,
    assembled_reads: PathBuf,
    out_dir: &PathBuf,
) -> Result<
    (
        ReadsToContigsOutputs,
        Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
        Vec<tokio::sync::oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // ================================================================
    // Durable outputs
    // ================================================================

    let reads_to_contigs_sam =
        out_dir.join("sample.reads_to_contigs.sam");

    let reads_to_contigs_bam =
        out_dir.join("sample.reads_to_contigs.bam");

    let reads_to_contigs_bai =
        out_dir.join("sample.reads_to_contigs.bam.bai");

    let reads_to_contigs_tsv =
        out_dir.join("reads_to_contigs.tsv");

    let non_contigs_fastq =
        out_dir.join("sample.non_contigs.fastq");

    let non_contigs_fasta =
        out_dir.join("sample.non_contigs.fasta");


    // ================================================================
    // minimap2
    //
    // Legacy:
    //
    // minimap2 -ax map-ont assembled_reads input_fastq
    //          -t 15 --secondary=no
    //
    // We preserve the map-ont preset and secondary=no.
    // ================================================================

    let minimap2_config = Minimap2Config {
        minimap2_index_path: assembled_reads.clone(),
        r1_path: Some(input_fastq),
        r2_path: None,
        option_fields: HashMap::from([
            (
                "-ax".to_string(),
                Some("map-ont".to_string()),
            ),
            (
                "--secondary=no".to_string(),
                None,
            ),
        ]),
        num_threads: Some(15),
    };

    let minimap2_args = generate_cli(
        MINIMAP2_TAG,
        &config,
        Some(&minimap2_config),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    info!(
        "RunReadsToContigs minimap2 args: {:?}",
        minimap2_args
    );

    let (mut minimap2_child, minimap2_err_task) =
        spawn_cmd(
            config.clone(),
            MINIMAP2_TAG,
            minimap2_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?;

    cleanup_tasks.push(minimap2_err_task);

    let minimap2_sam_rx = parse_child_output(
        &mut minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        &config,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;


    // ================================================================
    // samtools view -h -F 0x900
    //
    // Keep headers and discard secondary + supplementary alignments.
    // ================================================================

    let primary_view_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-h".to_string(), None),
            ("-F".to_string(), Some("0x900".to_string())),
            ("-".to_string(), None),
        ]),
    };

    let primary_view_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&primary_view_config),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (
        primary_child,
        primary_stream_task,
        primary_err_task,
    ) = stream_to_cmd(
        config.clone(),
        minimap2_sam_rx,
        SAMTOOLS_TAG,
        primary_view_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(primary_stream_task);
    cleanup_tasks.push(primary_err_task);

    let primary_sam_rx = {
        let mut guard = primary_child.lock().await;

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


    // ================================================================
    // filter_clipped_alignments.py equivalent
    //
    // Legacy invocation:
    //
    //     -m 20
    // ================================================================

    let (
        filtered_sam_stream,
        clipped_filter_task,
    ) = filter_clipped_alignments(
        ReceiverStream::new(primary_sam_rx),
        20.0,
        config.base_buffer_size,
    )
        .await?;

    cleanup_tasks.push(clipped_filter_task);


    // ================================================================
    // Fan out FINAL filtered SAM.
    //
    // Consumers:
    //
    //  1. durable SAM
    //  2. reads_to_contigs.tsv
    //  3. BAM pipeline
    //  4. non-contig FASTQ extraction
    // ================================================================

    let (filtered_sam_streams, filtered_sam_done_rx) =
        fanout_to_channels(
            filtered_sam_stream,
            4,
            "reads_to_contigs_filtered_sam",
            &config,
            StreamDataType::JustBytes,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(filtered_sam_done_rx);

    let mut filtered_sam_streams =
        filtered_sam_streams.into_iter();


    // ================================================================
    // Branch 1: durable filtered SAM
    // ================================================================

    let sam_file_stream = ReceiverStream::new(
        filtered_sam_streams
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let sam_write_task = write_byte_stream_to_file(
        &reads_to_contigs_sam,
        sam_file_stream,
        config.clone(),
        StreamDataType::JustBytes,
        "reads_to_contigs_sam",
        false,
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(sam_write_task);


    // ================================================================
    // Branch 2: reads_to_contigs.tsv
    //
    // Legacy:
    //
    // grep -v "^@" sample.reads_to_contigs.sam |
    // awk '{ print $1 "\t" $3 "\t" length($10)}'
    // ================================================================

    let mut tsv_input = ReceiverStream::new(
        filtered_sam_streams
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let tsv_path = reads_to_contigs_tsv.clone();

    let tsv_task = tokio::spawn(async move {
        use tokio::io::AsyncWriteExt;

        let file = tokio::fs::File::create(&tsv_path).await?;
        let mut writer = tokio::io::BufWriter::new(file);

        while let Some(item) = tsv_input.next().await {
            let ParseOutput::Bytes(bytes) = item else {
                continue;
            };

            if bytes.first() == Some(&b'@') {
                continue;
            }

            let line = std::str::from_utf8(&bytes)?
                .trim_end_matches(&['\r', '\n'][..]);

            let fields: Vec<&str> =
                line.split('\t').collect();

            if fields.len() < 10 {
                return Err(anyhow!(
                    "Malformed SAM record while generating reads_to_contigs.tsv"
                ));
            }

            let output = format!(
                "{}\t{}\t{}\n",
                fields[0],
                fields[2],
                fields[9].len(),
            );

            writer.write_all(output.as_bytes()).await?;
        }

        writer.flush().await?;

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(tsv_task);


    // ================================================================
    // Branch 3:
    //
    // samtools view -b -
    //      |
    // samtools sort -
    //
    // ================================================================

    let bam_input = filtered_sam_streams
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    let bam_view_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-b".to_string(), None),
            ("-".to_string(), None),
        ]),
    };

    let bam_view_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&bam_view_config),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (
        bam_view_child,
        bam_view_stream_task,
        bam_view_err_task,
    ) = stream_to_cmd(
        config.clone(),
        bam_input,
        SAMTOOLS_TAG,
        bam_view_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(bam_view_stream_task);
    cleanup_tasks.push(bam_view_err_task);

    let bam_view_rx = {
        let mut guard = bam_view_child.lock().await;

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


    // samtools sort

    let sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
        ]),
    };

    let sort_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&sort_config),
    )
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
        bam_view_rx,
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

    let sorted_bam_rx = {
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

    let bam_write_task = write_byte_stream_to_file(
        &reads_to_contigs_bam,
        ReceiverStream::new(sorted_bam_rx),
        config.clone(),
        StreamDataType::JustBytes,
        "reads_to_contigs_bam",
        false,
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;


    // ================================================================
    // Branch 4: non-contig reads
    //
    // samtools fastq -n -f 4 -
    // ================================================================

    let non_contig_sam_input =
        filtered_sam_streams
            .next()
            .ok_or(PipelineError::EmptyStream)?;

    let non_contig_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-f".to_string(), Some("4".to_string())),
            ("-".to_string(), None),
        ]),
    };

    let non_contig_fastq_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&non_contig_fastq_config),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (
        non_contig_child,
        non_contig_stream_task,
        non_contig_err_task,
    ) = stream_to_cmd(
        config.clone(),
        non_contig_sam_input,
        SAMTOOLS_TAG,
        non_contig_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(non_contig_stream_task);
    cleanup_tasks.push(non_contig_err_task);

    let non_contig_fastq_rx = {
        let mut guard = non_contig_child.lock().await;

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


    // ================================================================
    // Fan out non-contig FASTQ:
    //
    //   1. durable FASTQ
    //   2. FASTQ -> FASTA
    // ================================================================

    let (non_contig_streams, non_contig_done_rx) =
        fanout_to_channels(
            ReceiverStream::new(non_contig_fastq_rx),
            2,
            "non_contig_fastq",
            &config,
            StreamDataType::OntFastq,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(non_contig_done_rx);

    let mut non_contig_streams =
        non_contig_streams.into_iter();

    let non_contig_file_stream =
        ReceiverStream::new(
            non_contig_streams
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );

    let non_contig_fasta_input =
        ReceiverStream::new(
            non_contig_streams
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );


    // Durable FASTQ

    let non_contig_fastq_write_task =
        write_byte_stream_to_file(
            &non_contigs_fastq,
            non_contig_file_stream,
            config.clone(),
            StreamDataType::OntFastq,
            "non_contigs_fastq",
            false,
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(e.to_string())
            })?;

    cleanup_tasks.push(non_contig_fastq_write_task);


    // FASTQ -> FASTA

    let (
        non_contig_fasta_stream,
        non_contig_fasta_convert_task,
    ) = fastq_to_fasta_stream(
        non_contig_fasta_input,
        config.base_buffer_size,
    )
        .await?;

    cleanup_tasks.push(non_contig_fasta_convert_task);

    let non_contig_fasta_write_task =
        write_byte_stream_to_file(
            &non_contigs_fasta,
            non_contig_fasta_stream,
            config.clone(),
            StreamDataType::JustBytes,
            "non_contigs_fasta",
            false,
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(e.to_string())
            })?;

    cleanup_tasks.push(non_contig_fasta_write_task);


    // ================================================================
    // BAM must be completely materialized before samtools index starts.
    // ================================================================

    bam_write_task
        .await
        .map_err(|e| PipelineError::Other(e.into()))??;


    // ================================================================
    // samtools index
    // ================================================================

    let index_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Index,
        subcommand_fields: HashMap::from([
            (
                reads_to_contigs_bam
                    .to_string_lossy()
                    .to_string(),
                None,
            ),
            (
                reads_to_contigs_bai
                    .to_string_lossy()
                    .to_string(),
                None,
            ),
        ]),
    };

    let index_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&index_config),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut index_child, index_err_task) =
        spawn_cmd(
            config.clone(),
            SAMTOOLS_TAG,
            index_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

    let index_status =
        index_child
            .wait()
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

    index_err_task
        .await
        .map_err(|e| PipelineError::Other(e.into()))??;

    if !index_status.success() {
        return Err(PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: format!(
                "samtools index exited with status {}",
                index_status
            ),
        });
    }


    Ok((
        ReadsToContigsOutputs {
            reads_to_contigs_sam,
            reads_to_contigs_bam,
            reads_to_contigs_bai,
            reads_to_contigs_tsv,
            non_contigs_fastq,
            non_contigs_fasta,
        },
        cleanup_tasks,
        cleanup_receivers,
    ))
}



/// The main entry point for the long-read metagenomic NGS (mNGS) pipeline.
///
/// # Arguments
///
/// * `config`: the run configuration
///
/// # Returns
///
/// Result<()>: success or error
pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir()
        .map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();

    let efs_run_dir = config.efs_runs_dir.join(&config.run_id);
    let efs_work_dir = efs_run_dir.join("work");


    let mut cleanup_tasks: Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<tokio::sync::oneshot::Receiver<anyhow::Result<(), anyhow::Error>>> =
        Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let mut host_filter_temp_dirs: Vec<TempDir> = Vec::new();
    let mut final_temp_dirs: Vec<TempDir> = Vec::new();
    let temp_paths: Vec<PathBuf> = Vec::new();

    info!("Starting long read mNGS pipeline.");

    // *******************
    // Setup and Validation
    // *******************

    // External tools check

    let mut versions_vec = vec![
        FASTP_TAG,
    ];

    if !config.args.no_version_check {
        check_versions(versions_vec, &out_dir.clone(), &config)
            .await
            .map_err(|e| PipelineError::Other(e.into()))?;
    }


    let (input_fastq, _file2, sample_base_buf, sample_base) =
        validate_file_inputs(&config, &cwd, true).await?;


    // Input Validation
    let (
        val_out_stream,
        validate_cleanup_tasks,
        validate_cleanup_receivers,
        val_count_task,
    ) = validate_input(
        config.clone(),
        input_fastq.clone(),
        sample_base_buf.clone(),
        &out_dir,
    )
        .await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // QC with fastp
    let (qc_fastp_out_stream, qc_cleanup_tasks, qc_cleanup_receivers, qc_count_result_rx) =
        fastp_qc(config.clone(), val_out_stream).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);


    // *******************
    // Host Filtering Stage
    // *******************

    //Host filtering

    let (
        _host_ref_fasta_path,
        host_index_path,
        _host_ref_temp,
        _host_index_temp,
        host_temp_dir,
        _host_ref_tasks,
    ) = minimap2_index_prep(
        &config,
        &config.ram_temp_dir,
        config.args.host_sequence.clone(),
        config.args.host_mm2_index.clone(),
        "host",
    )
        .await?;

    debug!("Host index : {}", host_index_path.display());

    let (
        host_filter_out_stream,
        host_filter_stats,
        host_filter_cleanup_tasks,
        host_filter_cleanup_receivers,
    ) = minimap2_filter(
        config.clone(),
        qc_fastp_out_stream,
        host_index_path,
        out_dir.join("sample.hostfiltered.bam"),
    )
        .await?;

    cleanup_tasks.extend(host_filter_cleanup_tasks);
    cleanup_receivers.extend(host_filter_cleanup_receivers);
    final_temp_dirs.extend(host_temp_dir);

    //Human filtering

    let (
        _human_ref_fasta_path,
        human_index_path,
        _human_ref_temp,
        _human_index_temp,
        human_temp_dir,
        _human_ref_tasks,
    ) = minimap2_index_prep(
        &config,
        &config.ram_temp_dir,
        config.args.human_sequence.clone(),
        config.args.human_mm2_index.clone(),
        "host",
    )
        .await?;

    debug!("Human index : {}", human_index_path.display());

    let (
        human_filter_out_stream,
        human_filter_stats,
        human_filter_cleanup_tasks,
        human_filter_cleanup_receivers,
    ) = minimap2_filter(
        config.clone(),
        host_filter_out_stream,
        human_index_path,
        out_dir.join("sample.humanfiltered.bam"),
    )
        .await?;

    cleanup_tasks.extend(human_filter_cleanup_tasks);
    cleanup_receivers.extend(human_filter_cleanup_receivers);
    final_temp_dirs.extend(human_temp_dir);

    // Human filter output remains raw FASTQ bytes through both minimap2 stages.
    // Parse only once here because reservoir sampling operates on FASTQ records.
    let human_filtered_byte_reader =
        ChannelReader::new(human_filter_out_stream.into_inner());

    let human_filtered_fastq_rx = parse_fastq(
        human_filtered_byte_reader,
        &config,
        StreamDataType::OntFastq,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!(
    "Failed to parse human-filtered FASTQ stream: {}",
    e
)))?;

    let human_filtered_fastq_stream =
        ReceiverStream::new(human_filtered_fastq_rx);

    // *******************
    // Subsampling Stage
    // *******************
    //
    // Legacy RunSubsampling operates on the output of RunHumanFilter.
    // Long reads are single-end, so each FASTQ record is one sampling unit.

    let (
        subsampled_stream,
        subsample_count_rx,
        subsample_task,
    ) = subsample_uniform_reservoir(
        config.clone(),
        human_filtered_fastq_stream.into_inner(),
        config.args.max_subsample as u64,
        false,
    )
        .await?;

    cleanup_tasks.push(subsample_task);

    // Split the subsampled stream into:
    //
    //   1. pre-assembly FASTA
    //   2. materialized FASTQ for downstream file-based consumers
    //
    let (subsampled_streams, subsample_router_done_rx) = fanout_to_channels(
        subsampled_stream,
        2,
        "post_subsample",
        &config,
        StreamDataType::OntFastq,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(subsample_router_done_rx);

    let mut subsampled_streams_iter = subsampled_streams.into_iter();

    let subsampled_fasta_input_stream =
        ReceiverStream::new(
            subsampled_streams_iter
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );

    let subsampled_fastq_input_stream =
        ReceiverStream::new(
            subsampled_streams_iter
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );


    // *******************
    // Assembly Stage
    // *******************

    // Choose scratch for the materialized subsampled FASTQ and Flye work dir.
    let assembly_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        1,
        true,
    )
        .await?;

    let pre_assembly_fasta_path =
        out_dir.join("pre_assembly.fa");

    let flye_out_dir =
        assembly_temp_dir.path().join("temp_flye_out");

    let subsampled_fastq_path =
        assembly_temp_dir.path().join("subsampled.fastq");

    final_temp_dirs.push(assembly_temp_dir);

    //convert FASTQ to FASTA

    let (
        pre_assembly_fasta_stream,
        fastq_to_fasta_task,
    ) = fastq_to_fasta_stream(
        subsampled_fasta_input_stream,
        config.base_buffer_size,
    )
        .await?;

    cleanup_tasks.push(fastq_to_fasta_task);


    let pre_assembly_fasta_write_task = write_byte_stream_to_file(
        &pre_assembly_fasta_path,
        pre_assembly_fasta_stream,
        config.clone(),
        StreamDataType::OntFastq,
        "pre_assembly_fasta",
        false,
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(pre_assembly_fasta_write_task);


    // Materialize the subsampled FASTQ into scratch.
    let subsampled_fastq_write_task = write_byte_stream_to_file(
        &subsampled_fastq_path,
        subsampled_fastq_input_stream,
        config.clone(),
        StreamDataType::OntFastq,
        "pre_assembly_fastq",
        false,
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;


    // Flye needs the complete input file, so this is the synchronization point.
    subsampled_fastq_write_task
        .await
        .map_err(|e| PipelineError::Other(e.into()))??;



    let flye_config = FlyeConfig {
        input_fastq: subsampled_fastq_path.clone(),
        outdir_path: flye_out_dir.clone(),
        guppy_basecaller_setting: config.args.guppy_basecaller_setting.clone(),
        polishing_iterations: config.args.polishing_iterations,
    };

    let flye_args =
        generate_cli(
            FLYE_TAG,
            &config,
            Some(&flye_config),
        )
            .map_err(|e| PipelineError::ToolExecution {
                tool: FLYE_TAG.to_string(),
                error: e.to_string(),
            })?;

    info!("Starting Flye assembly: {:?}", flye_args);

    let (mut flye_child, flye_err_task) =
        crate::utils::streams::spawn_cmd(
            config.clone(),
            FLYE_TAG,
            flye_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FLYE_TAG.to_string(),
                error: e.to_string(),
            })?;

    // Flye must finish before anything can consume its assembly.
    let flye_status =
        flye_child
            .wait()
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FLYE_TAG.to_string(),
                error: e.to_string(),
            })?;

    if flye_status.success() {
        info!("Flye assembly completed successfully");
    } else {
    // Flye failure itself is NOT fatal.
        warn!(
        "Flye exited with status {}; continuing because assembly failure is non-fatal",
        flye_status
    );
    }

    flye_err_task
        .await
        .map_err(|e| PipelineError::Other(e.into()))??;


    // ------------------------------------------------------------------
    // Resolve Flye assembly output.
    //
    // Legacy WDL:
    //
    //     if [ -f temp_flye_out/assembly.fasta ]
    //     then
    //         cat temp_flye_out/assembly.fasta > sample.assembled_reads.fasta
    //     else
    //         touch sample.assembled_reads.fasta
    //     fi
    //
    //     mv sample.assembled_reads.fasta contigs.fasta
    //
    // Therefore the durable pipeline artifact is always contigs.fasta,
    // which is empty when Flye did not produce assembly.fasta.
    // ------------------------------------------------------------------

    let flye_assembly_path =
        flye_out_dir.join("assembly.fasta");

    let contigs_fasta_path =
        out_dir.join("contigs.fasta");

    if flye_assembly_path.exists() {
        tokio::fs::copy(
            &flye_assembly_path,
            &contigs_fasta_path,
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(format!(
                    "Failed to copy Flye assembly {} to {}: {}",
                    flye_assembly_path.display(),
                    contigs_fasta_path.display(),
                    e
                ))
            })?;

        info!(
        "Flye assembly copied to {}",
        contigs_fasta_path.display()
    );
    } else {
        warn!(
        "Flye did not produce {}; creating empty contigs.fasta",
        flye_assembly_path.display()
    );

        tokio::fs::File::create(&contigs_fasta_path)
            .await
            .map_err(|e| {
                PipelineError::IOError(format!(
                    "Failed to create empty {}: {}",
                    contigs_fasta_path.display(),
                    e
                ))
            })?;
    }

    // *******************
    // Write checkpoint
    // *******************


    // *******************
    // Temporary final stream drain for construction
    // *******************

    // let mut final_temp_drain_stream = subsampled_out_stream;
    //
    // let final_temp_drain_task = tokio::spawn(async move {
    //     while final_temp_drain_stream.next().await.is_some() {}
    //     Ok::<(), anyhow::Error>(())
    // });
    //
    // cleanup_tasks.push(final_temp_drain_task);


    // *******************
    // Results retrieval
    // *******************

    let validation_and_qc_fut = async move {

        let stats = join_with_error_handling(val_count_task).await?;
        info!(
            "Processed {} raw reads, {} raw bases, {} validated reads, {} validated bases, {} undersized, {} oversized reads",
            stats.raw_reads, stats.raw_bases, stats.validated, stats.validated_bases, stats.undersized, stats.oversized
        );

        Ok::<_, PipelineError>((stats))
    };

    let human_filter_stats_fut = async move {
        let stats = human_filter_stats
            .await
            .map_err(|e| PipelineError::Other(anyhow!("Human filter stats channel failed: {}", e)))?;

        info!(
        "Human filter: {} reads, {} bases",
        stats.reads,
        stats.bases
    );

        Ok::<_, PipelineError>(stats)
    };

    let subsample_counts_fut = async move {

        let subsample_count = subsample_count_rx
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Subsample count receiver failed: {}",
                e
            ))
            })?;


        info!(
        "Subsampled reads (pairs counted separately): {}",
        subsample_count
    );

        Ok::<_, PipelineError>((
            subsample_count,
        ))
    };


    let (_validation_stats, human_filtered_stats, _subsample_stats) = tokio::try_join!(
        validation_and_qc_fut,
        human_filter_stats_fut,
        subsample_counts_fut
    )?;

    info!("Human filtered stats: {:?}", human_filtered_stats);

    // *******************
    // Cleanup
    // *******************

    // Helper to handle downcast + NotFound tolerance in one place
    fn handle_cleanup_error(e: anyhow::Error, context: &str) -> anyhow::Result<(), PipelineError> {
        match e.downcast::<std::io::Error>() {
            Ok(io_err) => {
                if io_err.kind() == std::io::ErrorKind::NotFound {
                    warn!(
                        "Cleanup ignored missing file/dir (non-fatal) [{}]: {}",
                        context, io_err
                    );
                    Ok(())
                } else {
                    // Other io errors are fatal
                    Err(PipelineError::Other(anyhow::Error::new(io_err)))
                }
            }
            Err(remaining) => {
                // Not an io::Error → treat as fatal
                Err(PipelineError::Other(remaining))
            }
        }
    }

    // Wait for all cleanup tasks
    let cleanup_results = try_join_all(cleanup_tasks)
        .await
        .map_err(|join_err| PipelineError::Other(join_err.into()))?;

    // Process each task result
    for res in cleanup_results {
        if let Err(e) = res {
            handle_cleanup_error(e, "task")?;
        }
    }

    // Process cleanup receivers (oneshot channels)
    for rx in cleanup_receivers {
        match rx.await {
            Ok(Ok(_)) => {} // success

            Ok(Err(e)) => {
                handle_cleanup_error(e, "receiver")?;
            }

            Err(oneshot_err) => {
                // Channel broken (panic, drop, empty stream, etc.)
                // Usually safe to ignore if stream was empty
                warn!(
                    "Cleanup receiver channel failed (possibly empty): {}",
                    oneshot_err
                );
                // continue (non-fatal)
            }
        }
    }


    // Final temp resource drops — defensive, never fail the pipeline
    for tf in temp_files {
        if let Err(e) = tf.close() {
            if e.kind() == std::io::ErrorKind::NotFound {
                debug!("Temp file already gone during close: {}", e);
            } else {
                warn!("Temp file close failed (non-fatal): {}", e);
            }
        }
    }

    for td in final_temp_dirs {
        if let Err(e) = std::fs::remove_dir_all(td.path()) {
            if e.kind() == std::io::ErrorKind::NotFound {
                debug!("Temp dir already gone: {}", e);
            } else {
                warn!("Temp dir removal failed (non-fatal): {}", e);
            }
        }
    }

    drop(temp_paths);


    info!("Finished long read mNGS pipeline (cleanup completed, some non-fatal warnings may have occurred).");
    Ok(())
}