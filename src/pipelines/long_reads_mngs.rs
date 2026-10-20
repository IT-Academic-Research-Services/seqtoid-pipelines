//! Long Read Metagenomics Pipeline

use std::path::PathBuf;
use std::sync::Arc;


use anyhow::anyhow;
use futures::future::try_join_all;
use log::{debug, info, warn};
use tempfile::{NamedTempFile, TempDir};
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;

use crate::config::defs::{PipelineError, RunConfig, PairingMode, ReadStats, StreamDataType, BOWTIE2_TAG, MINIMAP2_TAG, KALLISTO_TAG, SPADES_TAG, MAKEBLASTDB_TAG, BLASTN_TAG, BLASTX_TAG, SORT_TAG, SAMTOOLS_TAG, FASTP_TAG, HISAT2_TAG};
use crate::cli::args::Technology;
use crate::utils::command::check_versions;
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file};
use crate::utils::fastx::{read_fastq};
use crate::utils::streams::{fanout_to_channels, join_with_error_handling, ParseOutput};

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



    let (_validation_and_qc, ) = tokio::try_join!(
        validation_and_qc_fut,
    )?;

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