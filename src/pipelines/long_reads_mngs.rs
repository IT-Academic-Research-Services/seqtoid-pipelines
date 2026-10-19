//! Long Read Metagenomics Pipeline

use std::path::PathBuf;
use std::sync::Arc;
use anyhow::anyhow;
use futures::future::try_join_all;
use log::{debug, info, warn};
use tempfile::{NamedTempFile, TempDir};
use tokio::task::JoinHandle;
use crate::config::defs::{PipelineError, RunConfig};
use crate::utils::file::validate_file_inputs;

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


    let (input_fastq, _file2, sample_base_buf, sample_base) =
        validate_file_inputs(&config, &cwd, true).await?;













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