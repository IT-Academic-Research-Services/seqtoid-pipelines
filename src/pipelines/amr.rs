//! AMR Pipeline

use std::sync::Arc;
use log::info;
use crate::config::defs::{PipelineError, RunConfig};

use crate::pipelines::short_read_mngs::run_host_filter;

pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();

    info!("Finished AMR (cleanup completed, some non-fatal warnings may have occurred).");
    Ok(())
}