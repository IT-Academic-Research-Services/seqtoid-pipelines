use std::sync::Arc;
use crate::config::defs::{PipelineError, RunConfig};

pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {

    println!("Finished short read mNGS.");
    Ok(())
}
