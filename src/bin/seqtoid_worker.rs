use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use clap::{Parser, ValueEnum};
use tokio::fs;

use seqtoid_pipelines::distributed::executor::{
    WorkerBackend,
    WorkerExecutor,
    WorkerExecutorConfig,
};
use seqtoid_pipelines::utils::work_units::WorkUnit;

#[derive(Debug, Clone, ValueEnum)]
enum Backend {
    MmseqsCpu,
    MmseqsGpu,
    Diamond,
}

impl From<Backend> for WorkerBackend {
    fn from(value: Backend) -> Self {
        match value {
            Backend::MmseqsCpu => WorkerBackend::MmseqsCpu,
            Backend::MmseqsGpu => WorkerBackend::MmseqsGpu,
            Backend::Diamond => WorkerBackend::Diamond,
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "seqtoid-worker",
    version,
    about = "Execute one distributed SeqToID NR work unit"
)]
struct WorkerArguments {
    /// Backend used to execute the work unit.
    #[arg(long, value_enum)]
    backend: Backend,

    /// JSON file containing one WorkUnit.
    #[arg(long)]
    work_unit: PathBuf,

    /// Worker-local scratch directory.
    #[arg(long, default_value = "/scratch")]
    scratch_dir: PathBuf,

    /// Durable directory in which completed results and metadata are written.
    #[arg(long)]
    results_dir: PathBuf,

    /// Number of threads assigned to the worker execution.
    #[arg(long, default_value_t = 112)]
    threads: usize,

    /// Maximum worker cores used by the shared NUMA policy.
    #[arg(long, default_value_t = 128)]
    max_cores: usize,

    /// Worker-local MMseqs database.
    ///
    /// Required for MMseqs CPU/GPU execution.
    #[arg(long)]
    mmseqs_db: Option<PathBuf>,

    /// Worker-local Diamond database prefix/path.
    ///
    /// Required for Diamond execution.
    #[arg(long)]
    diamond_db: Option<PathBuf>,

    /// Explicit worker ID.
    ///
    /// If omitted, the executable attempts to determine the EC2 instance ID.
    #[arg(long)]
    worker_id: Option<String>,
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::init();

    let args = WorkerArguments::parse();

    let worker_id = match args.worker_id {
        Some(id) => id,
        None => ec2_instance_id()
            .await
            .context("worker-id was not supplied and EC2 instance ID could not be determined")?,
    };

    log::info!(
        "Starting seqtoid-worker: worker_id={}, backend={:?}, work_unit={}",
        worker_id,
        args.backend,
        args.work_unit.display()
    );

    let work_unit = load_work_unit(&args.work_unit).await?;

    log::info!(
        "Loaded work unit {}: sample={}, chunk={}, paired_end={}, state={:?}, attempt={}",
        work_unit.id(),
        work_unit.sample_id,
        work_unit.chunk_id,
        work_unit.paired_end,
        work_unit.state,
        work_unit.attempt
    );

    let executor_config = WorkerExecutorConfig {
        worker_id,
        scratch_dir: args.scratch_dir,
        results_dir: args.results_dir,
        threads: args.threads,
        max_cores: args.max_cores,
        mmseqs_db: args.mmseqs_db,
        diamond_db: args.diamond_db,
    };

    let executor = WorkerExecutor::new(
        executor_config,
        args.backend.into(),
    );

    let mut work_unit = work_unit;

    executor
        .claim_and_execute(&mut work_unit)
        .await
        .context("worker execution failed")?;

    Ok(())
}

/// Load one serialized WorkUnit from disk.
async fn load_work_unit(path: &PathBuf) -> Result<WorkUnit> {
    let bytes = fs::read(path)
        .await
        .with_context(|| {
            format!(
                "failed to read work-unit file {}",
                path.display()
            )
        })?;

    serde_json::from_slice(&bytes)
        .with_context(|| {
            format!(
                "invalid WorkUnit JSON in {}",
                path.display()
            )
        })
}

/// Determine the EC2 instance ID using IMDSv2.
///
/// This is only used when --worker-id is not supplied. Supplying the worker
/// ID explicitly is useful for local development and unit/integration tests.
async fn ec2_instance_id() -> Result<String> {
    use tokio::process::Command;

    let token_output = Command::new("curl")
        .args([
            "-sS",
            "-X",
            "PUT",
            "-H",
            "X-aws-ec2-metadata-token-ttl-seconds: 300",
            "http://169.254.169.254/latest/api/token",
        ])
        .output()
        .await
        .context("failed to query EC2 metadata token")?;

    if !token_output.status.success() {
        return Err(anyhow!(
            "EC2 metadata token request failed with status {:?}",
            token_output.status.code()
        ));
    }

    let token = String::from_utf8(token_output.stdout)
        .context("EC2 metadata token was not valid UTF-8")?
        .trim()
        .to_string();

    if token.is_empty() {
        return Err(anyhow!("EC2 metadata token was empty"));
    }

    let instance_output = Command::new("curl")
        .args([
            "-sS",
            "-H",
            &format!("X-aws-ec2-metadata-token: {}", token),
            "http://169.254.169.254/latest/meta-data/instance-id",
        ])
        .output()
        .await
        .context("failed to query EC2 instance ID")?;

    if !instance_output.status.success() {
        return Err(anyhow!(
            "EC2 instance ID request failed with status {:?}",
            instance_output.status.code()
        ));
    }

    let instance_id = String::from_utf8(instance_output.stdout)
        .context("EC2 instance ID was not valid UTF-8")?
        .trim()
        .to_string();

    if instance_id.is_empty() {
        return Err(anyhow!("EC2 instance ID was empty"));
    }

    Ok(instance_id)
}