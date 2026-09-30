use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use clap::{Parser, ValueEnum};
use tokio::fs;

use seqtoid_pipelines::distributed::executor::{
    WorkUnitClaimed,
    WorkerBackend,
    WorkerExecutor,
    WorkerExecutorConfig,
};

use seqtoid_pipelines::utils::work_units::{
    WorkUnit,
    WorkUnitState,
};

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

    /// Directory containing WorkUnit JSON files for this run.
    #[arg(long)]
    work_dir: PathBuf,

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
        "Starting seqtoid-worker: worker_id={}, backend={:?}, work_dir={}",
        worker_id,
        args.backend,
        args.work_dir.display()
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

    loop {
        let work_unit_paths = list_work_unit_paths(&args.work_dir).await?;

        match find_first_available_work_unit(&work_unit_paths).await? {
            Some((path, mut work_unit)) => {
                log::info!(
                "Found AVAILABLE work unit: {} \
                 sample={}, chunk={}, attempt={}",
                path.display(),
                work_unit.sample_id,
                work_unit.chunk_id,
                work_unit.attempt,
            );

                match executor
                    .claim_and_execute(&path, &mut work_unit)
                    .await
                {
                    Ok(_) => {}

                    Err(err) if err.downcast_ref::<WorkUnitClaimed>().is_some() => {
                        log::debug!(
            "Work unit {} was claimed by another worker; trying another",
            path.display()
        );
                    }

                    Err(err) => {
                        return Err(err).context("worker execution failed");
                    }
                }
            }

            None => {
                if has_work_in_progress(&work_unit_paths).await? {
                    log::debug!(
            "No AVAILABLE work units currently available; other workers are still processing"
        );

                    tokio::time::sleep(
                        tokio::time::Duration::from_secs(2)
                    )
                        .await;

                    continue;
                }

                log::info!(
        "No AVAILABLE, CLAIMED, or RUNNING work units remain"
    );
                break;
            }


        }
    }


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


/// Find serialized WorkUnit files in a run's work directory.
async fn list_work_unit_paths(work_dir: &PathBuf) -> Result<Vec<PathBuf>> {
    let mut entries = fs::read_dir(work_dir)
        .await
        .with_context(|| {
            format!(
                "failed to read work directory {}",
                work_dir.display()
            )
        })?;

    let mut paths = Vec::new();

    while let Some(entry) = entries
        .next_entry()
        .await
        .with_context(|| {
            format!(
                "failed while reading work directory {}",
                work_dir.display()
            )
        })?
    {
        let path = entry.path();

        let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };

        if path.is_file()
            && file_name.starts_with("work_")
            && file_name.ends_with(".json")
        {
            paths.push(path);
        }
    }

    paths.sort();

    Ok(paths)
}

/// Find the first AVAILABLE WorkUnit in deterministic filename order.
async fn find_first_available_work_unit(
    work_unit_paths: &[PathBuf],
) -> Result<Option<(PathBuf, WorkUnit)>> {
    for path in work_unit_paths {
        let work_unit = load_work_unit(path).await?;

        if work_unit.state == WorkUnitState::Available {
            return Ok(Some((path.clone(), work_unit)));
        }
    }

    Ok(None)
}

/// Determine whether any work units are still in progress.
async fn has_work_in_progress(
    work_unit_paths: &[PathBuf],
) -> Result<bool> {
    for path in work_unit_paths {
        let work_unit = load_work_unit(path).await?;

        if matches!(
            work_unit.state,
            WorkUnitState::Claimed | WorkUnitState::Running
        ) {
            return Ok(true);
        }
    }

    Ok(false)
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