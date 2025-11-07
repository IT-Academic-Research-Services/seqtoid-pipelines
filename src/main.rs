mod pipelines;
mod utils;
mod config;
mod cli;

use sysinfo::CpuRefreshKind;
use std::time::{Instant, SystemTime};
use std::{env, fs};
use chrono::{DateTime};
use std::cmp::min;
use std::path::PathBuf;
use std::sync::Arc;
use sysinfo::{System, RefreshKind, MemoryRefreshKind};
use std::io::Write;

use anyhow::Result;
use log::{self, LevelFilter, debug, info, error, warn};
use env_logger::Builder;
use clap::Parser;
use rayon::ThreadPool;
use tokio::sync::Semaphore;
use rayon::ThreadPoolBuilder;
use tokio::process::Command as TokioCommand;
use tokio::time::{sleep, Duration};
use std::process::Output;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_core::{RngCore, OsRng};
use crate::cli::parse;
use crate::config::defs::{RunConfig, StreamDataType, PipelineError};
use crate::cli::args::Technology;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::r1r2_base;
use crate::utils::system::{detect_cores_and_load, compute_stream_threads, detect_ram, generate_rng, compute_base_buffer_size, get_ram_temp_dir};
use pipelines::consensus_genome;
use pipelines::short_read_mngs;
use pipelines::db;


#[tokio::main]
async fn main() -> Result<()> {
    let run_start = Instant::now();

    #[cfg(not(unix))]
    Err(anyhow!("Named pipes are not supported on non-Unix systems."));

    let args = parse();

    let log_level = if args.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    Builder::new()
        .filter_level(log_level)
        .format(|buf, record| {
            writeln!(
                buf,
                "[{}] {}: {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                record.level(),
                record.args()
            )
        })
        .init();

    println!("\n-------------\n SeqToID\n-------------\n");

    let dir = env::current_dir()?;
    info!("The current directory is {:?}\n", dir);

    let ram_temp_dir = get_ram_temp_dir();
    info!("The RAM temp directory is {:?}\n", ram_temp_dir);

    let (max_cores, cpu_load) = detect_cores_and_load(args.threads).await?;
    let stream_threads = compute_stream_threads(max_cores, cpu_load, args.threads);
    debug!("Detected {} physical cores; CPU load {}%; using {} threads for pool, {} for streams",
              max_cores, cpu_load, max_cores, stream_threads);

    let thread_pool = Arc::new(ThreadPoolBuilder::new()
        .num_threads(max_cores)
        .build()
        .expect("Failed to create thread pool"));

    let maximal_semaphore = Arc::new(Semaphore::new(2));

    let (total_ram, available_ram) = detect_ram()?;
    debug!("Available RAM: {} bytes (~{} GiB)", available_ram, available_ram / 1_073_741_824);
    debug!("Total RAM: {} bytes (~{} GiB)", total_ram, total_ram / 1_073_741_824);

    let input_size_mb = get_input_size_mb(&args.file1, &args.file2).unwrap_or(0);
    debug!("Total input file size: {} MB", input_size_mb);

    let rng = generate_rng(args.seed);

    let sdt = match args.technology {
        Technology::Illumina => StreamDataType::IlluminaFastq,
        Technology::ONT => StreamDataType::OntFastq,
    };
    let base_buffer_size = compute_base_buffer_size(available_ram, total_ram, sdt, stream_threads);

    let out_dir = setup_output_dir(&args, &dir)?;
    let module = args.module.clone();
    let run_config = Arc::new(RunConfig {
        cwd: dir,
        ram_temp_dir,
        out_dir,
        args,
        thread_pool,
        maximal_semaphore,
        base_buffer_size,
        input_size_mb,
        available_ram,
        rng,
        log_level
    });

    if let Err(e) = match module.as_str() {
        "consensus_genome" => consensus_genome_run(run_config).await,
        "short_read_mngs" => short_read_mngs_run(run_config).await,
        "build_taxid_lineages" => build_taxid_lineages_run(run_config).await,
        "build_accession2taxid" => build_accession2taxid_run(run_config).await,
        _ => Err(PipelineError::InvalidConfig(format!("Invalid module: {}", module))),
    } {
        error!("Pipeline failed: {} at {} milliseconds.", e, run_start.elapsed().as_millis());
        std::process::exit(1);
    }

    println!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
    Ok(())
}



async fn consensus_genome_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    consensus_genome::run(run_config).await
}

async fn short_read_mngs_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    short_read_mngs::run(run_config).await
}

async fn build_taxid_lineages_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    db::taxid_lineages_db(run_config).await
}

async fn build_accession2taxid_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    db::accession2taxid_db(run_config).await
}

/// Sets up output directory
/// If `out_dir` is specified from args, uses it;
/// otherwise, creates a directory named `<sample_base>_YYYYMMDD`.
/// Ensures the directory exists.
///
/// # Arguments
/// * `args` - The parsed command-line arguments.
/// * `cwd` - The current working directory.
/// # Returns
/// path to the output directory.
fn setup_output_dir(args: &cli::args::Arguments, cwd: &PathBuf) -> Result<PathBuf> {
    let out_dir = match &args.out_dir {
        Some(out) => {
            let path = PathBuf::from(out);
            if path.is_absolute() {
                path
            } else {
                cwd.join(path)
            }
        }
        None => {
            let file1_path: PathBuf = match &args.file1 {
                Some(file) => {
                    let file1_full_path = file_path_manipulator(&PathBuf::from(file), Some(cwd), None, None, "");
                    if file1_full_path.exists() {
                        file1_full_path
                    } else {
                        return Err(anyhow::anyhow!("Cannot find file 1 (-i)"));
                    }
                }
                None => return Err(anyhow::anyhow!("File1 path required")),
            };
            let file1_r1r2 = r1r2_base(&file1_path);

            let dir_base = match file1_r1r2.prefix {
                Some(prefix) => prefix,
                None => {
                    info!("No R1 tag found. Using bare file 1 stem as sample_base.");
                    file1_r1r2.file_name.unwrap_or_else(|| {
                        file1_path
                            .file_stem()
                            .map(|stem| stem.to_string_lossy().into_owned())
                            .unwrap_or_else(|| "default_sample".to_string())
                    })
                }
            };

            let timestamp = SystemTime::now()
                .duration_since(SystemTime::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .map(|secs| {
                    let dt = DateTime::from_timestamp(secs as i64, 0)
                        .unwrap_or_else(|| DateTime::from_timestamp(0, 0).unwrap());
                    dt.format("%Y%m%d").to_string()
                })
                .unwrap_or_else(|_| "19700101".to_string());
            cwd.join(format!("{}_{}", dir_base, timestamp))
        }
    };
    fs::create_dir_all(&out_dir)?;
    Ok(out_dir)
}


/// Creates a pool of threads for running sub-processes.
///
/// # Arguments
/// # 'max_cores' : Maximum cores allowed or discovered.
///
/// # Returns
/// rayon::ThreadPool
fn create_thread_pool(max_cores: usize) -> ThreadPool {
    ThreadPoolBuilder::new()
        .num_threads(max_cores)
        .build()
        .expect("Failed to create thread pool")
}


/// Gets the file size(s) of the input file(s) from metadata
///
/// # Arguments
/// # 'file1' : First (probably not optional) input.
/// # 'file2' : Second input.
///
/// # Returns
/// size in MB
fn get_input_size_mb(file1: &Option<String>, file2: &Option<String>) -> Result<u64> {
    let mut total_size = 0u64;
    if let Some(f1) = file1 {
        total_size += fs::metadata(f1)?.len();
    }
    if let Some(f2) = file2 {
        total_size += fs::metadata(f2)?.len();
    }
    Ok(total_size / 1_048_576) // Bytes -> MB
}