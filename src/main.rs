mod pipelines;
mod utils;
mod config;

use std::time::{Instant, SystemTime};
use std::{env, fs};
use std::cmp::min;
use std::path::PathBuf;
use std::sync::Arc;
use sysinfo::System;

use anyhow::Result;
use clap::Parser;
use num_cpus;
use rayon::ThreadPool;
use tokio::sync::Semaphore;
use tokio::runtime::Builder;
use rayon::ThreadPoolBuilder;

use crate::cli::parse;
use crate::config::defs::{RunConfig, StreamDataType, PipelineError};
use crate::cli::args::Technology;
use crate::utils::file::extension_remover;
use pipelines::consensus_genome;
use pipelines::db;

mod cli;

#[tokio::main]
async fn main() -> Result<()> {
    let run_start = Instant::now();
    println!("\n-------------\n SeqToID\n-------------\n");

    #[cfg(not(unix))]
    Err(anyhow!("Named pipes are not supported on non-Unix systems."));

    let dir = env::current_dir()?;
    println!("The current directory is {:?}\n", dir);

    let ram_temp_dir = get_ram_temp_dir();
    println!("The RAM temp directory is {:?}\n", ram_temp_dir);

    let args = parse();

    let max_cores = min(num_cpus::get(), args.threads);
    let thread_pool = Arc::new(create_thread_pool(max_cores));
    let maximal_semaphore = Arc::new(Semaphore::new(2));

    let mut system = System::new_all();
    system.refresh_memory();
    let available_ram = system.available_memory();

    let sdt = match args.technology {
        Technology::Illumina => StreamDataType::IlluminaFastq,
        Technology::ONT => StreamDataType::OntFastq,
    };
    let base_buffer_size = compute_base_buffer_size(available_ram, max_cores, sdt);

    let out_dir = setup_output_dir(&args, &dir)?;
    let module = args.module.clone();
    let run_config = Arc::new(RunConfig {
        cwd: dir,
        ram_temp_dir,
        out_dir,
        args,
        thread_pool,
        maximal_semaphore,
        base_buffer_size
    });

    if let Err(e) = match module.as_str() {
        "consensus_genome" => consensus_genome_run(run_config.clone()).await,
        "create_db" => create_db_run(&run_config).await,
        _ => Err(PipelineError::InvalidConfig(format!("Invalid module: {}", module))),
    } {
        eprintln!("Pipeline failed: {}", e);
        std::process::exit(1);
    }

    println!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
    Ok(())
}

async fn consensus_genome_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {

    // consensus_genome::old_run(&run_config)
    //     .await
    //     .map_err(|e| PipelineError::Other(e.into()))

    consensus_genome::run(run_config).await
}

async fn create_db_run(run_config: &RunConfig) -> Result<(), PipelineError> {
    db::create_db(run_config).await
        .map_err(|e| PipelineError::Other(e.into()))
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
            let file1_path = args.file1.as_ref().ok_or_else(|| anyhow::anyhow!("File1 path required"))?;
            let sample_base_buf = PathBuf::from(file1_path);
            let (no_ext_sample_base_buf, _) = extension_remover(&sample_base_buf);
            let no_ext_sample_base = no_ext_sample_base_buf.to_string_lossy().into_owned();
            let timestamp = SystemTime::now()
                .duration_since(SystemTime::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .map(|secs| {
                    let dt = chrono::NaiveDateTime::from_timestamp_opt(secs as i64, 0).unwrap();
                    dt.format("%Y%m%d").to_string()
                })
                .unwrap_or_else(|_| "19700101".to_string());
            cwd.join(format!("{}_{}", no_ext_sample_base, timestamp))
        }
    };
    fs::create_dir_all(&out_dir)?;
    Ok(out_dir)
}

/// Searches for a directory for RAM temp files.
/// Prefers /dev/shm (RAM disk) for linux, otherwise returns the standard temp dir.
///
/// # Returns
/// PathBuf: temp dir for RAM files.
fn get_ram_temp_dir() -> PathBuf {
    #[cfg(any(target_os = "linux", target_os = "macos"))]
    {
        if let Ok(metadata) = fs::metadata("/dev/shm") {
            if metadata.is_dir() {
                return PathBuf::from("/dev/shm");
            }
        }
        std::env::temp_dir()
    }

    #[cfg(target_os = "windows")]
    {
        std::env::temp_dir()
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos", target_os = "windows")))]
    {
        std::env::temp_dir()
    }
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

fn compute_base_buffer_size(available_ram: u64, num_cores: usize, data_type: StreamDataType) -> usize {
    let record_size = match data_type {
        StreamDataType::IlluminaFastq => 1_000,
        StreamDataType::OntFastq => 10_000,
        StreamDataType::JustBytes => 500,
    };
    let ram_fraction = 0.25;
    let min_buffer = 10_000;
    let max_buffer = 1_000_000;
    ((available_ram as f64 * ram_fraction / num_cores as f64) / record_size as f64)
        .clamp(min_buffer as f64, max_buffer as f64) as usize
}