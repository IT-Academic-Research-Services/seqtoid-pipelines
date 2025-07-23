mod pipelines;
mod utils;
mod config;

use std::time::Instant;
use std::{env, fs};
use std::cmp::min;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::Result;
use clap::Parser;
use num_cpus;
use rayon::ThreadPool;
use tokio::sync::Semaphore;
use tokio::runtime::Builder;
use rayon::ThreadPoolBuilder;

use crate::cli::parse;
use crate::config::defs::RunConfig;
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

    let module = args.module.clone();
    let run_config = RunConfig { cwd: dir, ram_temp_dir, args, thread_pool, maximal_semaphore};

    if let Err(e) = match module.as_str() {
        "consensus_genome" => consensus_genome_run(&run_config).await,
        "create_db" => create_db_run(&run_config).await,
        _ => Err(anyhow::anyhow!("Invalid module: {}", module)),
    } {
        eprintln!("Pipeline failed: {}", e);
        std::process::exit(1);
    }

    println!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
    Ok(())
}

async fn consensus_genome_run(run_config: &RunConfig) -> Result<()> {
    consensus_genome::run(run_config).await
}
async fn create_db_run(run_config: &RunConfig) -> Result<()> {
    db::create_db(run_config).await
}


//*****************
// Setup functions

/// Searches for a directory for RAM temp files.
/// Prefers /dev/shm (RAM disk) for linux, otherwise returns the standard temp dir.
/// # Arguments
///
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