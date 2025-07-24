mod pipelines;
mod utils;
mod config;

use std::time::Instant;
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
use crate::config::defs::{RunConfig, StreamDataType};
use crate::cli::args::Technology;
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
    let module = args.module.clone();
    let run_config = RunConfig { cwd: dir, ram_temp_dir, args, thread_pool, maximal_semaphore, base_buffer_size};

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