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
use crate::config::defs::{RunConfig, StreamDataType, PipelineError, NRAlignmentBackend, GpuDetection, GpuInfo};
use crate::cli::args::Technology;
use crate::utils::file::{file_path_manipulator, resolve_existing_input_path, derive_sample_base_from_file1};
use crate::utils::fastx::r1r2_base;
use crate::utils::system::{detect_cores_and_load, compute_stream_threads, detect_ram, generate_rng,
                           compute_base_buffer_size, get_ram_temp_dir, detect_gpus, detect_physical_cores,
                           detect_simd_level};
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

    let physical_cores = detect_physical_cores();
    let (max_cores, cpu_load) = detect_cores_and_load(args.threads, args.use_smt).await?;
    let stream_threads = compute_stream_threads(max_cores, cpu_load, args.threads, args.use_smt);

    let thread_pool = Arc::new(ThreadPoolBuilder::new()
        .num_threads(max_cores)
        .build()
        .expect("Failed to create thread pool"));

    let maximal_semaphore = Arc::new(Semaphore::new(2));

    let mut has_gpu = false;
    let mut use_mmseqs_gpu = false;
    let gpu_info = detect_gpus().unwrap_or(GpuDetection { count: 0, gpus: vec![] });

    if gpu_info.count > 0 {
        has_gpu = true;
        for gpu in &gpu_info.gpus {
            info!(
            "GPU {}: {} ({}{} MiB)",
            gpu.index,
            gpu.name,
            gpu.memory_mib.map_or("?".to_string(), |m| m.to_string()),
            if gpu.memory_mib.is_some() { "" } else { "unknown" }
        );
        }

        use_mmseqs_gpu = gpu_info.gpus.iter().any(|g| {
            g.name.contains("H100") || g.name.contains("L40") || g.name.contains("A10")
        });
    }

    let backend = if args.use_diamond {
        NRAlignmentBackend::Diamond
    } else if use_mmseqs_gpu {
        NRAlignmentBackend::MmseqsGpu
    } else {
        NRAlignmentBackend::MmseqsCpu
    };

    let (total_ram, available_ram) = detect_ram()?;
    debug!("Available RAM: {} bytes (~{} GiB)", available_ram, available_ram / 1_073_741_824);
    debug!("Total RAM: {} bytes (~{} GiB)", total_ram, total_ram / 1_073_741_824);

    let input_size = get_input_size(&args.file1, &args.file2).unwrap_or(0);
    debug!("Total input file size: {} Bytes", input_size);

    let rng = generate_rng(args.seed);

    let sdt = match args.technology {
        Technology::Illumina => StreamDataType::IlluminaFastq,
        Technology::ONT => StreamDataType::OntFastq,
    };
    let base_buffer_size = compute_base_buffer_size(input_size);

    let simd = detect_simd_level();

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
        input_size,
        physical_cores,
        max_cores,
        available_ram,
        rng,
        log_level,
        base_backpressure_pause: 1000, // NB: hardcoded for testing
        simd: simd,
        gpu_info: gpu_info,
        has_gpu: has_gpu,
        alignment_backend: backend,

    });

    if let Err(e) = match module.as_str() {
        "consensus_genome" => consensus_genome_run(run_config).await,
        "short_read_mngs" => short_read_mngs_run(run_config).await,
        "build_taxid_lineages" => build_taxid_lineages_run(run_config).await,
        "build_accession2taxid" => build_accession2taxid_run(run_config).await,
        "build_fasta_offset" => build_fasta_offset_db_run(run_config).await,
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

async fn build_fasta_offset_db_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    db::fasta_offset_db(run_config).await
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
            let file1 = match &args.file1 {
                Some(file) => resolve_existing_input_path(file, cwd)
                    .map_err(|e| anyhow::anyhow!(e.to_string()))?,
                None => return Err(anyhow::anyhow!("File1 path required")),
            };

            let sample_base = derive_sample_base_from_file1(&file1)
                .map_err(|e| anyhow::anyhow!(e.to_string()))?;

            let timestamp = chrono::Local::now().format("%Y%m%d").to_string();
            cwd.join(format!("{}_{}", sample_base.display(), timestamp))
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
fn get_input_size(file1: &Option<String>, file2: &Option<String>) -> Result<u64> {
    let mut total_size = 0u64;
    if let Some(f1) = file1 {
        total_size += fs::metadata(f1)?.len();
    }
    if let Some(f2) = file2 {
        total_size += fs::metadata(f2)?.len();
    }
    Ok(total_size) // Bytes -> MB
}