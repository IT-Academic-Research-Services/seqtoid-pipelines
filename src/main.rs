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

use anyhow::Result;
use clap::Parser;
use rayon::ThreadPool;
use tokio::sync::Semaphore;
use rayon::ThreadPoolBuilder;
use tokio::process::Command as TokioCommand;
use tokio::time::{sleep, Duration};
use std::process::Output;

use crate::cli::parse;
use crate::config::defs::{RunConfig, StreamDataType, PipelineError};
use crate::cli::args::Technology;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::r1r2_base;
use pipelines::consensus_genome;
use pipelines::short_read_mngs;


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

    let refresh_kind = RefreshKind::nothing().with_cpu(Default::default());
    let mut system = System::new_with_specifics(refresh_kind);
    system.refresh_cpu_all();
    let physical_cores = System::physical_core_count().unwrap_or(1);
    system.refresh_cpu_specifics(CpuRefreshKind::nothing().with_cpu_usage());
    tokio::time::sleep(std::time::Duration::from_millis(100)).await; // Allow usage to settle, bake until a golden brown
    let cpu_load = system.global_cpu_usage();

    // Dynamic thread scaling: I/O-bound -> 2x cores; cap for safety
    let stream_threads = if cfg!(target_os = "linux") && physical_cores > 50 {
        let max_threads = if cpu_load > 50.0 { physical_cores } else { physical_cores * 2 };
        max_threads.min(args.threads).min(256) // Cap for safety
    } else {
        20 // msmall systems
    };

    let max_cores = physical_cores.min(args.threads);
    eprintln!("Detected {} physical cores; CPU load {}%; using {} threads for pool, {} for streams",
              physical_cores, cpu_load, max_cores, stream_threads);

    // let max_cores = min(num_cpus::get(), args.threads);

    // Thread pool: NUMA NUMA hey
    let thread_pool = Arc::new(ThreadPoolBuilder::new()
        .num_threads(max_cores)
        .build()
        .expect("Failed to create thread pool"));

    let maximal_semaphore = Arc::new(Semaphore::new(2));

    let total_ram: u64;
    let available_ram: u64;

    #[cfg(target_os = "linux")]
    {
        let mut system = System::new_all();
        system.refresh_memory();
        available_ram = system.available_memory();
        total_ram = system.total_memory();
    }

    #[cfg(target_os = "macos")]
    {
        let refresh_kind = RefreshKind::nothing().with_memory(Default::default());
        let mut system = System::new_with_specifics(refresh_kind);
        system.refresh_memory_specifics(MemoryRefreshKind::everything());
        total_ram = system.total_memory();
        let used_ram = system.used_memory();
        available_ram = total_ram.saturating_sub(used_ram);
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        // Fallback: unknown platform — be conservative, maybe someone's trying to run this on Solaris or something....
        let mut system = System::new_all();
        system.refresh_memory();
        available_ram = system.available_memory();
        total_ram = system.total_memory();
        eprintln!("Unknown OS: using sysinfo defaults");
    }

    eprintln!("Available RAM: {} bytes (~{} GiB)", available_ram, available_ram / 1_073_741_824);
    eprintln!("Total RAM: {} bytes (~{} GiB)", total_ram, total_ram / 1_073_741_824);

    let input_size_mb = get_input_size_mb(&args.file1, &args.file2).unwrap_or(0);
    eprintln!("Total input file size: {} MB", input_size_mb);

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
        input_size_mb
    });

    // let io_monitor_handle = tokio::spawn(monitor_io_utilization("nvme0n1".to_string()));  // Replace with your NVMe device name, e.g., "nvme0n1"

    if let Err(e) = match module.as_str() {
        "consensus_genome" => consensus_genome_run(run_config).await,
        "short_read_mngs" => short_read_mngs_run(run_config).await,
        _ => Err(PipelineError::InvalidConfig(format!("Invalid module: {}", module))),
    } {
        eprintln!("Pipeline failed: {} at {} milliseconds.", e, run_start.elapsed().as_millis());
        std::process::exit(1);
    }

    // Await monitor task at end (optional, for clean shutdown)
    // io_monitor_handle.await.unwrap_or_else(|e| eprintln!("I/O monitor failed: {}", e));

    println!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
    Ok(())
}

async fn monitor_io_utilization(device: String) {
    loop {
        match run_iostat(&device).await {
            Ok(output) => {
                let stdout = String::from_utf8_lossy(&output.stdout);
                if let Some(util_line) = parse_iostat_util(&stdout, &device) {
                    eprintln!("NVMe I/O Util: {}% (device: {})", util_line, device);
                    if util_line.parse::<f32>().unwrap_or(0.0) > 80.0 {
                        eprintln!("WARNING: High I/O util (>80%) - Potential stall cause. Consider striping NVMe or more /dev/shm usage.");
                    }
                } else {
                    eprintln!("Failed to parse iostat output: {}", stdout);
                }
            }
            Err(e) => eprintln!("iostat error: {}", e),
        }
        sleep(Duration::from_secs(5)).await;  // Check every 5s; adjust for granularity
    }
}

// Helper: Run iostat -x -d <device> 1 1 (once, extended disk stats)
async fn run_iostat(device: &str) -> Result<Output> {
    TokioCommand::new("iostat")
        .args(&["-x", "-d", device, "1", "1"])  // Extended, device-specific, 1s interval, count=1
        .output()
        .await
        .map_err(|e| anyhow::anyhow!("Failed to run iostat: {}", e))
}

// Helper: Parse %util from iostat output (e.g., look for line with device and extract last column)
fn parse_iostat_util(output: &str, device: &str) -> Option<String> {
    for line in output.lines() {
        if line.contains(device) && !line.contains("Device") {  // Skip header
            let parts: Vec<&str> = line.split_whitespace().collect();
            if let Some(util) = parts.last() {
                return Some(util.to_string());
            }
        }
    }
    None
}



async fn consensus_genome_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    consensus_genome::run(run_config).await
}

async fn short_read_mngs_run(run_config: Arc<RunConfig>) -> Result<(), PipelineError> {
    short_read_mngs::run(run_config).await
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
        None => {  // Always places out dir in subdir of cwd regardless of location of file 1
            let file1_path: PathBuf = match &args.file1 {
                Some(file) => {
                    let file1_full_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd), None, None, "");
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
                    eprintln!("No R1 tag found. Using bare file 1 stem as sample_base.");
                    file1_r1r2.file_name.unwrap()
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
            cwd.join(format!("{}_{}", dir_base.clone(), timestamp))
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

fn compute_base_buffer_size(available_ram: u64, total_ram: u64, data_type: StreamDataType, stream_threads: usize) -> usize {
    let record_size = match data_type {
        StreamDataType::IlluminaFastq => 1_000,    // ~1KB per FASTQ record
        StreamDataType::OntFastq => 10_000,        // ~10KB for ONT
        StreamDataType::JustBytes => 500,          // Generic
    };

    let ram_fraction = if cfg!(target_os = "linux") { 0.5 } else { 0.3 };
    let min_buffer = 100_000;  // ~100 MB Illumina, ~1 GB ONT
    let max_buffer = 100_000_000;  // ~100 GB Illumina, ~1 TB ONT

    // linux trust available unless <10 GiB; macOS: fallback if <20% total
    let effective_ram = if cfg!(target_os = "linux") && available_ram < 10 * 1_073_741_824 {
        println!(
            "Warning: Critically low available RAM ({} GiB); using 50% total RAM",
            available_ram / 1_073_741_824
        );
        total_ram / 2
    } else if cfg!(target_os = "macos") && available_ram < total_ram / 5 {
        println!(
            "Warning: Low available RAM ({} GiB) vs total ({} GiB); using total RAM",
            available_ram / 1_073_741_824,
            total_ram / 1_073_741_824
        );
        total_ram
    } else {
        available_ram
    };

    let total_buffer_bytes = (effective_ram as f64 * ram_fraction) as usize;

    // OOM guard: Cap if buffers exceed input size or minimum threshold
    let low_ram_threshold = (min_buffer * record_size * stream_threads) as usize;
    if total_buffer_bytes < low_ram_threshold {
        println!(
            "Low RAM ({} GiB); capping at minimal buffer: {} records",
            effective_ram / 1_073_741_824,
            min_buffer
        );
        return min_buffer;
    }

    let max_records = (total_buffer_bytes / record_size) / stream_threads.max(1);
    let buffer_size = max_records.clamp(min_buffer, max_buffer);

    println!(
        "Computed base buffer size: {} records (~{} MB/channel, ~{} MB total est. for {} streams)",
        buffer_size,
        (buffer_size * record_size) / 1_048_576,
        (buffer_size * record_size * stream_threads) / 1_048_576,
        stream_threads
    );
    buffer_size
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