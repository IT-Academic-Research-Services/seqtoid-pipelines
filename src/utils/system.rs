// src/utils/system.rs: System functions

use std::cmp::min;
use std::fs;
use std::path::PathBuf;
use std::process::Output;
use std::time::Duration;

use sysinfo::{CpuRefreshKind, MemoryRefreshKind, RefreshKind, System};
use tokio::time::sleep;
use tokio::process::Command as TokioCommand;
use anyhow::{anyhow, Result};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_core::{OsRng, RngCore};

use crate::cli::args::{Arguments, Technology};
use crate::config::defs::{RunConfig, StreamDataType};


/// Determines number of cores that can be used for CPU based tasks
///
/// # Arguments
///
/// * `` -
///
/// # Returns
///
/// Result<usize, f32> maximum cores, current cpu usage
pub async fn detect_cores_and_load(args_threads: usize) -> Result<(usize, f32)> {
    let refresh_kind = RefreshKind::nothing().with_cpu(Default::default());
    let mut system = System::new_with_specifics(refresh_kind);
    system.refresh_cpu_all();
    let physical_cores = System::physical_core_count().unwrap_or(1);
    system.refresh_cpu_specifics(CpuRefreshKind::nothing().with_cpu_usage());
    sleep(Duration::from_millis(100)).await;
    let cpu_load = system.global_cpu_usage();
    let max_cores = physical_cores.min(args_threads);
    Ok((max_cores, cpu_load))
}


/// Computes the number of stream threads based on cores, load, and OS.
///
/// # Arguments
///
/// * `physical_cores` - Number of real cores on the system
/// * `cpu_load` - Estimate of load on CPU from detect_cores_and_load
///
/// # Returns
///
/// Result<usize, f32> maximum cores, current cpu usage
pub fn compute_stream_threads(physical_cores: usize, cpu_load: f32, args_threads: usize) -> usize {
    if cfg!(target_os = "linux") && physical_cores > 50 {
        let max_threads = if cpu_load > 50.0 { physical_cores } else { physical_cores * 2 };
        max_threads.min(args_threads).min(256)
    } else {
        20 // Lower bound for small systems
    }
}


/// Finds the amount of total and avai;lable RAM, keyed to OS
///
/// # Arguments
///
///
/// # Returns
///
/// Result<u64, u64> total ram, available ram
pub fn detect_ram() -> Result<(u64, u64)> {
    let (total_ram, available_ram) = if cfg!(target_os = "linux") {
        let mut system = System::new_all();
        system.refresh_memory();
        (system.total_memory(), system.available_memory())
    } else if cfg!(target_os = "macos") {
        let refresh_kind = RefreshKind::nothing().with_memory(Default::default());
        let mut system = System::new_with_specifics(refresh_kind);
        system.refresh_memory_specifics(MemoryRefreshKind::everything());
        let total = system.total_memory();
        let used = system.used_memory();
        (total, total.saturating_sub(used))
    } else {
        // Fallback for other OS
        let mut system = System::new_all();
        system.refresh_memory();
        let avail = system.available_memory();
        (system.total_memory(), avail)
    };

    if total_ram == 0 || available_ram == 0 {
        return Err(anyhow!("Failed to detect valid RAM values"));
    }

    Ok((total_ram, available_ram))
}


/// Creates a project-wide RNG from the system, using entropy pool. Optional seed for
/// reproducibility.
///
/// # Arguments
///
///  * `seed` - Seed number that allows reproducible results.
///
/// # Returns
///
/// A StdRng
pub fn generate_rng(seed: Option<u64>) -> StdRng {
    let seed = seed.unwrap_or_else(|| {
        let mut bytes = [0u8; 8];
        OsRng.fill_bytes(&mut bytes);
        u64::from_le_bytes(bytes)
    });
    StdRng::seed_from_u64(seed)
}



/// Computes the size of the RAM buffers used by system parameters.
///
///
/// # Arguments
///
///  * `available_ram` - Amoutn of usable RAM,
/// * `total_ram` - All physical ram on system
/// * `data_type` - StreamDataType (short or long read)
/// * `stream_threads` - Max Number of stream threads allowed
///
/// # Returns
///
/// Usize of basic buffer size per stream
pub fn compute_base_buffer_size(available_ram: u64, total_ram: u64, data_type: StreamDataType, stream_threads: usize) -> usize {
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


/// Searches for a directory for RAM temp files.
/// Prefers /dev/shm (RAM disk) for linux, otherwise returns the standard temp dir.
///
/// # Returns
/// PathBuf: temp dir for RAM files.
pub fn get_ram_temp_dir() -> PathBuf {
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



/// The iostat functions below are currently unused and just here for future use.
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
