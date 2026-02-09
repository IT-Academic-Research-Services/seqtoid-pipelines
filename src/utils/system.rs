// src/utils/system.rs: System functions

use std::cmp::min;
use std::fs;
use std::path::PathBuf;
use std::process::Output;
use std::time::Duration;
use std::collections::HashMap;

use log::{self, LevelFilter, debug, info, error, warn};
use sysinfo::{CpuRefreshKind, MemoryRefreshKind, RefreshKind, System};
use tokio::time::sleep;
use tokio::process::Command as TokioCommand;
use anyhow::{anyhow, Result};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_core::{OsRng, RngCore};

use crate::cli::args::{Arguments, Technology};
use crate::config::defs::{RunConfig, StreamDataType};




/// Detects physical cores (not logical) — fallback to lscpu if sysinfo doesn't distinguish
pub fn detect_physical_cores() -> usize {
    let physical = num_cpus::get_physical();

    if physical > 0 {
        info!("Detected {} physical cores", physical);
        physical
    } else {
        let logical = num_cpus::get();
        let fallback = logical / 2;
        warn!(
            "Failed to detect physical cores (num_cpus returned 0); falling back to {} (logical / 2)",
            fallback
        );
        fallback
    }
}

/// Determines number of cores that can be used for CPU based tasks
///
/// # Arguments
///
/// * `` -
///
/// # Returns
///
/// Result<usize, f32> maximum cores, current cpu usage
// In system.rs
pub async fn detect_cores_and_load(args_threads: usize, use_smt: bool) -> Result<(usize, f32)> {
    let refresh_kind = RefreshKind::nothing().with_cpu(Default::default());
    let mut system = System::new_with_specifics(refresh_kind);
    system.refresh_cpu_all();

    let logical_cores = system.cpus().len();
    let physical_cores = detect_physical_cores();

    let max_cores = if use_smt {
        logical_cores
    } else {
        physical_cores
    };

    let load = system.global_cpu_usage() / 100.0;

    let effective_cores = if args_threads > 0 {
        min(args_threads, max_cores)
    } else {
        max_cores
    };

    debug!("Detected {} logical cores ({} physical); CPU load {}%; using {} threads (use_smt={})", logical_cores, physical_cores, load * 100.0, effective_cores, use_smt);

    Ok((effective_cores, load))
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
pub fn compute_stream_threads(cores: usize, cpu_load: f32, args_threads: usize, use_smt: bool) -> usize {
    let base = if use_smt { cores * 2 } else { cores };  
    if cfg!(target_os = "linux") && base > 50 {
        let max_threads = if cpu_load > 50.0 { base } else { base * 2 };
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
pub fn compute_base_buffer_size(total_input_size_bytes: u64) -> usize {
    let input_gb = total_input_size_bytes as f64 / 1_000_000_000.0;

    // Base scaling: ~100k records per GB input
    let mut records_per_channel = (input_gb * 100_000.0) as usize;

    // Hard cap at 1M records/channel (~1 GB for typical reads) — prevents OOM/throttling
    records_per_channel = records_per_channel.min(1_000_000);

    // Minimum 10k records/channel — avoid tiny buffers on small inputs
    records_per_channel = records_per_channel.max(10_000);

    debug!(
        "Computed base buffer size: {} records/channel (input: {:.1} GB, capped at 1M)",
        records_per_channel, input_gb
    );

    records_per_channel
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
                    debug!("NVMe I/O Util: {}% (device: {})", util_line, device);
                    if util_line.parse::<f32>().unwrap_or(0.0) > 80.0 {
                        warn!("WARNING: High I/O util (>80%) - Potential stall cause. Consider striping NVMe or more /dev/shm usage.");
                    }
                } else {
                    error!("Failed to parse iostat output: {}", stdout);
                }
            }
            Err(e) => error!("iostat error: {}", e),
        }
        sleep(Duration::from_secs(5)).await;
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
