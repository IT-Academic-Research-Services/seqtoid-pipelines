// src/utils/system.rs: System functions

use std::cmp::min;
use std::fs;
use std::path::PathBuf;
use std::process::Output;
use std::time::Duration;
use std::process::Command;

use log::{self, debug, info, error, warn};
use sysinfo::{MemoryRefreshKind, RefreshKind, System};
use tokio::time::sleep;
use tokio::process::Command as TokioCommand;
use anyhow::{anyhow, Result};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_core::{OsRng, RngCore};

use crate::config::defs::{RunConfig, StreamDataType, SimdLevel, GpuDetection, GpuInfo};



/// Detects physical cores (not logical) — fallback to lscpu if sysinfo doesn't distinguish
pub fn detect_physical_cores() -> usize {
    let physical = num_cpus::get_physical();

    if physical > 0 {
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

/// Single source of truth for buffer sizes — fully dynamic.
pub fn compute_buffer_size(
    config: &RunConfig,
    phase: &str,
    data_type: StreamDataType,
    multiplier: f64,
) -> usize {
    let input_gb = config.input_size as f64 / 1_000_000_000.0;
    let avail_gib = config.available_ram as f64 / (1024.0 * 1024.0 * 1024.0);

    let mut records_per_channel = (input_gb * 100_000.0) as usize;
    records_per_channel = records_per_channel.clamp(10_000, 1_000_000);

    let record_size = match data_type {
        StreamDataType::IlluminaFastq => 1_000,
        StreamDataType::OntFastq => 10_000,
        StreamDataType::JustBytes => 512,
    };

    let mut bytes = (records_per_channel * record_size) as f64;

    // Hardware scaling
    if avail_gib >= 1000.0 {
        bytes *= 2.5;
        bytes = bytes.min(128.0 * 1024.0 * 1024.0);
    } else if avail_gib >= 500.0 {
        bytes *= 2.0;
        bytes = bytes.min(96.0 * 1024.0 * 1024.0);
    } else if avail_gib >= 128.0 {
        bytes *= 1.5;
        bytes = bytes.min(64.0 * 1024.0 * 1024.0);
    } else {
        bytes = bytes.min(16.0 * 1024.0 * 1024.0);
    }

    // Phase hot-path boost
    match phase {
        "stream_to_cmd" | "write_to_fifo" | "fanout" | "deinterleave" => {
            bytes *= 1.8;
            bytes = bytes.min(128.0 * 1024.0 * 1024.0);
        }
        "parse_fastq" | "parse_bytes" => {
            bytes *= 0.5;
            bytes = bytes.max(4.0 * 1024.0 * 1024.0);
        }
        "batch_rayon" => {
            bytes *= 2.0;
            bytes = bytes.min(64.0 * 1024.0 * 1024.0);
        }
        _ => {}
    }

    let final_bytes = (bytes * multiplier)
        .clamp(4.0 * 1024.0 * 1024.0, 128.0 * 1024.0 * 1024.0) as usize;

    debug!(
        "compute_buffer_size({}/{:?}) = {} MiB (input {:.1}GB, avail {:.1}GiB, mult {:.1})",
        phase, data_type, final_bytes / (1024*1024), input_gb, avail_gib, multiplier
    );

    final_bytes
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
#[allow(dead_code)]
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
#[allow(dead_code)]
async fn run_iostat(device: &str) -> Result<Output> {
    TokioCommand::new("iostat")
        .args(&["-x", "-d", device, "1", "1"])  // Extended, device-specific, 1s interval, count=1
        .output()
        .await
        .map_err(|e| anyhow::anyhow!("Failed to run iostat: {}", e))
}

// Helper: Parse %util from iostat output (e.g., look for line with device and extract last column)
#[allow(dead_code)]
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

pub fn detect_gpus() -> Result<GpuDetection> {
    #[cfg(target_os = "linux")]
    {
        // Try nvidia-smi first (most reliable on your cluster)
        if let Ok(nvidia) = detect_nvidia() {
            if !nvidia.gpus.is_empty() {
                info!("Detected {} NVIDIA GPU(s): {:?}", nvidia.count, nvidia.gpus);
                return Ok(nvidia);
            }
        }

        // Fallback: very basic lspci detection (just count + rough names)
        if let Ok(basic) = detect_lspci_basic() {
            if basic.count > 0 {
                info!("Detected {} GPU(s) via lspci (no nvidia-smi available)", basic.count);
                return Ok(basic);
            }
        }

        debug!("No GPUs detected on Linux");
        return Ok(GpuDetection { count: 0, gpus: vec![] });
    }

    #[cfg(target_os = "macos")]
    {
        if let Ok(macos) = detect_macos() {
            if macos.count > 0 {
                info!("Detected {} GPU(s) on macOS: {:?}", macos.count, macos.gpus);
                return Ok(macos);
            }
        }
        debug!("No GPUs detected on macOS");
        Ok(GpuDetection { count: 0, gpus: vec![] })
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        warn!("GPU detection not implemented for this platform");
        Ok(GpuDetection { count: 0, gpus: vec![] })
    }
}

#[allow(dead_code)]
fn detect_nvidia() -> Result<GpuDetection> {
    // Query name + memory + driver
    let output = Command::new("nvidia-smi")
        .args([
            "--query-gpu=index,name,memory.total,driver_version",
            "--format=csv,noheader,nounits",
        ])
        .output()
        .map_err(|e| anyhow!("nvidia-smi not found or failed: {}", e))?;

    if !output.status.success() {
        debug!("nvidia-smi command failed");
        return Ok(GpuDetection { count: 0, gpus: vec![] });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut gpus = vec![];

    for (i, line) in stdout.lines().enumerate() {
        let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if parts.len() < 4 {
            continue;
        }

        let index = i;
        let name = parts[1].to_string();
        let memory_str = parts[2].trim_end_matches(" MiB");
        let memory_mib = memory_str.parse::<u64>().ok();
        let driver = Some(parts[3].to_string());

        gpus.push(GpuInfo {
            index,
            name,
            memory_mib,
            is_discrete: true,           // NVIDIA cards are always discrete
            driver,
        });
    }

    Ok(GpuDetection {
        count: gpus.len(),
        gpus,
    })
}

fn detect_macos() -> Result<GpuDetection> {
    let output = Command::new("system_profiler")
        .args(["SPDisplaysDataType", "-json"])
        .output()
        .map_err(|e| anyhow!("system_profiler failed: {}", e))?;

    if !output.status.success() {
        return Ok(GpuDetection { count: 0, gpus: vec![] });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Very simple parsing — real parsing would use serde_json, but we avoid extra deps here
    let mut gpus = vec![];
    let mut index = 0;

    for line in stdout.lines() {
        if line.contains("\"Chipset Model\":") || line.contains("\"Device Name\":") {
            let name = line
                .split(':')
                .nth(1)
                .unwrap_or("")
                .trim()
                .trim_matches('"')
                .trim_matches(',')
                .to_string();

            gpus.push(GpuInfo {
                index,
                name: if name.is_empty() { "Apple GPU".to_string() } else { name },
                memory_mib: None,           // macOS doesn't expose easily via CLI
                is_discrete: false,         // Apple Silicon is integrated
                driver: None,
            });
            index += 1;
        }
    }

    // Fallback: assume at least 1 on modern Macs
    if gpus.is_empty() {
        gpus.push(GpuInfo {
            index: 0,
            name: "Apple integrated GPU".to_string(),
            memory_mib: None,
            is_discrete: false,
            driver: None,
        });
    }

    Ok(GpuDetection {
        count: gpus.len(),
        gpus,
    })
}

#[allow(dead_code)]
fn detect_lspci_basic() -> Result<GpuDetection> {
    let output = Command::new("lspci")
        .arg("-mm")           // machine-readable format, easier to parse
        .arg("-v")            // more verbose (includes device class)
        .output()
        .map_err(|e| anyhow!("Failed to run lspci: {}", e))?;

    if !output.status.success() {
        debug!("lspci command failed or not found");
        return Ok(GpuDetection { count: 0, gpus: vec![] });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut gpus = vec![];
    let mut index = 0;

    // Keywords that usually indicate a GPU
    let gpu_indicators = [
        "VGA compatible controller",
        "3D controller",
        "Display controller",
    ];

    // Common GPU vendors (for rough classification)
    let gpu_vendors = [
        "NVIDIA", "AMD", "ATI", "Intel", "Matrox", "Radeon", "GeForce",
        "Quadro", "Tesla", "RTX", "GTX", "FirePro", "Arc",
    ];

    for line in stdout.lines() {
        // Skip lines that are clearly not device descriptions
        if line.contains("Kernel driver") || line.contains("Subsystem") || line.is_empty() {
            continue;
        }

        let lower = line.to_lowercase();

        // Check if this line describes a graphics device
        if gpu_indicators.iter().any(|&kw| lower.contains(&kw.to_lowercase())) {
            // Try to extract something that looks like a device name/model
            let mut name_parts = vec![];

            // Split on ':' and take the part after the slot/class
            if let Some(after_class) = line.splitn(3, ':').nth(2) {
                let cleaned = after_class.trim();
                name_parts.push(cleaned.to_string());
            }

            // Look for vendor/model in the line
            for vendor in &gpu_vendors {
                if lower.contains(&vendor.to_lowercase()) {
                    name_parts.push(vendor.to_string());
                    break;
                }
            }

            let name = if name_parts.is_empty() {
                "Unknown GPU".to_string()
            } else {
                name_parts.join(" ")
            };

            gpus.push(GpuInfo {
                index,
                name,
                memory_mib: None,         // lspci doesn't show VRAM
                is_discrete: !lower.contains("intel") && !lower.contains("integrated"),
                driver: None,
            });

            index += 1;
        }
    }

    debug!("lspci basic detection found {} potential GPUs", gpus.len());

    Ok(GpuDetection {
        count: gpus.len(),
        gpus,
    })
}

/// Computes a safe concurrency level for a pipeline phase.
///
/// # Parameters
/// - `config`: RunConfig with max_cores and RAM info
/// - `phase_name`: For logging (e.g. "minimap2_nt", "call_hits_nr")
/// - `ram_per_thread_gb`: Estimated GB per thread (mm2 ~65, call_hits ~0.5–1)
/// - `cpu_divisor`: How many threads per core (mm2: 6, call_hits: 3–4)
/// - `max_cap`: Hard upper limit (phase-specific)
/// - `min_threads`: Minimum even if RAM/CPU say lower
pub fn compute_phase_concurrency(
    config: &RunConfig,
    phase_name: &str,
    ram_per_thread_gb: f64,
    cpu_divisor: f64,
    max_cap: usize,
    min_threads: usize,
) -> usize {
    let total_threads = config.max_cores as f64;
    let (total_bytes, avail_bytes) = detect_ram().unwrap_or((0, 0)); // fallback to 0 if fails

    let total_gib = total_bytes as f64 / (1024.0 * 1024.0 * 1024.0);
    let avail_gib = avail_bytes as f64 / (1024.0 * 1024.0 * 1024.0);

    let target_ram_gib = avail_gib * 0.80; // 80% safe usage
    let max_jobs_from_ram = (target_ram_gib / ram_per_thread_gb).floor() as usize;

    let max_jobs_from_cpu = (total_threads / cpu_divisor).max(1.0) as usize;

    let mut concurrency = max_jobs_from_ram
        .min(max_jobs_from_cpu)
        .min(max_cap)
        .max(min_threads);

    // Small machine overrides
    if total_gib < 128.0 {
        concurrency = concurrency.min(4);
    }
    if total_gib < 64.0 {
        concurrency = concurrency.min(2);
    }

    info!(
        "{} concurrency: {} jobs (RAM limit: {}, CPU limit: {}, cap: {}, min: {})",
        phase_name, concurrency, max_jobs_from_ram, max_jobs_from_cpu, max_cap, min_threads
    );

    concurrency
}

pub fn compute_batch_size(
    est_total_lines: Option<usize>,  // Optional: If known (e.g., from prior count); else None for defaults
    avg_line_bytes: usize,           // ~200
    target_batch_mb: usize,          // 100–500
    concurrency: usize,
) -> usize {
    let target_bytes = target_batch_mb * 1024 * 1024;
    let base_batch = target_bytes / avg_line_bytes;  // e.g., 100MB / 200B = 500k lines

    if let Some(total_lines) = est_total_lines {
        // Dynamic: ~10–20 batches per concurrent task for overlap
        let total_batches = concurrency * 15;
        let from_total = total_lines / total_batches;
        base_batch.min(from_total)
    } else {
        base_batch  // Fixed if unknown
    }
        .clamp(10_000, 500_000)  // Min/max
}

pub fn detect_simd_level() -> SimdLevel {
    #[cfg(not(target_arch = "x86_64"))]
    {
        info!("Non-x86_64 platform → forcing Scalar SIMD (MacBook Pro is safe here)");
        return SimdLevel::Scalar;
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx512f") &&
            is_x86_feature_detected!("avx512bw") &&
            is_x86_feature_detected!("avx512vl") &&
            is_x86_feature_detected!("avx512dq") {
            info!("AVX-512 (F+BW+VL+DQ) detected → enabling AVX512 path (r8id.48xlarge / EPYC)");
            SimdLevel::Avx512
        } else if is_x86_feature_detected!("avx2") {
            debug!("AVX2 detected (no AVX-512) → falling back to AVX2");
            SimdLevel::Avx2
        } else {
            info!("No vector extensions detected → Scalar fallback (safe on all hardware)");
            SimdLevel::Scalar
        }
    }
}

/// Configures Transparent Huge Pages (THP) safely for large EPYC nodes.
///
/// We deliberately use `madvise` instead of `always`. Forcing `always` globally
/// interacts very badly with GNU sort (when using large `-S` buffers) and some
/// other tools, causing severe memory fragmentation and massive slowdowns on
/// machines with ≥ 512 GiB RAM.
///
/// `madvise` is the safe default: it lets the application (and libraries) opt
/// into huge pages where it actually helps (large allocations, mmseqs, etc.)
/// without punishing tools like `sort`, `samtools sort`, etc.
pub async fn ensure_transparent_hugepages(config: &RunConfig) -> Result<()> {
    if !cfg!(target_os = "linux") || config.max_cores < 64 {
        debug!("THP configuration skipped (not a large Linux node)");
        return Ok(());
    }

    info!("Configuring Transparent Huge Pages (madvise mode) for large EPYC node...");

    // We only touch these two. defrag is left at the system default.
    let paths = [
        "/sys/kernel/mm/transparent_hugepage/enabled",
        "/sys/kernel/mm/transparent_hugepage/shmem_enabled",
    ];

    let target = "madvise";

    for path in &paths {
        let current = match tokio::fs::read_to_string(path).await {
            Ok(s) => s.trim().to_string(),
            Err(e) => {
                debug!("Could not read {}: {}", path, e);
                continue;
            }
        };

        // Already correct?
        if current.contains("[madvise]") || current == target {
            debug!("{} already set to madvise", path);
            continue;
        }

        if let Err(e) = tokio::fs::write(path, target).await {
            warn!("Failed to set {} to '{}': {}", path, target, e);
        } else {
            info!("Set {} → {}", path, target);
        }
    }

    info!("THP configuration complete (mode: madvise)");
    Ok(())
}