/// Analysis of data that does not fit elsewhere (i.e. not FASTX-based)
use std::path::Path;
use std::process::Command;
use anyhow::{Result, anyhow};



// For the show-coords executable from MUMmer suite
#[derive(Debug)]
pub struct AlignmentMetrics {
    pub total_aligned_length: u64,
    pub largest_alignment: u64,
    pub genome_fraction: f64,
    pub na50: u64,
    pub nga50: u64,
    pub na75: u64,
    pub nga75: u64,
    pub la50: u32,
    pub lga50: u32,
    pub la75: u32,
    pub lga75: u32,
    pub mismatches_per_100kbp: f64,
    pub indels_per_100kbp: f64,
    pub misassemblies: u32,
    pub unaligned_length: u64,
    pub duplication_ratio: f64,
}


/// Helper function to compute Nx (e.g., N50, N75)
///
/// # Arguments
///
/// - `lengths`: Lengths series from show-coords or alignment.delta
/// - `fraction`: The threashold i.e. the '50' in N50.
///
/// # Returns
///
/// An Nx (N50, N75) unsigned float.
pub fn compute_nx(lengths: &[u64], fraction: f64) -> u64 {
    let total: u64 = lengths.iter().sum();
    let target = (total as f64 * fraction).ceil() as u64;
    let mut cumsum = 0;
    for &len in lengths {
        cumsum += len;
        if cumsum >= target {
            return len;
        }
    }
    0
}

/// Helper function to compute Lx (e.g., L50, L75)
///
/// # Arguments
///
/// - `lengths`: Lengths series from show-coords or alignment.delta
/// - `fraction`: The threashold i.e. the '50' in N50.
///
/// # Returns
///
/// An Lx (L50, L75) unsigned float.
pub fn compute_lx(lengths: &[u64], fraction: f64) -> u32 {
    let total: u64 = lengths.iter().sum();
    let target = (total as f64 * fraction).ceil() as u64;
    let mut cumsum = 0;
    for (i, &len) in lengths.iter().enumerate() {
        cumsum += len;
        if cumsum >= target {
            return (i + 1) as u32;
        }
    }
    0
}
