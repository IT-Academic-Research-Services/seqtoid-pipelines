/// Analysis of data that does not fit elsewhere (i.e. not FASTX-based)
use std::path::Path;
use std::process::Command;
use anyhow::{Result, anyhow};
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;
use std::collections::HashMap;
use crate::utils::streams::ParseOutput;





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



/// Parses samtools stats lines.
///
/// # Arguments
///
/// * `rx` - samtools stats `ParseOutput::Bytes`
///
/// # Returns
///
/// Hamshmap <String, String>
pub async fn parse_samtools_stats(rx: mpsc::Receiver<ParseOutput>) -> Result<HashMap<String, String>> {
    let mut stats = HashMap::new();
    let mut stream = ReceiverStream::new(rx);

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(line_bytes) => {
                let line = String::from_utf8_lossy(&line_bytes).trim().to_string();
                if line.is_empty() {
                    continue;
                }

                // Parse lines starting with "SN" (Summary Numbers)
                if line.starts_with("SN") {
                    let parts: Vec<&str> = line.split('\t').collect();
                    if parts.len() >= 3 {
                        let key = parts[1].trim_end_matches(':').to_string();
                        let value = parts[2].to_string();
                        stats.insert(key, value);
                    }
                }
                // other sections (e.g., FFQ, COV) if needfed in future
            }
            _ => {
                return Err(anyhow!("Unexpected non-byte data in samtools stats stream"));
            }
        }
    }

    if stats.is_empty() {
        return Err(anyhow!("No valid samtools stats data parsed"));
    }

    Ok(stats)
}
