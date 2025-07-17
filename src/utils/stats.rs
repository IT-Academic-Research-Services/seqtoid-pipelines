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
/// - `fraction`: The threshold i.e. the '50' in N50.
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


/// Parses samtools depth lines.
///
/// # Arguments
///
/// * `rx` - samtools stats `ParseOutput::Bytes`
///
/// # Returns
///
/// Hamshmap <String, u32>
pub async fn parse_samtools_depth(rx: mpsc::Receiver<ParseOutput>) -> Result<HashMap<String, u32>> {
    let mut depth_map = HashMap::new();
    let mut stream = ReceiverStream::new(rx);
    let mut line_count = 0;

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(line_bytes) => {
                let line = String::from_utf8_lossy(&line_bytes).trim().to_string();
                if line.is_empty() {
                    continue;
                }

                line_count += 1;
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() != 3 {
                    return Err(anyhow!(
                        "Invalid samtools depth format at line {}: expected 3 fields, found {}",
                        line_count,
                        parts.len()
                    ));
                }

                let chr = parts[0].to_string();
                let pos: u64 = parts[1]
                    .parse()
                    .map_err(|e| anyhow!("Failed to parse position at line {}: {}", line_count, e))?;
                let depth: u32 = parts[2]
                    .parse()
                    .map_err(|e| anyhow!("Failed to parse depth at line {}: {}", line_count, e))?;

                let key = format!("{}:{}", chr, pos);
                depth_map.insert(key, depth);
            }
            _ => {
                return Err(anyhow!("Unexpected non-byte data in samtools depth stream at line {}", line_count));
            }
        }
    }

    if depth_map.is_empty() {
        return Err(anyhow!("No valid samtools depth data parsed"));
    }

    Ok(depth_map)
}


/// The is meant to be used by the results of parse_samtools_depth
///
/// # Arguments
///
/// * `depth_map` - Hamshmap <String, u32>`
///
/// # Returns
///
/// Hashmap <String, f64>
pub fn compute_depth_stats(depth_map: &HashMap<String, u32>) -> Result<HashMap<String, f64>> {
    if depth_map.is_empty() {
        return Err(anyhow!("InsufficientReadsError: There was insufficient coverage so a consensus genome could not be created."));
    }

    let mut depths: Vec<u32> = depth_map.values().copied().collect();
    let count = depths.len() as f64;

    let sum: u64 = depths.iter().map(|&x| x as u64).sum();
    let depth_avg = sum as f64 / count;

    // Sort depths for quantiles
    depths.sort_unstable();

    // Compute quantiles
    let depth_q25 = if count > 0.0 {
        let idx = ((count * 0.25).ceil() - 1.0) as usize;
        depths.get(idx).copied().unwrap_or(0) as f64
    } else {
        0.0
    };

    let depth_q5 = if count > 0.0 {
        let idx = ((count * 0.5).ceil() - 1.0) as usize;
        depths.get(idx).copied().unwrap_or(0) as f64
    } else {
        0.0
    };

    let depth_q75 = if count > 0.0 {
        let idx = ((count * 0.75).ceil() - 1.0) as usize;
        depths.get(idx).copied().unwrap_or(0) as f64
    } else {
        0.0
    };

    let depth_frac_above_10x = depths.iter().filter(|&&d| d >= 10).count() as f64 / count;
    let depth_frac_above_25x = depths.iter().filter(|&&d| d >= 25).count() as f64 / count;
    let depth_frac_above_50x = depths.iter().filter(|&&d| d >= 50).count() as f64 / count;
    let depth_frac_above_100x = depths.iter().filter(|&&d| d >= 100).count() as f64 / count;

    let mut stats = HashMap::new();
    stats.insert("depth_avg".to_string(), depth_avg);
    stats.insert("depth_q.25".to_string(), depth_q25);
    stats.insert("depth_q.5".to_string(), depth_q5);
    stats.insert("depth_q.75".to_string(), depth_q75);
    stats.insert("depth_frac_above_10x".to_string(), depth_frac_above_10x);
    stats.insert("depth_frac_above_25x".to_string(), depth_frac_above_25x);
    stats.insert("depth_frac_above_50x".to_string(), depth_frac_above_50x);
    stats.insert("depth_frac_above_100x".to_string(), depth_frac_above_100x);

    Ok(stats)
}


/// Pars seqkit stats output
///
/// # Arguments
///
/// * `rx` - A Receiver stream of ParseOutput::Lines containing seqkit stats output
///
/// # Returns
///
/// HashMap<String, String>
pub async fn parse_seqkit_stats(rx: mpsc::Receiver<ParseOutput>) -> Result<HashMap<String, String>> {
    let mut stats = HashMap::new();
    let mut stream = ReceiverStream::new(rx);
    let mut header_processed = false;

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(line_bytes) => {
                let line = String::from_utf8_lossy(&line_bytes).trim().to_string();
                if line.is_empty() {
                    continue;
                }

                let parts: Vec<&str> = line.split_whitespace().collect();

                if !header_processed {
                    // Skip header line
                    header_processed = true;
                    continue;
                }

                if parts.len() < 8 {
                    return Err(anyhow!("Invalid seqkit stats format: expected at least 8 fields, found {}", parts.len()));
                }

                // Extract fields
                stats.insert("file".to_string(), parts[0].to_string());
                stats.insert("format".to_string(), parts[1].to_string());
                stats.insert("type".to_string(), parts[2].to_string());
                stats.insert("num_seqs".to_string(), parts[3].replace(",", "").to_string());
                stats.insert("sum_len".to_string(), parts[4].replace(",", "").to_string());
                stats.insert("min_len".to_string(), parts[5].replace(",", "").to_string());
                stats.insert("avg_len".to_string(), parts[6].to_string());
                stats.insert("max_len".to_string(), parts[7].replace(",", "").to_string());
            }
            _ => {
                return Err(anyhow!("Unexpected non-byte data in seqkit stats stream"));
            }
        }
    }

    if stats.is_empty() {
        return Err(anyhow!("No valid seqkit stats data parsed"));
    }

    Ok(stats)
}


/// Parse  stats from ERCC stream
///
/// # Arguments
///
/// * `rx` - A Receiver stream of ParseOutput::Lines containing seqkit stats output
///
/// # Returns
///
/// HashMap<String, u64>
pub async fn parse_ercc_stats(rx: mpsc::Receiver<ParseOutput>) -> Result<HashMap<String, u64>> {
    let mut stats = HashMap::new();
    let mut stream = ReceiverStream::new(rx);
    let mut line_count = 0;

    while let Some(item) = stream.next().await {
        line_count += 1;
        match item {
            ParseOutput::Bytes(line_bytes) => {
                let line = String::from_utf8_lossy(&line_bytes).trim().to_string();
                if line.is_empty() {
                    continue;
                }

                if line.starts_with("SN") {
                    let parts: Vec<&str> = line.split('\t').collect();
                    if parts.len() >= 3 {
                        let key = parts[1].trim_end_matches(':').to_string();
                        let value_str = parts[2];
                        if let Ok(value) = value_str.parse::<u64>() {
                            match key.as_str() {
                                "reads mapped" => {
                                    stats.insert("ercc_mapped_reads".to_string(), value);
                                }
                                "reads mapped and paired" => {
                                    stats.insert("ercc_mapped_paired".to_string(), value);
                                }
                                _ => {}
                            }
                        } else {
                            eprintln!(
                                "Warning: Failed to parse value '{}' for key '{}' at line {}",
                                value_str, key, line_count
                            );
                        }
                    } else {
                        eprintln!(
                            "Warning: Invalid line format at line {}: {}",
                            line_count, line
                        );
                    }
                }
            }
            _ => {
                return Err(anyhow!(
                    "Unexpected non-byte data in ERCC stats stream at line {}",
                    line_count
                ));
            }
        }
    }

    if stats.is_empty() {
        eprintln!("Warning: No valid ERCC stats data parsed");
    }

    Ok(stats)
}