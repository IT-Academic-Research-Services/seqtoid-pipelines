/// Calculation of metrics from result files.
use anyhow::{anyhow, Result};
use crate::utils::stats::{compute_lx, compute_nx};
use crate::utils::streams::ParseOutput;
use tokio::sync::mpsc::Receiver;


// MUMmer
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

}

/// Calculation of metrics from the MUMmer show-coords program
///
/// # Arguments
///
/// - `lengths`: Lengths series from show-coords or alignment.delta
/// - `fraction`: The threashold i.e. the '50' in N50.
///
/// # Returns
///
/// An Nx (N50, N75) unsigned float.
pub async fn parse_show_coords(
    mut stdout_stream: Receiver<ParseOutput>,
    reference_length: u64,
    total_contig_length: u64,
) -> Result<AlignmentMetrics> {
    let mut alignments = Vec::new();

    // Process the stream
    while let Some(parse_output) = stdout_stream.recv().await {
        let line = match parse_output {
            ParseOutput::Bytes(bytes) => String::from_utf8(bytes)?,
            _ => continue, // Skip non-byte data
        };

        // Skip header and non-data lines
        if line.starts_with("NUCMER") || line.starts_with('[') || line.trim() == "(END)" {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 11 {
            continue; // Skip malformed lines
        }

        let alignment = (
            fields[0].parse::<u64>().unwrap_or(0), // S1
            fields[1].parse::<u64>().unwrap_or(0), // E1
            fields[3].parse::<u64>().unwrap_or(0), // S2
            fields[4].parse::<u64>().unwrap_or(0), // E2
            fields[6].parse::<u64>().unwrap_or(0), // LEN 2 (query length)
            fields[7].parse::<f64>().unwrap_or(0.0), // % IDY
            fields[9].parse::<f64>().unwrap_or(0.0), // COV R
            fields[10].parse::<f64>().unwrap_or(0.0), // COV Q
        );
        alignments.push(alignment);
    }

    if alignments.is_empty() {
        return Err(anyhow!("No valid alignments found in show-coords output"));
    }

    // Compute metrics
    let total_aligned_length = alignments.iter().map(|a| a.4).sum::<u64>();
    let largest_alignment = alignments.iter().map(|a| a.4).max().unwrap_or(0);
    let genome_fraction = (total_aligned_length as f64 / reference_length as f64) * 100.0;

    // Nx/Lx metrics
    let mut sorted_lengths: Vec<u64> = alignments.iter().map(|a| a.4).collect();
    sorted_lengths.sort_by(|a, b| b.cmp(a));
    let na50 = compute_nx(&sorted_lengths, 0.5);
    let nga50 = compute_nx(&sorted_lengths, 0.5); // Simplified; adjust for reference length if needed
    let na75 = compute_nx(&sorted_lengths, 0.75);
    let nga75 = compute_nx(&sorted_lengths, 0.75);
    let la50 = compute_lx(&sorted_lengths, 0.5);
    let lga50 = compute_lx(&sorted_lengths, 0.5);
    let la75 = compute_lx(&sorted_lengths, 0.75);
    let lga75 = compute_lx(&sorted_lengths, 0.75);


    Ok(AlignmentMetrics {
        total_aligned_length,
        largest_alignment,
        genome_fraction,
        na50,
        nga50,
        na75,
        nga75,
        la50,
        lga50,
        la75,
        lga75,
    })
}