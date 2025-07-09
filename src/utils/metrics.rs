/// Calculation of metrics from result files.
use std::collections::BTreeMap;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tokio::sync::mpsc;
use crate::utils::stats::{compute_lx, compute_nx};
use crate::utils::streams::ParseOutput;
use tokio::sync::mpsc::Receiver;
use crate::utils::fastx::{sequence_reader, SequenceReader, SequenceRecord};


/* Mummer-based metrics */

//Assembly Metrics
#[derive(Debug)]
pub struct AssemblyMetrics {
    pub num_contigs: usize,                    // Total number of contigs
    pub contig_counts: BTreeMap<usize, usize>, // Threshold -> Count of contigs >= threshold
    pub total_lengths: BTreeMap<usize, usize>, // Threshold -> Total length of contigs >= threshold
    pub total_length: usize,                   // Total length of all contigs
    pub largest_contig: usize,                 // Length of the largest contig
    pub n50: usize,                            // N50 statistic
    pub n75: usize,                            // N75 statistic
    pub l50: usize,                            // L50 statistic
    pub l75: usize,                            // L75 statistic
    pub gc_percent: f64,                       // GC percentage
}


/// Computes the QUAST-like metrics for a FASTA file
///
/// # Arguments
/// * `rx` - Receiver of SeqeunceRecord::FAsta
/// * `buffer_size` - Size of the output channel buffer.
///
/// # Returns
/// Result<AssemblyMetrics>
pub async fn compute_assembly_metrics(
    mut rx: Receiver<SequenceRecord>,
    thresholds: &[usize],
) -> Result<AssemblyMetrics> {
    let mut contig_lengths = Vec::new();
    let mut total_gc = 0;
    let mut total_bases = 0;

    // Collect lengths and GC content from the stream
    while let Some(record) = rx.recv().await {
        match record {
            SequenceRecord::Fasta { seq, .. } => {
                let len = seq.len();
                contig_lengths.push(len);
                total_bases += len;
                total_gc += seq.iter().filter(|&&b| b.to_ascii_uppercase() == b'G' || b.to_ascii_uppercase() == b'C').count();
            }
            SequenceRecord::Fastq { .. } => {
                return Err(anyhow!("Expected FASTA record, got FASTQ"));
            }
        }
    }

    // Compute metrics
    let num_contigs = contig_lengths.len();
    let total_length: usize = contig_lengths.iter().sum();
    let largest_contig = contig_lengths.iter().max().copied().unwrap_or(0);
    let gc_percent = if total_bases > 0 {
        (total_gc as f64 / total_bases as f64) * 100.0
    } else {
        0.0
    };

    // Compute counts and total lengths per threshold
    let mut contig_counts = BTreeMap::new();
    let mut total_lengths = BTreeMap::new();
    for &threshold in thresholds {
        let filtered_lengths: Vec<_> = contig_lengths.iter().filter(|&&len| len >= threshold).cloned().collect();
        contig_counts.insert(threshold, filtered_lengths.len());
        total_lengths.insert(threshold, filtered_lengths.iter().sum());
    }


    contig_lengths.sort_by(|a, b| b.cmp(a));
    let lengths_u64: Vec<u64> = contig_lengths.iter().map(|&len| len as u64).collect();

    let n50 = compute_nx(&lengths_u64, 0.5) as usize;
    let n75 = compute_nx(&lengths_u64, 0.75) as usize;
    let l50 = compute_lx(&lengths_u64, 0.5) as usize;
    let l75 = compute_lx(&lengths_u64, 0.75) as usize;

    Ok(AssemblyMetrics {
        num_contigs,
        contig_counts,
        total_lengths,
        total_length,
        largest_contig,
        n50,
        n75,
        l50,
        l75,
        gc_percent,
    })
}


//Reference Metrics

/// Computes the QUAST-like metrics for a reference FASTA file
///
/// # Arguments
/// * `path` - PAth to FASTA file. Internally, this is likely to be a RAM file in /dev/shm, but can be any file
///
/// # Returns
/// Result<(total_length, gc_percent)>
pub async fn compute_reference_metrics(path: &PathBuf) -> Result<(usize, f64)> {
    let mut reader = match sequence_reader(path)? {
        SequenceReader::Fasta(reader) => reader,
        _ => return Err(anyhow!("Expected FASTA file for reference")),
    };
    let mut total_length = 0;
    let mut total_gc = 0;
    for result in reader.into_records() {
        let record = result.map_err(|e| anyhow!("Error reading reference FASTA: {}", e))?;
        let seq = record.seq;
        total_length += seq.len();
        total_gc += seq.iter().filter(|&&b| b.to_ascii_uppercase() == b'G' || b.to_ascii_uppercase() == b'C').count();
    }
    let gc_percent = if total_length > 0 {
        (total_gc as f64 / total_length as f64) * 100.0
    } else {
        0.0
    };
    Ok((total_length, gc_percent))
}



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