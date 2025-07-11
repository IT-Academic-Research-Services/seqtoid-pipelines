use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader, AsyncWriteExt};
use std::path::Path;
use std::sync::Arc;
use tokio::sync::Mutex;
use crate::utils::stats::{compute_lx, compute_nx};
use crate::utils::streams::ParseOutput;
use crate::utils::fastx::SequenceRecord;
use tokio::sync::mpsc::{Receiver, Sender};
use crate::utils::fastx::{sequence_reader, SequenceReader};
use rust_htslib::{bam, bam::Read, bam::record::Record, bam::Header, bam::HeaderView};



#[derive(Debug)]
pub struct AssemblyMetrics {
    pub num_contigs: usize,
    pub contig_counts: BTreeMap<usize, usize>,
    pub total_lengths: BTreeMap<usize, usize>,
    pub total_length: usize,
    pub largest_contig: usize,
    pub n50: usize,
    pub n75: usize,
    pub l50: usize,
    pub l75: usize,
    pub gc_percent: f64,
    pub unaligned_length: usize,
    pub unaligned_contigs: usize,
    pub n_per_100kbp: f64,
}

#[derive(Debug)]
pub struct ReferenceMetrics {
    pub total_length: usize,
    pub gc_percent: f64,
}

#[derive(Debug)]
pub struct BAMMetrics {
    pub total_reads: usize,
    pub left_reads: usize,
    pub right_reads: usize,
    pub percent_mapped: f64,
    pub properly_paired: f64,
    pub avg_coverage_depth: usize,
    pub coverage_above_1x: f64,
    pub ref_mapped: f64,
    pub ref_properly_paired: f64,
    pub ref_avg_coverage_depth: usize,
    pub ref_coverage_above_1x: f64,
}

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
    pub duplication_ratio: f64,
    pub mismatches_per_100kbp: f64,
    pub indels_per_100kbp: f64,
    pub misassemblies: u32,
    pub misassembled_contigs: usize,
    pub misassembled_contigs_length: usize,
    pub local_misassemblies: u32,
    pub scaffold_gap_ext_mis: u32,
    pub scaffold_gap_loc_mis: u32,
    pub unaligned_mis_contigs: usize,
}

#[derive(Debug)]
pub struct AllMetrics {
    pub assembly: AssemblyMetrics,
    pub reference: ReferenceMetrics,
    pub bam: BAMMetrics,
    pub alignment: AlignmentMetrics,
}

pub async fn compute_assembly_metrics(
    mut rx: Receiver<SequenceRecord>,
    thresholds: &[usize],
    show_coords_stream: Receiver<ParseOutput>,
) -> Result<AssemblyMetrics> {
    let mut contig_lengths = Vec::new();
    let mut total_gc = 0;
    let mut total_bases = 0;
    let mut total_ns = 0;
    let mut record_count = 0;

    eprintln!("Starting compute_assembly_metrics");
    while let Some(record) = rx.recv().await {
        record_count += 1;
        // eprintln!("Received FASTA record {}: {:?}", record_count, record);
        match record {
            SequenceRecord::Fasta { seq, id, .. } => {
                let len = seq.len();
                eprintln!("FASTA record id: {}, length: {}", id, len);
                contig_lengths.push(len);
                total_bases += len;
                total_gc += seq.iter().filter(|&&b| b.to_ascii_uppercase() == b'G' || b.to_ascii_uppercase() == b'C').count();
                total_ns += seq.iter().filter(|&&b| b.to_ascii_uppercase() == b'N').count();
            }
            SequenceRecord::Fastq { .. } => {
                return Err(anyhow!("Expected FASTA record, got FASTQ at record {}", record_count));
            }
        }
    }
    eprintln!("Total FASTA records processed: {}, contig_lengths: {:?}", record_count, contig_lengths);

    if contig_lengths.is_empty() {
        return Err(anyhow!("No valid FASTA records found in consensus_eval_stream"));
    }

    let num_contigs = contig_lengths.len();
    let total_length: usize = contig_lengths.iter().sum();
    let largest_contig = contig_lengths.iter().max().copied().unwrap_or(0);
    let gc_percent = if total_bases > 0 {
        (total_gc as f64 / total_bases as f64) * 100.0
    } else {
        0.0
    };
    let n_per_100kbp = if total_bases > 0 {
        (total_ns as f64 / total_bases as f64) * 100_000.0
    } else {
        0.0
    };

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

    let (unaligned_length, unaligned_contigs) = compute_unaligned_metrics(show_coords_stream, total_length as u64, &contig_lengths).await?;

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
        unaligned_length: unaligned_length as usize,
        unaligned_contigs,
        n_per_100kbp,
    })
}

async fn compute_unaligned_metrics(
    mut stdout_stream: Receiver<ParseOutput>,
    total_contig_length: u64,
    contig_lengths: &[usize],
) -> Result<(u64, usize)> {
    let mut total_aligned_length = 0;
    let mut aligned_contigs = HashSet::new();
    eprintln!("Starting compute_unaligned_metrics");

    while let Some(parse_output) = stdout_stream.recv().await {
        let line = match parse_output {
            ParseOutput::Bytes(bytes) => String::from_utf8(bytes)
                .map_err(|e| anyhow!("Invalid UTF-8: {}", e))?,
            _ => {
                eprintln!("Skipping non-bytes output");
                continue;
            }
        };
        eprintln!("Received line: {}", line);
        if line.is_empty() || line.starts_with("NUCMER") || line.starts_with('[') || line.trim() == "(END)" {
            eprintln!("Skipping header, footer, or empty line: {}", line);
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        eprintln!("Fields (len={}): {:?}", fields.len(), fields);
        // Check for alignment line: must have at least 13 fields, numeric S1, LEN 2, and tags
        if fields.len() < 13 || !fields[0].parse::<u64>().is_ok() || !fields[6].parse::<u64>().is_ok() {
            eprintln!("Skipping non-alignment line: {}", line);
            continue;
        }
        // Find the query contig name (last field)
        if let Some(contig_name) = fields.last() {
            total_aligned_length += fields[6].parse::<u64>().map_err(|e| {
                anyhow!("Failed to parse LEN 2 from '{}': {}", fields[6], e)
            })?;
            aligned_contigs.insert(contig_name.to_string());
            eprintln!("Added query contig: {}", contig_name);
        } else {
            eprintln!("No query contig found in line: {}", line);
            continue;
        }
    }

    if total_aligned_length == 0 {
        return Err(anyhow!("compute_unaligned_metrics: No valid alignments found in show-coords output"));
    }

    if aligned_contigs.len() > contig_lengths.len() {
        eprintln!("Error: aligned_contigs.len() ({}) > contig_lengths.len() ({})",
                  aligned_contigs.len(), contig_lengths.len());
        eprintln!("Aligned contigs: {:?}", aligned_contigs);
        return Err(anyhow!("More aligned contigs than total contigs in assembly"));
    }

    let unaligned_contigs = contig_lengths.len() - aligned_contigs.len();
    let unaligned_length = total_contig_length.saturating_sub(total_aligned_length);

    eprintln!("Total aligned length: {}, Unaligned length: {}, Unaligned contigs: {}",
              total_aligned_length, unaligned_length, unaligned_contigs);

    Ok((unaligned_length, unaligned_contigs))
}

pub async fn compute_duplication_ratio(
    mut stdout_stream: Receiver<ParseOutput>,
    total_contig_length: u64,
) -> Result<f64> {
    let mut aligned_bases = Vec::new();

    while let Some(parse_output) = stdout_stream.recv().await {
        let line = match parse_output {
            ParseOutput::Bytes(bytes) => String::from_utf8(bytes)?,
            _ => continue,
        };
        if line.starts_with("NUCMER") || line.starts_with('[') || line.trim() == "(END)" {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 11 {
            continue;
        }
        let s2 = fields[3].parse::<u64>().unwrap_or(0);
        let e2 = fields[4].parse::<u64>().unwrap_or(0);
        aligned_bases.push((s2, e2));
    }

    if aligned_bases.is_empty() {
        return Err(anyhow!("compute_duplication_ratio: No valid alignments found in show-coords output"));
    }

    let total_aligned = aligned_bases.iter().map(|&(s, e)| e - s + 1).sum::<u64>();
    Ok(total_aligned as f64 / total_contig_length as f64)
}

pub async fn compute_mismatches_per_100kbp(delta_path: &Path) -> Result<f64> {
    let file = File::open(delta_path).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let mut mismatches = 0;
    let mut aligned_bases = 0;

    while let Some(line) = lines.next_line().await? {
        if line.starts_with('>') || line.starts_with("NUCMER") {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 7 {
            aligned_bases += fields[2].parse::<u64>().unwrap_or(0); // LEN 2
            if fields.len() > 7 {
                mismatches += fields.get(7).and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            }
        }
    }

    if aligned_bases == 0 {
        return Err(anyhow!("No aligned bases found in delta file"));
    }

    Ok((mismatches as f64 / aligned_bases as f64) * 100_000.0)
}

pub async fn compute_indels_per_100kbp(delta_path: &Path) -> Result<f64> {
    let file = File::open(delta_path).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let mut indels = 0;
    let mut aligned_bases = 0;

    while let Some(line) = lines.next_line().await? {
        if line.starts_with('>') || line.starts_with("NUCMER") {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 7 {
            aligned_bases += fields[2].parse::<u64>().unwrap_or(0); // LEN 2
            if fields.len() > 8 {
                indels += fields.get(8).and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            }
        }
    }

    if aligned_bases == 0 {
        return Err(anyhow!("No aligned bases found in delta file"));
    }

    Ok((indels as f64 / aligned_bases as f64) * 100_000.0)
}

pub async fn compute_misassemblies(
    delta_path: &Path,
    contig_lengths: &BTreeMap<String, usize>,
) -> Result<(u32, usize, usize, u32, u32, u32, usize)> {
    let file = File::open(delta_path).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let mut misassemblies = 0;
    let mut local_misassemblies = 0;
    let mut misassembled_contigs = HashSet::new();
    let mut prev_end = 0;
    let mut prev_contig = String::new();

    while let Some(line) = lines.next_line().await? {
        if line.starts_with('>') {
            prev_contig = line.split_whitespace().nth(1).unwrap_or("").to_string();
            prev_end = 0;
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 7 {
            let start = fields[3].parse::<u64>().unwrap_or(0); // S2
            let end = fields[4].parse::<u64>().unwrap_or(0); // E2
            let len2 = fields[2].parse::<u64>().unwrap_or(0); // LEN 2
            if !prev_contig.is_empty() {
                if start > prev_end + 1000 {
                    misassemblies += 1;
                    misassembled_contigs.insert(prev_contig.clone());
                } else if start > prev_end + 100 && len2 < 1000 {
                    local_misassemblies += 1;
                    misassembled_contigs.insert(prev_contig.clone());
                }
            }
            prev_end = end;
        }
    }

    let misassembled_contigs_length = misassembled_contigs
        .iter()
        .map(|contig| contig_lengths.get(contig).copied().unwrap_or(0))
        .sum::<usize>();

    Ok((
        misassemblies,
        misassembled_contigs.len(),
        misassembled_contigs_length,
        local_misassemblies,
        0, // scaffold_gap_ext_mis
        0, // scaffold_gap_loc_mis
        0, // unaligned_mis_contigs
    ))
}

pub async fn compute_reference_metrics(path: &PathBuf) -> Result<ReferenceMetrics> {
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
    Ok(ReferenceMetrics { total_length, gc_percent })
}

pub async fn parse_show_coords(
    mut stdout_stream: Receiver<ParseOutput>,
    reference_length: u64,
    total_contig_length: u64,
) -> Result<AlignmentMetrics> {
    let mut alignments = Vec::new();

    while let Some(parse_output) = stdout_stream.recv().await {
        let line = match parse_output {
            ParseOutput::Bytes(bytes) => String::from_utf8(bytes)?,
            _ => continue,
        };
        if line.is_empty() || line.starts_with("NUCMER") || line.starts_with('[') || line.trim() == "(END)" {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 13 || !fields[0].parse::<u64>().is_ok() {
            continue;
        }
        let alignment = (
            fields[0].parse::<u64>().unwrap_or(0), // S1
            fields[1].parse::<u64>().unwrap_or(0), // E1
            fields[3].parse::<u64>().unwrap_or(0), // S2
            fields[4].parse::<u64>().unwrap_or(0), // E2
            fields[6].parse::<u64>().unwrap_or(0), // LEN 2
            fields[7].parse::<f64>().unwrap_or(0.0), // % IDY
            fields[9].parse::<f64>().unwrap_or(0.0), // COV R
            fields[10].parse::<f64>().unwrap_or(0.0), // COV Q
        );
        alignments.push(alignment);
    }



    if alignments.is_empty() {
        return Err(anyhow!("parse_show_coords: No valid alignments found in show-coords output"));
    }

    let total_aligned_length = alignments.iter().map(|a| a.4).sum::<u64>();
    let largest_alignment = alignments.iter().map(|a| a.4).max().unwrap_or(0);
    let genome_fraction = if reference_length > 0 {
        (total_aligned_length as f64 / reference_length as f64) * 100.0
    } else {
        0.0
    };

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
        duplication_ratio: 0.0, // Computed separately
        mismatches_per_100kbp: 0.0, // Computed separately
        indels_per_100kbp: 0.0, // Computed separately
        misassemblies: 0, // Computed separately
        misassembled_contigs: 0, // Computed separately
        misassembled_contigs_length: 0, // Computed separately
        local_misassemblies: 0, // Computed separately
        scaffold_gap_ext_mis: 0, // Computed separately
        scaffold_gap_loc_mis: 0, // Computed separately
        unaligned_mis_contigs: 0, // Computed separately
    })
}


pub async fn compute_bam_metrics(
    mut bam_stream: Receiver<ParseOutput>,
    reference_length: u64,
) -> Result<BAMMetrics> {
    // Create a minimal header (adjust based on your reference if needed)
    let mut header = Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &"ref")
            .push_tag(b"LN", &(reference_length as i32)),
    );
    let header_view = HeaderView::from_header(&header);

    let mut total_reads = 0;
    let mut left_reads = 0;
    let mut right_reads = 0;
    let mut mapped_reads = 0;
    let mut properly_paired = 0;
    let mut total_depth = 0;
    let mut bases_covered = 0;
    let mut pileup_positions = Vec::new();

    while let Some(parse_output) = bam_stream.recv().await {
        let record_bytes = match parse_output {
            ParseOutput::Bytes(bytes) => bytes,
            _ => continue,
        };

        let sam_line = String::from_utf8(record_bytes)
            .map_err(|e| anyhow!("Invalid UTF-8 in SAM line: {}", e))?;
        let record = Record::from_sam(&header_view, sam_line.as_bytes())
            .map_err(|e| anyhow!("Failed to parse SAM record: {}", e))?;

        total_reads += 1;
        if record.is_first_in_template() {
            left_reads += 1;
        } else if record.is_last_in_template() {
            right_reads += 1;
        }
        if !record.is_unmapped() {
            mapped_reads += 1;
            let start = record.pos() as u64;
            let end_i64 = record.cigar().end_pos();
            if end_i64 < 0 {
                return Err(anyhow!("Invalid end position for record: {}", end_i64));
            }
            let end = end_i64 as u64;
            if end > reference_length {
                return Err(anyhow!("End position {} exceeds reference length {}", end, reference_length));
            }
            for pos in start..=end {
                pileup_positions.push(pos);
            }
        }
        if record.is_proper_pair() {
            properly_paired += 1;
        }
    }

    let mut coverage_counts = std::collections::HashMap::new();
    for pos in pileup_positions {
        *coverage_counts.entry(pos).or_insert(0) += 1;
    }
    total_depth = coverage_counts.values().sum::<u64>();
    bases_covered = coverage_counts.len() as u64;

    let percent_mapped = if total_reads > 0 {
        (mapped_reads as f64 / total_reads as f64) * 100.0
    } else {
        0.0
    };
    let properly_paired_percent = if total_reads > 0 {
        (properly_paired as f64 / total_reads as f64) * 100.0
    } else {
        0.0
    };
    let avg_coverage_depth = if reference_length > 0 {
        (total_depth / reference_length) as usize
    } else {
        0
    };
    let coverage_above_1x = if reference_length > 0 {
        (bases_covered as f64 / reference_length as f64) * 100.0
    } else {
        0.0
    };

    Ok(BAMMetrics {
        total_reads,
        left_reads,
        right_reads,
        percent_mapped,
        properly_paired: properly_paired_percent,
        avg_coverage_depth,
        coverage_above_1x,
        ref_mapped: percent_mapped,
        ref_properly_paired: properly_paired_percent,
        ref_avg_coverage_depth: avg_coverage_depth,
        ref_coverage_above_1x: coverage_above_1x,
    })
}



pub async fn compute_all_metrics(
    fasta_rx: Receiver<SequenceRecord>,
    ref_fasta_path: &PathBuf,
    bam_stream: Receiver<ParseOutput>,
    show_coords_stream: Receiver<ParseOutput>,
    delta_path: &Path,
    thresholds: &[usize],
) -> Result<AllMetrics> {
    // Split show-coords stream for three consumers
    let (tx1, rx1) = tokio::sync::mpsc::channel(100);
    let (tx2, rx2) = tokio::sync::mpsc::channel(100);
    let (tx3, rx3) = tokio::sync::mpsc::channel(100);
    tokio::spawn(async move {
        let mut rx = show_coords_stream;
        while let Some(output) = rx.recv().await {
            let _ = tx1.send(output.clone()).await;
            let _ = tx2.send(output.clone()).await;
            let _ = tx3.send(output).await;
        }
    });

    // Split fasta_rx to collect contig lengths and names
    let (tx_fasta, rx_fasta) = tokio::sync::mpsc::channel(100);
    let contig_length_map = Arc::new(Mutex::new(BTreeMap::new()));
    let contig_length_map_clone = Arc::clone(&contig_length_map);
    tokio::spawn(async move {
        let mut rx = fasta_rx;
        while let Some(record) = rx.recv().await {
            if let SequenceRecord::Fasta { id, seq, desc: _ } = &record {
                contig_length_map_clone.lock().await.insert(id.clone(), seq.len());
            }
            let _ = tx_fasta.send(record).await;
        }
        Result::<()>::Ok(())
    });

    // Compute reference metrics first
    let reference_metrics = compute_reference_metrics(ref_fasta_path).await?;

    // Compute assembly metrics
    let assembly_metrics = compute_assembly_metrics(rx_fasta, thresholds, rx1).await?;

    // Compute remaining metrics concurrently
    let (bam_metrics, alignment_metrics, duplication_ratio) = tokio::try_join!(
        compute_bam_metrics(bam_stream, reference_metrics.total_length as u64),
        parse_show_coords(rx2, reference_metrics.total_length as u64, assembly_metrics.total_length as u64),
        compute_duplication_ratio(rx3, assembly_metrics.total_length as u64),
    )?;

    let guard = contig_length_map.lock().await;
    let (misassemblies, misassembled_contigs, misassembled_contigs_length, local_misassemblies, scaffold_gap_ext_mis, scaffold_gap_loc_mis, unaligned_mis_contigs) =
        compute_misassemblies(delta_path, &*guard).await?;

    Ok(AllMetrics {
        assembly: AssemblyMetrics {
            num_contigs: assembly_metrics.num_contigs,
            contig_counts: assembly_metrics.contig_counts,
            total_lengths: assembly_metrics.total_lengths,
            total_length: assembly_metrics.total_length,
            largest_contig: assembly_metrics.largest_contig,
            n50: assembly_metrics.n50,
            n75: assembly_metrics.n75,
            l50: assembly_metrics.l50,
            l75: assembly_metrics.l75,
            gc_percent: assembly_metrics.gc_percent,
            unaligned_length: assembly_metrics.unaligned_length,
            unaligned_contigs: assembly_metrics.unaligned_contigs,
            n_per_100kbp: assembly_metrics.n_per_100kbp,
        },
        reference: ReferenceMetrics {
            total_length: reference_metrics.total_length,
            gc_percent: reference_metrics.gc_percent,
        },
        bam: bam_metrics,
        alignment: AlignmentMetrics {
            total_aligned_length: alignment_metrics.total_aligned_length,
            largest_alignment: alignment_metrics.largest_alignment,
            genome_fraction: alignment_metrics.genome_fraction,
            na50: alignment_metrics.na50,
            nga50: alignment_metrics.nga50,
            na75: alignment_metrics.na75,
            nga75: alignment_metrics.nga75,
            la50: alignment_metrics.la50,
            lga50: alignment_metrics.lga50,
            la75: alignment_metrics.la75,
            lga75: alignment_metrics.lga75,
            duplication_ratio,
            mismatches_per_100kbp: compute_mismatches_per_100kbp(delta_path).await?,
            indels_per_100kbp: compute_indels_per_100kbp(delta_path).await?,
            misassemblies,
            misassembled_contigs,
            misassembled_contigs_length,
            local_misassemblies,
            scaffold_gap_ext_mis,
            scaffold_gap_loc_mis,
            unaligned_mis_contigs,
        },
    })
}

pub async fn write_metrics_to_tsv(
    fasta_rx: Receiver<SequenceRecord>,
    ref_fasta_path: &PathBuf,
    bam_stream: Receiver<ParseOutput>,
    show_coords_stream: Receiver<ParseOutput>,
    delta_path: &PathBuf,
    thresholds: &[usize],
    output_path: &PathBuf,
) -> Result<()> {
    // Compute all metrics
    let metrics = compute_all_metrics(
        fasta_rx,
        ref_fasta_path,
        bam_stream,
        show_coords_stream,
        delta_path,
        thresholds,
    ).await?;

    // Build the TSV content
    let mut tsv_content = String::new();
    tsv_content.push_str("Assembly\tconsensus\n");

    // Assembly metrics with thresholds
    for &threshold in thresholds {
        tsv_content.push_str(&format!(
            "# contigs (>= {} bp)\t{}\n",
            threshold,
            metrics.assembly.contig_counts.get(&threshold).unwrap_or(&0)
        ));
    }
    for &threshold in thresholds {
        tsv_content.push_str(&format!(
            "Total length (>= {} bp)\t{}\n",
            threshold,
            metrics.assembly.total_lengths.get(&threshold).unwrap_or(&0)
        ));
    }
    tsv_content.push_str(&format!("# contigs\t{}\n", metrics.assembly.num_contigs));
    tsv_content.push_str(&format!("Largest contig\t{}\n", metrics.assembly.largest_contig));
    tsv_content.push_str(&format!("Total length\t{}\n", metrics.assembly.total_length));
    tsv_content.push_str(&format!("Reference length\t{}\n", metrics.reference.total_length));
    tsv_content.push_str(&format!("GC (%)\t{:.2}\n", metrics.assembly.gc_percent));
    tsv_content.push_str(&format!("Reference GC (%)\t{:.2}\n", metrics.reference.gc_percent));
    tsv_content.push_str(&format!("N50\t{}\n", metrics.assembly.n50));
    tsv_content.push_str(&format!("NG50\t{}\n", metrics.assembly.n50)); // Assuming NG50 same as N50
    tsv_content.push_str(&format!("N75\t{}\n", metrics.assembly.n75));
    tsv_content.push_str(&format!("NG75\t{}\n", metrics.assembly.n75)); // Assuming NG75 same as N75
    tsv_content.push_str(&format!("L50\t{}\n", metrics.assembly.l50));
    tsv_content.push_str(&format!("LG50\t{}\n", metrics.assembly.l50)); // Assuming LG50 same as L50
    tsv_content.push_str(&format!("L75\t{}\n", metrics.assembly.l75));
    tsv_content.push_str(&format!("LG75\t{}\n", metrics.assembly.l75)); // Assuming LG75 same as L75
    tsv_content.push_str(&format!("# total reads\t{}\n", metrics.bam.total_reads));
    tsv_content.push_str(&format!("# left\t{}\n", metrics.bam.left_reads));
    tsv_content.push_str(&format!("# right\t{}\n", metrics.bam.right_reads));
    tsv_content.push_str(&format!("Mapped (%)\t{:.1}\n", metrics.bam.percent_mapped));
    tsv_content.push_str(&format!("Reference mapped (%)\t{:.1}\n", metrics.bam.ref_mapped));
    tsv_content.push_str(&format!("Properly paired (%)\t{:.2}\n", metrics.bam.properly_paired));
    tsv_content.push_str(&format!("Reference properly paired (%)\t{:.2}\n", metrics.bam.ref_properly_paired));
    tsv_content.push_str(&format!("Avg. coverage depth\t{}\n", metrics.bam.avg_coverage_depth));
    tsv_content.push_str(&format!("Reference avg. coverage depth\t{}\n", metrics.bam.ref_avg_coverage_depth));
    tsv_content.push_str(&format!("Coverage >= 1x (%)\t{:.2}\n", metrics.bam.coverage_above_1x));
    tsv_content.push_str(&format!("Reference coverage >= 1x (%)\t{:.2}\n", metrics.bam.ref_coverage_above_1x));
    tsv_content.push_str(&format!("# misassemblies\t{}\n", metrics.alignment.misassemblies));
    tsv_content.push_str(&format!("# misassembled contigs\t{}\n", metrics.alignment.misassembled_contigs));
    tsv_content.push_str(&format!("Misassembled contigs length\t{}\n", metrics.alignment.misassembled_contigs_length));
    tsv_content.push_str(&format!("# local misassemblies\t{}\n", metrics.alignment.local_misassemblies));
    tsv_content.push_str(&format!("# scaffold gap ext. mis.\t{}\n", metrics.alignment.scaffold_gap_ext_mis));
    tsv_content.push_str(&format!("# scaffold gap loc. mis.\t{}\n", metrics.alignment.scaffold_gap_loc_mis));
    tsv_content.push_str(&format!("# unaligned mis. contigs\t{}\n", metrics.alignment.unaligned_mis_contigs));
    tsv_content.push_str(&format!("# unaligned contigs\t{} + {} part\n", metrics.assembly.unaligned_contigs, 0)); // Assuming no partial unaligned
    tsv_content.push_str(&format!("Unaligned length\t{}\n", metrics.assembly.unaligned_length));
    tsv_content.push_str(&format!("Genome fraction (%)\t{:.3}\n", metrics.alignment.genome_fraction));
    tsv_content.push_str(&format!("Duplication ratio\t{:.3}\n", metrics.alignment.duplication_ratio));
    tsv_content.push_str(&format!("# N's per 100 kbp\t{:.2}\n", metrics.assembly.n_per_100kbp));
    tsv_content.push_str(&format!("# mismatches per 100 kbp\t{:.2}\n", metrics.alignment.mismatches_per_100kbp));
    tsv_content.push_str(&format!("# indels per 100 kbp\t{:.2}\n", metrics.alignment.indels_per_100kbp));
    tsv_content.push_str(&format!("Largest alignment\t{}\n", metrics.alignment.largest_alignment));
    tsv_content.push_str(&format!("Total aligned length\t{}\n", metrics.alignment.total_aligned_length));
    tsv_content.push_str(&format!("NA50\t{}\n", metrics.alignment.na50));
    tsv_content.push_str(&format!("NGA50\t{}\n", metrics.alignment.nga50));
    tsv_content.push_str(&format!("NA75\t{}\n", metrics.alignment.na75));
    tsv_content.push_str(&format!("NGA75\t{}\n", metrics.alignment.nga75));
    tsv_content.push_str(&format!("LA50\t{}\n", metrics.alignment.la50));
    tsv_content.push_str(&format!("LGA50\t{}\n", metrics.alignment.lga50));
    tsv_content.push_str(&format!("LA75\t{}\n", metrics.alignment.la75));
    tsv_content.push_str(&format!("LGA75\t{}\n", metrics.alignment.lga75));

    // Write to file
    let mut file = File::create(output_path).await
        .map_err(|e| anyhow!("Failed to create output file {}: {}", output_path.display(), e))?;
    file.write_all(tsv_content.as_bytes()).await
        .map_err(|e| anyhow!("Failed to write to output file {}: {}", output_path.display(), e))?;
    file.flush().await
        .map_err(|e| anyhow!("Failed to flush output file {}: {}", output_path.display(), e))?;

    Ok(())
}