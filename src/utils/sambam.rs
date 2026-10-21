//! SAM/BAM alignment counter and insert size statistics.

use crate::config::defs::PipelineError;
use tokio::task::JoinHandle;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::sync::mpsc::Receiver as StdReceiver;

use anyhow::{anyhow, Context, Result};
use log::{self, debug, info, warn};
use dashmap::DashMap;
use futures::StreamExt;
use noodles::bam;
use noodles::bam::r#async::io::Reader as BamAsyncReader;
use noodles::bam::record::Record;
use noodles::bam::record::Record as BamRecord;
use tokio::fs::File;
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tokio_util::io::StreamReader;
use streampile::StreamingPileupBuilder;
use bytes::Bytes;
use serde_json::{json, Value};
use tokio::io::AsyncWriteExt;
use tokio::io::BufWriter;
use needletail::parse_fastx_file;
use noodles::sam::alignment::record::Flags;
use rayon::prelude::*;

use crate::config::defs::{ClusterInfo, RunConfig};
use crate::utils::streams::ParseOutput;

#[derive(Debug, Clone)]
pub struct InsertSizeStats {
    pub insert_sizes: Vec<(u32, u64)>, // (size, count) sorted by size
    pub total_proper_pairs: u64,
    pub mean: f64,
    pub median: f64,
    pub stddev: f64,
    pub mapped_proper_pairs: u64, // should == total_proper_pairs here
}

pub async fn stream_sam_alignment_counter(
    rx: mpsc::Receiver<ParseOutput>,
    early_exit: bool,
) -> Result<u64, anyhow::Error> {
    let mut stream = ReceiverStream::new(rx);
    let mut counter: u64 = 0;

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(chunk) => {
                let chunk_str = String::from_utf8_lossy(&chunk);
                for line in chunk_str.lines() {
                    if !line.starts_with('@') {
                        counter += 1;
                        if early_exit {
                            return Ok(1); // Early exit: Not empty
                        }
                    }
                }
            }
            _ => {
                return Err(anyhow!(
                    "Unexpected ParseOutput in SAM stream; expected Bytes"
                ))
            }
        }
    }

    Ok(counter)
}

/// Bridges the pipeline's asynchronous BAM byte stream into a synchronous
/// `std::io::Read` implementation for noodles + streampile.
///
/// This is used only inside `spawn_blocking`.
struct BlockingChannelReader {
    rx: mpsc::Receiver<ParseOutput>,
    runtime: tokio::runtime::Handle,
    buffer: Bytes,
    offset: usize,
}

impl BlockingChannelReader {
    fn new(
        rx: mpsc::Receiver<ParseOutput>,
        runtime: tokio::runtime::Handle,
    ) -> Self {
        Self {
            rx,
            runtime,
            buffer: Bytes::new(),
            offset: 0,
        }
    }
}

impl io::Read for BlockingChannelReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }

        loop {
            if self.offset < self.buffer.len() {
                let remaining = &self.buffer[self.offset..];
                let n = remaining.len().min(buf.len());

                buf[..n].copy_from_slice(&remaining[..n]);
                self.offset += n;

                if self.offset == self.buffer.len() {
                    self.buffer = Bytes::new();
                    self.offset = 0;
                }

                return Ok(n);
            }

            match self.runtime.block_on(self.rx.recv()) {
                Some(ParseOutput::Bytes(bytes)) => {
                    if bytes.is_empty() {
                        continue;
                    }

                    self.buffer = bytes;
                    self.offset = 0;
                }

                Some(_) => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "generate_assembly_coverage: expected \
                         ParseOutput::Bytes in BAM stream",
                    ));
                }

                None => {
                    return Ok(0);
                }
            }
        }
    }
}

fn anyhow_to_io(e: anyhow::Error) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e)
}

/// Generate read-to-contig mappings and unweighted contig read counts from
/// a name-sorted BAM stream.
///
/// This is the generic equivalent of the legacy
/// `generate_info_from_sam()` helper.
///
/// Legacy behavior:
///   * process alignments in input order
///   * skip SAM headers (handled by the BAM reader)
///   * keep only the first alignment encountered for each read ID
///   * count one unique read for its reported contig
///   * populate read2contig for mapped reads
///   * use "*" for reads with no reference sequence
///   * filter contigs below `min_contig_size`
///
/// Short-read duplicate-cluster weighting is intentionally NOT done here.
/// The caller that has duplicate-cluster information can transform the
/// unweighted contig statistics afterward.
pub async fn generate_info_from_bam_stream(
    rx: mpsc::Receiver<ParseOutput>,
    min_contig_size: u64,
) -> Result<(
    HashMap<String, String>, // read2contig
    HashMap<String, u64>,    // unweighted contig_stats
    HashMap<String, usize>,  // contig_unique_counts
)> {
    let byte_stream = ReceiverStream::new(rx).map(|item| match item {
        ParseOutput::Bytes(bytes) => Ok(bytes),
        _ => Err(anyhow_to_io(anyhow!(
            "BAM stream received non-Bytes variant"
        ))),
    });

    let stream_reader = StreamReader::new(byte_stream);
    let mut bam_reader = BamAsyncReader::new(stream_reader);

    let header = bam_reader
        .read_header()
        .await
        .map_err(|e| anyhow!("BAM header error: {}", e))?;

    let mut read2contig: HashMap<String, String> = HashMap::new();
    let mut contig_stats: HashMap<String, u64> = HashMap::new();
    let mut contig_unique_counts: HashMap<String, usize> = HashMap::new();
    let mut seen_reads: HashSet<String> = HashSet::new();

    let mut record = BamRecord::default();

    while bam_reader.read_record(&mut record).await? > 0 {
        // The legacy Python helper does not inspect SAM flags here.
        //
        // In particular, do NOT filter unmapped/secondary/supplementary
        // records at this generic layer. The Python code simply takes the
        // first alignment encountered for each read ID.
        let read_name = record
            .name()
            .and_then(|name| std::str::from_utf8(name.as_ref()).ok())
            .unwrap_or("*")
            .to_string();

        // Global, input-order deduplication.
        //
        // This must be one set for the entire BAM stream, not one per batch.
        if !seen_reads.insert(read_name.clone()) {
            continue;
        }

        let contig_name = match record.reference_sequence_id() {
            Some(Ok(rid)) => match header.reference_sequences().get_index(rid) {
                Some((name, _)) => name.to_string(),
                None => "*".to_string(),
            },
            _ => "*".to_string(),
        };

        // Generic unweighted statistics: one unique read contributes one.
        *contig_stats
            .entry(contig_name.clone())
            .or_insert(0) += 1;

        *contig_unique_counts
            .entry(contig_name.clone())
            .or_insert(0) += 1;

        // Match legacy generate_info_from_sam():
        // only mapped reads are entered into read2contig.
        if contig_name != "*" {
            read2contig.insert(read_name, contig_name);
        }
    }

    // Match the legacy post-pass filtering:
    // remove contigs whose unique-read count is below the threshold.
    contig_stats.retain(|contig, count| {
        *count >= min_contig_size
            && contig_unique_counts
            .get(contig)
            .copied()
            .unwrap_or(0) as u64
            >= min_contig_size
    });

    contig_unique_counts.retain(|contig, _| contig_stats.contains_key(contig));

    info!(
        "generate_info_from_bam_stream: {} unique reads, {} mapped read-to-contig entries, {} contigs retained",
        seen_reads.len(),
        read2contig.len(),
        contig_stats.len()
    );

    Ok((
        read2contig,
        contig_stats,
        contig_unique_counts,
    ))
}


/// Insert stats from an exisiting bam file
///
///
/// # Arguments
///
///  * `bam_path`
/// * `max_records_to_check` - limits records
/// * `_thread_pool - left in for possible future rayon based
///
/// # Returns
///
/// Rreult <InsertSizeStats>
pub async fn compute_insert_size_stats_from_bam(
    bam_path: PathBuf,
    max_records_to_check: Option<usize>,
) -> Result<InsertSizeStats> {
    let file = File::open(&bam_path)
        .await
        .map_err(|e| anyhow!("Failed to open BAM file {}: {}", bam_path.display(), e))?;

    let mut bam_reader = bam::r#async::io::Reader::new(file);

    // noodles-bam requires the header to be read first
    let _header = bam_reader
        .read_header()
        .await
        .map_err(|e| anyhow!("BAM header error: {}", e))?;

    let histogram = Arc::new(Mutex::new(HashMap::<u32, u64>::new()));
    let mut total_proper_pairs = 0u64;
    let mut processed = 0usize;

    const MAX_REASONABLE_INSERT: u32 = 10_000;
    const MIN_REASONABLE_INSERT: u32 = 10;

    let mut record = Record::default();

    loop {
        let bytes_read = match bam_reader.read_record(&mut record).await {
            Ok(n) => n,
            Err(e) => return Err(anyhow!("BAM record read error: {}", e)),
        };

        if bytes_read == 0 {
            break;
        }

        processed += 1;
        if let Some(max) = max_records_to_check {
            if processed > max {
                break;
            }
        }

        let flags = record.flags();

        // Standard proper-pair filter (matches Picard logic)
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || !flags.is_properly_segmented()
            || record.mate_reference_sequence_id().is_none()
            || record.template_length() == 0
        {
            continue;
        }

        // Count each proper pair exactly once (use first segment only)
        if !flags.is_first_segment() {
            continue;
        }

        let insert_size = record.template_length().abs() as u32;

        if insert_size < MIN_REASONABLE_INSERT || insert_size > MAX_REASONABLE_INSERT {
            continue;
        }

        {
            let mut hist = histogram.lock().unwrap();
            *hist.entry(insert_size).or_insert(0) += 1;
        }

        total_proper_pairs += 1;
    }

    // ─────────────────────────────────────────────────────────────
    // Build sorted histogram + stats (unchanged from your original)
    // ─────────────────────────────────────────────────────────────
    let mut insert_sizes: Vec<(u32, u64)> = {
        let hist = histogram.lock().unwrap();
        hist.iter().map(|(&size, &count)| (size, count)).collect()
    };

    insert_sizes.sort_by_key(|&(size, _)| size);

    let total_weight: u64 = insert_sizes.iter().map(|&(_, c)| c).sum();

    let mean = if total_weight > 0 {
        insert_sizes
            .iter()
            .map(|&(s, c)| s as f64 * c as f64)
            .sum::<f64>()
            / total_weight as f64
    } else {
        0.0
    };

    // Median
    let mut cumulative = 0u64;
    let mut median = 0.0;
    for &(size, count) in &insert_sizes {
        cumulative += count;
        if cumulative >= (total_weight + 1) / 2 {
            median = size as f64;
            break;
        }
    }

    let variance = if total_weight > 0 {
        insert_sizes
            .iter()
            .map(|&(s, c)| {
                let dev = s as f64 - mean;
                dev * dev * c as f64
            })
            .sum::<f64>()
            / total_weight as f64
    } else {
        0.0
    };

    let stddev = variance.sqrt();

    Ok(InsertSizeStats {
        insert_sizes,
        total_proper_pairs,
        mean,
        median,
        stddev,
        mapped_proper_pairs: total_proper_pairs,
    })
}

/// Generate assembly contig coverage statistics from a streamed BAM.
///
/// The BAM is consumed once through `streampile` using the existing streaming
/// pipeline. No intermediate BAM file is written.
///
/// Coverage is represented over the actual FASTA contig length, so positions
/// with no aligned reads remain explicit zeroes.
///
/// Pileup depth deliberately uses `unfiltered_depth()`, matching the legacy
/// Python:
///
///     pileup_column.get_num_aligned()
///
/// The alignment filter is explicitly configured to match the legacy pysam
/// pileup flag filter:
///
///     UNMAPPED | SECONDARY | QC_FAIL | DUPLICATE
///
/// Supplementary alignments are therefore retained.
///
/// The final per-contig statistics remain parallelized using the existing
/// Rayon thread pool.
///
/// JSON uses the legacy Python field names:
///
///     coverage, avg, p0, p100, p25, p50, p75,
///     avg2xcnt, cnt0, cnt1, cnt2, contig_len
///
/// CSV uses the legacy Python column names:
///
///     contig_name,avg,min,max,p25,p50,p75,avg2xcnt,cnt0,cnt1,cnt2
pub async fn generate_assembly_coverage(
    config: Arc<RunConfig>,
    assembly_bam_rx: ReceiverStream<ParseOutput>,
    contigs_fasta_path: &PathBuf,
    coverage_json_path: &PathBuf,
    coverage_csv_path: &PathBuf,
) -> Result<()> {
    const MIN_CONTIG_FILE_SIZE: u64 = 50;

    // ---------------------------------------------------------------------
    // Early exit for failed/empty assemblies.
    // ---------------------------------------------------------------------
    let meta = tokio::fs::metadata(contigs_fasta_path)
        .await
        .with_context(|| {
            format!(
                "Failed to stat contigs FASTA {}",
                contigs_fasta_path.display()
            )
        })?;

    if meta.len() < MIN_CONTIG_FILE_SIZE {
        tokio::fs::write(coverage_json_path, b"{}").await?;
        tokio::fs::write(coverage_csv_path, b"No Contigs\n").await?;
        return Ok(());
    }

    info!(
        "Streaming assembly coverage: {} -> {}",
        contigs_fasta_path.display(),
        coverage_json_path.display()
    );

    // ---------------------------------------------------------------------
    // Load contig order and actual FASTA lengths.
    //
    // These lengths are authoritative for the coverage vectors. A position
    // with no pileup remains zero.
    // ---------------------------------------------------------------------
    let (contig_order, contig_lengths): (Vec<String>, HashMap<String, usize>) = {
        let mut reader = parse_fastx_file(contigs_fasta_path)
            .map_err(|e| {
                anyhow!(
                    "Failed to open contigs FASTA {}: {}",
                    contigs_fasta_path.display(),
                    e
                )
            })?;

        let mut order = Vec::with_capacity(100_000);
        let mut lengths = HashMap::with_capacity(100_000);

        while let Some(result) = reader.next() {
            let record = result
                .map_err(|e| anyhow!("Failed to parse contigs FASTA: {}", e))?;

            let id = String::from_utf8_lossy(record.id()).into_owned();
            let len = record.seq().len();

            if len > 0 {
                order.push(id.clone());
                lengths.insert(id, len);
            }
        }

        if order.is_empty() {
            tokio::fs::write(coverage_json_path, b"{}").await?;
            tokio::fs::write(coverage_csv_path, b"No Contigs\n").await?;
            return Ok(());
        }

        (order, lengths)
    };

    // ---------------------------------------------------------------------
    // The BAM stream is consumed by a synchronous noodles/streampile worker.
    //
    // Keep copies of the FASTA metadata for the statistics phase after the
    // blocking pileup pass completes.
    // ---------------------------------------------------------------------
    let pileup_contig_lengths = contig_lengths.clone();

    let rx = assembly_bam_rx.into_inner();
    let runtime = tokio::runtime::Handle::current();

    let coverage_by_contig = tokio::task::spawn_blocking(move || -> Result<
        HashMap<String, Vec<u32>>
    > {
        let reader = BlockingChannelReader::new(rx, runtime);

        let mut bam_reader = bam::io::Reader::new(reader);

        let header = bam_reader
            .read_header()
            .map_err(|e| {
                anyhow!(
                    "Failed to read BAM header for assembly coverage: {}",
                    e
                )
            })?;

        let references = header.reference_sequences();

        info!(
            "Assembly coverage pileup: BAM header contains {} references",
            references.len()
        );

        // Match pysam pileup's legacy default flag filtering, while retaining
        // supplementary alignments.
        let exclude_flags =
            Flags::UNMAPPED
                | Flags::SECONDARY
                | Flags::QC_FAIL
                | Flags::DUPLICATE;

        let mut builder =
            StreamingPileupBuilder::new(bam_reader, &header)
                .map_err(|e| {
                    anyhow!(
                        "Failed to create StreamingPileupBuilder: {}",
                        e
                    )
                })?
                .exclude_flags(exclude_flags);

        let mut coverage_by_contig =
            HashMap::with_capacity(pileup_contig_lengths.len());

        let mut total_columns: u64 = 0;
        let mut total_depth: u64 = 0;

        // `StreamingPileupBuilder` is forward-only. Iterate references in BAM
        // header order rather than arbitrary HashMap/FASTA order.
        for (reference_index, (reference_name, reference_sequence)) in
            references.iter().enumerate()
        {
            let contig_name = reference_name.to_string();

            let Some(&contig_len) =
                pileup_contig_lengths.get(&contig_name)
            else {
                warn!(
                    "Assembly coverage: BAM reference '{}' is not present \
                     in contigs FASTA; skipping",
                    contig_name
                );
                continue;
            };

            let bam_len =
                reference_sequence.length().get() as usize;

            if bam_len != contig_len {
                warn!(
                    "Assembly coverage: length mismatch for {}: \
                     BAM={} FASTA={}",
                    contig_name,
                    bam_len,
                    contig_len
                );
            }

            info!(
                "Assembly coverage: processing reference {}/{}: {} ({} bp)",
                reference_index + 1,
                references.len(),
                contig_name,
                contig_len
            );

            // One depth value for every actual FASTA reference position.
            let mut coverage = vec![0u32; contig_len];

            let mut columns = builder
                .columns(&contig_name, 0, contig_len)
                .map_err(|e| {
                    anyhow!(
                        "Failed to create pileup columns for {}: {}",
                        contig_name,
                        e
                    )
                })?;

            let mut contig_columns: u64 = 0;
            let mut contig_depth: u64 = 0;

            // columns() sweeps every requested reference position, including
            // positions with zero coverage.
            while let Some(pileup) = columns.next_pileup() {
                let pileup = pileup.map_err(|e| {
                    anyhow!(
                        "Pileup failed for {} after {} columns: {}",
                        contig_name,
                        contig_columns,
                        e
                    )
                })?;

                let position = pileup.position();

                if position >= coverage.len() {
                    return Err(anyhow!(
                        "Pileup position {} is outside FASTA length {} \
                         for contig '{}'",
                        position,
                        coverage.len(),
                        contig_name
                    ));
                }

                // Direct equivalent of the legacy:
                //
                //     pileup_column.get_num_aligned()
                //
                // Do NOT use filtered_depth() here.
                let depth = pileup.unfiltered_depth() as u32;

                coverage[position] = depth;

                contig_columns += 1;
                contig_depth += depth as u64;
                total_columns += 1;
                total_depth += depth as u64;
            }

            info!(
                "Assembly coverage: completed {} — {} columns, \
                 cumulative unfiltered depth {}",
                contig_name,
                contig_columns,
                contig_depth
            );

            coverage_by_contig.insert(contig_name, coverage);
        }

        // Explicitly close so any remaining source-side error is surfaced.
        builder
            .close()
            .map_err(|e| {
                anyhow!(
                    "Failed to close assembly coverage pileup builder: {}",
                    e
                )
            })?;

        info!(
            "Assembly coverage pileup complete — {} columns, \
             cumulative unfiltered depth {}",
            total_columns,
            total_depth
        );

        Ok(coverage_by_contig)
    })
        .await
        .map_err(|e| {
            anyhow!(
            "Assembly coverage pileup blocking task failed: {}",
            e
        )
        })??;

    // ---------------------------------------------------------------------
    // Parallel per-contig statistics.
    // ---------------------------------------------------------------------
    let coverage_by_contig = Arc::new(coverage_by_contig);

    let stats: Vec<(String, Value)> = config.thread_pool.install(|| {
        contig_order
            .par_iter()
            .filter_map(|name| {
                let coverage = coverage_by_contig.get(name)?;
                let len = *contig_lengths.get(name)?;

                if len == 0 {
                    return None;
                }

                let sum: u64 =
                    coverage.iter().map(|&d| d as u64).sum();

                let avg = sum as f64 / len as f64;

                let mut sorted = coverage.clone();
                sorted.sort_unstable();

                let p25 =
                    sorted[((len as f64 * 0.25) as usize).min(len - 1)];

                let p50 =
                    sorted[((len as f64 * 0.50) as usize).min(len - 1)];

                let p75 =
                    sorted[((len as f64 * 0.75) as usize).min(len - 1)];

                let min = *sorted.first()?;
                let max = *sorted.last()?;

                let avg2xcnt = coverage
                    .iter()
                    .filter(|&&depth| depth as f64 > 2.0 * avg)
                    .count() as f64
                    / len as f64;

                let cnt0 = coverage
                    .iter()
                    .filter(|&&depth| depth == 0)
                    .count() as f64
                    / len as f64;

                let cnt1 = coverage
                    .iter()
                    .filter(|&&depth| depth == 1)
                    .count() as f64
                    / len as f64;

                let cnt2 = coverage
                    .iter()
                    .filter(|&&depth| depth == 2)
                    .count() as f64
                    / len as f64;

                let coverage_array: Vec<u64> =
                    coverage.iter().map(|&d| d as u64).collect();

                Some((
                    name.clone(),
                    json!({
                        "coverage": coverage_array,
                        "avg": avg,
                        "p0": min as u64,
                        "p100": max as u64,
                        "p25": p25 as u64,
                        "p50": p50 as u64,
                        "p75": p75 as u64,
                        "avg2xcnt": avg2xcnt,
                        "cnt0": cnt0,
                        "cnt1": cnt1,
                        "cnt2": cnt2,
                        "contig_len": len as u64,
                    }),
                ))
            })
            .collect()
    });

    // ---------------------------------------------------------------------
    // Write JSON and CSV asynchronously.
    // ---------------------------------------------------------------------
    let mut json_file =
        BufWriter::new(File::create(coverage_json_path).await?);

    let mut csv_file =
        BufWriter::new(File::create(coverage_csv_path).await?);

    csv_file
        .write_all(
            b"contig_name,avg,min,max,p25,p50,p75,avg2xcnt,cnt0,cnt1,cnt2\n",
        )
        .await?;

    json_file.write_all(b"{").await?;

    for (i, (name, stat)) in stats.iter().enumerate() {
        if i > 0 {
            json_file.write_all(b",").await?;
        }

        let key = serde_json::to_string(name)?;
        let value = serde_json::to_string(stat)?;

        json_file.write_all(key.as_bytes()).await?;
        json_file.write_all(b":").await?;
        json_file.write_all(value.as_bytes()).await?;
    }

    json_file.write_all(b"}").await?;
    json_file.flush().await?;

    for (name, stat) in &stats {
        let object = stat
            .as_object()
            .ok_or_else(|| {
                anyhow!(
                    "Coverage statistics for '{}' are not an object",
                    name
                )
            })?;

        let line = format!(
            "{},{:.6},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6}\n",
            name,
            object["avg"].as_f64().unwrap_or(0.0),
            object["p0"].as_u64().unwrap_or(0),
            object["p100"].as_u64().unwrap_or(0),
            object["p25"].as_u64().unwrap_or(0),
            object["p50"].as_u64().unwrap_or(0),
            object["p75"].as_u64().unwrap_or(0),
            object["avg2xcnt"].as_f64().unwrap_or(0.0),
            object["cnt0"].as_f64().unwrap_or(0.0),
            object["cnt1"].as_f64().unwrap_or(0.0),
            object["cnt2"].as_f64().unwrap_or(0.0),
        );

        csv_file.write_all(line.as_bytes()).await?;
    }

    csv_file.flush().await?;

    info!(
        "Assembly coverage complete — {} contigs written",
        stats.len()
    );

    Ok(())
}

/// Filter minimap2 SAM alignments with excessive terminal clipping.
///
/// This is the streaming Rust equivalent of filter_clipped_alignments.py.
///
/// Behavior intentionally matches the legacy implementation:
///
/// - SAM headers are passed through unchanged.
/// - Already-unmapped records (CIGAR == "*") are passed through unchanged.
/// - Only clipping at the beginning and end of the CIGAR is considered.
/// - Both soft clipping (`S`) and hard clipping (`H`) count toward clipping.
/// - Hard-clipped bases are added to SEQ length when calculating the
///   effective read length because they are absent from SEQ.
/// - Alignments with percent_clipped < max_percent are retained.
/// - Alignments with percent_clipped >= max_percent are rewritten as
///   unmapped SAM records while preserving QNAME, SEQ, and QUAL.
///
/// The output remains a SAM byte stream.
pub async fn filter_clipped_alignments(
    mut input_stream: ReceiverStream<ParseOutput>,
    max_percent: f64,
    buffer_size: usize,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        JoinHandle<anyhow::Result<(), anyhow::Error>>,
    ),
    PipelineError,
> {
    let (tx, rx) = tokio::sync::mpsc::channel(buffer_size);

    let task = tokio::spawn(async move {
        while let Some(item) = input_stream.next().await {
            let bytes = match item {
                ParseOutput::Bytes(bytes) => bytes,
                other => {
                    return Err(anyhow!(
                        "filter_clipped_alignments expected SAM bytes, got {:?}",
                        other
                    ));
                }
            };

            // ----------------------------------------------------------
            // SAM headers are preserved exactly.
            // ----------------------------------------------------------
            if bytes.first() == Some(&b'@') {
                tx.send(ParseOutput::Bytes(bytes))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "filter_clipped_alignments output receiver dropped"
                        )
                    })?;

                continue;
            }

            // ----------------------------------------------------------
            // Parse the mandatory SAM fields.
            //
            // We only need:
            //   0  QNAME
            //   5  CIGAR
            //   9  SEQ
            //   10 QUAL
            // ----------------------------------------------------------
            let line = std::str::from_utf8(&bytes)
                .map_err(|e| {
                    anyhow!(
                        "filter_clipped_alignments received invalid UTF-8 SAM record: {}",
                        e
                    )
                })?;

            let line = line.trim_end_matches(&['\r', '\n'][..]);

            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() < 11 {
                return Err(anyhow!(
                    "Malformed SAM record in filter_clipped_alignments: \
                     expected at least 11 fields, found {}: {}",
                    fields.len(),
                    line
                ));
            }

            let qname = fields[0];
            let cigar = fields[5];
            let seq = fields[9];
            let qual = fields[10];

            // ----------------------------------------------------------
            // Legacy Python:
            //
            //     elif row[5] != "*":
            //         ...
            //     else:
            //         output.writerow(row)
            //
            // So already-unmapped records pass through unchanged.
            // ----------------------------------------------------------
            if cigar == "*" {
                tx.send(ParseOutput::Bytes(bytes))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "filter_clipped_alignments output receiver dropped"
                        )
                    })?;

                continue;
            }

            // ----------------------------------------------------------
            // Calculate terminal clipping.
            //
            // Python patterns:
            //
            //     ^\d+[SH]
            //     \d+[SH]$
            //
            // Therefore ONLY the first and last CIGAR operations count.
            // Internal S/H operations, if present, do not.
            // ----------------------------------------------------------
            let (front_clipping, front_hard_clipping) =
                terminal_cigar_clipping(cigar, true)?;

            let (end_clipping, end_hard_clipping) =
                terminal_cigar_clipping(cigar, false)?;

            // Python starts with:
            //
            //     alength = len(row[9])
            //
            // and adds hard-clipped bases because they are not represented
            // in SEQ.
            let alignment_length = seq.len()
                + front_hard_clipping
                + end_hard_clipping;

            if alignment_length == 0 {
                return Err(anyhow!(
                    "SAM record {} has zero effective alignment length",
                    qname
                ));
            }

            let total_clipping =
                front_clipping + end_clipping;

            let percent_clipped =
                (total_clipping as f64 / alignment_length as f64) * 100.0;

            // ----------------------------------------------------------
            // Python uses:
            //
            //     if percent_clipped < max_percent:
            //
            // Therefore equality is filtered.
            // ----------------------------------------------------------
            if percent_clipped < max_percent {
                tx.send(ParseOutput::Bytes(bytes))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "filter_clipped_alignments output receiver dropped"
                        )
                    })?;

                continue;
            }

            // ----------------------------------------------------------
            // Excessively clipped alignment:
            //
            // Match the legacy UNMAPPED_READ template exactly:
            //
            // QNAME  FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
            // id     4    *     0   0    *     *     0     0    seq qual
            //
            // The original QNAME, SEQ, and QUAL are retained.
            // Everything else is replaced.
            // ----------------------------------------------------------
            let unmapped = format!(
                "{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\n",
                qname,
                seq,
                qual,
            );

            tx.send(ParseOutput::Bytes(bytes::Bytes::from(unmapped)))
                .await
                .map_err(|_| {
                    anyhow!(
                        "filter_clipped_alignments output receiver dropped"
                    )
                })?;
        }

        Ok(())
    });

    Ok((
        ReceiverStream::new(rx),
        task,
    ))
}


/// Return terminal clipping information for one end of a CIGAR string.
///
/// Returns:
///
///     (clipping_length, hard_clipping_length)
///
/// Examples:
///
///     100M          -> (0, 0)
///     20S100M       -> front: (20, 0)
///     20H100M       -> front: (20, 20)
///     100M30S       -> end:   (30, 0)
///     100M30H       -> end:   (30, 30)
///
/// This deliberately matches the legacy regular expressions:
///
///     ^\d+[SH]
///     \d+[SH]$
///
/// It does NOT combine adjacent clipping operations. For example,
/// `10H20S100M` matches only `10H` at the front, just as the Python
/// regular expression does.
fn terminal_cigar_clipping(
    cigar: &str,
    front: bool,
) -> anyhow::Result<(usize, usize)> {
    let cigar_bytes = cigar.as_bytes();

    if cigar_bytes.is_empty() {
        return Ok((0, 0));
    }

    if front {
        // --------------------------------------------------------------
        // Equivalent to:
        //
        //     re.compile(r"^\d+[SH]")
        // --------------------------------------------------------------
        let mut end = 0;

        while end < cigar_bytes.len()
            && cigar_bytes[end].is_ascii_digit()
        {
            end += 1;
        }

        // No leading number or no operation following it.
        if end == 0 || end >= cigar_bytes.len() {
            return Ok((0, 0));
        }

        let op = cigar_bytes[end];

        if op != b'S' && op != b'H' {
            return Ok((0, 0));
        }

        let clipping_length: usize = cigar[..end]
            .parse()
            .map_err(|e| {
                anyhow!(
                    "Invalid leading clipping length in CIGAR '{}': {}",
                    cigar,
                    e
                )
            })?;

        let hard_clipping_length = if op == b'H' {
            clipping_length
        } else {
            0
        };

        Ok((
            clipping_length,
            hard_clipping_length,
        ))
    } else {
        // --------------------------------------------------------------
        // Equivalent to:
        //
        //     re.compile(r"\d+[SH]$")
        // --------------------------------------------------------------
        let op = *cigar_bytes
            .last()
            .ok_or_else(|| anyhow!("Empty CIGAR"))?;

        if op != b'S' && op != b'H' {
            return Ok((0, 0));
        }

        // Start immediately before the final S/H and walk backwards
        // through its numeric length.
        let mut start = cigar_bytes.len() - 1;

        while start > 0
            && cigar_bytes[start - 1].is_ascii_digit()
        {
            start -= 1;
        }

        // A terminal S/H with no numeric length isn't a valid clipping
        // operation for purposes of the legacy regex.
        if start == cigar_bytes.len() - 1 {
            return Ok((0, 0));
        }

        let clipping_length: usize =
            cigar[start..cigar_bytes.len() - 1]
                .parse()
                .map_err(|e| {
                    anyhow!(
                        "Invalid trailing clipping length in CIGAR '{}': {}",
                        cigar,
                        e
                    )
                })?;

        let hard_clipping_length = if op == b'H' {
            clipping_length
        } else {
            0
        };

        Ok((
            clipping_length,
            hard_clipping_length,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tokio::sync::mpsc;

    #[tokio::test]
    async fn test_stream_sam_alignment_counter() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let sam_data = vec![
            ParseOutput::Bytes(b"@HD\tVN:1.6\n@SQ\tSN:ERCC-00002\tLN:1061\nread1\t77\t*\t0\t0\t*\t*\t0\t0\tATCG\tIIII\n".to_vec().into()),
        ];

        tokio::spawn(async move {
            for item in sam_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_sam_alignment_counter(rx, true).await?;
        assert_eq!(count, 1, "Should count 1 SAM alignment with early_exit");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_sam_alignment_counter_empty() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let sam_data = vec![ParseOutput::Bytes(
            b"@HD\tVN:1.6\n@SQ\tSN:ERCC-00002\tLN:1061\n"
                .to_vec()
                .into(),
        )];

        tokio::spawn(async move {
            for item in sam_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_sam_alignment_counter(rx, false).await?;
        assert_eq!(count, 0, "Should count 0 for header-only SAM");
        Ok(())
    }
}
