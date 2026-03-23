// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};
use std::env::args;
use std::sync::Arc;
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use log::{self, LevelFilter, debug, info, error, warn};
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;
use fst::Map;
use ahash::AHashMap;
use serde::{Deserialize, Serialize};
use tokio_stream::wrappers::ReceiverStream;
use tokio::fs::File as TokioFile;
use tokio::sync::mpsc;
use tokio::sync::mpsc::{channel, Sender};
use tokio::io::{AsyncBufReadExt, BufReader, BufWriter, AsyncWriteExt};
use tokio_stream::wrappers::BroadcastStream;
use tokio::sync::Semaphore;
use tokio::task::JoinHandle;
use dashmap::DashMap;
use tokio::time::{sleep, Duration, Instant};


use crate::config::defs::{Taxid, Lineage, NT_TAG, NR_TAG, RunConfig, ClusterInfo};
use crate::utils::streams::ParseOutput;
use crate::utils::taxonomy::validate_taxid_lineage;
use crate::utils::streams::ToBytes;
use crate::utils::file::write_byte_stream_to_file;
use crate::utils::system::{compute_batch_size, compute_phase_concurrency};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigSummaryEntry {
    contig_name: String,
    common_name: Option<String>,
    category_name: Option<String>,
    score: Option<f64>,
    db_type: String,
    reads: u64,
    bases: u64,
    species_taxid: i32,
    genus_taxid: i32,
    family_taxid: i32,
}

#[derive(Debug, Clone, Default)]
pub struct SpeciesAlignmentResults {
    pub contig: Option<Taxid>,
    pub read: Option<Taxid>,
}

/// Single BLAST m8 line
#[derive(Debug, Clone, Default)]
pub struct M8Record {
    pub qname: String,       // Read ID
    pub tname: String,       // Base Accession ID (no version suffix)
    pub pident: f64,         // Percent identity
    pub alen: u64,           // Alignment length
    pub mismatch: u64,       // Mismatches
    pub gapopen: u64,        // Gap openings
    pub qstart: u64,         // Query start
    pub qend: u64,           // Query end
    pub tstart: u64,         // Target start
    pub tend: u64,           // Target end
    pub evalue: f64,         // E-value
    pub bitscore: f64,       // Bitscore
    pub qlen: u64,   // Query length: From column in NT; 0 or contig len in NR (unused in NR logic)
    pub slen: u64,   // Subject length: From column in NT; 0 in NR (unused entirely)
}

impl M8Record {
    /// Parse 14-column NT (blastn) output — matches BlastnOutput6NTReader
    pub fn parse_line_nt(line: &str) -> Result<Self> {
        let line = line.trim_end();
        if line.is_empty() {
            return Err(anyhow!("empty line"));
        }

        let mut fields = line.split('\t');

        macro_rules! next {
            () => {
                fields.next().ok_or_else(|| anyhow!("missing field in NT m8 line"))?
            };
        }

        macro_rules! parse_float {
            () => {{
                let s = next!();
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!("invalid float '{}': {}", s, e))?
            }};
        }

        macro_rules! parse_u64 {
            () => {
                next!().parse::<u64>().map_err(|e| anyhow!("invalid u64: {}", e))?
            };
        }

        let qname = next!().to_string();
        let raw_accession = next!(); // e.g., "NC_123456.2"
        let tname = raw_accession
            .split('.')
            .next()
            .unwrap_or(raw_accession)
            .to_string(); // → "NC_123456"

        let pident = parse_float!();
        let alen = parse_u64!();
        let mismatch = parse_u64!();
        let gapopen = parse_u64!();
        let qstart = parse_u64!();
        let qend = parse_u64!();
        let tstart = parse_u64!();
        let tend = parse_u64!();
        let evalue = parse_float!();
        let bitscore = parse_float!();
        let qlen = parse_u64!();
        let slen = parse_u64!();

        // Defensive: Warn on extra columns
        if fields.next().is_some() {
            warn!("Extra columns in NT m8 line (expected 14): {}", line);
        }

        Ok(Self {
            qname,
            tname,
            pident,
            alen,
            mismatch,
            gapopen,
            qstart,
            qend,
            tstart,
            tend,
            evalue,
            bitscore,
            qlen,
            slen,
        })
    }

    /// Parse 12-column NR (blastx) output — matches BlastnOutput6Reader
    pub fn parse_line_nr(line: &str) -> Result<Self> {
        let line = line.trim_end();
        if line.is_empty() {
            return Err(anyhow!("empty line"));
        }

        let mut fields = line.split('\t');

        macro_rules! next {
            () => {
                fields.next().ok_or_else(|| anyhow!("missing field in NR m8 line"))?
            };
        }

        macro_rules! parse_float {
            () => {{
                let s = next!();
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!("invalid float '{}': {}", s, e))?
            }};
        }

        macro_rules! parse_u64 {
            () => {
                next!().parse::<u64>().map_err(|e| anyhow!("invalid u64: {}", e))?
            };
        }

        let qname = next!().to_string();
        let raw_accession = next!(); // e.g., "QIK02963.1"
        let tname = raw_accession
            .split('.')
            .next()
            .unwrap_or(raw_accession)
            .to_string(); // → "QIK02963"

        let pident = parse_float!();
        let alen = parse_u64!();
        let mismatch = parse_u64!();
        let gapopen = parse_u64!();
        let qstart = parse_u64!();
        let qend = parse_u64!();
        let tstart = parse_u64!();
        let tend = parse_u64!();
        let evalue = parse_float!();
        let bitscore = parse_float!();

        // Defensive: Warn on extra columns
        if fields.next().is_some() {
            warn!("Extra columns in NR m8 line (expected 12): {}", line);
        }

        Ok(Self {
            qname,
            tname,
            pident,
            alen,
            mismatch,
            gapopen,
            qstart,
            qend,
            tstart,
            tend,
            evalue,
            bitscore,
            qlen: 0,  // Safe default: Not present/used in NR logic
            slen: 0,  // Safe default: Not present/used anywhere
        })
    }


    pub fn to_tab_string(&self) -> String {
        format!(
            "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{}",
            self.qname,
            self.tname,
            self.pident,
            self.alen,
            self.mismatch,
            self.gapopen,
            self.qstart,
            self.qend,
            self.tstart,
            self.tend,
            self.evalue,
            self.bitscore
        )
    }
}

#[derive(Default)]
pub struct AggBucket {
    pub nonunique_count: u64,
    pub unique_count: u64,
    pub base_count: u64,
    pub sum_percent_identity: f64,
    pub sum_alignment_length: f64,
    pub sum_e_value: f64,
    pub source_count_type: HashSet<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonCount {
    pub tax_id: i32,
    pub tax_level: u8,
    pub genus_taxid: i32,
    pub family_taxid: i32,
    pub count: u64,
    pub nonunique_count: u64,
    pub unique_count: u64,
    pub dcr: f64,
    pub percent_identity: f64,
    pub alignment_length: f64,
    pub e_value: f64,
    pub count_type: String,
    pub base_count: u64,
    pub source_count_type: Option<Vec<String>>,
}


// *******************
// Taxonomy helper functions based on m8 records
// *******************

/// Build a negative taxid for “no-specific-call” at a given level.
///  -100 * level (level 1=species -100, 2=genus -200, 3=family -300)
fn negative_taxid(level: u8) -> i64 {
    -100 * (level as i64)
}


/// Computes the common lineage for a set of hits,
/// Assumes lineages are [species, genus, family] leaf-to-root with possible negatives.
///
/// # Arguments
///
///  * `hits` - list of m8 records
/// * `lineage_map` - derived from taxid DB (taxid -> [species, genus, family])
///
/// # Returns
///
/// (level, consensus_taxid, selected_hits) where level=1 species, 2 genus, 3 family, 0 none;
/// consensus_taxid is the taxid at that level.
pub fn consensus_level(
    hits: &[M8Record],
    lineage_map: &AHashMap<Taxid, Lineage>,
    acc2taxid_map: &Map<Vec<u8>>,
    should_keep: &Arc<impl Fn(&[i32]) -> bool + Send + Sync>,
) -> Result<(u8, i64, Vec<M8Record>)>   {
    let lineages: Vec<Lineage> = hits
        .iter()
        .filter_map(|r| {
            let accession = &r.tname;

            // 1. Try full accession (e.g., "ACC.1")
            let taxid_opt = acc2taxid_map.get(accession.as_bytes());

            // 2. Fallback: strip version → base accession (e.g., "ACC")
            let taxid_u64 = taxid_opt.or_else(|| {
                let base_acc = accession.split('.').next().unwrap_or(accession);
                if base_acc != accession {
                    acc2taxid_map.get(base_acc.as_bytes())
                } else {
                    None
                }
            })
                .ok_or_else(|| anyhow!("Accession not found in acc2taxid map: {}", accession))
                .ok()?;

            let taxid = taxid_u64 as i32;
            if taxid <= 0 {
                return None;
            }

            lineage_map.get(&taxid).cloned()
        })
        .collect();

    if lineages.is_empty() {
        return Ok((0, 0, Vec::new()));
    }

    let mut max_level = 0u8;
    let mut consensus_taxid = 0i64;

    for level in 0..3 {
        let taxids: HashSet<i32> = lineages
            .iter()
            .map(|lin| lin[level])
            .filter(|&t| t > 0)
            .collect();

        if taxids.len() == 1 {
            let agreed = *taxids.iter().next().unwrap() as i64;
            max_level = (level + 1) as u8;
            consensus_taxid = agreed;
        } else {
            break;
        }
    }

    if max_level > 0 {
        let rep_lineage = lineages[0];
        let validated = validate_taxid_lineage(&rep_lineage, consensus_taxid as i32, max_level);
        if !(should_keep)(&validated) {
            return Ok((0, 0, Vec::new()));
        }
    }

    Ok((max_level, consensus_taxid, hits.to_vec()))
}

pub fn parse_m8_batch_to_calls<F>(
    batch: Vec<u8>,
    lineage_map: &AHashMap<Taxid, Lineage>,
    acc2taxid_map: &Map<Vec<u8>>,
    should_keep_filter: &F,
    min_aln_len: u64,
) -> Vec<Vec<u8>>
where
    F: Fn(&[i32]) -> bool + Send + Sync,
{
    let batch_start = Instant::now();
    let mut read_groups: AHashMap<String, Vec<M8Record>> = AHashMap::with_capacity(1024);

    if let Ok(batch_str) = std::str::from_utf8(&batch) {
        for line in batch_str.lines() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            if let Ok(rec) = M8Record::parse_line_nr(line) {
                if rec.alen >= min_aln_len {
                    read_groups.entry(rec.qname.clone()).or_default().push(rec);
                }
            }
        }
    }

    let mut results = Vec::with_capacity(read_groups.len() * 2);

    for (read_id, hits) in read_groups {
        let mut valid_hits: Vec<&M8Record> = Vec::with_capacity(hits.len().min(64));

        for hit in &hits {
            if let Some(taxid_u64) = acc2taxid_map.get(hit.tname.as_bytes()) {
                let taxid = taxid_u64 as i32;
                if taxid <= 0 {
                    continue;
                }

                let lineage = lineage_map.get(&taxid).cloned().unwrap_or([-1i32; 3]);
                if !should_keep_filter(&lineage) {
                    continue;
                }
                valid_hits.push(hit);
            }
        }

        if valid_hits.is_empty() {
            continue;
        }

        let mut best: Option<&M8Record> = None;
        let mut max_bitscore = f64::NEG_INFINITY;

        for &h in &valid_hits {
            if h.bitscore > max_bitscore {
                max_bitscore = h.bitscore;
                best = Some(h);
            }
        }

        let best = match best {
            Some(b) => b,
            None => continue,
        };

        let dedup_line = format!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.3}\n",
            best.qname,
            best.tname,
            best.pident,
            best.alen,
            best.mismatch,
            best.gapopen,
            best.qstart,
            best.qend,
            best.tstart,
            best.tend,
            best.evalue,
            best.bitscore
        );

        let valid_owned: Vec<M8Record> = valid_hits.into_iter().cloned().collect();

        // Wrap should_keep_filter to match consensus_level's requirement
        let wrapped_filter = Arc::new(|l: &[i32]| should_keep_filter(l));
        let (tax_level, _cons_taxid, consensus_hits) =
            consensus_level(&valid_owned, lineage_map, acc2taxid_map, &wrapped_filter)
                .unwrap_or((0, 0, Vec::new()));

        let first_lineage = consensus_hits
            .first()
            .and_then(|h| acc2taxid_map.get(h.tname.as_bytes()))
            .map(|taxid_u64| taxid_u64 as i32)
            .and_then(|taxid| lineage_map.get(&taxid).cloned())
            .unwrap_or([-1i32; 3]);

        let summary_line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            read_id, best.tname, first_lineage[0], first_lineage[1], first_lineage[2], tax_level
        );

        results.push(dedup_line.into_bytes());
        results.push(summary_line.into_bytes());
    }

    results
}

fn parse_read_taxid_from_hitsummary(line: &str) -> Option<(String, i32)> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() >= 3 {
        Some((parts[0].to_string(), parts[2].parse().ok()?))
    } else {
        None
    }
}


pub async fn generate_taxon_count_json_from_m8(
    mut m8_stream_rx:           ReceiverStream<ParseOutput>,
    mut hit_summary_stream_rx:  ReceiverStream<ParseOutput>,
    db_type:                    &str,
    lineage_map:                Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter:         Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    duplicate_clusters:         Arc<DashMap<String, ClusterInfo>>,
    output_tx:                  mpsc::Sender<ParseOutput>,
    concurrency:                usize,
    batch_size_lines:           usize,
) -> Result<()> {
    if concurrency == 0 {
        return Err(anyhow::anyhow!("concurrency must be > 0"));
    }

    let db_type = db_type.to_string();

    // Bounded mpsc channel — hard concurrency cap + backpressure
    let (job_tx, job_rx) = mpsc::channel::<(Vec<String>, Vec<String>)>(concurrency);
    let shared_rx = Arc::new(tokio::sync::Mutex::new(job_rx));

    // Producer: build batches and send to bounded channel (exact same loop as before)
    let producer_handle = tokio::spawn({
        let mut m8_stream = m8_stream_rx;
        let mut hit_stream = hit_summary_stream_rx;
        async move {
            let mut batch_m8 = Vec::with_capacity(batch_size_lines);
            let mut batch_hit = Vec::with_capacity(batch_size_lines);

            loop {
                let m8_item = m8_stream.next().await;
                let hit_item = hit_stream.next().await;

                match (m8_item, hit_item) {
                    (Some(m8), Some(hit)) => {
                        if let Ok(bytes) = m8.to_bytes() {
                            let line = String::from_utf8_lossy(&bytes);
                            let trimmed = line.trim_end().to_string();
                            if !trimmed.is_empty() {
                                batch_m8.push(trimmed);
                            }
                        }
                        if let Ok(bytes) = hit.to_bytes() {
                            let line = String::from_utf8_lossy(&bytes);
                            let trimmed = line.trim_end().to_string();
                            if !trimmed.is_empty() {
                                batch_hit.push(trimmed);
                            }
                        }

                        if batch_m8.len() >= batch_size_lines && batch_hit.len() >= batch_size_lines {
                            if job_tx.send((std::mem::take(&mut batch_m8), std::mem::take(&mut batch_hit))).await.is_err() {
                                break;
                            }
                        }
                    }
                    _ => {
                        // final partial batch
                        if !batch_m8.is_empty() || !batch_hit.is_empty() {
                            let _ = job_tx.send((batch_m8, batch_hit)).await;
                        }
                        break;
                    }
                }
            }
            drop(job_tx);
            Ok::<(), anyhow::Error>(())
        }
    });

    // Fixed worker pool — exactly `concurrency` workers (hard cap)
    let mut workers = Vec::with_capacity(concurrency);
    for _ in 0..concurrency {
        let rx = shared_rx.clone();
        let output_tx = output_tx.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let db_type = db_type.clone();

        let worker = tokio::spawn(async move {
            loop {
                let job = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };

                if let Some((batch_m8, batch_hit)) = job {
                    let mut buckets: AHashMap<Taxid, AggBucket> = AHashMap::new();

                    for (m8_line, hit_line) in batch_m8.into_iter().zip(batch_hit) {
                        let hit_fields: Vec<&str> = hit_line.split('\t').collect();
                        if hit_fields.len() < 10 { continue; }

                        let read_id = hit_fields[0].to_string();
                        let level   = hit_fields.get(1).and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
                        let taxid   = hit_fields.get(2).and_then(|s| s.parse::<i32>().ok()).unwrap_or(0);

                        if taxid <= 0 || level == 0 { continue; }

                        if let Ok(m8) = M8Record::parse_line_nt(&m8_line)
                            .or_else(|_| M8Record::parse_line_nr(&m8_line))
                        {
                            let bucket = buckets.entry(taxid).or_default();

                            bucket.nonunique_count += 1;
                            bucket.unique_count += duplicate_clusters
                                .get(&read_id)
                                .map(|entry| entry.value().size)
                                .unwrap_or(1);

                            bucket.base_count += 1;
                            bucket.sum_percent_identity += m8.pident;
                            bucket.sum_alignment_length += m8.alen as f64;
                            bucket.sum_e_value += m8.evalue;
                            bucket.source_count_type.insert(db_type.clone());
                        }
                    }

                    // Emit results
                    for (taxid, bucket) in buckets {
                        if let Some(lineage) = lineage_map.get(&taxid) {
                            if !should_keep_filter(lineage) { continue; }

                            let dcr = if bucket.nonunique_count > 0 {
                                bucket.unique_count as f64 / bucket.nonunique_count as f64
                            } else { 0.0 };

                            let percent_identity = if bucket.base_count > 0 {
                                bucket.sum_percent_identity / bucket.base_count as f64
                            } else { 0.0 };

                            let alignment_length = if bucket.base_count > 0 {
                                bucket.sum_alignment_length / bucket.base_count as f64
                            } else { 0.0 };

                            let e_value = if bucket.base_count > 0 {
                                bucket.sum_e_value / bucket.base_count as f64
                            } else { 0.0 };

                            let count = TaxonCount {
                                tax_id: taxid,
                                tax_level: 1,
                                genus_taxid: lineage[1],
                                family_taxid: lineage[2],
                                count: bucket.unique_count,
                                nonunique_count: bucket.nonunique_count,
                                unique_count: bucket.unique_count,
                                dcr,
                                percent_identity,
                                alignment_length,
                                e_value,
                                count_type: db_type.clone(),
                                base_count: bucket.base_count,
                                source_count_type: Some(bucket.source_count_type.into_iter().collect()),
                            };

                            let json = serde_json::to_string(&count)? + "\n";
                            let _ = output_tx.send(ParseOutput::Bytes(Arc::new(json.into_bytes()))).await;
                        }
                    }
                } else {
                    break;
                }
            }
            Ok::<(), anyhow::Error>(())
        });
        workers.push(worker);
    }

    // Wait for producer + all workers
    producer_handle.await??;
    for w in workers {
        w.await??;
    }

    info!("Finished taxon count JSON generation for {} (bounded channel, {} workers)", db_type, concurrency);
    Ok(())
}


pub async fn compute_merged_taxon_counts(
    config: Arc<RunConfig>,
    nt_m8_stream: mpsc::Receiver<ParseOutput>,
    nt_hit_summary_stream: mpsc::Receiver<ParseOutput>,
    nt_contig_summary: Vec<ContigSummaryEntry>,

    nr_m8_stream: mpsc::Receiver<ParseOutput>,
    nr_hit_summary_stream: mpsc::Receiver<ParseOutput>,
    nr_contig_summary: Vec<ContigSummaryEntry>,

    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,

    merged_m8_path: PathBuf,
    merged_hitsummary_path: PathBuf,
    merged_taxon_counts_path: PathBuf,
    merged_contig_summary_path: PathBuf,
    nr_alignment_per_read: Arc<DashMap<String, SpeciesAlignmentResults>>,
) -> Result<()> {
    // Write merged m8 and hit summary to disk
    let mut merged_m8_file = BufWriter::new(TokioFile::create(&merged_m8_path).await?);
    let mut merged_hit_file = BufWriter::new(TokioFile::create(&merged_hitsummary_path).await?);

    // Channels for synchronized writing
    let (m8_write_tx, mut m8_write_rx) = mpsc::channel::<Arc<Vec<u8>>>(1024);
    let (hit_write_tx, mut hit_write_rx) = mpsc::channel::<Arc<Vec<u8>>>(1024);

    let m8_write_task = tokio::spawn(async move {
        while let Some(bytes) = m8_write_rx.recv().await {
            merged_m8_file.write_all(&bytes).await?;
            merged_m8_file.write_all(b"\n").await?;
        }
        merged_m8_file.flush().await?;
        Ok::<(), anyhow::Error>(())
    });

    let hit_write_task = tokio::spawn(async move {
        while let Some(bytes) = hit_write_rx.recv().await {
            merged_hit_file.write_all(&bytes).await?;
            merged_hit_file.write_all(b"\n").await?;
        }
        merged_hit_file.flush().await?;
        Ok::<(), anyhow::Error>(())
    });

    // NT pass task
    let nt_m8_write_tx = m8_write_tx.clone();
    let nt_hit_write_tx = hit_write_tx.clone();
    let nr_align_nt = nr_alignment_per_read.clone();
    let nt_task = tokio::spawn(async move {
        let mut nt_m8_stream = ReceiverStream::new(nt_m8_stream);
        let mut nt_hit_stream = ReceiverStream::new(nt_hit_summary_stream);

        loop {
            let m8_item = nt_m8_stream.next().await;
            let hit_item = nt_hit_stream.next().await;

            match (m8_item, hit_item) {
                (Some(m8_item), Some(hit_item)) => {
                    let m8_bytes = if let ParseOutput::Bytes(bytes) = m8_item {
                        bytes
                    } else {
                        return Err(anyhow!("Unexpected item type in NT m8 stream: expected Bytes"));
                    };
                    let hit_bytes = hit_item.to_bytes()?;
                    let hit_line = String::from_utf8_lossy(&hit_bytes);
                    let hit_trimmed = hit_line.trim_end();
                    if hit_trimmed.is_empty() {
                        continue;
                    }

                    let hit_fields: Vec<&str> = hit_trimmed.split('\t').collect();
                    if hit_fields.len() < 10 {
                        warn!("Malformed NT hit line: {}", hit_trimmed);
                        continue;
                    }

                    let read_id = hit_fields[0];

                    let nt_contig = hit_fields.get(9).and_then(|s| s.parse::<Taxid>().ok());
                    let nt_read = hit_fields.get(3).and_then(|s| s.parse::<Taxid>().ok());

                    let nr_align = nr_align_nt.get(read_id);
                    let has_nr_contig = nr_align.as_deref().map_or(false, |a| a.contig.is_some());
                    let has_nr_read = nr_align.as_deref().map_or(false, |a| a.read.is_some());

                    if nt_contig.is_some() || (!has_nr_contig && nt_read.is_some()) {
                        nt_m8_write_tx.send(m8_bytes).await?;

                        let mut hit_with_source = hit_fields.iter().copied().collect::<Vec<&str>>();
                        hit_with_source.push(NT_TAG);
                        nt_hit_write_tx.send(Arc::new(hit_with_source.join("\t").as_bytes().to_vec())).await?;

                        nr_align_nt.remove(read_id);
                    }
                }
                (None, None) => break,
                (Some(_), None) => return Err(anyhow!("NT m8 stream has more items than NT hit summary stream")),
                (None, Some(_)) => return Err(anyhow!("NT hit summary stream has more items than NT m8 stream")),
            }
        }
        Ok::<(), anyhow::Error>(())
    });

    // NR pass task
    let nr_m8_write_tx = m8_write_tx;
    let nr_hit_write_tx = hit_write_tx;
    let nr_align_nr = nr_alignment_per_read;
    let nr_task = tokio::spawn(async move {
        let mut nr_m8_stream = ReceiverStream::new(nr_m8_stream);
        let mut nr_hit_stream = ReceiverStream::new(nr_hit_summary_stream);

        loop {
            let m8_item = nr_m8_stream.next().await;
            let hit_item = nr_hit_stream.next().await;

            match (m8_item, hit_item) {
                (Some(m8_item), Some(hit_item)) => {
                    let hit_bytes = hit_item.to_bytes()?;
                    let hit_line = String::from_utf8_lossy(&hit_bytes);
                    let hit_trimmed = hit_line.trim_end();
                    if hit_trimmed.is_empty() {
                        continue;
                    }

                    let hit_fields: Vec<&str> = hit_trimmed.split('\t').collect();
                    let read_id = hit_fields[0];

                    if nr_align_nr.contains_key(read_id) {
                        let m8_bytes = if let ParseOutput::Bytes(bytes) = m8_item {
                            bytes
                        } else {
                            return Err(anyhow!("Unexpected item type in NR m8 stream: expected Bytes"));
                        };
                        nr_m8_write_tx.send(m8_bytes).await?;

                        let mut hit_with_source = hit_fields.iter().copied().collect::<Vec<&str>>();
                        hit_with_source.push(NR_TAG);
                        nr_hit_write_tx.send(Arc::new(hit_with_source.join("\t").as_bytes().to_vec())).await?;
                    }
                }
                (None, None) => break,
                (Some(_), None) => return Err(anyhow!("NR m8 stream has more items than NR hit summary stream")),
                (None, Some(_)) => return Err(anyhow!("NR hit summary stream has more items than NR m8 stream")),
            }
        }
        Ok::<(), anyhow::Error>(())
    });

    nt_task.await??;
    nr_task.await??;

    m8_write_task.await??;
    hit_write_task.await??;
    info!("Merged alignment files written to disk");

    // Stream merged files back for taxon counting
    let (m8_tx, m8_rx) = mpsc::channel(4096);
    let (hit_tx, hit_rx) = mpsc::channel(4096);

    let m8_path_clone = merged_m8_path.clone();
    let m8_task = tokio::spawn(async move {
        let file = TokioFile::open(m8_path_clone).await?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line).await? > 0 {
            if !line.trim_end().is_empty() {
                let _ = m8_tx.send(ParseOutput::Bytes(Arc::new(line.as_bytes().to_vec()))).await;
            }
            line.clear();
        }
        Ok::<(), anyhow::Error>(())
    });

    let hit_path_clone = merged_hitsummary_path.clone();
    let hit_task = tokio::spawn(async move {
        let file = TokioFile::open(hit_path_clone).await?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line).await? > 0 {
            if !line.trim_end().is_empty() {
                let _ = hit_tx.send(ParseOutput::Bytes(Arc::new(line.as_bytes().to_vec()))).await;
            }
            line.clear();
        }
        Ok::<(), anyhow::Error>(())
    });

    let merged_taxon_counts_path_clone = merged_taxon_counts_path.clone();
    let config_clone = config.clone();

    let (json_tx, json_rx) = mpsc::channel(4096);
    let write_task = tokio::spawn(async move {
        write_byte_stream_to_file(
            &merged_taxon_counts_path_clone,
            ReceiverStream::new(json_rx),
            Some(config_clone.base_buffer_size),
        ).await
    });

    let taxon_count_concurrency = compute_phase_concurrency(
        &config,
        "taxon_counting",
        0.4,      // very light
        4.0,
        64,
        4,
    );

    let taxon_count_batch_size = compute_batch_size(
        None,               // we don't know total lines yet
        220,                // rough avg bytes per m8 + hit line
        150,
        taxon_count_concurrency,
    );

    generate_taxon_count_json_from_m8(
        ReceiverStream::new(m8_rx),
        ReceiverStream::new(hit_rx),
        "merged_NT_NR",
        lineage_map.clone(),
        should_keep_filter.clone(),
        duplicate_clusters.clone(),
        json_tx,
        taxon_count_concurrency,
        taxon_count_batch_size
    )
        .await?;

    m8_task.await??;
    hit_task.await??;
    write_task.await??;
    info!("Merged taxon counts generated");

    // Contig summary merge
    let mut merged_contigs: HashMap<Taxid, Vec<ContigSummaryEntry>> = HashMap::new();
    let mut nt_contig_names: HashSet<String> = HashSet::new();

    for mut entry in nt_contig_summary {
        entry.db_type = "merged_NT_NR".to_string();
        nt_contig_names.insert(entry.contig_name.clone());
        merged_contigs.entry(entry.species_taxid).or_default().push(entry);
    }

    for mut entry in nr_contig_summary {
        if !nt_contig_names.contains(&entry.contig_name) {
            entry.db_type = "merged_NT_NR".to_string();
            merged_contigs.entry(entry.species_taxid).or_default().push(entry);
        }
    }

    let final_contigs: Vec<ContigSummaryEntry> = merged_contigs.into_values().flatten().collect();

    let json = serde_json::to_string_pretty(&final_contigs)?;
    let mut out = BufWriter::new(TokioFile::create(&merged_contig_summary_path).await?);
    out.write_all(json.as_bytes()).await?;
    out.flush().await?;
    info!("Merged contig summary written ({} entries)", final_contigs.len());

    info!("compute_merged_taxon_counts complete");
    Ok(())
}
