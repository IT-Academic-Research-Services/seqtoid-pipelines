use tokio_stream::StreamExt;
use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use std::fs::{self as StdFS};
use std::fs::File as StdFile;
use std::io::{BufReader as StdBufReader};

use std::sync::Arc;
use std::io::{BufWriter as StdBufWriter, Write};

use anyhow::{anyhow, Context, Result};
use log::{info, warn, debug};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use flate2::read::GzDecoder;
use csv::ReaderBuilder;
use sled::{Db, Tree, IVec};
use fst::{MapBuilder};
use needletail::{parse_fastx_file, FastxReader};
use bincode::{encode_into_std_write, decode_from_std_read};
use serde::{Serialize, Deserialize};
use ahash::AHashMap;
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::{mpsc, Semaphore, Mutex};
use crate::config::defs::{Taxid, Lineage, PipelineError, INVALID_CALL_BASE_ID, TaxonSeqLocation};
use crate::utils::blast::{M8Record, AggBucket, TaxonCount};
use crate::utils::streams::ParseOutput;

const SPECIES_NON_SPECIFIC: Taxid = -100;
const GENUS_NON_SPECIFIC: Taxid = -200;
const FAMILY_NON_SPECIFIC: Taxid = -300;

const NULL_TAXID: Taxid = -1;
const NULL_LINEAGE: Lineage = [NULL_TAXID; 3];

// *******************
// DB creation functions
// *******************

/// Builds a taxid-lineages map from taxdump.tar.gz files and serializes to bincode.
/// Typical location, sourced from NCBI:
/// ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
///
/// # Arguments
///
///  * `nodes_path` - path to the nodes.dmp file
/// * `merged_path` - path to the merged.dmp file
/// * `output_bincode_path` - output bincode path
///
/// # Returns
///
/// Result
pub async fn build_taxid_lineages_db(
    nodes_path: &PathBuf,
    merged_path: &PathBuf,
    output_bincode_path: &PathBuf,
) -> Result<()> {
    // Pre-allocate maps for efficiency (NCBI taxonomy ~2.7M entries)
    let mut parent_map: HashMap<Taxid, Taxid> = HashMap::with_capacity(3_000_000);
    let mut rank_map: HashMap<Taxid, String> = HashMap::with_capacity(3_000_000);

    let nodes_file = File::open(nodes_path).await?;
    let mut nodes_reader = BufReader::new(nodes_file).lines();
    while let Some(line) = nodes_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 3 { continue; }
        let Ok(taxid) = parts[0].parse::<Taxid>() else { continue; };
        let Ok(parent) = parts[1].parse::<Taxid>() else { continue; };
        let rank = parts[2].to_string();
        if taxid != 1 {
            parent_map.insert(taxid, parent);
            rank_map.insert(taxid, rank);
        }
    }
    info!("Loaded {} nodes", parent_map.len());

    let mut merged_map: HashMap<Taxid, Taxid> = HashMap::with_capacity(100_000);
    let merged_file = File::open(merged_path).await?;
    let mut merged_reader = BufReader::new(merged_file).lines();
    while let Some(line) = merged_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 2 { continue; }
        let Ok(old) = parts[0].parse::<Taxid>() else { continue; };
        let Ok(new) = parts[1].parse::<Taxid>() else { continue; };
        merged_map.insert(old, new);
    }
    info!("Loaded {} merges", merged_map.len());

    let mut resolved_parents = parent_map.clone();
    let mut resolved_ranks = rank_map.clone();

    for (&taxid, _) in parent_map.iter() {
        let mut current = taxid;
        while let Some(&merged) = merged_map.get(&current) {
            current = merged;
        }
        if current != taxid {
            if let Some(&parent) = parent_map.get(&current) {
                resolved_parents.insert(taxid, parent);
            }
            if let Some(rank) = rank_map.get(&current) {
                resolved_ranks.insert(taxid, rank.clone());
            }
        }
    }

    parent_map = resolved_parents;
    rank_map = resolved_ranks;

    let mut lineages: HashMap<Taxid, Lineage> = HashMap::with_capacity(parent_map.len());
    for (&taxid, _) in parent_map.iter() {
        let mut path = Vec::with_capacity(10);
        let mut seen = HashSet::new();
        let mut current = taxid;
        while current != 1 && seen.insert(current) {
            path.push(current);
            if let Some(&parent) = parent_map.get(&current) {
                current = parent;
            } else {
                break;
            }
        }
        path.reverse();

        let mut species = -1;
        let mut genus = -1;
        let mut family = -1;

        for &t in &path {
            let rank = rank_map.get(&t).map(|s| s.as_str()).unwrap_or("no_rank");
            if rank == "species" && species == -1 {
                species = t;
            } else if rank == "genus" && genus == -1 {
                genus = t;
            } else if rank == "family" && family == -1 {
                family = t;
            }
            if species != -1 && genus != -1 && family != -1 {
                break;
            }
        }

        lineages.insert(taxid, [species, genus, family]);
    }

    info!("Built {} lineages", lineages.len());

    let mut file = StdFile::create(output_bincode_path)?;
    encode_into_std_write(&lineages, &mut file, bincode::config::standard())
        .map_err(|e| anyhow!("Bincode encode failed: {}", e))?;

    info!("Saved to {}", output_bincode_path.display());
    Ok(())
}

pub async fn load_taxid_lineages_db(bincode_path: &PathBuf) -> Result<Arc<AHashMap<Taxid, Lineage>>> {
    let mut file = StdFile::open(bincode_path)?;
    let std_map: HashMap<Taxid, Lineage> = decode_from_std_read(&mut file, bincode::config::standard())
        .map_err(|e| anyhow!("Bincode decode failed: {}", e))?;
    let map = AHashMap::from_iter(std_map);
    Ok(Arc::new(map))
}



#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct AccTaxEntry {
    pub acc: String,
    pub taxid: Taxid,
}

/// Builds a accession2taxid FST DB from a known accession2taxid file
/// Skips accessions not present in NT/NR
/// Typical location, sourced from NCBI:
/// ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
/// ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
///
/// # Arguments
///
///  * `gz_path` - accession2taxid.gz
///  * `nt_file` - Path to NT FASTA file
///  * `nr_file` - Path to NR FASTA file
///  * `db_path` -  output db path
///
/// # Returns
/// Result
pub async fn build_accession2taxid_db(
    mapping_files: &[PathBuf],
    nt_file: Option<&PathBuf>,
    nr_file: Option<&PathBuf>,
    db_path: &PathBuf,
) -> Result<()> {
    let mut accession_bases: HashSet<String> = HashSet::new();

    let mut scan_tasks = Vec::new();

    if let Some(p) = nt_file.cloned() {
        scan_tasks.push(tokio::task::spawn_blocking(move || {
            let mut bases = HashSet::new();
            scan_fasta_bases(&p, &mut bases, "NT");
            bases
        }));
    }
    if let Some(p) = nr_file.cloned() {
        scan_tasks.push(tokio::task::spawn_blocking(move || {
            let mut bases = HashSet::new();
            scan_fasta_bases(&p, &mut bases, "NR");
            bases
        }));
    }

    for task in scan_tasks {
        let bases = task.await?;
        accession_bases.extend(bases);
    }

    let accession_bases = Arc::new(accession_bases);
    let (tx, mut rx) = mpsc::channel(100);
    let concurrency = std::cmp::max(1, num_cpus::get() / 2);
    let (job_tx, job_rx) = mpsc::channel::<PathBuf>(concurrency);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let mut worker_handles = Vec::new();

    for i in 0..concurrency {
        let job_rx = Arc::clone(&job_rx);
        let tx = tx.clone();
        let accession_bases = Arc::clone(&accession_bases);

        let handle = tokio::spawn(async move {
            loop {
                let mapping_file = {
                    let mut rx = job_rx.lock().await;
                    rx.recv().await
                };

                let mapping_file = match mapping_file {
                    Some(f) => f,
                    None => break,
                };

                let filename = mapping_file.file_name().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
                info!("Worker {} reading mapping file: {}", i, filename);

                let tx = tx.clone();
                let accession_bases = Arc::clone(&accession_bases);
                let res = tokio::task::spawn_blocking(move || -> Result<()> {
                    let file = StdFile::open(&mapping_file)?;
                    let gz = GzDecoder::new(StdBufReader::new(file));
                    let mut rdr = ReaderBuilder::new()
                        .delimiter(b'\t')
                        .has_headers(false)
                        .from_reader(gz);

                    let mut batch = Vec::with_capacity(1000);
                    for result in rdr.records() {
                        let record = result.map_err(|e| anyhow!("Parse error in {}: {}", filename, e))?;

                        let full_acc = record[0].to_string();
                        let acc_base = full_acc.split('.').next().unwrap_or(&full_acc).to_string();

                        if !accession_bases.contains(&acc_base) {
                            continue;
                        }

                        let taxid_idx = if record.len() == 2 { 1 } else { 2 };
                        let taxid: Taxid = record.get(taxid_idx)
                            .ok_or_else(|| anyhow!("Missing taxid column in {}", filename))?
                            .parse()
                            .map_err(|_| anyhow!("Invalid taxid in {}: {}", filename, record.get(taxid_idx).unwrap_or("<missing>")))?;

                        batch.push(AccTaxEntry { acc: acc_base, taxid });  // ← key = base only
                        if batch.len() >= 1000 {
                            tx.blocking_send(std::mem::replace(&mut batch, Vec::with_capacity(1000)))
                                .map_err(|_| anyhow!("Channel closed"))?;
                        }
                    }
                    if !batch.is_empty() {
                        tx.blocking_send(batch).map_err(|_| anyhow!("Channel closed"))?;
                    }
                    Ok(())
                }).await?;

                if let Err(e) = res {
                    return Err(e);
                }
            }
            Ok::<(), anyhow::Error>(())
        });
        worker_handles.push(handle);
    }

    // Producer for jobs
    let mapping_files_vec = mapping_files.to_vec();
    tokio::spawn(async move {
        for mapping_file in mapping_files_vec {
            if job_tx.send(mapping_file).await.is_err() {
                break;
            }
        }
    });

    drop(tx);

    let mut entries = Vec::new();
    while let Some(batch) = rx.recv().await {
        entries.extend(batch);
    }

    for handle in worker_handles {
        handle.await??;
    }

    // No need to sort — we use HashMap dedup below
    let mut acc_to_taxid: HashMap<String, Taxid> = HashMap::new();
    let mut total_mapped = 0;

    for entry in entries {
        // Dedup: keep first seen (or warn on conflict)
        if let Some(existing) = acc_to_taxid.insert(entry.acc.clone(), entry.taxid) {
            if existing != entry.taxid {
                warn!("Taxid conflict for base {}: {} vs {}", entry.acc, existing, entry.taxid);
            }
        } else {
            total_mapped += 1;
        }

        if total_mapped % 1_000_000 == 0 {
            info!("\t{}M base accessions mapped", total_mapped / 1_000_000);
        }
    }

    // Build FST with base accessions as keys
    let mut builder = MapBuilder::memory();
    let mut sorted_keys: Vec<String> = acc_to_taxid.keys().cloned().collect();
    sorted_keys.sort_unstable();  // deterministic order

    for acc_base in sorted_keys {
        if let Some(&taxid) = acc_to_taxid.get(&acc_base) {
            builder.insert(acc_base.as_bytes(), taxid as u64)?;
        }
    }

    let map_bytes = builder.into_inner()?;
    let mut out_file = StdFile::create(db_path)?;
    out_file.write_all(&map_bytes)?;
    info!("Built FST acc2taxid DB with {} entries at {}", total_mapped, db_path.display());

    Ok(())
}

/// scans a FASTA file for base accessions (i.e. no version at end) and inserts
pub fn scan_fasta_bases(
    path: &PathBuf,
    bases: &mut HashSet<String>,
    name: &str,
) {
    let mut reader = parse_fastx_file(path)
        .unwrap_or_else(|e| panic!("Failed to open {} FASTA {}: {}", name, path.display(), e));

    let mut count = 0;
    while let Some(record) = reader.next() {
        let rec = record.expect("FASTA parse error");
        if let Some(acc_base) = std::str::from_utf8(rec.id())
            .ok()
            .and_then(|id| id.split('.').next())
        {
            bases.insert(acc_base.to_string());
            count += 1;
        }
    }

    info!(
        "Extracted {} unique accession bases from {}: {}",
        bases.len(),
        name,
        path.display()
    );
}

/// Builds filter for taxa to keep from the three lists input
///
/// # Arguments
///
///  * `deuterostome_path` - path to deuterostome filter list
///  * `taxon_whitelist_path` - Always-keep list
///  * `taxon_blacklist_path` - Always-discard list
///
/// # Returns
/// RResult<Closure, anyhow::Error> – an Ok that contains a zero-cost, heap-allocated closure
/// which, when called with a slice of taxonomy IDs, tells you whether those IDs are allowed
/// to stay in the downstream aggregation step.
pub async fn build_should_keep_filter(
    deuterostome_path: Option<PathBuf>,
    taxon_whitelist_path: Option<PathBuf>,
    taxon_blacklist_path: Option<PathBuf>,
) -> Result<impl Fn(&[i32]) -> bool + Send + Sync + 'static> {
    let mut taxids_to_remove: HashSet<i32> = HashSet::from([9605, 9606]);

    if let Some(path) = taxon_blacklist_path {
        taxids_to_remove.extend(read_file_into_set(&path).await?);
    }
    if let Some(path) = deuterostome_path {
        taxids_to_remove.extend(read_file_into_set(&path).await?);
    }

    let taxids_to_keep: Option<HashSet<i32>> = if let Some(path) = taxon_whitelist_path {
        Some(read_file_into_set(&path).await?)
    } else {
        None
    };

    Ok(move |hits: &[i32]| {
        let is_blacklisted = hits.iter().any(|&t| t > 0 && taxids_to_remove.contains(&t));
        let is_whitelisted = if let Some(keep) = &taxids_to_keep {
            hits.iter().any(|&t| t > 0 && keep.contains(&t))
        } else {
            true
        };
        is_whitelisted && !is_blacklisted
    })
}


///Reading helper function for above filter list files
///
/// # Arguments
///
///  * `path` - path to  filter list
///
/// # Returns
/// RResult of hash set of taxon id numbers
pub async fn read_file_into_set(path: &PathBuf) -> Result<HashSet<i32>> {
    let file = File::open(path).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let mut set = HashSet::new();
    while let Some(line) = lines.next_line().await? {
        if let Ok(taxid) = line.trim().parse::<i32>() {
            set.insert(taxid);
        }
    }
    Ok(set)
}


/// Retrieves the top hit from an m8 stream
/// Ranking
/// 1. Higher AAI = pident * alen / qlen
/// 2. Lower evalue
/// 3. Higher bitscore
/// # Arguments
/// * `path` - path to hit sumamry file.
/// # Returns
/// hash map of HitSummaries
pub async fn get_top_m8_nt(
    mut input: ReceiverStream<ParseOutput>,
    output_tx: mpsc::Sender<ParseOutput>,
    concurrency: usize,
    batch_size: usize,   // Lines per batch; aim for 10-50MB/batch
) -> Result<()> {
    if concurrency == 0 {
        return Err(anyhow!("concurrency must be > 0"));
    }

    // Bounded mpsc channel — hard concurrency cap + backpressure
    let (job_tx, job_rx) = mpsc::channel::<Vec<String>>(concurrency);
    let shared_rx = Arc::new(tokio::sync::Mutex::new(job_rx));

    let mut worker_handles = Vec::with_capacity(concurrency);
    for i in 0..concurrency {
        let rx = shared_rx.clone();
        let out_tx = output_tx.clone();
        let handle = tokio::spawn(async move {
            loop {
                let batch = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };

                let batch_lines = match batch {
                    Some(b) => b,
                    None => break,
                };

                let mut hits_by_read: AHashMap<String, Vec<M8Record>> = AHashMap::new();
                for line in batch_lines {
                    match M8Record::parse_line_nt(&line) {
                        Ok(m8) => {
                            hits_by_read.entry(m8.qname.clone()).or_default().push(m8);
                        }
                        Err(e) => {
                            log::warn!("Worker {} failed to parse NT m8 line: {} — {}", i, e, line);
                        }
                    }
                }

                for (_read_id, mut hits) in hits_by_read {
                    if hits.is_empty() {
                        continue;
                    }
                    // Sort: Descending alen, ascending mismatch+gapopen, ascending evalue, descending bitscore
                    hits.sort_by(|a, b| {
                        b.alen.cmp(&a.alen)
                            .then_with(|| a.mismatch.cmp(&b.mismatch))
                            .then_with(|| a.gapopen.cmp(&b.gapopen))
                            .then_with(|| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal))
                            .then_with(|| b.bitscore.partial_cmp(&a.bitscore).unwrap_or(std::cmp::Ordering::Equal))
                    });
                    let best = hits.into_iter().next().expect("Non-empty vec after check");
                    let line = best.to_tab_string() + "\n";
                    if out_tx.send(ParseOutput::Bytes(Arc::new(line.into_bytes()))).await.is_err() {
                        return Err(anyhow!("Output send failed in worker {}", i));
                    }
                }
            }
            Ok::<(), anyhow::Error>(())
        });
        worker_handles.push(handle);
    }

    // Producer: build batches and send to bounded channel
    let mut batch = Vec::with_capacity(batch_size);
    while let Some(item) = input.next().await {
        if let ParseOutput::Bytes(bytes) = item {
            let line = String::from_utf8_lossy(&bytes).trim_end().to_string();
            if line.is_empty() {
                continue;
            }
            batch.push(line);

            if batch.len() >= batch_size {
                if job_tx.send(std::mem::take(&mut batch)).await.is_err() {
                    break;
                }
            }
        }
    }

    if !batch.is_empty() {
        let _ = job_tx.send(batch).await;
    }
    drop(job_tx);

    for handle in worker_handles {
        handle.await.map_err(|e| anyhow!("Worker join failed: {}", e))??;
    }

    Ok(())
}

#[derive(Default)]
struct MergedHsp {
    total_alen: u64,
    sum_pident: f64,
    evalue: f64,
    bitscore: f64,
    hsp_count: u64,
    representative: M8Record,
}


pub async fn get_top_m8_nr(
    mut input: ReceiverStream<ParseOutput>,
    output_tx: mpsc::Sender<ParseOutput>,
    concurrency: usize,
    batch_size: usize,   // Lines per batch; aim for 10-50MB/batch
) -> Result<()> {
    if concurrency == 0 {
        return Err(anyhow!("concurrency must be > 0"));
    }

    // Bounded mpsc channel — hard concurrency cap + backpressure
    let (job_tx, job_rx) = mpsc::channel::<Vec<String>>(concurrency);
    let shared_rx = Arc::new(tokio::sync::Mutex::new(job_rx));

    let mut worker_handles = Vec::with_capacity(concurrency);
    for i in 0..concurrency {
        let rx = shared_rx.clone();
        let out_tx = output_tx.clone();
        let handle = tokio::spawn(async move {
            loop {
                let batch = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };

                let batch_lines = match batch {
                    Some(b) => b,
                    None => break,
                };

                // Merge HSPs per (contig, subject)
                let mut merged: AHashMap<(String, String), MergedHsp> = AHashMap::new();
                for line in batch_lines {
                    match M8Record::parse_line_nr(&line) {
                        Ok(m8) => {
                            let key = (m8.qname.clone(), m8.tname.clone());
                            let entry = merged.entry(key).or_default();
                            entry.total_alen += m8.alen;
                            entry.sum_pident += m8.pident * m8.alen as f64;
                            entry.evalue = entry.evalue.min(m8.evalue);
                            entry.bitscore = entry.bitscore.max(m8.bitscore);
                            entry.hsp_count += 1;
                            // Update representative if better bitscore
                            if m8.bitscore > entry.representative.bitscore {
                                entry.representative = m8;
                            }
                        }
                        Err(e) => {
                            log::warn!("Worker {} failed to parse NR m8 line: {} — {}", i, e, line);
                        }
                    }
                }

                // Per contig, pick best subject
                let mut best_per_contig: AHashMap<String, MergedHsp> = AHashMap::new();
                for ((contig, _subject), merged_hsp) in merged {
                    let current = best_per_contig.entry(contig).or_default();
                    if merged_hsp.bitscore > current.bitscore ||
                        (merged_hsp.bitscore == current.bitscore && merged_hsp.evalue < current.evalue) {
                        *current = merged_hsp;
                    }
                }

                // Output reconstructed lines
                for (_contig, best) in best_per_contig {
                    let rep = &best.representative;
                    let avg_pident = if best.total_alen > 0 {
                        best.sum_pident / best.total_alen as f64
                    } else {
                        rep.pident
                    };
                    // Reconstruct line (12-column NR format: no qlen/slen)
                    let line = format!(
                        "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}\n",
                        rep.qname, rep.tname, avg_pident, best.total_alen,
                        rep.mismatch, rep.gapopen, rep.qstart, rep.qend,
                        rep.tstart, rep.tend, best.evalue, best.bitscore
                    );
                    if out_tx.send(ParseOutput::Bytes(Arc::new(line.into_bytes()))).await.is_err() {
                        return Err(anyhow!("Output send failed in worker {}", i));
                    }
                }
            }
            Ok::<(), anyhow::Error>(())
        });
        worker_handles.push(handle);
    }

    // Producer: build batches and send to bounded channel
    let mut batch = Vec::with_capacity(batch_size);
    while let Some(item) = input.next().await {
        if let ParseOutput::Bytes(bytes) = item {
            let line = String::from_utf8_lossy(&bytes).trim_end().to_string();
            if line.is_empty() {
                continue;
            }
            batch.push(line);

            if batch.len() >= batch_size {
                if job_tx.send(std::mem::take(&mut batch)).await.is_err() {
                    break;
                }
            }
        }
    }

    if !batch.is_empty() {
        let _ = job_tx.send(batch).await;
    }
    drop(job_tx);

    for handle in worker_handles {
        handle.await.map_err(|e| anyhow!("Worker join failed: {}", e))??;
    }

    Ok(())
}



// *******************
// DB access functions
// *******************

/// Validates and cleans a taxonomic lineage
/// /// # Arguments
/// * `lineage`  lineage from taxdump (species, genus, family, ... root). May contain 0s or gaps.
/// * `hit_taxid`   - The taxid at the level where consensus was reached
/// * `hit_level`   - 1: species, 2: genus, 3: family
///
///
/// # Returns
/// Validayte lineage vector with negtive taxids by converntion for no hit at that level
pub fn validate_taxid_lineage(
    lineage: &[i32], //
    hit_taxid: Taxid,
    hit_level: u8
) -> Vec<i32> {
    const INVALID_BASE: i32 = -2_000_000_000; //  base for artificial negative taxids

    let mut cleaned = lineage.to_vec();

    // Invalidate all levels below the consensus level
    // If hit_level == 0 → no consensus → invalidate everything
    let invalidate_up_to = if hit_level == 0 {
        cleaned.len()
    } else {
        hit_level.saturating_sub(1) as usize
    };

    for level in 0..invalidate_up_to {
        cleaned[level] = INVALID_BASE - (level as i32 + 1) * 100;
    }

    //  From the consensus level upward, fill missing/invalid entries
    // with artificial negative taxids based on parent
    if hit_level > 0 && hit_level <= cleaned.len() as u8 {
        let start_level = (hit_level as usize).saturating_sub(1); // inclusive
        let mut parent = hit_taxid;

        for level in start_level..cleaned.len() {
            if cleaned[level] <= 0 {
                // missing or invalid. so artificial negative ID
                cleaned[level] = INVALID_BASE - (level as i32 + 1) * 100 - parent;
            }
            parent = cleaned[level];
        }
    }

    cleaned
}


pub fn get_valid_lineage(
    hits_by_read_id: &AHashMap<String, (Taxid, u8)>,     // contig_id → (taxid, level)
    lineage_map: &Arc<AHashMap<Taxid, Lineage>>,
    read_id: &str,
) -> Lineage {
    let (hit_taxid, hit_level) = hits_by_read_id
        .get(read_id)
        .copied()
        .unwrap_or((-1, 255)); // 255 = invalid level

    if hit_taxid <= 0 {
        return NULL_LINEAGE;
    }

    if hit_level == 1 {

        lineage_map
            .get(&hit_taxid)
            .copied()
            .unwrap_or(NULL_LINEAGE)
    } else {

        [SPECIES_NON_SPECIFIC, GENUS_NON_SPECIFIC, FAMILY_NON_SPECIFIC]
    }
}

pub fn generate_locator_work(
    records: Arc<Vec<(String, String)>>, // (header_without_>, sequence)
    taxid_field: String,
    hit_type: String,
    output_fa: PathBuf,
    output_json: PathBuf,
) -> Result<()> {
    if records.is_empty() {
        std::fs::write(&output_fa, b"")?;
        std::fs::write(&output_json, b"[]")?;
        debug!("Empty FASTA/JSON written for taxid_field: {}", taxid_field);
        return Ok(());
    }

    // -----------------------------------------------------------------------
    // 1. Build sort indices using the same taxid extraction logic as Python
    // -----------------------------------------------------------------------
    let mut indices: Vec<usize> = (0..records.len()).collect();

    indices.sort_unstable_by_key(|&i| {
        get_taxid(&records[i].0, &taxid_field).unwrap_or(-1)
    });

    // -----------------------------------------------------------------------
    // 2. Prepare FASTA writer (large buffer, same as before)
    // -----------------------------------------------------------------------
    let std_file = StdFile::create(&output_fa)
        .map_err(|e| anyhow!("Failed to create FASTA file {}: {}", output_fa.display(), e))?;

    let mut writer = StdBufWriter::with_capacity(8 * 1024 * 1024, std_file);

    // -----------------------------------------------------------------------
    // 3. Write sorted FASTA + collect byte ranges
    // -----------------------------------------------------------------------
    let mut locations: Vec<TaxonSeqLocation> = Vec::with_capacity(100_000.min(records.len() / 10));
    let mut current_taxid: Taxid = -1;
    let mut first_byte: u64 = 0;
    let mut pos: u64 = 0;

    for &idx in &indices {
        let (header, seq) = &records[idx];
        let taxid = get_taxid(header, &taxid_field).unwrap_or(-1);

        // Close previous block if taxid changed
        if taxid != current_taxid && current_taxid != -1 {
            locations.push(TaxonSeqLocation {
                taxid: current_taxid,
                first_byte,
                last_byte: pos.saturating_sub(1),
                hit_type: hit_type.clone(),
            });
            first_byte = pos;
        }
        current_taxid = taxid;

        // Write header
        writer.write_all(b">")?;
        writer.write_all(header.as_bytes())?;
        writer.write_all(b"\n")?;
        pos += 1 + header.len() as u64 + 1;

        // Write sequence
        writer.write_all(seq.as_bytes())?;
        writer.write_all(b"\n")?;
        pos += seq.len() as u64 + 1;
    }

    // Final block
    if current_taxid != -1 {
        locations.push(TaxonSeqLocation {
            taxid: current_taxid,
            first_byte,
            last_byte: pos.saturating_sub(1),
            hit_type: hit_type.clone(),
        });
    }

    writer.flush()
        .map_err(|e| anyhow!("Failed to flush FASTA writer for {}: {}", output_fa.display(), e))?;

    // -----------------------------------------------------------------------
    // 4. Write JSON locator (same format as Python)
    // -----------------------------------------------------------------------
    let json_bytes = serde_json::to_vec(&locations)
        .map_err(|e| anyhow!("Failed to serialize JSON for {}: {}", output_json.display(), e))?;

    std::fs::write(&output_json, json_bytes)
        .map_err(|e| anyhow!("Failed to write JSON {}: {}", output_json.display(), e))?;

    debug!(
        "generate_locator_work complete: {} contigs → {} FASTA bytes, {} locations",
        records.len(),
        pos,
        locations.len()
    );

    Ok(())
}

/// Matches the original Python splitting logic as closely as possible
/// Returns -1 on failure (same fallback as your current version)
pub fn get_taxid(header: &str, field: &str) -> Result<i32> {
    // Remove leading '>' if present
    let cleaned = header.trim_start_matches('>');

    // Build delimiter exactly like Python: ":field:"
    let delimiter = format!(":{}:", field);

    // Split exactly on that delimiter (like Python split(f":{taxid_field}:"))
    let parts: Vec<&str> = cleaned.split(&delimiter).collect();

    if parts.len() <= 1 {
        // Field not found
        return Ok(-1);
    }

    // Take the part immediately after the delimiter
    let after = parts[1];

    // Take until next ':' (like Python parts[1].split(":")[0])
    let taxid_str = match after.split(':').next() {
        Some(s) => s.trim(),
        None => return Ok(-1),
    };

    // Parse (same as Python — fail silently to -1)
    let taxid = taxid_str
        .parse::<i32>()
        .unwrap_or_else(|e| {
            debug!("Failed to parse taxid '{}': {}", taxid_str, e);
            -1
        });

    Ok(taxid)
}


fn get_taxid_field_num(sample_headers: &[String], taxid_field: &str) -> i32 {
    for header in sample_headers.iter().take(100) {  // look at first 100 headers
        let cleaned = header.trim_start_matches('>');
        let parts: Vec<&str> = cleaned.split(':').collect();
        if let Some(idx) = parts.iter().position(|&p| p == taxid_field) {
            info!("Found taxid field '{}' at position {}", taxid_field, idx + 1);
            return (idx + 1) as i32;
        }
    }
    warn!("Taxid field '{}' not found in first 100 headers", taxid_field);
    -1
}


pub fn combine_taxon_loc_json(input_jsons: Vec<PathBuf>, output_json: PathBuf) -> Result<()> {
    let mut combined: Vec<TaxonSeqLocation> = Vec::new();
    for path in input_jsons {
        let data = StdFS::read(&path)?;
        let mut locs: Vec<TaxonSeqLocation> = serde_json::from_slice(&data)?;
        combined.append(&mut locs);
    }
    let json_data = serde_json::to_vec(&combined)?;
    StdFS::write(output_json, json_data)?;
    Ok(())
}

