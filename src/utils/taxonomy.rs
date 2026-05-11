use tokio_stream::StreamExt;
use std::path::PathBuf;
use std::collections::{HashMap, HashSet};
use std::fs::{self as StdFS};
use std::fs::File as StdFile;
use std::io::{BufReader as StdBufReader};

use std::sync::Arc;
use std::io::{BufWriter as StdBufWriter, Write};

use anyhow::{anyhow, Result};
use log::{info, warn, debug};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use flate2::read::GzDecoder;
use csv::ReaderBuilder;
use fst::{MapBuilder};
use needletail::parse_fastx_file;
use bincode::{encode_into_std_write, decode_from_std_read};
use ahash::AHashMap;
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::{mpsc, Mutex};
use crate::config::defs::{Taxid, Lineage, TaxonSeqLocation};
use crate::utils::blast::M8Record;
use crate::utils::streams::ParseOutput;
use bytes::Bytes;

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

    while let Some(record) = reader.next() {
        let rec = record.expect("FASTA parse error");
        if let Some(acc_base) = std::str::from_utf8(rec.id())
            .ok()
            .and_then(|id| id.split('.').next())
        {
            bases.insert(acc_base.to_string());
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
                    if out_tx.send(ParseOutput::Bytes(Bytes::from(line.into_bytes()))).await.is_err() {
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
                    if out_tx.send(ParseOutput::Bytes(Bytes::from(line.into_bytes()))).await.is_err() {
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
/// Matches the CZI lineage.validate_taxid_lineage  for 3-rank case
/// (species, genus, family).
///
/// - Levels below the consensus hit_level get artificial negative taxids.
/// - Levels at/above consensus keep real taxid if present, else get negative fallback
/// # Arguments
/// * `lineage`  lineage from taxdump (species, genus, family, ... root). May contain 0s or gaps.
/// * `hit_taxid`   - The taxid at the level where consensus was reached
/// * `hit_level`   - 1: species, 2: genus, 3: family
///
///
/// # Returns
/// Validayte lineage vector with negtive taxids by converntion for no hit at that level
pub fn validate_taxid_lineage(
    lineage: &[i32],        // full lineage from DB for the hit_taxid
    hit_taxid: Taxid,       // the taxid at which consensus was reached
    hit_level: u8,          // 1=species, 2=genus, 3=family, 0=no consensus
) -> Vec<i32> {
    let mut cleaned = lineage.to_vec();   // always length 3 in our usage

    if hit_level == 0 {
        // No consensus at any level → invalidate everything
        cleaned[0] = SPECIES_NON_SPECIFIC;   // -100
        cleaned[1] = GENUS_NON_SPECIFIC;     // -200
        cleaned[2] = FAMILY_NON_SPECIFIC;    // -300
        return cleaned;
    }

    // Invalidate levels below the consensus level
    // hit_level=1 → invalidate nothing (species is the consensus)
    // hit_level=2 → invalidate species (index 0)
    // hit_level=3 → invalidate species + genus (indices 0,1)
    for level in 0..(hit_level as usize).saturating_sub(1) {
        cleaned[level] = match level {
            0 => SPECIES_NON_SPECIFIC,
            1 => GENUS_NON_SPECIFIC,
            2 => FAMILY_NON_SPECIFIC,
            _ => -100 * (level as i32 + 1),
        };
    }

    // For the consensus level and above, keep real values or fill gaps with negatives
    // (in practice for 3-rank we rarely have gaps above consensus)
    let start = (hit_level as usize).saturating_sub(1);
    let mut parent = hit_taxid;

    for level in start..cleaned.len() {
        if cleaned[level] <= 0 {
            // missing or invalid → artificial negative based on parent
            cleaned[level] = -100 * (level as i32 + 1) - parent;
        }
        parent = cleaned[level];
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


#[allow(dead_code)]
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::Arc;
    use tempfile::tempdir;

    fn taxon_loc(taxid: i32, first_byte: u64, last_byte: u64, hit_type: &str) -> TaxonSeqLocation {
        TaxonSeqLocation {
            taxid,
            first_byte,
            last_byte,
            hit_type: hit_type.to_string(),
        }
    }

    #[test]
    fn test_validate_taxid_lineage_no_consensus() {
        let lineage = [10, 20, 30];
        let cleaned = validate_taxid_lineage(&lineage, 999, 0);
        assert_eq!(cleaned, [-100, -200, -300]);
    }

    #[test]
    fn test_validate_taxid_lineage_species_consensus_keeps_values() {
        let lineage = [10, 20, 30];
        let cleaned = validate_taxid_lineage(&lineage, 10, 1);
        assert_eq!(cleaned, [10, 20, 30]);
    }

    #[test]
    fn test_validate_taxid_lineage_genus_consensus_invalidates_species() {
        let lineage = [0, 20, 30];
        let cleaned = validate_taxid_lineage(&lineage, 20, 2);
        assert_eq!(cleaned, [-100, 20, 30]);
    }

    #[test]
    fn test_validate_taxid_lineage_family_consensus_invalidates_species_and_genus() {
        let lineage = [0, 0, 30];
        let cleaned = validate_taxid_lineage(&lineage, 30, 3);
        assert_eq!(cleaned, [-100, -200, 30]);
    }

    #[test]
    fn test_get_taxid_parses_field_and_ignores_missing_field() {
        let header = "sample:species_nt:123:readA";
        assert_eq!(get_taxid(header, "species_nt").unwrap(), 123);
        assert_eq!(get_taxid(header, "genus_nt").unwrap(), -1);
    }

    #[test]
    fn test_get_taxid_handles_leading_gt_and_bad_number() {
        let header = ">sample:species_nt:not_a_number:readA";
        assert_eq!(get_taxid(header, "species_nt").unwrap(), -1);
    }

    #[test]
    fn test_get_valid_lineage_species_hit_uses_map() {
        let mut hits = AHashMap::new();
        hits.insert("read1".to_string(), (42, 1));

        let mut lineage_map = AHashMap::new();
        lineage_map.insert(42, [11, 22, 33]);

        let lineage_map = Arc::new(lineage_map);
        let got = get_valid_lineage(&hits, &lineage_map, "read1");
        assert_eq!(got, [11, 22, 33]);
    }

    #[test]
    fn test_get_valid_lineage_non_species_hit_returns_non_specific_lineage() {
        let mut hits = AHashMap::new();
        hits.insert("read1".to_string(), (42, 2));

        let mut lineage_map = AHashMap::new();
        lineage_map.insert(42, [11, 22, 33]);

        let lineage_map = Arc::new(lineage_map);
        let got = get_valid_lineage(&hits, &lineage_map, "read1");
        assert_eq!(got, [-100, -200, -300]);
    }

    #[test]
    fn test_get_valid_lineage_missing_hit_returns_null_lineage() {
        let hits: AHashMap<String, (Taxid, u8)> = AHashMap::new();
        let lineage_map = Arc::new(AHashMap::new());

        let got = get_valid_lineage(&hits, &lineage_map, "missing");
        assert_eq!(got, [-1, -1, -1]);
    }

    #[tokio::test]
    async fn test_read_file_into_set_ignores_invalid_lines_and_dedups() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("taxids.txt");
        fs::write(&path, b"1\n2\nabc\n2\n  3  \n")?;

        let set = read_file_into_set(&path).await?;
        assert_eq!(set.len(), 3);
        assert!(set.contains(&1));
        assert!(set.contains(&2));
        assert!(set.contains(&3));
        Ok(())
    }

    #[tokio::test]
    async fn test_build_should_keep_filter_applies_white_black_lists() -> Result<()> {
        let dir = tempdir()?;

        let deuterostome = dir.path().join("deuterostome.txt");
        let whitelist = dir.path().join("whitelist.txt");
        let blacklist = dir.path().join("blacklist.txt");

        fs::write(&deuterostome, b"111\n")?;
        fs::write(&whitelist, b"333\n444\n")?;
        fs::write(&blacklist, b"222\n")?;

        let keep = build_should_keep_filter(
            Some(deuterostome.clone()),
            Some(whitelist.clone()),
            Some(blacklist.clone()),
        )
            .await?;

        assert!(keep(&[333]));
        assert!(keep(&[444, 999]));
        assert!(!keep(&[222]));
        assert!(!keep(&[111]));
        assert!(!keep(&[9605]));
        assert!(!keep(&[9606]));
        assert!(!keep(&[333, 222]));
        Ok(())
    }

    #[test]
    fn test_combine_taxon_loc_json_concatenates_inputs() -> Result<()> {
        let dir = tempdir()?;
        let in1 = dir.path().join("one.json");
        let in2 = dir.path().join("two.json");
        let out = dir.path().join("combined.json");

        let a = vec![
            taxon_loc(10, 0, 9, "NT"),
            taxon_loc(20, 10, 19, "NT"),
        ];
        let b = vec![taxon_loc(30, 0, 7, "NR")];

        fs::write(&in1, serde_json::to_vec(&a)?)?;
        fs::write(&in2, serde_json::to_vec(&b)?)?;

        combine_taxon_loc_json(vec![in1, in2], out.clone())?;

        let combined: Vec<TaxonSeqLocation> = serde_json::from_slice(&fs::read(out)?)?;
        assert_eq!(combined.len(), 3);
        assert_eq!(combined[0].taxid, 10);
        assert_eq!(combined[1].taxid, 20);
        assert_eq!(combined[2].taxid, 30);
        Ok(())
    }

    #[test]
    fn test_generate_locator_work_sorts_and_writes_locations() -> Result<()> {
        let dir = tempdir()?;
        let fasta_out = dir.path().join("locator.fasta");
        let json_out = dir.path().join("locator.json");

        let records = Arc::new(vec![
            ("contig_b:species_nt:20:extra".to_string(), "TT".to_string()),
            ("contig_a:species_nt:10:extra".to_string(), "AAAA".to_string()),
        ]);

        generate_locator_work(
            records,
            "species_nt".to_string(),
            "NT".to_string(),
            fasta_out.clone(),
            json_out.clone(),
        )?;

        let fasta = fs::read_to_string(&fasta_out)?;
        assert!(fasta.starts_with(">contig_a:species_nt:10:extra\nAAAA\n"));
        assert!(fasta.contains(">contig_b:species_nt:20:extra\nTT\n"));

        let locs: Vec<TaxonSeqLocation> = serde_json::from_slice(&fs::read(&json_out)?)?;
        assert_eq!(locs.len(), 2);
        assert_eq!(locs[0].taxid, 10);
        assert_eq!(locs[0].hit_type, "NT");
        assert_eq!(locs[1].taxid, 20);
        assert_eq!(locs[1].hit_type, "NT");
        Ok(())
    }
}

