use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use std::fs::File as StdFile;
use std::io::{BufReader as StdBufReader};
use std::sync::Arc;
use std::io::Write;

use anyhow::{Result, anyhow};
use log::{info, warn, debug};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use flate2::read::GzDecoder;
use csv::ReaderBuilder;
use sled::{Db, Tree, IVec};
use fst::{MapBuilder};
use needletail::{parse_fastx_file, FastxReader};
use rayon::prelude::*;
use bincode::{encode_into_std_write, decode_from_std_read};
use serde::{Serialize, Deserialize};
use ahash::AHashMap;

use crate::config::defs::{Taxid, Lineage, PipelineError, INVALID_CALL_BASE_ID};

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

/// Builds a accession2taxid sled DB from a known accession2taxid file
/// Skips accessions not present in NT/NR
/// Typical location, sourced from NCBI:
/// ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
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
    // Extract accession bases (no version) from FASTA
    let mut accession_bases: HashSet<String> = HashSet::new();

    for (path, name) in [(nt_file, "NT"), (nr_file, "NR")]
        .into_iter()
        .filter_map(|(opt_path, name)| opt_path.map(|p| (p, name)))
    {
        info!("Extracting accession bases from {}: {}", name, path.display());
        let file = tokio::fs::File::open(path).await?;
        let mut reader = tokio::io::BufReader::new(file).lines();
        let mut count = 0;

        while let Some(line) = reader.next_line().await? {
            if line.starts_with('>') {
                if let Some(acc_base) = line[1..].split_whitespace().next()
                    .and_then(|s| s.split('.').next())
                {
                    accession_bases.insert(acc_base.to_string());
                    count += 1;
                    if count % 1_000_000 == 0 {
                        info!("\tExtracted {}M accession bases from {}", count / 1_000_000, name);
                    }
                }
            }
        }
        info!("Extracted {} accession bases from {}", count, name);
    }


    let db = sled::open(db_path)?;
    let tree: Tree = db.open_tree("acc2taxid")?;

    let mut total_mapped = 0;
    let mut total_skipped = 0;

    //  Process each mapping file
    for mapping_file in mapping_files {
        let filename = mapping_file.file_name().and_then(|s| s.to_str()).unwrap_or("unknown");
        info!("Processing mapping file: {}", filename);

        let file = StdFile::open(mapping_file)?;
        let gz = GzDecoder::new(StdBufReader::new(file));
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(gz);

        let mut mapped = 0;
        let mut skipped = 0;

        for result in rdr.records() {
            let record = match result {
                Ok(r) => r,
                Err(e) => {
                    warn!("Parse error in {}: {}", filename, e);
                    skipped += 1;
                    continue;
                }
            };

            // andle prot.accession2taxid.FULL (rare): only 2 fields
            if record.len() == 2 {
                let acc_base = record[0].to_string();
                if !accession_bases.contains(&acc_base) {
                    skipped += 1;
                    continue;
                }
                let taxid: Taxid = match record[1].parse() {
                    Ok(t) => t,
                    Err(_) => { skipped += 1; continue; }
                };
                // Store full accession (no version in this file)
                tree.insert(acc_base.as_bytes(), taxid.to_le_bytes().as_slice())?;
                mapped += 1;
            }
            // Standard case: 3+ fields
            else if record.len() >= 3 {
                let full_acc = record[0].to_string(); // e.g., YP_009725295.1
                let acc_base = full_acc.split('.').next().unwrap_or(&full_acc).to_string();

                if !accession_bases.contains(&acc_base) {
                    skipped += 1;
                    continue;
                }

                let taxid: Taxid = match record[2].parse() {
                    Ok(t) => t,
                    Err(_) => { skipped += 1; continue; }
                };

                // Store FULL accession with version
                tree.insert(full_acc.as_bytes(), taxid.to_le_bytes().as_slice())?;
                mapped += 1;
            } else {
                skipped += 1;
            }

            if (mapped + skipped) % 1_000_000 == 0 {
                info!("\t{}M records processed (mapped: {}, skipped: {})",
                      (mapped + skipped) / 1_000_000, mapped, skipped);
            }
        }

        info!("Finished {}: mapped={}, skipped={}", filename, mapped, skipped);
        total_mapped += mapped;
        total_skipped += skipped;
    }

    db.flush_async().await?;
    info!("Built acc2taxid DB with {} entries. Skipped {} not in NT/NR.",
          total_mapped, total_skipped);

    Ok(())
}

pub async fn build_fst_acc2taxid(
    mapping_files: &[PathBuf],
    nt_file: Option<&PathBuf>,
    nr_file: Option<&PathBuf>,
    db_path: &PathBuf,
) -> Result<()> {

    let mut accession_bases: HashSet<String> = HashSet::new();

    if let Some(path) = nt_file {
        scan_fasta_bases(path, &mut accession_bases, "NT");
    }
    if let Some(path) = nr_file {
        scan_fasta_bases(path, &mut accession_bases, "NR");
    }

    // fst map (accession → taxid)
    let mut builder = MapBuilder::memory();
    let mut total_mapped = 0;

    for mapping_file in mapping_files {
        let filename = mapping_file.file_name().and_then(|s| s.to_str()).unwrap_or("unknown");
        info!("Processing mapping file: {}", filename);

        let file = StdFile::open(mapping_file)?;
        let gz = GzDecoder::new(StdBufReader::new(file));
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(gz);

        let mut mapped = 0;

        for result in rdr.records() {
            let record = match result {
                Ok(r) => r,
                Err(e) => {
                    warn!("Parse error in {}: {}", filename, e);
                    continue;
                }
            };

            let full_acc = record[0].to_string();
            let acc_base = full_acc.split('.').next().unwrap_or(&full_acc).to_string();

            if !accession_bases.contains(&acc_base) {
                continue;
            }

            let taxid_idx = if record.len() == 2 { 1 } else { 2 };
            let taxid: Taxid = match record[taxid_idx].parse() {
                Ok(t) => t,
                Err(_) => continue,
            };

            builder.insert(full_acc.as_bytes(), taxid as u64)?;
            mapped += 1;
            total_mapped += 1;

            if total_mapped % 1_000_000 == 0 {
                info!("\t{}M accessions mapped", total_mapped / 1_000_000);
            }
        }

        info!("Finished {}: mapped={}", filename, mapped);
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


// *******************
// DB access functions
// *******************


pub fn validate_taxid_lineage(lineage: &[i32], hit_taxid: Taxid, hit_level: u8) -> Vec<i32> {
    let mut cleaned = lineage.to_vec();
    for level in 0..(hit_level as usize - 1) {
        cleaned[level] = INVALID_CALL_BASE_ID - (level as i32 + 1) * 100;
    }
    let mut parent = hit_taxid;
    for level in (hit_level as usize - 1)..cleaned.len() {
        if cleaned[level] <= 0 {
            cleaned[level] = INVALID_CALL_BASE_ID - (level as i32 + 1) * 100 - parent;
        }
        parent = cleaned[level];
    }
    cleaned
}