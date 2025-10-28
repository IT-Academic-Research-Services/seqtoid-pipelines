use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use std::fs::File as StdFile;
use std::io::{BufReader as StdBufReader};

use anyhow::{Result, anyhow};
use log::{info, warn, debug};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use flate2::read::GzDecoder;
use csv::ReaderBuilder;
use sled::{Db, Tree, IVec};

use crate::config::defs::{Taxid, Lineage, PipelineError, INVALID_CALL_BASE_ID};

// *******************
// DB creation functions
// *******************

/// Builds a taxid-lineages sled DB from taxdump.tar.gz files
/// Typical location, sourced from NCBI:
/// ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
///
/// # Arguments
///
///  * `nodes_path` - path to the nodes.dmp file
/// * `merged_path` - path to the megred.dmp file
/// * `db_path` -  output db path
///
/// # Returns
///
/// Result
pub async fn build_taxid_lineages_db(
    nodes_path: &PathBuf,
    merged_path: &PathBuf,
    db_path: &PathBuf,
) -> Result<()> {
    // Pre-allocate maps for efficiency (NCBI taxonomy ~2.7M entries)
    let mut parent_map: HashMap<Taxid, Taxid> = HashMap::with_capacity(3_000_000);
    let mut rank_map: HashMap<Taxid, String> = HashMap::with_capacity(3_000_000);

    let nodes_file = File::open(nodes_path).await?;
    let mut nodes_reader = BufReader::new(nodes_file).lines();
    while let Some(line) = nodes_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 3 {
            warn!("Skipping invalid line in nodes.dmp: {}", line);
            continue;
        }
        let taxid: Taxid = match parts[0].parse() {
            Ok(t) => t,
            Err(e) => {
                warn!("Invalid taxid in nodes.dmp: {}, error: {}", parts[0], e);
                continue;
            }
        };
        let parent: Taxid = match parts[1].parse() {
            Ok(p) => p,
            Err(e) => {
                warn!("Invalid parent taxid in nodes.dmp: {}, error: {}", parts[1], e);
                continue;
            }
        };
        let rank = parts[2].to_string();
        if taxid != 1 { // Skip root
            parent_map.insert(taxid, parent);
            rank_map.insert(taxid, rank);
        }
    }
    info!("Loaded {} parents and ranks from nodes.dmp", parent_map.len());

    let mut merged_map: HashMap<Taxid, Taxid> = HashMap::with_capacity(100_000);
    let merged_file = File::open(merged_path).await?;
    let mut merged_reader = BufReader::new(merged_file).lines();
    while let Some(line) = merged_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 2 {
            warn!("Skipping invalid line in merged.dmp: {}", line);
            continue;
        }
        let old: Taxid = match parts[0].parse() {
            Ok(o) => o,
            Err(e) => {
                warn!("Invalid old taxid in merged.dmp: {}, error: {}", parts[0], e);
                continue;
            }
        };
        let new: Taxid = match parts[1].parse() {
            Ok(n) => n,
            Err(e) => {
                warn!("Invalid new taxid in merged.dmp: {}, error: {}", parts[1], e);
                continue;
            }
        };
        merged_map.insert(old, new);
    }
    info!("Loaded {} merges from merged.dmp", merged_map.len());

    let db = sled::open(db_path)?;
    let tree: Tree = db.open_tree("lineages")?;
    let mut seen = HashSet::with_capacity(parent_map.len());
    let mut missing_ranks = 0;

    for (&taxid, _) in &parent_map {
        if seen.contains(&taxid) {
            continue;
        }
        let mut species: Taxid = -100;
        let mut genus: Taxid = -200;
        let mut family: Taxid = -300;
        let mut current = taxid;
        let mut depth = 0;
        let max_depth = 50; // Increased for complex hierarchies

        while current != 1 && depth < max_depth {
            let rank = rank_map.get(&current).cloned().unwrap_or_default();
            if rank == "species" && species == -100 {
                species = current;
            } else if rank == "genus" && genus == -200 {
                genus = current;
            } else if rank == "family" && family == -300 {
                family = current;
            }
            if species != -100 && genus != -200 && family != -300 {
                break; // All ranks found
            }
            current = *parent_map.get(&current).unwrap_or(&1);
            if let Some(merged) = merged_map.get(&current) {
                current = *merged;
            }
            depth += 1;
        }

        if depth >= max_depth {
            warn!("Depth limit reached for taxid {}: species={}, genus={}, family={}", taxid, species, genus, family);
            missing_ranks += 1;
        }

        // Serialize lineage to bytes
        let mut bytes = Vec::with_capacity(12);
        bytes.extend_from_slice(&species.to_le_bytes());
        bytes.extend_from_slice(&genus.to_le_bytes());
        bytes.extend_from_slice(&family.to_le_bytes());

        tree.insert(taxid.to_le_bytes(), bytes)?;
        seen.insert(taxid);
    }

    db.flush_async().await?;
    info!("Built lineages DB with {} entries at {}. {} taxids reached depth limit.", seen.len(), db_path.display(), missing_ranks);
    Ok(())
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
    gz_path: &PathBuf,
    nt_file: Option<&PathBuf>,
    nr_file: Option<&PathBuf>,
    db_path: &PathBuf,
) -> Result<()> {
    // Collect accession IDs from NT/NR FASTA files
    let mut accessions = HashSet::with_capacity(50_000_000); // Typical NT/NR size
    for file_path in [nt_file, nr_file].into_iter().flatten() {
        let file = File::open(file_path).await?;
        let mut reader = BufReader::new(file).lines();
        let mut count = 0;
        while let Some(line) = reader.next_line().await? {
            if line.starts_with('>') {
                let accession = line[1..]
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| anyhow!("Invalid FASTA header: {}", line))?
                    .split('.')
                    .next()
                    .ok_or_else(|| anyhow!("Invalid accession in header: {}", line))?;
                accessions.insert(accession.to_string());
                count += 1;
                if count % 1_000_000 == 0 {
                    info!("Extracted {}M accession names from {}", count / 1_000_000, file_path.display());
                }
            }
        }
        info!("Extracted {} accession names from {}", count, file_path.display());
    }
    info!("Total unique accessions from NT/NR: {}", accessions.len());


    let db = sled::open(db_path)?;
    let tree: Tree = db.open_tree("acc2taxid")?;

    let file = StdFile::open(gz_path)?;
    let gz = GzDecoder::new(StdBufReader::new(file));
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(gz);

    let mut count = 0;
    let mut skipped = 0;
    let mut duplicates = 0;

    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                warn!("Failed to read record: {}", e);
                skipped += 1;
                continue;
            }
        };
        if record.len() < 3 {
            warn!("Skipping invalid record (too few fields): {:?}", record);
            skipped += 1;
            continue;
        }
        let acc = record[0].to_string();
        let acc_no_version = acc.split('.').next().unwrap_or(&acc);
        if !accessions.contains(acc_no_version) {
            skipped += 1;
            continue;
        }
        let taxid: Taxid = match record[2].parse() {
            Ok(t) => t,
            Err(e) => {
                warn!("Skipping record with invalid taxid: {:?}, error: {}", record, e);
                skipped += 1;
                continue;
            }
        };
        if tree.contains_key(acc.as_bytes())? {
            warn!("Duplicate accession found: {}", acc);
            duplicates += 1;
        }
        tree.insert(acc.as_bytes(), taxid.to_le_bytes().as_slice())?;
        count += 1;
        if count % 1_000_000 == 0 {
            info!("Processed {}M records", count / 1_000_000);
        }
    }

    db.flush_async().await?;
    info!(
        "Built acc2taxid DB with {} entries at {}. Skipped {} records ({} invalid, {} not in NT/NR), found {} duplicates.",
        count, db_path.display(), skipped, skipped - count, count, duplicates
    );
    Ok(())
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
) -> Result<impl Fn(&[i32]) -> bool> {
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

/// Unpacks the byte arrangement made in the sled DB build_taxid_lineages_db
/// Typical location, sourced from NCBI:
/// ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
///
/// # Arguments
///
///  * `gz_path` - accession2taxid.gz
/// * `db_path` -  output db path
///
/// # Returns
/// Result
pub fn unpack_lineage_bytes(ivec: &IVec) -> Result<Vec<Taxid>> {
    let data = ivec.as_ref();
    if data.len() != 12 { return Err(anyhow!("malformed lineage: expected 12 bytes")); }
    let mut lineage = Vec::with_capacity(3);
    for i in 0..3 {
        let start = 4 * i;
        lineage.push(i32::from_le_bytes(data[start..start+4].try_into().unwrap()));
    }
    Ok(lineage)
}


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