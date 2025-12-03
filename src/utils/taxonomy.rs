use tokio_stream::StreamExt;
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
use tokio_stream::wrappers::ReceiverStream;
use crate::config::defs::{Taxid, Lineage, PipelineError, INVALID_CALL_BASE_ID};
use crate::utils::blast::{M8Record, AggBucket, TaxonCount};
use crate::utils::streams::ParseOutput;
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
    if let Some(p) = nt_file {
        scan_fasta_bases(p, &mut accession_bases, "NT");
    }
    if let Some(p) = nr_file {
        scan_fasta_bases(p, &mut accession_bases, "NR");
    }

    let mut entries = Vec::new();

    for mapping_file in mapping_files {
        let filename = mapping_file.file_name().and_then(|s| s.to_str()).unwrap_or("unknown");
        info!("Reading mapping file: {}", filename);

        let file = StdFile::open(mapping_file)?;
        let gz = GzDecoder::new(StdBufReader::new(file));
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(gz);

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

            entries.push(AccTaxEntry { acc: full_acc, taxid });
        }
    }


    entries.sort_unstable_by(|a, b| a.acc.cmp(&b.acc));

    let mut builder = MapBuilder::memory();
    let mut total_mapped = 0;
    let mut last_acc: Option<String> = None;

    for entry in entries {
        // Skip duplicate keys (same accession)
        if last_acc.as_deref() == Some(&entry.acc) {
            continue;
        }

        last_acc = Some(entry.acc.clone());

        builder.insert(entry.acc.as_bytes(), entry.taxid as u64)?;
        total_mapped += 1;

        if total_mapped % 1_000_000 == 0 {
            info!("\t{}M accessions mapped", total_mapped / 1_000_000);
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
///
/// # Arguments
/// * `path` - path to hit sumamry file.
/// # Returns
/// hash map of HitSummaries
pub async fn get_top_m8_nt(
    mut input: ReceiverStream<ParseOutput>,
    mut output_tx: tokio::sync::mpsc::Sender<ParseOutput>,
) -> Result<()> {
    let mut best_per_contig: HashMap<String, M8Record> = HashMap::new();

    while let Some(item) = input.next().await {
        let bytes = match item {
            ParseOutput::Bytes(b) => b,
            _ => continue,
        };
        let line = String::from_utf8_lossy(&bytes);
        let line = line.trim_end();
        if line.is_empty() { continue; }

        let m8 = match M8Record::parse_line(line) {
            Ok(m8) => m8,
            Err(e) => {
                warn!("Failed to parse m8 line: {} — {}", e, line);
                continue;
            }
        };

        let contig_id = m8.qname.clone();
        if let Some(current) = best_per_contig.get_mut(&contig_id) {
            if m8.bitscore > current.bitscore ||
                (m8.bitscore == current.bitscore && m8.evalue < current.evalue) {
                *current = m8;
            }
        } else {
            best_per_contig.insert(contig_id, m8);
        }
    }

    for (_, best) in best_per_contig {
        let line = best.to_tab_string() + "\n";
        output_tx.send(ParseOutput::Bytes(Arc::new(line.into_bytes()))).await?;
    }

    Ok(())
}

pub async fn get_top_m8_nr(
    mut input: ReceiverStream<ParseOutput>,
    mut output_tx: tokio::sync::mpsc::Sender<ParseOutput>,
) -> Result<()> {
    // Identical to NT for now (adjust if NR ranking differs)
    get_top_m8_nt(input, output_tx).await
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