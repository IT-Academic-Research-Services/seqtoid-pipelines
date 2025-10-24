use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use std::fs::File as StdFile;
use std::io::{BufReader as StdBufReader};

use anyhow::{Result, anyhow};
use log::{info, warn, debug};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use csv::ReaderBuilder;
use sled::{Db, Tree, IVec};

use crate::config::defs::{Taxid, Lineage};

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

    // load parent map taxid → parent_taxid and rank map taxid → rank
    let mut parent_map: HashMap<Taxid, Taxid> = HashMap::new();
    let mut rank_map: HashMap<Taxid, String> = HashMap::new();
    let nodes_file = File::open(nodes_path).await?;
    let mut nodes_reader = BufReader::new(nodes_file).lines();
    while let Some(line) = nodes_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 3 { continue; }
        let taxid: Taxid = parts[0].parse()?;
        let parent: Taxid = parts[1].parse()?;
        let rank = parts[2].to_string();
        if taxid != 1 {  // Skip root
            parent_map.insert(taxid, parent);
            rank_map.insert(taxid, rank);
        }
    }
    info!("Loaded {} parents and ranks from nodes.dmp", parent_map.len());

    // load merged: old_taxid → new_taxid
    let mut merged_map: HashMap<Taxid, Taxid> = HashMap::new();
    let merged_file = File::open(merged_path).await?;
    let mut merged_reader = BufReader::new(merged_file).lines();
    while let Some(line) = merged_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 2 { continue; }
        let old: Taxid = parts[0].parse()?;
        let new: Taxid = parts[1].parse()?;
        merged_map.insert(old, new);
    }
    info!("Loaded {} merges from merged.dmp", merged_map.len());

    let db = sled::open(db_path)?;
    let tree: Tree = db.open_tree("lineages")?;

    // build 3-level lineages for all taxids (species, genus, family; leaf-to-root)
    let mut seen = HashSet::new();
    for (&taxid, _) in &parent_map {
        if seen.contains(&taxid) { continue; }
        let mut species: Taxid = -100;
        let mut genus: Taxid = -200;
        let mut family: Taxid = -300;
        let mut current = taxid;
        let mut depth = 0;
        while current != 1 && depth < 20 {  // Limit depth to prevent infinite loops
            let rank = rank_map.get(&current).cloned().unwrap_or_default();
            if rank == "species" && species == -100 {
                species = current;
            } else if rank == "genus" && genus == -200 {
                genus = current;
            } else if rank == "family" && family == -300 {
                family = current;
            }
            if species > -100 && genus > -200 && family > -300 {
                break;
            }
            current = *parent_map.get(&current).unwrap_or(&1);
            if let Some(merged) = merged_map.get(&current) {
                current = *merged;  // redirect deprecated using the merge map
            }
            depth += 1;
        }

        // Pack to bytes: three i32 taxids (little-endian), leaf-to-root: species, genus, family
        let mut bytes = Vec::with_capacity(12);
        bytes.extend_from_slice(&species.to_le_bytes());
        bytes.extend_from_slice(&genus.to_le_bytes());
        bytes.extend_from_slice(&family.to_le_bytes());

        tree.insert(taxid.to_le_bytes(), bytes)?;
        seen.insert(taxid);
    }

    db.flush_async().await?;
    info!("Built lineages DB with {} entries at {}", seen.len(), db_path.display());
    Ok(())
}



/// Builds a accession2taxid sled DB from a known accession2taxid file
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
pub async fn build_accession2taxid_db(
    gz_path: &PathBuf,
    db_path: &PathBuf,
) -> Result<()> {
    let db = sled::open(db_path)?;
    let tree: Tree = db.open_tree("acc2taxid")?;

    let file = StdFile::open(gz_path)?;
    let gz = flate2::read::GzDecoder::new(StdBufReader::new(file));
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(gz);

    let mut count = 0;
    for result in rdr.records() {
        let record = result?;
        if record.len() < 3 {
            warn!("Skipping invalid record: {:?}", record);
            continue;
        }
        let acc = record[0].as_bytes().to_vec(); // Use accession (first column)
        let taxid: Taxid = match record[2].parse() { // Use taxid (third column)
            Ok(t) => t,
            Err(e) => {
                warn!("Failed to parse taxid from record: {:?}", record);
                return Err(anyhow!("Failed to parse taxid: {}", e));
            }
        };
        tree.insert(acc, taxid.to_le_bytes().as_slice())?;
        count += 1;
    }

    db.flush_async().await?;
    info!("Built acc2taxid DB with {} entries at {}", count, db_path.display());
    Ok(())
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