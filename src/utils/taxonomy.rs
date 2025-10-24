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

    // load parent map taxid → parent_taxid
    let mut parent_map: HashMap<Taxid, Taxid> = HashMap::new();
    let nodes_file = File::open(nodes_path).await?;
    let mut nodes_reader = BufReader::new(nodes_file).lines();
    while let Some(line) = nodes_reader.next_line().await? {
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() < 2 { continue; }
        let taxid: Taxid = parts[0].parse()?;
        let parent: Taxid = parts[1].parse()?;
        if taxid != 1 {  // Skip root
            parent_map.insert(taxid, parent);
        }
    }
    info!("Loaded {} parents from nodes.dmp", parent_map.len());

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

    // build lineages for all taxids (traverse to root)
    let mut seen = HashSet::new();
    for (&taxid, _) in &parent_map {
        if seen.contains(&taxid) { continue; }
        let mut lineage = Vec::new();
        let mut current = taxid;
        while current != 1 {  // Stop at root
            lineage.push(current);
            current = *parent_map.get(&current).unwrap_or(&1);  // Fallback to root
            if let Some(merged) = merged_map.get(&current) {
                current = *merged;  // redirect deprecated using the merge map
            }
        }
        lineage.reverse();  // root going leafward

        // Pack to bytes: u32 count + u32 taxids (little-endian)
        let mut bytes = Vec::with_capacity(4 + 4 * lineage.len());
        bytes.extend_from_slice(&(lineage.len() as u32).to_le_bytes());
        for &t in &lineage {
            bytes.extend_from_slice(&t.to_le_bytes());
        }

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
    if data.len() < 4 { return Err(anyhow!("short lineage")); }
    let count = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    if data.len() != 4 + 4 * count { return Err(anyhow!("malformed lineage")); }
    let mut lineage = Vec::with_capacity(count);
    for i in 0..count {
        let start = 4 + 4 * i;
        lineage.push(u32::from_le_bytes(data[start..start+4].try_into().unwrap()));
    }
    Ok(lineage)
}