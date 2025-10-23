use std::path::Path;
use std::collections::{HashMap, HashSet};

use anyhow::{Result, anyhow};
use log::{info, warn, debug};
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use sled::{Db, Tree};

pub type Taxid = u32;  // NCBI taxids fit in u32

// *******************
// DB creation functions
// *******************
pub async fn build_taxid_lineages_db(
    nodes_path: &Path,  // nodes.dmp file
    merged_path: &Path, // merged.dmp
    db_path: &Path,     // output .db
) -> Result<Db> {

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
    info!("Built lineages DB with {} entries", seen.len());
    Ok(db)
}