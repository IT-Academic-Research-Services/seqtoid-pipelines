// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, Result};
use tokio::io::AsyncBufReadExt;
use tokio_stream::StreamExt;

use crate::config::defs::{Taxid, Lineage};

const INVALID_CALL_BASE_ID: i64 = -100;

/// Single BLAST m8 line
#[derive(Debug, Clone)]
pub struct M8Record {
    pub qname: String,       // Read ID
    pub tname: String,       // Accession ID
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
}

impl M8Record {
    pub fn parse_line(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        macro_rules! next {
            () => {
                fields.next().ok_or_else(|| anyhow!("missing field"))?
            };
        }
        Ok(Self {
            qname: next!().to_string(),
            tname: next!().to_string(),
            pident: next!().parse()?,
            alen: next!().parse()?,
            mismatch: next!().parse()?,
            gapopen: next!().parse()?,
            qstart: next!().parse()?,
            qend: next!().parse()?,
            tstart: next!().parse()?,
            tend: next!().parse()?,
            evalue: next!().parse()?,
            bitscore: next!().parse()?,
        })
    }
}


// *******************
// PTaxonomy helper functions based on m8 records
// *******************

/// Build a negative taxid for “no-specific-call” at a given level.
fn negative_taxid(level: u8, parent: i64) -> i64 {
    INVALID_CALL_BASE_ID - (level as i64) * 100 - parent
}


/// Walk up the lineage until we find a level where all hits agree.
///
/// # Arguments
///
///  * `hits` - list of m8 records
/// * `lineage_map` - derived from taxid DB
///
/// # Returns
///
/// Tuple ()
pub fn consensus_level(
    hits: &[M8Record],
    lineage_map: &HashMap<Taxid, Lineage>,
) -> (u8, i64, Vec<M8Record>) {
    let lineages: Vec<&Lineage> = hits
        .iter()
        // Map accession (tname) → taxid → lineage
        .filter_map(|r| {
            let taxid: Taxid = r.tname.parse().ok()?;
            lineage_map.get(&taxid)
        })
        .collect();

    if lineages.is_empty() {
        return (0, 0, Vec::new());
    }

    let max_depth = lineages.iter().map(|v| v.len()).max().unwrap_or(0);
    for depth in 0..max_depth {
        let taxids: HashSet<i64> = lineages
            .iter()
            .filter_map(|lin| lin.get(depth).copied().map(|t| t as i64))
            .collect();

        if taxids.len() == 1 {
            let consensus = *taxids.iter().next().unwrap();
            return (depth as u8 + 1, consensus, hits.to_vec());
        }
    }
    (0, 0, Vec::new())
}