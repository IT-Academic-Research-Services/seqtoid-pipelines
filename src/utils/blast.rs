// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, Result};
use tokio::io::AsyncBufReadExt;
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;

use crate::config::defs::{Taxid, Lineage};


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
            () => { fields.next().ok_or_else(|| anyhow!("missing field"))? };
        }
        macro_rules! parse_float {
            () => {{
                let s = next!();
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!("invalid float '{}': {}", s, e))?
            }};
        }

        Ok(Self {
            qname: next!().to_string(),
            tname: next!().to_string(),
            pident: parse_float!(),
            alen: next!().parse()?,
            mismatch: next!().parse()?,
            gapopen: next!().parse()?,
            qstart: next!().parse()?,
            qend: next!().parse()?,
            tstart: next!().parse()?,
            tend: next!().parse()?,
            evalue: parse_float!(),
            bitscore: parse_float!(),
        })
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

#[derive(Debug, Clone)]
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


/// Computes the common lineage for a set of hits, matching m8.py _common_lineage logic.
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
    lineage_map: &HashMap<Taxid, Lineage>,
) -> (u8, i64, Vec<M8Record>) {
    let lineages: Vec<Lineage> = hits
        .iter()
        // Map accession (tname) → taxid → lineage [species, genus, family]
        .filter_map(|r| {
            let taxid: Taxid = r.tname.parse().ok()?;
            lineage_map.get(&taxid).cloned()
        })
        .collect();

    if lineages.is_empty() {
        return (0, 0, Vec::new());
    }

    let mut max_level = 0u8;
    let mut consensus_taxid = 0i64;

    for level in 0..3 {  // 0=species, 1=genus, 2=family
        let taxids: HashSet<i32> = lineages
            .iter()
            .map(|lin| lin.get(level).copied().unwrap_or(-1))
            .filter(|&t| t > 0)
            .collect();

        if taxids.len() == 1 {
            let agreed = *taxids.iter().next().unwrap() as i64;
            max_level = (level + 1) as u8;  // Update to this more specific level
            consensus_taxid = agreed;
        } else {
            // No agreement at this level; stop as deeper levels won't agree
            break;
        }
    }

    (max_level, consensus_taxid, hits.to_vec())
}