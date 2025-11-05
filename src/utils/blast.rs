// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::{anyhow, Result};
use tokio::io::AsyncBufReadExt;
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;
use fst::Map;
use ahash::AHashMap;

use crate::config::defs::{Taxid, Lineage};
use crate::utils::taxonomy::validate_taxid_lineage;


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


/// Computes the common lineage for a set of hits,
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
pub fn consensus_level<D: AsRef<[u8]>>(
    hits: &[M8Record],
    lineage_map: &AHashMap<Taxid, Lineage>,
    acc2taxid_map: &Map<D>,
    should_keep: &Arc<impl Fn(&[i32]) -> bool + Send + Sync>,
) -> Result<(u8, i64, Vec<M8Record>)>   {
    let lineages: Vec<Lineage> = hits
        .iter()
        .filter_map(|r| {
            let accession = &r.tname;

            // 1. Try full accession (e.g., "ACC.1")
            let taxid_opt = acc2taxid_map.get(accession.as_bytes());

            // 2. Fallback: strip version → base accession (e.g., "ACC")
            let taxid_u64 = taxid_opt.or_else(|| {
                let base_acc = accession.split('.').next().unwrap_or(accession);
                if base_acc != accession {
                    acc2taxid_map.get(base_acc.as_bytes())
                } else {
                    None
                }
            })
                .ok_or_else(|| anyhow!("Accession not found in acc2taxid map: {}", accession))
                .ok()?; // ← This `?` is allowed here because we're in a `filter_map` returning `Option`

            let taxid = taxid_u64 as i32;
            if taxid <= 0 {
                return None;
            }

            lineage_map.get(&taxid).cloned()
        })
        .collect();

    if lineages.is_empty() {
        return Ok((0, 0, Vec::new()));
    }

    let mut max_level = 0u8;
    let mut consensus_taxid = 0i64;

    for level in 0..3 {
        let taxids: HashSet<i32> = lineages
            .iter()
            .map(|lin| lin[level])
            .filter(|&t| t > 0)
            .collect();

        if taxids.len() == 1 {
            let agreed = *taxids.iter().next().unwrap() as i64;
            max_level = (level + 1) as u8;
            consensus_taxid = agreed;
        } else {
            break;
        }
    }

    if max_level > 0 {
        let rep_lineage = lineages[0];
        let validated = validate_taxid_lineage(&rep_lineage, consensus_taxid as i32, max_level);
        if !(should_keep)(&validated) {
            return Ok((0, 0, Vec::new()));
        }
    }

    Ok((max_level, consensus_taxid, hits.to_vec()))
}