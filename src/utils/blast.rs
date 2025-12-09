// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::{anyhow, Result};
use log::{self, LevelFilter, debug, info, error, warn};
use tokio::io::AsyncBufReadExt;
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;
use fst::Map;
use ahash::AHashMap;
use serde::{Deserialize, Serialize};
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::mpsc::Sender;
use crate::config::defs::{Taxid, Lineage, NT_TAG};
use crate::utils::streams::ParseOutput;
use crate::utils::taxonomy::validate_taxid_lineage;
use crate::utils::streams::ToBytes;


/// Single BLAST m8 line
#[derive(Debug, Clone, Default)]
pub struct M8Record {
    pub qname: String,       // Read ID
    pub tname: String,       // Base Accession ID (no version suffix)
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
    pub qlen: u64,   // Query length: From column in NT; 0 or contig len in NR (unused in NR logic)
    pub slen: u64,   // Subject length: From column in NT; 0 in NR (unused entirely)
}

impl M8Record {
    /// Parse 14-column NT (blastn) output — matches BlastnOutput6NTReader
    pub fn parse_line_nt(line: &str) -> Result<Self> {
        let line = line.trim_end();
        if line.is_empty() {
            return Err(anyhow!("empty line"));
        }

        let mut fields = line.split('\t');

        macro_rules! next {
            () => {
                fields.next().ok_or_else(|| anyhow!("missing field in NT m8 line"))?
            };
        }

        macro_rules! parse_float {
            () => {{
                let s = next!();
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!("invalid float '{}': {}", s, e))?
            }};
        }

        macro_rules! parse_u64 {
            () => {
                next!().parse::<u64>().map_err(|e| anyhow!("invalid u64: {}", e))?
            };
        }

        let qname = next!().to_string();
        let raw_accession = next!(); // e.g., "NC_123456.2"
        let tname = raw_accession
            .split('.')
            .next()
            .unwrap_or(raw_accession)
            .to_string(); // → "NC_123456"

        let pident = parse_float!();
        let alen = parse_u64!();
        let mismatch = parse_u64!();
        let gapopen = parse_u64!();
        let qstart = parse_u64!();
        let qend = parse_u64!();
        let tstart = parse_u64!();
        let tend = parse_u64!();
        let evalue = parse_float!();
        let bitscore = parse_float!();
        let qlen = parse_u64!();
        let slen = parse_u64!();

        // Defensive: Warn on extra columns
        if fields.next().is_some() {
            warn!("Extra columns in NT m8 line (expected 14): {}", line);
        }

        Ok(Self {
            qname,
            tname,
            pident,
            alen,
            mismatch,
            gapopen,
            qstart,
            qend,
            tstart,
            tend,
            evalue,
            bitscore,
            qlen,
            slen,
        })
    }

    /// Parse 12-column NR (blastx) output — matches BlastnOutput6Reader
    pub fn parse_line_nr(line: &str) -> Result<Self> {
        let line = line.trim_end();
        if line.is_empty() {
            return Err(anyhow!("empty line"));
        }

        let mut fields = line.split('\t');

        macro_rules! next {
            () => {
                fields.next().ok_or_else(|| anyhow!("missing field in NR m8 line"))?
            };
        }

        macro_rules! parse_float {
            () => {{
                let s = next!();
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!("invalid float '{}': {}", s, e))?
            }};
        }

        macro_rules! parse_u64 {
            () => {
                next!().parse::<u64>().map_err(|e| anyhow!("invalid u64: {}", e))?
            };
        }

        let qname = next!().to_string();
        let raw_accession = next!(); // e.g., "QIK02963.1"
        let tname = raw_accession
            .split('.')
            .next()
            .unwrap_or(raw_accession)
            .to_string(); // → "QIK02963"

        let pident = parse_float!();
        let alen = parse_u64!();
        let mismatch = parse_u64!();
        let gapopen = parse_u64!();
        let qstart = parse_u64!();
        let qend = parse_u64!();
        let tstart = parse_u64!();
        let tend = parse_u64!();
        let evalue = parse_float!();
        let bitscore = parse_float!();

        // Defensive: Warn on extra columns
        if fields.next().is_some() {
            warn!("Extra columns in NR m8 line (expected 12): {}", line);
        }

        Ok(Self {
            qname,
            tname,
            pident,
            alen,
            mismatch,
            gapopen,
            qstart,
            qend,
            tstart,
            tend,
            evalue,
            bitscore,
            qlen: 0,  // Safe default: Not present/used in NR logic
            slen: 0,  // Safe default: Not present/used anywhere
        })
    }


    pub fn to_tab_string(&self) -> String {
        format!(
            "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{}",
            self.qname,
            self.tname,
            self.pident,
            self.alen,
            self.mismatch,
            self.gapopen,
            self.qstart,
            self.qend,
            self.tstart,
            self.tend,
            self.evalue,
            self.bitscore
        )
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

#[derive(Debug, Clone, Serialize, Deserialize)]
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
                .ok()?;

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

fn parse_read_taxid_from_hitsummary(line: &str) -> Option<(String, i32)> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() >= 3 {
        Some((parts[0].to_string(), parts[2].parse().ok()?))
    } else {
        None
    }
}

pub async fn generate_taxon_count_json_from_m8(
    mut refined_m8_stream: ReceiverStream<ParseOutput>,
    mut refined_hit_summary_stream: ReceiverStream<ParseOutput>,
    db_type: &str,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    duplicate_cluster_sizes: &HashMap<String, u64>,
    mut output_tx: Sender<ParseOutput>,
) -> Result<()> {
    let mut buckets: AHashMap<Taxid, AggBucket> = AHashMap::with_capacity(500_000);

    // 1. Read-to-taxid map from refined hit summary
    let mut read_to_taxid: AHashMap<String, Taxid> = AHashMap::with_capacity(10_000_000);
    while let Some(item) = refined_hit_summary_stream.next().await {
        let bytes = item.to_bytes()?;  // ← note the ? here
        let line = String::from_utf8_lossy(&bytes);

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }

        let read_id = fields[0].to_string();
        let taxid = fields[1].parse::<Taxid>().unwrap_or(0);
        read_to_taxid.insert(read_id, taxid);
    }


    while let Some(item) = refined_m8_stream.next().await {
        let bytes = item.to_bytes()?;  // ← again, ? to propagate any error
        let line = String::from_utf8_lossy(&bytes);

        let m8 = if db_type == NT_TAG {
            M8Record::parse_line_nt(&line)?
        } else {
            M8Record::parse_line_nr(&line)?
        };

        if let Some(&taxid) = read_to_taxid.get(&m8.qname) {
            if taxid <= 0 { continue; }

            let lineage = match lineage_map.get(&taxid) {
                Some(l) => l,
                None => continue,
            };

            if !(should_keep_filter)(lineage) { continue; }

            let bucket = buckets.entry(taxid).or_default();
            bucket.nonunique_count += 1;
            bucket.unique_count += duplicate_cluster_sizes.get(&m8.qname).cloned().unwrap_or(1);
            bucket.base_count += 1;
            bucket.sum_percent_identity += m8.pident;
            bucket.sum_alignment_length += m8.alen as f64;
            bucket.sum_e_value += m8.evalue;
            bucket.source_count_type.insert(db_type.to_string());
        }
    }


    for (taxid, bucket) in buckets {
        let lineage = match lineage_map.get(&taxid) {
            Some(l) => l,
            None => continue,
        };

        if !(should_keep_filter)(lineage) { continue; }

        let dcr = if bucket.nonunique_count > 0 {
            bucket.unique_count as f64 / bucket.nonunique_count as f64
        } else { 0.0 };

        let percent_identity = if bucket.base_count > 0 {
            bucket.sum_percent_identity / bucket.base_count as f64
        } else { 0.0 };

        let alignment_length = if bucket.base_count > 0 {
            bucket.sum_alignment_length / bucket.base_count as f64
        } else { 0.0 };

        let e_value = if bucket.base_count > 0 {
            bucket.sum_e_value / bucket.base_count as f64
        } else { 0.0 };

        let count = TaxonCount {
            tax_id: taxid,
            tax_level: 1, // always species-level aggregation in refined counts
            genus_taxid: lineage[1],
            family_taxid: lineage[2],
            count: bucket.unique_count,
            nonunique_count: bucket.nonunique_count,
            unique_count: bucket.unique_count,
            dcr,
            percent_identity,
            alignment_length,
            e_value,
            count_type: db_type.to_string(),
            base_count: bucket.base_count,
            source_count_type: Some(bucket.source_count_type.iter().cloned().collect()),
        };

        let json_line = serde_json::to_string(&count)? + "\n";
        output_tx
            .send(ParseOutput::Bytes(Arc::new(json_line.into_bytes())))
            .await
            .map_err(|_| anyhow!("taxon count receiver dropped"))?;
    }

    Ok(())
}