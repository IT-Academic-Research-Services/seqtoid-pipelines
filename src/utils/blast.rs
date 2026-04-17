// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};
use std::env::args;
use std::sync::Arc;
use std::path::PathBuf;
use std::hash::{BuildHasher, Hash, Hasher};
use rayon::prelude::*;
use memchr::memchr;

use anyhow::{anyhow, Context, Result};
use log::{self, LevelFilter, debug, info, error, warn};
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;
use fst::Map;
use ahash::{AHashMap, RandomState as AHashState};
use serde::{Deserialize, Serialize};
use tokio_stream::wrappers::ReceiverStream;
use tokio::fs::File as TokioFile;
use tokio::sync::mpsc;
use tokio::sync::mpsc::{channel, Sender};
use tokio::io::{AsyncBufReadExt, BufReader, BufWriter, AsyncWriteExt};
use tokio_stream::wrappers::BroadcastStream;
use tokio::sync::Semaphore;
use tokio::task::JoinHandle;
use dashmap::DashMap;
use tokio::time::{sleep, Duration, Instant};
use ahash::RandomState as AHashRandomState;
use bytes::Bytes;


use once_cell::sync::Lazy;

use crate::config::defs::{Taxid, Lineage, NT_TAG, NR_TAG, RunConfig, ClusterInfo, MIN_NORMAL_POSITIVE_DOUBLE, READ_COUNTING_MODE, ReadCountingMode, SIMD_LEVEL, SimdLevel};
use crate::utils::streams::ParseOutput;
use crate::utils::taxonomy::validate_taxid_lineage;
use crate::utils::streams::ToBytes;
use crate::utils::file::write_byte_stream_to_file;
use crate::utils::system::{compute_batch_size, compute_phase_concurrency};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigSummaryEntry {
    pub contig_name: String,
    pub common_name: Option<String>,
    pub category_name: Option<String>,
    pub score: Option<f64>,
    pub db_type: String,
    pub reads: u64,
    pub bases: u64,
    pub species_taxid: i32,
    pub genus_taxid: i32,
    pub family_taxid: i32,
}

#[derive(Debug, Clone, Default)]
pub struct SpeciesAlignmentResults {
    pub contig: Option<Taxid>,
    pub read: Option<Taxid>,
}

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

//for read-grouped bastching
#[derive(Debug, Clone)]
pub struct ReducedRead {
    pub seq: u64,
    pub dedup: Vec<u8>,
    pub summary: Vec<u8>,
    pub accession: String,
}

#[derive(Debug)]
pub struct PendingRead {
    pub read_id: String,
    pub hits: Vec<M8Record>,
}

#[derive(Debug)]
pub enum WorkerMsg {
    Line { seq: u64, line: Vec<u8> },
    Flush,
}

impl M8Record {
    // ── NT scalar ──────────────────────────────────────────────────────────

    /// Parse 14-column NT (blastn) output — scalar baseline.
    fn parse_line_nt_scalar(line: &str) -> Result<Self> {
        let line = line.trim_end();
        if line.is_empty() {
            return Err(anyhow!("empty line"));
        }
        let mut fields = line.split('\t');

        macro_rules! next {
            ($name:literal) => {
                fields.next().ok_or_else(|| anyhow!(concat!("missing ", $name, " in NT m8 line")))?
            };
        }
        macro_rules! parse_float {
            ($name:literal) => {{
                let s = next!($name);
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!(concat!("invalid float ", $name, ": {}"), e))?
            }};
        }
        macro_rules! parse_u64 {
            ($name:literal) => {{
                let s = next!($name);
                lexical_parse::<u64, _>(s.as_bytes())
                    .map_err(|e| anyhow!(concat!("invalid u64 ", $name, ": {}"), e))?
            }};
        }

        let qname           = next!("qname").to_string();
        let raw_accession   = next!("tname");
        let tname           = raw_accession.split('.').next().unwrap_or(raw_accession).to_string();
        let pident          = parse_float!("pident");
        let alen            = parse_u64!("alen");
        let mismatch        = parse_u64!("mismatch");
        let gapopen         = parse_u64!("gapopen");
        let qstart          = parse_u64!("qstart");
        let qend            = parse_u64!("qend");
        let tstart          = parse_u64!("tstart");
        let tend            = parse_u64!("tend");
        let evalue          = parse_float!("evalue");
        let bitscore        = parse_float!("bitscore");
        let qlen            = parse_u64!("qlen");
        let slen            = parse_u64!("slen");

        if fields.next().is_some() {
            warn!("Extra columns in NT m8 line (expected 14): {}", line);
        }
        Ok(Self { qname, tname, pident, alen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore, qlen, slen })
    }

    // ── NR scalar ──────────────────────────────────────────────────────────

    /// Parse 12-column NR (blastx) output — scalar baseline.
    fn parse_line_nr_scalar(line: &str) -> Result<Self> {
        let line = line.trim_end();
        if line.is_empty() {
            return Err(anyhow!("empty line"));
        }
        let mut fields = line.split('\t');

        macro_rules! next {
            ($name:literal) => {
                fields.next().ok_or_else(|| anyhow!(concat!("missing ", $name, " in NR m8 line")))?
            };
        }
        macro_rules! parse_float {
            ($name:literal) => {{
                let s = next!($name);
                lexical_parse::<f64, _>(s.as_bytes())
                    .map_err(|e| anyhow!(concat!("invalid float ", $name, ": {}"), e))?
            }};
        }
        macro_rules! parse_u64 {
            ($name:literal) => {{
                let s = next!($name);
                lexical_parse::<u64, _>(s.as_bytes())
                    .map_err(|e| anyhow!(concat!("invalid u64 ", $name, ": {}"), e))?
            }};
        }

        let qname           = next!("qname").to_string();
        let raw_accession   = next!("tname");
        let tname           = raw_accession.split('.').next().unwrap_or(raw_accession).to_string();
        let pident          = parse_float!("pident");
        let alen            = parse_u64!("alen");
        let mismatch        = parse_u64!("mismatch");
        let gapopen         = parse_u64!("gapopen");
        let qstart          = parse_u64!("qstart");
        let qend            = parse_u64!("qend");
        let tstart          = parse_u64!("tstart");
        let tend            = parse_u64!("tend");
        let evalue          = parse_float!("evalue");
        let bitscore        = parse_float!("bitscore");

        if fields.next().is_some() {
            warn!("Extra columns in NR m8 line (expected 12): {}", line);
        }
        Ok(Self { qname, tname, pident, alen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore, qlen: 0, slen: 0 })
    }

    // ── shared AVX-512 tab scanner ─────────────────────────────────────────

    /// Collect every tab position in `bytes` using a 64-byte AVX-512 sweep.
    /// Falls back to scalar for the tail (< 64 remaining bytes).
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512f,avx512bw")]
    unsafe fn collect_tab_positions(bytes: &[u8]) -> Vec<usize> {
        use std::arch::x86_64::*;
        use std::arch::x86_64::__m512i;
        let len = bytes.len();
        let mut tab_pos: Vec<usize> = Vec::with_capacity(16);
        let tab_splat = _mm512_set1_epi8(b'\t' as i8);
        let mut i = 0usize;
        while i + 64 <= len {
            let chunk = _mm512_loadu_si512(bytes.as_ptr().add(i).cast::<__m512i>());
            let mut mask: u64 = _mm512_cmpeq_epi8_mask(chunk, tab_splat);
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                tab_pos.push(i + bit);
                mask &= mask - 1;
            }
            i += 64;
        }
        while i < len {
            if bytes[i] == b'\t' { tab_pos.push(i); }
            i += 1;
        }
        tab_pos
    }

    // ── NT AVX-512 ─────────────────────────────────────────────────────────

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512f,avx512bw")]
    unsafe fn parse_line_nt_avx512_inner(line: &str) -> Result<Self> {
        let bytes = line.trim_end().as_bytes();
        if bytes.is_empty() { return Err(anyhow!("empty line")); }

        let tab_pos = Self::collect_tab_positions(bytes);
        let len = bytes.len();

        if tab_pos.len() < 13 {
            return Err(anyhow!("NT m8 line has {} tabs, need ≥13 for 14 fields", tab_pos.len()));
        }

        macro_rules! field {
            ($n:expr) => {{
                let start = if $n == 0 { 0 } else { tab_pos[$n - 1] + 1 };
                let end   = if $n < tab_pos.len() { tab_pos[$n] } else { len };
                &bytes[start..end]
            }};
        }
        macro_rules! parse_u64 {
            ($n:expr, $name:literal) => {
                lexical_parse::<u64, _>(field!($n))
                    .map_err(|e| anyhow!(concat!($name, ": {}"), e))?
            };
        }
        macro_rules! parse_f64 {
            ($n:expr, $name:literal) => {
                lexical_parse::<f64, _>(field!($n))
                    .map_err(|e| anyhow!(concat!($name, ": {}"), e))?
            };
        }

        let qname = std::str::from_utf8(field!(0)).map_err(|_| anyhow!("qname not UTF-8"))?.to_string();
        let raw   = std::str::from_utf8(field!(1)).map_err(|_| anyhow!("tname not UTF-8"))?;
        let tname = raw.split('.').next().unwrap_or(raw).to_string();

        let pident   = parse_f64!(2,  "pident");
        let alen     = parse_u64!(3,  "alen");
        let mismatch = parse_u64!(4,  "mismatch");
        let gapopen  = parse_u64!(5,  "gapopen");
        let qstart   = parse_u64!(6,  "qstart");
        let qend     = parse_u64!(7,  "qend");
        let tstart   = parse_u64!(8,  "tstart");
        let tend     = parse_u64!(9,  "tend");
        let evalue   = parse_f64!(10, "evalue");
        let bitscore = parse_f64!(11, "bitscore");
        let qlen     = parse_u64!(12, "qlen");
        let slen     = parse_u64!(13, "slen");

        if tab_pos.len() > 13 {
            warn!("Extra columns in NT m8 line (expected 14)");
        }
        Ok(Self { qname, tname, pident, alen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore, qlen, slen })
    }

    #[cfg(target_arch = "x86_64")]
    fn parse_line_nt_avx512(line: &str) -> Result<Self> {
        // Safety: only reachable when SIMD_LEVEL == Avx512, confirmed at startup.
        unsafe { Self::parse_line_nt_avx512_inner(line) }
    }

    // ── NR AVX-512 ─────────────────────────────────────────────────────────

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512f,avx512bw")]
    unsafe fn parse_line_nr_avx512_inner(line: &str) -> Result<Self> {
        let bytes = line.trim_end().as_bytes();
        if bytes.is_empty() { return Err(anyhow!("empty line")); }

        let tab_pos = Self::collect_tab_positions(bytes);
        let len = bytes.len();

        if tab_pos.len() < 11 {
            return Err(anyhow!("NR m8 line has {} tabs, need ≥11 for 12 fields", tab_pos.len()));
        }

        macro_rules! field {
            ($n:expr) => {{
                let start = if $n == 0 { 0 } else { tab_pos[$n - 1] + 1 };
                let end   = if $n < tab_pos.len() { tab_pos[$n] } else { len };
                &bytes[start..end]
            }};
        }
        macro_rules! parse_u64 {
            ($n:expr, $name:literal) => {
                lexical_parse::<u64, _>(field!($n))
                    .map_err(|e| anyhow!(concat!($name, ": {}"), e))?
            };
        }
        macro_rules! parse_f64 {
            ($n:expr, $name:literal) => {
                lexical_parse::<f64, _>(field!($n))
                    .map_err(|e| anyhow!(concat!($name, ": {}"), e))?
            };
        }

        let qname = std::str::from_utf8(field!(0)).map_err(|_| anyhow!("qname not UTF-8"))?.to_string();
        let raw   = std::str::from_utf8(field!(1)).map_err(|_| anyhow!("tname not UTF-8"))?;
        let tname = raw.split('.').next().unwrap_or(raw).to_string();

        let pident   = parse_f64!(2,  "pident");
        let alen     = parse_u64!(3,  "alen");
        let mismatch = parse_u64!(4,  "mismatch");
        let gapopen  = parse_u64!(5,  "gapopen");
        let qstart   = parse_u64!(6,  "qstart");
        let qend     = parse_u64!(7,  "qend");
        let tstart   = parse_u64!(8,  "tstart");
        let tend     = parse_u64!(9,  "tend");
        let evalue   = parse_f64!(10, "evalue");
        let bitscore = parse_f64!(11, "bitscore");

        if tab_pos.len() > 11 {
            warn!("Extra columns in NR m8 line (expected 12)");
        }
        Ok(Self { qname, tname, pident, alen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore, qlen: 0, slen: 0 })
    }

    #[cfg(target_arch = "x86_64")]
    fn parse_line_nr_avx512(line: &str) -> Result<Self> {
        // Safety: only reachable when SIMD_LEVEL == Avx512, confirmed at startup.
        unsafe { Self::parse_line_nr_avx512_inner(line) }
    }

    // ── public dispatch ────────────────────────────────────────────────────

    /// Parse 14-column NT (blastn) output — dispatches to best available path.
    pub fn parse_line_nt(line: &str) -> Result<Self> {
        PARSE_M8_NT(line)
    }

    /// Parse 12-column NR (blastx) output — dispatches to best available path.
    pub fn parse_line_nr(line: &str) -> Result<Self> {
        PARSE_M8_NR(line)
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

static PARSE_M8_NT: Lazy<fn(&str) -> Result<M8Record>> = Lazy::new(|| {
    #[cfg(target_arch = "x86_64")]
    if matches!(*SIMD_LEVEL, SimdLevel::Avx512) {
        debug!("M8 NT parser: using AVX-512 path");
        return M8Record::parse_line_nt_avx512;
    }
    debug!("M8 NT parser: using scalar path");
    M8Record::parse_line_nt_scalar
});

static PARSE_M8_NR: Lazy<fn(&str) -> Result<M8Record>> = Lazy::new(|| {
    #[cfg(target_arch = "x86_64")]
    if matches!(*SIMD_LEVEL, SimdLevel::Avx512) {
        debug!("M8 NR parser: using AVX-512 path");
        return M8Record::parse_line_nr_avx512;
    }
    debug!("M8 NR parser: using scalar path");
    M8Record::parse_line_nr_scalar
});

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



pub fn process_record_pair(
    m8_bytes: &[u8],
    hit_bytes: &[u8],
    agg: &mut AHashMap<Vec<i32>, AggBucket>,
    lineage_cache: &mut AHashMap<i32, Vec<i32>>,
    duplicate_clusters: &DashMap<String, ClusterInfo>,
    should_keep: &dyn Fn(&[i32]) -> bool,
    count_type: &str,
    source_count_type: Option<&str>,
) -> Result<(), anyhow::Error> {
    let hit_str = std::str::from_utf8(hit_bytes)?;
    let hit_fields: Vec<&str> = hit_str.trim_end().split('\t').collect();
    if hit_fields.len() < 7 {
        return Ok(()); // malformed hit summary line → skip (matches Python)
    }

    let read_id = hit_fields[0].to_string();
    let level: u8 = hit_fields.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
    let consensus_taxid: i32 = hit_fields.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);

    if consensus_taxid <= 0 || level == 0 {
        return Ok(());
    }

    let m8_str = std::str::from_utf8(m8_bytes)?;
    let m8 = M8Record::parse_line_nr(m8_str)?;

    let mut alen = m8.alen as f64;
    let pident = m8.pident;
    let mut ev = m8.evalue;
    if ev <= MIN_NORMAL_POSITIVE_DOUBLE {
        ev = MIN_NORMAL_POSITIVE_DOUBLE;
    }
    let ev_log10 = ev.log10();

    // Exact Python merged_NT_NR adjustment (still needed)
    if count_type == "merged_NT_NR" && source_count_type == Some("NR") {
        alen *= 3.0;
    }

    // Build raw lineage from hit summary columns (indices 4,5,6)
    let species = hit_fields.get(4).and_then(|s| s.parse().ok()).unwrap_or(-1);
    let genus   = hit_fields.get(5).and_then(|s| s.parse().ok()).unwrap_or(-1);
    let family  = hit_fields.get(6).and_then(|s| s.parse().ok()).unwrap_or(-1);
    let raw = vec![species, genus, family];

    // Lineage cache (exactly like summarize_hits)
    let cleaned = if let Some(cached) = lineage_cache.get(&consensus_taxid) {
        cached.clone()
    } else {
        let c = validate_taxid_lineage(&raw, consensus_taxid, level);
        lineage_cache.insert(consensus_taxid, c.clone());
        c
    };

    if !should_keep(&cleaned) {
        return Ok(());
    }

    let cluster_size = duplicate_clusters
        .get(&read_id)
        .map(|e| e.value().size)
        .unwrap_or(1);

    let mut agg_key = cleaned;
    while !agg_key.is_empty() {
        let bucket = agg.entry(agg_key.clone()).or_default();

        bucket.nonunique_count += cluster_size;
        bucket.unique_count += 1;
        bucket.base_count += (m8.alen as u64) * cluster_size;
        bucket.sum_percent_identity += pident;
        bucket.sum_alignment_length += alen;
        bucket.sum_e_value += ev_log10;

        if let Some(src) = source_count_type {
            bucket.source_count_type.insert(src.to_string());
        }

        agg_key = agg_key[1..].to_vec();
    }

    Ok(())
}

pub fn merge_aggregations(parts: Vec<AHashMap<Vec<i32>, AggBucket>>) -> AHashMap<Vec<i32>, AggBucket> {
    let mut merged = AHashMap::new();
    for part in parts {
        for (k, b) in part {
            let e: &mut AggBucket = merged.entry(k).or_default();

            e.nonunique_count += b.nonunique_count;
            e.unique_count += b.unique_count;
            e.base_count += b.base_count;
            e.sum_percent_identity += b.sum_percent_identity;
            e.sum_alignment_length += b.sum_alignment_length;
            e.sum_e_value += b.sum_e_value;

            if !b.source_count_type.is_empty() {
                e.source_count_type.extend(b.source_count_type);
            }
        }
    }
    merged
}

pub fn build_taxon_counts_list(agg: AHashMap<Vec<i32>, AggBucket>, count_type: &str) -> Vec<TaxonCount> {
    let mut rows = Vec::with_capacity(agg.len());
    for (agg_key, bucket) in agg {
        if bucket.unique_count == 0 { continue; }

        let len_key = agg_key.len();
        let tax_level = (3 - len_key + 1) as u8;   // exact Python: species=1, genus=2, family=3

        let genus_taxid = if tax_level <= 2 { agg_key[2 - tax_level as usize] } else { -200 };
        let family_taxid = if tax_level <= 3 { agg_key[3 - tax_level as usize] } else { -300 };

        let count = if READ_COUNTING_MODE == ReadCountingMode::CountAll {
            bucket.nonunique_count
        } else {
            bucket.unique_count
        };

        let source_vec = if bucket.source_count_type.is_empty() {
            None
        } else {
            let mut v: Vec<_> = bucket.source_count_type.into_iter().collect();
            v.sort_unstable();
            Some(v)
        };

        rows.push(TaxonCount {
            tax_id: agg_key[0],
            tax_level,
            genus_taxid,
            family_taxid,
            count,
            nonunique_count: bucket.nonunique_count,
            unique_count: bucket.unique_count,
            dcr: bucket.nonunique_count as f64 / bucket.unique_count as f64,
            percent_identity: bucket.sum_percent_identity / bucket.unique_count as f64,
            alignment_length: bucket.sum_alignment_length / bucket.unique_count as f64,
            e_value: bucket.sum_e_value / bucket.unique_count as f64,
            count_type: count_type.to_string(),
            base_count: bucket.base_count,
            source_count_type: source_vec,
        });
    }
    rows
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
pub fn consensus_level(
    hits: &[M8Record],
    lineage_map: &AHashMap<Taxid, Lineage>,
    acc2taxid_map: &Map<Vec<u8>>,
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



pub fn read_id_from_m8_line(line: &str) -> Option<&str> {
    let read_id = line.split('\t').next()?;
    if read_id.is_empty() { None } else { Some(read_id) }
}

pub fn shard_for_read_id(read_id: &str, workers: usize) -> usize {
    let workers = workers.max(1);
    let state = AHashState::with_seeds(0, 1, 2, 3);
    let mut hasher = state.build_hasher();
    read_id.hash(&mut hasher);
    (hasher.finish() as usize) % workers
}

pub fn summarize_m8_hits<F>(
    seq: u64,
    read_id: String,
    hits: &[M8Record],
    lineage_map: &AHashMap<Taxid, Lineage>,
    acc2taxid_map: &Map<Vec<u8>>,
    should_keep_filter: &F,
    min_aln_len: u64,
) -> Option<ReducedRead>
where
    F: Fn(&[i32]) -> bool + Send + Sync,
{
    // Preserve existing semantics: ignore short hits, invalid accessions, invalid lineages,
    // and only keep hits that pass should_keep_filter.
    let mut valid_hits: Vec<&M8Record> = Vec::with_capacity(hits.len().min(64));

    for hit in hits {
        if hit.alen < min_aln_len {
            continue;
        }

        if let Some(taxid_u64) = acc2taxid_map.get(hit.tname.as_bytes()) {
            let taxid = taxid_u64 as i32;
            if taxid <= 0 {
                continue;
            }

            let lineage = lineage_map.get(&taxid).cloned().unwrap_or([-1i32; 3]);
            if !should_keep_filter(&lineage) {
                continue;
            }

            valid_hits.push(hit);
        }
    }

    if valid_hits.is_empty() {
        return None;
    }

    // Best hit = max bitscore, first seen wins on ties.
    let mut best: Option<&M8Record> = None;
    let mut max_bitscore = f64::NEG_INFINITY;

    for &h in &valid_hits {
        if h.bitscore > max_bitscore {
            max_bitscore = h.bitscore;
            best = Some(h);
        }
    }

    let best = best?;

    let dedup_line = format!(
        "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.3}\n",
        best.qname,
        best.tname,
        best.pident,
        best.alen,
        best.mismatch,
        best.gapopen,
        best.qstart,
        best.qend,
        best.tstart,
        best.tend,
        best.evalue,
        best.bitscore
    )
        .into_bytes();

    // Keep the same consensus logic as the current code.
    let valid_owned: Vec<M8Record> = valid_hits.into_iter().cloned().collect();
    let wrapped_filter = Arc::new(|l: &[i32]| should_keep_filter(l));

    let (tax_level, _cons_taxid, consensus_hits) =
        consensus_level(&valid_owned, lineage_map, acc2taxid_map, &wrapped_filter)
            .unwrap_or((0, 0, Vec::new()));

    let first_lineage = consensus_hits
        .first()
        .and_then(|h| acc2taxid_map.get(h.tname.as_bytes()))
        .map(|taxid_u64| taxid_u64 as i32)
        .and_then(|taxid| lineage_map.get(&taxid).cloned())
        .unwrap_or([-1i32; 3]);

    let summary_line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\n",
        read_id,
        best.tname,
        first_lineage[0],
        first_lineage[1],
        first_lineage[2],
        tax_level
    )
        .into_bytes();

    Some(ReducedRead {
        seq,
        dedup: dedup_line,
        summary: summary_line,
        accession: best.tname.clone(),   
    })
}



pub fn parse_summary_batch(batch: Vec<u8>) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && line[0] != b'#')
        .flat_map(|line: &[u8]| vec![line.to_vec()])
        .collect()
}

pub fn parse_m8_metrics_batch(batch: Vec<u8>) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && line[0] != b'#')
        .flat_map(|line: &[u8]| vec![line.to_vec()])
        .collect()
}

pub fn parse_m8_acc_batch(batch: Vec<u8>) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && line[0] != b'#')
        .flat_map(|line_bytes: &[u8]| {
            // Strip trailing whitespace in bytes — no UTF-8 conversion needed
            let line_bytes = match line_bytes.iter().rposition(|&b| b > b' ') {
                Some(last) => &line_bytes[..=last],
                None => return vec![],
            };

            // Find the first tab → field 0 is read_id, field 1 starts after it
            let tab1 = match memchr(b'\t', line_bytes) {
                Some(pos) => pos,
                None => return vec![],
            };
            let read_id = &line_bytes[..tab1];

            // Find the second tab (or use end of line) → field 1 is accession
            let rest = &line_bytes[tab1 + 1..];
            let accession = match memchr(b'\t', rest) {
                Some(pos) => &rest[..pos],
                None => rest,
            };

            if read_id.is_empty() || accession.is_empty() {
                return vec![];
            }

            let mut out = Vec::with_capacity(read_id.len() + accession.len() + 2);
            out.extend_from_slice(read_id);
            out.push(b'\t');
            out.extend_from_slice(accession);
            out.push(b'\n');
            vec![out]
        })
        .collect()
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
    mut m8_stream_rx:           ReceiverStream<ParseOutput>,
    mut hit_summary_stream_rx:  ReceiverStream<ParseOutput>,
    db_type:                    &str,
    lineage_map:                Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter:         Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    duplicate_clusters:         Arc<DashMap<String, ClusterInfo>>,
    output_tx:                  mpsc::Sender<ParseOutput>,
    concurrency:                usize,
    batch_size_lines:           usize,
) -> Result<()> {
    if concurrency == 0 {
        return Err(anyhow::anyhow!("concurrency must be > 0"));
    }

    let db_type = db_type.to_string();

    // Bounded mpsc channel — hard concurrency cap + backpressure
    let (job_tx, job_rx) = mpsc::channel::<(Vec<String>, Vec<String>)>(concurrency);
    let shared_rx = Arc::new(tokio::sync::Mutex::new(job_rx));

    // Producer: build batches and send to bounded channel (exact same loop as before)
    let producer_handle = tokio::spawn({
        let mut m8_stream = m8_stream_rx;
        let mut hit_stream = hit_summary_stream_rx;
        async move {
            let mut batch_m8 = Vec::with_capacity(batch_size_lines);
            let mut batch_hit = Vec::with_capacity(batch_size_lines);

            loop {
                let m8_item = m8_stream.next().await;
                let hit_item = hit_stream.next().await;

                match (m8_item, hit_item) {
                    (Some(m8), Some(hit)) => {
                        if let Ok(bytes) = m8.to_bytes() {
                            let line = String::from_utf8_lossy(&bytes);
                            let trimmed = line.trim_end().to_string();
                            if !trimmed.is_empty() {
                                batch_m8.push(trimmed);
                            }
                        }
                        if let Ok(bytes) = hit.to_bytes() {
                            let line = String::from_utf8_lossy(&bytes);
                            let trimmed = line.trim_end().to_string();
                            if !trimmed.is_empty() {
                                batch_hit.push(trimmed);
                            }
                        }

                        if batch_m8.len() >= batch_size_lines && batch_hit.len() >= batch_size_lines {
                            if job_tx.send((std::mem::take(&mut batch_m8), std::mem::take(&mut batch_hit))).await.is_err() {
                                break;
                            }
                        }
                    }
                    _ => {
                        // final partial batch
                        if !batch_m8.is_empty() || !batch_hit.is_empty() {
                            let _ = job_tx.send((batch_m8, batch_hit)).await;
                        }
                        break;
                    }
                }
            }
            drop(job_tx);
            Ok::<(), anyhow::Error>(())
        }
    });

    // Fixed worker pool — exactly `concurrency` workers (hard cap)
    let mut workers = Vec::with_capacity(concurrency);
    for _ in 0..concurrency {
        let rx = shared_rx.clone();
        let output_tx = output_tx.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let db_type = db_type.clone();

        let worker = tokio::spawn(async move {
            loop {
                let job = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };

                if let Some((batch_m8, batch_hit)) = job {
                    let mut buckets: AHashMap<Taxid, AggBucket> = AHashMap::new();

                    for (m8_line, hit_line) in batch_m8.into_iter().zip(batch_hit) {
                        let hit_fields: Vec<&str> = hit_line.split('\t').collect();
                        if hit_fields.len() < 10 { continue; }

                        let read_id = hit_fields[0].to_string();
                        let level   = hit_fields.get(1).and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
                        let taxid   = hit_fields.get(2).and_then(|s| s.parse::<i32>().ok()).unwrap_or(0);

                        if taxid <= 0 || level == 0 { continue; }

                        if let Ok(m8) = M8Record::parse_line_nt(&m8_line)
                            .or_else(|_| M8Record::parse_line_nr(&m8_line))
                        {
                            let bucket = buckets.entry(taxid).or_default();

                            bucket.nonunique_count += 1;
                            bucket.unique_count += duplicate_clusters
                                .get(&read_id)
                                .map(|entry| entry.value().size)
                                .unwrap_or(1);

                            bucket.base_count += 1;
                            bucket.sum_percent_identity += m8.pident;
                            bucket.sum_alignment_length += m8.alen as f64;
                            bucket.sum_e_value += m8.evalue;
                            bucket.source_count_type.insert(db_type.clone());
                        }
                    }

                    // Emit results
                    for (taxid, bucket) in buckets {
                        if let Some(lineage) = lineage_map.get(&taxid) {
                            if !should_keep_filter(lineage) { continue; }

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
                                tax_level: 1,
                                genus_taxid: lineage[1],
                                family_taxid: lineage[2],
                                count: bucket.unique_count,
                                nonunique_count: bucket.nonunique_count,
                                unique_count: bucket.unique_count,
                                dcr,
                                percent_identity,
                                alignment_length,
                                e_value,
                                count_type: db_type.clone(),
                                base_count: bucket.base_count,
                                source_count_type: Some(bucket.source_count_type.into_iter().collect()),
                            };

                            let json = serde_json::to_string(&count)? + "\n";
                            let _ = output_tx
                                .send(ParseOutput::Bytes(Bytes::from(json)))
                                .await;
                        }
                    }
                } else {
                    break;
                }
            }
            Ok::<(), anyhow::Error>(())
        });
        workers.push(worker);
    }

    // Wait for producer + all workers
    producer_handle.await??;
    for w in workers {
        w.await??;
    }

    info!("Finished taxon count JSON generation for {} (bounded channel, {} workers)", db_type, concurrency);
    Ok(())
}


// pub async fn compute_merged_taxon_counts(
//     config: Arc<RunConfig>,
//     nt_m8_stream: mpsc::Receiver<ParseOutput>,
//     nt_hit_summary_stream: mpsc::Receiver<ParseOutput>,
//     nt_contig_summary: Vec<ContigSummaryEntry>,
//
//     nr_m8_stream: mpsc::Receiver<ParseOutput>,
//     nr_hit_summary_stream: mpsc::Receiver<ParseOutput>,
//     nr_contig_summary: Vec<ContigSummaryEntry>,
//
//     lineage_map: Arc<AHashMap<Taxid, Lineage>>,
//     should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
//     duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
//
//     merged_m8_path: PathBuf,
//     merged_hitsummary_path: PathBuf,
//     merged_taxon_counts_path: PathBuf,
//     merged_contig_summary_path: PathBuf,
//     nr_alignment_per_read: Arc<DashMap<String, SpeciesAlignmentResults, AHashRandomState>>,
// ) -> Result<()> {
//     // Write merged m8 and hit summary to disk
//     let mut merged_m8_file = BufWriter::new(TokioFile::create(&merged_m8_path).await?);
//     let mut merged_hit_file = BufWriter::new(TokioFile::create(&merged_hitsummary_path).await?);
//
//     // Channels for synchronized writing
//     let (m8_write_tx, mut m8_write_rx) = mpsc::channel::<Arc<Vec<u8>>>(1024);
//     let (hit_write_tx, mut hit_write_rx) = mpsc::channel::<Arc<Vec<u8>>>(1024);
//
//     let m8_write_task = tokio::spawn(async move {
//         while let Some(bytes) = m8_write_rx.recv().await {
//             merged_m8_file.write_all(&bytes).await?;
//             merged_m8_file.write_all(b"\n").await?;
//         }
//         merged_m8_file.flush().await?;
//         Ok::<(), anyhow::Error>(())
//     });
//
//     let hit_write_task = tokio::spawn(async move {
//         while let Some(bytes) = hit_write_rx.recv().await {
//             merged_hit_file.write_all(&bytes).await?;
//             merged_hit_file.write_all(b"\n").await?;
//         }
//         merged_hit_file.flush().await?;
//         Ok::<(), anyhow::Error>(())
//     });
//
//     // NT pass task
//     let nt_m8_write_tx = m8_write_tx.clone();
//     let nt_hit_write_tx = hit_write_tx.clone();
//     let nr_align_nt = nr_alignment_per_read.clone();
//     let nt_task = tokio::spawn(async move {
//         let mut nt_m8_stream = ReceiverStream::new(nt_m8_stream);
//         let mut nt_hit_stream = ReceiverStream::new(nt_hit_summary_stream);
//
//         loop {
//             let m8_item = nt_m8_stream.next().await;
//             let hit_item = nt_hit_stream.next().await;
//
//             match (m8_item, hit_item) {
//                 (Some(m8_item), Some(hit_item)) => {
//                     let m8_bytes = if let ParseOutput::Bytes(bytes) = m8_item {
//                         bytes
//                     } else {
//                         return Err(anyhow!("Unexpected item type in NT m8 stream: expected Bytes"));
//                     };
//                     let hit_bytes = hit_item.to_bytes()?;
//                     let hit_line = String::from_utf8_lossy(&hit_bytes);
//                     let hit_trimmed = hit_line.trim_end();
//                     if hit_trimmed.is_empty() {
//                         continue;
//                     }
//
//                     let hit_fields: Vec<&str> = hit_trimmed.split('\t').collect();
//                     if hit_fields.len() < 10 {
//                         warn!("Malformed NT hit line: {}", hit_trimmed);
//                         continue;
//                     }
//
//                     let read_id = hit_fields[0];
//
//                     let nt_contig = hit_fields.get(9).and_then(|s| s.parse::<Taxid>().ok());
//                     let nt_read = hit_fields.get(3).and_then(|s| s.parse::<Taxid>().ok());
//
//                     let nr_align = nr_align_nt.get(read_id);
//                     let has_nr_contig = nr_align.as_deref().map_or(false, |a| a.contig.is_some());
//                     let has_nr_read = nr_align.as_deref().map_or(false, |a| a.read.is_some());
//
//                     if nt_contig.is_some() || (!has_nr_contig && nt_read.is_some()) {
//                         nt_m8_write_tx.send(m8_bytes).await?;
//
//                         let mut hit_with_source = hit_fields.iter().copied().collect::<Vec<&str>>();
//                         hit_with_source.push(NT_TAG);
//                         nt_hit_write_tx.send(Arc::new(hit_with_source.join("\t").as_bytes().to_vec())).await?;
//
//                         nr_align_nt.remove(read_id);
//                     }
//                 }
//                 (None, None) => break,
//                 (Some(_), None) => return Err(anyhow!("NT m8 stream has more items than NT hit summary stream")),
//                 (None, Some(_)) => return Err(anyhow!("NT hit summary stream has more items than NT m8 stream")),
//             }
//         }
//         Ok::<(), anyhow::Error>(())
//     });
//
//     // NR pass task
//     let nr_m8_write_tx = m8_write_tx;
//     let nr_hit_write_tx = hit_write_tx;
//     let nr_align_nr = nr_alignment_per_read;
//     let nr_task = tokio::spawn(async move {
//         let mut nr_m8_stream = ReceiverStream::new(nr_m8_stream);
//         let mut nr_hit_stream = ReceiverStream::new(nr_hit_summary_stream);
//
//         loop {
//             let m8_item = nr_m8_stream.next().await;
//             let hit_item = nr_hit_stream.next().await;
//
//             match (m8_item, hit_item) {
//                 (Some(m8_item), Some(hit_item)) => {
//                     let hit_bytes = hit_item.to_bytes()?;
//                     let hit_line = String::from_utf8_lossy(&hit_bytes);
//                     let hit_trimmed = hit_line.trim_end();
//                     if hit_trimmed.is_empty() {
//                         continue;
//                     }
//
//                     let hit_fields: Vec<&str> = hit_trimmed.split('\t').collect();
//                     let read_id = hit_fields[0];
//
//                     if nr_align_nr.contains_key(read_id) {
//                         let m8_bytes = if let ParseOutput::Bytes(bytes) = m8_item {
//                             bytes
//                         } else {
//                             return Err(anyhow!("Unexpected item type in NR m8 stream: expected Bytes"));
//                         };
//                         nr_m8_write_tx.send(m8_bytes).await?;
//
//                         let mut hit_with_source = hit_fields.iter().copied().collect::<Vec<&str>>();
//                         hit_with_source.push(NR_TAG);
//                         nr_hit_write_tx.send(Arc::new(hit_with_source.join("\t").as_bytes().to_vec())).await?;
//                     }
//                 }
//                 (None, None) => break,
//                 (Some(_), None) => return Err(anyhow!("NR m8 stream has more items than NR hit summary stream")),
//                 (None, Some(_)) => return Err(anyhow!("NR hit summary stream has more items than NR m8 stream")),
//             }
//         }
//         Ok::<(), anyhow::Error>(())
//     });
//
//     nt_task.await??;
//     nr_task.await??;
//
//     m8_write_task.await??;
//     hit_write_task.await??;
//     info!("Merged alignment files written to disk");
//
//     // Stream merged files back for taxon counting
//     let (m8_tx, m8_rx) = mpsc::channel(4096);
//     let (hit_tx, hit_rx) = mpsc::channel(4096);
//
//     let m8_path_clone = merged_m8_path.clone();
//     let m8_task = tokio::spawn(async move {
//         let file = TokioFile::open(m8_path_clone).await?;
//         let mut reader = BufReader::new(file);
//         let mut line = String::new();
//         while reader.read_line(&mut line).await? > 0 {
//             if !line.trim_end().is_empty() {
//                 let _ = m8_tx.send(ParseOutput::Bytes(Arc::new(line.as_bytes().to_vec()))).await;
//             }
//             line.clear();
//         }
//         Ok::<(), anyhow::Error>(())
//     });
//
//     let hit_path_clone = merged_hitsummary_path.clone();
//     let hit_task = tokio::spawn(async move {
//         let file = TokioFile::open(hit_path_clone).await?;
//         let mut reader = BufReader::new(file);
//         let mut line = String::new();
//         while reader.read_line(&mut line).await? > 0 {
//             if !line.trim_end().is_empty() {
//                 let _ = hit_tx.send(ParseOutput::Bytes(Arc::new(line.as_bytes().to_vec()))).await;
//             }
//             line.clear();
//         }
//         Ok::<(), anyhow::Error>(())
//     });
//
//     let merged_taxon_counts_path_clone = merged_taxon_counts_path.clone();
//     let config_clone = config.clone();
//
//     let (json_tx, json_rx) = mpsc::channel(4096);
//     let write_task = tokio::spawn(async move {
//         write_byte_stream_to_file(
//             &merged_taxon_counts_path_clone,
//             ReceiverStream::new(json_rx),
//             Some(config_clone.base_buffer_size),
//             "compute_merged_taxon_counts"
//         ).await
//     });
//
//     let taxon_count_concurrency = compute_phase_concurrency(
//         &config,
//         "taxon_counting",
//         0.4,      // very light
//         4.0,
//         64,
//         4,
//     );
//
//     let taxon_count_batch_size = compute_batch_size(
//         None,               // we don't know total lines yet
//         220,                // rough avg bytes per m8 + hit line
//         150,
//         taxon_count_concurrency,
//     );
//
//     generate_taxon_count_json_from_m8(
//         ReceiverStream::new(m8_rx),
//         ReceiverStream::new(hit_rx),
//         "merged_NT_NR",
//         lineage_map.clone(),
//         should_keep_filter.clone(),
//         duplicate_clusters.clone(),
//         json_tx,
//         taxon_count_concurrency,
//         taxon_count_batch_size
//     )
//         .await?;
//
//     m8_task.await??;
//     hit_task.await??;
//     write_task.await??;
//     info!("Merged taxon counts generated");
//
//     // Contig summary merge
//     let mut merged_contigs: HashMap<Taxid, Vec<ContigSummaryEntry>> = HashMap::new();
//     let mut nt_contig_names: HashSet<String> = HashSet::new();
//
//     for mut entry in nt_contig_summary {
//         entry.db_type = "merged_NT_NR".to_string();
//         nt_contig_names.insert(entry.contig_name.clone());
//         merged_contigs.entry(entry.species_taxid).or_default().push(entry);
//     }
//
//     for mut entry in nr_contig_summary {
//         if !nt_contig_names.contains(&entry.contig_name) {
//             entry.db_type = "merged_NT_NR".to_string();
//             merged_contigs.entry(entry.species_taxid).or_default().push(entry);
//         }
//     }
//
//     let final_contigs: Vec<ContigSummaryEntry> = merged_contigs.into_values().flatten().collect();
//
//     let json = serde_json::to_string_pretty(&final_contigs)?;
//     let mut out = BufWriter::new(TokioFile::create(&merged_contig_summary_path).await?);
//     out.write_all(json.as_bytes()).await?;
//     out.flush().await?;
//     info!("Merged contig summary written ({} entries)", final_contigs.len());
//
//     info!("compute_merged_taxon_counts complete");
//     Ok(())
// }

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_m8_eq(a: &M8Record, b: &M8Record, ctx: &str) {
        assert_eq!(a.qname,    b.qname,    "{ctx}: qname");
        assert_eq!(a.tname,    b.tname,    "{ctx}: tname");
        assert!((a.pident - b.pident).abs() < 1e-9, "{ctx}: pident");
        assert_eq!(a.alen,     b.alen,     "{ctx}: alen");
        assert_eq!(a.mismatch, b.mismatch, "{ctx}: mismatch");
        assert_eq!(a.gapopen,  b.gapopen,  "{ctx}: gapopen");
        assert_eq!(a.qstart,   b.qstart,   "{ctx}: qstart");
        assert_eq!(a.qend,     b.qend,     "{ctx}: qend");
        assert_eq!(a.tstart,   b.tstart,   "{ctx}: tstart");
        assert_eq!(a.tend,     b.tend,     "{ctx}: tend");
        assert!((a.evalue   - b.evalue  ).abs() < 1e-30, "{ctx}: evalue");
        assert!((a.bitscore - b.bitscore).abs() < 1e-9,  "{ctx}: bitscore");
        assert_eq!(a.qlen,     b.qlen,     "{ctx}: qlen");
        assert_eq!(a.slen,     b.slen,     "{ctx}: slen");
    }

    fn compare_nt(line: &str) {
        let scalar = M8Record::parse_line_nt_scalar(line).expect("NT scalar failed");
        #[cfg(target_arch = "x86_64")]
        {
            let avx = M8Record::parse_line_nt_avx512(line).expect("NT avx512 failed");
            assert_m8_eq(&scalar, &avx, line);
        }
    }

    fn compare_nr(line: &str) {
        let scalar = M8Record::parse_line_nr_scalar(line).expect("NR scalar failed");
        #[cfg(target_arch = "x86_64")]
        {
            let avx = M8Record::parse_line_nr_avx512(line).expect("NR avx512 failed");
            assert_m8_eq(&scalar, &avx, line);
        }
    }

    // ── NT test lines (14 columns) ─────────────────────────────────────────

    const NT_BASIC: &str =
        "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903";

    const NT_ZERO_EVALUE: &str =
        "read2\tNC_000001.1\t100.000\t250\t0\t0\t1\t250\t1\t250\t0.0\t462.000\t250\t248956422";

    const NT_MISMATCHES: &str =
        "read3\tKJ660346.2\t95.238\t210\t10\t2\t5\t214\t500\t709\t4.56e-88\t330.000\t220\t18958";

    // Long enough to span two 64-byte SIMD chunks
    const NT_LONG_QNAME: &str =
        "a_very_long_read_identifier_that_exceeds_sixty_four_bytes_total\tNC_045512.2\t98.667\t150\t2\t0\t1\t150\t1\t150\t2.34e-70\t275.000\t150\t29903";

    const NT_TRAILING_WHITESPACE: &str =
        "read4\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903   ";

    // ── NR test lines (12 columns) ─────────────────────────────────────────

    const NR_BASIC: &str =
        "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000";

    const NR_ZERO_EVALUE: &str =
        "read2\tAAA12345.1\t100.000\t300\t0\t0\t1\t300\t1\t300\t0.0\t600.000";

    const NR_GAPS: &str =
        "read3\tXYZ99999.2\t78.431\t153\t33\t3\t10\t162\t5\t155\t8.9e-20\t89.700";

    // ── NT tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_nt_basic() { compare_nt(NT_BASIC); }

    #[test]
    fn test_nt_zero_evalue() { compare_nt(NT_ZERO_EVALUE); }

    #[test]
    fn test_nt_mismatches_gaps() { compare_nt(NT_MISMATCHES); }

    #[test]
    fn test_nt_long_qname_spans_two_chunks() { compare_nt(NT_LONG_QNAME); }

    #[test]
    fn test_nt_trailing_whitespace() {
        compare_nt(NT_TRAILING_WHITESPACE);
        let r = M8Record::parse_line_nt_scalar(NT_TRAILING_WHITESPACE).unwrap();
        assert_eq!(r.qname, "read4");
        assert_eq!(r.slen, 29903);
    }

    #[test]
    fn test_nt_accession_version_stripped() {
        let r = M8Record::parse_line_nt_scalar(NT_BASIC).unwrap();
        assert_eq!(r.tname, "NC_045512", "version suffix must be stripped");
        #[cfg(target_arch = "x86_64")]
        {
            let avx = M8Record::parse_line_nt_avx512(NT_BASIC).unwrap();
            assert_eq!(avx.tname, "NC_045512");
        }
    }

    #[test]
    fn test_nt_error_on_empty() {
        assert!(M8Record::parse_line_nt_scalar("").is_err());
        #[cfg(target_arch = "x86_64")]
        assert!(M8Record::parse_line_nt_avx512("").is_err());
    }

    #[test]
    fn test_nt_error_on_too_few_fields() {
        let short = "read1\tNC_045512.2\t99.333\t150";
        assert!(M8Record::parse_line_nt_scalar(short).is_err());
        #[cfg(target_arch = "x86_64")]
        assert!(M8Record::parse_line_nt_avx512(short).is_err());
    }

    // ── NR tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_nr_basic() { compare_nr(NR_BASIC); }

    #[test]
    fn test_nr_zero_evalue() { compare_nr(NR_ZERO_EVALUE); }

    #[test]
    fn test_nr_gaps() { compare_nr(NR_GAPS); }

    #[test]
    fn test_nr_qlen_slen_zero() {
        let r = M8Record::parse_line_nr_scalar(NR_BASIC).unwrap();
        assert_eq!(r.qlen, 0, "NR qlen must default to 0");
        assert_eq!(r.slen, 0, "NR slen must default to 0");
        #[cfg(target_arch = "x86_64")]
        {
            let avx = M8Record::parse_line_nr_avx512(NR_BASIC).unwrap();
            assert_eq!(avx.qlen, 0);
            assert_eq!(avx.slen, 0);
        }
    }

    #[test]
    fn test_nr_accession_version_stripped() {
        let r = M8Record::parse_line_nr_scalar(NR_BASIC).unwrap();
        assert_eq!(r.tname, "QIK02963");
        #[cfg(target_arch = "x86_64")]
        {
            let avx = M8Record::parse_line_nr_avx512(NR_BASIC).unwrap();
            assert_eq!(avx.tname, "QIK02963");
        }
    }

    #[test]
    fn test_nr_error_on_empty() {
        assert!(M8Record::parse_line_nr_scalar("").is_err());
        #[cfg(target_arch = "x86_64")]
        assert!(M8Record::parse_line_nr_avx512("").is_err());
    }

    #[test]
    fn test_nr_error_on_too_few_fields() {
        let short = "read1\tQIK02963.1\t87.500";
        assert!(M8Record::parse_line_nr_scalar(short).is_err());
        #[cfg(target_arch = "x86_64")]
        assert!(M8Record::parse_line_nr_avx512(short).is_err());
    }
}
