// BLAST-related file functions and structures
use std::collections::HashSet;
use std::sync::Arc;
use std::hash::{BuildHasher, Hash, Hasher};

use rayon::prelude::*;
use memchr::memchr;
use anyhow::{anyhow, Result};
use log::{self, debug, info, warn};
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;
use fst::Map;
use ahash::{AHashMap, RandomState as AHashState};
use serde::{Deserialize, Serialize};
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::mpsc;
use tokio::io::AsyncWriteExt;
use dashmap::DashMap;
use bytes::Bytes;
use once_cell::sync::Lazy;
use ahash::RandomState as AHashRandomState;

use crate::config::defs::{Taxid, Lineage, ClusterInfo, MIN_NORMAL_POSITIVE_DOUBLE, READ_COUNTING_MODE, ReadCountingMode, SIMD_LEVEL, SimdLevel, PipelineError, RunConfig, StreamDataType, SORT_TAG, NT_TAG, NR_TAG};
use crate::utils::streams::{ParseOutput, parse_lines, fanout_to_channels, ToBytes, spawn_cmd};
use crate::utils::taxonomy::{validate_taxid_lineage};
use crate::utils::file::{choose_temp_dir, write_byte_stream_to_file, rename_file_path};
use crate::utils::command::generate_cli;
use tokio::fs::File as TokioFile;
use tokio::time::Instant;
use tokio::sync::oneshot;
use std::path::PathBuf;
use std::collections::HashMap;
use futures::future::try_join_all;
use tokio::io::BufWriter;
use tokio::task::JoinHandle;
use crate::utils::system::compute_phase_concurrency;

/// Represents an entry in the contig summary, grouping metadata and hit stats.
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

/// Stores the best taxonomic alignment results for both a contig and its constituent read.
#[derive(Debug, Clone, Default)]
pub struct SpeciesAlignmentResults {
    pub contig: Option<Taxid>,
    pub read: Option<Taxid>,
}

/// Represents a single record in a BLAST m8 (tabular) output file.
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

/// Internal structure for grouped batch processing of read hits.
#[derive(Debug, Clone)]
pub struct ReducedRead {
    pub seq: u64,
    pub dedup: Vec<u8>,
    pub summary: Vec<u8>,
    pub accession: String,
}

/// Holds hits for a read that are waiting to be processed as a group.
#[derive(Debug)]
pub struct PendingRead {
    pub read_id: String,
    pub hits: Vec<M8Record>,
}

/// Message type for passing data to background workers.
#[derive(Debug)]
pub enum WorkerMsg {
    /// A single line of data with an associated sequence number.
    Line { seq: u64, line: Vec<u8> },
    /// Signals the worker to flush its current buffers.
    Flush,
}

impl M8Record {
    // ── NT scalar ──────────────────────────────────────────────────────────

    /// Parses a 14-column NT (blastn) output line using a scalar implementation.
    ///
    /// # Arguments
    ///
    /// * `line`: the m8 line string to parse
    ///
    /// # Returns
    ///
    /// Result<Self>: the parsed M8Record or an error
    pub fn parse_line_nt_scalar(line: &str) -> Result<Self> {
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

    /// Parses a 12-column NR (blastx) output line using a scalar implementation.
    ///
    /// # Arguments
    ///
    /// * `line`: the m8 line string to parse
    ///
    /// # Returns
    ///
    /// Result<Self>: the parsed M8Record or an error
    pub fn parse_line_nr_scalar(line: &str) -> Result<Self> {
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
        let trimmed = line.trim_end();
        let bytes = trimmed.as_bytes();
        if bytes.is_empty() {
            return Err(anyhow!("empty line"));
        }

        let tab_pos = Self::collect_tab_positions(bytes);
        let len = bytes.len();

        if tab_pos.len() < 13 {
            return Err(anyhow!(
            "NT m8 line has {} tabs, need ≥13 for 14 fields",
            tab_pos.len()
        ));
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
            warn!("Extra columns in NT m8 line (expected 14): {}", trimmed);
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
        let trimmed = line.trim_end();
        let bytes = trimmed.as_bytes();
        if bytes.is_empty() {
            return Err(anyhow!("empty line"));
        }

        let tab_pos = Self::collect_tab_positions(bytes);
        let len = bytes.len();

        if tab_pos.len() < 11 {
            return Err(anyhow!(
            "NR m8 line has {} tabs, need ≥11 for 12 fields",
            tab_pos.len()
        ));
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
            warn!("Extra columns in NT m8 line (expected 12): {}", trimmed);
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
            "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}\t{}\t{}",
            self.qname, self.tname, self.pident, self.alen,
            self.mismatch, self.gapopen, self.qstart, self.qend,
            self.tstart, self.tend, self.evalue, self.bitscore,
            self.qlen, self.slen
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



/// Processes a pair of m8 hit and hit summary records and updates the aggregation map.
///
/// # Arguments
///
/// * `m8_bytes`: raw bytes of the m8 alignment record
/// * `hit_bytes`: raw bytes of the corresponding hit summary record
/// * `agg`: shared map of lineage keys to aggregation buckets
/// * `lineage_cache`: cache for validated taxonomic lineages
/// * `duplicate_clusters`: map of cluster information for weighting
/// * `should_keep`: filter for taxonomic lineages to keep
/// * `count_type`: type of database being referenced (e.g., "NT")
/// * `source_count_type`: optional second source type (e.g., "NR")
///
/// # Returns
///
/// Result<()>: success or error
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

/// Merges multiple aggregation maps into a single one.
///
/// # Arguments
///
/// * `parts`: a vector of aggregation maps to merge
///
/// # Returns
///
/// AHashMap<Vec<i32>, AggBucket>: the merged aggregation map
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

/// Converts an aggregation map into a list of TaxonCount structures.
///
/// # Arguments
///
/// * `agg`: the aggregation map to convert
/// * `count_type`: the database type for the counts
///
/// # Returns
///
/// Vec<TaxonCount>: the list of generated taxon count records
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
#[allow(dead_code)]
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



/// Extract the read ID from a BLAST m8 line.
///
/// # Arguments
///
/// * `line`: the m8 line string
///
/// # Returns
///
/// Option<&str>: the extracted read ID if found
pub fn read_id_from_m8_line(line: &str) -> Option<&str> {
    let read_id = line.split('\t').next()?;
    if read_id.is_empty() { None } else { Some(read_id) }
}

/// Assigns a read ID to a specific shard for concurrent processing.
///
/// # Arguments
///
/// * `read_id`: the read ID string
/// * `workers`: total number of workers/shards
///
/// # Returns
///
/// usize: the assigned shard index
pub fn shard_for_read_id(read_id: &str, workers: usize) -> usize {
    let workers = workers.max(1);
    let state = AHashState::with_seeds(0, 1, 2, 3);
    let mut hasher = state.build_hasher();
    read_id.hash(&mut hasher);
    (hasher.finish() as usize) % workers
}

/// Sorts an BLAST m8 format byte stream by read ID (first column) using GNU sort.
/// RAM-first via choose_temp_dir, falls back to NVMe scratch.
/// Guarantees consecutive lines per read ID. Never drops data.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of unsorted m8 records
/// * `db_type`: database type (e.g., "nt", "nr") for tagging
///
/// # Returns
///
/// Result containing the sorted m8 record stream
pub async fn sort_m8_by_read_id(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    db_type: &str,
) -> Result<ReceiverStream<ParseOutput>, PipelineError> {
    let tag = format!("sort_m8_{}", db_type);
    let estimated_bytes = (config.input_size * 4).max(1_000_000_000);

    let temp_dir = choose_temp_dir(
        estimated_bytes,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    let unsorted_path = temp_dir.path().join(format!("{}_unsorted.m8", db_type));
    let sorted_path = temp_dir.path().join(format!("{}_sorted_by_read.m8", db_type));

    info!(
        "[{}] Starting m8 sort (est {} bytes) using RAM/NVMe scratch",
        tag, estimated_bytes
    );

    // 1) Write the unsorted stream to disk and wait for the file to be fully closed.
    let write_task = write_byte_stream_to_file(
        &unsorted_path,
        input_stream,
        config.clone(),
        StreamDataType::JustBytes,
        &tag,
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    write_task
        .await
        .map_err(|e| PipelineError::Other(anyhow!("{} writer task panicked: {}", tag, e)))?
        .map_err(|e| PipelineError::IOError(format!("{} writer task failed: {}", tag, e)))?;

    let sort_buffer = if config.available_ram > 256 * 1024 * 1024 * 1024 {
        "32G"
    } else {
        "8G"
    };

    // 2) Only now spawn GNU sort, so it reads a complete file.
    let sort_config = crate::utils::command::sort::SortConfig {
        key: "-k1,1".to_string(),
        parallel: None,
        buffer_size: Some(sort_buffer.to_string()),
        temp_dir: Some(temp_dir.path().to_string_lossy().into_owned()),
        output: Some(sorted_path.to_string_lossy().into_owned()),
        extra_fields: HashMap::new(),
        input: Some(unsorted_path.to_string_lossy().into_owned()),
    };

    let sort_args = generate_cli(
        crate::config::defs::SORT_TAG,
        &config,
        Some(&sort_config),
    )
        .map_err(|e| PipelineError::ToolExecution {
            tool: crate::config::defs::SORT_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut sort_child, sort_err_task) = spawn_cmd(
        config.clone(),
        crate::config::defs::SORT_TAG,
        sort_args,
        config.args.verbose,
        None
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: crate::config::defs::SORT_TAG.to_string(),
            error: e.to_string(),
        })?;

    let status = sort_child
        .wait()
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: crate::config::defs::SORT_TAG.to_string(),
            error: e.to_string(),
        })?;

    if !status.success() {
        let _ = sort_err_task.await;
        return Err(PipelineError::ToolExecution {
            tool: crate::config::defs::SORT_TAG.to_string(),
            error: format!("GNU sort exited with {:?}", status.code()),
        });
    }

    sort_err_task
        .await
        .map_err(|e| PipelineError::Other(anyhow!("sort stderr task panicked: {}", e)))?
        .map_err(|e| PipelineError::Other(anyhow!("sort stderr task failed: {}", e)))?;

    info!("[{}] m8 sort completed → {}", tag, sorted_path.display());

    // 3) Stream the sorted file back.
    let sorted_file = TokioFile::open(&sorted_path)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let meta = tokio::fs::metadata(&sorted_path)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    info!(
    "[{}] sorted output exists: {} bytes at {}",
    tag,
    meta.len(),
    sorted_path.display()
);


    let rx = parse_lines(sorted_file, &config, StreamDataType::JustBytes)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok(ReceiverStream::new(rx))
}

/// Calls taxonomic hits from a sorted m8 stream.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `m8_input`: stream of m8 records sorted by read ID
/// * `sample_base_buf`: base path for the sample
/// * `lineage_map`: map of taxids to lineages
/// * `acc2taxid_map`: map of accessions to taxids
/// * `should_keep_filter`: filter for taxids to keep
/// * `min_aln_len`: minimum alignment length to consider a hit
/// * `concurrency`: number of concurrent workers
/// * `tag`: tag for logging
///
/// # Returns
///
/// Result containing:
/// - stream of top m8 hits per read
/// - stream of hit summary records
/// - vector of cleanup tasks
/// - vector of cleanup receivers
pub async fn call_hits_m8(
    config: Arc<RunConfig>,
    mut m8_input: ReceiverStream<ParseOutput>, // MUST be sorted by read ID
    sample_base_buf: PathBuf,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    acc2taxid_map: Arc<Map<Vec<u8>>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    min_aln_len: u64,
    concurrency: usize,
    tag: String,
) -> Result<(
    ReceiverStream<ParseOutput>, // dedup m8 stream
    ReceiverStream<ParseOutput>, // summary stream
    Vec<tokio::task::JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
), PipelineError> {
    let worker_count = concurrency.max(1);
    let mut cleanup_tasks: Vec<tokio::task::JoinHandle<Result<()>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<Result<()>>> = Vec::new();

    info!(
        "[call_hits_m8:{}] STARTING STREAMING VERSION — workers={}, min_aln_len={}",
        tag, worker_count, min_aln_len
    );

    let (dedup_tx, dedup_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size * 256);
    let (summary_tx, summary_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size * 256);

    #[derive(Debug)]
    enum WorkerMsgLocal {
        ProcessRead {
            read_id: String,
            hits: Vec<M8Record>,
        },
        Flush,
    }

    // One sender per worker; reads are sharded by read_id.
    let mut worker_txs: Vec<mpsc::Sender<WorkerMsgLocal>> = Vec::with_capacity(worker_count);

    for worker_idx in 0..worker_count {
        let (worker_tx, mut worker_rx) = mpsc::channel::<WorkerMsgLocal>(config.base_buffer_size * 16);
        worker_txs.push(worker_tx);

        let lineage_map = lineage_map.clone();
        let acc2taxid_map = acc2taxid_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let dedup_tx = dedup_tx.clone();
        let summary_tx = summary_tx.clone();
        let worker_tag = tag.clone();

        let worker_handle = tokio::spawn(async move {
            let mut total_emitted = 0usize;
            let start = Instant::now();

            debug!(
                "[call_hits_m8:{}] worker {} started (streaming mode)",
                worker_tag, worker_idx
            );

            while let Some(msg) = worker_rx.recv().await {
                match msg {
                    WorkerMsgLocal::ProcessRead { read_id, hits } => {
                        let reduced = summarize_m8_hits(
                            0,
                            read_id,
                            &hits,
                            &lineage_map,
                            &acc2taxid_map,
                            &*should_keep_filter,
                            min_aln_len,
                        );

                        if let Some(reduced) = reduced {
                            total_emitted += 1;

                            dedup_tx
                                .send(ParseOutput::Bytes(Bytes::from(reduced.dedup)))
                                .await
                                .map_err(|e| anyhow!("dedup output channel closed: {}", e))?;

                            summary_tx
                                .send(ParseOutput::Bytes(Bytes::from(reduced.summary)))
                                .await
                                .map_err(|e| anyhow!("summary output channel closed: {}", e))?;
                        }
                    }
                    WorkerMsgLocal::Flush => break,
                }
            }

            info!(
                "[call_hits_m8:{}] worker {} finished — emitted {} reads in {:?}",
                worker_tag, worker_idx, total_emitted, start.elapsed()
            );
            Ok::<(), anyhow::Error>(())
        });

        cleanup_tasks.push(worker_handle);
    }

    // Coordinator: consume the sorted stream, accumulate one read at a time,
    // and dispatch finished read blocks to a shard.
    let coordinator_tag = tag.clone();
    let coordinator_handle = tokio::spawn(async move {
        let mut current_read_id: Option<String> = None;
        let mut current_hits: Vec<M8Record> = Vec::with_capacity(32);
        let mut total_lines = 0usize;
        let mut total_reads = 0usize;
        let start = Instant::now();

        #[cfg(debug_assertions)]
        let mut last_seen_read_id: Option<String> = None;

        debug!(
            "[call_hits_m8:{}] coordinator started (streaming group-by)",
            coordinator_tag
        );

        while let Some(item) = m8_input.next().await {
            let bytes = match item.to_bytes() {
                Ok(b) => b,
                Err(_) => continue,
            };

            let line_str = match std::str::from_utf8(&bytes) {
                Ok(s) => s.trim_end(),
                Err(_) => continue,
            };

            if line_str.is_empty() || line_str.starts_with('#') {
                continue;
            }

            total_lines += 1;

            let read_id = match read_id_from_m8_line(line_str) {
                Some(r) => r.to_string(),
                None => continue,
            };

            let rec = match M8Record::parse_line_nr(line_str) {
                Ok(r) => r,
                Err(_) => continue,
            };

            if total_lines == 1 {
                info!(
                    "[call_hits_m8:{}] first input line: read_id={} line={}",
                    coordinator_tag, read_id, line_str
                );
            }

            #[cfg(debug_assertions)]
            {
                if let Some(prev) = &last_seen_read_id {
                    if prev > &read_id {
                        warn!(
                            "[call_hits_m8:{}] read IDs appear to go backwards: prev={} current={}",
                            coordinator_tag, prev, read_id
                        );
                    }
                }
                last_seen_read_id = Some(read_id.clone());
            }

            // If the read ID changes, the previous read is complete.
            if current_read_id.as_deref() != Some(&read_id) && !current_hits.is_empty() {
                let prev_read_id = current_read_id.take().unwrap();
                let shard = shard_for_read_id(&prev_read_id, worker_count);

                worker_txs[shard]
                    .send(WorkerMsgLocal::ProcessRead {
                        read_id: prev_read_id,
                        hits: std::mem::take(&mut current_hits),
                    })
                    .await
                    .map_err(|e| anyhow!("failed to dispatch completed read to worker: {}", e))?;

                total_reads += 1;
            }

            if current_read_id.is_none() {
                current_read_id = Some(read_id.clone());
            }

            current_hits.push(rec);
        }

        // Final read.
        if let Some(read_id) = current_read_id {
            if !current_hits.is_empty() {
                let shard = shard_for_read_id(&read_id, worker_count);

                worker_txs[shard]
                    .send(WorkerMsgLocal::ProcessRead {
                        read_id,
                        hits: current_hits,
                    })
                    .await
                    .map_err(|e| anyhow!("failed to dispatch final read to worker: {}", e))?;

                total_reads += 1;
            }
        }

        // Flush all workers after all reads have been dispatched.
        for (_, tx) in worker_txs.iter().enumerate() {
            tx.send(WorkerMsgLocal::Flush)
                .await
                .map_err(|e| anyhow!("failed to flush worker: {}", e))?;
        }

        info!(
            "[call_hits_m8:{}] coordinator done — processed {} lines, {} reads in {:?}",
            coordinator_tag,
            total_lines,
            total_reads,
            start.elapsed()
        );

        if total_lines == 0 {
            warn!(
                "[call_hits_m8:{}] coordinator saw zero input lines; upstream sort/stream wiring is wrong",
                coordinator_tag
            );
        }

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(coordinator_handle);

    let dedup_stream = ReceiverStream::new(dedup_rx);
    let (dedup_branches, dedup_router_done_rx) = fanout_to_channels(
        dedup_stream,
        2,
        "call_hits_m8_dedup",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    // track router task (replaces cleanup_receivers)
    cleanup_receivers.push(dedup_router_done_rx);

    let mut dedup_branches = dedup_branches.into_iter();

    let dedup_main = ReceiverStream::new(
        dedup_branches.next().ok_or(PipelineError::EmptyStream)?
    );

    let dedup_file_stream = ReceiverStream::new(
        dedup_branches.next().ok_or(PipelineError::EmptyStream)?
    );

    let call_file_write_task = write_byte_stream_to_file(
        &config.out_dir.join(rename_file_path(
            &sample_base_buf,
            None,
            Some(&format!("{}.dedup.m8", tag)),
            "_",
        )),
        dedup_file_stream,
        config.clone(),
        StreamDataType::JustBytes,
        "call_hits_m8_dedup",
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    cleanup_tasks.push(call_file_write_task);

    let summary_stream = ReceiverStream::new(summary_rx);
    let (summary_branches, summary_router_done_rx) = fanout_to_channels(
        summary_stream,
        2,
        "call_hits_m8_summary",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    // track router task instead of oneshot receiver
    cleanup_receivers.push(summary_router_done_rx);

    let mut summary_branches = summary_branches.into_iter();

    let summary_main = ReceiverStream::new(
        summary_branches.next().ok_or(PipelineError::EmptyStream)?
    );

    let summary_file_stream = ReceiverStream::new(
        summary_branches.next().ok_or(PipelineError::EmptyStream)?
    );

    let summary_file_write_task = write_byte_stream_to_file(
        &config.out_dir.join(rename_file_path(
            &sample_base_buf,
            None,
            Some(&format!("{}.summary.txt", tag)),
            "_",
        )),
        summary_file_stream,
        config.clone(),
        StreamDataType::JustBytes,
        "call_hits_m8_summary",
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    cleanup_tasks.push(summary_file_write_task);

    info!("[call_hits_m8:{}] wiring complete — outputs ready", tag);

    Ok((
        dedup_main,
        summary_main,
        cleanup_tasks,
        cleanup_receivers,
    ))
}

/// Summarizes BLAST m8 hits for a single read into a ReducedRead structure.
///
/// # Arguments
///
/// * `seq`: sequence number of the read
/// * `read_id`: the read ID string
/// * `hits`: a slice of m8 records for the read
/// * `lineage_map`: map of taxids to taxonomic lineages
/// * `acc2taxid_map`: map of accessions to taxids
/// * `should_keep_filter`: filter for taxonomic lineages to keep
/// * `min_aln_len`: minimum alignment length to consider a hit
///
/// # Returns
///
/// Option<ReducedRead>: the summarized read info, or None if no hits pass filters
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



/// Parses a batch of summary records from raw bytes.
///
/// # Arguments
///
/// * `batch`: raw bytes containing summary lines
///
/// # Returns
///
/// Vec<Vec<u8>>: a vector of parsed summary line bytes
pub fn parse_summary_batch(batch: Vec<u8>) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && line[0] != b'#')
        .flat_map(|line: &[u8]| vec![line.to_vec()])
        .collect()
}

/// Parses a batch of m8 metrics from raw bytes.
///
/// # Arguments
///
/// * `batch`: raw bytes containing m8 metric lines
///
/// # Returns
///
/// Vec<Vec<u8>>: a vector of parsed metric line bytes
pub fn parse_m8_metrics_batch(batch: Vec<u8>) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && line[0] != b'#')
        .flat_map(|line: &[u8]| vec![line.to_vec()])
        .collect()
}

/// Parses a batch of m8 accession mappings from raw bytes.
///
/// # Arguments
///
/// * `batch`: raw bytes containing m8 records
///
/// # Returns
///
/// Vec<Vec<u8>>: a vector of (read_id, accession) line bytes
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

#[allow(dead_code)]
/// Extracts the read ID and taxid from a hit summary line.
///
/// # Arguments
///
/// * `line`: the hit summary line string
///
/// # Returns
///
/// Option<(String, i32)>: the extracted read ID and taxid if found
fn parse_read_taxid_from_hitsummary(line: &str) -> Option<(String, i32)> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() >= 3 {
        Some((parts[0].to_string(), parts[2].parse().ok()?))
    } else {
        None
    }
}


/// Generates taxonomic count JSON records from an m8 and hit summary stream.
///
/// # Arguments
///
/// * `m8_stream_rx`: receiver for m8 record stream
/// * `hit_summary_stream_rx`: receiver for hit summary record stream
/// * `db_type`: database type (e.g., "NT")
/// * `lineage_map`: map of taxids to taxonomic lineages
/// * `should_keep_filter`: filter for lineages to keep
/// * `duplicate_clusters`: map of cluster information for duplicate weighting
/// * `output_tx`: sender for the generated JSON count records
/// * `concurrency`: number of concurrent workers
/// * `batch_size_lines`: number of lines per processing batch
///
/// # Returns
///
/// Result<()>: success or error
pub async fn generate_taxon_count_json_from_m8(
    m8_stream_rx:           ReceiverStream<ParseOutput>,
    hit_summary_stream_rx:  ReceiverStream<ParseOutput>,
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



/// Generates taxonomic counts from a called m8 stream.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `m8_stream`: line-based m8 format input stream
/// * `summary_stream`: stream of summary information from m8
/// * `duplicate_clusters`: map of cluster information
/// * `should_keep_filter`: filter for taxids to keep
/// * `count_type`: type of database being referenced (e.g., "NT")
/// * `source_count_type`: optional second source type (e.g., "NR")
///
/// # Returns
///
/// Result containing a vector of taxon counts
pub async fn generate_taxon_counts(
    config: Arc<RunConfig>,
    m8_stream: ReceiverStream<ParseOutput>,
    summary_stream: ReceiverStream<ParseOutput>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    count_type: String,
    source_count_type: Option<String>, // kept for exact call-site compatibility
) -> Result<Vec<TaxonCount>, PipelineError> {
    use std::time::Duration;

    let num_workers = compute_phase_concurrency(
        &config,
        "generate_taxon_counts",
        0.6,   // ~600 MB transient per worker (lineage cache + agg map)
        2.5,
        128,   // EPYC 84c sweet spot
        8,     // still useful on MacBook Air
    );

    info!(
        "generate_taxon_counts({}) — {} workers, per-worker local AHashMap + lineage cache (exact Python semantics)",
        count_type, num_workers
    );

    let (worker_txs, worker_rxs): (Vec<_>, Vec<_>) = (0..num_workers)
        .map(|_| mpsc::channel::<(Vec<Bytes>, Vec<Bytes>)>(64))
        .unzip();

    // Spawn workers — each owns a completely private aggregation map + cache.
    // We also time each recv() so we can see whether workers are starved or blocked.
    let worker_handles: Vec<_> = worker_rxs
        .into_iter()
        .enumerate()
        .map(|(worker_idx, mut rx)| {
            let dup = duplicate_clusters.clone();
            let keep = should_keep_filter.clone();
            let src = source_count_type.clone();
            let ct = count_type.clone();

            tokio::spawn(async move {
                let mut agg: AHashMap<Vec<i32>, AggBucket> = AHashMap::new();
                let mut lineage_cache: AHashMap<i32, Vec<i32>> = AHashMap::new();

                let mut batches = 0usize;
                let mut total_pairs = 0usize;
                let mut max_recv_wait = Duration::from_secs(0);
                let mut max_batch_process = Duration::from_secs(0);
                let start = Instant::now();

                loop {
                    let recv_start = Instant::now();
                    let maybe_batch = rx.recv().await;
                    let recv_wait = recv_start.elapsed();
                    if recv_wait > max_recv_wait {
                        max_recv_wait = recv_wait;
                    }

                    let Some((m8_batch, hit_batch)) = maybe_batch else {
                        break;
                    };

                    batches += 1;
                    let batch_start = Instant::now();

                    for (m8_bytes, hit_bytes) in m8_batch.into_iter().zip(hit_batch) {
                        total_pairs += 1;
                        process_record_pair(
                            m8_bytes.as_ref(),
                            hit_bytes.as_ref(),
                            &mut agg,
                            &mut lineage_cache,
                            &dup,
                            &*keep,
                            &ct,
                            src.as_deref(),
                        )?;
                    }

                    let batch_elapsed = batch_start.elapsed();
                    if batch_elapsed > max_batch_process {
                        max_batch_process = batch_elapsed;
                    }

                    // if recv_wait >= Duration::from_millis(100)
                    //     || batch_elapsed >= Duration::from_millis(100)
                    // {
                    //     debug!(
                    //         "[generate_taxon_counts:{}] worker {} recv_wait={:?} batch_time={:?} batches={} pairs={} agg={} cache={}",
                    //         ct,
                    //         worker_idx,
                    //         recv_wait,
                    //         batch_elapsed,
                    //         batches,
                    //         total_pairs,
                    //         agg.len(),
                    //         lineage_cache.len(),
                    //     );
                    // }

                    // if last_log.elapsed() >= Duration::from_secs(10) {
                    //     info!(
                    //         "[generate_taxon_counts:{}] worker {} progress: batches={} pairs={} max_recv_wait={:?} max_batch_time={:?} agg={} cache={} elapsed={:?}",
                    //         ct,
                    //         worker_idx,
                    //         batches,
                    //         total_pairs,
                    //         max_recv_wait,
                    //         max_batch_process,
                    //         agg.len(),
                    //         lineage_cache.len(),
                    //         start.elapsed(),
                    //     );
                    //     last_log = Instant::now();
                    // }
                }

                info!(
                    "[generate_taxon_counts:{}] worker {} done: batches={} pairs={} max_recv_wait={:?} max_batch_time={:?} elapsed={:?}",
                    ct,
                    worker_idx,
                    batches,
                    total_pairs,
                    max_recv_wait,
                    max_batch_process,
                    start.elapsed(),
                );

                Ok::<AHashMap<Vec<i32>, AggBucket>, anyhow::Error>(agg)
            })
        })
        .collect();

    // Producer: pair the two streams and round-robin to workers.
    // Keep the original lockstep semantics, but add visibility at each await boundary.
    let mut m8_batch: Vec<Bytes> = Vec::with_capacity(8192);
    let mut hit_batch: Vec<Bytes> = Vec::with_capacity(8192);
    let mut batch_idx = 0usize;
    let mut total_pairs = 0usize;
    let mut m8_items_seen = 0usize;
    let mut hit_items_seen = 0usize;

    let mut m8_iter = m8_stream;
    let mut hit_iter = summary_stream;

    let mut max_m8_wait = Duration::from_secs(0);
    let mut max_hit_wait = Duration::from_secs(0);
    let mut max_send_wait = Duration::from_secs(0);
    let mut last_log = Instant::now();
    let producer_start = Instant::now();

    let mut m8_closed = false;
    let mut hit_closed = false;

    loop {
        let (m8_res, hit_res) = tokio::join!(
            async {
                let start = Instant::now();
                let item = m8_iter.next().await;
                (item, start.elapsed())
            },
            async {
                let start = Instant::now();
                let item = hit_iter.next().await;
                (item, start.elapsed())
            }
        );

        let (m8_item, m8_wait) = m8_res;
        let (hit_item, hit_wait) = hit_res;

        if m8_wait > max_m8_wait {
            max_m8_wait = m8_wait;
        }
        if hit_wait > max_hit_wait {
            max_hit_wait = hit_wait;
        }

        if m8_wait >= Duration::from_millis(50) || hit_wait >= Duration::from_millis(50) {
            let ratio = if hit_wait.as_nanos() > 0 {
                m8_wait.as_secs_f64() / hit_wait.as_secs_f64()
            } else {
                0.0
            };

            let classification = if m8_wait > hit_wait * 5 {
                "STARVING_ON_M8"
            } else if hit_wait > m8_wait * 5 {
                "STARVING_ON_SUMMARY"
            } else {
                "BOTH_SLOW"
            };

            warn!(
                "[generate_taxon_counts:{}] WAIT: m8={:?} hit={:?} ratio={:.2} class={} total_pairs={} batches={}",
                count_type,
                m8_wait,
                hit_wait,
                ratio,
                classification,
                total_pairs,
                batch_idx,
            );
        }

        let m8_present = m8_item.is_some();
        let hit_present = hit_item.is_some();

        let (Some(m8_item), Some(hit_item)) = (m8_item, hit_item) else {
            if !m8_present && !m8_closed {
                m8_closed = true;
                warn!(
                    "[generate_taxon_counts:{}] m8 stream closed at pairs={} batches={}",
                    count_type, total_pairs, batch_idx
                );
            }
            if !hit_present && !hit_closed {
                hit_closed = true;
                warn!(
                    "[generate_taxon_counts:{}] summary stream closed at pairs={} batches={}",
                    count_type, total_pairs, batch_idx
                );
            }

            warn!(
                "[generate_taxon_counts:{}] STREAM CLOSED: m8_present={} hit_present={} total_pairs={} batches={}",
                count_type,
                m8_present,
                hit_present,
                total_pairs,
                batch_idx
            );
            break;
        };

        if let ParseOutput::Bytes(b) = m8_item {
            if !b.is_empty() {
                m8_items_seen += 1;
                m8_batch.push(b);
                total_pairs += 1;
            }
        }

        if let ParseOutput::Bytes(b) = hit_item {
            if !b.is_empty() {
                hit_items_seen += 1;
                hit_batch.push(b);
            }
        }

        if m8_batch.len() >= 8192 && hit_batch.len() >= 8192 {
            let worker_idx = batch_idx % num_workers;

            let send_start = Instant::now();
            worker_txs[worker_idx]
                .send((std::mem::take(&mut m8_batch), std::mem::take(&mut hit_batch)))
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "generate_taxon_counts({}) worker send failed: {}",
                        count_type,
                        e
                    ))
                })?;
            let send_wait = send_start.elapsed();

            if send_wait > max_send_wait {
                max_send_wait = send_wait;
            }

            if m8_wait >= Duration::from_millis(100)
                || hit_wait >= Duration::from_millis(100)
                || send_wait >= Duration::from_millis(100)
            {
                info!(
                    "[generate_taxon_counts:{}] batch {} waits: m8={:?} hit={:?} send={:?} total_pairs={} worker={}",
                    count_type,
                    batch_idx,
                    m8_wait,
                    hit_wait,
                    send_wait,
                    total_pairs,
                    worker_idx,
                );
            }

            batch_idx += 1;
        }

        if last_log.elapsed() >= Duration::from_secs(10) {
            info!(
                "[generate_taxon_counts:{}] producer progress: pairs={} batches={} max_waits(m8={:?}, hit={:?}, send={:?}) elapsed={:?}",
                count_type,
                total_pairs,
                batch_idx,
                max_m8_wait,
                max_hit_wait,
                max_send_wait,
                producer_start.elapsed(),
            );
            last_log = Instant::now();
        }
    }

    // Final partial batch
    if !m8_batch.is_empty() || !hit_batch.is_empty() {
        let worker_idx = batch_idx % num_workers;
        let send_start = Instant::now();
        worker_txs[worker_idx]
            .send((m8_batch, hit_batch))
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "generate_taxon_counts({}) final worker send failed: {}",
                    count_type,
                    e
                ))
            })?;
        let send_wait = send_start.elapsed();
        if send_wait > max_send_wait {
            max_send_wait = send_wait;
        }

        info!(
            "[generate_taxon_counts:{}] final batch sent to worker {} in {:?}",
            count_type, worker_idx, send_wait
        );
    }

    drop(worker_txs); // close all workers

    // Merge the tiny partial maps from each worker.
    let partials: Vec<AHashMap<Vec<i32>, AggBucket>> = try_join_all(worker_handles)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;

    let global_agg = merge_aggregations(partials);

    // Final Python Loop 2 — identical output shape
    let result = build_taxon_counts_list(global_agg, &count_type);

    info!(
        "generate_taxon_counts({}) complete — {} taxa (max waits: m8={:?}, hit={:?}, send={:?})",
        count_type,
        result.len(),
        max_m8_wait,
        max_hit_wait,
        max_send_wait,
    );

    Ok(result)
}


/// Computes merged taxonomic counts from both NT and NR databases.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `nt_m8_stream`: receiver for NT m8 records
/// * `nt_hit_summary_stream`: receiver for NT hit summary records
/// * `nt_contig_summary`: summary entries for NT contigs
/// * `nr_m8_stream`: receiver for NR m8 records
/// * `nr_hit_summary_stream`: receiver for NR hit summary records
/// * `nr_contig_summary`: summary entries for NR contigs
/// * `_lineage_map`: map of taxids to lineages (unused)
/// * `should_keep_filter`: filter for taxids to keep
/// * `duplicate_clusters`: map of cluster information
/// * `merged_m8_path`: output path for merged m8 file
/// * `merged_hitsummary_path`: output path for merged hit summary file
/// * `merged_taxon_counts_path`: output path for merged taxon counts file
/// * `merged_contig_summary_path`: output path for merged contig summary file
/// * `_nr_alignment_per_read`: map of NR alignments per read (unused)
///
/// # Returns
///
/// Result containing a vector of cleanup tasks
pub async fn compute_merged_taxon_counts(
    config: Arc<RunConfig>,
    nt_m8_stream: mpsc::Receiver<ParseOutput>,
    nt_hit_summary_stream: mpsc::Receiver<ParseOutput>,
    nt_contig_summary: Vec<ContigSummaryEntry>,

    nr_m8_stream: mpsc::Receiver<ParseOutput>,
    nr_hit_summary_stream: mpsc::Receiver<ParseOutput>,
    nr_contig_summary: Vec<ContigSummaryEntry>,

    _lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,

    merged_m8_path: PathBuf,
    merged_hitsummary_path: PathBuf,
    merged_taxon_counts_path: PathBuf,
    merged_contig_summary_path: PathBuf,
    _nr_alignment_per_read: Arc<DashMap<String, SpeciesAlignmentResults, AHashRandomState>>,
) -> Result<Vec<JoinHandle<Result<(), anyhow::Error>>>, PipelineError> {
    let mut cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();

    // 1. Parallel NR preload
    let (nr_alignment_map, nr_hit_lines) = tokio::task::spawn_blocking(move || {
        let mut alignment_map: AHashMap<String, SpeciesAlignmentResults> = AHashMap::new();
        let mut hit_lines: AHashMap<String, Bytes> = AHashMap::new();

        let mut nr_hit_stream = ReceiverStream::new(nr_hit_summary_stream);

        while let Some(item) = futures::executor::block_on(nr_hit_stream.next()) {
            let bytes = match item.to_bytes() {
                Ok(b) => b,
                Err(_) => continue,
            };
            let line = String::from_utf8_lossy(&bytes);
            let trimmed = line.trim_end();
            if trimmed.is_empty() { continue; }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 10 { continue; }

            let read_id = fields[0].to_string();
            let contig_species = fields.get(9).and_then(|s| s.parse::<Taxid>().ok());
            let read_species   = fields.get(3).and_then(|s| s.parse::<Taxid>().ok());

            alignment_map.insert(read_id.clone(), SpeciesAlignmentResults {
                contig: contig_species,
                read: read_species,
            });
            hit_lines.insert(read_id, bytes);
        }
        (alignment_map, hit_lines)
    }).await.map_err(|e| PipelineError::Other(e.into()))?;

    // 2. Parallel file writers (background cleanup)
    let mut merged_m8_file = BufWriter::new(TokioFile::create(&merged_m8_path).await?);
    let mut merged_hit_file = BufWriter::new(TokioFile::create(&merged_hitsummary_path).await?);

    let (m8_write_tx, mut m8_write_rx) = mpsc::channel::<Bytes>(4096);
    let (hit_write_tx, mut hit_write_rx) = mpsc::channel::<Bytes>(4096);

    cleanup_tasks.push(tokio::spawn(async move {
        while let Some(bytes) = m8_write_rx.recv().await {
            let _ = merged_m8_file.write_all(&bytes).await;
            let _ = merged_m8_file.write_all(b"\n").await;
        }
        let _ = merged_m8_file.flush().await;
        Ok::<(), anyhow::Error>(())
    }));

    cleanup_tasks.push(tokio::spawn(async move {
        while let Some(bytes) = hit_write_rx.recv().await {
            let _ = merged_hit_file.write_all(&bytes).await;
            let _ = merged_hit_file.write_all(b"\n").await;
        }
        let _ = merged_hit_file.flush().await;
        Ok::<(), anyhow::Error>(())
    }));

    // 3. Parallel selection passes (NT first, then NR)
    let (merged_m8_tx, merged_m8_rx) = mpsc::channel::<ParseOutput>(4096);
    let (merged_hit_tx, merged_hit_rx) = mpsc::channel::<ParseOutput>(4096);

    // NT pass
    let nt_m8_write_tx = m8_write_tx.clone();
    let nt_hit_write_tx = hit_write_tx.clone();
    let merged_m8_tx_nt = merged_m8_tx.clone();
    let merged_hit_tx_nt = merged_hit_tx.clone();

    let nt_task = tokio::spawn({
        let nr_alignment_map = nr_alignment_map.clone();
        async move {
            let mut nt_m8 = ReceiverStream::new(nt_m8_stream);
            let mut nt_hit = ReceiverStream::new(nt_hit_summary_stream);

            while let (Some(m8_item), Some(hit_item)) = (nt_m8.next().await, nt_hit.next().await) {
                let m8_bytes = match m8_item { ParseOutput::Bytes(b) => b, _ => continue };
                let hit_bytes = hit_item.to_bytes()?;

                let hit_line = String::from_utf8_lossy(&hit_bytes);
                let trimmed = hit_line.trim_end();
                if trimmed.is_empty() { continue; }

                let hit_fields: Vec<&str> = trimmed.split('\t').collect();
                if hit_fields.len() < 10 { continue; }

                let read_id = hit_fields[0];

                let nt_contig = hit_fields.get(9).and_then(|s| s.parse::<Taxid>().ok());
                let nt_read   = hit_fields.get(3).and_then(|s| s.parse::<Taxid>().ok());

                let nr_align = nr_alignment_map.get(read_id);
                let has_nr_contig = nr_align.map_or(false, |a| a.contig.is_some());
                let _has_nr_read   = nr_align.map_or(false, |a| a.read.is_some());

                if nt_contig.is_some() || (!has_nr_contig && nt_read.is_some()) {
                    let _ = nt_m8_write_tx.send(m8_bytes.clone()).await;
                    let mut hit_with_source = hit_fields.to_vec();
                    hit_with_source.push(NT_TAG);
                    let hit_bytes_out = Bytes::from(hit_with_source.join("\t"));
                    let _ = nt_hit_write_tx.send(hit_bytes_out.clone()).await;

                    let _ = merged_m8_tx_nt.send(ParseOutput::Bytes(m8_bytes)).await;
                    let _ = merged_hit_tx_nt.send(ParseOutput::Bytes(hit_bytes_out)).await;
                }
            }
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(nt_task);

    // NR pass (remaining reads)
    let nr_m8_write_tx = m8_write_tx;
    let nr_hit_write_tx = hit_write_tx;
    let merged_m8_tx_nr = merged_m8_tx;
    let merged_hit_tx_nr = merged_hit_tx;

    let nr_task = tokio::spawn({
        let nr_hit_lines = nr_hit_lines.clone();
        async move {
            let mut nr_m8 = ReceiverStream::new(nr_m8_stream);

            while let Some(m8_item) = nr_m8.next().await {
                let m8_bytes = match m8_item { ParseOutput::Bytes(b) => b, _ => continue };

                let m8_str = String::from_utf8_lossy(&m8_bytes);
                let trimmed_m8 = m8_str.trim_end();
                if trimmed_m8.is_empty() { continue; }

                let m8_fields: Vec<&str> = trimmed_m8.split('\t').collect();
                let read_id = m8_fields[0];

                if let Some(hit_bytes) = nr_hit_lines.get(read_id) {
                    let _ = nr_m8_write_tx.send(m8_bytes.clone()).await;

                    let mut hit_with_source = String::from_utf8_lossy(hit_bytes)
                        .trim_end()
                        .split('\t')
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>();
                    hit_with_source.push(NR_TAG.to_string());

                    let hit_bytes_out = Bytes::from(hit_with_source.join("\t"));
                    let _ = nr_hit_write_tx.send(hit_bytes_out.clone()).await;

                    let _ = merged_m8_tx_nr.send(ParseOutput::Bytes(m8_bytes)).await;
                    let _ = merged_hit_tx_nr.send(ParseOutput::Bytes(hit_bytes_out)).await;
                }
            }
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(nr_task);

    // 4. Taxon counting (already highly parallel from earlier work)
    let merged_taxon_counts = generate_taxon_counts(
        config.clone(),
        ReceiverStream::new(merged_m8_rx),
        ReceiverStream::new(merged_hit_rx),
        duplicate_clusters.clone(),
        should_keep_filter.clone(),
        "merged_NT_NR".to_string(),
        None,
    ).await?;

    let json = serde_json::to_string_pretty(&merged_taxon_counts)?;
    tokio::fs::write(&merged_taxon_counts_path, json).await?;
    info!("Merged taxon counts generated (fully streaming)");

    // 5. Parallel contig merge (exact Python _merge_contigs)
    let final_contigs = tokio::task::spawn_blocking(move || {
        let mut merged_contigs: HashMap<Taxid, Vec<ContigSummaryEntry>> = HashMap::new();
        let mut nt_contig_names: HashSet<String> = HashSet::new();

        for mut entry in nt_contig_summary {
            entry.db_type = "merged_NT_NR".to_string();
            nt_contig_names.insert(entry.contig_name.clone());
            merged_contigs.entry(entry.species_taxid).or_default().push(entry);
        }

        for mut entry in nr_contig_summary {
            if !nt_contig_names.contains(&entry.contig_name) {
                entry.db_type = "merged_NT_NR".to_string();
                merged_contigs.entry(entry.species_taxid).or_default().push(entry);
            }
        }

        merged_contigs.into_values().flatten().collect::<Vec<_>>()
    }).await.map_err(|e| PipelineError::Other(e.into()))?;

    let json = serde_json::to_string_pretty(&final_contigs)?;
    tokio::fs::write(&merged_contig_summary_path, json).await?;
    info!("Merged contig summary written ({} entries)", final_contigs.len());

    info!("compute_merged_taxon_counts complete — exact Python logic, fully streaming");
    Ok(cleanup_tasks)
}



#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // ── Strategy for generating tag fields (valid + malformed) ─────────────
    fn tag_field_strategy() -> impl Strategy<Value = String> {
        prop_oneof![
            // Valid tag: key:type:value
            ("[A-Za-z0-9_]{1,10}", "[A-Za-z]:", "[^\\t\\n]{0,40}")
                .prop_map(|(k, t, v)| format!("{}:{}{}", k, t, v)),

            // Malformed: missing type (only one colon)
            ("[A-Za-z0-9_]{1,10}", "[^:\\t\\n]{0,30}")
                .prop_map(|(k, v)| format!("{}:{}", k, v)),

            // Malformed: empty key
            (":[A-Za-z]:[^\\t\\n]{0,30}").prop_map(|s| s.to_string()),

            // Malformed: empty value
            ("[A-Za-z0-9_]{1,10}:[A-Za-z]:").prop_map(|s| s.to_string()),

            // Malformed: garbage / no colons
            "[A-Za-z0-9_:\\-]{1,40}".prop_map(|s| s),

            // Malformed: trailing colon only
            "[A-Za-z0-9_]{1,10}:".prop_map(|s| s.to_string()),

            // Empty tag field (consecutive tabs)
            Just("".to_string()),
        ]
    }

    // ── Proptests with random extra tags (valid + malformed) ───────────────
    proptest! {
        #[test]
        fn test_parse_line_nt_equivalence_with_malformed_tags(
            qname in "[a-zA-Z0-9._-]{1,100}",
            tname in "[a-zA-Z0-9._-]{1,50}",
            pident in 0.0..100.0,
            alen in 0u64..10000,
            mismatch in 0u64..1000,
            gapopen in 0u64..100,
            qstart in 1u64..10000,
            qend in 1u64..10000,
            tstart in 1u64..1000000,
            tend in 1u64..1000000,
            evalue in 0.0..1.0,
            bitscore in 0.0..1000.0,
            qlen in 1u64..10000,
            slen in 1u64..1000000,
            extra_tags in prop::collection::vec(tag_field_strategy(), 0..35)
        ) {
            let mut line = format!(
                "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{:.3}\t{}\t{}",
                qname, tname, pident, alen, mismatch, gapopen, qstart, qend,
                tstart, tend, evalue, bitscore, qlen, slen
            );

            for tag in extra_tags {
                line.push_str(&format!("\t{}", tag));
            }

            let scalar_res = M8Record::parse_line_nt_scalar(&line);
            let dispatched_res = M8Record::parse_line_nt(&line);

            match (scalar_res, dispatched_res) {
                (Ok(s), Ok(d)) => assert_m8_eq(&s, &d, &line),
                (Err(_), Err(_)) => {}
                _ => panic!("Scalar and AVX-512 disagreed on success/failure for line: {}", line),
            }
        }

        #[test]
        fn test_parse_line_nr_equivalence_with_malformed_tags(
            qname in "[a-zA-Z0-9._-]{1,100}",
            tname in "[a-zA-Z0-9._-]{1,50}",
            pident in 0.0..100.0,
            alen in 0u64..10000,
            mismatch in 0u64..1000,
            gapopen in 0u64..100,
            qstart in 1u64..10000,
            qend in 1u64..10000,
            tstart in 1u64..1000000,
            tend in 1u64..1000000,
            evalue in 0.0..1.0,
            bitscore in 0.0..1000.0,
            extra_tags in prop::collection::vec(tag_field_strategy(), 0..35)
        ) {
            let mut line = format!(
                "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{:.3}",
                qname, tname, pident, alen, mismatch, gapopen, qstart, qend,
                tstart, tend, evalue, bitscore
            );

            for tag in extra_tags {
                line.push_str(&format!("\t{}", tag));
            }

            let scalar_res = M8Record::parse_line_nr_scalar(&line);
            let dispatched_res = M8Record::parse_line_nr(&line);

            match (scalar_res, dispatched_res) {
                (Ok(s), Ok(d)) => assert_m8_eq(&s, &d, &line),
                (Err(_), Err(_)) => {}
                _ => panic!("Scalar and AVX-512 disagreed on success/failure for line: {}", line),
            }
        }
    }

    // ── helpers ────────────────────────────────────────────────────────────

    #[allow(dead_code)]
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
        let dispatched = M8Record::parse_line_nt(line).expect("NT dispatch failed");
        assert_m8_eq(&scalar, &dispatched, line);
    }

    fn compare_nr(line: &str) {
        let scalar = M8Record::parse_line_nr_scalar(line).expect("NR scalar failed");
        let dispatched = M8Record::parse_line_nr(line).expect("NR dispatch failed");
        assert_m8_eq(&scalar, &dispatched, line);
    }

    // ── NT test constants ──────────────────────────────────────────────────
    const NT_BASIC: &str =
        "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903";

    const NT_ZERO_EVALUE: &str =
        "read2\tNC_000001.1\t100.000\t250\t0\t0\t1\t250\t1\t250\t0.0\t462.000\t250\t248956422";

    const NT_MISMATCHES: &str =
        "read3\tKJ660346.2\t95.238\t210\t10\t2\t5\t214\t500\t709\t4.56e-88\t330.000\t220\t18958";

    const NT_LONG_QNAME: &str =
        "a_very_long_read_identifier_that_exceeds_sixty_four_bytes_total\tNC_045512.2\t98.667\t150\t2\t0\t1\t150\t1\t150\t2.34e-70\t275.000\t150\t29903";

    const NT_TRAILING_WHITESPACE: &str =
        "read4\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903   ";

    const NT_EXTRA_COLUMNS: &str =
        "read5\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903\textra1\textra2";

    // ── NR test constants ──────────────────────────────────────────────────
    const NR_BASIC: &str =
        "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000";

    const NR_ZERO_EVALUE: &str =
        "read2\tAAA12345.1\t100.000\t300\t0\t0\t1\t300\t1\t300\t0.0\t600.000";

    const NR_GAPS: &str =
        "read3\tXYZ99999.2\t78.431\t153\t33\t3\t10\t162\t5\t155\t8.9e-20\t89.700";

    const NR_TRAILING_WHITESPACE: &str =
        "read4\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000   ";

    const NR_EXTRA_COLUMNS: &str =
        "read5\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000\textra1\textra2";

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
    fn test_nt_extra_columns() {
        compare_nt(NT_EXTRA_COLUMNS);
        let r = M8Record::parse_line_nt(NT_EXTRA_COLUMNS).unwrap();
        assert_eq!(r.tname, "NC_045512");
    }

    #[test]
    fn test_nt_accession_version_stripped() {
        let r = M8Record::parse_line_nt_scalar(NT_BASIC).unwrap();
        assert_eq!(r.tname, "NC_045512", "version suffix must be stripped");
    }

    #[test]
    fn test_nt_error_on_empty() {
        assert!(M8Record::parse_line_nt_scalar("").is_err());
        assert!(M8Record::parse_line_nt("").is_err());
    }

    #[test]
    fn test_nt_error_on_too_few_fields() {
        let short = "read1\tNC_045512.2\t99.333\t150";
        assert!(M8Record::parse_line_nt_scalar(short).is_err());
        assert!(M8Record::parse_line_nt(short).is_err());
    }

    // ── NR tests ──────────────────────────────────────────────────────────
    #[test]
    fn test_nr_basic() { compare_nr(NR_BASIC); }

    #[test]
    fn test_nr_zero_evalue() { compare_nr(NR_ZERO_EVALUE); }

    #[test]
    fn test_nr_gaps() { compare_nr(NR_GAPS); }

    #[test]
    fn test_nr_trailing_whitespace() {
        compare_nr(NR_TRAILING_WHITESPACE);
        let r = M8Record::parse_line_nr_scalar(NR_TRAILING_WHITESPACE).unwrap();
        assert_eq!(r.qname, "read4");
    }

    #[test]
    fn test_nr_extra_columns() {
        compare_nr(NR_EXTRA_COLUMNS);
        let r = M8Record::parse_line_nr(NR_EXTRA_COLUMNS).unwrap();
        assert_eq!(r.tname, "QIK02963");
    }

    #[test]
    fn test_nr_qlen_slen_zero() {
        let r = M8Record::parse_line_nr_scalar(NR_BASIC).unwrap();
        assert_eq!(r.qlen, 0);
        assert_eq!(r.slen, 0);
    }

    #[test]
    fn test_nr_accession_version_stripped() {
        let r = M8Record::parse_line_nr_scalar(NR_BASIC).unwrap();
        assert_eq!(r.tname, "QIK02963");
    }

    #[test]
    fn test_nr_error_on_empty() {
        assert!(M8Record::parse_line_nr_scalar("").is_err());
        assert!(M8Record::parse_line_nr("").is_err());
    }

    #[test]
    fn test_nr_error_on_too_few_fields() {
        let short = "read1\tQIK02963.1\t87.500";
        assert!(M8Record::parse_line_nr_scalar(short).is_err());
        assert!(M8Record::parse_line_nr(short).is_err());
    }

    // ── Additional strong edge-case tests ──────────────────────────────────

    #[test]
    fn test_nt_malformed_and_extra_columns() {
        let cases = vec![
            "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000",
            "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903\textra1\textra2\textra3",
            "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903\t",
            "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t\t150\t29903",
        ];

        for line in cases {
            let scalar_res = M8Record::parse_line_nt_scalar(line);
            let dispatched_res = M8Record::parse_line_nt(line);

            match (scalar_res, dispatched_res) {
                (Ok(s), Ok(d)) => assert_m8_eq(&s, &d, line),
                (Err(_), Err(_)) => {}
                _ => panic!("Scalar vs dispatched disagreement on: {}", line),
            }
        }
    }

    #[test]
    fn test_nr_malformed_and_extra_columns() {
        let cases = vec![
            "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38",
            "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000\textra1\textra2",
            "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000\t",
            "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t\t96\t1.45e-38\t152.000",
        ];

        for line in cases {
            let scalar_res = M8Record::parse_line_nr_scalar(line);
            let dispatched_res = M8Record::parse_line_nr(line);

            match (scalar_res, dispatched_res) {
                (Ok(s), Ok(d)) => assert_m8_eq(&s, &d, line),
                (Err(_), Err(_)) => {}
                _ => panic!("Scalar vs dispatched disagreement on: {}", line),
            }
        }
    }

    #[test]
    fn test_nt_long_line_many_tags() {
        let mut line = NT_BASIC.to_string();
        for i in 0..60 {
            line.push_str(&format!("\ttag{}:Z:value{}", i, i));
        }
        compare_nt(&line);
    }

    #[test]
    fn test_nr_long_line_many_tags() {
        let mut line = NR_BASIC.to_string();
        for i in 0..60 {
            line.push_str(&format!("\ttag{}:Z:value{}", i, i));
        }
        compare_nr(&line);
    }

    // ── Existing consensus / aggregation / summarize tests ────────────────

    use fst::MapBuilder;

    fn lineage(species: i32, genus: i32, family: i32) -> Lineage {
        [species, genus, family]
    }

    fn build_lineage_map() -> AHashMap<Taxid, Lineage> {
        let mut map = AHashMap::default();
        map.insert(1 as Taxid, lineage(1, 10, 100));
        map.insert(2 as Taxid, lineage(2, 20, 200));
        map.insert(3 as Taxid, lineage(3, 30, 300));
        map.insert(4 as Taxid, lineage(1, 10, 999));
        map
    }

    fn build_acc_map(entries: &[(&str, u64)]) -> Map<Vec<u8>> {
        let mut builder = MapBuilder::memory();
        for (k, v) in entries {
            builder.insert(k, *v).unwrap();
        }
        builder.into_map()
    }

    fn nr_record(qname: &str, tname: &str, alen: u64, bitscore: f64) -> M8Record {
        M8Record {
            qname: qname.to_string(),
            tname: tname.to_string(),
            pident: 99.0,
            alen,
            mismatch: 0,
            gapopen: 0,
            qstart: 1,
            qend: alen,
            tstart: 1,
            tend: alen,
            evalue: 1e-40,
            bitscore,
            qlen: 0,
            slen: 0,
        }
    }

    #[test]
    fn test_consensus_level_single_hit_returns_family_level() {
        let hits = vec![nr_record("read1", "ACC1", 100, 200.0)];
        let lineage_map = build_lineage_map();
        let acc2taxid_map = build_acc_map(&[("ACC1", 1)]);
        let should_keep = Arc::new(|_: &[i32]| true);

        let (level, taxid, kept) = consensus_level(&hits, &lineage_map, &acc2taxid_map, &should_keep)
            .expect("consensus_level should succeed");

        assert_eq!(level, 3);
        assert_eq!(taxid, 100);
        assert_eq!(kept.len(), 1);
    }

    #[test]
    fn test_consensus_level_species_and_genus_agree() {
        let mut lineage_map = AHashMap::default();
        lineage_map.insert(1 as Taxid, lineage(1, 10, 100));
        lineage_map.insert(2 as Taxid, lineage(1, 10, 200));

        let hits = vec![
            nr_record("read1", "ACC1", 100, 200.0),
            nr_record("read1", "ACC2", 100, 199.0),
        ];
        let acc2taxid_map = build_acc_map(&[("ACC1", 1), ("ACC2", 2)]);
        let should_keep = Arc::new(|_: &[i32]| true);

        let (level, taxid, kept) = consensus_level(&hits, &lineage_map, &acc2taxid_map, &should_keep)
            .expect("consensus_level should succeed");

        assert_eq!(level, 2);
        assert_eq!(taxid, 10);
        assert_eq!(kept.len(), 2);
    }

    #[test]
    fn test_consensus_level_rejected_by_filter() {
        let hits = vec![nr_record("read1", "ACC1", 100, 200.0)];
        let lineage_map = build_lineage_map();
        let acc2taxid_map = build_acc_map(&[("ACC1", 1)]);
        let should_keep = Arc::new(|_: &[i32]| false);

        let (level, taxid, kept) = consensus_level(&hits, &lineage_map, &acc2taxid_map, &should_keep)
            .expect("consensus_level should return ok");

        assert_eq!((level, taxid), (0, 0));
        assert!(kept.is_empty());
    }

    #[test]
    fn test_merge_aggregations_is_additive() {
        let mut part1: AHashMap<Vec<i32>, AggBucket> = AHashMap::default();
        let mut part2: AHashMap<Vec<i32>, AggBucket> = AHashMap::default();

        let mut b1 = AggBucket::default();
        b1.nonunique_count = 2;
        b1.unique_count = 1;
        b1.base_count = 100;
        b1.sum_percent_identity = 99.0;
        b1.sum_alignment_length = 100.0;
        b1.sum_e_value = -40.0;
        b1.source_count_type.insert("NT".to_string());

        let mut b2 = AggBucket::default();
        b2.nonunique_count = 3;
        b2.unique_count = 4;
        b2.base_count = 250;
        b2.sum_percent_identity = 198.0;
        b2.sum_alignment_length = 250.0;
        b2.sum_e_value = -80.0;
        b2.source_count_type.insert("NR".to_string());

        part1.insert(vec![1, 10, 100], b1);
        part2.insert(vec![1, 10, 100], b2);

        let merged = merge_aggregations(vec![part1, part2]);
        let bucket = merged.get(&vec![1, 10, 100]).unwrap();

        assert_eq!(bucket.nonunique_count, 5);
        assert_eq!(bucket.unique_count, 5);
        assert_eq!(bucket.base_count, 350);
        assert_eq!(bucket.sum_percent_identity, 297.0);
        assert_eq!(bucket.sum_alignment_length, 350.0);
        assert_eq!(bucket.sum_e_value, -120.0);
        assert!(bucket.source_count_type.contains("NT"));
        assert!(bucket.source_count_type.contains("NR"));
    }

    #[test]
    fn test_build_taxon_counts_list_preserves_taxon_fields_and_sources() {
        let mut agg: AHashMap<Vec<i32>, AggBucket> = AHashMap::default();

        let mut bucket = AggBucket::default();
        bucket.nonunique_count = 5;
        bucket.unique_count = 2;
        bucket.base_count = 250;
        bucket.sum_percent_identity = 198.0;
        bucket.sum_alignment_length = 210.0;
        bucket.sum_e_value = -80.0;
        bucket.source_count_type.insert("NR".to_string());
        bucket.source_count_type.insert("NT".to_string());

        agg.insert(vec![11, 22, 33], bucket);

        let rows = build_taxon_counts_list(agg, "merged_NT_NR");
        assert_eq!(rows.len(), 1);

        let row = &rows[0];
        assert_eq!(row.tax_id, 11);
        assert_eq!(row.tax_level, 1);
        assert_eq!(row.genus_taxid, 22);
        assert_eq!(row.family_taxid, 33);
        assert_eq!(row.nonunique_count, 5);
        assert_eq!(row.unique_count, 2);
        assert_eq!(row.base_count, 250);
        assert_eq!(row.count_type, "merged_NT_NR");
        assert_eq!(
            row.source_count_type.as_ref().unwrap(),
            &vec!["NR".to_string(), "NT".to_string()]
        );
        assert_eq!(row.dcr, 2.5);
    }

    #[test]
    fn test_process_record_pair_skips_malformed_or_invalid_hits() {
        let mut agg: AHashMap<Vec<i32>, AggBucket> = AHashMap::default();
        let mut lineage_cache: AHashMap<i32, Vec<i32>> = AHashMap::default();
        let duplicate_clusters: DashMap<String, ClusterInfo> = DashMap::new();
        let should_keep = |_: &[i32]| true;

        let bad_hit = b"read1\t1\t0";
        let m8 = b"read1\tACC1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200";

        process_record_pair(
            m8,
            bad_hit,
            &mut agg,
            &mut lineage_cache,
            &duplicate_clusters,
            &should_keep,
            "NT",
            None,
        )
            .expect("process_record_pair should not fail");

        assert!(agg.is_empty());
    }

    #[test]
    fn test_process_record_pair_populates_aggregations() {
        let mut agg: AHashMap<Vec<i32>, AggBucket> = AHashMap::default();
        let mut lineage_cache: AHashMap<i32, Vec<i32>> = AHashMap::default();
        let duplicate_clusters: DashMap<String, ClusterInfo> = DashMap::new();
        let should_keep = |_: &[i32]| true;

        let hit = b"read1\t1\t1\t0\t1\t10\t100";
        let m8 = b"read1\tACC1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200";

        process_record_pair(
            m8,
            hit,
            &mut agg,
            &mut lineage_cache,
            &duplicate_clusters,
            &should_keep,
            "NT",
            None,
        )
            .expect("process_record_pair should succeed");

        assert!(!agg.is_empty());
        let bucket = agg.values().next().unwrap();
        assert_eq!(bucket.unique_count, 1);
        assert_eq!(bucket.nonunique_count, 1);
        assert_eq!(bucket.base_count, 100);
        assert!((bucket.sum_alignment_length - 100.0).abs() < 1e-9);
    }

    #[test]
    fn test_process_record_pair_multiplies_nr_alignment_length_for_merged_nt_nr() {
        let mut agg: AHashMap<Vec<i32>, AggBucket> = AHashMap::default();
        let mut lineage_cache: AHashMap<i32, Vec<i32>> = AHashMap::default();
        let duplicate_clusters: DashMap<String, ClusterInfo> = DashMap::new();
        let should_keep = |_: &[i32]| true;

        let hit = b"read1\t1\t1\t0\t1\t10\t100";
        let m8 = b"read1\tACC1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200";

        process_record_pair(
            m8,
            hit,
            &mut agg,
            &mut lineage_cache,
            &duplicate_clusters,
            &should_keep,
            "merged_NT_NR",
            Some("NR"),
        )
            .expect("process_record_pair should succeed");

        let bucket = agg.values().next().unwrap();
        assert_eq!(bucket.base_count, 100);
        assert!((bucket.sum_alignment_length - 300.0).abs() < 1e-9);
    }

    #[test]
    fn test_summarize_m8_hits_filters_short_alignments_and_picks_best_bitscore() {
        let lineage_map = build_lineage_map();
        let acc2taxid_map = build_acc_map(&[("ACC1", 1), ("ACC2", 1), ("SHORT", 1)]);
        let should_keep = |_: &[i32]| true;

        let hits = vec![
            nr_record("read1", "SHORT", 20, 500.0),
            nr_record("read1", "ACC2", 100, 150.0),
            nr_record("read1", "ACC1", 100, 150.0),
        ];

        let reduced = summarize_m8_hits(
            7,
            "read1".to_string(),
            &hits,
            &lineage_map,
            &acc2taxid_map,
            &should_keep,
            50,
        )
            .expect("expected a reduced read");

        assert_eq!(reduced.seq, 7);
        assert_eq!(reduced.accession, "ACC2");
        assert!(std::str::from_utf8(&reduced.dedup).unwrap().starts_with("read1\tACC2\t"));
        assert!(std::str::from_utf8(&reduced.summary).unwrap().starts_with("read1\tACC2\t"));
    }

    #[test]
    fn test_summarize_m8_hits_returns_none_when_all_hits_filtered_out() {
        let lineage_map = build_lineage_map();
        let acc2taxid_map = build_acc_map(&[("SHORT", 1)]);
        let should_keep = |_: &[i32]| true;

        let hits = vec![nr_record("read1", "SHORT", 20, 500.0)];

        let reduced = summarize_m8_hits(
            7,
            "read1".to_string(),
            &hits,
            &lineage_map,
            &acc2taxid_map,
            &should_keep,
            50,
        );

        assert!(reduced.is_none());
    }

    #[test]
    fn test_read_id_from_m8_line_basic() {
        assert_eq!(read_id_from_m8_line("read1\tACC1\tx"), Some("read1"));
    }

    #[test]
    fn test_read_id_from_m8_line_empty_first_field() {
        assert_eq!(read_id_from_m8_line("\tACC1\tx"), None);
    }

    #[test]
    fn test_shard_for_read_id_is_deterministic() {
        let a = shard_for_read_id("read1", 8);
        let b = shard_for_read_id("read1", 8);
        assert_eq!(a, b);
    }

    #[test]
    fn test_shard_for_read_id_handles_zero_workers() {
        let shard = shard_for_read_id("read1", 0);
        assert_eq!(shard, 0);
    }

    #[test]
    fn test_parse_summary_batch_skips_comments_and_blanks() {
        let batch = b"read1\tfoo\n#comment\n\nread2\tbar\n".to_vec();
        let lines = parse_summary_batch(batch);
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], b"read1\tfoo".to_vec());
        assert_eq!(lines[1], b"read2\tbar".to_vec());
    }

    #[test]
    fn test_parse_m8_metrics_batch_skips_comments_and_blanks() {
        let batch = b"read1\t100\n#comment\n\nread2\t200\n".to_vec();
        let lines = parse_m8_metrics_batch(batch);
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], b"read1\t100".to_vec());
        assert_eq!(lines[1], b"read2\t200".to_vec());
    }

    #[test]
    fn test_parse_m8_acc_batch_parses_read_and_accession() {
        let batch = b"read1\tACC1\textra\n#comment\nread2\tACC2   \n\tACC3\nbadline\n".to_vec();
        let lines = parse_m8_acc_batch(batch);
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], b"read1\tACC1\n".to_vec());
        assert_eq!(lines[1], b"read2\tACC2\n".to_vec());
    }
}