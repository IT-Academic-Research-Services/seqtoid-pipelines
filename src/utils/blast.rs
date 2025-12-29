// BLAST-related file functions and structures
use std::collections::{HashMap, HashSet};
use std::env::args;
use std::sync::Arc;
use std::path::PathBuf;

use anyhow::{anyhow, Result};
use log::{self, LevelFilter, debug, info, error, warn};
use tokio_stream::StreamExt;
use lexical::parse as lexical_parse;
use fst::Map;
use ahash::AHashMap;
use serde::{Deserialize, Serialize};
use tokio_stream::wrappers::ReceiverStream;
use tokio::fs::File as TokioFile;
use tokio::sync::mpsc;
use tokio::sync::mpsc::{channel, Sender};
use tokio::io::{AsyncBufReadExt, BufReader, BufWriter, AsyncWriteExt};
use tokio_stream::wrappers::BroadcastStream;
use crate::config::defs::{Taxid, Lineage, NT_TAG, NR_TAG, RunConfig, ClusterInfo};
use crate::utils::streams::ParseOutput;
use crate::utils::taxonomy::validate_taxid_lineage;
use crate::utils::streams::ToBytes;
use crate::utils::file::write_byte_stream_to_file;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigSummaryEntry {
    contig_name: String,
    common_name: Option<String>,
    category_name: Option<String>,
    score: Option<f64>,
    db_type: String,
    reads: u64,
    bases: u64,
    species_taxid: i32,
    genus_taxid: i32,
    family_taxid: i32,
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
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    mut output_tx: Sender<ParseOutput>,
) -> Result<()> {
    let mut buckets: AHashMap<Taxid, AggBucket> = AHashMap::with_capacity(500_000);


    // 0: read_id
    // 1: accession
    // 2: hit_taxid (often isolate-level)
    // 3: genus_taxid
    // 4: family_taxid
    // 5: consensus_level (1=genus, 2=family, 3=species)
    let mut read_to_taxid: AHashMap<String, Taxid> = AHashMap::with_capacity(10_000_000);
    while let Some(item) = refined_hit_summary_stream.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            warn!(
                "Malformed refined hit_summary line (expected 6 fields, got {}): {}",
                fields.len(),
                line
            );
            continue;
        }

        let read_id = fields[0].to_string();

        let taxid: Taxid = match fields[2].parse() {
            Ok(t) if t > 0 => t,
            Ok(t) => {
                debug!("Read {} has non-positive taxid {} → skipping", read_id, t);
                continue;
            }
            Err(e) => {
                warn!("Invalid taxid '{}' in refined hit_summary (parse error: {}): {}", fields[2], e, line);
                continue;
            }
        };

        read_to_taxid.insert(read_id, taxid);
    }

    info!("Loaded {} read → taxid mappings from refined hit summary", read_to_taxid.len());

    while let Some(item) = refined_m8_stream.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);

        // Parse m8 line — defensive with warnings
        let m8 = if db_type == NT_TAG {
            match M8Record::parse_line_nt(&line) {
                Ok(m) => m,
                Err(e) => {
                    warn!("Failed to parse NT m8 line: {} — {}", e, line);
                    continue;
                }
            }
        } else {
            match M8Record::parse_line_nr(&line) {
                Ok(m) => m,
                Err(e) => {
                    warn!("Failed to parse NR m8 line: {} — {}", e, line);
                    continue;
                }
            }
        };

        let Some(&taxid) = read_to_taxid.get(&m8.qname) else {
            // This is expected for unmapped reads — not an error
            continue;
        };

        if taxid <= 0 {
            debug!("Read {} has invalid refined taxid {} — skipping m8 hit", m8.qname, taxid);
            continue;
        }

        let lineage = match lineage_map.get(&taxid) {
            Some(l) => l,
            None => {
                warn!("Missing lineage for taxid {} (read: {}) — skipping hit", taxid, m8.qname);
                continue;
            }
        };

        if !(should_keep_filter)(lineage) {
            debug!("Read {} (taxid {}) filtered out by should_keep_filter", m8.qname, taxid);
            continue;
        }

        let bucket = buckets.entry(taxid).or_default();
        bucket.nonunique_count += 1;
        bucket.unique_count +=  duplicate_clusters
            .get(&m8.qname)
            .map(|cluster| cluster.size)
            .unwrap_or(1u64);
        bucket.base_count += 1;
        bucket.sum_percent_identity += m8.pident;
        bucket.sum_alignment_length += m8.alen as f64;
        bucket.sum_e_value += m8.evalue;
        bucket.source_count_type.insert(db_type.to_string());
    }


    for (taxid, bucket) in buckets {
        let lineage = match lineage_map.get(&taxid) {
            Some(l) => l,
            None => {
                warn!("No lineage for final taxid {} — skipping output", taxid);
                continue;
            }
        };

        if !(should_keep_filter)(lineage) {
            continue;
        }

        let dcr = if bucket.nonunique_count > 0 {
            bucket.unique_count as f64 / bucket.nonunique_count as f64
        } else {
            0.0
        };

        let percent_identity = if bucket.base_count > 0 {
            bucket.sum_percent_identity / bucket.base_count as f64
        } else {
            0.0
        };

        let alignment_length = if bucket.base_count > 0 {
            bucket.sum_alignment_length / bucket.base_count as f64
        } else {
            0.0
        };

        let e_value = if bucket.base_count > 0 {
            bucket.sum_e_value / bucket.base_count as f64
        } else {
            0.0
        };

        let count = TaxonCount {
            tax_id: taxid,
            tax_level: 1, // species-level aggregation
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
            .map_err(|_| anyhow!("taxon count output channel dropped"))?;
    }

    Ok(())
}


pub async fn compute_merged_taxon_counts(
    config: Arc<RunConfig>,
    nt_m8_stream: mpsc::Receiver<ParseOutput>,
    nt_hit_summary_stream: mpsc::Receiver<ParseOutput>,
    nt_contig_summary: Vec<ContigSummaryEntry>,

    nr_m8_stream: mpsc::Receiver<ParseOutput>,
    nr_hit_summary_stream: mpsc::Receiver<ParseOutput>,
    nr_contig_summary: Vec<ContigSummaryEntry>,

    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,

    merged_m8_path: PathBuf,
    merged_hitsummary_path: PathBuf,
    merged_taxon_counts_path: PathBuf,
    merged_contig_summary_path: PathBuf,

    mut nr_alignment_per_read: HashMap<String, SpeciesAlignmentResults>,
) -> Result<()> {
    // Write merged m8 and hit summary to disk
    let mut merged_m8_file = BufWriter::new(TokioFile::create(&merged_m8_path).await?);
    let mut merged_hit_file = BufWriter::new(TokioFile::create(&merged_hitsummary_path).await?);

    // NT pass
    let mut nt_m8_stream = ReceiverStream::new(nt_m8_stream);
    let mut nt_hit_stream = ReceiverStream::new(nt_hit_summary_stream);

    while let (Some(m8_item), Some(hit_item)) = (nt_m8_stream.next().await, nt_hit_stream.next().await) {
        let m8_bytes = m8_item.to_bytes()?;
        let hit_bytes = hit_item.to_bytes()?;
        let hit_line = String::from_utf8_lossy(&hit_bytes);
        let hit_trimmed = hit_line.trim_end();
        if hit_trimmed.is_empty() {
            continue;
        }

        let hit_fields: Vec<&str> = hit_trimmed.split('\t').collect();
        if hit_fields.len() < 10 {
            warn!("Malformed NT hit line: {}", hit_trimmed);
            continue;
        }

        let read_id = hit_fields[0];

        let nt_contig = hit_fields.get(9).and_then(|s| s.parse::<Taxid>().ok());
        let nt_read = hit_fields.get(3).and_then(|s| s.parse::<Taxid>().ok());

        let nr_align = nr_alignment_per_read.get(read_id);
        let has_nr_contig = nr_align.map_or(false, |a| a.contig.is_some());
        let has_nr_read = nr_align.map_or(false, |a| a.read.is_some());

        if nt_contig.is_some() || (!has_nr_contig && nt_read.is_some()) {
            merged_m8_file.write_all(&m8_bytes).await?;
            merged_m8_file.write_all(b"\n").await?;

            let mut hit_with_source = hit_fields.iter().copied().collect::<Vec<&str>>();
            hit_with_source.push(NT_TAG);
            merged_hit_file
                .write_all(hit_with_source.join("\t").as_bytes())
                .await?;
            merged_hit_file.write_all(b"\n").await?;

            nr_alignment_per_read.remove(read_id);
        }
    }

    // Remaining NR pass
    let mut nr_m8_stream = ReceiverStream::new(nr_m8_stream);
    let mut nr_hit_stream = ReceiverStream::new(nr_hit_summary_stream);

    while let (Some(m8_item), Some(hit_item)) = (nr_m8_stream.next().await, nr_hit_stream.next().await) {
        let hit_bytes = hit_item.to_bytes()?;
        let hit_line = String::from_utf8_lossy(&hit_bytes);
        let hit_trimmed = hit_line.trim_end();
        if hit_trimmed.is_empty() {
            continue;
        }

        let hit_fields: Vec<&str> = hit_trimmed.split('\t').collect();
        let read_id = hit_fields[0];

        if nr_alignment_per_read.contains_key(read_id) {
            let m8_bytes = m8_item.to_bytes()?;
            merged_m8_file.write_all(&m8_bytes).await?;
            merged_m8_file.write_all(b"\n").await?;

            let mut hit_with_source = hit_fields.iter().copied().collect::<Vec<&str>>();
            hit_with_source.push(NR_TAG);
            merged_hit_file
                .write_all(hit_with_source.join("\t").as_bytes())
                .await?;
            merged_hit_file.write_all(b"\n").await?;
        }
    }

    merged_m8_file.flush().await?;
    merged_hit_file.flush().await?;
    drop(merged_m8_file);
    drop(merged_hit_file);
    info!("Merged alignment files written to disk");

    // Stream merged files back for taxon counting
    let (m8_tx, m8_rx) = mpsc::channel(4096);
    let (hit_tx, hit_rx) = mpsc::channel(4096);

    let m8_path_clone = merged_m8_path.clone();
    let m8_task = tokio::spawn(async move {
        let file = TokioFile::open(m8_path_clone).await?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line).await? > 0 {
            if !line.trim_end().is_empty() {
                let _ = m8_tx.send(ParseOutput::Bytes(Arc::new(line.as_bytes().to_vec()))).await;
            }
            line.clear();
        }
        Ok::<(), anyhow::Error>(())
    });

    let hit_path_clone = merged_hitsummary_path.clone();
    let hit_task = tokio::spawn(async move {
        let file = TokioFile::open(hit_path_clone).await?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line).await? > 0 {
            if !line.trim_end().is_empty() {
                let _ = hit_tx.send(ParseOutput::Bytes(Arc::new(line.as_bytes().to_vec()))).await;
            }
            line.clear();
        }
        Ok::<(), anyhow::Error>(())
    });

    let merged_taxon_counts_path_clone = merged_taxon_counts_path.clone();
    let config_clone = config.clone();

    let (json_tx, json_rx) = mpsc::channel(1024);
    let write_task = tokio::spawn(async move {
        write_byte_stream_to_file(
            &merged_taxon_counts_path_clone,
            ReceiverStream::new(json_rx),
            Some(config_clone.base_buffer_size),
        ).await
    });

    generate_taxon_count_json_from_m8(
        ReceiverStream::new(m8_rx),
        ReceiverStream::new(hit_rx),
        "merged_NT_NR",
        lineage_map.clone(),
        should_keep_filter.clone(),
        duplicate_clusters.clone(),
        json_tx,
    )
        .await?;

    m8_task.await??;
    hit_task.await??;
    write_task.await??;
    info!("Merged taxon counts generated");

    // Contig summary merge
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

    let final_contigs: Vec<ContigSummaryEntry> = merged_contigs.into_values().flatten().collect();

    let json = serde_json::to_string_pretty(&final_contigs)?;
    let mut out = BufWriter::new(TokioFile::create(&merged_contig_summary_path).await?);
    out.write_all(json.as_bytes()).await?;
    out.flush().await?;
    info!("Merged contig summary written ({} entries)", final_contigs.len());

    info!("compute_merged_taxon_counts complete");
    Ok(())
}
