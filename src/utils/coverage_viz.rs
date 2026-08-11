//! Coverage visualization utilities.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, Context, Result};
use futures::stream::{FuturesUnordered, StreamExt};
use log::{self, warn, info, debug};
use serde::Serialize;
use serde_json::{json, Value};
use tokio::fs::create_dir_all;
use tokio_stream::wrappers::{ReceiverStream, UnboundedReceiverStream};

use crate::utils::blast::M8Record;
use crate::utils::streams::{ParseOutput, ToBytes};

const MAX_NUM_BINS_COVERAGE: usize = 500;
const NUM_ACCESSIONS_PER_TAXON: usize = 10;
const MIN_CONTIG_SIZE: u64 = 4;

#[derive(Debug, Clone, Default)]
struct TaxonData {
    accessions: HashSet<String>,
    num_total_accessions: usize,
    best_accessions: Vec<String>,
}

#[derive(Debug, Clone, Default)]
struct AccessionData {
    reads: Vec<String>,
    contigs: Vec<String>,
    name: String,
    total_length: u64,
    score: f64,
}

#[derive(Debug, Clone)]
struct ContigHit {
    accession: String,
    subject_start: u64,
    subject_end: u64,
    query_start: u64,
    query_end: u64,
    alignment_length: u64,
    percent_id: f64,
    num_mismatches: u64,
    num_gaps: u64,
    total_length: u64,       // bp — set in augment_contig_data_with_coverage
    coverage: Vec<f64>,
    prop_mismatch: f64,
    num_reads: u64,          // from valid_contigs
    byterange: Option<(u64, u64)>, // [offset, length] in contigs.fasta; set by byterange augment
}

#[derive(Debug, Clone)]
struct ReadHit {
    accession: String,
    subject_start: u64,
    subject_end: u64,
    alignment_length: u64,
    percent_id: f64,
    num_mismatches: u64,
    num_gaps: u64,
    prop_mismatch: f64,
}

// Output JSON structures
#[derive(Serialize)]
struct CoverageVizData {
    total_length: u64,
    name: String,
    hit_groups: Value,
    coverage: Vec<[f64; 5]>, // [bin_index, depth, breadth%, num_contigs, num_reads]
    coverage_bin_size: f64,
    max_aligned_length: u64,
    coverage_depth: f64,
    coverage_breadth: f64,
    avg_prop_mismatch: f64,
}

#[derive(Serialize)]
struct CoverageVizSummary {
    #[serde(flatten)]
    taxons: HashMap<String, TaxonSummary>,
}

#[derive(Serialize)]
struct TaxonSummary {
    best_accessions: Vec<AccessionSummary>,
    num_accessions: usize,
}

#[derive(Serialize)]
struct AccessionSummary {
    id: String,
    name: String,
    num_contigs: usize,
    num_reads: usize,
    score: f64,
    coverage_breadth: f64,
    coverage_depth: f64,
}

#[derive(Debug, Clone, Default)]
struct CoverageBin {
    depth: f64,
    endpoints: Vec<(f64, i8)>,
    num_contigs: u64,
    num_reads: u64,
}

pub async fn generate_coverage_viz(
    refined_gsnap_hitsummary2_tab: ReceiverStream<ParseOutput>,
    refined_gsnap_blast_top_m8: UnboundedReceiverStream<ParseOutput>,
    contig_coverage_json: PathBuf,
    contig_stats_json: PathBuf,
    contigs_fasta: PathBuf,
    gsnap_deduped_m8: ReceiverStream<ParseOutput>,
    nt_info_db: PathBuf,
    output_dir: PathBuf,
    max_num_bins_coverage: Option<usize>,
    num_accessions_per_taxon: Option<usize>,
    min_contig_size: Option<u64>,
    keep_taxons_with_no_contigs: bool,
) -> Result<()> {
    let max_bins = max_num_bins_coverage.unwrap_or(MAX_NUM_BINS_COVERAGE);
    let acc_per_taxon = num_accessions_per_taxon.unwrap_or(NUM_ACCESSIONS_PER_TAXON);
    let min_size = min_contig_size.unwrap_or(MIN_CONTIG_SIZE);

    let info_dict = load_nt_info_db(&nt_info_db)?;

    let (taxon_data, accession_data, contig_data, read_data) = prepare_data(
        refined_gsnap_hitsummary2_tab,
        refined_gsnap_blast_top_m8,
        &contig_coverage_json,
        &contig_stats_json,
        &contigs_fasta,
        gsnap_deduped_m8,
        info_dict,
        min_size,
        acc_per_taxon,
        keep_taxons_with_no_contigs,
    )
    .await?;

    let coverage_viz_data =
        generate_coverage_viz_data(&accession_data, &contig_data, &read_data, max_bins).await?;

    let summary_data =
        generate_coverage_viz_summary_data(&taxon_data, &accession_data, &coverage_viz_data);

    let summary_path = output_dir.join("coverage_viz_summary.json");
    let mut summary_file = BufWriter::new(File::create(&summary_path)?);
    serde_json::to_writer_pretty(&mut summary_file, &summary_data)?;

    let viz_dir = output_dir.join("coverage_viz");
    create_dir_all(&viz_dir).await?;

    let mut tasks = FuturesUnordered::new();
    for (acc_id, data) in coverage_viz_data {
        let viz_dir = viz_dir.clone();
        tasks.push(tokio::task::spawn_blocking(move || {
            let path = viz_dir.join(format!("{}_coverage_viz.json", acc_id));
            let file = File::create(path)?;
            let mut writer = BufWriter::new(file);
            serde_json::to_writer_pretty(&mut writer, &data)
                .map_err(|e| anyhow!("Failed writing {}: {}", acc_id, e))
        }));
    }

    while let Some(result) = tasks.next().await {
        let result: Result<Result<()>, tokio::task::JoinError> = result;
        result.context("Task join error")??;
    }

    Ok(())
}

async fn prepare_data(
    hit_summary: ReceiverStream<ParseOutput>,
    blast_top_m8: UnboundedReceiverStream<ParseOutput>,
    contig_coverage_json: &Path,
    contig_stats_json: &Path,
    contigs_fasta: &Path,
    gsnap_deduped_m8: ReceiverStream<ParseOutput>,
    info_dict: HashMap<String, (String, u64)>,
    min_contig_size: u64,
    num_accessions_per_taxon: usize,
    keep_taxons_with_no_contigs: bool,
) -> Result<(
    HashMap<String, TaxonData>,
    HashMap<String, AccessionData>,
    HashMap<String, Vec<ContigHit>>,
    HashMap<String, Vec<ReadHit>>,
)> {
    let valid_contigs = get_valid_contigs_with_read_counts(contig_stats_json, min_contig_size)?;

    let (mut accession_data, mut taxon_data) =
        generate_accession_data(hit_summary, &valid_contigs).await?;

    if !keep_taxons_with_no_contigs {
        remove_taxons_with_no_contigs(&mut accession_data, &mut taxon_data);
    }

    augment_accession_data_with_info(&info_dict, &mut accession_data);

    let assigned_reads = get_unassigned_reads_set(&accession_data); // actually assigned, per Python

    let mut contig_data = generate_contig_data(blast_top_m8, &valid_contigs).await?;

    let read_data = generate_read_data(gsnap_deduped_m8, &assigned_reads).await?;

    augment_contig_data_with_coverage(contig_coverage_json, &mut contig_data)?;
    augment_contig_data_with_byteranges(contigs_fasta, &mut contig_data)?;

    let (taxon_data, accession_data) = select_best_accessions_per_taxon(
        taxon_data,
        accession_data,
        &contig_data,
        num_accessions_per_taxon,
    );

    info!(
  "[coverage prepare] valid_contigs={} accession_data={} taxon_data={} contig_data={} read_data={}",
  valid_contigs.len(),
  accession_data.len(),
  taxon_data.len(),
  contig_data.len(),
  read_data.len()
);
    let n_contig_acc = accession_data.values().filter(|a| !a.contigs.is_empty()).count();
    info!("[coverage prepare] accessions_with_contigs={}", n_contig_acc);

    Ok((taxon_data, accession_data, contig_data, read_data))
}

fn load_nt_info_db(path: &Path) -> Result<HashMap<String, (String, u64)>> {
    let file = File::open(path)
        .with_context(|| format!("open nt_info {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    let mut bad = 0u64;

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            bad += 1;
            continue;
        }
        let len_str = parts[parts.len() - 1];
        let Ok(len) = len_str.parse::<u64>() else {
            bad += 1;
            continue;
        };
        let acc_raw = parts[0];
        let name = if parts.len() == 3 {
            parts[1].to_string()
        } else {
            parts[1..parts.len() - 1].join(" ")
        };

        // Prefer file already unversioned; still accept versioned keys.
        let acc = acc_raw.split('.').next().unwrap_or(acc_raw).to_string();
        map.entry(acc).or_insert((name, len));
    }

    info!(
        "[load_nt_info_db] path={} entries={} bad_lines={}",
        path.display(),
        map.len(),
        bad
    );
    Ok(map)
}

fn get_valid_contigs_with_read_counts(path: &Path, min_size: u64) -> Result<HashMap<String, u64>> {
    let file = File::open(path)?;
    let mut contents = String::new();
    std::io::Read::read_to_string(&mut BufReader::new(file), &mut contents)?;
    let counts: HashMap<String, u64> = serde_json::from_str(&contents)?;
    Ok(counts
        .into_iter()
        .filter(|(c, cnt)| c != "*" && *cnt >= min_size)
        .collect())
}

async fn generate_accession_data(
    mut hit_summary: ReceiverStream<ParseOutput>,
    valid_contigs: &HashMap<String, u64>,
) -> Result<(HashMap<String, AccessionData>, HashMap<String, TaxonData>)> {
    let mut accession_data: HashMap<String, AccessionData> = HashMap::new();
    let mut taxon_data: HashMap<String, TaxonData> = HashMap::new();
    // Contigs deduped per accession (Python uses a set)
    let mut contig_sets: HashMap<String, HashSet<String>> = HashMap::new();

    let mut line_count = 0u64;

    while let Some(item) = hit_summary.next().await {
        line_count += 1;
        if line_count % 100_000 == 0 {
            log::info!("Processed {} hit_summary lines", line_count);
        }

        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.is_empty() || (fields.len() == 1 && fields[0].is_empty()) {
            continue;
        }

        if fields.len() >= 12 && valid_contigs.contains_key(fields[7]) {
            let taxon = fields[9].to_string();
            let acc = fields[8].to_string();
            let contig = fields[7].to_string();

            taxon_data
                .entry(taxon)
                .or_default()
                .accessions
                .insert(acc.clone());
            contig_sets.entry(acc.clone()).or_default().insert(contig);
            accession_data.entry(acc).or_default();
        } else if fields.len() >= 5 {
            // Includes 12-col rows whose contig is not in valid_contigs (Python else branch)
            let taxon = fields[4].to_string();
            let acc = fields[3].to_string();
            let read = fields[0].to_string();

            taxon_data
                .entry(taxon)
                .or_default()
                .accessions
                .insert(acc.clone());
            accession_data.entry(acc).or_default().reads.push(read);
        }
    }

    for (acc, set) in contig_sets {
        if let Some(ad) = accession_data.get_mut(&acc) {
            ad.contigs = set.into_iter().collect();
        }
    }

    for td in taxon_data.values_mut() {
        td.num_total_accessions = td.accessions.len();
    }

    Ok((accession_data, taxon_data))
}

fn remove_taxons_with_no_contigs(
    accession_data: &mut HashMap<String, AccessionData>,
    taxon_data: &mut HashMap<String, TaxonData>,
) {
    let mut to_remove = Vec::new();
    for (taxon, td) in taxon_data.iter() {
        let contig_cnt: usize = td
            .accessions
            .iter()
            .filter_map(|acc| accession_data.get(acc))
            .map(|ad| ad.contigs.len())
            .sum();
        if contig_cnt == 0 {
            to_remove.push(taxon.clone());
        }
    }
    for taxon in to_remove {
        if let Some(td) = taxon_data.remove(&taxon) {
            for acc in td.accessions {
                accession_data.remove(&acc);
            }
        }
    }
}

fn augment_accession_data_with_info(
    info_dict: &HashMap<String, (String, u64)>,
    accession_data: &mut HashMap<String, AccessionData>,
) {
    for (acc_id, ad) in accession_data.iter_mut() {
        if let Some((name, len)) = info_dict.get(acc_id) {
            ad.name = name.clone();
            ad.total_length = *len;
        } else {
            ad.name = "Unknown accession".to_string();
            ad.total_length = 0;
        }
    }
}

fn get_unassigned_reads_set(accession_data: &HashMap<String, AccessionData>) -> HashSet<String> {
    let mut set = HashSet::new();
    for ad in accession_data.values() {
        set.extend(ad.reads.iter().cloned());
    }
    set
}

async fn generate_contig_data(
    mut blast_top_m8: UnboundedReceiverStream<ParseOutput>,
    valid_contigs: &HashMap<String, u64>,
) -> Result<HashMap<String, Vec<ContigHit>>> {
    let mut contig_data: HashMap<String, Vec<ContigHit>> = HashMap::new();
    let mut line_count = 0u64;
    let mut parse_err = 0u64;

    while let Some(item) = blast_top_m8.next().await {
        line_count += 1;
        if line_count % 100_000 == 0 {
            info!("Processed {} blast_top_m8 lines", line_count);
        }

        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let line_trim = line.trim_end();
        if line_trim.is_empty() {
            continue;
        }

        let m8 = match M8Record::parse_line_nt(line_trim)
            .or_else(|_| M8Record::parse_line_nr(line_trim))
        {
            Ok(record) => record,
            Err(e) => {
                parse_err += 1;
                if parse_err <= 20 {
                    warn!("Failed to parse blast_top_m8 line: {} — {}", line_trim, e);
                }
                continue;
            }
        };

        let contig_name = m8.qname.clone();
        if !valid_contigs.contains_key(&contig_name) {
            continue;
        }

        // total_length is contig length in bp — filled later from coverage JSON
        // (valid_contigs values are read counts, not lengths)
        let prop_mismatch = if m8.alen > 0 {
            m8.mismatch as f64 / m8.alen as f64
        } else {
            0.0
        };

        contig_data.entry(contig_name.clone()).or_default().push(ContigHit {
            accession: m8.tname.clone(),
            subject_start: m8.tstart,
            subject_end: m8.tend,
            query_start: m8.qstart,
            query_end: m8.qend,
            alignment_length: m8.alen,
            percent_id: m8.pident,
            num_mismatches: m8.mismatch,
            num_gaps: m8.gapopen,
            total_length: 0, // bp — set in augment_contig_data_with_coverage
            coverage: vec![],
            prop_mismatch,
            num_reads: *valid_contigs.get(&contig_name).unwrap_or(&0),
            byterange: None, // set in augment_contig_data_with_byteranges
        });
    }

    info!(
        "[generate_contig_data] done: lines={} parse_err={} contigs={}",
        line_count,
        parse_err,
        contig_data.len()
    );
    Ok(contig_data)
}

async fn generate_read_data(
    mut gsnap_deduped_m8: ReceiverStream<ParseOutput>,
    assigned_reads: &HashSet<String>,
) -> Result<HashMap<String, Vec<ReadHit>>> {
    let mut read_data: HashMap<String, Vec<ReadHit>> = HashMap::new();
    let mut line_count = 0u64;
    let mut parse_err = 0u64;

    while let Some(item) = gsnap_deduped_m8.next().await {
        line_count += 1;
        if line_count % 100_000 == 0 {
            info!("Processed {} gsnap_deduped_m8 lines", line_count);
        }

        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let line_trim = line.trim_end();
        if line_trim.is_empty() {
            continue;
        }
        
        let m8 = match M8Record::parse_line_nt(line_trim)
            .or_else(|_| M8Record::parse_line_nr(line_trim))
        {
            Ok(record) => record,
            Err(e) => {
                parse_err += 1;
                if parse_err <= 20 {
                    warn!("Failed to parse deduped m8 line: {} — {}", line_trim, e);
                }
                continue;
            }
        };

        let read_name = m8.qname.clone();
        if !assigned_reads.contains(&read_name) {
            continue; // not a loose/unassigned read for coverage viz
        }

        let prop_mismatch = if m8.alen > 0 {
            m8.mismatch as f64 / m8.alen as f64
        } else {
            0.0
        };

        read_data.entry(read_name).or_default().push(ReadHit {
            accession: m8.tname.clone(),
            subject_start: m8.tstart,
            subject_end: m8.tend,
            alignment_length: m8.alen,
            percent_id: m8.pident,
            num_mismatches: m8.mismatch,
            num_gaps: m8.gapopen,
            prop_mismatch,
        });
    }

    info!(
        "[generate_read_data] done: lines={} parse_err={} reads={}",
        line_count,
        parse_err,
        read_data.len()
    );
    Ok(read_data)
}
fn augment_contig_data_with_coverage(
    path: &Path,
    contig_data: &mut HashMap<String, Vec<ContigHit>>,
) -> Result<()> {
    let file = File::open(path)?;
    let root: Value = serde_json::from_reader(BufReader::new(file))
        .map_err(|e| anyhow!("Failed to parse contig_coverage_json: {}", e))?;

    let obj = root
        .as_object()
        .ok_or_else(|| anyhow!("contig_coverage_json root is not an object"))?;

    for (name, v) in obj {
        let (depths, contig_len) = match v {
            Value::Array(a) => {
                let depths: Vec<f64> = a.iter().filter_map(|x| x.as_f64()).collect();
                let len = depths.len() as u64; // 1 depth per base
                (depths, len)
            }
            Value::Object(o) => {
                let arr = o
                    .get("coverage")
                    .or_else(|| o.get("depths"))
                    .or_else(|| o.get("depth"))
                    .and_then(|x| x.as_array())
                    .ok_or_else(|| {
                        anyhow!(
                            "contig_coverage entry '{}' missing coverage[]",
                            name
                        )
                    })?;
                let depths: Vec<f64> = arr.iter().filter_map(|x| x.as_f64()).collect();
                // Prefer explicit contig_len; fall back to coverage array length
                let len = o
                    .get("contig_len")
                    .or_else(|| o.get("total_length"))
                    .or_else(|| o.get("length"))
                    .and_then(|x| x.as_u64())
                    .unwrap_or(depths.len() as u64);
                (depths, len)
            }
            _ => {
                warn!("Skipping contig_coverage entry '{}': unexpected JSON type", name);
                continue;
            }
        };

        if let Some(hits) = contig_data.get_mut(name) {
            for hit in hits.iter_mut() {
                hit.coverage = depths.clone();
                hit.total_length = contig_len;
            }
        }
        // contigs present only in coverage JSON (no blast hit) are intentionally ignored
    }

    // Contigs that had blast hits but no coverage entry stay at total_length=0;
    // calculate_accession_coverage already guards on empty coverage / zero length.
    Ok(())
}

fn select_best_accessions_per_taxon(
    mut taxon_data: HashMap<String, TaxonData>,
    mut accession_data: HashMap<String, AccessionData>,
    contig_data: &HashMap<String, Vec<ContigHit>>,
    num_per_taxon: usize,
) -> (HashMap<String, TaxonData>, HashMap<String, AccessionData>) {
    // Score = max_contig_alen + sum_contig_alen + num_reads  (Python get_score)
    // Do not zero out score when total_length == 0 — contig accessions still rank by alen.
    for ad in accession_data.values_mut() {
        let mut max_alen = 0u64;
        let mut sum_alen = 0u64;
        for cname in &ad.contigs {
            if let Some(hits) = contig_data.get(cname) {
                for h in hits {
                    max_alen = max_alen.max(h.alignment_length);
                    sum_alen += h.alignment_length;
                }
            }
        }
        ad.score = max_alen as f64 + sum_alen as f64 + ad.reads.len() as f64;
    }

    let mut filtered_accessions: HashMap<String, AccessionData> = HashMap::new();

    for td in taxon_data.values_mut() {
        let mut sorted: Vec<String> = td.accessions.iter().cloned().collect();
        sorted.sort_by(|a, b| {
            let sa = accession_data
                .get(a)
                .map(|x| x.score)
                .unwrap_or(f64::NEG_INFINITY);
            let sb = accession_data
                .get(b)
                .map(|x| x.score)
                .unwrap_or(f64::NEG_INFINITY);
            sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
        });

        // Python: len(contigs) >= 1 — no total_length check
        let with_contigs: Vec<String> = sorted
            .iter()
            .filter(|a| {
                accession_data
                    .get(*a)
                    .map(|x| !x.contigs.is_empty())
                    .unwrap_or(false)
            })
            .cloned()
            .collect();

        let best = if with_contigs.len() >= num_per_taxon {
            with_contigs
        } else {
            // Python: fill with accessions that have zero contigs
            let without: Vec<String> = sorted
                .iter()
                .filter(|a| {
                    accession_data
                        .get(*a)
                        .map(|x| x.contigs.is_empty())
                        .unwrap_or(false)
                })
                .cloned()
                .collect();
            let need = num_per_taxon - with_contigs.len();
            let mut v = with_contigs;
            v.extend(without.into_iter().take(need));
            v.sort_by(|a, b| {
                let sa = accession_data
                    .get(a)
                    .map(|x| x.score)
                    .unwrap_or(f64::NEG_INFINITY);
                let sb = accession_data
                    .get(b)
                    .map(|x| x.score)
                    .unwrap_or(f64::NEG_INFINITY);
                sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
            });
            v
        };

        for a in &best {
            if let Some(ad) = accession_data.remove(a) {
                filtered_accessions.insert(a.clone(), ad);
            }
        }
        td.best_accessions = best;
    }

    taxon_data.retain(|_, td| !td.best_accessions.is_empty());

    (taxon_data, filtered_accessions)
}

async fn generate_coverage_viz_data(
    accession_data: &HashMap<String, AccessionData>,
    contig_data: &HashMap<String, Vec<ContigHit>>,
    read_data: &HashMap<String, Vec<ReadHit>>,
    max_num_bins: usize,
) -> Result<HashMap<String, CoverageVizData>> {
    let contig_data = Arc::new(contig_data.clone());
    let read_data = Arc::new(read_data.clone());

    let mut tasks = FuturesUnordered::new();

    for (acc_id, acc_obj) in accession_data {
        let acc_id = acc_id.clone();
        let acc_obj = acc_obj.clone();
        let contig_data = Arc::clone(&contig_data);
        let read_data = Arc::clone(&read_data);

        if acc_obj.total_length == 0 {
            warn!(
            "Skipping zero-length / unknown accession {} (no nt_info entry)",
            acc_id
        );
            continue;
        }

        tasks.push(tokio::task::spawn_blocking(move || {
            let total_len = acc_obj.total_length as f64;

            let num_bins = max_num_bins.min(total_len as usize);
            let bin_size = total_len / num_bins as f64;

            let (coverage, _) = calculate_accession_coverage(
                &acc_id,
                &acc_obj,
                &contig_data,
                &read_data,
                num_bins,
                bin_size,
            )?;

            let stats = calculate_accession_stats(&acc_obj, &contig_data, &read_data, total_len)?;

            let hit_groups = generate_hit_group_json(
                &acc_obj,
                &acc_id,
                &contig_data,
                &read_data,
                num_bins,
            );

            let viz = CoverageVizData {
                total_length: acc_obj.total_length,
                name: acc_obj.name.clone(),
                hit_groups: hit_groups,
                coverage,
                coverage_bin_size: bin_size,
                max_aligned_length: stats.max_aligned_length,
                coverage_depth: format_number(stats.coverage_depth),
                coverage_breadth: format_percent(stats.coverage_breadth),
                avg_prop_mismatch: format_percent(stats.avg_prop_mismatch),
            };

            Ok((acc_id, viz))
        }));
    }

    let mut result = HashMap::with_capacity(accession_data.len());
    while let Some(task_result) = tasks.next().await {
        let task_result: Result<Result<(String, CoverageVizData)>, tokio::task::JoinError> =
            task_result;
        let (acc_id, viz) = task_result.context("Task join error")??;
        result.insert(acc_id, viz);
    }

    Ok(result)
}

fn calculate_accession_coverage(
    acc_id: &str,
    acc_obj: &AccessionData,
    contig_data: &HashMap<String, Vec<ContigHit>>,
    read_data: &HashMap<String, Vec<ReadHit>>,
    num_bins: usize,
    bin_size: f64,
) -> Result<(Vec<[f64; 5]>, f64)> {
    let mut coverage = vec![CoverageBin::default(); num_bins];

    for contig_name in &acc_obj.contigs {
        if let Some(hits) = contig_data.get(contig_name) {
            for hit in hits {
                if hit.accession != acc_id {
                    continue;
                }
                let (s_start, s_end) =
                    decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
                let (bin_start, bin_end) = align_interval(s_start / bin_size, s_end / bin_size);

                let bin_start_i = floor_with_min(bin_start, 0);
                let bin_end_i = ceil_with_max(bin_end, num_bins as i64);

                for i in bin_start_i..bin_end_i {
                    let i = i as usize;
                    let acc_interval = [
                        bin_size * bin_start.max(i as f64),
                        bin_size * bin_end.min((i + 1) as f64),
                    ];

                    let contig_interval = transform_interval(
                        &[acc_interval[0], acc_interval[1]],
                        s_start,
                        s_end,
                        hit.query_start as f64,
                        hit.query_end as f64,
                    );

                    let coverage_interval = if hit.total_length as usize == hit.coverage.len() {
                        align_interval(contig_interval.0, contig_interval.1)
                    } else {
                        let tmp = transform_interval(
                            &[contig_interval.0, contig_interval.1],
                            0.0,
                            hit.total_length as f64,
                            0.0,
                            hit.coverage.len() as f64,
                        );
                        align_interval(tmp.0, tmp.1)
                    };

                    let cov_start = floor_with_min(coverage_interval.0, 0) as usize;
                    let cov_end =
                        ceil_with_max(coverage_interval.1, hit.coverage.len() as i64) as usize;

                    if cov_end > cov_start {
                        let avg_cov = hit.coverage[cov_start..cov_end].iter().sum::<f64>()
                            / (cov_end - cov_start) as f64;
                        let proportion = (acc_interval[1] - acc_interval[0]) / bin_size;
                        let contrib = avg_cov * proportion;

                        coverage[i].depth += contrib;
                        coverage[i].endpoints.push((acc_interval[0], 1));
                        coverage[i].endpoints.push((acc_interval[1], -1));
                        coverage[i].num_contigs += 1;
                    }
                }
            }
        }
    }

    for read_name in &acc_obj.reads {
        if let Some(hits) = read_data.get(read_name) {
            for hit in hits {
                if hit.accession != acc_id {
                    continue;
                }
                let (s_start, s_end) =
                    decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
                let (bin_start, bin_end) = align_interval(s_start / bin_size, s_end / bin_size);

                let bin_start_i = floor_with_min(bin_start, 0);
                let bin_end_i = ceil_with_max(bin_end, num_bins as i64);

                for i in bin_start_i..bin_end_i {
                    let i = i as usize;
                    let acc_range = [
                        bin_size * bin_start.max(i as f64),
                        bin_size * bin_end.min((i + 1) as f64),
                    ];
                    let contrib = (acc_range[1] - acc_range[0]) / bin_size;

                    coverage[i].depth += contrib;
                    coverage[i].endpoints.push((acc_range[0], 1));
                    coverage[i].endpoints.push((acc_range[1], -1));
                    coverage[i].num_reads += 1;
                }
            }
        }
    }

    let mut final_coverage = Vec::new();
    for (i, bin) in coverage.iter().enumerate() {
        if bin.depth > 0.0 {
            let breadth = calculate_covered_length(&bin.endpoints) / bin_size;
            final_coverage.push([
                i as f64,
                format_number(bin.depth),
                format_percent(breadth),
                bin.num_contigs as f64,
                bin.num_reads as f64,
            ]);
        }
    }

    Ok((final_coverage, bin_size))
}

fn calculate_covered_length(endpoints: &[(f64, i8)]) -> f64 {
    let mut sorted = endpoints.to_vec();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut total = 0.0;
    let mut depth = 0;
    let mut last_start = 0.0;

    for &(pos, change) in &sorted {
        depth += change;
        if change == 1 && depth == 1 {
            last_start = pos;
        } else if change == -1 && depth == 0 {
            total += pos - last_start;
        }
    }
    total
}

struct AccessionStats {
    max_aligned_length: u64,
    coverage_depth: f64,
    coverage_breadth: f64,
    avg_prop_mismatch: f64,
}

fn calculate_accession_stats(
    acc_obj: &AccessionData,
    contig_data: &HashMap<String, Vec<ContigHit>>,
    read_data: &HashMap<String, Vec<ReadHit>>,
    total_len: f64,
) -> Result<AccessionStats> {
    let mut max_len = 0u64;
    let mut cov_sum = 0.0;
    let mut mismatch_sum = 0.0;
    let mut endpoints = Vec::new();

    for contig in &acc_obj.contigs {
        if let Some(hits) = contig_data.get(contig) {
            for hit in hits {
                let (dec_start, dec_end) =
                    decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
                let (s_start, s_end) = align_interval(dec_start, dec_end);
                let aligned_len = (s_end - s_start) as u64;
                max_len = max_len.max(aligned_len);

                let (q_dec_start, q_dec_end) =
                    decrement_lower_bound(hit.query_start as f64, hit.query_end as f64);
                let (q_start, q_end) = align_interval(q_dec_start, q_dec_end);
                if !hit.coverage.is_empty() {
                    let lo = (q_start as usize).min(hit.coverage.len());
                    let hi = (q_end as usize).min(hit.coverage.len()).max(lo);
                    cov_sum += hit.coverage[lo..hi].iter().sum::<f64>();
                }

                mismatch_sum += hit.prop_mismatch;
                endpoints.push((s_start, 1));
                endpoints.push((s_end, -1));
            }
        }
    }

    for read in &acc_obj.reads {
        if let Some(hits) = read_data.get(read) {
            for hit in hits {
                let (dec_start, dec_end) =
                    decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
                let (s_start, s_end) = align_interval(dec_start, dec_end);
                let aligned_len = (s_end - s_start) as u64;
                max_len = max_len.max(aligned_len);

                cov_sum += aligned_len as f64;
                mismatch_sum += hit.prop_mismatch;
                endpoints.push((s_start, 1));
                endpoints.push((s_end, -1));
            }
        }
    }

    let breadth = if total_len > 0.0 {
        calculate_covered_length(&endpoints) / total_len
    } else {
        0.0
    };

    //prop_total_mismatch / (len(contigs) + len(reads))
    let n_names = acc_obj.contigs.len() + acc_obj.reads.len();

    Ok(AccessionStats {
        max_aligned_length: max_len,
        coverage_depth: if total_len > 0.0 {
            cov_sum / total_len
        } else {
            0.0
        },
        coverage_breadth: breadth,
        avg_prop_mismatch: if n_names > 0 {
            mismatch_sum / n_names as f64
        } else {
            0.0
        },
    })
}

fn generate_coverage_viz_summary_data(
    taxon_data: &HashMap<String, TaxonData>,
    accession_data: &HashMap<String, AccessionData>,
    coverage_viz: &HashMap<String, CoverageVizData>,
) -> CoverageVizSummary {
    let mut taxons = HashMap::new();

    for (taxon_id, td) in taxon_data {
        let best_acc: Vec<AccessionSummary> = td
            .best_accessions
            .iter()
            .filter_map(|acc_id| {
                let ad = accession_data.get(acc_id)?;
                let cv = coverage_viz.get(acc_id)?;
                Some(AccessionSummary {
                    id: acc_id.clone(),
                    name: ad.name.clone(),
                    num_contigs: ad.contigs.len(),
                    num_reads: ad.reads.len(),
                    score: format_number(ad.score),
                    coverage_breadth: cv.coverage_breadth,
                    coverage_depth: cv.coverage_depth,
                })
            })
            .collect();

        taxons.insert(
            taxon_id.clone(),
            TaxonSummary {
                best_accessions: best_acc,
                num_accessions: td.num_total_accessions,
            },
        );
    }

    CoverageVizSummary { taxons }
}

/// Aggregate reads/contigs into hit groups for one accession.
/// Layout of each group matches Python get_hit_group_json:
/// [num_contigs, num_reads, contig_r, start, end, avg_alen, avg_pident,
///  avg_mismatch, avg_gaps, bin_index, contig_byteranges]
fn generate_hit_group_json(
    accession_obj: &AccessionData,
    accession_id: &str,
    contig_data: &HashMap<String, Vec<ContigHit>>,
    read_data: &HashMap<String, Vec<ReadHit>>,
    num_bins: usize,
) -> Value {
    if accession_obj.total_length == 0 || num_bins == 0 {
        return json!([]);
    }

    let bin_size = accession_obj.total_length as f64 / num_bins as f64;

    // Individual hits (span >= bin_size) vs binned small hits
    let mut individual_reads: Vec<(String, usize)> = Vec::new();
    let mut individual_contigs: Vec<(String, usize)> = Vec::new();
    let mut read_bins: Vec<Vec<(String, usize)>> = vec![Vec::new(); num_bins];
    let mut contig_bins: Vec<Vec<(String, usize)>> = vec![Vec::new(); num_bins];

    // ── process one hit ──────────────────────────────────────────────
    let mut process_contig = |name: &str| {
        let Some(hits) = contig_data.get(name) else {
            warn!("Could not find contig in map: {}", name);
            return;
        };
        for (ind, hit) in hits.iter().enumerate() {
            if hit.accession != accession_id {
                warn!(
                    "Mismatched accession for {}: {} (hit) versus {} (hitsummary)",
                    name, hit.accession, accession_id
                );
                continue;
            }
            let (start, end) =
                align_interval_u64(hit.subject_start, hit.subject_end);
            let (dec_s, dec_e) = decrement_lower_bound(start as f64, end as f64);
            let (acc_s, acc_e) = align_interval(dec_s, dec_e);
            let span = acc_e - acc_s;

            if span >= bin_size {
                individual_contigs.push((name.to_string(), ind));
            } else {
                let mid = (acc_s + acc_e) / 2.0;
                let bin_idx = floor_with_min(mid / bin_size, 0) as usize;
                let bin_idx = bin_idx.min(num_bins.saturating_sub(1));
                contig_bins[bin_idx].push((name.to_string(), ind));
            }
        }
    };

    let mut process_read = |name: &str| {
        let Some(hits) = read_data.get(name) else {
            warn!("Could not find read in map: {}", name);
            return;
        };
        for (ind, hit) in hits.iter().enumerate() {
            if hit.accession != accession_id {
                continue;
            }
            let (start, end) =
                align_interval_u64(hit.subject_start, hit.subject_end);
            let (dec_s, dec_e) = decrement_lower_bound(start as f64, end as f64);
            let (acc_s, acc_e) = align_interval(dec_s, dec_e);
            let span = acc_e - acc_s;

            if span >= bin_size {
                individual_reads.push((name.to_string(), ind));
            } else {
                let mid = (acc_s + acc_e) / 2.0;
                let bin_idx = floor_with_min(mid / bin_size, 0) as usize;
                let bin_idx = bin_idx.min(num_bins.saturating_sub(1));
                read_bins[bin_idx].push((name.to_string(), ind));
            }
        }
    };

    for r in &accession_obj.reads {
        process_read(r);
    }
    for c in &accession_obj.contigs {
        process_contig(c);
    }

    // ── emit groups ──────────────────────────────────────────────────
    let mut hit_groups: Vec<Value> = Vec::new();

    for (name, ind) in &individual_reads {
        if let Some(hits) = read_data.get(name) {
            if let Some(h) = hits.get(*ind) {
                hit_groups.push(get_hit_group_json(&[], &[h], bin_size));
            }
        }
    }
    for (name, ind) in &individual_contigs {
        if let Some(hits) = contig_data.get(name) {
            if let Some(h) = hits.get(*ind) {
                hit_groups.push(get_hit_group_json(&[h], &[], bin_size));
            }
        }
    }

    for i in 0..num_bins {
        let read_refs: Vec<&ReadHit> = read_bins[i]
            .iter()
            .filter_map(|(n, ind)| read_data.get(n).and_then(|v| v.get(*ind)))
            .collect();
        let contig_refs: Vec<&ContigHit> = contig_bins[i]
            .iter()
            .filter_map(|(n, ind)| contig_data.get(n).and_then(|v| v.get(*ind)))
            .collect();

        if read_refs.is_empty() && contig_refs.is_empty() {
            continue;
        }
        hit_groups.push(get_hit_group_json(&contig_refs, &read_refs, bin_size));
    }

    json!(hit_groups)
}

/// One hit-group array
fn get_hit_group_json(
    contig_objs: &[&ContigHit],
    read_objs: &[&ReadHit],
    bin_size: f64,
) -> Value {
    let num_contigs = contig_objs.len();
    let num_reads = read_objs.len();
    let num_hits = num_contigs + num_reads;
    if num_hits == 0 {
        return json!([]);
    }

    // Unique byteranges → sum num_reads once per contig sequence
    let mut seen = HashSet::new();
    let mut contig_r: u64 = 0;
    let mut contig_byteranges: Vec<Value> = Vec::new();
    for c in contig_objs {
        let key = c.byterange.unwrap_or((0, 0));
        if seen.insert(key) {
            contig_r += c.num_reads;
            if let Some((off, len)) = c.byterange {
                contig_byteranges.push(json!([off, len]));
            }
        }
    }

    // Bounds + averages over all hits
    let mut endpoints: Vec<u64> = Vec::with_capacity(num_hits * 2);
    let mut sum_alen = 0.0;
    let mut sum_pident = 0.0;
    let mut sum_mismatch = 0.0;
    let mut sum_gaps = 0.0;

    for c in contig_objs {
        endpoints.push(c.subject_start);
        endpoints.push(c.subject_end);
        sum_alen += c.alignment_length as f64;
        sum_pident += c.percent_id;
        sum_mismatch += c.num_mismatches as f64;
        sum_gaps += c.num_gaps as f64;
    }
    for r in read_objs {
        endpoints.push(r.subject_start);
        endpoints.push(r.subject_end);
        sum_alen += r.alignment_length as f64;
        sum_pident += r.percent_id;
        sum_mismatch += r.num_mismatches as f64;
        sum_gaps += r.num_gaps as f64;
    }

    let hit_group_start = *endpoints.iter().min().unwrap_or(&0);
    let hit_group_end = *endpoints.iter().max().unwrap_or(&0);
    let mid = ((hit_group_start as f64 - 1.0) + hit_group_end as f64) / 2.0;
    let bin_index = floor_with_min(mid / bin_size, 0);

    let n = num_hits as f64;
    json!([
        num_contigs,
        num_reads,
        contig_r,
        hit_group_start,
        hit_group_end,
        format_number(sum_alen / n),
        format_percent(sum_pident / n / 100.0), //  avg(percent_id)/100
        format_number(sum_mismatch / n),
        format_number(sum_gaps / n),
        bin_index,
        contig_byteranges,
    ])
}

fn align_interval_u64(a: u64, b: u64) -> (u64, u64) {
    (a.min(b), a.max(b))
}



fn augment_contig_data_with_byteranges(
    contigs_fasta: &Path,
    contig_data: &mut HashMap<String, Vec<ContigHit>>,
) -> Result<()> {
    let file = File::open(contigs_fasta)
        .with_context(|| format!("open contigs_fasta {}", contigs_fasta.display()))?;
    let reader = BufReader::new(file);

    let mut seq_offset: u64 = 0;
    let mut seq_len: u64 = 0;
    let mut contig_name = String::new();

    for line in reader.lines() {
        let line = line?;
        // Re-add the newline that lines() strips — file offsets must match on-disk bytes.
        let line_bytes = (line.len() + 1) as u64; // +1 for '\n'

        if line.starts_with('>') {
            // Finalize previous record
            if seq_len > 0 {
                if let Some(hits) = contig_data.get_mut(&contig_name) {
                    for hit in hits.iter_mut() {
                        hit.byterange = Some((seq_offset, seq_len));
                    }
                }
                seq_offset += seq_len;
            }
            seq_len = line_bytes;
            contig_name = line[1..].trim().to_string();
        } else {
            seq_len += line_bytes;
        }
    }

    // Last contig in the file
    if seq_len > 0 {
        if let Some(hits) = contig_data.get_mut(&contig_name) {
            for hit in hits.iter_mut() {
                hit.byterange = Some((seq_offset, seq_len));
            }
        }
    }

    Ok(())
}

// Interval utilities — exact ports from Python
fn decrement_lower_bound(start: f64, end: f64) -> (f64, f64) {
    if start < end {
        (start - 1.0, end)
    } else {
        (start, end - 1.0)
    }
}

fn align_interval(a: f64, b: f64) -> (f64, f64) {
    (a.min(b), a.max(b))
}

fn transform_interval(
    interval: &[f64; 2],
    f_start: f64,
    f_end: f64,
    s_start: f64,
    s_end: f64,
) -> (f64, f64) {
    let pos1 = (interval[0] - f_start) / (f_end - f_start);
    let pos2 = (interval[1] - f_start) / (f_end - f_start);
    (
        pos1 * (s_end - s_start) + s_start,
        pos2 * (s_end - s_start) + s_start,
    )
}

fn floor_with_min(n: f64, min_v: i64) -> i64 {
    let f = n.floor() as i64;
    if f < min_v {
        min_v
    } else {
        f
    }
}

fn ceil_with_max(n: f64, max_v: i64) -> i64 {
    let c = n.ceil() as i64;
    if c > max_v {
        max_v
    } else {
        c
    }
}

fn format_number(n: f64) -> f64 {
    if n < 0.1 {
        (n * 1000.0).round() / 1000.0
    } else if n < 1.0 {
        (n * 100.0).round() / 100.0
    } else {
        (n * 10.0).round() / 10.0
    }
}

fn format_percent(n: f64) -> f64 {
    (n * 1000.0).round() / 1000.0
}
