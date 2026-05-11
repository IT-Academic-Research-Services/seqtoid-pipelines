use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use log::{self, warn};

use anyhow::{anyhow, Context, Result};
use futures::stream::{FuturesUnordered, StreamExt};
use serde::Serialize;
use serde_json::{json, Value};
use tokio::fs::create_dir_all;
use tokio_stream::wrappers::ReceiverStream;

use crate::utils::blast::M8Record;
use crate::utils::streams::{ParseOutput, ToBytes};

const MAX_NUM_BINS_COVERAGE: usize = 500;
const NUM_ACCESSIONS_PER_TAXON: usize = 10;
const MIN_CONTIG_SIZE: u64 = 500;


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
    subject_start: u64,
    subject_end: u64,
    query_start: u64,
    query_end: u64,
    total_length: u64,
    coverage: Vec<f64>,
    prop_mismatch: f64,
}

#[derive(Debug, Clone)]
struct ReadHit {
    subject_start: u64,
    subject_end: u64,
    prop_mismatch: f64,
    accession: String,
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
    refined_gsnap_blast_top_m8: ReceiverStream<ParseOutput>,
    contig_coverage_json: PathBuf,
    contig_stats_json: PathBuf,
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
        gsnap_deduped_m8,
        info_dict,
        min_size,
        acc_per_taxon,
        keep_taxons_with_no_contigs,
    )
        .await?;

    let coverage_viz_data = generate_coverage_viz_data(&accession_data, &contig_data, &read_data, max_bins).await?;

    let summary_data = generate_coverage_viz_summary_data(
        &taxon_data,
        &accession_data,
        &coverage_viz_data,
    );


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
    blast_top_m8: ReceiverStream<ParseOutput>,
    contig_coverage_json: &Path,
    contig_stats_json: &Path,
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

    let assigned_reads = get_unassigned_reads_set(&accession_data);  // actually assigned, per Python

    let mut contig_data = generate_contig_data(blast_top_m8, &valid_contigs).await?;

    let read_data = generate_read_data(gsnap_deduped_m8, &assigned_reads).await?;

    augment_contig_data_with_coverage(contig_coverage_json, &mut contig_data)?;

    let (taxon_data, accession_data) = select_best_accessions_per_taxon(
        taxon_data,
        accession_data,
        num_accessions_per_taxon,
    );

    Ok((taxon_data, accession_data, contig_data, read_data))
}

fn load_nt_info_db(path: &Path) -> Result<HashMap<String, (String, u64)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let acc = parts[0].to_string();
            let name = parts[1].to_string();
            let len: u64 = parts[2].parse()?;
            map.insert(acc, (name, len));
        }
    }
    Ok(map)
}

fn get_valid_contigs_with_read_counts(
    path: &Path,
    min_size: u64,
) -> Result<HashMap<String, u64>> {
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

    let mut line_count = 0;

    while let Some(item) = hit_summary.next().await {
        line_count += 1;
        if line_count % 100_000 == 0 {
            log::info!("Processed {} hit_summary lines", line_count);
        }

        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() >= 12 && valid_contigs.contains_key(fields[7]) {
            let taxon = fields[9].to_string();
            let acc = fields[8].to_string();
            let contig = fields[7].to_string();

            taxon_data.entry(taxon).or_default().accessions.insert(acc.clone());
            accession_data.entry(acc).or_default().contigs.push(contig);
        } else if fields.len() >= 5 {
            let taxon = fields[4].to_string();
            let acc = fields[3].to_string();
            let read = fields[0].to_string();

            taxon_data.entry(taxon).or_default().accessions.insert(acc.clone());
            accession_data.entry(acc).or_default().reads.push(read);
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
    mut blast_top_m8: ReceiverStream<ParseOutput>,
    valid_contigs: &HashMap<String, u64>,
) -> Result<HashMap<String, Vec<ContigHit>>> {
    let mut contig_data: HashMap<String, Vec<ContigHit>> = HashMap::new();

    let mut line_count = 0;

    while let Some(item) = blast_top_m8.next().await {
        line_count += 1;
        if line_count % 100_000 == 0 {
            log::info!("Processed {} blast_top_m8 lines", line_count);
        }

        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let line_trim = line.trim_end();
        if line_trim.is_empty() {
            continue;
        }

        let m8 = match M8Record::parse_line_nt(line_trim) {
            Ok(record) => record,
            Err(e) => {
                warn!("Failed to parse m8 line: {} - {}", line_trim, e);
                continue;
            }
        };

        let contig_name = m8.qname.clone();
        if !valid_contigs.contains_key(&contig_name) {
            continue;
        }

        let total_length = *valid_contigs.get(&contig_name).unwrap_or(&0);
        let prop_mismatch = m8.mismatch as f64 / m8.alen as f64;

        let hit = ContigHit {
            subject_start: m8.tstart,
            subject_end: m8.tend,
            query_start: m8.qstart,
            query_end: m8.qend,
            total_length,
            coverage: vec![],
            prop_mismatch,
        };

        contig_data.entry(contig_name).or_default().push(hit);
    }

    Ok(contig_data)
}

async fn generate_read_data(
    mut gsnap_deduped_m8: ReceiverStream<ParseOutput>,
    assigned_reads: &HashSet<String>,
) -> Result<HashMap<String, Vec<ReadHit>>> {
    let mut read_data: HashMap<String, Vec<ReadHit>> = HashMap::new();

    let mut line_count = 0;

    while let Some(item) = gsnap_deduped_m8.next().await {
        line_count += 1;
        if line_count % 100_000 == 0 {
            log::info!("Processed {} gsnap_deduped_m8 lines", line_count);
        }

        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let line_trim = line.trim_end();
        if line_trim.is_empty() {
            continue;
        }

        let m8 = match M8Record::parse_line_nt(line_trim) {
            Ok(record) => record,
            Err(e) => {
                warn!("Failed to parse m8 line: {} - {}", line_trim, e);
                continue;
            }
        };

        let read_name = m8.qname.clone();
        if assigned_reads.contains(&read_name) {
            continue;
        }

        let accession_id = m8.tname.clone();
        let prop_mismatch = m8.mismatch as f64 / m8.alen as f64;

        let hit = ReadHit {
            subject_start: m8.tstart,
            subject_end: m8.tend,
            prop_mismatch,
            accession: accession_id,
        };

        read_data.entry(read_name).or_default().push(hit);
    }

    Ok(read_data)
}

fn augment_contig_data_with_coverage(
    path: &Path,
    contig_data: &mut HashMap<String, Vec<ContigHit>>,
) -> Result<()> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let contig_coverage: HashMap<String, Vec<f64>> = serde_json::from_reader(&mut reader)
        .map_err(|e| anyhow!("Failed to parse contig_coverage_json: {}", e))?;

    for (contig_name, hits) in contig_data.iter_mut() {
        if let Some(cov) = contig_coverage.get(contig_name) {
            for hit in hits.iter_mut() {
                hit.coverage = cov.clone();
            }
        } else {
            warn!("No coverage data for contig: {}", contig_name);
        }
    }

    Ok(())
}

fn select_best_accessions_per_taxon(
    mut taxon_data: HashMap<String, TaxonData>,
    mut accession_data: HashMap<String, AccessionData>,
    num_per_taxon: usize,
) -> (HashMap<String, TaxonData>, HashMap<String, AccessionData>) {
    for ad in accession_data.values_mut() {
        ad.score = ad.contigs.len() as f64 * 1000.0 + ad.reads.len() as f64;
    }

    for td in taxon_data.values_mut() {
        let mut scored: Vec<_> = td
            .accessions
            .iter()
            .filter_map(|acc| accession_data.get(acc).map(|ad| (acc.clone(), ad.score)))
            .collect();
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        td.best_accessions = scored.into_iter().take(num_per_taxon).map(|(acc, _)| acc).collect();
    }

    (taxon_data, accession_data)
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

        tasks.push(tokio::task::spawn_blocking(move || {
            let total_len = acc_obj.total_length as f64;
            if total_len == 0.0 {
                return Err(anyhow!("Zero length accession: {}", acc_id));
            }

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

            let viz = CoverageVizData {
                total_length: acc_obj.total_length,
                name: acc_obj.name.clone(),
                hit_groups: json!([]),
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
        let task_result: Result<Result<(String, CoverageVizData)>, tokio::task::JoinError> = task_result;
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
                let (s_start, s_end) = decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
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
                    let cov_end = ceil_with_max(coverage_interval.1, hit.coverage.len() as i64) as usize;

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
                let (s_start, s_end) = decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
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
    let mut hit_count = 0usize;
    let mut endpoints = Vec::new();

    for contig in &acc_obj.contigs {
        if let Some(hits) = contig_data.get(contig) {
            for hit in hits {
                let (dec_start, dec_end) = decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
                let (s_start, s_end) = align_interval(dec_start, dec_end);
                let aligned_len = (s_end - s_start) as u64;
                max_len = max_len.max(aligned_len);

                let (q_dec_start, q_dec_end) = decrement_lower_bound(hit.query_start as f64, hit.query_end as f64);
                let (q_start, q_end) = align_interval(q_dec_start, q_dec_end);
                cov_sum += hit.coverage[q_start as usize..q_end as usize].iter().sum::<f64>();

                mismatch_sum += hit.prop_mismatch;
                hit_count += 1;

                endpoints.push((s_start, 1));
                endpoints.push((s_end, -1));
            }
        }
    }

    for read in &acc_obj.reads {
        if let Some(hits) = read_data.get(read) {
            for hit in hits {
                let (dec_start, dec_end) = decrement_lower_bound(hit.subject_start as f64, hit.subject_end as f64);
                let (s_start, s_end) = align_interval(dec_start, dec_end);
                let aligned_len = (s_end - s_start) as u64;
                max_len = max_len.max(aligned_len);

                cov_sum += aligned_len as f64;

                mismatch_sum += hit.prop_mismatch;
                hit_count += 1;

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

    Ok(AccessionStats {
        max_aligned_length: max_len,
        coverage_depth: if total_len > 0.0 { cov_sum / total_len } else { 0.0 },
        coverage_breadth: breadth,
        avg_prop_mismatch: if hit_count > 0 { mismatch_sum / hit_count as f64 } else { 0.0 },
    })
}

fn generate_coverage_viz_summary_data(
    taxon_data: &HashMap<String, TaxonData>,
    accession_data: &HashMap<String, AccessionData>,
    coverage_viz: &HashMap<String, CoverageVizData>,
) -> CoverageVizSummary {
    let mut taxons = HashMap::new();

    for (taxon_id, td) in taxon_data {
        let best_acc: Vec<AccessionSummary> = td.best_accessions.iter().map(|acc_id| {
            let ad = &accession_data[acc_id];
            let cv = &coverage_viz[acc_id];
            AccessionSummary {
                id: acc_id.clone(),
                name: ad.name.clone(),
                num_contigs: ad.contigs.len(),
                num_reads: ad.reads.len(),
                score: format_number(ad.score),
                coverage_breadth: cv.coverage_breadth,
                coverage_depth: cv.coverage_depth,
            }
        }).collect();

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
    if f < min_v { min_v } else { f }
}

fn ceil_with_max(n: f64, max_v: i64) -> i64 {
    let c = n.ceil() as i64;
    if c > max_v { max_v } else { c }
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