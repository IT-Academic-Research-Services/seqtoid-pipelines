//! AMR Pipeline

use crate::pipelines::short_read_mngs::HostFilterOutputs;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;
use serde::ser::{Serialize, SerializeMap, Serializer};

use anyhow::anyhow;
use log::{error, info, warn};
use tempfile::TempDir;
use bytes::Bytes;
use dashmap::DashMap;
use tokio::sync::{mpsc, oneshot};
use tokio::task::JoinHandle;
use tokio::io::AsyncBufReadExt;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;

use crate::config::defs::{ClusterInfo, PairingMode, PipelineError, RunConfig, StreamDataType,
                          SPADES_TAG, KMA_TAG, SAMTOOLS_TAG, BAMTOOLS_TAG, BEDTOOLS_TAG, PIGZ_TAG,
                          BedtoolsSubcommand, SamtoolsSubcommand, BamtoolsSubcommand};
use crate::pipelines::short_read_mngs::{
    kallisto_results,
    run_host_filter,
    HostFilterCompletion,
};
use crate::utils::command::kma::KmaConfig;
use crate::utils::command::pigz::PigzConfig;
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::bamtools::BamtoolsConfig;
use crate::utils::command::bedtools::BedtoolsConfig;
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::generate_cli;
use crate::utils::fastx::{parse_byte_stream_to_fastq, read_fasta, read_fastq, write_fasta_stream_to_file, SequenceRecord};
use crate::utils::file::{choose_temp_dir, resolve_existing_input_path, validate_file_inputs, write_count_file, write_byte_stream_to_file};
use crate::utils::sambam::parse_sam_record;
use crate::utils::sequence::translate_rgi_kma_consensus;
use crate::utils::streams::{
    bytes_to_lines,
    fanout_to_channels,
    interleave_fastq_streams,
    join_with_error_handling,
    parse_bytes,
    parse_child_output,
    spawn_cmd,
    stream_to_file,
    stream_to_cmd,
    ChildStream,
    ParseMode,
    ParseOutput,
};
use crate::utils::subsample::load_duplicate_clusters_csv;

const RGI_MIN_CONSENSUS_DEPTH: f64 = 5.0;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AmrInputMode {
    Raw,
    Filtered,
}

impl AmrInputMode {
    fn from_config(config: &RunConfig) -> Result<Self, PipelineError> {
        match config.args.amr_input_mode.as_str() {
            "raw" => Ok(Self::Raw),
            "filtered" => Ok(Self::Filtered),
            other => Err(PipelineError::InvalidConfig(format!(
                "Invalid AMR input mode '{}'; expected 'raw' or 'filtered'",
                other
            ))),
        }
    }
}

struct RunRedupHandle {
    /// One live FASTA stream per input end, in WDL `redups_1.fa`,
    /// `redups_2.fa` order.
    redups_streams: Vec<ReceiverStream<ParseOutput>>,

    /// Completes after every per-end producer has reached EOF.
    /// Callers must install consumers for `redups_streams` before awaiting it.
    completion_task: JoinHandle<Result<(), anyhow::Error>>,
}

#[derive(Debug)]
struct RunKmaHandle {
    /// Complete line-oriented KMA SAM for the native AMR parser.
    sam_stream: ReceiverStream<ParseOutput>,

    /// Raw KMA SAM bytes, unchanged, for direct streaming into samtools fixmate.
    samtools_stream: mpsc::Receiver<ParseOutput>,

    /// Completes when every raw KMA stdout byte has reached both branches.
    sam_fanout_done: oneshot::Receiver<Result<(), anyhow::Error>>,

    /// Completes when the parser branch has converted raw byte chunks to lines.
    sam_lines_task: JoinHandle<Result<(), anyhow::Error>>,

    completion_task: JoinHandle<Result<(), anyhow::Error>>,

    // Keep KMA sidecars alive for downstream AMR post-processing.
    output_temp_dir: TempDir,
    output_prefix: PathBuf,
}

#[derive(Debug, Clone)]
struct KmaReference {
    name: String,
    length: u64,
}

#[derive(Debug, Clone)]
struct KmaAlignment {
    qname: String,
    flag: u16,
    rname: String,
    pos: u64,
    mapq: u8,
    cigar: String,
    rnext: String,
    seq: String,
}

#[derive(Debug)]
struct KmaSamResults {
    references: HashMap<String, KmaReference>,
    alignments: Vec<KmaAlignment>,

    /// Number of non-header SAM records emitted by KMA.
    total_records: u64,

    /// KMA records that are not usable alignments despite potentially
    /// misleading FLAG values, e.g. RNAME="*", POS=0, CIGAR="*".
    invalid_alignment_records: u64,
}

#[derive(Debug, Clone, serde::Serialize)]
struct RgiReadCounts {
    reference_length: u64,
    mapped: u64,
    unmapped: u64,
    all: u64,
}

/// Per-reference breadth of coverage from the BEDTools genomecov histogram.
/// This is not KMA's read depth, which is consumed separately from its .res file.
#[derive(Debug, Clone, Default)]
struct RgiCoverage {
    reference_length: u64,
    covered_bases: u64,
    uncovered_bases: u64,
    percent_coverage: f64,
}

#[derive(Debug, Clone, serde::Serialize)]
struct RgiAlignment {
    qname: String,
    flag: u16,
    rname: String,
    pos: u64,
    mapq: u8,
    mrnm: String,
}

#[derive(Debug, Clone)]
struct CardModel {
    model_id: String,
    aro_accession: String,
    model_name: String,
    model_type: String,

    // BWT.py:
    // category_aro_class_name -> unique category_aro_name values
    categories: HashMap<String, Vec<String>>,

    // Unique first-two-word NCBI taxonomy names, preserving encounter order.
    taxon: Vec<String>,
}

#[derive(Debug)]
struct CardModels {
    by_model_id: HashMap<String, CardModel>,
    by_aro_accession: HashMap<String, CardModel>,
}


#[derive(Debug, Clone, serde::Serialize)]
struct RgiPercentCoverage {
    covered: String,
    uncovered: String,
}

#[derive(Debug, Clone, serde::Serialize)]
struct RgiLengthCoverage {
    covered: String,
    uncovered: String,
}

#[derive(Debug, Clone, serde::Serialize)]
struct RgiReference {
    sequence_length: String,
}

#[derive(Debug, Clone, serde::Serialize)]
struct RgiAlleleSummary {
    id: String,
    cvterm_name: String,
    aro_accession: String,
    model_type: String,
    database: String,
    reference_allele_source: String,
    observed_in_genomes: String,
    observed_in_plasmids: String,
    observed_in_pathogens: Vec<String>,
    range_of_reference_allele_source: f64,
    reads: RgiReadCounts,
    alignments: Vec<RgiAlignment>,
    mapq_average: String,
    number_of_mapped_baits: u64,
    number_of_mapped_baits_with_reads: u64,
    average_bait_coverage: u64,
    bait_coverage_coefficient_of_variation: u64,
    mate_pair: Vec<String>,
    percent_coverage: RgiPercentCoverage,
    length_coverage: RgiLengthCoverage,
    reference: RgiReference,
    mutation: String,
    resistomes: HashMap<String, Vec<String>>,
    predicted_pathogen: String,
    depth: String,
    snps: String,
    consensus_sequence_dna: String,
    consensus_sequence_protein: String,
}

/// One result from BWT.py summary(). Python has three observable outcomes:
/// a populated dict, `{}` for KMA consensus-block failures, or implicit None
/// for failures caught by the outer summary() exception handler.
#[derive(Debug, Clone)]
enum RgiAlleleSummaryResult {
    Summary(RgiAlleleSummary),
    Empty,
    Null,
}

impl Serialize for RgiAlleleSummaryResult {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            Self::Summary(summary) => summary.serialize(serializer),
            Self::Empty => {
                let map = serializer.serialize_map(Some(0))?;
                map.end()
            }
            Self::Null => serializer.serialize_none(),
        }
    }
}

async fn write_rgi_allele_json(
    results: &[RgiAlleleSummaryResult],
    out_dir: &std::path::Path,
) -> anyhow::Result<PathBuf> {
    let path = out_dir.join("sr_amr_report.allele_mapping_data.json");
    let json = serde_json::to_vec(results)?;
    tokio::fs::write(&path, json).await?;
    info!(
        "[AMR] RGI allele JSON written: results={} path={}",
        results.len(),
        path.display()
    );
    Ok(path)
}

fn build_rgi_allele_summaries(
    reads: &HashMap<String, RgiReadCounts>,
    coverage: &HashMap<String, RgiCoverage>,
    alignments: &HashMap<String, Vec<RgiAlignment>>,
    sidecars: &KmaSidecars,
    models: &CardModels,
) -> Vec<RgiAlleleSummaryResult> {
    // BWT.py's ThreadPoolExecutor makes final JSON array order scheduling-
    // dependent. Keep Rust deterministic; parity comparisons should compare
    // contents rather than array position.
    let mut references: Vec<&String> = reads.keys().collect();
    references.sort();

    let mut results = Vec::with_capacity(references.len());
    let mut summaries = 0usize;
    let mut empty_results = 0usize;
    let mut null_results = 0usize;
    let mut skipped_no_consensus = 0usize;
    let mut skipped_translation_error = 0usize;

    for alignment_hit in references {
        let Some(model) = resolve_card_model(alignment_hit, models) else {
            // get_model_id() can yield an empty model_id; models[model_id]
            // then raises inside summary()'s outer try and Python returns None.
            warn!(
                "[AMR] missing model for alignment {}; preserving RGI null summary result",
                alignment_hit
            );
            results.push(RgiAlleleSummaryResult::Null);
            null_results += 1;
            continue;
        };

        // Everything below corresponds to BWT.py's outer summary() try.
        // Missing coverage/model-derived evidence therefore becomes Python
        // None rather than a fatal pipeline error.
        let Some(cov) = coverage.get(alignment_hit) else {
            warn!(
                "[AMR] missing genomecov evidence for {}; preserving RGI null summary result",
                alignment_hit
            );
            results.push(RgiAlleleSummaryResult::Null);
            null_results += 1;
            continue;
        };

        let counts = reads
            .get(alignment_hit)
            .expect("reference came from reads keys");
        let hit_alignments = alignments.get(alignment_hit).cloned().unwrap_or_default();

        let mapq_average = if hit_alignments.is_empty() {
            0.0
        } else {
            hit_alignments.iter().map(|a| a.mapq as f64).sum::<f64>()
                / hit_alignments.len() as f64
        };

        // BWT.py keeps unique MRNM values in encounter order, excluding '=' and
        // mates whose reference string contains this model's ARO accession.
        let aro_marker = format!("ARO:{}", model.aro_accession);
        let mut mate_pair = Vec::<String>::new();
        for alignment in &hit_alignments {
            if alignment.mrnm != "="
                && !alignment.mrnm.contains(&aro_marker)
                && !mate_pair.iter().any(|m| m == &alignment.mrnm)
            {
                mate_pair.push(alignment.mrnm.clone());
            }
        }

        // This pipeline currently uses canonical CARD references, not the
        // optional Resistomes & Variants / wildcard input handled by BWT.py.
        let database = "CARD".to_string();
        let reference_allele_source = "CARD curation".to_string();
        let observed_in_genomes = "no data".to_string();
        let observed_in_plasmids = "no data".to_string();
        let observed_in_pathogens = model.taxon.clone();
        let percent_identity = 100.0;

        let mut resistomes = model.categories.clone();
        resistomes.entry("AMR Gene Family".to_string()).or_default();
        resistomes.entry("Drug Class".to_string()).or_default();
        resistomes.entry("Resistance Mechanism".to_string()).or_default();

        // BWT.py's KMA block is its own inner try/except. Any exception here
        // returns {}, which is observably different from the outer None path.
        let Some(consensus_sequence_dna) = sidecars.consensus_dna.get(alignment_hit) else {
            skipped_no_consensus += 1;
            empty_results += 1;
            warn!(
                "[AMR] model with id: {}, has few mapped reads to make consensus sequence skipping: missing consensus for {}",
                model.model_id,
                alignment_hit
            );
            results.push(RgiAlleleSummaryResult::Empty);
            continue;
        };

        let mut consensus_sequence_protein =
            match translate_rgi_kma_consensus(consensus_sequence_dna) {
                Ok(protein) => protein,
                Err(e) => {
                    skipped_translation_error += 1;
                    empty_results += 1;
                    warn!(
                        "[AMR] model with id: {}, has few mapped reads to make consensus sequence skipping: {}",
                        model.model_id,
                        e
                    );
                    results.push(RgiAlleleSummaryResult::Empty);
                    continue;
                }
            };

        let depth = sidecars
            .read_depth
            .get(alignment_hit)
            .cloned()
            .unwrap_or_default();
        let mut snps = sidecars
            .mutations
            .get(alignment_hit)
            .map(|items| items.join("; "))
            .unwrap_or_default();
        let mut consensus_sequence_dna = consensus_sequence_dna.clone();

        // BWT.py performs float(read_coverage_depth.strip() or 0) outside the
        // inner KMA try but inside the outer summary() try. A malformed nonempty
        // depth therefore yields Python None, not {} and not a pipeline failure.
        let parsed_depth = if depth.trim().is_empty() {
            0.0
        } else {
            match depth.trim().parse::<f64>() {
                Ok(value) => value,
                Err(e) => {
                    warn!(
                        "[AMR] invalid KMA depth {:?} for model {} / {}; preserving RGI null summary result: {}",
                        depth,
                        model.model_id,
                        alignment_hit,
                        e
                    );
                    results.push(RgiAlleleSummaryResult::Null);
                    null_results += 1;
                    continue;
                }
            }
        };

        if parsed_depth < RGI_MIN_CONSENSUS_DEPTH {
            snps.clear();
            consensus_sequence_dna.clear();
            consensus_sequence_protein.clear();
        }

        results.push(RgiAlleleSummaryResult::Summary(RgiAlleleSummary {
            id: alignment_hit.clone(),
            cvterm_name: model.model_name.clone(),
            aro_accession: model.aro_accession.clone(),
            model_type: model.model_type.clone(),
            database,
            reference_allele_source,
            observed_in_genomes,
            observed_in_plasmids,
            observed_in_pathogens,
            range_of_reference_allele_source: percent_identity,
            reads: counts.clone(),
            alignments: hit_alignments,
            mapq_average: format!("{:.2}", mapq_average),
            // include_baits is false in the WDL invocation; BWT.py returns four zeros.
            number_of_mapped_baits: 0,
            number_of_mapped_baits_with_reads: 0,
            average_bait_coverage: 0,
            bait_coverage_coefficient_of_variation: 0,
            mate_pair,
            percent_coverage: RgiPercentCoverage {
                covered: format!("{:.2}", cov.percent_coverage),
                uncovered: format!("{:.2}", 100.0 - cov.percent_coverage),
            },
            length_coverage: RgiLengthCoverage {
                covered: cov.covered_bases.to_string(),
                uncovered: cov.uncovered_bases.to_string(),
            },
            reference: RgiReference {
                sequence_length: cov.reference_length.to_string(),
            },
            mutation: "N/A".to_string(),
            resistomes,
            predicted_pathogen: "N/A".to_string(),
            depth,
            snps,
            consensus_sequence_dna,
            consensus_sequence_protein,
        }));
        summaries += 1;
    }

    info!(
        "[AMR] RGI allele summaries built: reportable_references={} results={} summaries={} empty_results={} null_results={} skipped_no_consensus={} skipped_translation_error={}",
        reads.len(),
        results.len(),
        summaries,
        empty_results,
        null_results,
        skipped_no_consensus,
        skipped_translation_error
    );

    results
}

/// Native streaming equivalent of the WDL RunRedup task.
///
/// WDL semantics:
///   1. `subsampled_reads` form the prefix of each `redups_N.fa`.
///   2. Only a duplicate cluster whose representative survived subsampling is
///      eligible for re-expansion.
///   3. The representative itself is already present in the subsampled FASTA;
///      append only its non-representative cluster members.
///   4. Appended records are selected from the corresponding host-filtered,
///      pre-dedup FASTQ and converted to FASTA.
///
/// Rust keeps the cluster relationship in memory, so the WDL-only
/// `cluster_sizes` / `clusters.csv` parsing and temporary ID files are not
/// reproduced here.
async fn run_redup(
    config: Arc<RunConfig>,
    mut host_filtered_stream: ReceiverStream<ParseOutput>,
    mut subsampled_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    duplicate_clusters_task: JoinHandle<
        Result<Arc<DashMap<String, ClusterInfo>>, anyhow::Error>
    >,
) -> Result<RunRedupHandle, PipelineError> {
    let part_count = if paired { 2 } else { 1 };

    let mut output_streams = Vec::with_capacity(part_count);
    let mut output_txs = Vec::with_capacity(part_count);

    for _ in 0..part_count {
        let (tx, rx) =
            mpsc::channel::<ParseOutput>(config.base_buffer_size);

        output_txs.push(tx);
        output_streams.push(ReceiverStream::new(rx));
    }

    // Drain the pre-dedup stream immediately. It is upstream of dedup and
    // therefore cannot be left unread while dedup/subsampling finalize.
    //
    // Keep records separated by end so duplicate members can later be
    // appended to the corresponding redups_N.fa output.
    let host_buffer_task = tokio::spawn(async move {
        let mut records_by_part =
            (0..part_count)
                .map(|_| Vec::<SequenceRecord>::new())
                .collect::<Vec<_>>();

        let mut record_index = 0usize;

        while let Some(item) = host_filtered_stream.next().await {
            let ParseOutput::Fastq(record) = item else {
                return Err(anyhow!(
                    "RunRedup expected FASTQ on interleaved host-filtered stream"
                ));
            };

            let part_index = if paired {
                record_index % 2
            } else {
                0
            };

            records_by_part[part_index].push(record);
            record_index += 1;
        }

        if paired && record_index % 2 != 0 {
            return Err(anyhow!(
                "RunRedup received incomplete paired pre-dedup stream: {} records",
                record_index
            ));
        }

        Ok::<Vec<Vec<SequenceRecord>>, anyhow::Error>(records_by_part)
    });

    // Forward the complete subsampled input immediately as FASTA while
    // recording which deduplicated representatives survived subsampling.
    let prefix_task = tokio::spawn(async move {
        let mut passed_representatives = HashSet::<String>::new();
        let mut record_index = 0usize;

        while let Some(item) = subsampled_stream.next().await {
            let ParseOutput::Fastq(record) = item else {
                return Err(anyhow!(
                    "RunRedup expected FASTQ on interleaved subsampled stream"
                ));
            };

            let part_index = if paired {
                record_index % 2
            } else {
                0
            };

            passed_representatives.insert(record.id().to_string());

            let fasta = SequenceRecord::Fasta {
                id: record.id().to_string(),
                desc: None,
                seq: record.seq(),
            };

            output_txs[part_index]
                .send(ParseOutput::Fasta(fasta))
                .await
                .map_err(|_| {
                    anyhow!(
                        "RunRedup downstream receiver dropped for part {}",
                        part_index + 1
                    )
                })?;

            record_index += 1;
        }

        if paired && record_index % 2 != 0 {
            return Err(anyhow!(
                "RunRedup received incomplete paired subsampled stream: {} records",
                record_index
            ));
        }

        Ok::<_, anyhow::Error>((
            passed_representatives,
            output_txs,
        ))
    });

    let completion_task = tokio::spawn(async move {
        // Both live inputs already have consumers installed above.
        let duplicate_clusters = duplicate_clusters_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunRedup duplicate-cluster task panicked: {}",
                    e
                )
            })??;

        let (passed_representatives, output_txs) = prefix_task
            .await
            .map_err(|e| {
                anyhow!("RunRedup prefix task panicked: {}", e)
            })??;

        let host_records_by_part = host_buffer_task
            .await
            .map_err(|e| {
                anyhow!("RunRedup host-buffer task panicked: {}", e)
            })??;

        let mut duplicate_member_ids =
            HashSet::<String>::new();

        for entry in duplicate_clusters.iter() {
            let rep_id = entry.key();
            let info = entry.value();

            if info.size <= 1
                || !passed_representatives.contains(rep_id.as_str())
            {
                continue;
            }

            for member_id in &info.members {
                if member_id != rep_id {
                    duplicate_member_ids.insert(member_id.clone());
                }
            }
        }

        for (part_index, host_records) in
            host_records_by_part.into_iter().enumerate()
        {
            let tx = &output_txs[part_index];

            // Preserve original pre-dedup order within each end.
            for record in host_records {
                if !duplicate_member_ids.contains(record.id()) {
                    continue;
                }

                let fasta = SequenceRecord::Fasta {
                    id: record.id().to_string(),
                    desc: None,
                    seq: record.seq(),
                };

                tx.send(ParseOutput::Fasta(fasta))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "RunRedup downstream receiver dropped for part {}",
                            part_index + 1
                        )
                    })?;
            }
        }

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunRedupHandle {
        redups_streams: output_streams,
        completion_task,
    })
}

#[derive(Debug)]
struct RunSpadesHandle {
    contigs_path: PathBuf,
    scaffolds_path: PathBuf,
    completion_task: JoinHandle<Result<(), anyhow::Error>>,
}

/// Native AMR equivalent of the WDL RunSpades task.
///
/// RunRedup produces live FASTA streams, but SPAdes inspects and reopens its
/// read inputs during a run. Materializing the inputs is therefore a genuine
/// external-tool boundary rather than a pipeline checkpoint. The per-end
/// writers are installed before this function returns so RunRedup is never
/// left blocked on an unconsumed bounded output.
async fn run_spades_amr(
    config: Arc<RunConfig>,
    reduplicated_streams: Vec<ReceiverStream<ParseOutput>>,
) -> Result<RunSpadesHandle, PipelineError> {
    if reduplicated_streams.is_empty() || reduplicated_streams.len() > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunSpades expected one or two reduplicated FASTA streams, got {}",
            reduplicated_streams.len()
        )));
    }

    let input_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        2,
        false,
    )
        .await?;

    let spades_out_dir = config.out_dir.join("spades");
    let contigs_path = spades_out_dir.join("contigs.fasta");
    let scaffolds_path = spades_out_dir.join("scaffolds.fasta");

    let mut input_paths = Vec::with_capacity(reduplicated_streams.len());
    let mut input_write_tasks = Vec::with_capacity(reduplicated_streams.len());

    for (idx, stream) in reduplicated_streams.into_iter().enumerate() {
        let path = input_temp_dir
            .path()
            .join(format!("redups_{}.fa", idx + 1));

        input_write_tasks.push(write_fasta_stream_to_file(
            stream,
            path.clone(),
            config.clone(),
            StreamDataType::JustBytes,
            &format!("AMR RunSpades input part {}", idx + 1),
        ));
        input_paths.push(path);
    }

    let completion_config = config.clone();
    let completion_contigs = contigs_path.clone();
    let completion_scaffolds = scaffolds_path.clone();

    let completion_task = tokio::spawn(async move {
        // Keep the scratch directory alive for the entire SPAdes invocation.
        let _input_temp_dir = input_temp_dir;

        for (idx, task) in input_write_tasks.into_iter().enumerate() {
            task.await
                .map_err(|e| anyhow!(
                    "RunSpades input writer {} panicked: {}",
                    idx + 1,
                    e
                ))??;
        }

        // A RunRedup/input-adapter failure is an upstream pipeline failure.
        // From this point onward, match WDL RunSpades: assembly failures are
        // non-fatal and produce compatibility marker files.
        async fn write_failure_outputs(
            contigs_path: &PathBuf,
            scaffolds_path: &PathBuf,
        ) -> Result<(), anyhow::Error> {
            if let Some(parent) = contigs_path.parent() {
                tokio::fs::create_dir_all(parent).await?;
            }
            tokio::fs::write(contigs_path, b";ASSEMBLY FAILED\n").await?;
            tokio::fs::write(scaffolds_path, b";ASSEMBLY FAILED\n").await?;
            Ok(())
        }

        if spades_out_dir.exists() {
            tokio::fs::remove_dir_all(&spades_out_dir).await?;
        }

        let mut options = HashMap::new();
        options.insert("--only-assembler".to_string(), None);

        let spades_config = SpadesConfig {
            r1_path: input_paths[0].clone(),
            r2_path_opt: input_paths.get(1).cloned(),
            outdir_path: spades_out_dir.clone(),
            tempdir_path: Some(_input_temp_dir.path().to_path_buf()),
            option_fields: options,
        };

        let spades_args = match generate_cli(
            SPADES_TAG,
            &completion_config,
            Some(&spades_config),
        ) {
            Ok(args) => args,
            Err(e) => {
                warn!("[AMR] RunSpades command generation failed (non-fatal): {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        };

        info!("[AMR] RunSpades args: {:?}", spades_args);

        let (mut child, stderr_task) = match spawn_cmd(
            completion_config.clone(),
            SPADES_TAG,
            spades_args,
            completion_config.args.verbose,
            None,
        )
            .await
        {
            Ok(v) => v,
            Err(e) => {
                warn!("[AMR] RunSpades spawn failed (non-fatal): {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        };

        // Drain stdout concurrently so a verbose SPAdes process cannot block on
        // its pipe. Preserve it as a useful sidecar log.
        let stdout_task = match parse_child_output(
            &mut child,
            ChildStream::Stdout,
            ParseMode::Lines,
            &completion_config,
        )
            .await
        {
            Ok(stdout_rx) => Some(tokio::spawn(stream_to_file(
                stdout_rx,
                completion_config.out_dir.join("amr_spades_stdout.log"),
            ))),
            Err(e) => {
                warn!("[AMR] RunSpades stdout capture failed (non-fatal): {}", e);
                None
            }
        };

        let status = child.wait().await;
        let stderr_result = stderr_task.await;
        if let Some(task) = stdout_task {
            let _ = task.await;
        }

        let spades_ok = matches!(status, Ok(ref s) if s.success())
            && matches!(stderr_result, Ok(Ok(())));

        if !spades_ok || !completion_contigs.is_file() {
            warn!(
                "[AMR] RunSpades assembly failed or contigs.fasta is missing; writing WDL-compatible marker"
            );
            write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
            return Ok(());
        }

        let metadata = tokio::fs::metadata(&completion_contigs).await?;
        if metadata.len() == 0 {
            warn!("[AMR] RunSpades produced empty contigs.fasta; writing failure marker");
            write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
            return Ok(());
        }

        // Match `seqtk seq -L min_contig_length`: read the completed SPAdes
        // FASTA and rewrite contigs.fasta with records at or above the configured
        // minimum length. Use a separate path to avoid truncating the input.
        let unfiltered_contigs = spades_out_dir.join("contigs_unfiltered.fasta");
        if let Err(e) = tokio::fs::rename(&completion_contigs, &unfiltered_contigs).await {
            warn!("[AMR] RunSpades could not stage contigs for filtering: {}", e);
            write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
            return Ok(());
        }

        let filtered_rx = match read_fasta(
            unfiltered_contigs.clone(),
            u64::MAX,
            Some(completion_config.args.min_contig_length),
            None,
            &completion_config,
        ) {
            Ok(rx) => rx,
            Err(e) => {
                warn!("[AMR] RunSpades contig filtering read failed: {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        };

        let filter_write_task = write_fasta_stream_to_file(
            ReceiverStream::new(filtered_rx),
            completion_contigs.clone(),
            completion_config.clone(),
            StreamDataType::JustBytes,
            "AMR RunSpades contig filtering",
        );

        match filter_write_task.await {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                warn!("[AMR] RunSpades contig filtering write failed: {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
            Err(e) => {
                warn!("[AMR] RunSpades contig filtering writer panicked: {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        }

        let _ = tokio::fs::remove_file(unfiltered_contigs).await;

        info!(
            "[AMR] RunSpades complete: contigs={}, scaffolds={}",
            completion_contigs.display(),
            completion_scaffolds.display()
        );

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunSpadesHandle {
        contigs_path,
        scaffolds_path,
        completion_task,
    })
}




async fn fastq_file_to_record_stream(
    path: PathBuf,
    config: Arc<RunConfig>,
    tag: &str,
) -> Result<(ReceiverStream<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>), PipelineError> {
    let (byte_rx, read_task) = read_fastq(
        path,
        None,
        PairingMode::Strict,
        None,
        u64::MAX,
        None,
        None,
        tag,
        &config,
    )
        .map_err(PipelineError::Other)?;

    let (record_rx, parse_task) = parse_byte_stream_to_fastq(
        byte_rx,
        config.base_buffer_size,
        config.args.stall_threshold,
    )
        .await
        .map_err(PipelineError::Other)?;

    let completion_task = tokio::spawn(async move {
        read_task
            .await
            .map_err(|e| anyhow!("FASTQ file reader task panicked: {}", e))??;
        parse_task
            .await
            .map_err(|e| anyhow!("FASTQ byte parser task panicked: {}", e))??;
        Ok(())
    });

    Ok((ReceiverStream::new(record_rx), completion_task))
}

async fn run_kma_amr(
    config: Arc<RunConfig>,
    reduplicated_streams: Vec<ReceiverStream<ParseOutput>>,
    kma_db: PathBuf,
) -> Result<RunKmaHandle, PipelineError> {
    let stream_count = reduplicated_streams.len();

    if stream_count == 0 || stream_count > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunKMA expected one or two reduplicated FASTA streams, got {}",
            stream_count
        )));
    }

    let paired = stream_count == 2;

    // ------------------------------------------------------------
    // Adapt RunRedup output to KMA stdin.
    //
    // Single-end:
    //
    //     R1 -> KMA -i --
    //
    // Paired-end:
    //
    //     R1 --\
    //            -> R1,R2,R1,R2,... -> KMA -int --
    //     R2 --/
    //
    // RunRedup emits ParseOutput::Fasta records. Keep them structured
    // until stream_to_cmd serializes them to KMA stdin.
    // ------------------------------------------------------------

    let (kma_input_rx, input_adapter_task) = if paired {
        let mut streams = reduplicated_streams.into_iter();

        let mut r1 = streams
            .next()
            .expect("paired RunKMA input missing R1 stream");

        let mut r2 = streams
            .next()
            .expect("paired RunKMA input missing R2 stream");

        let (tx, rx) =
            mpsc::channel::<ParseOutput>(config.base_buffer_size);

        let task = tokio::spawn(async move {
            let mut pairs = 0u64;

            loop {
                let left = r1.next().await;
                let right = r2.next().await;

                match (left, right) {
                    (
                        Some(ParseOutput::Fasta(r1_record)),
                        Some(ParseOutput::Fasta(r2_record)),
                    ) => {
                        tx.send(ParseOutput::Fasta(r1_record))
                            .await
                            .map_err(|_| {
                                anyhow!(
                                    "RunKMA paired input receiver dropped on R1"
                                )
                            })?;

                        tx.send(ParseOutput::Fasta(r2_record))
                            .await
                            .map_err(|_| {
                                anyhow!(
                                    "RunKMA paired input receiver dropped on R2"
                                )
                            })?;

                        pairs += 1;
                    }

                    (None, None) => {
                        break;
                    }

                    (None, Some(_)) | (Some(_), None) => {
                        return Err(anyhow!(
                            "RunKMA R1/R2 streams have unequal record counts \
                             after {} complete pairs",
                            pairs
                        ));
                    }

                    (Some(left), Some(right)) => {
                        return Err(anyhow!(
                            "RunKMA expected paired FASTA records, got {:?} and {:?}",
                            left,
                            right
                        ));
                    }
                }
            }

            info!(
                "[AMR] RunKMA prepared {} interleaved read pairs",
                pairs
            );

            Ok::<(), anyhow::Error>(())
        });

        (rx, task)
    } else {
        let mut streams = reduplicated_streams.into_iter();

        let mut r1 = streams
            .next()
            .expect("single-end RunKMA input missing R1 stream");

        let (tx, rx) =
            mpsc::channel::<ParseOutput>(config.base_buffer_size);

        let task = tokio::spawn(async move {
            let mut reads = 0u64;

            while let Some(item) = r1.next().await {
                let ParseOutput::Fasta(record) = item else {
                    return Err(anyhow!(
                        "RunKMA expected single-end FASTA record, got {:?}",
                        item
                    ));
                };

                tx.send(ParseOutput::Fasta(record))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "RunKMA single-end input receiver dropped"
                        )
                    })?;

                reads += 1;
            }

            info!(
                "[AMR] RunKMA prepared {} single-end reads",
                reads
            );

            Ok::<(), anyhow::Error>(())
        });

        (rx, task)
    };

    // ------------------------------------------------------------
    // KMA requires an output prefix even though SAM is emitted on
    // stdout.
    //
    // Keep the TempDir in RunKmaHandle so KMA sidecars remain alive
    // for downstream RunRgiBwtKma processing.
    // ------------------------------------------------------------

    let kma_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        2,
        false,
    )
        .await?;

    let kma_output_prefix =
        kma_temp_dir.path().join("read_amr");

    let kma_config = KmaConfig::amr_stdin(
        kma_db,
        kma_output_prefix.clone(),
        paired,
    );

    let kma_args = generate_cli(
        KMA_TAG,
        &config,
        Some(&kma_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunKMA mode={}, args={:?}",
        if paired {
            "paired-interleaved-stdin"
        } else {
            "single-end-stdin"
        },
        kma_args
    );

    let stderr_path =
        config.out_dir.join("amr_kma_stderr.log");

    // ------------------------------------------------------------
    // Feed structured FASTA directly to KMA.
    // ------------------------------------------------------------

    let (child, stdin_task, stderr_task) =
        stream_to_cmd(
            config.clone(),
            kma_input_rx,
            KMA_TAG,
            kma_args,
            StreamDataType::IlluminaFastq,
            config.args.verbose,
            Some(stderr_path),
        )
            .await
            .map_err(PipelineError::Other)?;

    // ------------------------------------------------------------
    // Take KMA stdout without retaining the child mutex.
    //
    // stream_to_cmd's stdin task ultimately waits on the child.
    // Holding the child lock while consuming stdout would deadlock.
    // ------------------------------------------------------------

    let stdout = {
        let mut child_guard = child.lock().await;

        child_guard
            .stdout
            .take()
            .ok_or_else(|| {
                PipelineError::Other(anyhow!(
                    "RunKMA child stdout not available"
                ))
            })?
    };

    // ------------------------------------------------------------
    // KMA `-sam` emits serialized SAM on stdout.
    //
    // IMPORTANT:
    //
    // Keep this as RAW BYTES before fanout. This preserves the exact
    // SAM byte stream, including line endings, for samtools.
    //
    //                         ┌-> bytes_to_lines()
    //                         │      -> native SAM parser
    // KMA -> parse_bytes() -> fanout
    //                         │
    //                         └-> samtools fixmate stdin
    //
    // There is deliberately no raw-SAM file boundary here.
    // ------------------------------------------------------------

    let sam_bytes_rx = parse_bytes(
        stdout,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(PipelineError::Other)?;

    let (mut sam_branches, sam_fanout_done) =
        fanout_to_channels(
            ReceiverStream::new(sam_bytes_rx),
            2,
            "AMR KMA SAM bytes",
            &config,
            StreamDataType::JustBytes,
        )
            .await?;

    // ------------------------------------------------------------
    // Branch 1: native AMR parser.
    //
    // KMA stdout was chunked arbitrarily by parse_bytes(), so convert
    // the byte stream back into complete SAM lines before handing it
    // to consume_kma_sam().
    //
    // bytes_to_lines() preserves the newline in each emitted line.
    // consume_kma_sam() therefore needs to trim the line ending before
    // parsing (see small change below).
    // ------------------------------------------------------------

    let native_sam_bytes_rx = sam_branches.remove(0);

    let (sam_lines_rx, sam_lines_task) =
        bytes_to_lines(
            native_sam_bytes_rx,
            config.base_buffer_size,
        )
            .await
            .map_err(PipelineError::Other)?;

    let sam_stream =
        ReceiverStream::new(sam_lines_rx);

    // ------------------------------------------------------------
    // Branch 2: exact raw SAM bytes for samtools fixmate.
    //
    // Do NOT parse, rewrite, materialize, or otherwise modify this
    // branch. stream_to_cmd() can feed these ParseOutput::Bytes
    // directly into samtools stdin.
    // ------------------------------------------------------------

    let samtools_stream =
        sam_branches.remove(0);

    // ------------------------------------------------------------
    // Own KMA-side producer/process completion in one task.
    //
    // kma_temp_dir is deliberately NOT moved into this task. The
    // returned RunKmaHandle owns it so KMA sidecars survive process
    // completion.
    //
    // Neither SAM branch is consumed here. The caller must install
    // both downstream consumers immediately:
    //
    //   sam_stream      -> consume_kma_sam()
    //   samtools_stream -> samtools fixmate
    //
    // This ensures KMA stdout cannot backpressure indefinitely.
    // ------------------------------------------------------------

    let completion_task = tokio::spawn(async move {
        input_adapter_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunKMA input adapter task panicked: {}",
                    e
                )
            })??;

        stdin_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunKMA stdin task panicked: {}",
                    e
                )
            })??;

        stderr_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunKMA stderr task panicked: {}",
                    e
                )
            })??;

        info!(
            "[AMR] RunKMA process completed ({})",
            if paired {
                "paired-end"
            } else {
                "single-end"
            }
        );

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunKmaHandle {
        sam_stream,
        samtools_stream,
        sam_fanout_done,
        sam_lines_task,
        completion_task,
        output_temp_dir: kma_temp_dir,
        output_prefix: kma_output_prefix,
    })
}

fn parse_kma_sq_header(line: &str) -> anyhow::Result<Option<KmaReference>> {
    if !line.starts_with("@SQ\t") {
        return Ok(None);
    }

    let mut name: Option<&str> = None;
    let mut length: Option<u64> = None;

    for field in line.split('\t').skip(1) {
        if let Some(value) = field.strip_prefix("SN:") {
            name = Some(value);
        } else if let Some(value) = field.strip_prefix("LN:") {
            length = Some(
                value
                    .parse::<u64>()
                    .map_err(|e| {
                        anyhow!(
                            "Invalid @SQ LN value '{}' in KMA SAM header '{}': {}",
                            value,
                            line,
                            e
                        )
                    })?,
            );
        }
    }

    let name = name.ok_or_else(|| {
        anyhow!("KMA @SQ header missing SN field: {}", line)
    })?;

    let length = length.ok_or_else(|| {
        anyhow!("KMA @SQ header missing LN field: {}", line)
    })?;

    Ok(Some(KmaReference {
        name: name.to_string(),
        length,
    }))
}

fn is_valid_kma_alignment(record: &crate::utils::sambam::SamRecordRef<'_>) -> bool {
    !record.flag_is_unmapped()
        && record.rname != "*"
        && record.pos != 0
        && record.cigar != "*"
}

async fn consume_kma_sam(
    mut sam_stream: ReceiverStream<ParseOutput>,
) -> anyhow::Result<KmaSamResults> {
    let mut references = HashMap::<String, KmaReference>::new();
    let mut alignments = Vec::<KmaAlignment>::new();

    let mut total_records = 0u64;
    let mut invalid_alignment_records = 0u64;

    while let Some(item) = sam_stream.next().await {
        let ParseOutput::Bytes(line_bytes) = item else {
            return Err(anyhow!(
                "RunKMA expected line-oriented Bytes output"
            ));
        };

        let line = std::str::from_utf8(&line_bytes)
            .map_err(|e| {
                anyhow!(
            "RunKMA emitted non-UTF8 SAM line: {}",
            e
        )
            })?
            .trim_end_matches(['\r', '\n']);

        if line.starts_with('@') {
            if let Some(reference) = parse_kma_sq_header(line)? {
                references.insert(
                    reference.name.clone(),
                    reference,
                );
            }

            continue;
        }

        total_records += 1;

        let record = parse_sam_record(line)?;

        if !is_valid_kma_alignment(&record) {
            invalid_alignment_records += 1;
            continue;
        }

        alignments.push(KmaAlignment {
            qname: record.qname.to_string(),
            flag: record.flag,
            rname: record.rname.to_string(),
            pos: record.pos,
            mapq: record.mapq,
            cigar: record.cigar.to_string(),
            rnext: record.rnext.to_string(),
            seq: record.seq.to_string(),
        });
    }

    info!(
        "[AMR] RunKMA SAM parsed: {} records, {} valid mapped alignments, \
         {} non-alignment records, {} reference sequences",
        total_records,
        alignments.len(),
        invalid_alignment_records,
        references.len(),
    );

    Ok(KmaSamResults {
        references,
        alignments,
        total_records,
        invalid_alignment_records,
    })
}


/// Parse one `samtools view` SAM alignment record into the fields retained
/// by RGI's `get_qname_rname_sequence()` / `preload_alignments()`.
///
/// RGI retains SAM fields:
///
///   1  QNAME
///   2  FLAG
///   3  RNAME
///   4  POS
///   5  MAPQ
///   7  RNEXT / MRNM
///
/// `samtools view` emits the complete SAM alignment record, so this parser
/// selects those fields directly instead of reproducing RGI's:
///
///   cut -f 1,2,3,4,5,7 | sort -s -n -k 1,1
///
/// The sort is unnecessary for the Rust representation because alignments
/// will be grouped directly by RNAME.
fn parse_rgi_alignment(line: &str) -> anyhow::Result<RgiAlignment> {
    let fields: Vec<&str> = line.split('\t').collect();

    // We need through SAM field 7 (RNEXT), i.e. zero-based index 6.
    if fields.len() < 7 {
        return Err(anyhow!(
            "RunRgiBwtKma samtools view emitted malformed SAM record: \
             expected at least 7 tab-delimited fields, found {}: {:?}",
            fields.len(),
            line
        ));
    }

    let flag = fields[1].parse::<u16>().map_err(|e| {
        anyhow!(
            "RunRgiBwtKma invalid SAM FLAG {:?} for QNAME {:?}: {}",
            fields[1],
            fields[0],
            e
        )
    })?;

    let pos = fields[3].parse::<u64>().map_err(|e| {
        anyhow!(
            "RunRgiBwtKma invalid SAM POS {:?} for QNAME {:?}: {}",
            fields[3],
            fields[0],
            e
        )
    })?;

    let mapq = fields[4].parse::<u8>().map_err(|e| {
        anyhow!(
            "RunRgiBwtKma invalid SAM MAPQ {:?} for QNAME {:?}: {}",
            fields[4],
            fields[0],
            e
        )
    })?;

    Ok(RgiAlignment {
        qname: fields[0].to_string(),
        flag,
        rname: fields[2].to_string(),
        pos,
        mapq,
        mrnm: fields[6].to_string(),
    })
}

/// Evidence from BWT.py get_read_coverage/get_consensus_sequence/get_mutation_details.
/// Depth stays textual; mutation order and duplicates are intentionally retained.
#[derive(Debug, Default)]
struct KmaSidecars {
    read_depth: HashMap<String, String>,
    consensus_dna: HashMap<String, String>,
    mutations: HashMap<String, Vec<String>>,
}

fn kma_sidecar_path(prefix: &std::path::Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

fn parse_kma_res_line(line: &str, depths: &mut HashMap<String, String>) -> anyhow::Result<()> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields[0] == "#Template" {
        return Ok(());
    }
    if fields.len() < 9 {
        return Err(anyhow!("KMA .res row has {} fields; expected at least 9", fields.len()));
    }
    // BWT.py uses row[8].strip(), with last duplicate template winning.
    depths.insert(fields[0].to_string(), fields[8].trim().to_string());
    Ok(())
}

fn parse_kma_vcf_line(line: &str, mutations: &mut HashMap<String, Vec<String>>) -> anyhow::Result<()> {
    let fields: Vec<&str> = line.split('\t').collect();
    // Match BWT.py's contains test, not just starts_with('#').
    if fields[0].contains('#') {
        return Ok(());
    }
    if fields.len() < 5 {
        return Err(anyhow!("KMA VCF row has {} fields; expected at least 5", fields.len()));
    }
    // No QUAL/FILTER/ALT filtering, normalization, sorting, or deduplication.
    mutations.entry(fields[0].to_string()).or_default().push(format!(
        "{}{}{}", fields[3].trim(), fields[1].trim(), fields[4].trim()
    ));
    Ok(())
}

async fn read_kma_res(path: &std::path::Path) -> anyhow::Result<HashMap<String, String>> {
    let file = tokio::fs::File::open(path).await
        .map_err(|e| anyhow!("Open KMA sidecar {}: {}", path.display(), e))?;
    let mut lines = tokio::io::BufReader::new(file).lines();
    let mut depths = HashMap::new();
    let mut row = 0u64;
    while let Some(line) = lines.next_line().await? {
        row += 1;
        parse_kma_res_line(&line, &mut depths)
            .map_err(|e| anyhow!("{} row {}: {}", path.display(), row, e))?;
    }
    Ok(depths)
}

async fn read_kma_fsa(path: &std::path::Path) -> anyhow::Result<HashMap<String, String>> {
    let file = tokio::fs::File::open(path).await
        .map_err(|e| anyhow!("Open KMA sidecar {}: {}", path.display(), e))?;
    let mut lines = tokio::io::BufReader::new(file).lines();
    let mut sequences = HashMap::new();
    let mut id: Option<String> = None;
    let mut sequence = String::new();
    while let Some(line) = lines.next_line().await? {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.take() {
                if sequences.insert(previous.clone(), std::mem::take(&mut sequence)).is_some() {
                    return Err(anyhow!("Duplicate FASTA ID {:?} in {}", previous, path.display()));
                }
            }
            // SeqIO record.id is the first whitespace-separated header token.
            id = Some(header.split_whitespace().next()
                .ok_or_else(|| anyhow!("Empty FASTA header in {}", path.display()))?.to_string());
        } else if id.is_some() {
            // Biopython plain FASTA removes spaces and tabs, preserves case/gaps.
            sequence.extend(line.chars().filter(|c| *c != ' ' && *c != '\t'));
        }
    }
    if let Some(previous) = id {
        if sequences.insert(previous.clone(), sequence).is_some() {
            return Err(anyhow!("Duplicate FASTA ID {:?} in {}", previous, path.display()));
        }
    }
    Ok(sequences)
}

async fn read_kma_vcf(
    config: Arc<RunConfig>,
    path: PathBuf,
) -> anyhow::Result<HashMap<String, Vec<String>>> {
    let pigz_config = PigzConfig::decompress_stdout(path.clone());
    let args = generate_cli(PIGZ_TAG, &config, Some(&pigz_config))?;
    info!("[AMR] KMA VCF pigz args: {:?}", args);
    let (mut child, stderr_task) = spawn_cmd(
        config.clone(), PIGZ_TAG, args, config.args.verbose, None,
    ).await?;

    // Read stdout directly so I/O errors propagate (parse_lines only logs them).
    // The process is still launched exclusively through the shared launcher.
    let parse_result = async {
        let stdout = child.stdout.take().ok_or_else(|| anyhow!("KMA VCF pigz stdout unavailable"))?;
        let mut lines = tokio::io::BufReader::new(stdout).lines();
        let mut mutations = HashMap::new();
        let mut row = 0u64;
        while let Some(line) = lines.next_line().await? {
            row += 1;
            parse_kma_vcf_line(&line, &mut mutations)
                .map_err(|e| anyhow!("{} row {}: {}", path.display(), row, e))?;
        }
        Ok::<_, anyhow::Error>(mutations)
    }.await;
    // On parse/read failure, close stdout and terminate/reap the producer.
    if parse_result.is_err() {
        let _ = child.kill().await;
    }
    let status = child.wait().await;
    let stderr_result = stderr_task.await;
    let mutations = parse_result?;
    let status = status?;
    stderr_result.map_err(|e| anyhow!("KMA VCF pigz stderr task panicked: {}", e))??;
    if !status.success() {
        return Err(anyhow!("KMA VCF pigz failed for {}: {}", path.display(), status));
    }
    Ok(mutations)
}

async fn read_kma_sidecars(config: Arc<RunConfig>, prefix: &std::path::Path) -> anyhow::Result<KmaSidecars> {
    let read_depth = read_kma_res(&kma_sidecar_path(prefix, ".res")).await?;
    let consensus_dna = read_kma_fsa(&kma_sidecar_path(prefix, ".fsa")).await?;
    let mutations = read_kma_vcf(config, kma_sidecar_path(prefix, ".vcf.gz")).await?;
    info!(
        "[AMR] KMA sidecars parsed: res_references={} fsa_references={} vcf_references={} mutation_records={}",
        read_depth.len(), consensus_dna.len(), mutations.len(),
        mutations.values().map(Vec::len).sum::<usize>()
    );
    Ok(KmaSidecars { read_depth, consensus_dna, mutations })
}

async fn load_card_models(
    path: &std::path::Path,
) -> anyhow::Result<CardModels> {
    let bytes = tokio::fs::read(path)
        .await
        .map_err(|e| anyhow!(
            "Failed to read CARD metadata {}: {}",
            path.display(),
            e
        ))?;

    let root: serde_json::Value =
        serde_json::from_slice(&bytes)
            .map_err(|e| anyhow!(
                "Failed to parse CARD metadata {}: {}",
                path.display(),
                e
            ))?;

    let object = root.as_object().ok_or_else(|| {
        anyhow!(
            "CARD metadata root is not a JSON object: {}",
            path.display()
        )
    })?;

    let mut by_model_id = HashMap::new();
    let mut by_aro_accession = HashMap::new();

    for (key, value) in object {
        // Match BWT.py: if i.isdigit()
        if key.is_empty() || !key.bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }

        let model = value.as_object().ok_or_else(|| {
            anyhow!("CARD model {} is not a JSON object", key)
        })?;

        let required_string =
            |field: &str| -> anyhow::Result<String> {
                model
                    .get(field)
                    .and_then(serde_json::Value::as_str)
                    .map(str::to_owned)
                    .ok_or_else(|| {
                        anyhow!(
                            "CARD model {} missing string field {:?}",
                            key,
                            field
                        )
                    })
            };

        let model_id = required_string("model_id")?;
        let aro_accession = required_string("ARO_accession")?;
        let model_name = required_string("model_name")?;
        let model_type = required_string("model_type")?;

        let mut categories =
            HashMap::<String, Vec<String>>::new();

        if let Some(category_object) =
            model.get("ARO_category")
                .and_then(serde_json::Value::as_object)
        {
            for category in category_object.values() {
                let Some(category) = category.as_object() else {
                    continue;
                };

                let Some(class_name) = category
                    .get("category_aro_class_name")
                    .and_then(serde_json::Value::as_str)
                else {
                    continue;
                };

                let Some(category_name) = category
                    .get("category_aro_name")
                    .and_then(serde_json::Value::as_str)
                else {
                    continue;
                };

                let values =
                    categories.entry(class_name.to_string())
                        .or_default();

                // Match BWT.py uniqueness behavior.
                if !values.iter().any(|v| v == category_name) {
                    values.push(category_name.to_string());
                }
            }
        }

        let mut taxon = Vec::<String>::new();

        if let Some(sequences) = model
            .get("model_sequences")
            .and_then(|v| v.get("sequence"))
            .and_then(serde_json::Value::as_object)
        {
            for sequence in sequences.values() {
                let Some(name) = sequence
                    .get("NCBI_taxonomy")
                    .and_then(|v| v.get("NCBI_taxonomy_name"))
                    .and_then(serde_json::Value::as_str)
                else {
                    continue;
                };

                // Exact BWT.py:
                // " ".join(name.split()[:2])
                let taxa = name
                    .split_whitespace()
                    .take(2)
                    .collect::<Vec<_>>()
                    .join(" ");

                if !taxa.is_empty()
                    && !taxon.iter().any(|existing| existing == &taxa)
                {
                    taxon.push(taxa);
                }
            }
        }

        let card_model = CardModel {
            model_id: model_id.clone(),
            aro_accession: aro_accession.clone(),
            model_name,
            model_type,
            categories,
            taxon,
        };

        by_model_id.insert(model_id, card_model.clone());
        by_aro_accession.insert(aro_accession, card_model);
    }

    info!(
        "[AMR] CARD metadata loaded: models={} aro_accessions={}",
        by_model_id.len(),
        by_aro_accession.len()
    );

    Ok(CardModels {
        by_model_id,
        by_aro_accession,
    })
}

fn value_after_colon(field: &str) -> Option<&str> {
    field.split_once(':').map(|(_, value)| value)
}

fn resolve_card_model<'a>(
    alignment_hit: &str,
    models: &'a CardModels,
) -> Option<&'a CardModel> {
    let fields: Vec<&str> = alignment_hit.split('|').collect();

    if alignment_hit.starts_with("Prevalence_Sequence_ID")
        || alignment_hit.starts_with("ARO:")
    {
        let model_id =
            fields.get(1)
                .and_then(|field| value_after_colon(field));

        return model_id
            .and_then(|id| models.by_model_id.get(id));
    }

    let accession =
        fields.get(4)
            .and_then(|field| value_after_colon(field));

    accession
        .and_then(|aro| models.by_aro_accession.get(aro))
}

pub async fn run(
    config: Arc<RunConfig>,
) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir()
        .map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();
    info!("Starting AMR pipeline.");

    // WDL branch boundary: only input preparation differs by input mode.
    // Return live reads BEFORE awaiting preparation: bounded producers need
    // the shared KMA consumer below to be installed before they can finish.
    let (kma_streams, preparation_task) = match AmrInputMode::from_config(&config)? {
        AmrInputMode::Raw => {
            info!(
                "[AMR] raw input selected; running shared host filtering"
            );

            // ------------------------------------------------------------
            // Raw-input setup
            // ------------------------------------------------------------

            let (
                file1_path,
                file2_path,
                sample_base_buf,
                _sample_base,
            ) = validate_file_inputs(
                &config,
                &cwd,
                true,
            )
                .await?;

            let paired = file2_path.is_some();

            let mut final_temp_dirs: Vec<TempDir> = Vec::new();

            let mut cleanup_tasks:
                Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> =
                Vec::new();

            let mut cleanup_receivers:
                Vec<
                    oneshot::Receiver<
                        anyhow::Result<(), anyhow::Error>
                    >
                > =
                Vec::new();

            // ------------------------------------------------------------
            // Shared host filtering
            //
            // AMR retains the exact host-filtered stream immediately
            // before deduplication for later AMR processing.
            // ------------------------------------------------------------

            let host_filter = run_host_filter(
                config.clone(),
                &cwd,
                &out_dir,
                &file1_path,
                file2_path.clone(),
                &sample_base_buf,
                paired,
                true,
                true,
                &mut final_temp_dirs,
                &mut cleanup_tasks,
                &mut cleanup_receivers,
            )
                .await?;




            let HostFilterCompletion {
                outputs,

                val_count_task,
                qc_count_result_rx,
                ercc_count_rx,

                ercc_bt2_bam_write_handle,
                host_bt2_bam_write_handle,
                optional_human_bam_write_handle,

                host_bt2_count_rx,

                ercc_bt2_bam_path: _,
                host_hisat2_count_rx,
                host_bt2_bam_path: _,

                kallisto_exit_task,
                kallisto_ercc_tx,
                kallisto_ercc_rx,

                dedup_count_rx,
                subsample_count_rx,
            } = host_filter;

            // ------------------------------------------------------------
            // Split host-filter outputs by ownership.
            //
            // Raw AMR consumes:
            //
            //   * pre_dedup_interleaved
            //   * subsampled_interleaved
            //   * duplicate_clusters
            //
            // run_host_filter() also creates the normal short-read
            // deinterleaved non_host_r1/non_host_r2 outputs. AMR does not
            // use those outputs, but they MUST still be drained. Otherwise
            // their bounded channels can fill and backpressure the shared
            // subsampled fanout.
            // ------------------------------------------------------------

            let HostFilterOutputs {
                non_host_r1,
                non_host_r2,
                subsampled_interleaved,
                duplicate_clusters,
                pre_dedup_interleaved,
            } = outputs;

            // ------------------------------------------------------------
            // Drain the unused normal short-read non-host outputs.
            //
            // Install these consumers immediately, before awaiting any
            // host-filter completion work.
            // ------------------------------------------------------------

            let mut unused_non_host_r1 = non_host_r1;

            cleanup_tasks.push(tokio::spawn(async move {
                while unused_non_host_r1.next().await.is_some() {}
                Ok::<(), anyhow::Error>(())
            }));

            let has_non_host_r2 = non_host_r2.is_some();

            if let Some(mut unused_non_host_r2) = non_host_r2 {
                cleanup_tasks.push(tokio::spawn(async move {
                    while unused_non_host_r2.next().await.is_some() {}
                    Ok::<(), anyhow::Error>(())
                }));
            }

            // ------------------------------------------------------------
            // Validate and take the AMR-specific retained streams.
            // ------------------------------------------------------------

            let pre_dedup_stream =
                pre_dedup_interleaved.ok_or_else(|| {
                    PipelineError::Other(anyhow!(
                        "AMR raw RunRedup missing retained interleaved pre-dedup stream"
                    ))
                })?;

            let subsampled_stream =
                subsampled_interleaved.ok_or_else(|| {
                    PipelineError::Other(anyhow!(
                        "AMR raw RunRedup missing retained interleaved subsampled stream"
                    ))
                })?;

            info!(
                "[AMR] host-filter streaming outputs ready: \
                 pre_dedup_interleaved=true, \
                 subsampled_interleaved=true, \
                 non_host_r1=draining, non_host_r2={}, \
                 duplicate_clusters=pending",
                has_non_host_r2,
            );

            // ------------------------------------------------------------
            // Install RunRedup consumers immediately.
            //
            // Both retained streams are structured interleaved FASTQ:
            //
            //   pre_dedup_stream
            //       exact host-filtered records immediately before dedup
            //
            //   subsampled_stream
            //       deduplicated + subsampled representatives
            //
            // duplicate_clusters completes only after dedup reaches EOF.
            // run_redup() installs consumers for both live streams before
            // awaiting that task.
            // ------------------------------------------------------------

            let redup = run_redup(
                config.clone(),
                pre_dedup_stream,
                subsampled_stream,
                paired,
                duplicate_clusters,
            )
                .await?;

            // ------------------------------------------------------------
            // Raw WDL topology:
            //
            //     RunRedupRaw -> RunSpades
            //
            // Install the SPAdes consumer immediately so every bounded
            // RunRedup output has a downstream reader before any
            // host-filter completion work is awaited.
            // ------------------------------------------------------------


            let mut spades_streams = Vec::new();
            let mut downstream_streams = Vec::new();
            let mut redup_fanout_tasks = Vec::new();

            for (idx, redup_stream) in
                redup.redups_streams.into_iter().enumerate()
            {
                let (mut branches, done_rx) = fanout_to_channels(
                    redup_stream,
                    2,
                    &format!("AMR redups part {}", idx + 1),
                    &config,
                    StreamDataType::IlluminaFastq,
                )
                    .await?;

                spades_streams.push(
                    ReceiverStream::new(branches.remove(0))
                );

                downstream_streams.push(
                    ReceiverStream::new(branches.remove(0))
                );

                redup_fanout_tasks.push(done_rx);
            }



            let spades = run_spades_amr(
                config.clone(),
                spades_streams,
            )
                .await?;

            let redup_completion_task = redup.completion_task;
            let preparation_out_dir = out_dir.clone();
            let preparation_task = tokio::spawn(async move {
                let out_dir = preparation_out_dir;
                // ------------------------------------------------------------
                // Drain validation / QC results.
                // ------------------------------------------------------------

                let validation_and_qc_fut = async move {
                    let stats =
                        join_with_error_handling(val_count_task).await?;

                    info!(
                        "[AMR] validation: raw_reads={}, raw_bases={}, \
                         validated_reads={}, validated_bases={}, \
                         undersized={}, oversized={}",
                        stats.raw_reads,
                        stats.raw_bases,
                        stats.validated,
                        stats.validated_bases,
                        stats.undersized,
                        stats.oversized,
                    );

                    let qc_fastp_read_count =
                        match qc_count_result_rx.await {
                            Ok(Ok(count)) => {
                                info!(
                                    "[AMR] fastp output reads: {}",
                                    count
                                );
                                count
                            }

                            Ok(Err(e)) => {
                                error!(
                                    "[AMR] error getting fastp \
                                     read count: {}",
                                    e
                                );
                                0
                            }

                            Err(e) => {
                                error!(
                                    "[AMR] fastp count receiver \
                                     dropped: {}",
                                    e
                                );
                                0
                            }
                        };

                    let ercc_filtered_count =
                        ercc_count_rx
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "ERCC count receiver failed: {}",
                                    e
                                ))
                            })?;

                    Ok::<_, PipelineError>((
                        stats,
                        ercc_filtered_count,
                        qc_fastp_read_count,
                    ))
                };

                // ------------------------------------------------------------
                // Drain BAM writers and host counts.
                // ------------------------------------------------------------

                let host_counts_fut = async move {
                    ercc_bt2_bam_write_handle
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "ERCC BAM writer task panicked: {}",
                                e
                            ))
                        })?
                        .map_err(PipelineError::Other)?;

                    host_bt2_bam_write_handle
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Host Bowtie2 BAM writer task \
                                 panicked: {}",
                                e
                            ))
                        })?
                        .map_err(PipelineError::Other)?;

                    if let Some(handle) =
                        optional_human_bam_write_handle
                    {
                        handle
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Human BAM writer task \
                                     panicked: {}",
                                    e
                                ))
                            })?
                            .map_err(PipelineError::Other)?;
                    }

                    let host_bt2_count =
                        host_bt2_count_rx
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Host Bowtie2 count receiver \
                                     failed: {}",
                                    e
                                ))
                            })?;

                    let host_hisat2_count =
                        host_hisat2_count_rx
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Host HISAT2 count receiver \
                                     failed: {}",
                                    e
                                ))
                            })?;

                    info!(
                        "[AMR] host filtering counts: \
                         bowtie2={}, hisat2={}",
                        host_bt2_count,
                        host_hisat2_count
                    );

                    Ok::<_, PipelineError>((
                        host_bt2_count,
                        host_hisat2_count,
                    ))
                };

                // ------------------------------------------------------------
                // Drain Kallisto.
                // ------------------------------------------------------------

                let out_dir_kallisto = out_dir.clone();

                let kallisto_fut = async move {
                    let kallisto_exit =
                        join_with_error_handling(
                            kallisto_exit_task
                        )
                            .await;

                    match kallisto_exit {
                        Ok(_) => {
                            info!(
                                "[AMR] Kallisto completed successfully"
                            );
                        }

                        Err(e) => {
                            warn!(
                                "[AMR] Kallisto failed \
                                 (non-fatal, possibly no ERCC \
                                 spike-in): {}",
                                e
                            );
                        }
                    }

                    let results_task =
                        kallisto_results(
                            out_dir_kallisto.join("kallisto"),
                            kallisto_ercc_tx,
                        )
                            .await?;

                    join_with_error_handling(results_task).await?;

                    let results =
                        kallisto_ercc_rx
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Kallisto results receiver \
                                     failed: {}",
                                    e
                                ))
                            })?;

                    info!(
                        "[AMR] Kallisto results drained: \
                         ercc_entries={}, transcript_mappings={}",
                        results.ercc_counts.len(),
                        results.transcript_to_gene.len(),
                    );

                    Ok::<(), PipelineError>(())
                };

                // ------------------------------------------------------------
                // Drain dedup / subsampling counts.
                // ------------------------------------------------------------

                let dedup_and_subsample_fut = async move {
                    let unique_units =
                        dedup_count_rx
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Dedup count receiver failed: {}",
                                    e
                                ))
                            })?;

                    let unique_reads =
                        unique_units * if paired { 2 } else { 1 };

                    let subsample_units =
                        subsample_count_rx
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Subsample count receiver \
                                     failed: {}",
                                    e
                                ))
                            })?;

                    let subsample_reads =
                        subsample_units * if paired { 2 } else { 1 };

                    info!(
                        "[AMR] dedup/subsample: \
                         unique_reads={}, subsampled_reads={}",
                        unique_reads,
                        subsample_reads
                    );

                    Ok::<_, PipelineError>((
                        unique_reads,
                        subsample_reads,
                    ))
                };

                let (validation_and_qc, host_counts, (), dedup_and_subsample) =
                    tokio::try_join!(
                        validation_and_qc_fut,
                        host_counts_fut,
                        kallisto_fut,
                        dedup_and_subsample_fut,
                    )?;
                let (
                    validation_stats,
                    ercc_filtered_count,
                    qc_fastp_read_count,
                ) = validation_and_qc;

                let (
                    host_bt2_filtered_count,
                    host_hisat2_filtered_count,
                ) = host_counts;

                let (
                    unique_reads,
                    subsample_reads,
                ) = dedup_and_subsample;

                // ------------------------------------------------------------
                // Write host-filter result artifacts.
                //
                // These exist only for raw AMR input because only raw input
                // executes the host-filter stage.
                // ------------------------------------------------------------

                write_count_file(
                    &out_dir,
                    "fastqs.count",
                    validation_stats.raw_reads,
                )
                    .await?;

                write_count_file(
                    &out_dir,
                    "validate_input_out.count",
                    validation_stats.validated,
                )
                    .await?;

                let validation_summary = serde_json::json!({
                    "<50": validation_stats.validation_too_short,
                    "50-500": validation_stats.validation_normal,
                    "500-10000": validation_stats.validation_long,
                    "10000+": validation_stats.validation_too_long,
                });

                let validation_summary_path =
                    out_dir.join("validate_input_summary.json");

                tokio::fs::write(
                    &validation_summary_path,
                    serde_json::to_vec(
                        &validation_summary
                    )
                        .map_err(|e| {
                            PipelineError::Other(e.into())
                        })?,
                )
                    .await
                    .map_err(|e| {
                        PipelineError::IOError(format!(
                            "Failed to write {}: {}",
                            validation_summary_path.display(),
                            e
                        ))
                    })?;

                write_count_file(
                    &out_dir,
                    "fastp_out.count",
                    qc_fastp_read_count,
                )
                    .await?;

                write_count_file(
                    &out_dir,
                    "bowtie2_ercc_filtered_out.count",
                    ercc_filtered_count,
                )
                    .await?;

                write_count_file(
                    &out_dir,
                    "bowtie2_host_filtered_out.count",
                    host_bt2_filtered_count,
                )
                    .await?;

                write_count_file(
                    &out_dir,
                    "hisat2_host_filtered_out.count",
                    host_hisat2_filtered_count,
                )
                    .await?;

                write_count_file(
                    &out_dir,
                    "czid_dedup_out.count",
                    unique_reads,
                )
                    .await?;

                write_count_file(
                    &out_dir,
                    "subsampled_out.count",
                    subsample_reads,
                )
                    .await?;

                redup_completion_task
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                "RunRedup task panicked: {}",
                e
            ))
                    })?
                    .map_err(PipelineError::Other)?;

                // Every RunRedup item must have reached both SPAdes and KMA.
                for (idx, done_rx) in redup_fanout_tasks.into_iter().enumerate() {
                    done_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                    "RunRedup fanout {} completion receiver dropped: {}",
                    idx + 1,
                    e
                ))
                        })?
                        .map_err(PipelineError::Other)?;
                }

                spades.completion_task
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                "RunSpades task panicked: {}",
                e
            ))
                    })?
                    .map_err(PipelineError::Other)?;

                info!(
                    "[AMR] raw contigs ready: {} (scaffolds: {})",
                    spades.contigs_path.display(),
                    spades.scaffolds_path.display()
                );

                // ------------------------------------------------------------
                // Drain cleanup owned by the host-filter stage.
                // ------------------------------------------------------------

                info!(
                    "[AMR] draining {} cleanup tasks and \
                     {} cleanup receivers",
                    cleanup_tasks.len(),
                    cleanup_receivers.len()
                );

                for task in cleanup_tasks {
                    task
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "AMR host-filter cleanup task \
                                 panicked: {}",
                                e
                            ))
                        })?
                        .map_err(PipelineError::Other)?;
                }

                for rx in cleanup_receivers {
                    rx.await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "AMR host-filter cleanup receiver \
                                 dropped: {}",
                                e
                            ))
                        })?
                        .map_err(PipelineError::Other)?;
                }

                // Keep host-filter temp directories alive until all
                // host-filter outputs and cleanup work have completed.
                info!(
                    "[AMR] raw-input host filtering complete: \
                     {} retained temp directories can now be released",
                    final_temp_dirs.len()
                );
                Ok::<(), PipelineError>(())
            });
            (downstream_streams, preparation_task)
        }
        AmrInputMode::Filtered => {
            info!(
                "[AMR] filtered input selected; \
                 skipping host filtering"
            );


            let (file1_path, file2_path, _, _) =
                validate_file_inputs(&config, &cwd, true).await?;

            let paired = file2_path.is_some();

            let pre_dedup_r1 = config
                .args
                .amr_pre_dedup_file1
                .as_deref()
                .ok_or_else(|| PipelineError::InvalidConfig(
                    "--amr-pre-dedup-file1 is required for --amr-input-mode filtered".to_string(),
                ))?;
            let pre_dedup_r1 = resolve_existing_input_path(pre_dedup_r1, &cwd, true)?;

            let pre_dedup_r2 = match config.args.amr_pre_dedup_file2.as_deref() {
                Some(path) => Some(resolve_existing_input_path(path, &cwd, true)?),
                None => None,
            };

            if file2_path.is_some() != pre_dedup_r2.is_some() {
                return Err(PipelineError::InvalidConfig(
                    "Filtered AMR RunRedup requires file2 and amr-pre-dedup-file2 to either both be present or both be absent".to_string(),
                ));
            }

            let clusters_path = config
                .args
                .amr_duplicate_clusters
                .as_deref()
                .ok_or_else(|| PipelineError::InvalidConfig(
                    "--amr-duplicate-clusters is required for --amr-input-mode filtered".to_string(),
                ))?;
            let clusters_path = resolve_existing_input_path(clusters_path, &cwd, true)?;

            let mut input_tasks = Vec::new();

            let (subsampled_r1, task) =
                fastq_file_to_record_stream(file1_path, config.clone(), "amr_filtered_subsampled_r1").await?;
            input_tasks.push(task);

            let subsampled_stream = if let Some(path) = file2_path {
                let (subsampled_r2, task) =
                    fastq_file_to_record_stream(path, config.clone(), "amr_filtered_subsampled_r2").await?;
                input_tasks.push(task);

                let (interleaved_rx, interleave_task) = interleave_fastq_streams(
                    subsampled_r1.into_inner(),
                    subsampled_r2.into_inner(),
                    config.base_buffer_size * 4,
                )
                    .await?;
                input_tasks.push(interleave_task);
                ReceiverStream::new(interleaved_rx)
            } else {
                subsampled_r1
            };

            let (pre_r1, task) =
                fastq_file_to_record_stream(
                    pre_dedup_r1,
                    config.clone(),
                    "amr_filtered_pre_dedup_r1",
                )
                    .await?;

            input_tasks.push(task);

            let pre_dedup_stream = if let Some(path) = pre_dedup_r2 {
                let (pre_r2, task) =
                    fastq_file_to_record_stream(
                        path,
                        config.clone(),
                        "amr_filtered_pre_dedup_r2",
                    )
                        .await?;

                input_tasks.push(task);

                let (interleaved_rx, interleave_task) =
                    interleave_fastq_streams(
                        pre_r1.into_inner(),
                        pre_r2.into_inner(),
                        config.base_buffer_size * 4,
                    )
                        .await?;

                input_tasks.push(interleave_task);

                ReceiverStream::new(interleaved_rx)
            } else {
                pre_r1
            };

            let redup = run_redup(
                config.clone(),
                pre_dedup_stream,
                subsampled_stream,
                paired,
                load_duplicate_clusters_csv(clusters_path),
            )
                .await?;



            let downstream_streams = redup.redups_streams;
            let preparation_task = tokio::spawn(async move {
                redup.completion_task.await
                    .map_err(|e| PipelineError::Other(anyhow!(
                        "Filtered RunRedup task panicked: {}", e
                    )))?
                    .map_err(PipelineError::Other)?;
                for task in input_tasks {
                    task.await
                        .map_err(|e| PipelineError::Other(anyhow!(
                            "Filtered FASTQ adapter task panicked: {}", e
                        )))?
                        .map_err(PipelineError::Other)?;
                }
                Ok::<(), PipelineError>(())
            });
            (downstream_streams, preparation_task)
        }
    };

    // Shared WDL RunRgiBwtKma: BOTH input modes execute everything below.
    // These scratch directories belong to shared processing, not host filtering.
    let mut final_temp_dirs: Vec<TempDir> = Vec::new();


    let kma_db = config
        .args
        .amr_kma_db
        .as_deref()
        .ok_or_else(|| {
            PipelineError::InvalidConfig(
                "--amr-kma-db is required for AMR".to_string()
            )
        })?;

    let kma = run_kma_amr(
        config.clone(),
        kma_streams,
        PathBuf::from(kma_db),
    )
        .await?;


    let RunKmaHandle {
        sam_stream: kma_sam_stream,
        samtools_stream: kma_samtools_stream,
        sam_fanout_done,
        sam_lines_task,
        completion_task: kma_completion_task,
        output_temp_dir: kma_output_temp_dir,
        output_prefix: kma_output_prefix,
    } = kma;
    // Keep KMA sidecars alive through downstream AMR post-processing.
    final_temp_dirs.push(kma_output_temp_dir);


    // ------------------------------------------------------------
    // Parse KMA SAM stdout directly into the native AMR
    // representation.
    //
    // This consumer is installed immediately after RunKMA returns so
    // KMA stdout can never block on the bounded line channel.
    // ------------------------------------------------------------

    let kma_sam_task = tokio::spawn(
        consume_kma_sam(kma_sam_stream)
    );

    // ------------------------------------------------------------
    // Stream raw KMA SAM directly through the same samtools fixmate
    // transformation used by RGI.
    //
    // RGI:
    //
    //   samtools fixmate
    //       --verbosity 0
    //       -z all
    //       --threads N
    //       -m
    //       input.sam
    //       -O bam
    //       output.bam
    //
    // Here both file operands are "-" so the transformation remains
    // fully streaming:
    //
    //   KMA stdout -> samtools fixmate stdin -> BAM stdout
    // ------------------------------------------------------------

    let fixmate_threads =
        config.thread_allocation(SAMTOOLS_TAG, Some("fixmate"));

    let fixmate_config =
        SamtoolsConfig::new(SamtoolsSubcommand::Fixmate)
            .option("--verbosity", "0")
            .option("-z", "all")
            .option("--threads", fixmate_threads.to_string())
            .flag("-m")
            .option("-O", "bam")
            .positional("-")
            .positional("-");

    let fixmate_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&fixmate_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma samtools fixmate args: {:?}",
        fixmate_args
    );

    let (
        fixmate_child,
        fixmate_stdin_task,
        fixmate_stderr_task,
    ) = stream_to_cmd(
        config.clone(),
        kma_samtools_stream,
        SAMTOOLS_TAG,
        fixmate_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(PipelineError::Other)?;

    // ------------------------------------------------------------
    // Stream fixmate BAM stdout directly into samtools sort.
    //
    // RGI:
    //
    //   samtools sort
    //       --threads N
    //       -T <prefix>.sorted
    //       -o <sorted.bam>
    //       <unsorted.bam>
    //
    // We preserve the sort operation but keep the pipeline streaming:
    //
    //   fixmate stdout -> samtools sort stdin -> sorted BAM stdout
    //
    // No intermediate unsorted BAM is materialized.
    // ------------------------------------------------------------

    let fixmate_stdout = {
        let mut child_guard = fixmate_child.lock().await;

        child_guard
            .stdout
            .take()
            .ok_or_else(|| {
                PipelineError::Other(anyhow!(
        "RunRgiBwtKma samtools fixmate stdout not available"
    ))
            })?
    };

    let fixmate_bam_rx = parse_bytes(
        fixmate_stdout,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(PipelineError::Other)?;

    let sort_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    let sort_temp_prefix =
        sort_temp_dir.path().join("sr_amr_report.sorted");

    let sort_threads =
        config.thread_allocation(SAMTOOLS_TAG, Some("sort"));

    let sort_config =
        SamtoolsConfig::new(SamtoolsSubcommand::Sort)
            .option("--threads", sort_threads.to_string())
            .option(
                "-T",
                sort_temp_prefix.to_string_lossy().to_string(),
            )
            .option("-O", "bam")
            .positional("-");

    let sort_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&sort_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma samtools sort args: {:?}",
        sort_args
    );

    let (
        sort_child,
        sort_stdin_task,
        sort_stderr_task,
    ) = stream_to_cmd(
        config.clone(),
        fixmate_bam_rx,
        SAMTOOLS_TAG,
        sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(PipelineError::Other)?;



    // ------------------------------------------------------------
    // Stream the coordinate-sorted BAM through the same predicate-free
    // bamtools filter transformation used by RGI, then re-emit it with
    // samtools before materializing the final length_100 BAM:
    //
    //   sort_stdout
    //       -> parse_bytes
    //       -> stream_to_cmd(BAMTOOLS_TAG, ["filter"])
    //       -> parse_bytes
    //       -> stream_to_cmd(SAMTOOLS_TAG, ["view", "-b", "-"])
    //       -> parse_bytes
    //       -> write_byte_stream_to_file(length_100.bam)
    //
    // BamTools 2.5.2 accepts BAM on stdin and emits the same decoded
    // alignments as its file-based filter invocation, but its stdout BAM
    // omits the BGZF EOF marker. `samtools view -b -` re-emits a properly
    // finalized BAM while preserving the alignment stream. This lets us
    // avoid both RGI's sorted.temp.bam materialization and its unused
    // intermediate BAI; the only BAM file boundary is the durable
    // sr_amr_report.sorted.length_100.bam required downstream.
    // ------------------------------------------------------------

    let sort_stdout = {
        let mut child_guard = sort_child.lock().await;

        child_guard
            .stdout
            .take()
            .ok_or_else(|| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools sort stdout not available"
                ))
            })?
    };

    let sorted_bam_rx = parse_bytes(
        sort_stdout,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(PipelineError::Other)?;

    let bamtools_config =
        BamtoolsConfig::new(BamtoolsSubcommand::Filter);

    let bamtools_args = generate_cli(
        BAMTOOLS_TAG,
        &config,
        Some(&bamtools_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma bamtools filter args: {:?}",
        bamtools_args
    );

    let (
        bamtools_child,
        bamtools_stdin_task,
        bamtools_stderr_task,
    ) = stream_to_cmd(
        config.clone(),
        sorted_bam_rx,
        BAMTOOLS_TAG,
        bamtools_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(PipelineError::Other)?;

    let bamtools_stdout = {
        let mut child_guard = bamtools_child.lock().await;

        child_guard
            .stdout
            .take()
            .ok_or_else(|| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma bamtools filter stdout not available"
                ))
            })?
    };

    let bamtools_bam_rx = parse_bytes(
        bamtools_stdout,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(PipelineError::Other)?;

    let finalize_config =
        SamtoolsConfig::new(SamtoolsSubcommand::View)
            .flag("-b")
            .positional("-");

    let finalize_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&finalize_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma samtools BAM finalizer args: {:?}",
        finalize_args
    );

    let (
        finalize_child,
        finalize_stdin_task,
        finalize_stderr_task,
    ) = stream_to_cmd(
        config.clone(),
        bamtools_bam_rx,
        SAMTOOLS_TAG,
        finalize_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(PipelineError::Other)?;

    let finalize_stdout = {
        let mut child_guard = finalize_child.lock().await;

        child_guard
            .stdout
            .take()
            .ok_or_else(|| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools BAM finalizer stdout not available"
                ))
            })?
    };

    let final_bam_rx = parse_bytes(
        finalize_stdout,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(PipelineError::Other)?;

    let sorted_bam_path =
        out_dir.join("sr_amr_report.sorted.length_100.bam");

    let sorted_bam_write_task = write_byte_stream_to_file(
        &sorted_bam_path,
        ReceiverStream::new(final_bam_rx),
        config.clone(),
        StreamDataType::JustBytes,
        "RunRgiBwtKma length_100 BAM",
        false,
    )
        .await
        .map_err(PipelineError::Other)?;

    // Keep samtools sort scratch alive until sort has completed.
    final_temp_dirs.push(sort_temp_dir);


    // ------------------------------------------------------------
    // Complete the streaming
    // KMA -> fixmate -> sort -> bamtools -> samtools view -> BAM chain.
    //
    // All stdout consumers have already been installed above, so it
    // is now safe to await producer/process completion. The final BAM
    // is not considered ready until every stdin/process task has
    // completed successfully and samtools view stdout has been fully
    // materialized.
    // ------------------------------------------------------------

    let rgi_bwt_bam_fut = async move {
        let (
            fixmate_stdin_result,
            sort_stdin_result,
            bamtools_stdin_result,
            finalize_stdin_result,
            sorted_bam_write_result,
        ) = tokio::join!(
            fixmate_stdin_task,
            sort_stdin_task,
            bamtools_stdin_task,
            finalize_stdin_task,
            sorted_bam_write_task,
        );

        fixmate_stdin_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools fixmate stdin task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        sort_stdin_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools sort stdin task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        bamtools_stdin_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma bamtools filter stdin task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        finalize_stdin_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools BAM finalizer stdin task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        sorted_bam_write_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma length_100 BAM writer task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        let (
            fixmate_stderr_result,
            sort_stderr_result,
            bamtools_stderr_result,
            finalize_stderr_result,
        ) = tokio::join!(
            fixmate_stderr_task,
            sort_stderr_task,
            bamtools_stderr_task,
            finalize_stderr_task,
        );

        fixmate_stderr_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools fixmate stderr task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        sort_stderr_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools sort stderr task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        bamtools_stderr_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma bamtools filter stderr task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        finalize_stderr_result
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools BAM finalizer stderr task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        info!(
            "[AMR] RunRgiBwtKma length_100 BAM ready: {}",
            sorted_bam_path.display()
        );

        Ok::<PathBuf, PipelineError>(sorted_bam_path)
    };

    let sorted_bam_path = rgi_bwt_bam_fut.await?;

    // ------------------------------------------------------------
    // Index the completed coordinate-sorted AMR BAM.
    //
    // rgi_bwt_bam_fut does not return until the final
    // sr_amr_report.sorted.length_100.bam writer has reached EOF.
    // samtools index therefore starts only after the durable BAM
    // file boundary is complete.
    // ------------------------------------------------------------

    let sorted_bai_path =
        PathBuf::from(format!("{}.bai", sorted_bam_path.display()));

    let index_threads =
        config.thread_allocation(SAMTOOLS_TAG, Some("index"));

    let index_config =
        SamtoolsConfig::new(SamtoolsSubcommand::Index)
            .option("-@", index_threads.to_string())
            .flag("-b")
            .positional(
                sorted_bam_path
                    .to_string_lossy()
                    .to_string(),
            )
            .positional(
                sorted_bai_path
                    .to_string_lossy()
                    .to_string(),
            );

    let index_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&index_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma samtools index args: {:?}",
        index_args
    );

    let (mut index_child, index_stderr_task) =
        spawn_cmd(
            config.clone(),
            SAMTOOLS_TAG,
            index_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(PipelineError::Other)?;

    let index_status =
        index_child
            .wait()
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
        "RunRgiBwtKma samtools index wait failed: {}",
        e
    ))
            })?;

    index_stderr_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
    "RunRgiBwtKma samtools index stderr task panicked: {}",
    e
            ))
        })?
        .map_err(PipelineError::Other)?;

    if !index_status.success() {
        return Err(PipelineError::Other(anyhow!(
            "RunRgiBwtKma samtools index failed with exit status {}",
            index_status
        )));
    }

    if !sorted_bai_path.is_file() {
        return Err(PipelineError::Other(anyhow!(
            "RunRgiBwtKma samtools index completed but BAI was not created: {}",
            sorted_bai_path.display()
        )));
    }

    info!(
        "[AMR] RunRgiBwtKma BAM index ready: {}",
        sorted_bai_path.display()
    );

    // ------------------------------------------------------------
    // Pull per-reference alignment counts with samtools idxstats.
    //
    // RGI materializes this stdout as sr_amr_report.temp.txt and
    // get_reads_count() later keeps only rows with mapped > 0.
    // There is no external-tool requirement for that text-file
    // boundary, so consume idxstats stdout directly and preserve the
    // same mapped-reference semantics in memory.
    // ------------------------------------------------------------

    let idxstats_config =
        SamtoolsConfig::new(SamtoolsSubcommand::Idxstats)
            .positional(
                sorted_bam_path
                    .to_string_lossy()
                    .to_string(),
            );

    let idxstats_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&idxstats_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma samtools idxstats args: {:?}",
        idxstats_args
    );

    let (mut idxstats_child, idxstats_stderr_task) =
        spawn_cmd(
            config.clone(),
            SAMTOOLS_TAG,
            idxstats_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(PipelineError::Other)?;

    let mut idxstats_stdout = parse_child_output(
        &mut idxstats_child,
        ChildStream::Stdout,
        ParseMode::Lines,
        &config,
    )
        .await
        .map_err(PipelineError::Other)?;

    let mut rgi_read_counts = HashMap::<String, RgiReadCounts>::new();

    while let Some(item) = idxstats_stdout.recv().await {
        let ParseOutput::Bytes(line_bytes) = item else {
            return Err(PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools idxstats expected line-oriented Bytes output"
            )));
        };

        let line = std::str::from_utf8(&line_bytes)
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools idxstats emitted non-UTF8 output: {}",
                    e
                ))
            })?
            .trim_end_matches(['\r', '\n']);

        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() != 4 {
            return Err(PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools idxstats expected 4 tab-separated fields, got {}: {}",
                fields.len(),
                line
            )));
        }

        let reference_length = fields[1]
            .parse::<u64>()
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma invalid idxstats reference length '{}': {}",
                    fields[1],
                    e
                ))
            })?;

        let mapped = fields[2]
            .parse::<u64>()
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma invalid idxstats mapped count '{}': {}",
                    fields[2],
                    e
                ))
            })?;

        let unmapped = fields[3]
            .parse::<u64>()
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma invalid idxstats unmapped count '{}': {}",
                    fields[3],
                    e
                ))
            })?;

        // Match RGI get_reads_count(): retain only references with at
        // least one mapped read. The '*' idxstats row therefore drops
        // out naturally because its mapped count is zero.
        if mapped > 0 {
            rgi_read_counts.insert(
                fields[0].to_string(),
                RgiReadCounts {
                    reference_length,
                    mapped,
                    unmapped,
                    all: mapped + unmapped,
                },
            );
        }
    }

    let idxstats_status = idxstats_child
        .wait()
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools idxstats wait failed: {}",
                e
            ))
        })?;

    idxstats_stderr_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools idxstats stderr task panicked: {}",
                e
            ))
        })?
        .map_err(PipelineError::Other)?;

    if !idxstats_status.success() {
        return Err(PipelineError::Other(anyhow!(
            "RunRgiBwtKma samtools idxstats failed with exit status {}",
            idxstats_status
        )));
    }

    info!(
        "[AMR] RunRgiBwtKma idxstats parsed: {} references with mapped reads",
        rgi_read_counts.len()
    );

    for (reference, counts) in &rgi_read_counts {
        info!(
            "[AMR] RGI idxstats: reference={} length={} mapped={} unmapped={} all={}",
            reference,
            counts.reference_length,
            counts.mapped,
            counts.unmapped,
            counts.all,
        );
    }



    // ------------------------------------------------------------
    // Pull alignment records from the completed AMR BAM and build
    // the in-memory alignment representation used by RGI reporting.
    //
    // RGI:
    //
    //   samtools view --threads <threads> <bam>
    //       | cut -f 1,2,3,4,5,7
    //       | sort -s -n -k 1,1
    //
    // followed by preload_alignments(), which groups the selected
    // SAM fields by RNAME.
    //
    // We eliminate both cut and sort. samtools view stdout is parsed
    // one SAM line at a time and the six fields retained by RGI are
    // extracted directly:
    //
    //   QNAME, FLAG, RNAME, POS, MAPQ, RNEXT
    //
    // Alignments are grouped directly by RNAME, so no temporary
    // sr_amr_report.seqs.temp.txt file is required.
    // ------------------------------------------------------------

    let view_threads =
        config.thread_allocation(SAMTOOLS_TAG, Some("view"));

    let view_config =
        SamtoolsConfig::new(SamtoolsSubcommand::View)
            .option("-@", view_threads.to_string())
            .positional(
                sorted_bam_path
                    .to_string_lossy()
                    .to_string(),
            );

    let view_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&view_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma samtools view args: {:?}",
        view_args
    );

    let (mut view_child, view_stderr_task) =
        spawn_cmd(
            config.clone(),
            SAMTOOLS_TAG,
            view_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(PipelineError::Other)?;

    // ParseMode::Lines emits one complete SAM line per
    // ParseOutput::Bytes, with trailing CR/LF removed.
    let mut view_rx =
        parse_child_output(
            &mut view_child,
            ChildStream::Stdout,
            ParseMode::Lines,
            &config,
        )
            .await
            .map_err(PipelineError::Other)?;

    let mut rgi_alignments:
        HashMap<String, Vec<RgiAlignment>> = HashMap::new();

    let mut alignment_count: u64 = 0;

    while let Some(item) = view_rx.recv().await {
        let bytes = match item {
            ParseOutput::Bytes(bytes) => bytes,

            other => {
                return Err(PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools view produced unexpected \
                     parse output: {:?}",
                    other
                )));
            }
        };

        let line =
            std::str::from_utf8(&bytes)
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "RunRgiBwtKma samtools view emitted \
                         non-UTF8 SAM line: {}",
                        e
                    ))
                })?;

        let alignment =
            parse_rgi_alignment(line)
                .map_err(PipelineError::Other)?;

        alignment_count += 1;

        rgi_alignments
            .entry(alignment.rname.clone())
            .or_default()
            .push(alignment);
    }

    // stdout is completely drained at this point, so samtools can
    // safely be awaited without risking a blocked stdout pipe.
    let view_status =
        view_child
            .wait()
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools view wait failed: {}",
                    e
                ))
            })?;

    view_stderr_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools view stderr task panicked: {}",
                e
            ))
        })?
        .map_err(PipelineError::Other)?;

    if !view_status.success() {
        return Err(PipelineError::Other(anyhow!(
            "RunRgiBwtKma samtools view failed with exit status {}",
            view_status
        )));
    }

    info!(
        "[AMR] RunRgiBwtKma samtools view parsed: \
         {} alignment records across {} reference names",
        alignment_count,
        rgi_alignments.len()
    );

    for (rname, alignments) in &rgi_alignments {
        info!(
            "[AMR] RGI alignments: reference={} records={}",
            rname,
            alignments.len()
        );
    }

    // ------------------------------------------------------------
    // Compute the coverage histogram used by RGI reporting.
    //
    // IMPORTANT FILE BOUNDARY:
    //
    // sr_amr_report.sorted.length_100.bam is fully materialized
    // before execution can reach this point. rgi_bwt_bam_fut awaited
    // sorted_bam_write_task to EOF and returned sorted_bam_path only
    // after that writer completed successfully.
    //
    // The BAM has also already been successfully indexed and consumed
    // by samtools idxstats and samtools view above.
    //
    // RGI:
    //
    //   bedtools genomecov
    //       -ibam sr_amr_report.sorted.length_100.bam
    //
    // RGI writes genomecov stdout to
    // sr_amr_report.coverage_all_positions.temp.txt. We do not need
    // that temporary file boundary, so keep stdout as a line stream.
    // Parsing/aggregation of the histogram is the following stage.
    // ------------------------------------------------------------

    if !sorted_bam_path.is_file() {
        return Err(PipelineError::Other(anyhow!(
            "RunRgiBwtKma BEDTools genomecov input BAM is missing: {}",
            sorted_bam_path.display()
        )));
    }

    let genomecov_config =
        BedtoolsConfig::new(BedtoolsSubcommand::Genomecov)
            .option(
                "-ibam",
                sorted_bam_path
                    .to_string_lossy()
                    .to_string(),
            );

    let genomecov_args = generate_cli(
        BEDTOOLS_TAG,
        &config,
        Some(&genomecov_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunRgiBwtKma bedtools genomecov args: {:?}",
        genomecov_args
    );

    let (mut genomecov_child, genomecov_stderr_task) =
        spawn_cmd(
            config.clone(),
            BEDTOOLS_TAG,
            genomecov_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(PipelineError::Other)?;

    // genomecov emits a tab-delimited coverage histogram on stdout.
    // Keep it line-oriented so the next stage can parse each histogram
    // record directly without materializing RGI's temporary text file.
    let mut genomecov_rx =
        parse_child_output(
            &mut genomecov_child,
            ChildStream::Stdout,
            ParseMode::Lines,
            &config,
        )
            .await
            .map_err(PipelineError::Other)?;

    // Sum histogram base counts, not the rounded fraction column.
    // Keep only idxstats references with mapped reads, as RGI reporting
    // does. BEDTools' aggregate "genome" rows are not reference rows.
    let mut rgi_coverage = HashMap::<String, RgiCoverage>::new();
    while let Some(item) = genomecov_rx.recv().await {
        let ParseOutput::Bytes(bytes) = item else {
            return Err(PipelineError::Other(anyhow!(
                "RunRgiBwtKma genomecov expected line-oriented Bytes output"
            )));
        };
        let line = std::str::from_utf8(&bytes)
            .map_err(|e| PipelineError::Other(anyhow!(
                "RunRgiBwtKma genomecov emitted non-UTF8 output: {}", e
            )))?
            .trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 5 {
            return Err(PipelineError::Other(anyhow!(
                "RunRgiBwtKma genomecov expected 5 fields: {:?}", line
            )));
        }
        if fields[0] == "genome" {
            continue;
        }
        let Some(counts) = rgi_read_counts.get(fields[0]) else {
            continue;
        };
        let parse_integer = |value: &str| -> Result<u64, PipelineError> {
            value.parse::<u64>().map_err(|e| PipelineError::Other(anyhow!(
                "RunRgiBwtKma invalid genomecov integer {:?} in {:?}: {}",
                value, line, e
            )))
        };
        let depth = parse_integer(fields[1])?;
        let bases = parse_integer(fields[2])?;
        let length = parse_integer(fields[3])?;
        if length == 0 || length != counts.reference_length {
            return Err(PipelineError::Other(anyhow!(
                "RunRgiBwtKma genomecov length mismatch for {}: histogram={}, idxstats={}",
                fields[0], length, counts.reference_length
            )));
        }
        let coverage = rgi_coverage.entry(fields[0].to_string())
            .or_insert_with(|| RgiCoverage {
                reference_length: length,
                ..Default::default()
            });
        let bucket = if depth == 0 {
            &mut coverage.uncovered_bases
        } else {
            &mut coverage.covered_bases
        };
        *bucket = bucket.checked_add(bases).ok_or_else(|| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma genomecov base-count overflow for {}", fields[0]
            ))
        })?;
    }

    // Drain stdout before waiting, so BEDTools cannot block on its pipe.
    let genomecov_status = genomecov_child.wait().await
        .map_err(|e| PipelineError::Other(anyhow!(
            "RunRgiBwtKma genomecov wait failed: {}", e
        )))?;
    genomecov_stderr_task.await
        .map_err(|e| PipelineError::Other(anyhow!(
            "RunRgiBwtKma genomecov stderr task panicked: {}", e
        )))?
        .map_err(PipelineError::Other)?;
    if !genomecov_status.success() {
        return Err(PipelineError::Other(anyhow!(
            "RunRgiBwtKma genomecov failed with exit status {}", genomecov_status
        )));
    }

    // Require a complete histogram for every reportable reference.
    // A fully covered reference legitimately has no depth-zero row.
    let mut coverage_references: Vec<&String> = rgi_read_counts.keys().collect();
    coverage_references.sort();
    for reference in coverage_references {
        let coverage = rgi_coverage.get_mut(reference).ok_or_else(|| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma genomecov missing mapped reference {}", reference
            ))
        })?;
        if coverage.covered_bases.checked_add(coverage.uncovered_bases)
            != Some(coverage.reference_length)
        {
            return Err(PipelineError::Other(anyhow!(
                "RunRgiBwtKma incomplete genomecov histogram for {}: covered={}, uncovered={}, length={}",
                reference, coverage.covered_bases, coverage.uncovered_bases,
                coverage.reference_length
            )));
        }
        coverage.percent_coverage =
            100.0 * coverage.covered_bases as f64 / coverage.reference_length as f64;
        info!(
            "[AMR] RGI coverage: reference={} length={} covered={} uncovered={} coverage={:.2}%",
            reference, coverage.reference_length, coverage.covered_bases,
            coverage.uncovered_bases, coverage.percent_coverage
        );
    }
    info!(
        "[AMR] RunRgiBwtKma genomecov complete: {} reportable references",
        rgi_coverage.len()
    );

    // KMA stdout must be drained concurrently with KMA execution.
    // At this point wait for both sides to finish.
    let kma_sam_results = kma_sam_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
    "RunKMA SAM parser task panicked: {}",
    e
            ))
        })?
        .map_err(PipelineError::Other)?;

    sam_lines_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
    "RunKMA SAM line adapter task panicked: {}",
    e
            ))
        })?
        .map_err(PipelineError::Other)?;

    sam_fanout_done
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
    "RunKMA SAM fanout completion receiver dropped: {}",
    e
            ))
        })?
        .map_err(PipelineError::Other)?;



    info!(
            "[AMR] RunKMA native results: records={}, mapped={}, references={}",
            kma_sam_results.total_records,
            kma_sam_results.alignments.len(),
            kma_sam_results.references.len(),
        );

    for alignment in &kma_sam_results.alignments {
        info!(
            "[AMR] KMA mapped: qname={} flag={} rname={} pos={} mapq={} cigar={}",
            alignment.qname,
            alignment.flag,
            alignment.rname,
            alignment.pos,
            alignment.mapq,
            alignment.cigar,
        );
    }



    kma_completion_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
    "RunKMA completion task panicked: {}",
    e
            ))
        })?
        .map_err(PipelineError::Other)?;


    // KMA has completed successfully; shared scratch still owns all sidecars.
    // Keep these maps in memory for subsequent allele-summary construction.
    let kma_sidecars = read_kma_sidecars(config.clone(), &kma_output_prefix)
        .await.map_err(PipelineError::Other)?;
    let missing_depth = rgi_read_counts.keys()
        .filter(|r| !kma_sidecars.read_depth.contains_key(*r)).count();
    let missing_consensus = rgi_read_counts.keys()
        .filter(|r| !kma_sidecars.consensus_dna.contains_key(*r)).count();
    info!(
        "[AMR] KMA sidecar mapped-reference check: mapped_references={} missing_depth={} missing_consensus={} mutation_references={}",
        rgi_read_counts.len(), missing_depth, missing_consensus, kma_sidecars.mutations.len()
    );
    // Missing per-reference evidence is not fatal here: BWT.py handles it
    // during summary construction.

    let card_json = config
        .args
        .amr_card_json
        .as_deref()
        .ok_or_else(|| {
            PipelineError::InvalidConfig(
                "--amr-card-json is required for AMR".to_string()
            )
        })?;

    let card_json =
        resolve_existing_input_path(card_json, &cwd, true)?;

    let card_models =
        load_card_models(&card_json)
            .await
            .map_err(PipelineError::Other)?;

    let mut resolved = 0usize;
    let mut unresolved = Vec::new();

    for reference in rgi_read_counts.keys() {
        if resolve_card_model(reference, &card_models).is_some() {
            resolved += 1;
        } else {
            unresolved.push(reference.clone());
        }
    }

    info!(
    "[AMR] CARD model resolution: mapped_references={} resolved={} unresolved={}",
    rgi_read_counts.len(),
    resolved,
    unresolved.len()
);

    for reference in &unresolved {
        warn!(
        "[AMR] CARD model unresolved: reference={}",
        reference
    );
    }

    // Reproduce BWT.py summary() for the canonical CARD/KMA evidence available
    // at this checkpoint, including consensus translation and depth suppression.
    let rgi_allele_summaries = build_rgi_allele_summaries(
        &rgi_read_counts,
        &rgi_coverage,
        &rgi_alignments,
        &kma_sidecars,
        &card_models,
    );

    let _allele_json_path = write_rgi_allele_json(
        &rgi_allele_summaries,
        &out_dir,
    )
        .await
        .map_err(PipelineError::Other)?;

    for result in &rgi_allele_summaries {
        let RgiAlleleSummaryResult::Summary(summary) = result else {
            continue;
        };

        info!(
            "[AMR] RGI allele summary: reference={} model={} aro={} mapped={} coverage={} mapq={} depth={} snps={} consensus_bp={} mate_pairs={}",
            summary.id,
            summary.cvterm_name,
            summary.aro_accession,
            summary.reads.mapped,
            summary.percent_coverage.covered,
            summary.mapq_average,
            summary.depth,
            if summary.snps.is_empty() { 0 } else { summary.snps.split("; ").count() },
            summary.consensus_sequence_dna.len(),
            summary.mate_pair.len(),
        );
    }

    preparation_task.await
        .map_err(|e| PipelineError::Other(anyhow!(
            "AMR input preparation task panicked: {}", e
        )))??;

    info!("Finished AMR processing through shared genomecov and KMA sidecar evidence collection.");
    Ok(())
}




#[cfg(test)]
mod kma_sidecar_tests {
    use super::*;

    #[test]
    fn res_preserves_depth_text_and_last_duplicate() {
        let mut depths = HashMap::new();
        parse_kma_res_line("#Template\tignored", &mut depths).unwrap();
        parse_kma_res_line("ref\t1\t2\t3\t4\t5\t6\t7\t 0.50 \t9", &mut depths).unwrap();
        assert_eq!(depths["ref"], "0.50");
        parse_kma_res_line("ref\t1\t2\t3\t4\t5\t6\t7\t5.00", &mut depths).unwrap();
        assert_eq!(depths["ref"], "5.00");
        assert!(parse_kma_res_line("bad", &mut depths).is_err());
    }

    #[test]
    fn vcf_preserves_order_duplicates_and_unfiltered_alt() {
        let mut mutations = HashMap::new();
        parse_kma_vcf_line("##fileformat=VCFv4.2", &mut mutations).unwrap();
        for _ in 0..2 {
            parse_kma_vcf_line("ref\t 12 \t.\t A \t C,G \t0\tFAIL", &mut mutations).unwrap();
        }
        parse_kma_vcf_line("ref\t13\t.\tT\t.", &mut mutations).unwrap();
        assert_eq!(mutations["ref"], vec!["A12C,G", "A12C,G", "T13."]);
        assert!(parse_kma_vcf_line("bad", &mut mutations).is_err());
    }

    #[tokio::test]
    async fn fsa_multiline_case_gaps_and_duplicate_ids() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fsa");
        tokio::fs::write(&path, b">ref description\nAc g\nT\t-N\n>empty\n").await.unwrap();
        let sequences = read_kma_fsa(&path).await.unwrap();
        assert_eq!(sequences["ref"], "AcgT-N");
        assert_eq!(sequences["empty"], "");
        tokio::fs::write(&path, b">ref one\nA\n>ref two\nC\n").await.unwrap();
        assert!(read_kma_fsa(&path).await.is_err());
    }
}

