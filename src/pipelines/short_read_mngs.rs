//! Short Read Metagenomics Pipeline

use std::cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader};
use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};

use ahash::AHashMap;
use ahash::RandomState as AHashRandomState;
use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use dashmap::DashMap;
use fst::Map;
use futures::future::try_join_all;
use log::{self, debug, error, info, warn};
use needletail::parse_fastx_file;
use noodles::bam::r#async::io::Reader as BamAsyncReader;
use noodles::bam::record::Record;
use noodles::sam::alignment::record::cigar::op::Kind as OpKind;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand_core::RngCore;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use tempfile::NamedTempFile;
use tempfile::TempDir;
use tokio::fs;
use tokio::fs::File as TokioFile;
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader as TokioBufReader, BufWriter};
use tokio::process::Child;
use tokio::process::Command;
use tokio::sync::mpsc;
use tokio::sync::mpsc::Sender;
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio::time::{sleep, Duration, Instant};
use tokio_stream::wrappers::{ReceiverStream, UnboundedReceiverStream};
use tokio_stream::StreamExt;
use tokio_util::io::StreamReader;
use twox_hash::XxHash64;
use fst::Streamer;

use crate::config::defs::{
    ClusterInfo, DiamondSubcommand, DuplicateClusters, KallistoSubcommand, Lineage,
    NRAlignmentBackend, PairingMode, PipelineError, ReadCountingMode, ReadStats, RunConfig,
    SamtoolsSubcommand, StreamDataType, Taxid, BLASTN_TAG, BLASTX_TAG, BOWTIE2_TAG, DIAMOND_TAG,
    FASTP_TAG, HISAT2_TAG, KALLISTO_TAG, MAKEBLASTDB_TAG, MINIMAP2_TAG, MMSEQS_TAG, NR_TAG, NT_TAG,
    READ_COUNTING_MODE, SAMTOOLS_TAG, SORT_TAG, SPADES_TAG, ReducedRead, ExecutionMode
};
use crate::utils::align::{
    sharded_minimap2_align,
    Minimap2Preset,
    ShardedMinimap2Config,
    run_diamond_single_file,
};
use crate::utils::blast::{
    build_taxon_counts_list, call_hits_m8, compute_merged_taxon_counts,
    generate_taxon_count_json_from_m8, generate_taxon_counts, merge_aggregations,
    process_record_pair, read_id_from_m8_line, shard_for_read_id, sort_m8_by_read_id,
    summarize_m8_hits, combine_taxon_counts, AggBucket, ContigSummaryEntry, M8Record, SpeciesAlignmentResults,
    TaxonCount,
};
use crate::utils::command::blastn::BlastnConfig;
use crate::utils::command::blastx::BlastxConfig;
use crate::utils::command::bowtie2::{bowtie2_index_prep, Bowtie2Config};
use crate::utils::command::diamond::{
    compute_optimal_block_size, compute_diamond_index_chunks, diamond_index_prep, DiamondConfig,
};
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::hisat2::{hisat2_index_prep, Hisat2Config};
use crate::utils::command::kallisto::KallistoConfig;
use crate::utils::command::makeblastdb::MakeblastdbConfig;
use crate::utils::command::minimap2::{minimap2_index_prep, Minimap2Config};
use crate::utils::command::mmseqs::{MmseqsBackend, MmseqsConfig, MmseqsSubcommand, build_mmseqs_cpu_search_config, build_mmseqs_gpu_search_config};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::coverage_viz::generate_coverage_viz;
use crate::utils::fastx::{
    generate_taxid_fasta, generate_taxid_locator, parse_byte_stream_to_fastq, raw_read_count,
    read_fasta, read_fastq, stream_record_counter, write_combined_fastq,
    write_fasta_stream_to_file, SequenceRecord, chunk_paired_fastq, generate_annotated_fasta_stream,
    UNMAPPED_HEADER_PREFIX, StreamRecordStats, collect_hit_summary_to_accession_map_concurrent
};
use crate::utils::file::{
    choose_temp_dir, file_path_manipulator, file_size, rename_file_path, resolve_optional_path,
    validate_file_inputs, write_byte_stream_to_file, write_parse_output_to_file,
    write_vecu8_to_file, copy_file_streaming,
};
use crate::utils::paf::paf_to_m8;
use crate::utils::plotting::plot_insert_sizes;
use crate::utils::sambam::{compute_insert_size_stats_from_bam, generate_info_from_bam_stream, generate_assembly_coverage};
use crate::utils::streams::deinterleave_fastq_stream_to_fifos;
use crate::utils::streams::{
    batch_rayon_process, deinterleave_fastq_stream, fanout_to_channels, join_with_error_handling,
    parse_bytes, parse_child_output, parse_fastq, parse_lines,
    read_child_output_to_vec, spawn_cmd, stream_to_cmd, stream_to_file, ChannelReader, ChildStream,
    ParseMode, ParseOutput, ToBytes,
};
use crate::utils::subsample::subsample_uniform_reservoir;
use crate::utils::system::{compute_batch_size, compute_phase_concurrency};
use crate::utils::taxonomy::{
    build_should_keep_filter, get_top_m8_nr, get_top_m8_nt, load_lineage_and_acc2tax_maps,
};

const MAX_SPADES_WORK_DIR: u64 = 500_000_000;

const MAX_NUM_BINS_COVERAGE: usize = 500;
const NUM_ACCESSIONS_PER_TAXON: usize = 10;
const MIN_CONTIG_SIZE: u64 = 4;

const MIN_REF_FASTA_SIZE: u64 = 25;
const MIN_ASSEMBLED_CONTIG_SIZE: u64 = 25;

const REFERENCE_METADATA_BUCKET: &str = "seqtoid-public-references";

const REFERENCE_METADATA_PREFIX: &str = "phase2/refs/metadata";

const DISTRIBUTED_WORKER_STALE_SECS: u64 = 30;

/// Holds the results from a Kallisto quantification run.
#[derive(Debug)]
pub struct KallistoResults {
    /// ERCC counts as a vector of (target_id, est_counts).
    pub ercc_counts: Vec<(String, f64)>,
    /// Mapping from transcript_id to gene_id.
    pub transcript_to_gene: Vec<(String, String)>,
}

/// Represents an item in a reservoir sample.
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct SampleItem {
    /// Key used for A-Res sampling.
    key: f64,
    /// Sequence records (1 for single-end, 2 for paired-end).
    records: Vec<SequenceRecord>,
    /// Number of duplicate reads this item represents.
    weight: u64,
    /// Optional hash key for disk-backed sampling.
    hash_key: Option<u64>,
}

impl PartialEq for SampleItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}
impl Eq for SampleItem {}
impl Ord for SampleItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.partial_cmp(&other.key).unwrap_or(Ordering::Equal)
    }
}
impl PartialOrd for SampleItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Combined taxonomic counts from both NT and NR databases.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CombinedTaxonCount {
    pub tax_id: i32,
    pub tax_level: u8,
    pub genus_taxid: i32,
    pub family_taxid: i32,
    pub nt: Option<TaxonMetrics>,
    pub nr: Option<TaxonMetrics>,
}

/// Detailed metrics for a taxon hit in a specific database.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TaxonMetrics {
    /// Adjusted read count (with duplicate compression).
    pub count: u64,
    pub nonunique_count: u64,
    pub unique_count: u64,
    /// Duplicate compression ratio.
    pub dcr: f64,
    pub percent_identity: f64,
    pub alignment_length: f64,
    pub e_value: f64,
    pub base_count: u64,
    pub source_count_type: Option<Vec<String>>,
}

/// Handle for managing an ongoing Spades assembly task.
#[derive(Debug)]
pub struct AssemblyHandle {
    pub spades_task: JoinHandle<Result<(), anyhow::Error>>,
    pub work_dir: PathBuf,
    pub out_dir: PathBuf,
}

/// Paths and metadata for assembly coverage analysis outputs.
#[derive(Debug)]
pub struct CoverageOutputs {
    pub contigs_fasta: PathBuf,
    pub contigs_all_fasta: PathBuf,
    pub scaffolds_fasta: PathBuf,
    pub bam_path: PathBuf,
    pub contig_stats_json: PathBuf,
    pub contig_stats: Arc<HashMap<String, u64>>,
    pub read2contig: Arc<HashMap<String, String>>,
}

/// Information about a single read's taxonomic hit.
#[derive(Debug, Clone)]
pub struct ReadHit {
    pub level: u8,
    pub taxid: i32,
    pub accession_id: String,
    pub species_taxid: i32,
    pub genus_taxid: i32,
    pub family_taxid: i32,
    pub contig_id: Option<String>,
    pub contig_accession_id: Option<String>,
    pub contig_species_taxid: i32,
    pub contig_genus_taxid: i32,
    pub contig_family_taxid: i32,
    pub from_assembly: bool,
    pub source_count_type: Option<String>,
}

impl ReadHit {
    /// Serializes the read hit to a tab-separated string.
    ///
    /// # Returns
    ///
    /// String: tab-separated representation of the hit
    pub fn to_tab_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.level,
            self.taxid,
            &self.accession_id,
            "", // alignment length
            "", // percent identity
            "", // bitscore
            "", // evalue
            self.species_taxid,
            self.genus_taxid,
            self.family_taxid,
            "", // name
            "", // common_name
            self.contig_id.as_deref().unwrap_or(""),
            self.contig_accession_id.as_deref().unwrap_or(""),
            self.contig_species_taxid,
            self.contig_genus_taxid,
            self.contig_family_taxid,
            "", // contig_length
            "", // contig_name
            "", // contig_accession
            "", // contig_taxid
            "", // contig_species_taxid duplicate
            "", // contig_genus_taxid duplicate
            "", // contig_family_taxid duplicate
            if self.from_assembly { "1" } else { "0" },
            "1", // count
            self.source_count_type.as_deref().unwrap_or("")
        )
    }

    /// Serializes the read hit to a tab-separated string with an explicit read ID.
    ///
    /// # Arguments
    ///
    /// * `read_id`: the ID of the read
    ///
    /// # Returns
    ///
    /// String: tab-separated representation of the hit
    /// HitSummary line: read_id, level, taxid, accession, species, genus, family
    /// Optional contig fields appended when present (hitsummary2-compatible).
    pub fn to_full_tab_line(&self, read_id: &str) -> String {
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            read_id,
            self.level,
            self.taxid,
            self.accession_id,
            self.species_taxid,
            self.genus_taxid,
            self.family_taxid,
        );
        if self.contig_id.is_some() || self.from_assembly {
            line.push_str(&format!(
                "\t{}\t{}\t{}\t{}\t{}",
                self.contig_id.as_deref().unwrap_or(""),
                self.contig_accession_id.as_deref().unwrap_or(""),
                self.contig_species_taxid,
                self.contig_genus_taxid,
                self.contig_family_taxid,
            ));
        }
        line
    }

    /// Parses a read hit from a tab-separated line.
    ///
    /// # Arguments
    ///
    /// * `line`: the tab-separated line to parse
    ///
    /// # Returns
    ///
    /// Result<(String, Self)>: the read ID and the parsed ReadHit
    pub fn from_tab_line(line: &str) -> Result<(String, Self)> {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 28 {
            return Err(anyhow!(
                "Invalid hit summary line: expected 28 fields, got {}",
                fields.len()
            ));
        }

        let read_id = fields[0].to_string();
        let hit = Self {
            level: fields[1].parse()?,
            taxid: fields[2].parse()?,
            accession_id: fields[3].to_string(),
            species_taxid: fields[8].parse()?,
            genus_taxid: fields[9].parse()?,
            family_taxid: fields[10].parse()?,
            contig_id: if fields[13].is_empty() {
                None
            } else {
                Some(fields[13].to_string())
            },
            contig_accession_id: if fields[14].is_empty() {
                None
            } else {
                Some(fields[14].to_string())
            },
            contig_species_taxid: fields[15].parse()?,
            contig_genus_taxid: fields[16].parse()?,
            contig_family_taxid: fields[17].parse()?,
            from_assembly: fields[25] == "1",
            source_count_type: if fields[27].is_empty() {
                None
            } else {
                Some(fields[27].to_string())
            },
        };

        Ok((read_id, hit))
    }
}

#[derive(Default, Clone, Debug)]
pub struct AccessionHit {
    pub taxid: i32, // i.e. species_taxid
    pub genus_taxid: i32,
    pub family_taxid: i32,
    pub count: u64,
}

/// Represents a sequence cluster from a deduplication step.
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct Cluster {
    /// ID of the representative read for this cluster.
    rep_id: String,
    /// Number of reads in the cluster.
    size: u64,
    /// IDs of all reads belonging to this cluster.
    members: Vec<String>,
}


/// Validates the input FASTQ files and returns an interleaved stream of sequence records.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `file1_path`: path to the R1 (or single-end) FASTQ file
/// * `file2_path`: optional path to the R2 FASTQ file
/// * `sample_base_buf`: base path for the sample
/// * `out_dir`: directory where output files should be saved
///
/// # Returns
///
/// Result containing:
/// - interleaved FASTQ stream
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - handle for the raw read count task
/// - handle for the validated read count task
async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> anyhow::Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
        JoinHandle<anyhow::Result<ReadStats, anyhow::Error>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated.fq.gz"),
        "_",
    );


    let (rx, val_count_task) = read_fastq(
        file1_path,
        file2_path,
        PairingMode::Strict,
        Some(config.args.technology.clone()),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        "validate_input",
        &config,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);

    // Split the byte stream for fastp and write
    let (val_streams, val_router_done_rx) = fanout_to_channels(
        val_rx_stream,
        2,
        "validate_input",
        &config,
        StreamDataType::IlluminaFastq,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(val_router_done_rx);

    let mut val_streams_iter = val_streams.into_iter();

    let val_ercc_bowtie2_filter_out_stream =
        ReceiverStream::new(val_streams_iter.next().ok_or(PipelineError::EmptyStream)?);

    let val_file_stream =
        ReceiverStream::new(val_streams_iter.next().ok_or(PipelineError::EmptyStream)?);

    let val_file_write_task = write_byte_stream_to_file(
        &validated_interleaved_file_path,
        val_file_stream,
        config.clone(),
        StreamDataType::IlluminaFastq,
        "validate_input",
        false
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(val_file_write_task);

    Ok((
        val_ercc_bowtie2_filter_out_stream,
        cleanup_tasks,
        cleanup_receivers,
        val_count_task,
    ))
}

/// Aligns reads using Bowtie2 and sorts the output BAM stream.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `bt2_index_path`: path to the Bowtie2 index
/// * `paired`: whether the input is paired-end
/// * `bowtie2_options`: additional command-line options for Bowtie2
/// * `output_bam_path`: path where the sorted BAM file will be written
///
/// # Returns
///
/// Result containing:
/// - stream of BAM records
/// - BAM write task handle
/// - receiver for the count of mapped reads
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - temporary directory used for sorting
async fn bowtie2_align_and_sort_stream(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    output_bam_path: PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>, // BAM → for fastq extraction
        JoinHandle<Result<()>>,      // BAM write task
        oneshot::Receiver<u64>,      // count rx
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        TempDir,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    // ─────────────────────────────
    // Bowtie2
    // ─────────────────────────────

    let mut bt2_config = Bowtie2Config::new()
        .flag("--very-sensitive-local")
        .option(
            "-x",
            bt2_index_path.to_string_lossy().to_string(),
        )
        .option(
            "-p",
            config
                .thread_allocation(BOWTIE2_TAG, None)
                .to_string(),
        );

    bt2_config = if paired {
        bt2_config
            .flag("--interleaved")
            .positional("-")
    } else {
        bt2_config.option("-U", "-")
    };

    let bt2_args =
        generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config))?;

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config))?;

    let (bt2_child, bt2_stream_task, bt2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        BOWTIE2_TAG,
        bt2_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await?;

    cleanup_tasks.push(bt2_stream_task);
    cleanup_tasks.push(bt2_err_task);

    let bt2_out = {
        let mut guard = bt2_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config).await?
    };

    // ─────────────────────────────
    // samtools sort -n, needs name sorting before conversion
    // ─────────────────────────────

    let sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("sort"))
                .to_string(),
        )
        .option("-m", "2G")
        .flag("-n") // name-sorted
        .flag("-u") // uncompressed
        .option("-O", "bam")
        .option("-T", temp_dir.path().to_string_lossy().to_string())
        .positional("-");

    let sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&sort_config))?;

    let (sort_child, sort_task, sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out,
        SAMTOOLS_TAG,
        sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await?;

    cleanup_tasks.push(sort_task);
    cleanup_tasks.push(sort_err_task);

    let sorted_bam_stream = {
        let mut guard = sort_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config).await?
    };

    // ─────────────────────────────
    // split stream
    // ─────────────────────────────

    let (streams, router_done_rx) = fanout_to_channels(
        ReceiverStream::new(sorted_bam_stream),
        3,
        "bt2_split",
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    // track the task instead of a oneshot receiver
    cleanup_receivers.push(router_done_rx);

    let mut it = streams.into_iter();

    let fastq_stream = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
    let count_stream = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
    let bam_write_stream = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);

    // ─────────────────────────────
    // BAM write
    // ─────────────────────────────

    let bam_write_task = tokio::spawn({
        let config = config.clone();
        async move {
            write_byte_stream_to_file(
                &output_bam_path,
                bam_write_stream,
                config.clone(),
                StreamDataType::JustBytes,
                "bt2_bam_write",
                false
            )
                .await?;
            Ok(())
        }
    });

    // ─────────────────────────────
    // count
    // ─────────────────────────────

    let (count_tx, count_rx) = oneshot::channel();
    let count_task = tokio::spawn(async move {
        let mut s = count_stream;
        let mut lines = 0u64;
        while let Some(item) = s.next().await {
            if let ParseOutput::Bytes(b) = item {
                if !b.is_empty() {
                    lines += 1;
                }
            }
        }
        let _ = count_tx.send(lines);
        Ok(())
    });

    cleanup_tasks.push(count_task);

    Ok((
        fastq_stream,
        bam_write_task,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        temp_dir,
    ))
}

/// Filters reads using Bowtie2 and returns a stream of unmapped FASTQ records.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `bt2_index_path`: path to the Bowtie2 index
/// * `paired`: whether the input is paired-end
/// * `bowtie2_options`: additional command-line options for Bowtie2
/// * `output_bam_path`: path where the sorted BAM file will be written
///
/// # Returns
///
/// Result containing:
/// - stream of unmapped FASTQ records
/// - receiver for the count of mapped reads
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - BAM write task handle
/// - path to the generated BAM file
/// - vector of temporary directories used
async fn bowtie2_filter_stream(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    output_bam_path: PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        oneshot::Receiver<u64>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        JoinHandle<Result<()>>,
        PathBuf,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let (
        fastq_input_stream,
        bam_write_task,
        count_rx,
        mut cleanup_tasks,
        cleanup_receivers,
        temp_dir,
    ) = bowtie2_align_and_sort_stream(
        config.clone(),
        input_stream,
        bt2_index_path,
        paired,
        output_bam_path.clone(),
    )
        .await?;

    // samtools fastq → stdout
    let fastq_flag = if paired { "-f13" } else { "-f4" };

    let config_fastq = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("fastq"))
                .to_string(),
        )
        .option("-c", "6")
        .flag(fastq_flag)
        .flag("-N")
        .positional("-");

    let args =
        generate_cli(SAMTOOLS_TAG, &config, Some(&config_fastq))?;

    let (child, stdin_task, err_task) = stream_to_cmd(
        config.clone(),
        fastq_input_stream.into_inner(),
        SAMTOOLS_TAG,
        args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await?;

    cleanup_tasks.push(stdin_task);
    cleanup_tasks.push(err_task);

    let out_stream = {
        let mut guard = child.lock().await;

        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            &config,
        )
            .await?
    };

    Ok((
        ReceiverStream::new(out_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        bam_write_task,
        output_bam_path,
        vec![temp_dir],
    ))
}

/// Filters reads from disk using Bowtie2 and returns paths to the unmapped FASTQ files.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `bt2_index_path`: path to the Bowtie2 index
/// * `paired`: whether the input is paired-end
/// * `bowtie2_options`: additional command-line options for Bowtie2
/// * `output_bam_path`: path where the sorted BAM file will be written
///
/// # Returns
///
/// Result containing:
/// - path to the R1 (or single-end) unmapped FASTQ file
/// - optional path to the R2 unmapped FASTQ file
/// - receiver for the count of mapped reads
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - BAM write task handle
/// - path to the generated BAM file
/// - vector of temporary directories used
async fn bowtie2_filter_files(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    output_bam_path: PathBuf,
) -> Result<
    (
        PathBuf,
        Option<PathBuf>,
        oneshot::Receiver<u64>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        JoinHandle<Result<()>>,
        PathBuf,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let (
        fastq_input_stream,
        bam_write_task,
        count_rx,
        mut cleanup_tasks,
        cleanup_receivers,
        temp_dir,
    ) = bowtie2_align_and_sort_stream(
        config.clone(),
        input_stream,
        bt2_index_path,
        paired,
        output_bam_path.clone(),
    )
        .await?;

    let fastq_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        2,
        false,
    )
        .await?;

    let fq1 = fastq_temp_dir.path().join("bt2_R1.fastq");
    let fq2 = if paired {
        Some(fastq_temp_dir.path().join("bt2_R2.fastq"))
    } else {
        None
    };

    if paired {
        let config_fastq = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
            .option(
                "-@",
                config
                    .thread_allocation(SAMTOOLS_TAG, Some("fastq"))
                    .to_string(),
            )
            .option("-c", "6")
            .option("-f", "13")
            .option("-1", fq1.to_string_lossy().to_string())
            .option(
                "-2",
                fq2.as_ref().unwrap().to_string_lossy().to_string(),
            )
            .flag("-N")
            .option("-0", "/dev/null")
            .option("-s", "/dev/null")
            .positional("-");

        let args = generate_cli(SAMTOOLS_TAG, &config, Some(&config_fastq))?;

        let (child, stdin_task, err_task) = stream_to_cmd(
            config.clone(),
            fastq_input_stream.into_inner(),
            SAMTOOLS_TAG,
            args,
            StreamDataType::JustBytes,
            config.args.verbose,
            None,
        )
            .await?;

        cleanup_tasks.push(stdin_task);
        cleanup_tasks.push(err_task);

        {
            let mut guard = child.lock().await;
            guard.wait().await?;
        }
    } else {
        // streaming fastq (no -o flag)
        let fastq_config = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
            .option(
                "-@",
                config
                    .thread_allocation(SAMTOOLS_TAG, Some("fastq"))
                    .to_string(),
            )
            .option("-c", "6")
            .flag("-f4")
            .flag("-N")
            .positional("-");

        let fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&fastq_config))?;

        let (fastq_child, fastq_task, fastq_err_task) = stream_to_cmd(
            config.clone(),
            fastq_input_stream.into_inner(),
            SAMTOOLS_TAG,
            fastq_args,
            StreamDataType::JustBytes,
            config.args.verbose,
            None,
        )
            .await?;

        cleanup_tasks.extend([fastq_task, fastq_err_task]);

        let fastq_stream = {
            let mut guard = fastq_child.lock().await;
            parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Fastq, &config).await?
        };

        let write_task = write_byte_stream_to_file(
            &fq1,
            ReceiverStream::new(fastq_stream),
            config.clone(),
            StreamDataType::JustBytes,
            "bt2_fastq_se",
            false
        )
            .await?;
        cleanup_tasks.push(write_task);
    }

    Ok((
        fq1,
        fq2,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        bam_write_task,
        output_bam_path,
        vec![temp_dir, fastq_temp_dir],
    ))
}

/// Filters reads using HISAT2 and returns a stream of unmapped FASTQ records.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `r1_path`: path to the R1 (or single-end) FASTQ file
/// * `r2_path_opt`: optional path to the R2 FASTQ file
/// * `hisat2_index_path`: path to the HISAT2 index
/// * `paired`: whether the input is paired-end
/// * `hisat2_options`: additional command-line options for HISAT2
/// * `output_bam_path`: optional path to save the aligned BAM file
/// * `headroom`: estimated memory headroom available for sorting
///
/// # Returns
///
/// Result containing:
/// - stream of unmapped FASTQ records
/// - receiver for the count of mapped reads
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - vector of temporary directories used
async fn hisat2_filter(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    hisat2_index_path: PathBuf,
    paired: bool,
    hisat2_options: Vec<(String, Option<String>)>,
    output_bam_path: Option<PathBuf>,
    headroom: u64,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        oneshot::Receiver<u64>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // PAIRED-END PATH
    if paired {
        let temp_dir = choose_temp_dir(
            config.input_size * headroom,
            &config.ram_temp_dir,
            &config.args.nvme_scratch,
            4,
            false,
        )
            .await?;

        // 1. HISAT2 → SAM
        let sam_path = temp_dir.path().join("hisat2.sam");
        let r2_path = r2_path_opt.as_ref().ok_or_else(|| {
            PipelineError::InvalidConfig(
                "Paired HISAT2 alignment requested but R2 path is missing".to_string(),
            )
        })?;

        let mut hisat2_config = Hisat2Config::new()
            .option("-x", hisat2_index_path.to_string_lossy().to_string())
            .option(
                "-p",
                config
                    .thread_allocation(HISAT2_TAG, None)
                    .to_string(),
            );

        for (key, value) in &hisat2_options {
            hisat2_config = match value {
                Some(value) => hisat2_config.option(key.clone(), value.clone()),
                None => hisat2_config.flag(key.clone()),
            };
        }

        hisat2_config = hisat2_config
            .option("-1", r1_path.to_string_lossy().to_string())
            .option("-2", r2_path.to_string_lossy().to_string())
            .option("-S", sam_path.to_string_lossy().to_string());

        let hisat2_args =
            generate_cli(HISAT2_TAG, &config, Some(&hisat2_config)).map_err(|e| {
                PipelineError::ToolExecution {
                    tool: HISAT2_TAG.to_string(),
                    error: e.to_string(),
                }
            })?;

        let (mut hisat2_child, hisat2_err_task) = spawn_cmd(
            config.clone(),
            HISAT2_TAG,
            hisat2_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: HISAT2_TAG.to_string(),
                error: e.to_string(),
            })?;

        cleanup_tasks.push(hisat2_err_task);

        let status = hisat2_child
            .wait()
            .await
            .map_err(|e| PipelineError::Other(e.into()))?;
        if !status.success() {
            return Err(PipelineError::ToolExecution {
                tool: HISAT2_TAG.to_string(),
                error: "HISAT2 execution failed".to_string(),
            });
        }

        // 2. sort -n → BAM
        let bam_path = temp_dir.path().join("hisat2_sorted.bam");
        let sort_threads = config.thread_allocation(SAMTOOLS_TAG, Some("sort"));

        let sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
            .option("-@", sort_threads.to_string())
            .option("-m", "2G")
            .flag("-n")
            .option("-o", bam_path.to_string_lossy().to_string())
            .option("-T", temp_dir.path().to_string_lossy().to_string())
            .positional(sam_path.to_string_lossy().to_string());

        let sort_args =
            generate_cli(SAMTOOLS_TAG, &config, Some(&sort_config)).map_err(|e| {
                PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                }
            })?;

        let (mut sort_child, sort_err_task) = spawn_cmd(
            config.clone(),
            SAMTOOLS_TAG,
            sort_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

        cleanup_tasks.push(sort_err_task);
        sort_child
            .wait()
            .await
            .map_err(|e| PipelineError::Other(e.into()))?;

        if let Some(out) = output_bam_path {
            tokio::fs::copy(&bam_path, &out)
                .await
                .map_err(|e| PipelineError::IOError(e.to_string()))?;
        }

        // 3. samtools fastq → R1/R2 files
        let fq1_path = temp_dir.path().join("hisat2_filtered_R1.fastq");
        let fq2_path = temp_dir.path().join("hisat2_filtered_R2.fastq");

        let fastq_config = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
            .option(
                "-@",
                config
                    .thread_allocation(SAMTOOLS_TAG, Some("fastq"))
                    .to_string(),
            )
            .option("-c", "6")
            .option("-f", "13")
            .option("-1", fq1_path.to_string_lossy().to_string())
            .option("-2", fq2_path.to_string_lossy().to_string())
            .option("-0", "/dev/null")
            .option("-s", "/dev/null")
            .positional(bam_path.to_string_lossy().to_string());

        let fastq_args =
            generate_cli(SAMTOOLS_TAG, &config, Some(&fastq_config)).map_err(|e| {
                PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                }
            })?;

        let (mut fastq_child, fastq_err_task) = spawn_cmd(
            config.clone(),
            SAMTOOLS_TAG,
            fastq_args,
            config.args.verbose,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?;

        cleanup_tasks.push(fastq_err_task);
        fastq_child
            .wait()
            .await
            .map_err(|e| PipelineError::Other(e.into()))?;

        // 4. convert to stream
        let (rx, read_task) = read_fastq(
            fq1_path.clone(),
            Some(fq2_path.clone()),
            PairingMode::Strict,
            None,
            u64::MAX,
            None,
            None,
            "hisat2_filter",
            &config,
        )
            .map_err(|e| PipelineError::Other(e))?;

        cleanup_tasks.push(tokio::spawn(async move {
            read_task.await??;
            Ok(())
        }));

        let out_stream = ReceiverStream::new(rx);

        // 5. mapped count
        let (count_tx, count_rx) = oneshot::channel();

        let count_config = SamtoolsConfig::new(SamtoolsSubcommand::View)
            .option(
                "-@",
                config
                    .thread_allocation(SAMTOOLS_TAG, Some("view"))
                    .to_string(),
            )
            .flag("--no-PG")
            .flag("-c")
            .flag("-F13")
            .positional(bam_path.to_string_lossy().to_string());

        let count_args =
            generate_cli(SAMTOOLS_TAG, &config, Some(&count_config))
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;

        let count_config_clone = config.clone();

        let count_task = tokio::spawn(async move {
            let (mut child, err_task) = spawn_cmd(
                count_config_clone.clone(),
                SAMTOOLS_TAG,
                count_args,
                count_config_clone.args.verbose,
                None,
            )
                .await?;

            let lines = read_child_output_to_vec(
                &mut child,
                ChildStream::Stdout,
                &count_config_clone,
            )
                .await?;

            let status = child.wait().await?;

            err_task.await??;

            if !status.success() {
                return Err(anyhow!(
            "samtools view count exited with status {}",
            status
        ));
            }

            let count = lines
                .first()
                .map(|line| line.trim())
                .unwrap_or("0")
                .parse::<u64>()
                .unwrap_or(0);

            let _ = count_tx.send(count);

            Ok::<(), anyhow::Error>(())
        });

        cleanup_tasks.push(count_task);

        return Ok((
            out_stream,
            count_rx,
            cleanup_tasks,
            cleanup_receivers,
            vec![temp_dir],
        ));
    }

    let temp_dir = choose_temp_dir(
        config.input_size * headroom,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    // 1. HISAT2 → SAM (streaming)
    let mut hisat2_config = Hisat2Config::new()
        .option("-x", hisat2_index_path.to_string_lossy().to_string())
        .option(
            "-p",
            config
                .thread_allocation(HISAT2_TAG, None)
                .to_string(),
        );

    for (key, value) in &hisat2_options {
        hisat2_config = match value {
            Some(value) => hisat2_config.option(key.clone(), value.clone()),
            None => hisat2_config.flag(key.clone()),
        };
    }

    hisat2_config = hisat2_config
        .option("-U", r1_path.to_string_lossy().to_string())
        .option("-S", "/dev/stdout");

    let hisat2_args =
        generate_cli(HISAT2_TAG, &config, Some(&hisat2_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: HISAT2_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (mut hisat2_child, hisat2_err_task) = spawn_cmd(
        config.clone(),
        HISAT2_TAG,
        hisat2_args,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: HISAT2_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(hisat2_err_task);

    let hisat2_sam_stream = {
        let mut guard = hisat2_child;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config).await?
    };

    // 2. samtools sort -n
    let sort_threads = config.thread_allocation(SAMTOOLS_TAG, Some("sort"));

    let sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option("-@", sort_threads.to_string())
        .option("-m", "2G")
        .flag("-n")
        .flag("-u")
        .option("-O", "bam")
        .option("-T", temp_dir.path().to_string_lossy().to_string())
        .positional("-");

    let sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&sort_config))?;

    let (sort_child, sort_task, sort_err_task) = stream_to_cmd(
        config.clone(),
        hisat2_sam_stream,
        SAMTOOLS_TAG,
        sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await?;

    cleanup_tasks.extend([sort_task, sort_err_task]);

    let sorted_bam_stream = {
        let mut guard = sort_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config).await?
    };

    // 3. Always fanout for count + fastq (reuse existing fanout + splits pattern).
    //    Write tee only if output_bam_path requested. Keeps streaming when None.
    let num_tees = if output_bam_path.is_some() { 3 } else { 2 };
    let (streams, hisat2_bam_split_done_rx) = fanout_to_channels(
        ReceiverStream::new(sorted_bam_stream),
        num_tees,
        "hisat2_bam_split",
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    let mut splits = streams.into_iter();
    let count_stream = splits.next().ok_or(PipelineError::EmptyStream)?;
    let fastq_stream = splits.next().ok_or(PipelineError::EmptyStream)?;

    if let Some(bam_path) = output_bam_path {
        let bam_write_stream = splits.next().ok_or(PipelineError::EmptyStream)?;
        let write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(bam_write_stream),
            config.clone(),
            StreamDataType::JustBytes,
            "hisat2_bam_write",
            false
        )
            .await?;
        cleanup_tasks.push(write_task);
    }
    cleanup_receivers.push(hisat2_bam_split_done_rx);

    let fastq_input_stream = ReceiverStream::new(fastq_stream);

    // 4. samtools fastq -f 4 → stdout (streaming)
    let fastq_threads = config.thread_allocation(SAMTOOLS_TAG, Some("fastq"));

    let fastq_config = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
        .option("-@", fastq_threads.to_string())
        .option("-c", "6")
        .flag("-f4")
        .flag("-N")
        .positional("-");

    let fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&fastq_config))?;

    let (fastq_child, fastq_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        fastq_input_stream.into_inner(),
        SAMTOOLS_TAG,
        fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await?;

    cleanup_tasks.extend([fastq_task, fastq_err_task]);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Fastq, &config).await?
    };

    let out_stream = ReceiverStream::new(unmapped_fastq_stream);

    // 5. mapped count via streaming view -c -F4 on its tee (correct SE count, no data drop)
    let (count_tx, count_rx) = oneshot::channel();
    let samtools_count_config = SamtoolsConfig::new(SamtoolsSubcommand::View)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("view"))
                .to_string(),
        )
        .flag("--no-PG")
        .flag("-c")
        .flag("-F4")
        .positional("-");
    let count_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_count_config))?;

    let (count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        count_stream,
        SAMTOOLS_TAG,
        count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await?;

    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let count_future = tokio::spawn({
        let config = config.clone();
        async move {
            let mut guard = count_child_arc.lock().await;
            let lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout, &config).await?;
            let count_str = lines.first().map(|s| s.trim()).unwrap_or("0");
            let count: u64 = count_str.parse().unwrap_or(0);
            let _ = count_tx.send(count);
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(count_future);

    Ok((
        out_stream,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        vec![temp_dir],
    ))
}


/// Performs quality control on the input FASTQ stream using fastp.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `paired`: whether the input is paired-end
/// * `input_stream`: stream of raw FASTQ records
///
/// # Returns
///
/// Result containing:
/// - stream of quality-controlled FASTQ records
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - receiver for the count of records that passed QC
async fn fastp_qc(
    config: Arc<RunConfig>,
    paired: bool,
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        oneshot::Receiver<Result<u64, anyhow::Error>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let mut qc_fastp_config_view = FastpConfig::new()
        .flag("--stdin")
        .flag("--stdout")
        .option(
            "-j",
            config
                .out_dir
                .join("fastp.json")
                .to_string_lossy()
                .to_string(),
        )
        .option(
            "-h",
            config
                .out_dir
                .join("fastp.html")
                .to_string_lossy()
                .to_string(),
        )
        .option(
            "-w",
            config
                .thread_allocation(FASTP_TAG, None)
                .to_string(),
        )
        // These default QC thresholds are loosely based on the
        // pre-2022 CZI pipeline using PriceSeq & LZW.
        .flag("--interleaved_in")
        .flag("--dont_eval_duplication")
        .option("--length_required", "35")
        .option("--qualified_quality_phred", "17")
        .option("--unqualified_percent_limit", "15")
        .option("--n_base_limit", "15")
        .option("--complexity_threshold", "60");
    // .flag("--sdust_complexity_filter");

    if let Some(adapter_fasta) = &config.args.adapter_fasta {
        let cwd = std::env::current_dir()?;

        let adapter_path = file_path_manipulator(
            &PathBuf::from(adapter_fasta),
            Some(&cwd),
            None,
            None,
            "",
        );

        if !adapter_path.exists() {
            return Err(PipelineError::InvalidConfig(format!(
                "Adapter FASTA file does not exist: {}",
                adapter_path.display()
            )));
        }

        qc_fastp_config_view = qc_fastp_config_view.option(
            "--adapter_fasta",
            adapter_path.to_string_lossy().into_owned(),
        );

        if paired {
            qc_fastp_config_view =
                qc_fastp_config_view.flag("--detect_adapter_for_pe");
        }
    }

    let qc_fastp_args =
        generate_cli(
            FASTP_TAG,
            &config,
            Some(&qc_fastp_config_view),
        )
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?;
    let (qc_fastp_child, qc_fastp_stream_task, qc_fastp_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        FASTP_TAG,
        qc_fastp_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(qc_fastp_stream_task);
    cleanup_tasks.push(qc_fastp_err_task);

    let qc_fastp_out_stream = {
        let mut guard = qc_fastp_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Tee the byte stream for counting
    let tee_count_input = ReceiverStream::new(qc_fastp_out_stream);

    let (tee_count_streams, tee_count_router_done_rx) = fanout_to_channels(
        tee_count_input,
        2,
        "fastp_count",
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "fanout_to_channels".to_string(),
            error: e.to_string(),
        })?;

    // track the task instead of oneshot receiver
    cleanup_receivers.push(tee_count_router_done_rx);

    let mut tee_count_streams_iter = tee_count_streams.into_iter();

    let qc_out_stream = ReceiverStream::new(
        tee_count_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let qc_count_stream = ReceiverStream::new(
        tee_count_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let (count_tx, count_result_rx) = oneshot::channel::<Result<u64, anyhow::Error>>();

    let count_task = tokio::spawn(async move {
        match stream_record_counter(qc_count_stream.into_inner(), false).await {
            Ok(stats) => {
                debug!("Fastp output reads: {}", stats.reads);
                let _ = count_tx.send(Ok(stats.reads));
                Ok(())
            }
            Err(e) => {
                let _ = count_tx.send(Err(e));
                Err(anyhow!("Failed to count fastp output reads"))
            }
        }
    });

    cleanup_tasks.push(count_task);

    Ok((
        qc_out_stream,
        cleanup_tasks,
        cleanup_receivers,
        count_result_rx,
    ))
}

/// Runs Kallisto quantification on the input FASTQ stream.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `out_dir`: directory where Kallisto outputs will be written
/// * `paired`: whether the input is paired-end
/// * `sample_base_buf`: base path for the sample
///
/// # Returns
///
/// Result containing:
/// - sender for ERCC results
/// - receiver for ERCC results
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - handle for the Kallisto task
async fn kallisto_quant(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    out_dir: PathBuf,
    paired: bool,
    sample_base_buf: PathBuf,
) -> Result<
    (
        oneshot::Sender<KallistoResults>,
        oneshot::Receiver<KallistoResults>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        JoinHandle<Result<(), anyhow::Error>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    // Convert sample_base_buf to &str
    let sample_base = sample_base_buf
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| {
            PipelineError::IOError(format!(
                "Invalid sample base path: {}",
                sample_base_buf.display()
            ))
        })?
        .to_string();

    // Create Kallisto output dir (idempotent, fast on NVMe)
    let kallisto_out_dir = out_dir.join("kallisto");
    fs::create_dir_all(&kallisto_out_dir).await.map_err(|e| {
        PipelineError::IOError(format!("Failed to create Kallisto output directory: {}", e))
    })?;
    debug!(
        "Kallisto output dir created/confirmed: {}",
        kallisto_out_dir.display()
    );

    // Convert raw byte stream to FASTQ
    let byte_rx = input_stream.into_inner();
    let byte_reader = ChannelReader::new(byte_rx);
    let fastq_rx = parse_fastq(byte_reader, &config, StreamDataType::IlluminaFastq)
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "parse_fastq".to_string(),
            error: e.to_string(),
        })?;
    let fastq_stream = ReceiverStream::new(fastq_rx);

    debug!("Deinterleaving FASTQ stream for sample: {}", sample_base);
    let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) =
        deinterleave_fastq_stream_to_fifos(config.clone(), fastq_stream, &sample_base, paired)
            .await?;

    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle_opt {
        cleanup_tasks.push(r2_handle);
    }
    let kallisto_index = config
        .args
        .kallisto_index
        .clone()
        .ok_or_else(|| {
            PipelineError::MissingArgument(
                "Must provide kallisto index with --kallisto-index arg.".to_string(),
            )
        })?;

    let kallisto_threads = config
        .thread_allocation(KALLISTO_TAG, None)
        .to_string();

    let kallisto_config = if paired {
        let r2_fifo = r2_fifo_opt.as_ref().ok_or_else(|| {
            PipelineError::InvalidConfig(
                "Paired mode requested but R2 FIFO not created".to_string(),
            )
        })?;

        KallistoConfig::new(KallistoSubcommand::Quant)
            .option("-i", kallisto_index.clone())
            .option(
                "-o",
                kallisto_out_dir.to_string_lossy().to_string(),
            )
            .flag("--plaintext")
            .option("--threads", kallisto_threads.clone())
            .positional(r1_fifo.to_string_lossy().to_string())
            .positional(r2_fifo.to_string_lossy().to_string())
    } else {
        KallistoConfig::new(KallistoSubcommand::Quant)
            .option("-i", kallisto_index)
            .option(
                "-o",
                kallisto_out_dir.to_string_lossy().to_string(),
            )
            .flag("--plaintext")
            .option("--threads", kallisto_threads)
            .flag("--single")
            .option("-l", "200")
            .option("-s", "20")
            .positional(r1_fifo.to_string_lossy().to_string())
    };

    let kallisto_args =
        generate_cli(KALLISTO_TAG, &config, Some(&kallisto_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: KALLISTO_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    debug!("Spawning Kallisto with args: {:?}", kallisto_args);

    let (mut kallisto_child, kallisto_err_task) = spawn_cmd(
        config.clone(),
        KALLISTO_TAG,
        kallisto_args,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(kallisto_err_task);

    // Create Kallisto exit task with safe FIFO cleanup
    let r1_fifo_clone = r1_fifo.clone();
    let r2_fifo_opt_clone = r2_fifo_opt.clone();
    let kallisto_exit_task = tokio::spawn(async move {
        let status = kallisto_child
            .wait()
            .await
            .map_err(|e| anyhow!("Kallisto wait failed: {}", e))?;
        if !status.success() {
            return Err(anyhow!("Kallisto exited with code: {:?}", status.code()));
        }

        // Clean up R1 FIFO
        tokio::fs::remove_file(&r1_fifo_clone).await.map_err(|e| {
            anyhow!(
                "Failed to remove R1 FIFO {}: {}",
                r1_fifo_clone.display(),
                e
            )
        })?;

        // Clean up R2 FIFO only if it exists
        if let Some(r2_fifo) = r2_fifo_opt_clone {
            tokio::fs::remove_file(&r2_fifo)
                .await
                .map_err(|e| anyhow!("Failed to remove R2 FIFO {}: {}", r2_fifo.display(), e))?;
        }

        Ok(())
    });

    // Create oneshot channel for results
    let (ercc_tx, ercc_rx) = oneshot::channel::<KallistoResults>();

    Ok((
        ercc_tx,
        ercc_rx,
        cleanup_tasks,
        cleanup_receivers,
        kallisto_exit_task,
    ))
}

/// Processes Kallisto output files to extract ERCC quantification results.
///
/// # Arguments
///
/// * `kallisto_out_dir`: directory containing Kallisto results
/// * `ercc_tx`: sender to transmit the extracted results
///
/// # Returns
///
/// Result containing a handle to the processing task.
async fn kallisto_results(
    kallisto_out_dir: PathBuf,
    ercc_tx: oneshot::Sender<KallistoResults>,
) -> Result<JoinHandle<Result<(), anyhow::Error>>, PipelineError> {
    let abundance_path = kallisto_out_dir.join("abundance.tsv");

    let parse_task = tokio::spawn(async move {
        // Wait for non-empty file
        let mut attempts = 5;
        let file = loop {
            match fs::metadata(&abundance_path).await {
                Ok(meta) if meta.len() > 0 => break TokioFile::open(&abundance_path).await?,
                _ if attempts == 0 => {
                    return Err(anyhow!("kallisto abundance.tsv never appeared or is empty"));
                }
                _ => {
                    attempts -= 1;
                    warn!(
                        "Waiting for kallisto abundance.tsv ({} attempts left)...",
                        attempts + 1
                    );
                    sleep(Duration::from_secs(2)).await;
                }
            }
        };

        let mut lines = TokioBufReader::new(file).lines();

        // Skip header
        if lines.next_line().await?.is_none() {
            warn!("kallisto abundance.tsv has no header — assuming empty");
            let _ = ercc_tx.send(KallistoResults {
                ercc_counts: vec![],
                transcript_to_gene: vec![],
            });
            return Ok(());
        }

        let mut ercc_counts = Vec::new();
        let mut transcript_to_gene = Vec::new();
        let mut parsed = 0;
        let mut skipped = 0;

        while let Some(line_result) = lines.next_line().await.transpose() {
            let line = match line_result {
                Ok(l) => l,
                Err(e) => {
                    warn!("I/O error reading kallisto abundance.tsv: {}", e);
                    skipped += 1;
                    continue;
                }
            };

            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                warn!(
                    "Malformed kallisto line ({} fields): {}",
                    fields.len(),
                    line
                );
                skipped += 1;
                continue;
            }

            let target_id = fields[0].to_string();

            let est_counts: f64 = match fields[3].parse() {
                Ok(v) => v,
                Err(e) => {
                    warn!(
                        "Invalid est_counts '{}' → {} (line: {})",
                        fields[3], e, line
                    );
                    skipped += 1;
                    continue;
                }
            };

            // Optional: parse TPM too
            let _tpm: f64 = fields[4].parse().unwrap_or(0.0);

            // Keep ERCCs only
            if target_id.starts_with("ERCC-") {
                ercc_counts.push((target_id.clone(), est_counts));
            }

            transcript_to_gene.push((target_id, format!("gene_{}", fields[0])));
            parsed += 1;
        }

        info!(
            "Kallisto parsing finished — {} lines parsed, {} skipped, {} ERCCs found",
            parsed,
            skipped,
            ercc_counts.len()
        );

        ercc_tx
            .send(KallistoResults {
                ercc_counts,
                transcript_to_gene,
            })
            .map_err(|_| anyhow!("Kallisto receiver dropped"))?;

        Ok(())
    });

    Ok(parse_task)
}

/// Filters reads using minimap2 and returns a stream of unmapped FASTQ records.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `paired`: whether the input is paired-end
/// * `output_bam_path`: optional path to save the aligned BAM file
///
/// # Returns
///
/// Result containing:
/// - stream of unmapped FASTQ records
/// - receiver for the count of mapped reads
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - vector of temporary files created
/// - vector of temporary directories used
#[allow(dead_code)]
async fn minimap2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    output_bam_path: Option<PathBuf>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        oneshot::Receiver<u64>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<NamedTempFile>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();
    let mut temp_dirs: Vec<TempDir> = Vec::new();
    let mut named_temp_files: Vec<NamedTempFile> = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    let (
        _host_ref_fasta_path,
        host_index_path,
        host_ref_temp,
        host_index_temp,
        host_temp_dir,
        host_ref_tasks,
    ) = minimap2_index_prep(
        &config,
        &config.ram_temp_dir,
        config.args.host_sequence.clone(),
        config.args.host_mm2_index.clone(),
        "host",
    )
        .await?;
    try_join_all(host_ref_tasks).await?; // Have to wait on this to make sure index is in place ahead of running minimap2
    debug!("Host index : {}", host_index_path.display());

    if let Some(temp) = host_ref_temp {
        named_temp_files.push(temp);
    }
    if let Some(temp) = host_index_temp {
        named_temp_files.push(temp);
    }
    if let Some(temp_dir) = host_temp_dir {
        temp_dirs.push(temp_dir);
    }

    let minimap2_config = Minimap2Config::new()
        .option(
            "-t",
            config
                .thread_allocation(MINIMAP2_TAG, None)
                .to_string(),
        )
        .option("-ax", "sr") // Short-read preset
        .option("-k", "15") // Default k-mer
        .option("-w", "10") // Default window
        .option("-s", "30") // Match HISAT2's score threshold
        .option("-A", "2") // Match score
        .option("-B", "4") // Mismatch penalty
        .option("-O", "4,24") // Gap open
        .option("-E", "2,1") // Gap extension
        .positional(host_index_path.to_string_lossy().to_string())
        .positional("-");

    let minimap2_args =
        generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            }
        })?;

    let (minimap2_child, minimap2_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::IlluminaFastq, // Assuming interleaved FASTQ
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("sort"))
                .to_string(),
        )
        .option("-m", "2G")
        .flag("-n")
        .flag("-u")
        .option("-O", "bam")
        .option("-T", temp_dir.path().to_string_lossy().to_string())
        .positional("-");

    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    use tokio_stream::wrappers::ReceiverStream;

    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_router_done_rx) = fanout_to_channels(
        bam_rx_stream,
        num_tees,
        "minimap2_bam_split",
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    // track task instead of oneshot receiver
    cleanup_receivers.push(bam_router_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream =
            ReceiverStream::new(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?);

        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            stream,
            config.clone(),
            StreamDataType::JustBytes,
            "mm2_filter_bam",
            false
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream =
        ReceiverStream::new(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?);
    let mapped_flag = if paired {
        "-F12".to_string()
    } else {
        "-F4".to_string()
    };
    let samtools_count_config = SamtoolsConfig::new(SamtoolsSubcommand::View)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("view"))
                .to_string(),
        )
        .flag("--no-PG")
        .flag("-c")
        .flag(mapped_flag)
        .positional("-");
    let samtools_count_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_count_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        bam_count_stream.into_inner(),
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (_count_tx, count_rx) = oneshot::channel::<u64>();
    let count_future = tokio::spawn({
        let config = config.clone();
        async move {
            let mut guard = count_child_arc.lock().await;
            let _count_lines =
                read_child_output_to_vec(&mut guard, ChildStream::Stdout, &config).await?;
            Ok(())
        }
    });
    cleanup_tasks.push(count_future);

    // Extract unmapped FASTQ
    let unmapped_stream =
        ReceiverStream::new(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?);

    let unmapped_flag = if paired {
        "-f13".to_string()
    } else {
        "-f4".to_string()
    };
    let samtools_fastq_config = SamtoolsConfig::new(SamtoolsSubcommand::Fastq)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("fastq"))
                .to_string(),
        )
        .option("-c", "6")
        .flag(unmapped_flag)
        .positional("-");
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream.into_inner(),
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Fastq, &config)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    temp_dirs.push(temp_dir);
    Ok((
        ReceiverStream::new(unmapped_fastq_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        named_temp_files,
        temp_dirs,
    ))
}

/// Deduplicates the input FASTQ stream.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_stream`: stream of interleaved FASTQ records
/// * `paired`: whether the input is paired-end
/// * `prefix_len`: length of the prefix to use for deduplication (if any)
/// * `out_dir`: directory where deduplication stats will be written
///
/// # Returns
///
/// Result containing:
/// - stream of unique FASTQ records
/// - receiver for the count of unique reads
/// - map of cluster information (representative to members)
/// - vector of cleanup tasks
/// - vector of cleanup receivers
async fn dedup(
    config: Arc<RunConfig>,
    input_stream: mpsc::Receiver<ParseOutput>,
    paired: bool,
    prefix_len: Option<usize>,
    out_dir: PathBuf,
) -> Result<
    (
        ReceiverStream<ParseOutput>, // uniques stream
        oneshot::Receiver<u64>,      // unique count rx

        // Completes after all dedup workers have finished and returns the
        // complete duplicate-cluster map.
        JoinHandle<Result<Arc<DashMap<String, ClusterInfo>>, anyhow::Error>>,

        Vec<JoinHandle<Result<(), anyhow::Error>>>,        // cleanup tasks
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, // cleanup rx
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    // CPU-bound hashing, moderate memory
    let num_workers = compute_phase_concurrency(
        &config,
        "dedup",
        1.0, // ~1GB per worker
        1.0, // 1 thread per core
        config.max_cores,
        4, // min 4 workers
    );

    let num_shards = (num_workers * 16).max(1024);

    let shards: Vec<
        Arc<Mutex<ahash::AHashMap<u64, (String, u64, Vec<String>)>>>,
    > = (0..num_shards)
        .map(|_| Arc::new(Mutex::new(ahash::AHashMap::new())))
        .collect();

    let total_count = Arc::new(AtomicU64::new(0));
    let unique_count_atomic = Arc::new(AtomicU64::new(0));

    // -----------------------------------------------------------------
    // clusters.csv
    //
    // This consumer is installed immediately. Workers may therefore
    // write cluster membership records without waiting for dedup() to
    // return.
    // -----------------------------------------------------------------

    let (csv_tx, csv_rx) =
        mpsc::channel::<(String, String)>(
            config.base_buffer_size * num_workers,
        );

    let csv_path = out_dir.join("clusters.csv");

    let csv_task = tokio::spawn(async move {
        let mut file = BufWriter::new(
            TokioFile::create(&csv_path)
                .await
                .context("Failed to create clusters.csv")?,
        );

        file.write_all(b"representative read id,read id\n")
            .await
            .context("CSV header write failed")?;

        let mut rx = csv_rx;

        while let Some((rep, member)) = rx.recv().await {
            file.write_all(
                format!("{},{}\n", rep, member).as_bytes(),
            )
                .await
                .context("CSV write failed")?;
        }

        file.flush()
            .await
            .context("CSV flush failed")?;

        Ok(())
    });

    cleanup_tasks.push(csv_task);

    // -----------------------------------------------------------------
    // Unique-read output.
    //
    // This is the bounded computational channel that must be returned
    // while its producers are still running.
    // -----------------------------------------------------------------

    let (uniques_tx, uniques_rx) =
        mpsc::channel::<ParseOutput>(
            config.base_buffer_size * num_workers,
        );

    let uniques_tx =
        Arc::new(tokio::sync::Mutex::new(uniques_tx));

    // -----------------------------------------------------------------
    // Dedup workers.
    // -----------------------------------------------------------------

    let mut worker_handles:
        Vec<JoinHandle<Result<(), anyhow::Error>>> =
        Vec::with_capacity(num_workers);

    let worker_txs: Vec<mpsc::Sender<ParseOutput>> =
        (0..num_workers)
            .map(|_| {
                let (tx, rx) =
                    mpsc::channel(config.base_buffer_size);

                let csv_tx_clone = csv_tx.clone();
                let uniques_tx_clone = uniques_tx.clone();

                let handle = tokio::spawn(dedup_worker(
                    rx,
                    shards.clone(),
                    uniques_tx_clone,
                    csv_tx_clone,
                    paired,
                    prefix_len,
                    total_count.clone(),
                    unique_count_atomic.clone(),
                ));

                worker_handles.push(handle);

                tx
            })
            .collect();

    // The original uniques sender is no longer needed. Each worker owns
    // its own Arc clone. Once all workers finish, the final sender is
    // dropped and uniques_rx receives clean EOF.
    drop(uniques_tx);

    let (count_tx, count_rx) =
        oneshot::channel::<u64>();

    // -----------------------------------------------------------------
    // Coordinator / finalizer.
    //
    // CRITICAL LIFECYCLE RULE:
    //
    // Do not perform distribution or worker joins synchronously inside
    // dedup(). Both can ultimately block behind the bounded uniques
    // channel.
    //
    // Returning uniques_rx immediately allows the caller to install
    // subsample_uniform_reservoir() as its concurrent consumer.
    // -----------------------------------------------------------------

    let dedup_finalize_task = tokio::spawn(async move {
        let mut stream =
            ReceiverStream::new(input_stream);

        let mut i = 0usize;

        // -------------------------------------------------------------
        // Distributor: round-robin complete pairs to workers.
        //
        // The input remains interleaved:
        //
        //   R1, R2, R1, R2, ...
        //
        // Validate each pair before sending either record to a worker.
        // -------------------------------------------------------------

        if paired {
            while let Some(r1_item) = stream.next().await {
                let r1 = match r1_item {
                    ParseOutput::Fastq(rec) => rec,

                    _ => {
                        return Err(anyhow!(
                            "Unexpected non-FASTQ item in paired dedup input"
                        ));
                    }
                };

                let r2_item =
                    stream.next().await.ok_or_else(|| {
                        anyhow!(
                            "Missing R2 in paired dedup stream"
                        )
                    })?;

                let r2 = match r2_item {
                    ParseOutput::Fastq(rec) => rec,

                    _ => {
                        return Err(anyhow!(
                            "Unexpected non-FASTQ R2 item in paired dedup input"
                        ));
                    }
                };

                if r1.id() != r2.id() {
                    return Err(anyhow!(
                        "Mismatched pair IDs: R1={}, R2={}",
                        r1.id(),
                        r2.id()
                    ));
                }

                let worker_tx =
                    &worker_txs[i % num_workers];

                worker_tx
                    .send(ParseOutput::Fastq(r1))
                    .await
                    .map_err(|e| {
                        anyhow!(
                            "Send R1 to dedup worker failed: {}",
                            e
                        )
                    })?;

                worker_tx
                    .send(ParseOutput::Fastq(r2))
                    .await
                    .map_err(|e| {
                        anyhow!(
                            "Send R2 to dedup worker failed: {}",
                            e
                        )
                    })?;

                i += 1;
            }
        } else {
            while let Some(item) = stream.next().await {
                let worker_tx =
                    &worker_txs[i % num_workers];

                worker_tx
                    .send(item)
                    .await
                    .map_err(|e| {
                        anyhow!(
                            "Send to dedup worker failed: {}",
                            e
                        )
                    })?;

                i += 1;
            }
        }

        // Input is exhausted. Close every worker input channel.
        drop(worker_txs);

        // -------------------------------------------------------------
        // Workers can now finish.
        //
        // This wait is safe because uniques_rx has already been returned
        // from dedup() and the caller has installed its consumer.
        // -------------------------------------------------------------

        for handle in worker_handles {
            handle
                .await
                .map_err(|e| {
                    anyhow!(
                        "Dedup worker join failed: {}",
                        e
                    )
                })?
                .context("Dedup worker error")?;
        }

        // All worker-owned csv_tx clones are gone now. Drop the
        // coordinator's original sender so clusters.csv receives EOF.
        drop(csv_tx);

        // -------------------------------------------------------------
        // Build the complete duplicate-cluster map.
        //
        // Workers are finished, so the shards are now immutable.
        // -------------------------------------------------------------

        let duplicate_clusters =
            Arc::new(
                DashMap::with_capacity(num_shards),
            );

        for shard in &shards {
            let guard = shard
                .lock()
                .map_err(|e| {
                    anyhow!(
                        "Shard lock poisoned: {}",
                        e
                    )
                })?;

            for (_, (rep_id, size, members)) in guard.iter() {
                duplicate_clusters.insert(
                    rep_id.clone(),
                    ClusterInfo {
                        size: *size,
                        members: members.clone(),
                    },
                );
            }
        }

        // -------------------------------------------------------------
        // duplicate_cluster_sizes.tsv
        //
        // Preserve the existing artifact, but write it here after the
        // cluster map is complete.
        // -------------------------------------------------------------

        let tsv_path =
            out_dir.join("duplicate_cluster_sizes.tsv");

        let mut file = BufWriter::new(
            TokioFile::create(&tsv_path)
                .await
                .context(
                    "Failed to create duplicate_cluster_sizes.tsv",
                )?,
        );

        file.write_all(
            b"representative read id\tcluster size\n",
        )
            .await
            .context("TSV header write failed")?;

        for entry in duplicate_clusters.iter() {
            let (rep_id, info) = entry.pair();

            file.write_all(
                format!(
                    "{}\t{}\n",
                    rep_id,
                    info.size
                )
                    .as_bytes(),
            )
                .await
                .context("TSV write failed")?;
        }

        file.flush()
            .await
            .context("TSV flush failed")?;

        // -------------------------------------------------------------
        // Publish the final unique count.
        // -------------------------------------------------------------

        let unique_count =
            unique_count_atomic.load(
                AtomicOrdering::Relaxed,
            );

        count_tx
            .send(unique_count)
            .map_err(|_| {
                anyhow!(
                    "Unique count receiver dropped"
                )
            })?;

        info!(
            "Dedup complete: {} total records, {} unique clusters",
            total_count.load(AtomicOrdering::Relaxed),
            unique_count
        );

        Ok::<
            Arc<DashMap<String, ClusterInfo>>,
            anyhow::Error,
        >(duplicate_clusters)
    });

    // -----------------------------------------------------------------
    // Return immediately.
    //
    // The caller can now consume uniques_rx concurrently with the
    // distributor/workers above.
    // -----------------------------------------------------------------

    Ok((
        ReceiverStream::new(uniques_rx),
        count_rx,
        dedup_finalize_task,
        cleanup_tasks,
        cleanup_receivers,
    ))
}


/// For paired-end input, this worker always treats two consecutive records
/// as one logical pair. Unique pairs are emitted while holding the shared
/// output-sender mutex across BOTH sends, guaranteeing:
///
/// ```text
/// R1(pair A), R2(pair A), R1(pair B), R2(pair B), ...
/// ```
///
/// even when multiple dedup workers are running concurrently.
async fn dedup_worker(
    mut rx: mpsc::Receiver<ParseOutput>,
    shards: Vec<Arc<Mutex<ahash::AHashMap<u64, (String, u64, Vec<String>)>>>>,
    uniques_tx: Arc<tokio::sync::Mutex<mpsc::Sender<ParseOutput>>>,
    csv_tx: mpsc::Sender<(String, String)>,
    paired: bool,
    prefix_len: Option<usize>,
    total_count: Arc<AtomicU64>,
    unique_count: Arc<AtomicU64>,
) -> Result<(), anyhow::Error> {
    let mut orphan_r1: Option<SequenceRecord> = None;

    while let Some(item) = rx.recv().await {
        let record = match item {
            ParseOutput::Fastq(rec) => rec,
            _ => {
                return Err(anyhow!(
                    "dedup_worker received unexpected non-FASTQ item"
                ));
            }
        };

        if paired {
            if orphan_r1.is_none() {
                orphan_r1 = Some(record);
                continue;
            }

            let r1 = orphan_r1.take().unwrap();
            let r2 = record;

            if r1.id() != r2.id() {
                return Err(anyhow!(
                    "Mismatched pair IDs in worker: {} != {}",
                    r1.id(),
                    r2.id()
                ));
            }

            let id = r1.id().to_string();

            let r1_prefix_len = prefix_len
                .map_or(r1.seq().len(), |l| l.min(r1.seq().len()));
            let r2_prefix_len = prefix_len
                .map_or(r2.seq().len(), |l| l.min(r2.seq().len()));

            let r1_prefix = &r1.seq()[..r1_prefix_len];
            let r2_prefix = &r2.seq()[..r2_prefix_len];

            let mut hasher = XxHash64::default();
            hasher.write(r1_prefix);
            hasher.write_u8(0);
            hasher.write(r2_prefix);

            let hash_key = hasher.finish();
            let shard_idx = (hash_key as usize) % shards.len();

            // Only hold the shard lock for the map operation.
            let (rep_id, is_new) = {
                let mut shard = shards[shard_idx]
                    .lock()
                    .map_err(|e| anyhow!("Shard lock poisoned: {}", e))?;

                let entry = shard.entry(hash_key).or_insert_with(|| {
                    unique_count.fetch_add(1, AtomicOrdering::Relaxed);
                    (id.clone(), 0, vec![])
                });

                entry.1 += 1;
                entry.2.push(id.clone());

                let rep_id = entry.0.clone();
                let is_new = entry.1 == 1;

                (rep_id, is_new)
            };

            total_count.fetch_add(1, AtomicOrdering::Relaxed);

            csv_tx
                .send((rep_id, id))
                .await
                .context("CSV send failed")?;

            if is_new {
                // CRITICAL:
                // Hold the sender lock across BOTH sends. This makes the pair
                // atomic at the logical output-stream level.
                let tx = uniques_tx.lock().await;

                tx.send(ParseOutput::Fastq(r1))
                    .await
                    .context("Send R1 unique failed")?;

                tx.send(ParseOutput::Fastq(r2))
                    .await
                    .context("Send R2 unique failed")?;

                // tx is dropped here, allowing another worker to emit.
            }
        } else {
            let prefix_len =
                prefix_len.map_or(record.seq().len(), |l| l.min(record.seq().len()));

            let prefix = &record.seq()[..prefix_len];

            let mut hasher = XxHash64::default();
            hasher.write(prefix);

            let hash_key = hasher.finish();
            let shard_idx = (hash_key as usize) % shards.len();

            let (rep_id, is_new) = {
                let mut shard = shards[shard_idx]
                    .lock()
                    .map_err(|e| anyhow!("Shard lock poisoned: {}", e))?;

                let id = record.id().to_string();

                let entry = shard.entry(hash_key).or_insert_with(|| {
                    unique_count.fetch_add(1, AtomicOrdering::Relaxed);
                    (id.clone(), 0, vec![])
                });

                entry.1 += 1;
                entry.2.push(id.clone());

                let rep_id = entry.0.clone();
                let is_new = entry.1 == 1;

                (rep_id, is_new)
            };

            total_count.fetch_add(1, AtomicOrdering::Relaxed);

            csv_tx
                .send((rep_id, record.id().to_string()))
                .await
                .context("CSV send failed")?;

            if is_new {
                let tx = uniques_tx.lock().await;

                tx.send(ParseOutput::Fastq(record))
                    .await
                    .context("Send unique failed")?;
            }
        }
    }

    if orphan_r1.is_some() {
        return Err(anyhow!("Orphan R1 in worker at end"));
    }

    Ok(())
}


/// Alignment of deduped/subsampled stream to a non-host DB (canoncially NT)
/// for output in PAF format, allowing downstream metagenomic calling.
///
/// - One task per chunk (concurrency already limited by estimation)
/// - Direct byte streaming from stdout
///  - Workers fire-and-forget into cleanup_tasks (no blocking here)
/// - Merged stream returns immediately so paf_to_m8 can start consuming right away
///
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream (paired or single-end).
///
/// # Returns
/// Tuple:
/// - Interleaved FASTQ stream : unmapped reads.
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
/// - Optional temporary FASTA file for non-host reference.
/// - Optional temporary minimap2 index file.
/// Aligns unmapped reads against NT index chunks using minimap2.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `r1_path`: path to the R1 unmapped FASTQ file
/// * `r2_path_opt`: optional path to the R2 unmapped FASTQ file
///
/// # Returns
///
/// Result containing:
/// - stream of alignment results (PAF format)
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// Short-read NT non-host alignment.
///
/// Short-read-specific wrapper around the shared sharded minimap2 aligner.
async fn nt_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let alignment = ShardedMinimap2Config {
        index_dir: PathBuf::from(&config.args.nt_split_dir),
        preset: Minimap2Preset::ShortRead,
        query1: r1_path,
        query2: r2_path_opt,
    };

    sharded_minimap2_align(config, alignment).await
}

/// Converts a PAF stream to an BLAST m8 format stream and writes it to a file.
///
/// # Arguments
///



async fn nr_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    sample_id: String,
) -> Result<
    (
        mpsc::Receiver<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    match config.execution_mode {
        ExecutionMode::Single => match config.alignment_backend {
            NRAlignmentBackend::Diamond => {
                diamond_non_host_align(
                    config,
                    r1_path,
                    r2_path_opt,
                )
                    .await
            }

            NRAlignmentBackend::MmseqsCpu => {
                mmseqs_non_host_align(
                    config,
                    r1_path,
                    r2_path_opt,
                    MmseqsBackend::Cpu,
                )
                    .await
            }

            NRAlignmentBackend::MmseqsGpu => {
                mmseqs_non_host_align(
                    config,
                    r1_path,
                    r2_path_opt,
                    MmseqsBackend::Gpu,
                )
                    .await
            }
        },

        ExecutionMode::Distributed => {
            nr_distributed_non_host_align(
                config,
                r1_path,
                r2_path_opt,
                sample_id,
            )
                .await
        }
    }
}



async fn get_reference_version(
    backend: NRAlignmentBackend,
) -> Result<String, PipelineError> {
    let key = match backend {
        NRAlignmentBackend::MmseqsCpu => {
            format!("{}/mmseqs-cpu.version", REFERENCE_METADATA_PREFIX)
        }

        NRAlignmentBackend::MmseqsGpu => {
            format!("{}/mmseqs-gpu.version", REFERENCE_METADATA_PREFIX)
        }

        NRAlignmentBackend::Diamond => {
            format!("{}/diamond.version", REFERENCE_METADATA_PREFIX)
        }
    };

    let sdk_config =
        aws_config::load_defaults(
            aws_config::BehaviorVersion::latest(),
        )
            .await;

    let s3_client =
        aws_sdk_s3::Client::new(&sdk_config);

    let response = s3_client
        .get_object()
        .bucket(REFERENCE_METADATA_BUCKET)
        .key(&key)
        .send()
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to retrieve reference version from s3://{}/{}: {}",
                REFERENCE_METADATA_BUCKET,
                key,
                e
            ))
        })?;

    let bytes = response
        .body
        .collect()
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed to read reference version object s3://{}/{}: {}",
                REFERENCE_METADATA_BUCKET,
                key,
                e
            ))
        })?
        .into_bytes();

    let version = String::from_utf8(bytes.to_vec())
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Reference version object s3://{}/{} is not valid UTF-8: {}",
                REFERENCE_METADATA_BUCKET,
                key,
                e
            ))
        })?
        .trim()
        .to_string();

    if version.is_empty() {
        return Err(PipelineError::Other(anyhow!(
            "Reference version object s3://{}/{} is empty",
            REFERENCE_METADATA_BUCKET,
            key
        )));
    }

    info!(
        "Distributed NR reference version: backend={:?}, version={}",
        backend,
        version
    );

    Ok(version)
}

/// Runs the distributed non-host alignment preparation path.
///
/// Copies the non-host FASTQs to EFS, creates paired FASTQ chunks, creates
/// one AVAILABLE WorkUnit per chunk, and leaves those work units on EFS for
/// workers to claim.
///
/// Worker discovery is best-effort and does not prevent work-unit creation.
async fn nr_distributed_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    sample_id: String,
) -> Result<
    (
        mpsc::Receiver<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let requested_workers =
        config.distributed_workers;

    let efs_base =
        config.efs_runs_dir.join(&config.run_id);

    let non_host_r1_efs =
        efs_base.join("nonhost_R1.fastq");

    let non_host_r2_efs =
        r2_path_opt
            .as_ref()
            .map(|_| {
                efs_base.join("nonhost_R2.fastq")
            });

    // ------------------------------------------------------------------
    // 1. Copy non-host FASTQs to EFS.
    // ------------------------------------------------------------------

    info!(
        "Copying non-host R1 to EFS: {}",
        non_host_r1_efs.display()
    );

    copy_file_streaming(
        &r1_path,
        &non_host_r1_efs,
        Some(config.clone()),
        None,
    )
        .await
        .map_err(PipelineError::Other)?
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "R1 copy task join failed: {e}"
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "R1 EFS copy failed: {e}"
            ))
        })?;

    if let (
        Some(local_r2),
        Some(efs_r2),
    ) = (
        &r2_path_opt,
        &non_host_r2_efs,
    ) {
        info!(
            "Copying non-host R2 to EFS: {}",
            efs_r2.display()
        );

        copy_file_streaming(
            local_r2,
            efs_r2,
            Some(config.clone()),
            None,
        )
            .await
            .map_err(PipelineError::Other)?
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "R2 copy task join failed: {e}"
                ))
            })?
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "R2 EFS copy failed: {e}"
                ))
            })?;
    }

    info!(
        "Non-host FASTQs on EFS under {}",
        efs_base.display()
    );

    // ------------------------------------------------------------------
    // 2. Chunk FASTQs.
    //
    // Chunk count is determined by non-host FASTQ volume rather than
    // worker count. Worker count determines execution capacity; chunk
    // count determines the amount of Diamond work in each WorkUnit.
    //
    // Empirical target:
    //   ~75 MiB of R1 + ~75 MiB of R2 per paired-end WorkUnit.
    //
    // This is deliberately aggressive because Diamond has a substantial
    // per-invocation database-processing cost. We want enough query
    // material in each invocation to amortize that cost while keeping
    // individual Diamond runs in the approximate 1-1.5 hour target range.
    // ------------------------------------------------------------------

    const TARGET_CHUNK_BYTES_PER_END: u64 = 75 * 1024 * 1024;

    let total_records =
        raw_read_count(
            non_host_r1_efs.clone(),
            non_host_r2_efs.clone(),
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
            "Non-host FASTQ count task join failed: {e}"
        ))
            })?
            .map_err(|e| {
                PipelineError::Other(anyhow!(
            "Failed to count non-host FASTQ records: {e}"
        ))
            })?;

    if total_records == 0 {
        return Err(PipelineError::EmptyStream);
    }

    let r1_bytes = tokio::fs::metadata(&non_host_r1_efs)
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
            "Failed to stat non-host R1 FASTQ {}: {e}",
            non_host_r1_efs.display()
        ))
        })?
        .len();

    let r2_bytes = match &non_host_r2_efs {
        Some(path) => tokio::fs::metadata(path)
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Failed to stat non-host R2 FASTQ {}: {e}",
                path.display()
            ))
            })?
            .len(),
        None => 0,
    };

    // Use the larger end for conservative paired-end chunk sizing.
    let bytes_per_end = r1_bytes.max(r2_bytes);

    let target_chunks = bytes_per_end
        .div_ceil(TARGET_CHUNK_BYTES_PER_END)
        .max(1);

    info!(
    "Distributed NR size-based chunking: R1={} bytes, R2={} bytes, \
     {} requested workers, target {} chunks, target ~{} MiB/end/chunk",
    r1_bytes,
    r2_bytes,
    requested_workers,
    target_chunks,
    TARGET_CHUNK_BYTES_PER_END / (1024 * 1024)
);

    let chunks_dir = efs_base.join("chunks");

    let chunk_work_units:
        Vec<(
            u64,
            u64,
            PathBuf,
            Option<PathBuf>,
        )>;

    if let Some(non_host_r2_efs) =
        non_host_r2_efs.clone()
    {
        // --------------------------------------------------------------
        // PAIRED-END
        // --------------------------------------------------------------

        if total_records % 2 != 0 {
            return Err(
                PipelineError::Other(anyhow!(
                    "Paired non-host FASTQ record count is not even: {}",
                    total_records
                )),
            );
        }

        let total_pairs =
            total_records / 2;

        let target_chunks =
            target_chunks.min(total_pairs.max(1));

        let pairs_per_chunk =
            (total_pairs + target_chunks - 1)
                / target_chunks;

        info!(
        "Distributed NR chunking: {} pairs, {} requested workers, \
         size-based target {} chunks, {} pairs/chunk",
        total_pairs,
        requested_workers,
        target_chunks,
        pairs_per_chunk
    );

        let chunk_summary =
            chunk_paired_fastq(
                non_host_r1_efs.clone(),
                non_host_r2_efs,
                chunks_dir.clone(),
                pairs_per_chunk,
            )
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to chunk paired non-host FASTQs: {e}"
                    ))
                })?;

        info!(
            "Distributed NR chunking complete: {} pairs, {} R1 records, \
             {} R2 records, {} chunks",
            chunk_summary.total_pairs,
            chunk_summary.total_r1_records,
            chunk_summary.total_r2_records,
            chunk_summary.chunks.len()
        );

        if chunk_summary.total_pairs
            != total_pairs
        {
            return Err(
                PipelineError::Other(anyhow!(
                    "Distributed chunk reconciliation failed: \
                     source pairs={}, chunked pairs={}",
                    total_pairs,
                    chunk_summary.total_pairs
                )),
            );
        }

        chunk_work_units =
            chunk_summary
                .chunks
                .iter()
                .map(|chunk| {
                    (
                        chunk.chunk_id,
                        chunk.pair_count,
                        chunk.r1_path.clone(),
                        Some(chunk.r2_path.clone()),
                    )
                })
                .collect();
    } else {
        // --------------------------------------------------------------
        // SINGLE-END
        // --------------------------------------------------------------

        let total_reads =
            total_records;

        let target_chunks =
            (target_chunks as u64)
                .min(total_reads)
                .max(1);

        let reads_per_chunk =
            (total_reads
                + target_chunks
                - 1)
                / target_chunks;

        info!(
        "Distributed NR single-end chunking: {} reads, {} requested workers, \
         target {} chunks, {} reads/chunk",
        total_reads,
        requested_workers,
        target_chunks,
        reads_per_chunk
    );

        tokio::fs::create_dir_all(
            &chunks_dir,
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Failed to create distributed chunks directory {}: {}",
                chunks_dir.display(),
                e
            ))
            })?;

        let (rx, read_task) =
            read_fastq(
                non_host_r1_efs.clone(),
                None,
                PairingMode::Strict,
                None,
                u64::MAX,
                None,
                None,
                "distributed_non_host_align_single_end_chunking",
                &config,
            )
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                "Failed to start single-end FASTQ reader: {e}"
            ))
                })?;

        let byte_reader =
            ChannelReader::new(
                rx
            );

        let fastq_rx =
            parse_fastq(
                byte_reader,
                &config,
                StreamDataType::IlluminaFastq,
            )
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                    "Failed to start single-end FASTQ parser: {e}"
                ))
                })?;

        let mut stream =
            ReceiverStream::new(
                fastq_rx
            );

        let mut chunks =
            Vec::<(
                u64,
                u64,
                PathBuf,
            )>::new();

        let mut chunk_id =
            0u64;

        let mut chunk_reads =
            0u64;

        let mut total_reads_written =
            0u64;

        let mut current_chunk_path =
            chunks_dir.join(format!(
                "chunk_{:08}_R1.fastq",
                chunk_id
            ));

        let mut writer =
            Some(
                BufWriter::new(
                    tokio::fs::File::create(
                        &current_chunk_path
                    )
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                        "Failed to create single-end chunk {}: {}",
                        current_chunk_path.display(),
                        e
                    ))
                        })?,
                )
            );

        while let Some(item) =
            stream.next().await
        {
            let record =
                match item {
                    ParseOutput::Fastq(record) => record,

                    ParseOutput::Fasta(_) => {
                        return Err(
                            PipelineError::InvalidFastqFormat(
                                "Single-end distributed chunking received FASTA output"
                                    .to_string(),
                            )
                        );
                    }

                    ParseOutput::Bytes(_) => {
                        return Err(
                            PipelineError::InvalidFastqFormat(
                                "Single-end distributed chunking received byte output"
                                    .to_string(),
                            )
                        );
                    }
                };

            let bytes =
                record
                    .to_bytes()
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                        "Failed to serialize single-end FASTQ record: {e}"
                    ))
                    })?;

            writer
                .as_mut()
                .expect("single-end chunk writer missing")
                .write_all(&bytes)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                    "Failed writing single-end chunk {}: {}",
                    current_chunk_path.display(),
                    e
                ))
                })?;

            chunk_reads += 1;
            total_reads_written += 1;

            if chunk_reads == reads_per_chunk {
                let mut completed_writer =
                    writer
                        .take()
                        .expect("single-end chunk writer missing");

                completed_writer
                    .flush()
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                        "Failed flushing single-end chunk {}: {}",
                        current_chunk_path.display(),
                        e
                    ))
                    })?;

                drop(completed_writer);

                chunks.push((
                    chunk_id,
                    chunk_reads,
                    current_chunk_path.clone(),
                ));

                chunk_id += 1;
                chunk_reads = 0;

                if total_reads_written < total_reads {
                    current_chunk_path =
                        chunks_dir.join(format!(
                            "chunk_{:08}_R1.fastq",
                            chunk_id
                        ));

                    writer =
                        Some(
                            BufWriter::new(
                                tokio::fs::File::create(
                                    &current_chunk_path
                                )
                                    .await
                                    .map_err(|e| {
                                        PipelineError::Other(anyhow!(
                                    "Failed to create single-end chunk {}: {}",
                                    current_chunk_path.display(),
                                    e
                                ))
                                    })?,
                            )
                        );
                }
            }
        }

        if chunk_reads > 0 {
            let mut completed_writer =
                writer
                    .take()
                    .expect("single-end chunk writer missing");

            completed_writer
                .flush()
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                    "Failed flushing final single-end chunk {}: {}",
                    current_chunk_path.display(),
                    e
                ))
                })?;

            drop(completed_writer);

            chunks.push((
                chunk_id,
                chunk_reads,
                current_chunk_path.clone(),
            ));
        } else {
            // The last exact-size chunk already closed its writer.
            drop(writer.take());
        }

        read_task
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Single-end FASTQ reader task join failed: {e}"
            ))
            })?
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Single-end FASTQ reader failed during chunking: {e}"
            ))
            })?;

        if total_reads_written
            != total_reads
        {
            return Err(
                PipelineError::Other(anyhow!(
                "Distributed single-end chunk reconciliation failed: \
                 source reads={}, chunked reads={}",
                total_reads,
                total_reads_written
            )),
            );
        }

        chunk_work_units =
            chunks
                .into_iter()
                .map(|(
                          chunk_id,
                          read_count,
                          r1_path,
                      )| {
                    (
                        chunk_id,
                        read_count,
                        r1_path,
                        None,
                    )
                })
                .collect();

        info!(
            "Distributed NR single-end chunking complete: {} reads, {} chunks",
            total_reads_written,
            chunk_work_units.len()
        );
    }

    // ------------------------------------------------------------------
    // 3. Retrieve the required backend-specific reference version
    //    from the canonical S3 metadata object.
    // ------------------------------------------------------------------

    let reference_version =
        get_reference_version(
            config.alignment_backend,
        )
            .await?;

    // ------------------------------------------------------------------
    // 4. Create one AVAILABLE WorkUnit per chunk.
    // ------------------------------------------------------------------

    let work_dir =
        config.efs_runs_dir
            .join(&config.run_id)
            .join("work");

    for (
        chunk_id,
        expected_input_count,
        r1_chunk_path,
        r2_chunk_path,
    ) in &chunk_work_units
    {
        let work_unit =
            crate::utils::work_units::WorkUnit::new(
                config.run_id.clone(),
                sample_id.clone(),
                *chunk_id,
                r2_chunk_path.is_some(),
                *expected_input_count,
                reference_version.clone(),
                r1_chunk_path.clone(),
                r2_chunk_path.clone(),
            );

        let payload =
            serde_json::to_vec_pretty(
                &work_unit,
            )
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to serialize work unit {}: {}",
                        work_unit.id(),
                        e
                    ))
                })?;

        let final_path =
            work_dir.join(
                format!(
                    "work_{:08}.json",
                    chunk_id
                ),
            );

        let tmp_path =
            work_dir.join(
                format!(
                    "work_{:08}.json.tmp",
                    chunk_id
                ),
            );

        // Write the complete WorkUnit before making its
        // final filename visible to workers.
        tokio::fs::write(
            &tmp_path,
            payload,
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "Failed to write work unit {}: {}",
                    tmp_path.display(),
                    e
                ))
            })?;

        tokio::fs::rename(
            &tmp_path,
            &final_path,
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "Failed to publish work unit {} -> {}: {}",
                    tmp_path.display(),
                    final_path.display(),
                    e
                ))
            })?;

        info!(
            "Published distributed work unit {}: {} {} -> {}",
            work_unit.id(),
            expected_input_count,
            if r2_chunk_path.is_some() {
                "pairs"
            } else {
                "reads"
            },
            final_path.display()
        );
    }

    info!(
        "Distributed NR: {} AVAILABLE work units created under {}",
        chunk_work_units.len(),
        work_dir.display()
    );

    // ------------------------------------------------------------------
    // 5. Wait for the requested number of READY workers.
    //
    // Workers are started externally for now. Each worker is launched with
    // this run's work_dir and polls the shared EFS work directory.
    // ------------------------------------------------------------------

    let worker_manager =
        crate::utils::workers::WorkerManager::new(
            config.efs_base_dir.clone().into(),
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Distributed NR: could not initialize worker manager: {}",
                e
            ))
            })?;

    match worker_manager
        .discover_ready_workers(
            config.alignment_backend,
        )
        .await
    {
        Ok(workers) => {
            info!(
            "Distributed NR: discovered {} READY workers; \
             {} requested",
            workers.len(),
            requested_workers
        );

            for worker in &workers {
                info!(
                "Distributed NR worker: instance_id={}, \
                 private_ip={}, instance_type={}, az={:?}",
                worker.instance_id,
                worker.private_ip,
                worker.instance_type,
                worker.availability_zone,
            );
            }
        }

        Err(e) => {
            warn!(
            "Distributed NR: READY worker discovery failed: {}; \
             continuing to scheduler",
            e
        );
        }
    }

    // ------------------------------------------------------------------
    // 6. Monitor WorkUnits, detect stale workers, and retry failed work.
    //
    // The scheduler loop remains the component that drives progress.
    //
    // CLAIMED:
    //   A worker may die after claiming but before RUNNING. Because no
    //   heartbeat exists yet, use the WorkUnit JSON mtime as the claim age.
    //
    // RUNNING:
    //   Use last_heartbeat to determine whether the worker is still alive.
    //
    // FAILED:
    //   Retry only when failure.retryable=true.
    //
    // DONE:
    //   Never modify.
    //
    // Recovery uses the same per-work-unit lock used by workers.
    // ------------------------------------------------------------------

    loop {
        let mut available = 0usize;
        let mut claimed = 0usize;
        let mut running = 0usize;
        let mut done = 0usize;
        let mut failed = 0usize;
        let mut recovered = 0usize;

        let mut entries =
            tokio::fs::read_dir(&work_dir)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to read distributed work directory {}: {}",
                        work_dir.display(),
                        e
                    ))
                })?;

        while let Some(entry) =
            entries
                .next_entry()
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed reading distributed work directory {}: {}",
                        work_dir.display(),
                        e
                    ))
                })?
        {
            let path =
                entry.path();

            let Some(name) =
                path.file_name()
                    .and_then(|v| v.to_str())
            else {
                continue;
            };

            if !name.starts_with("work_")
                || !name.ends_with(".json")
            {
                continue;
            }

            let bytes =
                tokio::fs::read(&path)
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "Failed to read distributed work unit {}: {}",
                            path.display(),
                            e
                        ))
                    })?;

            let work_unit:
                crate::utils::work_units::WorkUnit =
                serde_json::from_slice(&bytes)
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "Failed to parse distributed work unit {}: {}",
                            path.display(),
                            e
                        ))
                    })?;

            match work_unit.state {
                crate::utils::work_units::WorkUnitState::Available => {
                    available += 1;
                }

                crate::utils::work_units::WorkUnitState::Claimed
                | crate::utils::work_units::WorkUnitState::Running => {
                    let now =
                        SystemTime::now()
                            .duration_since(UNIX_EPOCH)
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "System clock error while checking work-unit {}: {}",
                                    path.display(),
                                    e
                                ))
                            })?
                            .as_secs();

                    let stale =
                        match work_unit.state {
                            crate::utils::work_units::WorkUnitState::Claimed => {
                                let metadata =
                                    tokio::fs::metadata(&path)
                                        .await
                                        .map_err(|e| {
                                            PipelineError::Other(anyhow!(
                                                "Failed to stat claimed work unit {}: {}",
                                                path.display(),
                                                e
                                            ))
                                        })?;

                                let modified =
                                    metadata
                                        .modified()
                                        .map_err(|e| {
                                            PipelineError::Other(anyhow!(
                                                "Failed to read mtime for claimed work unit {}: {}",
                                                path.display(),
                                                e
                                            ))
                                        })?;

                                let modified_secs =
                                    modified
                                        .duration_since(UNIX_EPOCH)
                                        .map_err(|e| {
                                            PipelineError::Other(anyhow!(
                                                "Invalid mtime for claimed work unit {}: {}",
                                                path.display(),
                                                e
                                            ))
                                        })?
                                        .as_secs();

                                now.saturating_sub(modified_secs)
                                    > DISTRIBUTED_WORKER_STALE_SECS
                            }

                            crate::utils::work_units::WorkUnitState::Running => {
                                match work_unit.last_heartbeat {
                                    Some(last_heartbeat) => {
                                        now.saturating_sub(last_heartbeat)
                                            > DISTRIBUTED_WORKER_STALE_SECS
                                    }

                                    None => {
                                        let metadata =
                                            tokio::fs::metadata(&path)
                                                .await
                                                .map_err(|e| {
                                                    PipelineError::Other(anyhow!(
                                                        "Failed to stat running work unit {}: {}",
                                                        path.display(),
                                                        e
                                                    ))
                                                })?;

                                        let modified =
                                            metadata
                                                .modified()
                                                .map_err(|e| {
                                                    PipelineError::Other(anyhow!(
                                                        "Failed to read mtime for running work unit {}: {}",
                                                        path.display(),
                                                        e
                                                    ))
                                                })?;

                                        let modified_secs =
                                            modified
                                                .duration_since(UNIX_EPOCH)
                                                .map_err(|e| {
                                                    PipelineError::Other(anyhow!(
                                                        "Invalid mtime for running work unit {}: {}",
                                                        path.display(),
                                                        e
                                                    ))
                                                })?
                                                .as_secs();

                                        now.saturating_sub(modified_secs)
                                            > DISTRIBUTED_WORKER_STALE_SECS
                                    }
                                }
                            }

                            _ => false,
                        };

                    if !stale {
                        match work_unit.state {
                            crate::utils::work_units::WorkUnitState::Claimed => {
                                claimed += 1;
                            }

                            crate::utils::work_units::WorkUnitState::Running => {
                                running += 1;
                            }

                            _ => {}
                        }

                        continue;
                    }

                    // Acquire the same lock used by workers before modifying
                    // the work unit. Re-read the WorkUnit after acquiring the
                    // lock because the state may have changed since the first
                    // read above.
                    let lock_path =
                        path.with_extension("json.lock");

                    let lock_file =
                        match tokio::fs::OpenOptions::new()
                            .write(true)
                            .create_new(true)
                            .open(&lock_path)
                            .await
                        {
                            Ok(file) => file,

                            Err(err)
                            if err.kind()
                                == std::io::ErrorKind::AlreadyExists =>
                                {
                                    match work_unit.state {
                                        crate::utils::work_units::WorkUnitState::Claimed => {
                                            claimed += 1;
                                        }

                                        crate::utils::work_units::WorkUnitState::Running => {
                                            running += 1;
                                        }

                                        _ => {}
                                    }

                                    continue;
                                }

                            Err(err) => {
                                return Err(
                                    PipelineError::Other(anyhow!(
                                        "Failed to acquire recovery lock {}: {}",
                                        lock_path.display(),
                                        err
                                    ))
                                );
                            }
                        };

                    let recovery_result: Result<bool, PipelineError> = async {
                        let current_bytes =
                            tokio::fs::read(&path)
                                .await
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "Failed to reread work unit {} during recovery: {}",
                                        path.display(),
                                        e
                                    ))
                                })?;

                        let mut current_work_unit:
                            crate::utils::work_units::WorkUnit =
                            serde_json::from_slice(&current_bytes)
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "Failed to parse work unit {} during recovery: {}",
                                        path.display(),
                                        e
                                    ))
                                })?;

                        let now =
                            SystemTime::now()
                                .duration_since(UNIX_EPOCH)
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "System clock error during recovery: {}",
                                        e
                                    ))
                                })?
                                .as_secs();

                        let current_stale =
                            match current_work_unit.state {
                                crate::utils::work_units::WorkUnitState::Claimed => {
                                    let metadata =
                                        tokio::fs::metadata(&path)
                                            .await
                                            .map_err(|e| {
                                                PipelineError::Other(anyhow!(
                                                    "Failed to stat claimed work unit {} during recovery: {}",
                                                    path.display(),
                                                    e
                                                ))
                                            })?;

                                    let modified =
                                        metadata
                                            .modified()
                                            .map_err(|e| {
                                                PipelineError::Other(anyhow!(
                                                    "Failed to read claimed work-unit mtime {}: {}",
                                                    path.display(),
                                                    e
                                                ))
                                            })?;

                                    let modified_secs =
                                        modified
                                            .duration_since(UNIX_EPOCH)
                                            .map_err(|e| {
                                                PipelineError::Other(anyhow!(
                                                    "Invalid claimed work-unit mtime {}: {}",
                                                    path.display(),
                                                    e
                                                ))
                                            })?
                                            .as_secs();

                                    now.saturating_sub(modified_secs)
                                        > DISTRIBUTED_WORKER_STALE_SECS
                                }

                                crate::utils::work_units::WorkUnitState::Running => {
                                    match current_work_unit.last_heartbeat {
                                        Some(last_heartbeat) => {
                                            now.saturating_sub(last_heartbeat)
                                                > DISTRIBUTED_WORKER_STALE_SECS
                                        }

                                        None => false,
                                    }
                                }

                                _ => false,
                            };

                        if !current_stale {
                            return Ok(false);
                        }

                        let worker_id =
                            current_work_unit
                                .claimed_by
                                .clone()
                                .unwrap_or_else(|| {
                                    "unknown-worker".to_string()
                                });

                        let attempt =
                            current_work_unit.attempt;

                        let state =
                            current_work_unit.state;

                        let reason =
                            match state {
                                crate::utils::work_units::WorkUnitState::Claimed => {
                                    format!(
                                        "worker {} became stale while CLAIMED; \
                                         no RUNNING transition within {} seconds",
                                        worker_id,
                                        DISTRIBUTED_WORKER_STALE_SECS
                                    )
                                }

                                crate::utils::work_units::WorkUnitState::Running => {
                                    match current_work_unit.last_heartbeat {
                                        Some(last_heartbeat) => format!(
                                            "worker {} became stale while RUNNING; \
                                             last heartbeat was {} seconds ago",
                                            worker_id,
                                            now.saturating_sub(last_heartbeat)
                                        ),

                                        None => format!(
                                            "worker {} became stale while RUNNING; \
                                             no heartbeat recorded",
                                            worker_id
                                        ),
                                    }
                                }

                                _ => {
                                    return Ok(false);
                                }
                            };

                        warn!(
                            "Distributed NR: worker failure detected: \
                             work unit={} worker={} attempt={} state={:?}: {}",
                            current_work_unit.id(),
                            worker_id,
                            attempt,
                            state,
                            reason
                        );

                        current_work_unit
                            .fail(
                                &worker_id,
                                attempt,
                                reason,
                                true,
                            )
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Failed to mark stale work unit {} as FAILED: {}",
                                    current_work_unit.id(),
                                    e
                                ))
                            })?;

                        current_work_unit
                            .retry()
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Failed to requeue stale work unit {}: {}",
                                    current_work_unit.id(),
                                    e
                                ))
                            })?;

                        let payload =
                            serde_json::to_vec_pretty(
                                &current_work_unit,
                            )
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "Failed to serialize requeued work unit {}: {}",
                                        current_work_unit.id(),
                                        e
                                    ))
                                })?;

                        let tmp_path =
                            path.with_extension("json.tmp");

                        tokio::fs::write(
                            &tmp_path,
                            payload,
                        )
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Failed to write requeued work unit {}: {}",
                                    tmp_path.display(),
                                    e
                                ))
                            })?;

                        tokio::fs::rename(
                            &tmp_path,
                            &path,
                        )
                            .await
                            .map_err(|e| {
                                PipelineError::Other(anyhow!(
                                    "Failed to publish requeued work unit {}: {}",
                                    path.display(),
                                    e
                                ))
                            })?;

                        info!(
                            "Distributed NR: requeued work unit {} \
                             after stale worker {} attempt {}",
                            current_work_unit.id(),
                            worker_id,
                            attempt
                        );

                        Ok(true)
                    }
                        .await;

                    drop(lock_file);

                    if let Err(err) =
                        tokio::fs::remove_file(&lock_path).await
                    {
                        warn!(
                            "Distributed NR: failed to remove recovery lock {}: {}",
                            lock_path.display(),
                            err
                        );
                    }

                    if recovery_result? {
                        recovered += 1;
                        available += 1;
                    } else {
                        match work_unit.state {
                            crate::utils::work_units::WorkUnitState::Claimed => {
                                claimed += 1;
                            }

                            crate::utils::work_units::WorkUnitState::Running => {
                                running += 1;
                            }

                            _ => {}
                        }
                    }
                }

                crate::utils::work_units::WorkUnitState::Done => {
                    done += 1;
                }

                crate::utils::work_units::WorkUnitState::Failed => {
                    if work_unit
                        .failure
                        .as_ref()
                        .map(|failure| failure.retryable)
                        .unwrap_or(false)
                    {
                        let lock_path =
                            path.with_extension("json.lock");

                        let lock_file =
                            match tokio::fs::OpenOptions::new()
                                .write(true)
                                .create_new(true)
                                .open(&lock_path)
                                .await
                            {
                                Ok(file) => file,

                                Err(err)
                                if err.kind()
                                    == std::io::ErrorKind::AlreadyExists =>
                                    {
                                        failed += 1;
                                        continue;
                                    }

                                Err(err) => {
                                    return Err(
                                        PipelineError::Other(anyhow!(
                                            "Failed to acquire retry lock {}: {}",
                                            lock_path.display(),
                                            err
                                        ))
                                    );
                                }
                            };

                        let retry_result: Result<bool, PipelineError> = async {
                            let current_bytes =
                                tokio::fs::read(&path)
                                    .await
                                    .map_err(|e| {
                                        PipelineError::Other(anyhow!(
                                            "Failed to reread FAILED work unit {}: {}",
                                            path.display(),
                                            e
                                        ))
                                    })?;

                            let mut current_work_unit:
                                crate::utils::work_units::WorkUnit =
                                serde_json::from_slice(&current_bytes)
                                    .map_err(|e| {
                                        PipelineError::Other(anyhow!(
                                            "Failed to parse FAILED work unit {}: {}",
                                            path.display(),
                                            e
                                        ))
                                    })?;

                            if current_work_unit.state
                                != crate::utils::work_units::WorkUnitState::Failed
                            {
                                return Ok(false);
                            }

                            let retryable =
                                current_work_unit
                                    .failure
                                    .as_ref()
                                    .map(|failure| failure.retryable)
                                    .unwrap_or(false);

                            if !retryable {
                                return Ok(false);
                            }

                            let previous_attempt =
                                current_work_unit.attempt;

                            current_work_unit
                                .retry()
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "Failed to retry work unit {}: {}",
                                        current_work_unit.id(),
                                        e
                                    ))
                                })?;

                            let payload =
                                serde_json::to_vec_pretty(
                                    &current_work_unit,
                                )
                                    .map_err(|e| {
                                        PipelineError::Other(anyhow!(
                                            "Failed to serialize requeued work unit {}: {}",
                                            current_work_unit.id(),
                                            e
                                        ))
                                    })?;

                            let tmp_path =
                                path.with_extension("json.tmp");

                            tokio::fs::write(
                                &tmp_path,
                                payload,
                            )
                                .await
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "Failed to write requeued work unit {}: {}",
                                        tmp_path.display(),
                                        e
                                    ))
                                })?;

                            tokio::fs::rename(
                                &tmp_path,
                                &path
                            )
                                .await
                                .map_err(|e| {
                                    PipelineError::Other(anyhow!(
                                        "Failed to publish requeued work unit {}: {}",
                                        path.display(),
                                        e
                                    ))
                                })?;

                            info!(
                                "Distributed NR: requeued retryable FAILED work unit {} \
                                 from attempt {}",
                                current_work_unit.id(),
                                previous_attempt
                            );

                            Ok(true)
                        }
                            .await;

                        drop(lock_file);

                        if let Err(err) =
                            tokio::fs::remove_file(&lock_path).await
                        {
                            warn!(
                                "Distributed NR: failed to remove retry lock {}: {}",
                                lock_path.display(),
                                err
                            );
                        }

                        if retry_result? {
                            recovered += 1;
                            available += 1;
                        } else {
                            failed += 1;
                        }
                    } else {
                        failed += 1;
                    }
                }
            }
        }

        info!(
            "Distributed NR scheduler: AVAILABLE={} CLAIMED={} RUNNING={} DONE={} FAILED={} RECOVERED={}",
            available, claimed, running, done, failed, recovered
        );

        // Only non-retryable FAILED work units are terminal failures.
        if failed > 0 {
            return Err(
                PipelineError::Other(anyhow!(
                    "Distributed NR execution failed: {} work units are FAILED",
                    failed
                ))
            );
        }

        if done == chunk_work_units.len() {
            info!(
                "Distributed NR scheduling complete: all {} work units DONE",
                done
            );
            break;
        }

        tokio::time::sleep(
            tokio::time::Duration::from_secs(2)
        )
            .await;
    }

    // ------------------------------------------------------------------
    // 7. Collect completed result m8 files in deterministic chunk order.
    //
    // Each DONE WorkUnit contains the authoritative result metadata,
    // including the attempt number and durable result_path. Use that
    // metadata rather than discovering result files by filename alone,
    // because retries may leave multiple attempt files in results/.
    // ------------------------------------------------------------------

    let mut completed_work_units:
        Vec<crate::utils::work_units::WorkUnit> =
        Vec::with_capacity(
            chunk_work_units.len()
        );

    let mut entries =
        tokio::fs::read_dir(
            &work_dir
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "Failed to read distributed work directory {}: {}",
                    work_dir.display(),
                    e
                ))
            })?;

    while let Some(entry) =
        entries
            .next_entry()
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "Failed reading distributed work directory {}: {}",
                    work_dir.display(),
                    e
                ))
            })?
    {
        let path =
            entry.path();

        let Some(name) =
            path.file_name()
                .and_then(|v| v.to_str())
        else {
            continue;
        };

        if !name.starts_with("work_")
            || !name.ends_with(".json")
        {
            continue;
        }

        let bytes =
            tokio::fs::read(&path)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to read completed work unit {}: {}",
                        path.display(),
                        e
                    ))
                })?;

        let work_unit:
            crate::utils::work_units::WorkUnit =
            serde_json::from_slice(&bytes)
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to parse completed work unit {}: {}",
                        path.display(),
                        e
                    ))
                })?;

        if work_unit.state
            != crate::utils::work_units::WorkUnitState::Done
        {
            return Err(
                PipelineError::Other(anyhow!(
                    "Expected work unit {} to be DONE after scheduler completion, found {:?}",
                    work_unit.id(),
                    work_unit.state
                ))
            );
        }

        let result =
            work_unit
                .result
                .as_ref()
                .ok_or_else(|| {
                    PipelineError::Other(anyhow!(
                        "DONE work unit {} has no result metadata",
                        work_unit.id()
                    ))
                })?;

        if result.attempt
            != work_unit.attempt
        {
            return Err(
                PipelineError::Other(anyhow!(
                    "Result attempt mismatch for {}: work unit attempt={}, result attempt={}",
                    work_unit.id(),
                    work_unit.attempt,
                    result.attempt
                ))
            );
        }

        let result_path =
            &result.result_path;

        let metadata =
            tokio::fs::metadata(
                result_path
            )
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to stat result for {}: {}: {}",
                        work_unit.id(),
                        result_path.display(),
                        e
                    ))
                })?;

        if metadata.len()
            != result.result_bytes
        {
            return Err(
                PipelineError::Other(anyhow!(
                    "Result size mismatch for {} attempt {}: \
                     metadata says {} bytes, work unit records {} bytes",
                    work_unit.id(),
                    result.attempt,
                    metadata.len(),
                    result.result_bytes
                ))
            );
        }

        completed_work_units.push(
            work_unit
        );
    }

    completed_work_units.sort_by_key(
        |work_unit| work_unit.chunk_id
    );

    if completed_work_units.len()
        != chunk_work_units.len()
    {
        return Err(
            PipelineError::Other(anyhow!(
                "Distributed NR result collection found {} DONE work units, expected {}",
                completed_work_units.len(),
                chunk_work_units.len()
            ))
        );
    }

    for (expected_idx, work_unit)
    in completed_work_units.iter().enumerate()
    {
        let expected_chunk_id =
            chunk_work_units[expected_idx].0;

        if work_unit.chunk_id
            != expected_chunk_id
        {
            return Err(
                PipelineError::Other(anyhow!(
                    "Distributed NR result collection chunk mismatch: \
                     expected chunk {}, found {}",
                    expected_chunk_id,
                    work_unit.chunk_id
                ))
            );
        }
    }

    let result_temp_dir =
        choose_temp_dir(
            config.input_size,
            &config.ram_temp_dir,
            &config.args.nvme_scratch,
            2,
            false,
        )
            .await?;

    let merged_m8 =
        result_temp_dir
            .path()
            .join(
                "distributed_nr_merged.m8"
            );

    let mut merged =
        TokioFile::create(
            &merged_m8
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(
                    format!(
                        "Failed to create merged distributed NR m8 {}: {}",
                        merged_m8.display(),
                        e
                    )
                )
            })?;

    let mut total_result_rows =
        0u64;

    let mut total_result_bytes =
        0u64;

    for work_unit
    in &completed_work_units
    {
        let result =
            work_unit
                .result
                .as_ref()
                .expect(
                    "DONE work unit result checked above"
                );

        info!(
            "Distributed NR result collection: chunk={} attempt={} rows={} bytes={} path={}",
            work_unit.chunk_id,
            result.attempt,
            result.result_rows,
            result.result_bytes,
            result.result_path.display()
        );

        let mut input =
            TokioFile::open(
                &result.result_path
            )
                .await
                .map_err(|e| {
                    PipelineError::IOError(
                        format!(
                            "Failed to open distributed NR result {}: {}",
                            result.result_path.display(),
                            e
                        )
                    )
                })?;

        tokio::io::copy(
            &mut input,
            &mut merged,
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(
                    format!(
                        "Failed to append distributed NR result {} to {}: {}",
                        result.result_path.display(),
                        merged_m8.display(),
                        e
                    )
                )
            })?;

        total_result_rows +=
            result.result_rows;

        total_result_bytes +=
            result.result_bytes;
    }

    merged
        .flush()
        .await
        .map_err(|e| {
            PipelineError::IOError(
                format!(
                    "Failed to flush merged distributed NR m8 {}: {}",
                    merged_m8.display(),
                    e
                )
            )
        })?;

    drop(merged);

    let merged_metadata =
        TokioFile::open(
            &merged_m8
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(
                    format!(
                        "Failed to reopen merged distributed NR m8 {}: {}",
                        merged_m8.display(),
                        e
                    )
                )
            })?;

    let merged_size =
        tokio::fs::metadata(
            &merged_m8
        )
            .await
            .map_err(|e| {
                PipelineError::IOError(
                    format!(
                        "Failed to stat merged distributed NR m8 {}: {}",
                        merged_m8.display(),
                        e
                    )
                )
            })?
            .len();

    if merged_size
        != total_result_bytes
    {
        return Err(
            PipelineError::Other(anyhow!(
                "Distributed NR merged result size mismatch: \
                 merged={} bytes, sum of chunk results={} bytes",
                merged_size,
                total_result_bytes
            ))
        );
    }

    info!(
        "Distributed NR result collection complete: \
         {} chunks, {} rows, {} bytes -> {}",
        completed_work_units.len(),
        total_result_rows,
        total_result_bytes,
        merged_m8.display()
    );

    let rx =
        parse_lines(
            merged_metadata,
            &config,
            StreamDataType::JustBytes,
        )
            .await
            .map_err(|e| {
                PipelineError::Other(
                    e.into()
                )
            })?;

    info!(
        "Distributed NR execution complete: {} chunks processed",
        chunk_work_units.len()
    );

    Ok((
        rx,
        Vec::new(),
        Vec::new(),
        vec![result_temp_dir],
    ))
}





/// Aligns unmapped reads against NR database using Diamond.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `r1_path`: path to the R1 unmapped FASTQ file
/// * `r2_path_opt`: optional path to the R2 unmapped FASTQ file
///
/// # Returns
///
/// Result containing:
/// - receiver for alignment results (m8 format)
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - vector of temporary directories used
async fn diamond_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
) -> Result<
    (
        mpsc::Receiver<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    let temp_dir = choose_temp_dir(
        100_000_000_000,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        true,
    )
        .await?;

    let diamond_db = config
        .args
        .diamond_db
        .clone()
        .ok_or(PipelineError::MissingArgument(
            "diamond_db required for NR alignment".into(),
        ))?;

    let (db_prefix, prep_tasks) =
        diamond_index_prep(Some(diamond_db), "non_host").await?;

    cleanup_tasks.extend(prep_tasks);

    // Preserve the existing single-node thread-allocation policy.
    let available_ram_gb =
        config.available_ram as f64 / 1_000_000_000.0;

    let threads = if available_ram_gb < 64.0 {
        config
            .thread_allocation(DIAMOND_TAG, None)
            .min(6)
            .max(1)
    } else {
        config.thread_allocation(DIAMOND_TAG, None)
    };

    info!(
        "Diamond non-host: threads={}, tmpdir={}",
        threads,
        temp_dir.path().display()
    );

    // Build the complete Diamond configuration through the shared
    // Diamond module. All Diamond-specific parameters are constructed
    // there using this machine's resources.
    let r1_config =
        crate::utils::command::diamond::build_diamond_config(
            db_prefix.clone(),
            r1_path,
            temp_dir.path().to_path_buf(),
            threads,
            config.available_ram,
            None,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: DIAMOND_TAG.to_string(),
                error: e.to_string(),
            })?;

    // Run R1.
    let r1_m8 = run_diamond_single_file(
        config.clone(),
        r1_config,
        &temp_dir,
        "nr_r1",
    )
        .await?;

    // Run R2 (if present), using the same shared Diamond parameter
    // construction and the same machine-local resource values.
    if let Some(r2_path) = r2_path_opt {
        let r2_config =
            crate::utils::command::diamond::build_diamond_config(
                db_prefix,
                r2_path,
                temp_dir.path().to_path_buf(),
                threads,
                config.available_ram,
                None,
            )
                .await
                .map_err(|e| PipelineError::ToolExecution {
                    tool: DIAMOND_TAG.to_string(),
                    error: e.to_string(),
                })?;

        let _r2_m8 = run_diamond_single_file(
            config.clone(),
            r2_config,
            &temp_dir,
            "nr_r2",
        )
            .await?;
    }

    // Merge results.
    let merged_m8 =
        temp_dir.path().join("diamond_merged_nr.m8");

    let mut merged = TokioFile::create(&merged_m8)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let mut f1 = TokioFile::open(&r1_m8).await?;

    tokio::io::copy(&mut f1, &mut merged).await?;

    merged.flush().await?;

    let m8_file = TokioFile::open(&merged_m8)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let m8_rx = parse_lines(
        m8_file,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok((
        m8_rx,
        cleanup_tasks,
        cleanup_receivers,
        vec![temp_dir],
    ))
}





/// Streaing parser if fasta to annotated
///
/// # Arguments
/// * `config` - RunConfig struct
/// * `duplicate_clusters` - Dashmap of duploicates
/// * `nt_map` - small map of only contig_id → best NT accession
/// * `n2_map` - small map of only contig_id → best NR accession
///
/// # Returns
/// Result of hashmap of id:accession
///
fn parse_fasta_batch_to_annotated_streaming(
    batch: Vec<u8>,
    duplicate_clusters: &Arc<DashMap<String, ClusterInfo>>,
    nt_map: &AHashMap<String, String>, // tiny contig_id → best NT accession
    nr_map: &AHashMap<String, String>, // tiny contig_id → best NR accession
) -> Vec<Vec<u8>> {
    let mut results = Vec::with_capacity(batch.len() / 512);
    let mut reader = match needletail::parse_fastx_reader(std::io::Cursor::new(batch)) {
        Ok(r) => r,
        Err(e) => {
            error!(
                "parse_fasta_batch_to_annotated_streaming: failed to create reader: {}",
                e
            );
            return results;
        }
    };

    while let Some(record) = reader.next() {
        let rec = match record {
            Ok(r) => r,
            Err(e) => {
                warn!(
                    "parse_fasta_batch_to_annotated_streaming: bad record: {}",
                    e
                );
                continue;
            }
        };

        let record: SequenceRecord = rec.into();
        let rep_id = record.id().to_string();

        let nr_acc = nr_map.get(&rep_id).cloned().unwrap_or_default();
        let nt_acc = nt_map.get(&rep_id).cloned().unwrap_or_default();
        let is_unidentified = nr_acc.is_empty() && nt_acc.is_empty();

        let dup_count = duplicate_clusters
            .get(&rep_id)
            .map(|r| r.size)
            .unwrap_or(1u64);

        if !is_unidentified {
            let id = format!("NR:{}:NT:{}:{}", nr_acc, nt_acc, rep_id);
            let mapped_rec = SequenceRecord::Fasta {
                id,
                desc: None,
                seq: record.seq(),
            };
            if let Ok(bytes) = mapped_rec.to_bytes() {
                let mut tagged = vec![0u8];
                tagged.extend_from_slice(&bytes);
                results.push(tagged);
            }
        } else {
            let fasta_rep = SequenceRecord::Fasta {
                id: rep_id.clone(),
                desc: None,
                seq: record.seq(),
            };
            if let Ok(bytes) = fasta_rep.to_bytes() {
                let mut tagged_unid = vec![1u8];
                tagged_unid.extend_from_slice(&bytes);
                results.push(tagged_unid.clone());

                if dup_count > 1 {
                    for _ in 1..dup_count {
                        results.push(tagged_unid.clone());
                    }
                }

                let mut tagged_uniq = vec![2u8];
                tagged_uniq.extend_from_slice(&bytes);
                results.push(tagged_uniq);
            }
        }
    }
    results
}


/// Short-read-specific annotated FASTA generation.
///
/// The pipeline-independent work is performed by
/// `utils::fastx::generate_annotated_fasta_stream`.
///
/// Existing short-read Rust behavior is preserved:
///
///   mapped:
///       NR:<nr_accessession>:NT:<nt_accession>:<read_id>
///
///   unidentified:
///       <read_id>
///
/// For each unidentified deduplicated representative:
///
///   - `unique_rx` receives the representative once.
///   - `expanded_rx` receives the representative `ClusterInfo.size` times.
///
///
/// The only behavioral change from the previous Rust implementation is that
/// an unidentified representative missing from `duplicate_clusters` is now
/// treated as an invariant violation instead of silently being treated as a
/// singleton.
///
async fn generate_annotated_fasta(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    nt_map: AHashMap<String, String>,
    nr_map: AHashMap<String, String>,
    _concurrency: usize,
) -> Result<(
    mpsc::Receiver<ParseOutput>, // mapped, accession-annotated
    mpsc::Receiver<ParseOutput>, // unidentified, duplicate-expanded
    mpsc::Receiver<ParseOutput>, // unidentified, unique representatives
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
)> {
    // Generic portion shared by short- and long-read pipelines.
    //
    // For short reads this preserves the previous behavior:
    //   mapped       -> accession-annotated ID
    //   unidentified -> original representative ID
    let (
        mapped_rx,
        generic_unidentified_rx,
        mut cleanup_tasks,
        mut cleanup_receivers,
    ) = generate_annotated_fasta_stream(
        config.clone(),
        input_stream,
        nt_map,
        nr_map,
    )
        .await?;

    let channel_cap = (config.base_buffer_size / 4).max(10_000);

    let (expanded_tx, expanded_rx) = mpsc::channel(channel_cap);
    let (unique_tx, unique_rx) = mpsc::channel(channel_cap);

    // This is tokio::sync::oneshot because `oneshot` is already imported
    // from tokio::sync at the top of short_read_mngs.rs.
    let (done_tx, done_rx) = oneshot::channel::<Result<()>>();

    let expansion_task = tokio::spawn(async move {
        let mut unidentified_stream =
            ReceiverStream::new(generic_unidentified_rx);

        while let Some(item) = unidentified_stream.next().await {
            let record = match item {
                ParseOutput::Fasta(record) => record,
                other => {
                    return Err(anyhow!(
                        "generate_annotated_fasta: expected FASTA record in \
                         unidentified stream, got {:?}",
                        other
                    ));
                }
            };

            let rep_id = record.id().to_string();

            // Every record reaching this point is a deduplicated
            // representative, so it must exist in duplicate_clusters.
            //
            // Copy the size out before any await so the DashMap guard
            // is not held across an asynchronous send.
            let dup_count = {
                let cluster = duplicate_clusters
                    .get(&rep_id)
                    .ok_or_else(|| {
                        anyhow!(
                            "generate_annotated_fasta: invariant violation: \
                             unidentified representative '{}' is missing \
                             from duplicate_clusters",
                            rep_id
                        )
                    })?;

                cluster.size
            };

            // One copy of each deduplicated representative goes to the
            // unique-unidentified output.
            unique_tx
                .send(ParseOutput::Fasta(record.clone()))
                .await
                .map_err(|_| {
                    anyhow!(
                        "generate_annotated_fasta: unique unidentified \
                         receiver dropped"
                    )
                })?;

            // Preserve the previous Rust behavior exactly: expand the
            // unidentified representative according to its duplicate-cluster
            // size, retaining the representative ID for every copy.
            for _ in 0..dup_count {
                expanded_tx
                    .send(ParseOutput::Fasta(record.clone()))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "generate_annotated_fasta: expanded unidentified \
                             receiver dropped"
                        )
                    })?;
            }
        }

        let _ = done_tx.send(Ok(()));

        info!(
            "Short-read unidentified duplicate expansion completed"
        );

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(expansion_task);
    cleanup_receivers.push(done_rx);

    Ok((
        mapped_rx,
        expanded_rx,
        unique_rx,
        cleanup_tasks,
        cleanup_receivers,
    ))
}



#[allow(dead_code)]
async fn write_dummy_assembly_files(assembly_dir: &PathBuf) -> Result<(), anyhow::Error> {
    let contigs = assembly_dir.join("contigs.fasta");
    let contigs_all = assembly_dir.join("contigs_all.fasta");
    let scaffolds = assembly_dir.join("scaffolds.fasta");
    let bowtie_sam = assembly_dir.join("read-contig.sam"); // or wherever you place it
    let stats_json = assembly_dir.join("contig_stats.json");

    let failed_marker = b";ASSEMBLY FAILED\n";
    let no_info = b"@NO INFO\n";
    let empty_json = b"{}";

    tokio::fs::write(&contigs, failed_marker).await?;
    tokio::fs::write(&contigs_all, failed_marker).await?;
    tokio::fs::write(&scaffolds, failed_marker).await?;
    tokio::fs::write(&bowtie_sam, no_info).await?;
    tokio::fs::write(&stats_json, empty_json).await?;

    info!("Wrote dummy assembly files after SPAdes failure");
    Ok(())
}

/// Runs SPAdes assembly and generates read-to-contig mapping.
///
/// This intentionally matches the failure semantics of the original Python
/// RunAssembly implementation:
///
/// - SPAdes failure is non-fatal to the overall pipeline.
/// - Missing expected SPAdes outputs are assembly failures.
/// - Contig filtering failures are assembly failures.
/// - Bowtie2/read-to-contig mapping failures are assembly failures.
/// - Contig-stat generation failures are assembly failures.
/// - On ANY such assembly failure, the assembly outputs are replaced with
///   dummy files and the overall pipeline continues.
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `assembly_out_dir` - Final <out_dir>/assembly directory.
/// * `r1_path` - Non-host R1 FASTQ.
/// * `r2_path_opt` - Optional non-host R2 FASTQ.
/// * `duplicate_clusters` - Duplicate cluster metadata.
/// * `paired` - Whether input is paired-end.
/// * `assembly_headroom` - Temp-space headroom multiplier.
///
/// # Returns
async fn process_assembly(
    config: Arc<RunConfig>,
    assembly_out_dir: &PathBuf,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    paired: bool,
    assembly_headroom: u64,
) -> Result<
    (
        CoverageOutputs,
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<()>>>,
        Vec<oneshot::Receiver<Result<()>>>,
        Vec<NamedTempFile>,
        TempDir,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();
    let temp_files: Vec<NamedTempFile> = Vec::new();

    let est_temp_bytes = config.input_size + MAX_SPADES_WORK_DIR;

    let temp_dir = choose_temp_dir(
        est_temp_bytes,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        assembly_headroom,
        true,
    )
        .await?;

    // ---------------------------------------------------------------------
    // Temp paths used internally by the Rust implementation
    // ---------------------------------------------------------------------

    let contigs_path = temp_dir.path().join("contigs.fasta");
    let contigs_all_path = temp_dir.path().join("contigs_all.fasta");
    let scaffolds_path = temp_dir.path().join("scaffolds.fasta");
    let bam_path = temp_dir.path().join("read-contig.bam");
    let stats_json_path = temp_dir.path().join("contig_stats.json");

    // ---------------------------------------------------------------------
    // Durable assembly outputs
    // ---------------------------------------------------------------------

    let final_contigs_path = assembly_out_dir.join("contigs.fasta");
    let final_contigs_all_path = assembly_out_dir.join("contigs_all.fasta");
    let final_scaffolds_path = assembly_out_dir.join("scaffolds.fasta");
    let final_bam_path = assembly_out_dir.join("read-contig.bam");
    let final_stats_json_path = assembly_out_dir.join("contig_stats.json");
    let stdout_log_path = assembly_out_dir.join("spades_stdout.log");

    // ---------------------------------------------------------------------
    // Helper result for the assembly-failure path.
    //
    // This is the Rust equivalent of the original Python except: block:
    //
    //   ;ASSEMBLY FAILED
    //   ;ASSEMBLY FAILED
    //   ;ASSEMBLY FAILED
    //   @NO INFO
    //   {}
    //
    // ---------------------------------------------------------------------

    macro_rules! return_assembly_failure {
        ($reason:expr) => {{
            let reason = $reason;

            error!(
                "[assembly] assembly failed; writing compatibility dummy outputs: {}",
                reason
            );

            let failed_marker = b";ASSEMBLY FAILED\n";
            let no_info = b"@NO INFO\n";
            let empty_json = b"{}";

            tokio::fs::write(&final_contigs_path, failed_marker)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to write dummy contigs.fasta after assembly failure: {}",
                        e
                    ))
                })?;

            tokio::fs::write(&final_contigs_all_path, failed_marker)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to write dummy contigs_all.fasta after assembly failure: {}",
                        e
                    ))
                })?;

            tokio::fs::write(&final_scaffolds_path, failed_marker)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to write dummy scaffolds.fasta after assembly failure: {}",
                        e
                    ))
                })?;

            tokio::fs::write(&final_bam_path, no_info)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to write dummy read-contig.bam after assembly failure: {}",
                        e
                    ))
                })?;

            tokio::fs::write(&final_stats_json_path, empty_json)
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Failed to write dummy contig_stats.json after assembly failure: {}",
                        e
                    ))
                })?;

            info!(
                "[assembly] compatibility dummy outputs written: \
                 contigs={}, contigs_all={}, scaffolds={}, bam={}, stats={}",
                final_contigs_path.display(),
                final_contigs_all_path.display(),
                final_scaffolds_path.display(),
                final_bam_path.display(),
                final_stats_json_path.display(),
            );

            let (empty_tx, empty_rx) = mpsc::channel(1);
            drop(empty_tx);

            return Ok((
                CoverageOutputs {
                    contigs_fasta: final_contigs_path,
                    contigs_all_fasta: final_contigs_all_path,
                    scaffolds_fasta: final_scaffolds_path,
                    bam_path: final_bam_path,
                    contig_stats_json: final_stats_json_path,
                    contig_stats: Arc::new(HashMap::new()),
                    read2contig: Arc::new(HashMap::new()),
                },
                ReceiverStream::new(empty_rx),
                cleanup_tasks,
                cleanup_receivers,
                temp_files,
                temp_dir,
            ));
        }};
    }

    // =====================================================================
    // SPAdes
    // =====================================================================

    info!(
        "[spades] starting: r1={}, r2={:?}, final_assembly_dir={}, temp_dir={}",
        r1_path.display(),
        r2_path_opt
            .as_ref()
            .map(|p| p.display().to_string()),
        assembly_out_dir.display(),
        temp_dir.path().display(),
    );

    let mut options = HashMap::new();
    options.insert("--only-assembler".to_string(), None);

    let spades_config = SpadesConfig {
        r1_path: r1_path.clone(),
        r2_path_opt: r2_path_opt.clone(),
        outdir_path: temp_dir.path().to_path_buf(),
        tempdir_path: Some(temp_dir.path().to_path_buf()),
        option_fields: options,
    };

    let spades_args = match generate_cli(
        SPADES_TAG,
        &config,
        Some(&spades_config),
    ) {
        Ok(args) => args,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed to generate SPAdes command: {}",
                e
            ));
        }
    };

    info!("[spades] command args: {:?}", spades_args);

    let (mut spades_child, spades_stderr_task) = match spawn_cmd(
        config.clone(),
        SPADES_TAG,
        spades_args,
        config.args.verbose,
        None,
    )
        .await
    {
        Ok(v) => v,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed to spawn SPAdes: {}",
                e
            ));
        }
    };

    info!(
        "[spades] child spawned: stdout_present={}, stderr_present={}",
        spades_child.stdout.is_some(),
        spades_child.stderr.is_some()
    );

    let spades_out_stream = match parse_child_output(
        &mut spades_child,
        ChildStream::Stdout,
        ParseMode::Lines,
        &config,
    )
        .await
    {
        Ok(stream) => stream,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed to capture SPAdes stdout: {}",
                e
            ));
        }
    };

    let spades_write_task =
        tokio::spawn(stream_to_file(spades_out_stream, stdout_log_path.clone()));

    cleanup_tasks.push(spades_write_task);

    // Match the original broad assembly exception behavior:
    // stderr handling problems also classify the assembly as failed.
    match spades_stderr_task.await {
        Ok(Ok(())) => {}
        Ok(Err(e)) => {
            return_assembly_failure!(format!(
                "SPAdes stderr task failed: {}",
                e
            ));
        }
        Err(e) => {
            return_assembly_failure!(format!(
                "SPAdes stderr task panicked: {}",
                e
            ));
        }
    }

    let spades_exit = match spades_child.wait().await {
        Ok(status) => status,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed waiting for SPAdes: {}",
                e
            ));
        }
    };

    info!(
        "[spades] exit status={:?} success={} stdout_log={}",
        spades_exit.code(),
        spades_exit.success(),
        stdout_log_path.display()
    );

    match fs::metadata(&stdout_log_path).await {
        Ok(meta) => info!(
            "[spades] stdout log written: {} bytes at {}",
            meta.len(),
            stdout_log_path.display()
        ),
        Err(e) => warn!(
            "[spades] stdout log missing or unreadable at {}: {}",
            stdout_log_path.display(),
            e
        ),
    }

    if !spades_exit.success() {
        return_assembly_failure!(format!(
            "SPAdes exited with status {:?}",
            spades_exit.code()
        ));
    }

    info!("[spades] process exited successfully");

    // =====================================================================
    // Validate required SPAdes outputs.
    //
    // IMPORTANT:
    // The original Python pipeline required BOTH contigs.fasta and
    // scaffolds.fasta. If either was missing, its broad except block marked
    // the entire assembly as failed.
    //
    // Therefore scaffolds.fasta is intentionally NOT optional here.
    // =====================================================================

    if !contigs_path.is_file() {
        return_assembly_failure!(format!(
            "SPAdes exited successfully but required output is missing: {}",
            contigs_path.display()
        ));
    }

    if !scaffolds_path.is_file() {
        return_assembly_failure!(format!(
            "SPAdes exited successfully but required output is missing: {}",
            scaffolds_path.display()
        ));
    }

    info!(
        "[spades] required outputs present: contigs={}, scaffolds={}",
        contigs_path.display(),
        scaffolds_path.display()
    );

    // =====================================================================
    // Preserve original contigs
    // =====================================================================

    if let Err(e) = fs::copy(&contigs_path, &contigs_all_path).await {
        return_assembly_failure!(format!(
            "failed to preserve SPAdes contigs {} -> {}: {}",
            contigs_path.display(),
            contigs_all_path.display(),
            e
        ));
    }

    // =====================================================================
    // Contig length filtering
    // =====================================================================

    let rx = match read_fasta(
        contigs_all_path.clone(),
        u64::MAX,
        Some(config.args.min_contig_length),
        None,
        &config,
    ) {
        Ok(rx) => rx,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed reading SPAdes contigs for length filtering: {}",
                e
            ));
        }
    };

    let write_handle = write_fasta_stream_to_file(
        ReceiverStream::new(rx),
        contigs_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "process_assembly_contigs_filtering",
    );

    match write_handle.await {
        Ok(Ok(())) => {}
        Ok(Err(e)) => {
            return_assembly_failure!(format!(
                "contig filtering writer failed: {}",
                e
            ));
        }
        Err(e) => {
            return_assembly_failure!(format!(
                "contig filtering writer task panicked: {}",
                e
            ));
        }
    }

    // =====================================================================
    // Durable copies
    // =====================================================================

    if let Err(e) = tokio::fs::copy(
        &contigs_path,
        &final_contigs_path,
    )
        .await
    {
        return_assembly_failure!(format!(
            "failed copying filtered contigs {} -> {}: {}",
            contigs_path.display(),
            final_contigs_path.display(),
            e
        ));
    }

    if let Err(e) = tokio::fs::copy(
        &contigs_all_path,
        &final_contigs_all_path,
    )
        .await
    {
        return_assembly_failure!(format!(
            "failed copying original contigs {} -> {}: {}",
            contigs_all_path.display(),
            final_contigs_all_path.display(),
            e
        ));
    }

    if let Err(e) = tokio::fs::copy(
        &scaffolds_path,
        &final_scaffolds_path,
    )
        .await
    {
        return_assembly_failure!(format!(
            "failed copying scaffolds {} -> {}: {}",
            scaffolds_path.display(),
            final_scaffolds_path.display(),
            e
        ));
    }

    info!(
        "[spades] assembly FASTA outputs copied successfully: \
         contigs={}, contigs_all={}, scaffolds={}",
        final_contigs_path.display(),
        final_contigs_all_path.display(),
        final_scaffolds_path.display()
    );

    // =====================================================================
    // Bowtie2 index
    //
    // This is ALSO inside the original Python try/except.
    // Failure here must therefore become ASSEMBLY FAILED rather than killing
    // the entire short-read pipeline.
    // =====================================================================

    let index_dir = temp_dir.path().join("bowtie_index");

    if let Err(e) = fs::create_dir_all(&index_dir).await {
        return_assembly_failure!(format!(
            "failed creating Bowtie2 index directory {}: {}",
            index_dir.display(),
            e
        ));
    }

    let index_prefix = index_dir.join("contigs");

    let num_index_cores =
        RunConfig::thread_allocation(&*config, BOWTIE2_TAG, None);

    info!(
        "[assembly] building bowtie2 index: prefix={}, threads={}",
        index_prefix.display(),
        num_index_cores
    );

    let threads_string = num_index_cores.to_string();

    let contigs_str = match contigs_path.to_str() {
        Some(v) => v,
        None => {
            return_assembly_failure!(format!(
                "contigs path is not valid UTF-8: {}",
                contigs_path.display()
            ));
        }
    };

    let index_prefix_str = match index_prefix.to_str() {
        Some(v) => v,
        None => {
            return_assembly_failure!(format!(
                "Bowtie2 index path is not valid UTF-8: {}",
                index_prefix.display()
            ));
        }
    };

    let build_status = match Command::new("bowtie2-build")
        .args([
            "--threads",
            threads_string.as_str(),
            "--quiet",
            contigs_str,
            index_prefix_str,
        ])
        .status()
        .await
    {
        Ok(status) => status,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed to spawn bowtie2-build: {}",
                e
            ));
        }
    };

    if !build_status.success() {
        return_assembly_failure!(format!(
            "bowtie2-build exited with status {}",
            build_status
        ));
    }

    info!("[assembly] bowtie2-build completed successfully");

    // =====================================================================
    // Bowtie2 read-to-contig alignment
    // =====================================================================

    let mut bt2_config_view = Bowtie2Config::new()
        .flag("--very-sensitive")
        .option(
            "-x",
            index_prefix.to_string_lossy().to_string(),
        )
        .option(
            "-p",
            config
                .thread_allocation(BOWTIE2_TAG, None)
                .to_string(),
        );

    bt2_config_view = if paired {
        let r2_path = match r2_path_opt.as_ref() {
            Some(path) => path,
            None => {
                return_assembly_failure!(
                "Paired Bowtie2 alignment requested but R2 path is missing"
                    .to_string()
            );
            }
        };

        bt2_config_view
            .option(
                "-1",
                r1_path.to_string_lossy().to_string(),
            )
            .option(
                "-2",
                r2_path.to_string_lossy().to_string(),
            )
    } else {
        bt2_config_view.option(
            "-U",
            r1_path.to_string_lossy().to_string(),
        )
    };

    let bt2_args = match generate_cli(
        BOWTIE2_TAG,
        &config,
        Some(&bt2_config_view),
    ) {
        Ok(args) => args,
        Err(e) => {
            return_assembly_failure!(format!(
            "failed generating Bowtie2 arguments: {}",
            e
        ));
        }
    };
    info!(
        "[assembly] bowtie2 args generated: {:?}",
        bt2_args
    );

    let (mut bt2_child, bt2_err_task) = match spawn_cmd(
        config.clone(),
        BOWTIE2_TAG,
        bt2_args,
        config.args.verbose,
        None,
    )
        .await
    {
        Ok(v) => v,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed spawning Bowtie2: {}",
                e
            ));
        }
    };

    cleanup_tasks.push(bt2_err_task);

    let bt2_out_stream = match parse_child_output(
        &mut bt2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        &config,
    )
        .await
    {
        Ok(stream) => stream,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed reading Bowtie2 output: {}",
                e
            ));
        }
    };

    // =====================================================================
    // Fan out Bowtie2 BAM stream into three independent branches:
    //
    // 1. durable BAM output: ordering is irrelevant
    // 2. name-sorted BAM: used for read-to-contig statistics
    // 3. coordinate-sorted BAM: returned for downstream assembly coverage
    //
    // The two samtools sorts are started independently and therefore run
    // concurrently. Neither sort waits for the other.
    // =====================================================================

    let (assembly_bam_streams, assembly_bam_fanout_done) =
        match fanout_to_channels(
            ReceiverStream::new(bt2_out_stream),
            3,
            "process_assembly_bam",
            &config,
            StreamDataType::JustBytes,
        )
            .await
        {
            Ok(v) => v,
            Err(e) => {
                return_assembly_failure!(format!(
                    "failed splitting assembly BAM stream: {}",
                    e
                ));
            }
        };

    cleanup_receivers.push(assembly_bam_fanout_done);

    let mut assembly_bam_streams = assembly_bam_streams.into_iter();

    // ---------------------------------------------------------------------
    // Branch 1: durable BAM output.
    //
    // This branch is deliberately unsorted. It is an artifact only; no
    // downstream computation depends on its record ordering.
    // ---------------------------------------------------------------------

    let bam_for_file = match assembly_bam_streams.next() {
        Some(rx) => ReceiverStream::new(rx),
        None => {
            return_assembly_failure!(
                "assembly BAM file branch was not created".to_string()
            );
        }
    };

    let write_bam_task = match write_byte_stream_to_file(
        &final_bam_path,
        bam_for_file,
        config.clone(),
        StreamDataType::JustBytes,
        "process_assembly",
        false,
    )
        .await
    {
        Ok(task) => task,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed starting assembly BAM writer: {}",
                e
            ));
        }
    };

    cleanup_tasks.push(write_bam_task);

    // ---------------------------------------------------------------------
    // Branch 2: query-name sort for read-to-contig statistics.
    //
    // This preserves the existing contract of generate_info_from_bam_stream
    // and the legacy Python generate_info_from_sam() input ordering.
    // ---------------------------------------------------------------------

    let bam_for_name_sort = match assembly_bam_streams.next() {
        Some(rx) => rx,
        None => {
            return_assembly_failure!(
                "assembly BAM name-sort branch was not created".to_string()
            );
        }
    };

    let name_sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("sort"))
                .to_string(),
        )
        .option("-m", "2G")
        .flag("-n")
        .flag("-u")
        .option("-O", "bam")
        .option("-T", temp_dir.path().to_string_lossy().to_string())
        .positional("-");

    let name_sort_args = match generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&name_sort_config),
    ) {
        Ok(args) => args,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed generating samtools name-sort arguments: {}",
                e
            ));
        }
    };

    info!(
        "[assembly] starting query-name sort for read/contig statistics: {:?}",
        name_sort_args
    );

    let (
        name_sort_child,
        name_sort_task,
        name_sort_err_task,
    ) = match stream_to_cmd(
        config.clone(),
        bam_for_name_sort,
        SAMTOOLS_TAG,
        name_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
    {
        Ok(v) => v,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed starting assembly BAM name sort: {}",
                e
            ));
        }
    };

    cleanup_tasks.push(name_sort_task);
    cleanup_tasks.push(name_sort_err_task);

    let name_sort_rx = {
        let mut guard = name_sort_child.lock().await;

        match parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
        {
            Ok(rx) => rx,
            Err(e) => {
                return_assembly_failure!(format!(
                    "failed reading coordinate-independent name-sort output: {}",
                    e
                ));
            }
        }
    };

    // ---------------------------------------------------------------------
    // Branch 3: coordinate sort for streampile coverage.
    //
    // This is intentionally launched independently of the name sort above.
    // The two samtools processes now consume separate fanout branches in
    // parallel.
    // ---------------------------------------------------------------------

    let bam_for_coordinate_sort = match assembly_bam_streams.next() {
        Some(rx) => rx,
        None => {
            return_assembly_failure!(
                "assembly BAM coordinate-sort branch was not created".to_string()
            );
        }
    };

    let coordinate_sort_config = SamtoolsConfig::new(SamtoolsSubcommand::Sort)
        .option(
            "-@",
            config
                .thread_allocation(SAMTOOLS_TAG, Some("sort"))
                .to_string(),
        )
        .option("-m", "2G")
        .flag("-u")
        .option("-O", "bam")
        .option("-T", temp_dir.path().to_string_lossy().to_string())
        .positional("-");

    let coordinate_sort_args = match generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&coordinate_sort_config),
    ) {
        Ok(args) => args,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed generating samtools coordinate-sort arguments: {}",
                e
            ));
        }
    };

    info!(
        "[assembly] starting coordinate sort for assembly coverage: {:?}",
        coordinate_sort_args
    );

    let (
        coordinate_sort_child,
        coordinate_sort_task,
        coordinate_sort_err_task,
    ) = match stream_to_cmd(
        config.clone(),
        bam_for_coordinate_sort,
        SAMTOOLS_TAG,
        coordinate_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None,
    )
        .await
    {
        Ok(v) => v,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed starting assembly BAM coordinate sort: {}",
                e
            ));
        }
    };

    cleanup_tasks.push(coordinate_sort_task);
    cleanup_tasks.push(coordinate_sort_err_task);

    let coordinate_sort_rx = {
        let mut guard = coordinate_sort_child.lock().await;

        match parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
        {
            Ok(rx) => rx,
            Err(e) => {
                return_assembly_failure!(format!(
                    "failed reading coordinate-sort output: {}",
                    e
                ));
            }
        }
    };

    // This is the stream that will ultimately be returned from
    // process_assembly() and consumed by generate_assembly_coverage().
    let bam_for_output = ReceiverStream::new(coordinate_sort_rx);

    // =====================================================================
    // Contig statistics
    //
    // name_sort_rx is already being drained by parse_child_output() while
    // the coordinate sort above is running independently.
    // =====================================================================

    let (read2contig, contig_stats, _contig_uniques) =
        match generate_short_read_info_from_bam_stream(
            name_sort_rx,
            &duplicate_clusters,
            MIN_CONTIG_SIZE,
        )
            .await
        {
            Ok(v) => v,
            Err(e) => {
                return_assembly_failure!(format!(
                "failed generating read-to-contig mapping/statistics: {}",
                e
            ));
            }
        };

    info!(
    "[assembly] contig statistics generated: contigs={}, read2contig={}",
    contig_stats.len(),
    read2contig.len()
);

    let stats_json = match serde_json::to_string_pretty(&contig_stats) {
        Ok(json) => json,
        Err(e) => {
            return_assembly_failure!(format!(
                "failed serializing contig statistics: {}",
                e
            ));
        }
    };

    if let Err(e) = std::fs::write(
        &stats_json_path,
        &stats_json,
    ) {
        return_assembly_failure!(format!(
            "failed writing temporary contig_stats.json {}: {}",
            stats_json_path.display(),
            e
        ));
    }

    if let Err(e) = tokio::fs::write(
        &final_stats_json_path,
        &stats_json,
    )
        .await
    {
        return_assembly_failure!(format!(
            "failed writing final contig_stats.json {}: {}",
            final_stats_json_path.display(),
            e
        ));
    }

    info!(
        "[assembly] completed successfully"
    );

    Ok((
        CoverageOutputs {
            contigs_fasta: contigs_path,
            contigs_all_fasta: contigs_all_path,
            scaffolds_fasta: scaffolds_path,
            bam_path,
            contig_stats_json: stats_json_path,
            contig_stats: Arc::new(contig_stats),
            read2contig: Arc::new(read2contig),
        },
        bam_for_output,
        cleanup_tasks,
        cleanup_receivers,
        temp_files,
        temp_dir,
    ))
}


/// Short-read-specific wrapper around the generic BAM read-to-contig
/// bookkeeping.
///
/// The generic function performs the common work:
///   * global first-seen read-ID deduplication
///   * read_id -> contig mapping
///   * unique read counts per contig
///   * unweighted contig statistics
///   * minimum-contig-size filtering
///
/// Short reads additionally weight each representative read by the size of
/// its duplicate cluster, reproducing the legacy
/// `generate_info_from_sam(..., duplicate_cluster_sizes_path=...)` behavior.
pub async fn generate_short_read_info_from_bam_stream(
    rx: mpsc::Receiver<ParseOutput>,
    duplicate_clusters: &Arc<DashMap<String, ClusterInfo>>,
    min_contig_size: u64,
) -> Result<(
    HashMap<String, String>, // read2contig
    HashMap<String, u64>,    // cluster-adjusted contig_stats
    HashMap<String, usize>,  // contig_unique_counts
)> {
    let (
        read2contig,
        unweighted_contig_stats,
        contig_unique_counts,
    ) = generate_info_from_bam_stream(
        rx,
        min_contig_size,
    )
        .await?;

    // Start from zero and reconstruct the cluster-weighted contig counts
    // from the generic read-to-contig mapping.
    //
    // Only retained contigs are included. The generic function has already
    // applied min_contig_size to the unique-read counts.
    let retained_contigs: HashSet<String> =
        unweighted_contig_stats.keys().cloned().collect();

    let mut contig_stats: HashMap<String, u64> =
        unweighted_contig_stats
            .keys()
            .map(|contig| (contig.clone(), 0u64))
            .collect();

    for (read_id, contig_name) in &read2contig {
        if !retained_contigs.contains(contig_name) {
            continue;
        }

        let cluster_size = match duplicate_clusters.get(read_id) {
            Some(entry) => entry.value().size,

            None => {
                // Match legacy get_read_cluster_size():
                // if the exact read ID is absent, strip the final two
                // characters and try that key. The result must exist, and
                // the stripped suffix must be /1 or /2.
                if read_id.len() < 2 {
                    return Err(anyhow!(
                        "Read ID not found in duplicate cluster sizes dict: {}",
                        read_id
                    ));
                }

                let suffix = &read_id[read_id.len() - 2..];
                let prefix = &read_id[..read_id.len() - 2];

                if suffix != "/1" && suffix != "/2" {
                    return Err(anyhow!(
                        "Read ID not found in duplicate cluster sizes dict: {}",
                        read_id
                    ));
                }

                duplicate_clusters
                    .get(prefix)
                    .map(|entry| entry.value().size)
                    .ok_or_else(|| {
                        anyhow!(
                            "Read ID not found in duplicate cluster sizes dict: {}",
                            read_id
                        )
                    })?
            }
        };

        *contig_stats
            .entry(contig_name.clone())
            .or_insert(0) += cluster_size;
    }

    info!(
        "generate_short_read_info_from_bam_stream: {} mapped reads, {} retained contigs",
        read2contig.len(),
        contig_stats.len()
    );

    Ok((
        read2contig,
        contig_stats,
        contig_unique_counts,
    ))
}

/// Fixes FASTA header by removing leading comma in description if present.
/// apparently there are a few broken accessions out there
/// Handles multi-headers separated by \x01 (CTRL-A).
///
/// # Arguments
/// * `header` - header
/// # Returns
/// header as String
fn fix_header(header_line: &str) -> String {
    if !header_line.starts_with('>') {
        return header_line.to_owned();
    }

    // Split on \x01 (multi-accession entries), fix commas in each part, rejoin
    let parts: Vec<&str> = header_line
        .split(|c: char| {
            c == ','
                && header_line
                .as_bytes()
                .get(header_line.find(',').unwrap_or(0) + 1)
                .map_or(false, |&b| b != b' ')
        })
        .collect();

    format!(">{}", parts.join("\x01"))
}

/// Extracts accessions from a FASTA, then uses an offset FST db to locatre
/// offsets in a FASTa. Then extracts those accessions + seqs and
/// writes them to a file, ideally to RAM
/// a
///
/// # Arguments
/// * `config` - RunConfig struct from main.
/// * `selelcted_genera` - Hashmap of
/// * `fasta_path` - path to fulle fasta (e.g. NT.fa) where all accession/seqs stored
/// * `index_path` - path to FST offset index associated witht erh abover fasta_path
/// * `out_path
/// # Returns
/// Result of hashset of the accessions
async fn build_reference_fasta_from_selected_genera(
    config: Arc<RunConfig>,
    selected_genera: &HashMap<i32, Vec<String>>, // genus_taxid → [accessions]
    db_type: &str,
    source_fasta_path: &PathBuf, // full NT/NR FASTA (e.g. nt.fa)
    source_index_path: &PathBuf, // FST index for fast lookup
) -> Result<(PathBuf, TempDir)> {
    // ← Change return to include TempDir
    // Flatten and deduplicate accessions
    let mut accessions: HashSet<String> = HashSet::new();
    for acc_list in selected_genera.values() {
        accessions.extend(acc_list.iter().cloned());
    }

    let temp_dir = choose_temp_dir(
        1024,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?; // ← Create early

    if accessions.is_empty() {
        warn!("selected_genera is empty — writing empty reference FASTA");
        let empty_path = temp_dir.path().join(format!("{}_empty_ref.fasta", db_type));

        // Create and write minimal FASTA
        let file = TokioFile::create(&empty_path).await.with_context(|| {
            format!("Failed to create empty ref FASTA: {}", empty_path.display())
        })?;
        let mut writer = BufWriter::new(file);
        writer.write_all(b">empty\nN\n").await.with_context(|| {
            format!(
                "Failed to write to empty ref FASTA: {}",
                empty_path.display()
            )
        })?;
        writer.flush().await?;
        drop(writer);

        // Verify
        if !empty_path.exists() {
            return Err(anyhow!(
                "Empty ref FASTA {} was not created after write",
                empty_path.display()
            ));
        }

        info!(
            "Wrote empty ref FASTA to {} (verified exists)",
            empty_path.display()
        );
        return Ok((empty_path, temp_dir)); // ← Return dir to keep alive
    }

    let mut acc_vec: Vec<String> = accessions.into_iter().collect();
    acc_vec.sort(); // deterministic order

    info!(
        "Building reference FASTA from {} filtered accessions (selected_genera)",
        acc_vec.len()
    );

    // Load FST index
    let index_data = std::fs::read(source_index_path)
        .with_context(|| format!("Failed to read FST index: {}", source_index_path.display()))?;
    let index = fst::Map::new(index_data)
        .with_context(|| format!("Failed to load FST index: {}", source_index_path.display()))?;

    // Open source FASTA
    let file = File::open(source_fasta_path).with_context(|| {
        format!(
            "Failed to open source FASTA: {}",
            source_fasta_path.display()
        )
    })?;
    let mut reader = BufReader::new(file);

    // Estimate size (use existing temp_dir)

    let output_path = temp_dir
        .path()
        .join(format!("{}_ref_from_selected_genera.fasta", db_type));
    let mut writer = std::io::BufWriter::new(
        File::create(&output_path)
            .with_context(|| format!("Failed to create ref FASTA: {}", output_path.display()))?,
    );

    let mut missing = 0;
    let mut skipped_long = 0;

    for acc in &acc_vec {
        let key = acc.as_bytes();
        if let Some(offset) = index.get(key) {
            reader.seek(SeekFrom::Start(offset))?;

            let mut buf = String::with_capacity(1024);
            buf.clear();
            if reader.read_line(&mut buf)? == 0 {
                warn!("EOF while reading header for {}", acc);
                missing += 1;
                continue;
            }

            let header = fix_header(&buf);
            let mut seq_len: u64 = header.len() as u64;
            writer.write_all(header.as_bytes())?;
            writer.write_all(b"\n")?;

            buf.clear();
            loop {
                let bytes_read = reader.read_line(&mut buf)?;
                if bytes_read == 0 || buf.starts_with('>') {
                    // Rewind so next accession can read its header
                    let rewind = buf.len() as i64;
                    reader.seek(SeekFrom::Current(-rewind))?;
                    break;
                }

                seq_len += bytes_read as u64;
                if seq_len > 100_000 {
                    skipped_long += 1;
                    info!("Skipping very long sequence: {} ({} bp)", acc, seq_len);
                    // Drain rest of sequence
                    while reader.read_line(&mut buf)? > 0 && !buf.starts_with('>') {
                        buf.clear();
                    }
                    break;
                }

                writer.write_all(buf.as_bytes())?;
                buf.clear();
            }
        } else {
            missing += 1;
            // This is normal — selected_genera may include accessions not in current DB
        }
    }

    writer.flush()?;
    drop(writer);

    if missing > 0 {
        warn!(
            "{} accessions from selected_genera not found in source FASTA",
            missing
        );
    }
    if skipped_long > 0 {
        warn!("{} sequences skipped (too long)", skipped_long);
    }

    info!(
        "Reference FASTA built: {} → {} accessions extracted",
        output_path.display(),
        acc_vec.len() - missing
    );

    Ok((output_path, temp_dir)) // ← Keep dir alive
}

/// Summarizes an M8 stream as a
/// # Arguments
/// * `stream` - stream in m8 format
/// * `min_reads_per_genus` - Minimum threshold of reads per genus to be included
/// # Returns
/// read dicttionary of read id's -> Readhit
///  accesions -> Accesiohit
async fn summarize_hits(
    config: Arc<RunConfig>,
    mut stream: ReceiverStream<ReducedRead>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    min_reads_per_genus: usize,
) -> Result<(
    AHashMap<String, Arc<ReadHit>>,
    AHashMap<String, AccessionHit>,
    HashMap<i32, Vec<String>>,
    usize,
)> {
    let concurrency = compute_phase_concurrency(&config, "summarize_hits", 0.35, 4.0, 96, 8);

    let (batch_tx, batch_rx) = mpsc::channel::<Vec<ReducedRead>>(concurrency * 2);
    let batch_rx = Arc::new(tokio::sync::Mutex::new(batch_rx));

    let producer = tokio::spawn(async move {
        let mut batch = Vec::with_capacity(8192);
        while let Some(item) = stream.next().await {
            batch.push(item);
            if batch.len() >= 8192 {
                batch_tx
                    .send(std::mem::take(&mut batch))
                    .await
                    .map_err(|_| anyhow!("summarize_hits: batch channel closed"))?;
            }
        }
        if !batch.is_empty() {
            batch_tx
                .send(batch)
                .await
                .map_err(|_| anyhow!("summarize_hits: batch channel closed on flush"))?;
        }
        Ok::<(), anyhow::Error>(())
    });

    let final_read_dict = Arc::new(DashMap::with_capacity(8_000_000));
    let final_accession_dict = Arc::new(DashMap::with_capacity(2_000_000));
    let final_genus_read = Arc::new(DashMap::<i32, usize>::new());
    let final_genus_species = Arc::new(DashMap::<i32, HashSet<i32>>::new());
    let final_genus_accessions = Arc::new(DashMap::<i32, HashSet<String>>::new());

    let total_reads = Arc::new(AtomicUsize::new(0));

    let mut worker_handles = Vec::with_capacity(concurrency);
    for _ in 0..concurrency {
        let rx = batch_rx.clone();
        let dup_clusters = duplicate_clusters.clone();

        let f_read = final_read_dict.clone();
        let f_acc = final_accession_dict.clone();
        let f_gr = final_genus_read.clone();
        let f_gs = final_genus_species.clone();
        let f_ga = final_genus_accessions.clone();
        let tr = total_reads.clone();

        worker_handles.push(tokio::spawn(async move {
            let mut local_read_dict: AHashMap<String, Arc<ReadHit>> = AHashMap::with_capacity(8192);
            let mut local_accession_dict: AHashMap<String, AccessionHit> =
                AHashMap::with_capacity(4096);
            let mut local_genus_read: AHashMap<i32, usize> = AHashMap::with_capacity(2048);
            let mut local_genus_species: AHashMap<i32, HashSet<i32>> =
                AHashMap::with_capacity(1024);
            let mut local_genus_accessions: AHashMap<i32, HashSet<String>> =
                AHashMap::with_capacity(1024);

            loop {
                let batch = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };

                let batch = match batch {
                    Some(b) => b,
                    None => break,
                };

                // HitSummary from ReducedRead.summary:
                // read_id, level, taxid, accession, species, genus, family
                let parsed: Vec<(String, u8, i32, String, i32, i32, i32)> =
                    tokio::task::spawn_blocking({
                        let batch = batch;
                        move || {
                            batch
                                .par_iter()
                                .filter_map(|item| {
                                    let line =
                                        std::str::from_utf8(&item.summary).ok()?.trim_end();
                                    if line.is_empty() {
                                        return None;
                                    }

                                    let parts: Vec<&str> = line.split('\t').collect();
                                    if parts.len() < 7 {
                                        return None;
                                    }

                                    let read_id = parts[0].to_string();
                                    let level: u8 = parts[1].parse().unwrap_or(0);
                                    let taxid: i32 = parts[2].parse().unwrap_or(0);
                                    let accession_id = if parts[3].is_empty()
                                        || parts[3] == "-"
                                        || parts[3] == "None"
                                    {
                                        "-".to_string()
                                    } else {
                                        parts[3].to_string()
                                    };
                                    let species_taxid: i32 = parts[4].parse().unwrap_or(-1);
                                    let genus_taxid: i32 = parts[5].parse().unwrap_or(-1);
                                    let family_taxid: i32 = parts[6].parse().unwrap_or(-1);

                                    Some((
                                        read_id,
                                        level,
                                        taxid,
                                        accession_id,
                                        species_taxid,
                                        genus_taxid,
                                        family_taxid,
                                    ))
                                })
                                .collect()
                        }
                    })
                        .await
                        .map_err(|e| anyhow!("summarize_hits batch parsing panicked: {}", e))?;

                for (
                    read_id,
                    level,
                    taxid,
                    accession_id,
                    species_taxid,
                    genus_taxid,
                    family_taxid,
                ) in parsed
                {
                    tr.fetch_add(1, AtomicOrdering::Relaxed);

                    let cluster_size: u64 = dup_clusters
                        .get(&read_id)
                        .map(|e| e.value().size)
                        .unwrap_or(1);

                    local_read_dict.insert(
                        read_id.clone(),
                        Arc::new(ReadHit {
                            level,
                            taxid,
                            accession_id: accession_id.clone(),
                            species_taxid,
                            genus_taxid,
                            family_taxid,
                            contig_id: None,
                            contig_accession_id: None,
                            contig_species_taxid: 0,
                            contig_genus_taxid: 0,
                            contig_family_taxid: 0,
                            from_assembly: false,
                            source_count_type: None,
                        }),
                    );

                    if accession_id != "-" && genus_taxid > 0 {
                        local_accession_dict
                            .entry(accession_id.clone())
                            .and_modify(|e| e.count += cluster_size)
                            .or_insert(AccessionHit {
                                taxid: species_taxid,
                                genus_taxid,
                                family_taxid,
                                count: cluster_size,
                            });

                        *local_genus_read.entry(genus_taxid).or_insert(0) += 1;
                        local_genus_species
                            .entry(genus_taxid)
                            .or_default()
                            .insert(species_taxid);
                        local_genus_accessions
                            .entry(genus_taxid)
                            .or_default()
                            .insert(accession_id);
                    }
                }
            }

            for (k, v) in local_read_dict {
                f_read.insert(k, v);
            }
            for (k, v) in local_accession_dict {
                f_acc
                    .entry(k)
                    .and_modify(|e: &mut AccessionHit| e.count += v.count)
                    .or_insert(v);
            }
            for (g, c) in local_genus_read {
                *f_gr.entry(g).or_insert(0) += c;
            }
            for (g, s) in local_genus_species {
                f_gs.entry(g).or_default().extend(s);
            }
            for (g, a) in local_genus_accessions {
                f_ga.entry(g).or_default().extend(a);
            }

            Ok::<(), anyhow::Error>(())
        }));
    }

    producer.await??;
    for h in worker_handles {
        h.await??;
    }

    let total_reads_val = total_reads.load(AtomicOrdering::SeqCst);

    info!(
        "[summarize_hits] processed {} reads, {} unique accessions ({} workers)",
        total_reads_val,
        final_accession_dict.len(),
        concurrency
    );

    let mut selected_genera: HashMap<i32, Vec<String>> = HashMap::new();
    for entry in final_genus_read.iter() {
        let genus_taxid = *entry.key();
        let read_count = *entry.value();

        if read_count < min_reads_per_genus {
            continue;
        }

        let species_count = final_genus_species
            .get(&genus_taxid)
            .map(|s| s.len())
            .unwrap_or(0);
        if species_count <= 1 {
            continue;
        }

        if let Some(acc_set) = final_genus_accessions.get(&genus_taxid) {
            let accessions: Vec<String> = acc_set.iter().cloned().collect();
            selected_genera.insert(genus_taxid, accessions);
        }
    }

    let read_dict_final: AHashMap<String, Arc<ReadHit>> = Arc::try_unwrap(final_read_dict)
        .map_err(|_| anyhow!("Failed to unwrap final_read_dict"))?
        .into_iter()
        .collect();

    let accession_dict_final: AHashMap<String, AccessionHit> =
        Arc::try_unwrap(final_accession_dict)
            .map_err(|_| anyhow!("Failed to unwrap final_accession_dict"))?
            .into_iter()
            .collect();

    Ok((
        read_dict_final,
        accession_dict_final,
        selected_genera,
        total_reads_val,
    ))
}


/// Updates the read dictionary based on top BLAST m8 results for contigs.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `read2contig`: map from read ID to contig ID
/// * `top_m8_stream`: stream of top m8 results for contigs
/// * `read_dict`: shared map of read IDs to their taxonomic hits
/// * `lineage_map`: map of taxids to lineages
/// * `accession_map`: map of accessions to hit metadata
/// * `should_keep`: filter for taxids to keep
/// * `db_type`: database type (e.g., "nt", "nr")
/// * `contig2lineage_tx`: sender for contig-to-lineage mapping records
/// * `read2blastm8_tx`: sender for refined m8 records
/// * `updated_tx`: sender for records of reads with updated hits
/// * `added_tx`: sender for records of reads with new hits
///
/// # Returns
///
/// Result<()>: success or error
pub async fn update_read_dict(
    config: Arc<RunConfig>,
    read2contig: Arc<HashMap<String, String>>,
    mut top_m8_stream: ReceiverStream<ParseOutput>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    _lineage_map: Arc<AHashMap<Taxid, Lineage>>, // unused: Python uses accession_dict only
    accession_map: Arc<AHashMap<String, AccessionHit>>,
    _should_keep: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>, // not applied here (Python)
    db_type: &str,
    contig2lineage_tx: Sender<ParseOutput>,
    read2blastm8_tx: Sender<ParseOutput>,
    updated_tx: Sender<ParseOutput>,
    added_tx: Sender<ParseOutput>,
) -> Result<()> {
    let _ = compute_phase_concurrency(&config, "update_read_dict", 0.05, 8.0, 128, 4);

    fn merge_m8_rows(curr: &mut M8Record, row: &M8Record) {
        let curr_len = curr.alen as f64;
        let row_len = row.alen as f64;
        let denom = curr_len + row_len;
        if denom > 0.0 {
            curr.pident = (curr.pident * curr_len + row.pident * row_len) / denom;
            curr.evalue = (curr.evalue * curr_len + row.evalue * row_len) / denom;
        }
        curr.mismatch += row.mismatch;
        curr.gapopen += row.gapopen;
        curr.bitscore += row.bitscore;
        curr.alen += row.alen;
        curr.qstart = curr.qstart.min(row.qstart);
        curr.qend = curr.qend.max(row.qend);
        curr.tstart = curr.tstart.min(row.tstart);
        curr.tend = curr.tend.max(row.tend);
    }

    // Python accession_dict keys are hit_summary accessions; BLAST sseqid may be versioned.
    fn lookup_acc<'a>(
        map: &'a AHashMap<String, AccessionHit>,
        acc: &str,
    ) -> Option<(&'a String, &'a AccessionHit)> {
        if let Some((k, h)) = map.get_key_value(acc) {
            return Some((k, h));
        }
        if let Some(base) = acc.split('.').next() {
            if base != acc {
                if let Some((k, h)) = map.get_key_value(base) {
                    return Some((k, h));
                }
            }
        }
        None
    }

    let db_type_upper = db_type.to_uppercase();

    // ── loop 1: blast top → contig2accession + contig2lineage ──────
    let mut contig2accession: AHashMap<String, (String, M8Record)> = AHashMap::new();
    let mut contig2lineage: AHashMap<String, [i32; 3]> = AHashMap::new();
    let mut tops_seen = 0u64;
    let mut acc_miss = 0u64;
    const ACC_MISS_LOG_CAP: u64 = 20; // first 20 only

    while let Some(item) = top_m8_stream.next().await {
        let ParseOutput::Bytes(bytes) = item else { continue };
        let line = String::from_utf8_lossy(&bytes);
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }

        let m8 = if db_type.eq_ignore_ascii_case("nt") {
            M8Record::parse_line_nt(line)?
        } else {
            M8Record::parse_line_nr(line)?
        };
        tops_seen += 1;

        let contig_id = m8.qname.clone();
        let accession_id = m8.tname.clone();

        // contig2lineage[contig] = accession_dict[accession]  (must exist)
        let Some((_key, acc_hit)) = lookup_acc(&accession_map, &accession_id) else {
            acc_miss += 1;
            if acc_miss <= ACC_MISS_LOG_CAP {
                warn!(
                "[update_read_dict:{}] acc_miss contig={} acc={} (raw tname from top m8)",
                db_type, contig_id, accession_id
            );
            }
            // Do not insert contig without lineage — Python would KeyError
            continue;
        };
        let lineage = [acc_hit.taxid, acc_hit.genus_taxid, acc_hit.family_taxid];
        contig2lineage.insert(contig_id.clone(), lineage);

        if let Some((_, curr)) = contig2accession.get_mut(&contig_id) {
            merge_m8_rows(curr, &m8);
        } else {
            contig2accession.insert(contig_id, (accession_id, m8));
        }

        if tops_seen % 5_000 == 0 {
            info!(
            "[update_read_dict:{}] tops={} contigs={} lineage={} acc_miss={}",
            db_type,
            tops_seen,
            contig2accession.len(),
            contig2lineage.len(),
            acc_miss
        );
        }
    }

    info!(
    "[update_read_dict:{}] pass1 done: tops={} contigs={} lineage={} acc_miss={}",
    db_type,
    tops_seen,
    contig2accession.len(),
    contig2lineage.len(),
    acc_miss
);
    if acc_miss > ACC_MISS_LOG_CAP {
        warn!(
        "[update_read_dict:{}] acc_miss total={} (only first {} logged)",
        db_type, acc_miss, ACC_MISS_LOG_CAP
    );
    }

    // ── loop 2: for read_id, contig_id in read2contig ──────────────
    let existing_reads: AHashMap<String, Arc<ReadHit>> = {
        let g = read_dict
            .lock()
            .map_err(|e| anyhow!("read_dict lock poisoned: {e}"))?;
        g.clone()
    };

    let mut updated_entries: Vec<(String, Arc<ReadHit>)> = Vec::new();
    let mut added_entries: Vec<(String, Arc<ReadHit>)> = Vec::new();
    let mut read2blastm8_lines: Vec<Vec<u8>> = Vec::new();
    let mut reads_with_contig_hit = 0u64;
    let mut reads_acc_ok = 0u64;

    for (read_id, contig_id) in read2contig.iter() {
        let Some((accession, m8_row)) = contig2accession.get(contig_id) else {
            continue; // contig has no BLAST top
        };
        reads_with_contig_hit += 1;

        // Python: if accession and accession in accession_dict
        if let Some((_k, acc_hit)) = lookup_acc(&accession_map, accession) {
            reads_acc_ok += 1;
            let species_taxid = acc_hit.taxid;
            let genus_taxid = acc_hit.genus_taxid;
            let family_taxid = acc_hit.family_taxid;

            if let Some(existing) = existing_reads.get(read_id) {
                let mut hit: ReadHit = (**existing).clone();
                hit.taxid = species_taxid;
                hit.contig_id = Some(contig_id.clone());
                hit.contig_accession_id = Some(accession.clone());
                hit.contig_species_taxid = species_taxid;
                hit.contig_genus_taxid = genus_taxid;
                hit.contig_family_taxid = family_taxid;
                hit.from_assembly = true;
                hit.source_count_type = Some(db_type_upper.clone());
                updated_entries.push((read_id.clone(), Arc::new(hit)));
            } else {
                added_entries.push((
                    read_id.clone(),
                    Arc::new(ReadHit {
                        level: 1,
                        taxid: species_taxid,
                        accession_id: accession.clone(),
                        species_taxid,
                        genus_taxid,
                        family_taxid,
                        contig_id: Some(contig_id.clone()),
                        contig_accession_id: Some(accession.clone()),
                        contig_species_taxid: species_taxid,
                        contig_genus_taxid: genus_taxid,
                        contig_family_taxid: family_taxid,
                        from_assembly: true,
                        source_count_type: Some(db_type_upper.clone()),
                    }),
                ));
            }
        }

        // Python: if m8_row: read2blastm8[read_id] = m8_row  (qseqid forced to read_id later)
        let mut m8_for_read = m8_row.clone();
        m8_for_read.qname = read_id.clone();
        let mut m8_line = m8_for_read.to_tab_string();
        m8_line.push('\n');
        read2blastm8_lines.push(m8_line.into_bytes());
    }

    info!(
        "[update_read_dict:{}] pass2: reads_with_contig_hit={} acc_ok={} updated={} added={} blast_m8={}",
        db_type,
        reads_with_contig_hit,
        reads_acc_ok,
        updated_entries.len(),
        added_entries.len(),
        read2blastm8_lines.len()
    );

    // Emit contig2lineage
    for (contig_id, lin) in &contig2lineage {
        let line = format!("{}\t{}\t{}\t{}\n", contig_id, lin[0], lin[1], lin[2]);
        contig2lineage_tx
            .send(ParseOutput::Bytes(Bytes::from(line)))
            .await
            .map_err(|_| anyhow!("contig2lineage_tx closed"))?;
    }

    {
        let mut dict = read_dict
            .lock()
            .map_err(|e| anyhow!("read_dict lock poisoned: {e}"))?;
        for (id, hit) in &updated_entries {
            dict.insert(id.clone(), hit.clone());
        }
        for (id, hit) in &added_entries {
            dict.insert(id.clone(), hit.clone());
        }
    }

    for line in read2blastm8_lines {
        read2blastm8_tx
            .send(ParseOutput::Bytes(Bytes::from(line)))
            .await
            .map_err(|_| anyhow!("read2blastm8_tx closed"))?;
    }
    for (read_id, hit) in updated_entries {
        let mut line = hit.to_full_tab_line(&read_id);
        line.push('\n');
        updated_tx
            .send(ParseOutput::Bytes(Bytes::from(line)))
            .await
            .map_err(|_| anyhow!("updated_tx closed"))?;
    }
    for (read_id, hit) in added_entries {
        let mut line = hit.to_full_tab_line(&read_id);
        line.push('\n');
        added_tx
            .send(ParseOutput::Bytes(Bytes::from(line)))
            .await
            .map_err(|_| anyhow!("added_tx closed"))?;
    }

    info!(
        "[update_read_dict:{}] emitted contig2lineage={}",
        db_type,
        contig2lineage.len()
    );
    Ok(())
}

/// Generates a JSON summary of contig taxonomic assignments.
///
/// # Arguments
///
/// * `read2contig`: map from read ID to contig ID
/// * `contig2lineage`: map from contig ID to validated taxonomic lineage
/// * `read_dict`: shared map of read IDs to their taxonomic hits
/// * `db_type`: database type (e.g., "nt", "nr")
/// * `duplicate_clusters`: map of cluster information for duplicate weighting
/// * `min_contig_size`: minimum number of reads a contig must have to be included
/// * `output_tx`: sender for the generated JSON lines
///
/// # Returns
///
/// Result<()>: success or error
async fn generate_contig_summary_json(
    read2contig: Arc<HashMap<String, String>>,
    contig2lineage: AHashMap<String, [i32; 3]>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    db_type: &str,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    should_keep: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    min_contig_size: u64, // Python MIN_CONTIG_SIZE = 4
    output_tx: Sender<ParseOutput>,
) -> Result<()> {
    // taxid → contig → [unique_reads, nonunique_weighted]
    let mut genus_summary: AHashMap<i32, AHashMap<String, [u64; 2]>> = AHashMap::new();
    let mut species_summary: AHashMap<i32, AHashMap<String, [u64; 2]>> = AHashMap::new();

    let db_type_upper = db_type.to_uppercase();

    {
        let dict = read_dict
            .lock()
            .map_err(|e| anyhow!("read_dict lock poisoned: {}", e))?;

        for (read_id, hit_arc) in dict.iter() {
            let hit: &ReadHit = hit_arc.as_ref();

            let contig = read2contig
                .get(read_id)
                .cloned()
                .unwrap_or_else(|| "*".to_string());

            let (species_taxid, genus_taxid, family_taxid) = if contig != "*" {
                if let Some(lineage) = contig2lineage.get(&contig) {
                    (lineage[0], lineage[1], lineage[2])
                } else {
                    (hit.species_taxid, hit.genus_taxid, hit.family_taxid)
                }
            } else {
                (hit.species_taxid, hit.genus_taxid, hit.family_taxid)
            };

            if !should_keep(&[species_taxid, genus_taxid, family_taxid]) {
                continue;
            }

            let cluster_size = duplicate_clusters
                .get(read_id.as_str())
                .map(|e| e.value().size)
                .unwrap_or(1u64);

            let sp = species_summary
                .entry(species_taxid)
                .or_default()
                .entry(contig.clone())
                .or_insert([0, 0]);
            sp[0] += 1;
            sp[1] += cluster_size;

            let gn = genus_summary
                .entry(genus_taxid)
                .or_default()
                .entry(contig)
                .or_insert([0, 0]);
            gn[0] += 1;
            gn[1] += cluster_size;
        }
    }

    for (tax_level, summary_map) in [(1u8, species_summary), (2u8, genus_summary)] {
        for (taxid, contig_counts) in summary_map {
            if taxid <= 0 {
                continue;
            }

            let filtered: HashMap<String, u64> = contig_counts
                .into_iter()
                .filter(|(_, [unique, _])| *unique >= min_contig_size)
                .map(|(contig, [_, weighted])| (contig, weighted))
                .collect();

            if filtered.is_empty() {
                continue;
            }

            let entry = ContigSummaryEntry {
                taxid,
                tax_level,
                count_type: db_type_upper.clone(),
                contig_counts: filtered,
            };
            let line = serde_json::to_string(&entry)? + "\n";
            output_tx
                .send(ParseOutput::Bytes(Bytes::from(line)))
                .await
                .map_err(|_| anyhow!("contig_summary_tx dropped"))?;
        }
    }

    Ok(())
}

/// Extracts the first field from a tab-separated byte slice.
///
/// # Arguments
///
/// * `bytes`: the tab-separated line bytes
///
/// # Returns
///
/// String: the extracted first field
fn first_field(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes)
        .split('\t')
        .next()
        .unwrap_or("")
        .trim_end()
        .to_string()
}

/// Collects a stream of bytes into a map keyed by the first tab-delimited field.
///
/// # Arguments
///
/// * `stream`: stream of bytes items
///
/// # Returns
///
/// Result containing the map and the order in which keys were first seen
async fn collect_line_map(
    mut stream: ReceiverStream<ParseOutput>,
) -> Result<(AHashMap<String, Bytes>, Vec<String>)> {
    let mut map: AHashMap<String, Bytes> = AHashMap::new();
    let mut order: Vec<String> = Vec::new();

    while let Some(item) = stream.next().await {
        let ParseOutput::Bytes(bytes) = item else {
            continue;
        };

        let key = first_field(&bytes);
        if !map.contains_key(&key) {
            order.push(key.clone());
        }
        map.insert(key, bytes);
    }

    Ok((map, order))
}

/// Collects a stream of bytes into an ordered vector of rows.
///
/// # Arguments
///
/// * `stream`: stream of bytes items
///
/// # Returns
///
/// Result containing the vector of row bytes
#[allow(dead_code)]
async fn collect_rows(mut stream: ReceiverStream<ParseOutput>) -> Result<Vec<Bytes>> {
    let mut rows = Vec::new();

    while let Some(item) = stream.next().await {
        let ParseOutput::Bytes(bytes) = item else {
            continue;
        };
        rows.push(bytes);
    }

    Ok(rows)
}

/// Collects a stream of bytes into keyed rows, preserving first-seen key order.
///
/// # Arguments
///
/// * `stream`: stream of bytes items
///
/// # Returns
///
/// Result containing the vector of (key, bytes) pairs and the key order
async fn collect_keyed_rows(
    mut stream: ReceiverStream<ParseOutput>,
) -> Result<(Vec<(String, Bytes)>, Vec<String>)> {
    let mut rows: Vec<(String, Bytes)> = Vec::new();
    let mut order: Vec<String> = Vec::new();

    while let Some(item) = stream.next().await {
        let ParseOutput::Bytes(bytes) = item else {
            continue;
        };

        let key = first_field(bytes.as_ref());
        order.push(key.clone());
        rows.push((key, bytes));
    }

    Ok((rows, order))
}

/// Merges refined hit information into existing m8 and hit summary streams.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `updated_reads_stream`: stream of reads that were updated by contig hits
/// * `added_reads_stream`: stream of reads that were newly identified by contig hits
/// * `blast_hits_stream`: stream of raw BLAST hits for contigs
/// * `original_hit_summary_stream`: stream of original hit summary records
/// * `original_deduped_m8_stream`: stream of original deduped m8 records
/// * `refined_m8_tx`: sender for the refined m8 stream
/// * `refined_hit_summary_tx`: sender for the refined hit summary stream
///
/// # Returns
///
/// Result<()>: success or error
async fn generate_m8_and_hit_summary(
    config: Arc<RunConfig>,
    updated_reads_stream: ReceiverStream<ParseOutput>,
    added_reads_stream: ReceiverStream<ParseOutput>,
    blast_hits_stream: ReceiverStream<ParseOutput>,
    original_hit_summary_stream: ReceiverStream<ParseOutput>,
    original_deduped_m8_stream: ReceiverStream<ParseOutput>,
    refined_m8_tx: Sender<ParseOutput>,
    refined_hit_summary_tx: Sender<ParseOutput>,
) -> Result<()> {
    info!("[generate_m8_and_hit_summary] collecting inputs (5 streams)");

    let updated_fut = collect_line_map(updated_reads_stream);
    let added_fut = collect_line_map(added_reads_stream);
    let blast_hits_fut = collect_line_map(blast_hits_stream);
    let hit_rows_fut = collect_keyed_rows(original_hit_summary_stream);
    let m8_rows_fut = collect_keyed_rows(original_deduped_m8_stream);

    let (
        (updated_reads, _updated_order),
        (added_reads, added_order),
        (blast_hits, _blast_order),
        (original_hit_rows, original_hit_order),
        (original_m8_rows, original_m8_order),
    ) = tokio::try_join!(
        updated_fut,
        added_fut,
        blast_hits_fut,
        hit_rows_fut,
        m8_rows_fut
    )?;

    let original_hit_map: AHashMap<String, Bytes> = original_hit_rows.into_iter().collect();
    let original_m8_map: AHashMap<String, Bytes> = original_m8_rows.into_iter().collect();

    info!(
        "[generate_m8_and_hit_summary] collected: updated={} added={} blast_hits={} orig_hit={} orig_m8={}",
        updated_reads.len(),
        added_reads.len(),
        blast_hits.len(),
        original_hit_map.len(),
        original_m8_map.len()
    );

    let original_hit_ids: HashSet<String> = original_hit_order.iter().cloned().collect();
    let original_m8_ids: HashSet<String> = original_m8_order.iter().cloned().collect();

    let thread_pool = config.thread_pool.clone();

    let (rewritten_m8, rewritten_hit_summary) = thread_pool.install(|| {
        rayon::join(
            || {
                original_m8_order
                    .par_iter()
                    .filter_map(|read_id| {
                        blast_hits
                            .get(read_id)
                            .cloned()
                            .or_else(|| original_m8_map.get(read_id).cloned())
                    })
                    .collect::<Vec<Bytes>>()
            },
            || {
                original_hit_order
                    .par_iter()
                    .filter_map(|read_id| {
                        updated_reads
                            .get(read_id)
                            .or_else(|| added_reads.get(read_id))
                            .cloned()
                            .or_else(|| original_hit_map.get(read_id).cloned())
                    })
                    .collect::<Vec<Bytes>>()
            },
        )
    });

    info!(
        "[generate_m8_and_hit_summary] rewrite done: m8={} hit={}",
        rewritten_m8.len(),
        rewritten_hit_summary.len()
    );

    // m8 append only for added_reads (new_read_ids)
    let m8_tx = refined_m8_tx;
    let blast_hits_m8 = blast_hits.clone();
    let added_order_m8 = added_order.clone();
    let original_m8_ids_m8 = original_m8_ids;
    let m8_send = async move {
        let mut emitted = 0usize;
        for row in rewritten_m8 {
            m8_tx
                .send(ParseOutput::Bytes(row))
                .await
                .map_err(|_| anyhow!("refined_m8_tx dropped"))?;
            emitted += 1;
            if emitted % 100_000 == 0 {
                info!("[generate_m8_and_hit_summary] m8 emit progress {}", emitted);
            }
        }

        let mut appended = 0usize;
        for read_id in &added_order_m8 {
            if !original_m8_ids_m8.contains(read_id) {
                if let Some(row) = blast_hits_m8.get(read_id) {
                    m8_tx
                        .send(ParseOutput::Bytes(row.clone()))
                        .await
                        .map_err(|_| anyhow!("refined_m8_tx dropped"))?;
                    appended += 1;
                }
            }
        }

        info!(
            "[generate_m8_and_hit_summary] m8 emit done: rewritten_path={} appended={}",
            emitted, appended
        );
        drop(m8_tx);
        Ok::<(), anyhow::Error>(())
    };

    // hit append only for added_reads
    let hit_tx = refined_hit_summary_tx;
    let added_reads_hit = added_reads;
    let added_order_hit = added_order;
    let hit_send = async move {
        let mut emitted = 0usize;
        for row in rewritten_hit_summary {
            hit_tx
                .send(ParseOutput::Bytes(row))
                .await
                .map_err(|_| anyhow!("refined_hit_summary_tx dropped"))?;
            emitted += 1;
            if emitted % 100_000 == 0 {
                info!(
                    "[generate_m8_and_hit_summary] hit emit progress {}",
                    emitted
                );
            }
        }

        let mut appended = 0usize;
        for read_id in &added_order_hit {
            if !original_hit_ids.contains(read_id) {
                if let Some(row) = added_reads_hit.get(read_id) {
                    hit_tx
                        .send(ParseOutput::Bytes(row.clone()))
                        .await
                        .map_err(|_| anyhow!("refined_hit_summary_tx dropped"))?;
                    appended += 1;
                }
            }
        }

        info!(
            "[generate_m8_and_hit_summary] hit emit done: rewritten_path={} appended={}",
            emitted, appended
        );
        drop(hit_tx);
        Ok::<(), anyhow::Error>(())
    };

    tokio::try_join!(m8_send, hit_send)?;
    info!("[generate_m8_and_hit_summary] complete");
    Ok(())
}

/// Handles a graceful early exit for BLAST contigs when no usable results are found.
///
/// # Arguments
///
/// * `db_type`: database type (e.g., "nt", "nr")
/// * `blast_m8`: path for the BLAST m8 file
/// * `blast_top_m8`: path for the BLAST top m8 file
/// * `refined_m8`: path for the refined m8 file
/// * `refined_hit_summary`: path for the refined hit summary file
/// * `refined_counts`: path for the refined counts file
/// * `contig_summary`: path for the contig summary file
/// * `deduped_m8`: optional path to the input deduped m8 file
/// * `hit_summary`: optional path to the input hit summary file
/// * `orig_counts`: optional path to the input original counts file
/// * `read_dict`: shared map of read IDs to their taxonomic hits
/// * `taxon_counts`: existing taxon counts
/// * `fn_start`: start time of the blast_contigs function
///
/// # Returns
///
/// Result containing empty result sets and placeholders
async fn early_blast_exit(
    db_type: &str,
    blast_m8: &PathBuf,
    blast_top_m8: &PathBuf,
    refined_m8: &PathBuf,
    refined_hit_summary: &PathBuf,
    refined_counts: &PathBuf,
    contig_summary: &PathBuf,
    deduped_m8: Option<&PathBuf>,
    hit_summary: Option<&PathBuf>,
    orig_counts: Option<&PathBuf>,
    read_dict: &Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    taxon_counts: Vec<TaxonCount>,
    fn_start: tokio::time::Instant,
) -> Result<(
    AHashMap<String, Arc<ReadHit>>,
    Vec<TaxonCount>,
    Vec<ContigSummaryEntry>,
    mpsc::Receiver<ReducedRead>,
    mpsc::UnboundedReceiver<ParseOutput>,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
)> {
    warn!(
        "[blast_contigs:{}] graceful early exit (no usable BLAST results)",
        db_type
    );

    if let (Some(d), Some(h), Some(c)) = (deduped_m8, hit_summary, orig_counts) {
        // Preserve previous read-level results (preferred path)
        tokio::fs::write(blast_m8, b" ").await?;
        tokio::fs::write(blast_top_m8, b" ").await?;
        tokio::fs::copy(d, refined_m8).await?;
        tokio::fs::copy(h, refined_hit_summary).await?;
        tokio::fs::copy(c, refined_counts).await?;
        tokio::fs::write(contig_summary, b"[]").await?;
    } else {
        // Fallback when we don't have the source files (current streaming situation)
        tokio::fs::write(blast_m8, b" ").await?;
        tokio::fs::write(blast_top_m8, b" ").await?;
        tokio::fs::write(refined_m8, b" ").await?;
        tokio::fs::write(refined_hit_summary, b" ").await?;
        tokio::fs::write(refined_counts, b" ").await?;
        tokio::fs::write(contig_summary, b"[]").await?;
    }

    let final_read_dict = read_dict.lock().unwrap().clone();
    let (_pairs_tx, pairs_rx) = mpsc::channel::<ReducedRead>(1);
    let (_top_tx, top_rx) = mpsc::unbounded_channel::<ParseOutput>();

    info!(
        "[blast_contigs:{}] empty-output fast path complete after {:?}",
        db_type,
        fn_start.elapsed()
    );


    Ok((
        final_read_dict,
        taxon_counts,
        vec![],
        pairs_rx,
        top_rx,
        vec![],
        vec![],
    ))
}

/// Refines taxonomic assignments of reads by aligning their parent contigs against NT/NR.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `db_type`: database type ("nt" or "nr")
/// * `deduped_m8_stream`: stream of original m8 hits for reads
/// * `hit_summary_stream`: stream of original hit summary records
/// * `read_dict`: shared map of read IDs to their taxonomic hits
/// * `accession_map`: map of accessions to hit metadata
/// * `taxon_counts`: existing taxon counts for reads
/// * `assembled_contig_fasta`: path to the assembled contigs FASTA file
/// * `read2contig`: map from read ID to contig ID
/// * `reference_fasta`: path to the reference sequences for alignment
/// * `duplicate_clusters`: map of cluster information for duplicate weighting
/// * `lineage_map`: map of taxids to lineages
/// * `should_keep_filter`: filter for taxids to keep
/// * `blast_headroom`: memory headroom for BLAST/Diamond
/// * `concurrency`: number of concurrent workers
///
/// # Returns
///
/// Result containing refined hits and summary data
async fn blast_contigs(
    config: Arc<RunConfig>,
    db_type: &'static str,
    sample_base: &String,
    mut pairs: ReceiverStream<ReducedRead>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    accession_map: Arc<AHashMap<String, AccessionHit>>,
    taxon_counts: Vec<TaxonCount>,
    assembled_contig_fasta: &PathBuf,
    read2contig: Arc<HashMap<String, String>>,
    reference_fasta: &PathBuf,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    blast_headroom: u64,
    concurrency: usize,
) -> Result<(
    AHashMap<String, Arc<ReadHit>>,
    Vec<TaxonCount>,
    Vec<ContigSummaryEntry>,
    mpsc::Receiver<ReducedRead>,          // refined pairs (m8 + hit summary)
    mpsc::UnboundedReceiver<ParseOutput>, // top contig m8 for coverage viz
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
    Option<TempDir>,
)> {
    fn blastdb_ready(prefix: &PathBuf, db_type: &str) -> bool {
        let (a, b, c, alias) = if db_type == NT_TAG {
            ("nhr", "nin", "nsq", "nal")
        } else {
            ("phr", "pin", "psq", "pal")
        };
        let p = |ext: &str| prefix.with_extension(ext).exists();
        p(alias) || (p(a) && p(b) && p(c))
    }

    fn read_id_from_bytes(b: &[u8]) -> String {
        let s = String::from_utf8_lossy(b);
        s.split('\t').next().unwrap_or("").trim().to_string()
    }

    use tokio::time::Instant;

    let fn_start = Instant::now();
    info!(
        "[blast_contigs:{}] START — concurrency={}, blast_headroom={}, out_dir={}",
        db_type,
        concurrency,
        blast_headroom,
        config.out_dir.display()
    );

    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    let out_dir = config.out_dir.join("assembly");
    tokio::fs::create_dir_all(&out_dir).await?;

    let (
        blast_m8_name,
        blast_top_m8_name,
        refined_m8_name,
        refined_hit_summary_name,
        refined_counts_name,
        contig_summary_name,
    ) = match db_type {
        "nt" => (
            "gsnap.blast.m8",
            "gsnap.blast.top.m8",
            "gsnap.reassigned.m8",
            "gsnap.hitsummary2.tab",
            "refined_gsnap_counts_with_dcr.json",
            "gsnap_contig_summary.json",
        ),
        "nr" => (
            "rapsearch2.blast.m8",
            "rapsearch2.blast.top.m8",
            "rapsearch2.reassigned.m8",
            "rapsearch2.hitsummary2.tab",
            "refined_rapsearch2_counts_with_dcr.json",
            "rapsearch2_contig_summary.json",
        ),
        other => {
            return Err(anyhow!(
        "blast_contigs: unsupported db_type '{}'",
        other
        ));
        }
    };

    let blast_m8_path = out_dir.join(blast_m8_name);
    let blast_top_m8_path = out_dir.join(blast_top_m8_name);
    let refined_m8_path = out_dir.join(refined_m8_name);
    let refined_hit_summary_path =
        out_dir.join(refined_hit_summary_name);
    let refined_counts_path = out_dir.join(refined_counts_name);
    let contig_summary_path = out_dir.join(contig_summary_name);

    let contig_size = file_size(assembled_contig_fasta).await?;
    let ref_size = file_size(reference_fasta).await?;
    info!(
        "[blast_contigs:{}] input sizes — contig={} ref={}",
        db_type, contig_size, ref_size
    );

    if contig_size < MIN_ASSEMBLED_CONTIG_SIZE || ref_size < MIN_REF_FASTA_SIZE {
        while pairs.next().await.is_some() {}

        let (
            dict,
            counts,
            contigs,
            pairs_rx,
            top_rx,
            tasks,
            recvs,
        ) = early_blast_exit(
            db_type,
            &blast_m8_path,
            &blast_top_m8_path,
            &refined_m8_path,
            &refined_hit_summary_path,
            &refined_counts_path,
            &contig_summary_path,
            None,
            None,
            None,
            &read_dict,
            taxon_counts,
            fn_start,
        )
            .await?;

        // No makeblastdb → no index TempDir
        return Ok((
            dict, counts, contigs, pairs_rx, top_rx, tasks, recvs, None,
        ));
    }

    // Split ReducedRead → original m8 + hit for generate_m8_and_hit_summary
    // (that function keys by read_id; pairing is preserved in the maps)
    let (orig_m8_tx, orig_m8_rx) = mpsc::channel::<ParseOutput>(131_072);
    let (orig_hit_tx, orig_hit_rx) = mpsc::channel::<ParseOutput>(131_072);
    let pairs_split_handle = tokio::spawn({
        let db_type = db_type;
        async move {
            let mut n = 0u64;
            while let Some(pair) = pairs.next().await {
                if pair.dedup.is_empty() && pair.summary.is_empty() {
                    continue;
                }
                if orig_m8_tx
                    .send(ParseOutput::Bytes(Bytes::from(pair.dedup)))
                    .await
                    .is_err()
                {
                    return Err(anyhow!(
                        "[blast_contigs:{}] orig_m8 consumer closed after {} pairs",
                        db_type,
                        n
                    ));
                }
                if orig_hit_tx
                    .send(ParseOutput::Bytes(Bytes::from(pair.summary)))
                    .await
                    .is_err()
                {
                    return Err(anyhow!(
                        "[blast_contigs:{}] orig_hit consumer closed after {} pairs",
                        db_type,
                        n
                    ));
                }
                n += 1;
            }
            info!(
                "[blast_contigs:{}] pairs split complete: {} pairs",
                db_type, n
            );
            Ok::<(), anyhow::Error>(())
        }
    });
    let temp_dir = choose_temp_dir(
        ref_size + contig_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        blast_headroom,
        true,
    )
        .await?;

    let blastdb_path = temp_dir.path().join(format!("{}_blastindex", db_type));

    let makeblastdb_config = MakeblastdbConfig::new()
        .option(
            "-in",
            reference_fasta.to_string_lossy().to_string(),
        )
        .option(
            "-dbtype",
            if db_type == NT_TAG {
                "nucl"
            } else {
                "prot"
            },
        )
        .option(
            "-out",
            blastdb_path.to_string_lossy().to_string(),
        );

    let makeblastdb_args =
        generate_cli(MAKEBLASTDB_TAG, &config, Some(&makeblastdb_config))
            .map_err(|e| PipelineError::ToolExecution {
                tool: MAKEBLASTDB_TAG.to_string(),
                error: e.to_string(),
            })?;

    info!("[blast_contigs:{}] launching makeblastdb", db_type);

    let (mut makeblastdb_child, makeblastdb_err_task) = spawn_cmd(
        config.clone(),
        MAKEBLASTDB_TAG,
        makeblastdb_args,
        config.args.verbose,
        None,
    )
        .await?;

    let makeblastdb_status = makeblastdb_child
        .wait()
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAKEBLASTDB_TAG.to_string(),
            error: format!("makeblastdb wait failed: {}", e),
        })?;

    makeblastdb_err_task.await??;

    if !makeblastdb_status.success() {
        return Err(anyhow!(
        "makeblastdb exited with {:?} (db_type={}, prefix={})",
        makeblastdb_status.code(),
        db_type,
        blastdb_path.display()
    )
            .into());
    }

    if !blastdb_ready(&blastdb_path, db_type) {
        return Err(anyhow!(
        "makeblastdb succeeded but index not ready at prefix={}",
        blastdb_path.display()
    )
            .into());
    }

    info!(
    "[blast_contigs:{}] blastdb ready at prefix={}",
    db_type,
    blastdb_path.display()
);

    let blast_command = if db_type == NT_TAG {
        BLASTN_TAG
    } else {
        BLASTX_TAG
    };

    let blast_args = if db_type == NT_TAG {
        let blastn_config = BlastnConfig::new()
            .option(
                "-query",
                assembled_contig_fasta.to_string_lossy().to_string(),
            )
            .option(
                "-db",
                blastdb_path.to_string_lossy().to_string(),
            )
            .option(
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
            )
            .option("-evalue", 1e-10.to_string())
            .option("-max_target_seqs", "5000")
            .option(
                "-num_threads",
                config
                    .thread_allocation(BLASTN_TAG, None)
                    .to_string(),
            );

        generate_cli(
            BLASTN_TAG,
            &config,
            Some(&blastn_config),
        )?
    } else {
        let blastx_config = BlastxConfig::new()
            .option(
                "-query",
                assembled_contig_fasta.to_string_lossy().to_string(),
            )
            .option(
                "-db",
                blastdb_path.to_string_lossy().to_string(),
            )
            .option("-outfmt", "6")
            .option("-evalue", 1e-10.to_string())
            .option("-num_alignments", "5")
            .option(
                "-num_threads",
                config
                    .thread_allocation(BLASTX_TAG, None)
                    .to_string(),
            );

        generate_cli(
            BLASTX_TAG,
            &config,
            Some(&blastx_config),
        )?
    };

    info!(
    "[blast_contigs:{}] launching {} query={} db={}",
    db_type,
    blast_command,
    assembled_contig_fasta.display(),
    blastdb_path.display()
);

    let (mut blast_child, err_task) = spawn_cmd(
        config.clone(),
        blast_command,
        blast_args,
        config.args.verbose,
        None,
    )
        .await?;
    cleanup_tasks.push(err_task);

    let blast_out_stream = parse_child_output(
        &mut blast_child,
        ChildStream::Stdout,
        ParseMode::Lines,
        &config,
    )
        .await?;

    // ── Top hits ──────────────────────────────────────────────────────────
    let (top_internal_tx, top_internal_rx) = mpsc::channel(131_072);
    let (top_for_update_tx, top_for_update_rx) = mpsc::channel(131_072);
    let (top_for_caller_tx, top_for_caller_rx) = mpsc::unbounded_channel();

    let top_handle = if db_type == NT_TAG {
        let c = compute_phase_concurrency(&config, "top_m8_nt", 0.25, 1.0, 128, 2);
        let bs = compute_batch_size(
            None,
            200,
            if config.available_ram < 64 { 100 } else { 250 },
            c,
        );
        let t0 = Instant::now();
        tokio::spawn(async move {
            info!("[blast_contigs:{}] get_top_m8_nt started", db_type);
            let res =
                get_top_m8_nt(ReceiverStream::new(blast_out_stream), top_internal_tx, c, bs).await;
            info!(
                "[blast_contigs:{}] get_top_m8_nt finished after {:?}",
                db_type,
                t0.elapsed()
            );
            res
        })
    } else {
        let c = compute_phase_concurrency(&config, "top_m8_nr", 0.5, 1.5, 128, 2);
        let bs = compute_batch_size(
            None,
            200,
            if config.available_ram < 64 { 100 } else { 250 },
            c,
        );
        let t0 = Instant::now();
        tokio::spawn(async move {
            info!("[blast_contigs:{}] get_top_m8_nr started", db_type);
            let res =
                get_top_m8_nr(ReceiverStream::new(blast_out_stream), top_internal_tx, c, bs).await;
            info!(
                "[blast_contigs:{}] get_top_m8_nr finished after {:?}",
                db_type,
                t0.elapsed()
            );
            res
        })
    };

    let forward_handle = tokio::spawn({
        let mut rx = top_internal_rx;
        let tx_update = top_for_update_tx;
        let tx_caller = top_for_caller_tx;
        async move {
            let mut total = 0usize;
            while let Some(item) = rx.recv().await {
                total += 1;
                if tx_update.send(item.clone()).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] top forward: update dropped at {}",
                        db_type, total
                    );
                    break;
                }
                if tx_caller.send(item).is_err() {
                    warn!(
                        "[blast_contigs:{}] top forward: caller dropped at {}",
                        db_type, total
                    );
                    break;
                }
            }
            info!(
                "[blast_contigs:{}] top forward complete: total={}",
                db_type, total
            );
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(forward_handle);

    let (contig2lineage_tx, contig2lineage_rx) = mpsc::channel(64_000);
    let (read2blastm8_tx, read2blastm8_rx) = mpsc::channel(64_000);
    let (updated_tx, updated_rx) = mpsc::channel(64_000);
    let (added_tx, added_rx) = mpsc::channel(64_000);

    let update_handle = tokio::spawn({
        let t0 = Instant::now();
        let fut = update_read_dict(
            config.clone(),
            read2contig.clone(),
            ReceiverStream::new(top_for_update_rx),
            read_dict.clone(),
            lineage_map.clone(),
            accession_map.clone(),
            should_keep_filter.clone(),
            db_type,
            contig2lineage_tx,
            read2blastm8_tx,
            updated_tx,
            added_tx,
        );
        async move {
            let res = fut.await;
            info!(
                "[blast_contigs:{}] update_read_dict finished after {:?}",
                db_type,
                t0.elapsed()
            );
            res
        }
    });

    let contig2lineage_task = tokio::spawn(async move {
        let mut contig2lineage: AHashMap<String, [i32; 3]> = AHashMap::new();
        let mut stream = ReceiverStream::new(contig2lineage_rx);
        while let Some(item) = stream.next().await {
            let bytes = item.to_bytes()?;
            let line = String::from_utf8_lossy(&bytes);
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                contig2lineage.insert(
                    fields[0].to_string(),
                    [
                        fields[1].parse().unwrap_or(0),
                        fields[2].parse().unwrap_or(0),
                        fields[3].parse().unwrap_or(0),
                    ],
                );
            }
        }
        Ok::<_, anyhow::Error>(contig2lineage)
    });

    // Internal dual streams for generate_m8 + counts (paired again before export)
    let (refined_m8_local_tx, refined_m8_local_rx) = mpsc::channel(2_000_000);
    let (refined_hit_local_tx, refined_hit_local_rx) = mpsc::channel(2_000_000);

    let generate_m8_handle = tokio::spawn({
        let t0 = Instant::now();
        let fut = generate_m8_and_hit_summary(
            config.clone(),
            ReceiverStream::new(updated_rx),
            ReceiverStream::new(added_rx),
            ReceiverStream::new(read2blastm8_rx),
            ReceiverStream::new(orig_hit_rx),
            ReceiverStream::new(orig_m8_rx),
            refined_m8_local_tx,
            refined_hit_local_tx,
        );
        async move {
            let res = fut.await;
            info!(
                "[blast_contigs:{}] generate_m8_and_hit_summary finished after {:?}",
                db_type,
                t0.elapsed()
            );
            res
        }
    });

    let (refined_m8_for_counts_tx, refined_m8_for_counts_rx) = mpsc::channel(2_000_000);
    let (refined_hit_for_counts_tx, refined_hit_for_counts_rx) = mpsc::channel(2_000_000);

    let (m8_acc_tx, mut m8_acc_rx) = mpsc::unbounded_channel::<ParseOutput>();
    let (hit_acc_tx, mut hit_acc_rx) = mpsc::unbounded_channel::<ParseOutput>();

    let m8_forward_handle = tokio::spawn({
        let mut rx = refined_m8_local_rx;
        let counts_tx = refined_m8_for_counts_tx;
        async move {
            let mut total = 0usize;
            while let Some(item) = rx.recv().await {
                total += 1;
                if counts_tx.send(item.clone()).await.is_err() {
                    warn!("[blast_contigs] m8 counts branch closed at {}", total);
                    break;
                }
                if m8_acc_tx.send(item).is_err() {
                    break;
                }
            }
            info!("[blast_contigs] m8 forward complete: total={}", total);
            Ok::<(), anyhow::Error>(())
        }
    });

    let hit_forward_handle = tokio::spawn({
        let mut rx = refined_hit_local_rx;
        let counts_tx = refined_hit_for_counts_tx;
        async move {
            let mut total = 0usize;
            while let Some(item) = rx.recv().await {
                total += 1;
                if counts_tx.send(item.clone()).await.is_err() {
                    warn!("[blast_contigs] hit counts branch closed at {}", total);
                    break;
                }
                if hit_acc_tx.send(item).is_err() {
                    break;
                }
            }
            info!("[blast_contigs] hit forward complete: total={}", total);
            Ok::<(), anyhow::Error>(())
        }
    });

    let m8_collect_handle = tokio::spawn(async move {
        let mut v = Vec::with_capacity(1_000_000);
        while let Some(item) = m8_acc_rx.recv().await {
            v.push(item);
        }
        v
    });
    let hit_collect_handle = tokio::spawn(async move {
        let mut v = Vec::with_capacity(1_000_000);
        while let Some(item) = hit_acc_rx.recv().await {
            v.push(item);
        }
        v
    });

    let taxon_count_concurrency =
        compute_phase_concurrency(&config, "taxon_counting", 0.4, 4.0, 64, 4);
    let taxon_count_batch_size = compute_batch_size(None, 220, 150, taxon_count_concurrency);
    let (refined_counts_tx, refined_counts_rx) = mpsc::unbounded_channel();

    let counts_handle = tokio::spawn({
        let t0 = Instant::now();
        let fut = generate_taxon_count_json_from_m8(
            ReceiverStream::new(refined_m8_for_counts_rx),
            ReceiverStream::new(refined_hit_for_counts_rx),
            db_type,
            lineage_map.clone(),
            should_keep_filter.clone(),
            duplicate_clusters.clone(),
            refined_counts_tx,
            taxon_count_concurrency,
            taxon_count_batch_size,
        );
        async move {
            let res = fut.await;
            info!(
                "[blast_contigs:{}] generate_taxon_count_json_from_m8 finished after {:?}",
                db_type,
                t0.elapsed()
            );
            res
        }
    });

    // ── Ordered joins ─────────────────────────────────────────────────────
    top_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("top_handle panicked: {e}")))?
        .map_err(|e| PipelineError::Other(anyhow!("get_top_m8 failed: {e}")))?;
    info!(
        "[blast_contigs:{}] top stage joined after {:?}",
        db_type,
        fn_start.elapsed()
    );

    update_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("update_handle panicked: {e}")))?
        .map_err(|e| PipelineError::Other(anyhow!("update_read_dict failed: {e}")))?;
    info!(
        "[blast_contigs:{}] update_read_dict joined after {:?}",
        db_type,
        fn_start.elapsed()
    );

    let contig2lineage = contig2lineage_task
        .await
        .map_err(|e| PipelineError::Other(anyhow!("contig2lineage panicked: {e}")))??;
    info!(
        "[blast_contigs:{}] contig2lineage map: {} entries",
        db_type,
        contig2lineage.len()
    );

    let (contig_summary_tx, contig_summary_rx) = mpsc::channel(64_000);
    let contig_summary_handle = tokio::spawn({
        let t0 = Instant::now();
        let fut = generate_contig_summary_json(
            read2contig.clone(),
            contig2lineage,
            read_dict.clone(),
            db_type,
            duplicate_clusters.clone(),
            should_keep_filter.clone(),
            MIN_CONTIG_SIZE,
            contig_summary_tx,
        );
        async move {
            let res = fut.await;
            info!(
                "[blast_contigs:{}] generate_contig_summary_json finished after {:?}",
                db_type,
                t0.elapsed()
            );
            res
        }
    });

    let (gen_res, counts_res, contig_res, m8_fwd_res, hit_fwd_res, split_res) = tokio::try_join!(
        generate_m8_handle,
        counts_handle,
        contig_summary_handle,
        m8_forward_handle,
        hit_forward_handle,
        pairs_split_handle,
    )
        .map_err(|e| PipelineError::Other(anyhow!("contig task panicked: {e}")))?;

    gen_res.map_err(|e| PipelineError::Other(anyhow!("generate_m8: {e}")))?;
    counts_res.map_err(|e| PipelineError::Other(anyhow!("counts: {e}")))?;
    contig_res.map_err(|e| PipelineError::Other(anyhow!("contig_summary: {e}")))?;
    m8_fwd_res.map_err(|e| PipelineError::Other(anyhow!("m8 forward: {e}")))?;
    hit_fwd_res.map_err(|e| PipelineError::Other(anyhow!("hit forward: {e}")))?;
    split_res.map_err(|e| PipelineError::Other(anyhow!("pairs split: {e}")))?;

    let refined_m8_items = m8_collect_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("m8 collect panicked: {e}")))?;
    let refined_hit_items = hit_collect_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("hit collect panicked: {e}")))?;

    info!(
        "[blast_contigs:{}] refined buffered: m8={} hit={}",
        db_type,
        refined_m8_items.len(),
        refined_hit_items.len()
    );

    // ── Pair by read_id → ReducedRead export (never by index) ─────────────
    let mut m8_by_id: AHashMap<String, Bytes> =
        AHashMap::with_capacity(refined_m8_items.len());
    for item in refined_m8_items {
        let ParseOutput::Bytes(b) = item else { continue };
        let id = read_id_from_bytes(&b);
        if !id.is_empty() {
            m8_by_id.insert(id, b);
        }
    }

    let mut pairs_out: Vec<ReducedRead> = Vec::with_capacity(refined_hit_items.len());
    let mut seq = 0u64;
    let mut orphan_hits = 0u64;
    for item in refined_hit_items {
        let ParseOutput::Bytes(summary) = item else { continue };
        let id = read_id_from_bytes(&summary);
        let Some(dedup) = m8_by_id.remove(&id) else {
            orphan_hits += 1;
            continue;
        };
        let accession = {
            let s = String::from_utf8_lossy(&summary);
            s.split('\t').nth(3).unwrap_or("").trim().to_string()
        };
        seq += 1;
        pairs_out.push(ReducedRead {
            seq,
            dedup: dedup.to_vec(),
            summary: summary.to_vec(),
            accession,
        });
    }
    if !m8_by_id.is_empty() {
        warn!(
            "[blast_contigs:{}] {} refined m8 rows with no matching hit summary",
            db_type,
            m8_by_id.len()
        );
    }
    if orphan_hits > 0 {
        warn!(
            "[blast_contigs:{}] {} refined hit rows with no matching m8",
            db_type, orphan_hits
        );
    }
    info!(
        "[blast_contigs:{}] refined pairs: {}",
        db_type,
        pairs_out.len()
    );


    // Durable hitsummary2 + reassigned m8 (CZID: gsnap.hitsummary2.tab / gsnap.reassigned.m8)
    // Written from the same pairs_out that feeds merge / taxid / coverage.
    let hitsummary2_path = config
        .out_dir
        .join(format!("{}_{}.hitsummary2.tab", sample_base, db_type));
    let reassigned_m8_path = config
        .out_dir
        .join(format!("{}_{}.reassigned.m8", sample_base, db_type));

    {
        use tokio::io::AsyncWriteExt;
        let mut hit_w =
            tokio::io::BufWriter::new(tokio::fs::File::create(&hitsummary2_path).await?);
        let mut m8_w =
            tokio::io::BufWriter::new(tokio::fs::File::create(&reassigned_m8_path).await?);

        let mut long_rows = 0u64;
        for p in &pairs_out {
            let nf = p.summary.split(|&b| b == b'\t').count();
            if nf >= 12 {
                long_rows += 1;
            }
            hit_w.write_all(&p.summary).await?;
            if !p.summary.ends_with(b"\n") {
                hit_w.write_all(b"\n").await?;
            }
            m8_w.write_all(&p.dedup).await?;
            if !p.dedup.ends_with(b"\n") {
                m8_w.write_all(b"\n").await?;
            }
        }
        hit_w.flush().await?;
        m8_w.flush().await?;

        info!(
            "[blast_contigs:{}] wrote hitsummary2={} lines={} long_rows(>=12)={} reassigned_m8={}",
            db_type,
            hitsummary2_path.display(),
            pairs_out.len(),
            long_rows,
            reassigned_m8_path.display()
        );
    }



    let (refined_pairs_tx, refined_pairs_rx) = mpsc::channel(pairs_out.len().max(1));
    tokio::spawn(async move {
        for p in pairs_out {
            if refined_pairs_tx.send(p).await.is_err() {
                break;
            }
        }
    });

    let blast_status = blast_child
        .wait()
        .await
        .map_err(|e| anyhow!("blast wait failed: {}", e))?;
    if !blast_status.success() {
        return Err(anyhow!(
            "{} exited with {:?} (db_type={})",
            blast_command,
            blast_status.code(),
            db_type
        )
            .into());
    }

    let mut refined_counts = Vec::new();
    let mut rx = UnboundedReceiverStream::new(refined_counts_rx);
    while let Some(item) = rx.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        refined_counts.push(serde_json::from_str(&line)?);
    }
    info!(
        "[blast_contigs:{}] refined_counts: {}",
        db_type,
        refined_counts.len()
    );

    let mut contig_summary: Vec<ContigSummaryEntry> = Vec::new();
    let mut rx = ReceiverStream::new(contig_summary_rx);
    while let Some(item) = rx.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        contig_summary.push(
            serde_json::from_str(line)
                .map_err(|e| anyhow!("contig_summary JSON parse failed: {e}; line={line}"))?,
        );
    }
    info!(
        "[blast_contigs:{}] contig_summary: {}",
        db_type,
        contig_summary.len()
    );

    let final_read_dict = read_dict.lock().unwrap().clone();
    info!(
        "[blast_contigs:{}] DONE after {:?} — counts={} contig_summary={} read_dict={}",
        db_type,
        fn_start.elapsed(),
        refined_counts.len(),
        contig_summary.len(),
        final_read_dict.len()
    );

    Ok((
        final_read_dict,
        refined_counts,
        contig_summary,
        refined_pairs_rx,
        top_for_caller_rx,
        cleanup_tasks,
        cleanup_receivers,
        Some(temp_dir),
    ))
}



/// Resolve the production MMseqs database path once and reuse it everywhere.
fn resolve_mmseqs_db_path(config: &RunConfig) -> Result<PathBuf, PipelineError> {
    let db = config
        .args
        .mmseqs_db
        .as_ref()
        .ok_or_else(|| PipelineError::MissingArgument("mmseqs_db required".to_string()))?;

    let db_path = PathBuf::from(db);
    if !db_path.exists() {
        return Err(PipelineError::IOError(format!(
            "MMseqs DB path does not exist: {}",
            db_path.display()
        )));
    }

    Ok(db_path)
}

/// Runs a single MMseqs2 command and monitors its progress.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `args`: command-line arguments for MMseqs2
/// * `label`: label for logging and identification
///
/// # Returns
///
/// Result<()>: success or error
async fn run_mmseqs_step(
    config: Arc<RunConfig>,
    args: Vec<String>,
    label: &'static str,
) -> Result<(), PipelineError> {
    let start = Instant::now();

    let mmseqs_stderr_log = config
        .out_dir
        .join("mmseqs")
        .join(format!("{}.stderr.log", label));

    let (mut child, stderr_task) = spawn_cmd(
        config.clone(),
        MMSEQS_TAG,
        args,
        config.args.verbose,
        Some(mmseqs_stderr_log),
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MMSEQS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let wait_task = tokio::spawn(async move {
        let status = child.wait().await?;
        if !status.success() {
            return Err(anyhow!(
                "mmseqs {} exited with code {:?}",
                label,
                status.code()
            ));
        }
        Ok::<(), anyhow::Error>(())
    });

    wait_task.await.map_err(|e| {
        PipelineError::Other(anyhow!("mmseqs {} wait task panicked: {}", label, e))
    })??;

    stderr_task.await.map_err(|e| {
        PipelineError::Other(anyhow!("mmseqs {} stderr task panicked: {}", label, e))
    })??;

    let elapsed = start.elapsed();
    info!("[mmseqs:{}] completed in {:?}", label, elapsed);

    Ok(())
}

/// Starts the MMseqs2 GPU server for a specific database.
///
/// # Arguments
///
/// * `db`: path to the MMseqs2 database
/// * `verbose`: whether to inherit stderr for verbose logging
///
/// # Returns
///
/// Result containing the child process handle for the GPU server
async fn start_gpuserver(db: &PathBuf, verbose: bool) -> Result<Child, PipelineError> {
    let mut cmd = Command::new(MMSEQS_TAG);
    cmd.arg("gpuserver")
        .arg(db)
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(if verbose {
            std::process::Stdio::inherit()
        } else {
            std::process::Stdio::null()
        });

    let child = cmd.spawn().map_err(|e| PipelineError::ToolExecution {
        tool: "mmseqs gpuserver".to_string(),
        error: e.to_string(),
    })?;

    info!("Started MMseqs GPU server for {}", db.display());
    Ok(child)
}

/// Stops the MMseqs2 GPU server.
///
/// # Arguments
///
/// * `child`: handle to the GPU server child process
///
/// # Returns
///
/// Result<()>: success or error
async fn stop_gpuserver(child: &mut Child) -> Result<(), PipelineError> {
    match child.try_wait() {
        Ok(Some(_)) => return Ok(()),
        Ok(None) => {}
        Err(e) => {
            return Err(PipelineError::Other(anyhow!(
                "Failed to query gpuserver state: {}",
                e
            )));
        }
    }

    child
        .kill()
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "mmseqs gpuserver".to_string(),
            error: e.to_string(),
        })?;

    let _ = child.wait().await;
    Ok(())
}

/// Runs the MMseqs2 workflow to align a FASTQ file against a target database.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `input_fastq`: path to the input FASTQ file
/// * `temp_dir`: temporary directory for MMseqs2 work
/// * `label`: label for logging
/// * `backend`: the alignment backend to use (CPU or GPU)
///
/// # Returns
///
/// Result containing the path to the output m8 file
async fn mmseqs_fastq_to_m8_file(
    config: Arc<RunConfig>,
    input_fastq: PathBuf,
    temp_dir: &TempDir,
    label: &str,
    backend: MmseqsBackend,
) -> Result<PathBuf, PipelineError> {
    let tmp_dir = temp_dir.path().join(format!("{}.tmp", label));
    let query_db = temp_dir.path().join(format!("{}.queryDB", label));
    let result_db = temp_dir.path().join(format!("{}.resultDB", label));
    let m8_path = temp_dir.path().join(format!("{}.m8", label));

    fs::create_dir_all(&tmp_dir)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let target_db = resolve_mmseqs_db_path(&config)?;

    // 1) FASTQ -> query DB
    let createdb_cfg = MmseqsConfig {
        subcommand: MmseqsSubcommand::Createdb,
        backend: MmseqsBackend::Cpu,
        input: Some(input_fastq),
        target_db: None,
        result_db: None,
        output: Some(query_db.clone()),
        tmp_dir: None,
        threads: None,
        sensitivity: None,
        search_type: None,
        max_seqs: None,
        prefilter_mode: None,
        db_load_mode: None,
        alignment_mode: None,
        index_subset: None,
        format_output: None,
        cuda_visible_devices: None,
        option_fields: HashMap::from([
            ("--dbtype".to_string(), Some("2".to_string())), // nucleotide
        ]),
        gpu_server: false,
    };

    let createdb_args = generate_cli(MMSEQS_TAG, &config, Some(&createdb_cfg)).map_err(|e| {
        PipelineError::ToolExecution {
            tool: MMSEQS_TAG.to_string(),
            error: e.to_string(),
        }
    })?;

    info!("[mmseqs:{}] createdb args: {:?}", label, createdb_args);
    run_mmseqs_step(config.clone(), createdb_args, "createdb").await?;

    if !query_db.exists() {
        return Err(PipelineError::Other(anyhow!(
            "mmseqs createdb did not produce queryDB: {}",
            query_db.display()
        )));
    }
    let meta = fs::metadata(&query_db)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    if meta.len() < 10_000 {
        return Err(PipelineError::Other(anyhow!(
            "queryDB looks too small ({} bytes) — createdb probably failed",
            meta.len()
        )));
    }

    info!(
        "Query DB created successfully: {} ({} bytes)",
        query_db.display(),
        meta.len()
    );

    // 2) GPU server, only for GPU backend
    let mut gpu_server: Option<Child> = None;

    if backend == MmseqsBackend::Gpu {
        gpu_server = Some(start_gpuserver(&target_db, config.args.verbose).await?);

        for _ in 0..20 {
            if let Some(server) = gpu_server.as_mut() {
                if let Ok(Some(status)) = server.try_wait() {
                    return Err(PipelineError::Other(anyhow!(
                        "MMseqs GPU server exited during warmup for {} with status {:?}",
                        target_db.display(),
                        status
                    )));
                }
            }
            tokio::time::sleep(Duration::from_secs(1)).await;
        }
    }

    // 3) search -> result DB
    let search_cfg = match backend {
        MmseqsBackend::Cpu => build_mmseqs_cpu_search_config(
            query_db.clone(),
            target_db.clone(),
            result_db.clone(),
            tmp_dir.clone(),
            config.mmseqs_threads(),
        ),

        MmseqsBackend::Gpu => build_mmseqs_gpu_search_config(
            query_db.clone(),
            target_db.clone(),
            result_db.clone(),
            tmp_dir.clone(),
            config.mmseqs_threads(),
        ),
    };

    // 3) search -> result DB
    let search_args = generate_cli(MMSEQS_TAG, &config, Some(&search_cfg))?;
    info!("[mmseqs:{}] search args: {:?}", label, search_args);

    // NO manual prepend here — spawn_cmd handles NUMA automatically
    let search_res = run_mmseqs_step(config.clone(), search_args, "search").await;

    if let Some(mut server) = gpu_server {
        info!("Stopping MMseqs GPU server...");
        stop_gpuserver(&mut server).await?;
    }

    search_res?;

    // 4) convertalis -> m8
    let convert_cfg = MmseqsConfig {
        subcommand: MmseqsSubcommand::ConvertAlis,
        backend: MmseqsBackend::Cpu,
        input: Some(query_db.clone()),
        target_db: Some(target_db.clone()),
        result_db: Some(result_db.clone()),
        output: Some(m8_path.clone()),
        tmp_dir: None,
        threads: None,
        sensitivity: None,
        search_type: None,
        max_seqs: None,
        prefilter_mode: None,
        db_load_mode: None,
        alignment_mode: None,
        index_subset: None,
        format_output: Some(
            "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
                .to_string(),
        ),
        cuda_visible_devices: None,
        option_fields: HashMap::new(),
        gpu_server: false,
    };

    let convert_args = generate_cli(MMSEQS_TAG, &config, Some(&convert_cfg)).map_err(|e| {
        PipelineError::ToolExecution {
            tool: MMSEQS_TAG.to_string(),
            error: e.to_string(),
        }
    })?;

    info!("[mmseqs:{}] convertalis args: {:?}", label, convert_args);

    run_mmseqs_step(config.clone(), convert_args, "convertalis").await?;

    let meta = fs::metadata(&m8_path)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    if meta.len() == 0 {
        return Err(PipelineError::Other(anyhow!(
            "mmseqs produced an empty m8 file: {}",
            m8_path.display()
        )));
    }

    Ok(m8_path)
}

/// Aligns unmapped reads against NR database using MMseqs2.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `r1_path`: path to the R1 unmapped FASTQ file
/// * `r2_path_opt`: optional path to the R2 unmapped FASTQ file
/// * `backend`: the alignment backend to use (CPU or GPU)
///
/// # Returns
///
/// Result containing:
/// - receiver for alignment results (m8 format)
/// - vector of cleanup tasks
/// - vector of cleanup receivers
/// - vector of temporary directories used
async fn mmseqs_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    backend: MmseqsBackend,
) -> Result<
    (
        tokio::sync::mpsc::Receiver<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<tokio::sync::oneshot::Receiver<Result<(), anyhow::Error>>>,
        Vec<TempDir>,
    ),
    PipelineError,
> {
    let cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    info!("Running mmseqs non-host align (backend: {:?})", backend);

    let temp_dir = choose_temp_dir(
        100_000_000_000,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        true,
    )
        .await?;

    let combined_fastq = temp_dir.path().join("mmseqs_non_host_combined.fastq");

    write_combined_fastq(r1_path, r2_path_opt.clone(), &combined_fastq).await?;

    let merged_m8 = mmseqs_fastq_to_m8_file(
        config.clone(),
        combined_fastq,
        &temp_dir,
        "nr_merged",
        backend,
    )
        .await?;

    // Hard fail before streaming if convertalis produced nothing.
    let meta = fs::metadata(&merged_m8).await.map_err(|e| {
        PipelineError::IOError(format!(
            "mmseqs m8 metadata failed for {}: {}",
            merged_m8.display(),
            e
        ))
    })?;

    let m8_size = meta.len();
    info!(
        "[mmseqs NR] m8 file on disk: {} ({} bytes)",
        merged_m8.display(),
        m8_size
    );

    if m8_size == 0 {
        return Err(PipelineError::Other(anyhow!(
            "mmseqs produced an empty m8 file: {}",
            merged_m8.display()
        )));
    }

    let m8_file = fs::File::open(&merged_m8)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let m8_rx = parse_lines(m8_file, &config, StreamDataType::JustBytes)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok((m8_rx, cleanup_tasks, cleanup_receivers, vec![temp_dir]))
}

async fn collect_m8_accession_map(
    mut stream: ReceiverStream<ReducedRead>,
) -> Result<AHashMap<String, String>> {
    let mut map = AHashMap::with_capacity(1_000_000);
    while let Some(pair) = stream.next().await {
        let line = String::from_utf8_lossy(&pair.dedup);
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let Some(q) = fields.next() else { continue };
        let Some(s) = fields.next() else { continue };
        if !s.is_empty() && s != "-" {
            // Prefer already-parsed accession when present
            let acc = if !pair.accession.is_empty() {
                pair.accession.as_str()
            } else {
                s
            };
            map.insert(q.to_string(), acc.to_string());
        }
    }
    info!("collect_m8_accession_map → {} entries", map.len());
    Ok(map)
}

/// The main entry point for the short-read metagenomic NGS (mNGS) pipeline.
///
/// # Arguments
///
/// * `config`: the run configuration
///
/// # Returns
///
/// Result<()>: success or error
pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();

    let efs_run_dir = config.efs_runs_dir.join(&config.run_id);
    let efs_work_dir = efs_run_dir.join("work");

    if matches!(config.execution_mode, ExecutionMode::Distributed) {
        tokio::fs::create_dir_all(&efs_work_dir)
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Failed to create distributed EFS work directory {}: {}",
                efs_work_dir.display(),
                e
            ))
            })?;

        info!(
        "[run] Distributed NR EFS work directory: {}",
        efs_work_dir.display()
    );
    }

    let mut cleanup_tasks: Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>> =
        Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let mut host_filter_temp_dirs: Vec<TempDir> = Vec::new();
    let mut final_temp_dirs: Vec<TempDir> = Vec::new();
    let temp_paths: Vec<PathBuf> = Vec::new();

    info!("Starting short read mNGS pipeline.");

    // *******************
    // Setup and Validation
    // *******************

    // External tools check

    let mut versions_vec = vec![
        BOWTIE2_TAG,
        MINIMAP2_TAG,
        KALLISTO_TAG,
        SPADES_TAG,
        MAKEBLASTDB_TAG,
        BLASTN_TAG,
        BLASTX_TAG,
        SORT_TAG,
        SAMTOOLS_TAG,
        FASTP_TAG,
        HISAT2_TAG,
    ];

    match config.alignment_backend {
        NRAlignmentBackend::Diamond => {
            versions_vec.push(DIAMOND_TAG);
        }
        _ => {
            versions_vec.push(MMSEQS_TAG);
        }
    }

    if !config.args.no_version_check {
        check_versions(versions_vec, &out_dir.clone(), &config)
            .await
            .map_err(|e| PipelineError::Other(e.into()))?;
    }

    // Check required files
    let host_bowtie2_index: String = config.args.host_bowtie2_index.clone().ok_or_else(|| {
        PipelineError::MissingArgument("host_bowtie2_index is required".to_string())
    })?;

    let host_hisat2_index: String = config.args.host_hisat2_index.clone().ok_or_else(|| {
        PipelineError::MissingArgument("host_hisat2_index is required".to_string())
    })?;

    let (file1_path, file2_path, sample_base_buf, sample_base) =
        validate_file_inputs(&config, &cwd, true).await?;
    let paired = file2_path.is_some();

    debug!(
        "file1 path: {:?}  sample base buf: {:?}  sample base: {:?}",
        file1_path.display(),
        sample_base_buf.display(),
        sample_base
    );

    let taxonomy_handle = tokio::spawn(load_lineage_and_acc2tax_maps(config.clone()));

    // Input Validation
    let (
        val_out_stream,
        validate_cleanup_tasks,
        validate_cleanup_receivers,
        val_count_task,
    ) = validate_input(
        config.clone(),
        file1_path.clone(),
        file2_path.clone(),
        sample_base_buf.clone(),
        &out_dir,
    )
        .await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // *******************
    // Host Filtering Stage
    // *******************

    // ERCC bt2 filtering and count
    let ercc_bt2_index_path = bowtie2_index_prep(&config.args.ercc_bowtie2_index, &cwd)?;
    let (
        ercc_bt2_out_stream,
        ercc_count_rx,
        ercc_bt2_cleanup_tasks,
        ercc_bt2_cleanup_receivers,
        ercc_bt2_bam_write_handle,
        _ercc_bt2_bam_path,
        ercc_bt2_temp_dirs,
    ) = bowtie2_filter_stream(
        config.clone(),
        val_out_stream,
        ercc_bt2_index_path,
        paired,
        out_dir.join("ercc_bt2_aligned_sorted.bam"),
    )
        .await?;
    cleanup_tasks.extend(ercc_bt2_cleanup_tasks);
    final_temp_dirs.extend(ercc_bt2_temp_dirs);

    // QC with fastp
    let (qc_fastp_out_stream, qc_cleanup_tasks, qc_cleanup_receivers, qc_count_result_rx) =
        fastp_qc(config.clone(), paired, ercc_bt2_out_stream).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);

    // Split for Kallisto and Bowtie2
    let (kallisto_streams, kallisto_done) = fanout_to_channels(
        qc_fastp_out_stream,
        2,
        "kallisto_split",
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    let mut kallisto_streams_iter = kallisto_streams.into_iter();

    let kallisto_bypass_stream = ReceiverStream::new(
        kallisto_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let kallisto_stream = ReceiverStream::new(
        kallisto_streams_iter
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );

    let (
        kallisto_ercc_tx,
        kallisto_ercc_rx,
        kallisto_cleanup_tasks,
        kallisto_cleanup_receivers,
        kallisto_exit_task,
    ) = kallisto_quant(
        config.clone(),
        kallisto_stream,
        out_dir.clone(),
        paired,
        sample_base_buf.clone(),
    )
        .await?;
    cleanup_tasks.extend(kallisto_cleanup_tasks);
    cleanup_receivers.extend(kallisto_cleanup_receivers);

    // Host filtering: bt2
    let host_bt2_index_path = bowtie2_index_prep(host_bowtie2_index, &cwd)?;
    let (
        host_bt2_r1,
        host_bt2_r2,
        host_bt2_count_rx,
        host_bt2_cleanup_tasks,
        host_bt2_cleanup_receivers,
        host_bt2_bam_write_handle,
        host_bt2_bam_path,
        host_bt2_temp_dirs,
    ) = bowtie2_filter_files(
        config.clone(),
        kallisto_bypass_stream,
        host_bt2_index_path,
        paired,
        out_dir.join("host_bt2_aligned_sorted.bam"),
    )
        .await?;
    cleanup_tasks.extend(host_bt2_cleanup_tasks);
    host_filter_temp_dirs.extend(host_bt2_temp_dirs);

    // HISAT2 forced checkpoint cleanup block
    // ERCC bt2 fanout is also done by this point in the pipeline
    for rx in ercc_bt2_cleanup_receivers {
        rx.await??;
    }

    for rx in qc_cleanup_receivers {
        rx.await??;
    }

    // because the stream is forced to checkpoint by hisat2, we can await the kallisto receiver from the fanout here
    kallisto_done.await?;

    // Host bt2 fanout is done (its fastq branch fed into hisat2)
    for rx in host_bt2_cleanup_receivers {
        rx.await??;
    }

    //host filtering hisat2
    // Host filtering HISAT2
    let host_hisat2_index_path =
        hisat2_index_prep(host_hisat2_index, &cwd)?;

    let (
        host_hisat2_out_stream,
        host_hisat2_count_rx,
        mut host_hisat2_cleanup_tasks,
        mut host_hisat2_cleanup_receivers,
        hisat_temp_dirs,
    ) = hisat2_filter(
        config.clone(),
        host_bt2_r1,
        host_bt2_r2,
        host_hisat2_index_path,
        paired,
        Vec::new(),
        None,
        4,
    )
        .await?;
    cleanup_tasks.append(&mut host_hisat2_cleanup_tasks);
    host_filter_temp_dirs.extend(hisat_temp_dirs);

    // If host is no huma, run an additional filter stage using a human reference
    let mut optional_human_bam_write_handle: Option<JoinHandle<Result<()>>> = None;
    let mut optional_human_bam_path: Option<PathBuf> = None;

    let mut hisat2_cleanup_receivers = host_hisat2_cleanup_receivers;

    let post_filter_stream = if config.args.human_host {
        host_hisat2_out_stream
    } else {
        let human_bowtie2_index: String = config.args.human_bowtie2_index.clone();
        let human_hisat2_index: String = config.args.human_hisat2_index.clone();

        let human_bt2_index_path =
            bowtie2_index_prep(human_bowtie2_index, &cwd)?;

        let (
            human_bt2_r1,
            human_bt2_r2,
            _human_bt2_count_rx,
            human_bt2_cleanup_tasks,
            human_bt2_cleanup_receivers,
            human_bt2_bam_write_handle,
            human_bt2_bam_path,
            human_bt2_temp_dirs,
        ) = bowtie2_filter_files(
            config.clone(),
            host_hisat2_out_stream,
            human_bt2_index_path,
            paired,
            out_dir.join("human_bt2_aligned_sorted.bam"),
        )
            .await?;

        cleanup_tasks.extend(human_bt2_cleanup_tasks);
        host_filter_temp_dirs.extend(human_bt2_temp_dirs);

        for rx in human_bt2_cleanup_receivers {
            rx.await??;
        }

        // Save the BAM write handle and path for later awaiting / insert stats
        optional_human_bam_write_handle = Some(human_bt2_bam_write_handle);
        optional_human_bam_path = Some(human_bt2_bam_path);

        let human_hisat2_index_path = hisat2_index_prep(human_hisat2_index, &cwd)?;

        let (
            human_hisat2_out_stream,
            _human_hisat2_count_rx,
            mut human_hisat2_cleanup_tasks,
            mut human_hisat2_cleanup_receivers,
            human_hisat_temp_dirs,
        ) = hisat2_filter(
            config.clone(),
            human_bt2_r1,
            human_bt2_r2,
            human_hisat2_index_path,
            paired,
            Vec::new(),
            None,
            4,
        )
            .await?;

        cleanup_tasks.append(&mut human_hisat2_cleanup_tasks);
        hisat2_cleanup_receivers.extend(human_hisat2_cleanup_receivers);
        host_filter_temp_dirs.extend(human_hisat_temp_dirs);

        human_hisat2_out_stream
    };


    let pre_dedup_parsed_stream =
        if paired {
            let (stream, parse_task) =
                parse_byte_stream_to_fastq(
                    post_filter_stream.into_inner(),
                    config.base_buffer_size,
                    config.args.stall_threshold,
                )
                    .await?;

            cleanup_tasks.push(parse_task);

            stream
        } else {
            post_filter_stream.into_inner()
        };

    let (
        dedup_stream,
        dedup_count_rx,
        dedup_finalize_task,
        mut dedup_cleanup_tasks,
        mut dedup_cleanup_receivers,
    ) = dedup(
        config.clone(),
        pre_dedup_parsed_stream,
        paired,
        Some(70),
        out_dir.clone(),
    )
        .await?;

    cleanup_tasks.append(&mut dedup_cleanup_tasks);
    cleanup_receivers.append(&mut dedup_cleanup_receivers);


    // Uniform subsampling without replacement.
    // In paired mode, max_subsample and the returned count are logical pairs.
    let (
        subsampled_stream,
        subsample_count_rx,
        _subsample_base_rx,
        subsample_send_task,
    ) = subsample_uniform_reservoir(
        config.clone(),
        dedup_stream.into_inner(),
        config.args.max_subsample as u64,
        paired,
    )
        .await?;

    cleanup_tasks.push(subsample_send_task);

    // *******************
    // Non host Alignment
    // *******************

    // Build filter
    let deuterostome_path = resolve_optional_path(&config.args.deuterostome_list, &config.cwd)?;
    let taxon_whitelist_path = resolve_optional_path(&config.args.taxon_whitelist, &config.cwd)?;
    let taxon_blacklist_path = resolve_optional_path(&config.args.taxon_blacklist, &config.cwd)?;

    let should_keep_filter = Arc::new(
        build_should_keep_filter(
            deuterostome_path,
            taxon_whitelist_path,
            taxon_blacklist_path,
        )
            .await
            .map_err(|e| PipelineError::Other(e))?,
    );

    //write de-interleaved non-host file. easier for the scatter-gather
    let (r1_rx, r2_rx_opt, deinterleave_handle) =
        deinterleave_fastq_stream(subsampled_stream, paired, config.base_buffer_size * 4)
            .await
            .map_err(|e| PipelineError::Other(e))?;

    let non_host_temp_dir = choose_temp_dir(
        config.input_size * 2,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        true,
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Compute paths after temp_dir but before moving rx (only borrows with as_ref())
    let non_host_r1_path = non_host_temp_dir.path().join("nonhost_R1.fastq");
    let non_host_r2_path_opt = r2_rx_opt
        .as_ref()
        .map(|_| non_host_temp_dir.path().join("nonhost_R2.fastq"));

    info!("Checkpoint: deinterleaving and writing non-host R1/R2 to temp (forces upstream completion). Writing to: {:?}", non_host_temp_dir);
    final_temp_dirs.push(non_host_temp_dir);

    let r1_write_task_fut = write_parse_output_to_file(
        &non_host_r1_path,
        ReceiverStream::new(r1_rx),
        Some(config.base_buffer_size * 4),
    );

    let r2_write_task_fut_opt = r2_rx_opt.map(|r2_rx| {
        let r2_path = non_host_r2_path_opt
            .as_ref()
            .expect("R2 path should exist when r2_rx_opt is Some");

        write_parse_output_to_file(
            r2_path,
            ReceiverStream::new(r2_rx),
            Some(config.base_buffer_size * 4),
        )
    });

    // start both writer tasks concurrently
    let (r1_write_task, r2_write_task_opt) = tokio::try_join!(r1_write_task_fut, async {
        match r2_write_task_fut_opt {
            Some(fut) => fut.await.map(Some),
            None => Ok(None),
        }
    })?;

    // then wait for both background writes and deinterleave together
    tokio::try_join!(
        async {
            r1_write_task.await?;
            Ok::<(), anyhow::Error>(())
        },
        async {
            if let Some(task) = r2_write_task_opt {
                task.await?;
            }
            Ok::<(), anyhow::Error>(())
        },
        async {
            deinterleave_handle.await?;
            Ok::<(), anyhow::Error>(())
        },
    )?;

    let duplicate_clusters = dedup_finalize_task
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
            "Dedup finalization task panicked: {}",
            e
        ))
        })?
        .map_err(PipelineError::Other)?;

    for rx in hisat2_cleanup_receivers {
        rx.await??;
    }

    info!(
        "Checkpoint complete. Files: r1:{:?}   r2:{:?}",
        non_host_r1_path, non_host_r2_path_opt
    );

    // Early cleanup upstream temp dirs (host filter, etc.)
    for td in &host_filter_temp_dirs {
        if let Err(e) = std::fs::remove_dir_all(td.path()) {
            if e.kind() != std::io::ErrorKind::NotFound {
                warn!("Early temp dir cleanup failed (non-fatal): {}", e);
            } else {
                debug!("Temp dir already gone: {}", td.path().display());
            }
        } else {
            info!("Early cleanup: removed temp dir {}", td.path().display());
        }
    }
    host_filter_temp_dirs.clear();

    // Clone paths **before** moving them into async copy tasks
    let non_host_r1_path_clone = non_host_r1_path.clone();
    let non_host_r2_path_opt_clone = non_host_r2_path_opt.clone();

    // Async copy to out_dir (non-blocking)
    let out_dir_clone = out_dir.clone();
    let r1_copy_task = tokio::spawn(async move {
        let dest_r1 = out_dir_clone.join("nonhost_R1.fastq");
        if let Err(e) = tokio::fs::copy(&non_host_r1_path_clone, &dest_r1).await {
            warn!("Failed to copy non-host R1 to out_dir: {}", e);
        } else {
            info!("Copied non-host R1 to {}", dest_r1.display());
        }
        Ok(())
    });
    cleanup_tasks.push(r1_copy_task);

    if let Some(r2_path) = non_host_r2_path_opt_clone {
        let r2_path_clone = r2_path.clone();
        let out_dir_clone = out_dir.clone();
        let r2_copy_task = tokio::spawn(async move {
            let dest_r2 = out_dir_clone.join("nonhost_R2.fastq");
            if let Err(e) = tokio::fs::copy(&r2_path_clone, &dest_r2).await {
                warn!("Failed to copy non-host R2 to out_dir: {}", e);
            } else {
                info!("Copied non-host R2 to {}", dest_r2.display());
            }
            Ok(())
        });
        cleanup_tasks.push(r2_copy_task);
    }


    // =====================================================================
    // Non-host Alignment
    //
    // The non-host FASTQ files are now fully materialized above, so NT and
    // NR have independent, stable inputs.
    //
    // Execution topology:
    //
    // Single:
    //     NT non-host alignment
    //       -> paf_to_m8
    //       -> sort NT m8
    //       -> NT call_hits / fanout
    //       -> NR non-host alignment
    //       -> sort NR m8
    //       -> NR call_hits / fanout
    //
    // Distributed:
    //     NT non-host alignment
    //       -> paf_to_m8
    //       -> sort NT m8
    //
    //                 || concurrent ||
    //
    //     NR distributed alignment
    //       -> NR m8
    //       -> sort NR m8
    //
    // We intentionally overlap ONLY the two non-host alignment branches.
    // Downstream NT/NR hit processing remains unchanged.
    // =====================================================================

    // ---------------------------------------------------------------------
    // NT non-host alignment future
    //
    // This includes:
    //   nt_non_host_align
    //   paf_to_m8
    //   sort_m8_by_read_id(nt)
    //
    // The future does not execute until awaited below. That lets distributed
    // mode start this branch concurrently with the NR branch.
    // ---------------------------------------------------------------------

    let nt_non_host_align_fut = async {
        info!("[run] NT non-host alignment phase started");

        let (
            non_host_nt_out_stream,
            mut non_host_nt_cleanup_tasks,
            mut non_host_nt_cleanup_receivers,
        ) = nt_non_host_align(
            config.clone(),
            non_host_r1_path.clone(),
            non_host_r2_path_opt.clone(),
        )
            .await?;

        let nt_m8_file_path = out_dir.join("gsnap.m8");

        let (
            m8_stream,
            nt_m8_write_task,
            mut m8_cleanup_tasks,
            m8_cleanup_receivers,
        ) = paf_to_m8(
            config.clone(),
            non_host_nt_out_stream,
            nt_m8_file_path,
        )
            .await?;

        // PAF -> m8 tasks belong to the NT branch.
        non_host_nt_cleanup_tasks.append(&mut m8_cleanup_tasks);
        non_host_nt_cleanup_receivers.extend(m8_cleanup_receivers);

        // -------------------------------------------------------------
        // Sort NT m8 by read ID before call_hits_m8.
        //
        // This guarantees consecutive lines per read and therefore
        // preserves the streaming group-by requirement downstream.
        // -------------------------------------------------------------

        let m8_sorted_start = Instant::now();

        let m8_sorted = sort_m8_by_read_id(
            config.clone(),
            m8_stream,
            "nt",
        )
            .await?;

        info!(
            "[run] sort_m8_by_read_id(nt) completed after {:?}",
            m8_sorted_start.elapsed()
        );

        // The NT branch
        // does not complete until its PAF -> m8 completion receivers have
        // completed successfully.
        for rx in non_host_nt_cleanup_receivers.drain(..) {
            rx.await??;
        }

        info!("[run] NT non-host alignment phase complete");

        Ok::<
            (
                ReceiverStream<ParseOutput>,
                JoinHandle<anyhow::Result<(), anyhow::Error>>,
                Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
                Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
            ),
            PipelineError,
        >((
            m8_sorted,
            nt_m8_write_task,
            non_host_nt_cleanup_tasks,
            non_host_nt_cleanup_receivers,
        ))
    };

    let distributed_nr_run_dir =
        config.efs_runs_dir.join(&config.run_id);

    if matches!(config.execution_mode, ExecutionMode::Distributed) {
        info!(
        "[run] Distributed NR EFS run directory: {}",
        distributed_nr_run_dir.display()
    );
    }

    // ---------------------------------------------------------------------
    // NR non-host alignment future
    //
    // This includes:
    //   nr_non_host_align
    //   NR m8 fanout
    //   NR m8 file writer
    //   sort_m8_by_read_id(nr)
    //
    // In distributed mode, nr_non_host_align() contains the distributed
    // scheduler / worker execution, so this future remains pending while
    // that remote work proceeds.
    // ---------------------------------------------------------------------

    let nr_non_host_align_fut = async {
        info!("[run] NR non-host alignment phase started");

        let (
            nr_non_host_m8_stream,
            mut nr_non_host_cleanup_tasks,
            nr_non_host_cleanup_receivers,
            nr_non_host_align_temp_dirs,
        ) = nr_non_host_align(
            config.clone(),
            non_host_r1_path.clone(),
            non_host_r2_path_opt.clone(),
            sample_base.clone(),
        )
            .await?;

        // Fan out the completed/streaming NR m8:
        //
        //   branch 0 -> sort -> call_hits_m8
        //   branch 1 -> durable nr.m8 file
        //
        let (nr_m8_streams, paf_to_m8_stream_done_rx) = fanout_to_channels(
            ReceiverStream::new(nr_non_host_m8_stream),
            2,
            "nr_m8_stream",
            &config,
            StreamDataType::JustBytes,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?;

        let mut nr_m8_streams_it = nr_m8_streams.into_iter();

        let nr_m8_stream = ReceiverStream::new(
            nr_m8_streams_it
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );

        let nr_m8_file_stream = ReceiverStream::new(
            nr_m8_streams_it
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );

        let nr_m8_file_path = out_dir.join("rapsearch2.m8");

        let write_task = write_byte_stream_to_file(
            &nr_m8_file_path,
            nr_m8_file_stream,
            config.clone(),
            StreamDataType::JustBytes,
            "nr_m8_file",
            true,
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

        nr_non_host_cleanup_tasks.push(write_task);

        // The fanout completion receiver is also owned by this branch.
        let mut nr_non_host_cleanup_receivers = nr_non_host_cleanup_receivers;
        nr_non_host_cleanup_receivers.push(paf_to_m8_stream_done_rx);

        // -------------------------------------------------------------
        // Sort NR m8 by read ID before call_hits_m8.
        //
        // Guarantees consecutive lines per read and therefore enables
        // the downstream streaming group-by.
        // -------------------------------------------------------------

        let nr_sort_start = Instant::now();

        let nr_m8_sorted = sort_m8_by_read_id(
            config.clone(),
            nr_m8_stream,
            "nr",
        )
            .await?;

        info!(
            "[run] sort_m8_by_read_id(nr) completed after {:?}",
            nr_sort_start.elapsed()
        );

        info!("[run] NR non-host alignment phase complete");

        Ok::<
            (
                ReceiverStream<ParseOutput>,
                Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>,
                Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>,
                Vec<TempDir>,
            ),
            PipelineError,
        >((
            nr_m8_sorted,
            nr_non_host_cleanup_tasks,
            nr_non_host_cleanup_receivers,
            nr_non_host_align_temp_dirs,
        ))
    };

    // ---------------------------------------------------------------------
    //
    // SINGLE:
    //     Await NT completely before starting NR.
    //
    // DISTRIBUTED:
    //     Drive both futures concurrently.
    //
    // ---------------------------------------------------------------------

    let (
        (
            m8_sorted,
            nt_m8_write_task,
            mut nt_non_host_cleanup_tasks,
            mut nt_non_host_cleanup_receivers,
        ),
        (
            nr_m8_sorted,
            mut nr_non_host_cleanup_tasks,
            mut nr_non_host_cleanup_receivers,
            nr_non_host_align_temp_dirs,
        ),
    ) = match config.execution_mode {
        ExecutionMode::Single => {
            info!(
                "[run] Non-host alignment mode=single: \
                 running NT then NR serially"
            );

            let nt_result = nt_non_host_align_fut.await?;

            let nr_result = nr_non_host_align_fut.await?;

            (nt_result, nr_result)
        }

        ExecutionMode::Distributed => {
            info!(
                "[run] Non-host alignment mode=distributed: \
                 running NT and NR concurrently"
            );

            tokio::try_join!(
                nt_non_host_align_fut,
                nr_non_host_align_fut,
            )?
        }
    };

    // ---------------------------------------------------------------------
    // Transfer branch-owned cleanup resources into the run-level cleanup
    // collections.
    // ---------------------------------------------------------------------

    cleanup_tasks.append(&mut nt_non_host_cleanup_tasks);
    cleanup_receivers.append(&mut nt_non_host_cleanup_receivers);

    cleanup_tasks.append(&mut nr_non_host_cleanup_tasks);
    cleanup_receivers.append(&mut nr_non_host_cleanup_receivers);

    final_temp_dirs.extend(nr_non_host_align_temp_dirs);

    // =====================================================================
    // NT downstream processing
    // =====================================================================

    let (lineage_map, acc2taxid_map) = taxonomy_handle.await??;

    let nt_concurrency =
        compute_phase_concurrency(
            &config,
            "call_hits_nt",
            1.0,
            3.5,
            64,
            16,
        );

    info!(
        "call hits nt concurrency {}",
        nt_concurrency
    );

    let nt_call_hits_start = Instant::now();

    let (
        nt_pairs,
        mut nt_call_cleanup_tasks,
        mut nt_call_cleanup_receivers,
    ) = call_hits_m8(
        config.clone(),
        m8_sorted,
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        36,
        nt_concurrency,
        "nt".to_string(),
    )
        .await?;

    info!(
        "[run] call_hits_m8(nt) returned after {:?}",
        nt_call_hits_start.elapsed()
    );

    cleanup_tasks.append(&mut nt_call_cleanup_tasks);
    cleanup_receivers.append(&mut nt_call_cleanup_receivers);

    // ---------------------------------------------------------------------
    // One fanout of paired NT hits.
    // ---------------------------------------------------------------------

    let nt_split_start = Instant::now();

    info!(
        "[run] fanout_to_channels(nt_pairs) x5 ReducedRead"
    );

    let (nt_pair_rxs, nt_pairs_done_rx) = fanout_to_channels(
        nt_pairs,
        6,
        "nt_pairs",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(nt_pairs_done_rx);

    info!(
        "[run] fanout_to_channels(nt_pairs) ready after {:?}",
        nt_split_start.elapsed()
    );

    let mut nt_pair_rxs = nt_pair_rxs.into_iter();
    let nt_pairs_taxon = ReceiverStream::new(
        nt_pair_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );
    let nt_pairs_summarize = ReceiverStream::new(
        nt_pair_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );
    let nt_pairs_initial = ReceiverStream::new(
        nt_pair_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );
    let nt_pairs_blast = ReceiverStream::new(
        nt_pair_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );
    let nt_pairs_taxid_initial = ReceiverStream::new(
        nt_pair_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );
    let nt_pairs_artifact = ReceiverStream::new(
        nt_pair_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?,
    );


    let (nt_deduped_tx, nt_deduped_rx) =
        mpsc::channel::<ParseOutput>(config.base_buffer_size);

    let (nt_hitsummary_tx, nt_hitsummary_rx) =
        mpsc::channel::<ParseOutput>(config.base_buffer_size);

    let nt_artifact_split_task = tokio::spawn(async move {
        let mut pairs = nt_pairs_artifact;

        while let Some(pair) = pairs.next().await {
            nt_deduped_tx
                .send(ParseOutput::Bytes(Bytes::from(pair.dedup)))
                .await
                .map_err(|_| anyhow!(
                "gsnap.deduped.m8 artifact receiver dropped"
            ))?;

            nt_hitsummary_tx
                .send(ParseOutput::Bytes(Bytes::from(pair.summary)))
                .await
                .map_err(|_| anyhow!(
                "gsnap.hitsummary.tab artifact receiver dropped"
            ))?;
        }

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(nt_artifact_split_task);

    let nt_deduped_path =
        out_dir.join("gsnap.deduped.m8");

    let nt_hitsummary_path =
        out_dir.join("gsnap.hitsummary.tab");

    let nt_deduped_write_task =
        write_byte_stream_to_file(
            &nt_deduped_path,
            ReceiverStream::new(nt_deduped_rx),
            config.clone(),
            StreamDataType::JustBytes,
            "gsnap_deduped_m8",
            true,
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let nt_hitsummary_write_task =
        write_byte_stream_to_file(
            &nt_hitsummary_path,
            ReceiverStream::new(nt_hitsummary_rx),
            config.clone(),
            StreamDataType::JustBytes,
            "gsnap_hitsummary",
            true,
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(nt_deduped_write_task);
    cleanup_tasks.push(nt_hitsummary_write_task);

    info!("[run] NT pair wiring complete");
    info!("[run] nt_pairs_taxon      -> generate_taxon_counts");
    info!("[run] nt_pairs_summarize  -> summarize_hits");
    info!("[run] nt_pairs_initial    -> collect_hit_summary / initial annotate maps");
    info!("[run] nt_pairs_blast      -> blast_contigs");
    info!("[run] nt_pairs_taxid      -> generate_taxid_fasta (experimental)");
    info!("[run] nt_pairs_artifact   -> gsnap.deduped.m8 + gsnap.hitsummary.tab");

    let nt_hit_summary_handle = tokio::spawn({
        let config = config.clone();
        let nt_pairs_summarize = nt_pairs_summarize;
        let duplicate_clusters = duplicate_clusters.clone();

        async move {
            let start = Instant::now();

            info!("[run] summarize_hits(nt) started");

            let res = summarize_hits(
                config.clone(),
                nt_pairs_summarize,
                duplicate_clusters,
                0,
            )
                .await;

            info!(
                "[run] summarize_hits(nt) finished after {:?}",
                start.elapsed()
            );

            res
        }
    });

    let nt_counts_task = tokio::spawn({
        let config = config.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let nt_pairs_taxon = nt_pairs_taxon;

        async move {
            let start = Instant::now();

            info!(
                "[run] generate_taxon_counts(NT) started"
            );

            let res = generate_taxon_counts(
                config,
                nt_pairs_taxon,
                duplicate_clusters,
                should_keep_filter,
                "NT".to_string(),
                None,
            )
                .await;

            info!(
                "[run] generate_taxon_counts(NT) finished after {:?}",
                start.elapsed()
            );

            res
        }
    });

    // =====================================================================
    // NR downstream processing
    // =====================================================================

    let nr_concurrency = compute_phase_concurrency(
        &config,
        "call_hits_nr",
        1.0, // ~1 GB per thread max (mostly transient)
        3.5,
        64, // higher cap for CPU-bound phase
        16, // min for meaningful parallelism
    );

    info!(
        "call hits nr concurrency {}",
        nr_concurrency
    );

    let nr_call_hits_start = Instant::now();

    let (
        nr_pairs,
        mut nr_call_cleanup_tasks,
        mut nr_call_cleanup_receivers,
    ) = call_hits_m8(
        config.clone(),
        nr_m8_sorted,
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        0,
        nr_concurrency,
        "nr".to_string(),
    )
        .await?;

    info!(
        "[run] call_hits_m8(nr) returned after {:?}",
        nr_call_hits_start.elapsed()
    );

    cleanup_tasks.append(&mut nr_call_cleanup_tasks);
    cleanup_receivers.append(&mut nr_call_cleanup_receivers);

    // ---------------------------------------------------------------------
    // Single 6-way fanout of ReducedRead.
    //
    // Each logical pair remains together for every downstream consumer.
    // ---------------------------------------------------------------------

    let nr_split_start = Instant::now();

    info!(
        "[run] starting fanout_to_channels for NR pairs (5 private channels)"
    );

    let (nr_pair_rxs, nr_pairs_done_rx) = fanout_to_channels(
        nr_pairs,
        6,
        "nr_pairs",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    cleanup_receivers.push(nr_pairs_done_rx);

    info!(
        "[run] fanout_to_channels(nr_pairs) ready after {:?} with 5 private channels",
        nr_split_start.elapsed()
    );

    let mut nr_pair_rxs_iter = nr_pair_rxs.into_iter();
    let nr_pairs_taxon =
        ReceiverStream::new(nr_pair_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_pairs_summarize =
        ReceiverStream::new(nr_pair_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_pairs_initial =
        ReceiverStream::new(nr_pair_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_pairs_blast =
        ReceiverStream::new(nr_pair_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_pairs_taxid_initial =
        ReceiverStream::new(nr_pair_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_pairs_artifact =
        ReceiverStream::new(
            nr_pair_rxs_iter
                .next()
                .ok_or(PipelineError::EmptyStream)?,
        );

    info!("[run] NR pair wiring complete");
    info!("[run] nr_pairs_taxon     -> generate_taxon_counts");
    info!("[run] nr_pairs_summarize -> summarize_hits");
    info!("[run] nr_pairs_initial   -> initial accession map");
    info!("[run] nr_pairs_blast     -> blast_contigs");
    info!("[run] nr_pairs_taxid     -> generate_taxid_fasta");
    info!("[run] nr_pairs_artifact     -> dedup and hitsummary ");


    let (nr_deduped_tx, nr_deduped_rx) =
        mpsc::channel::<ParseOutput>(config.base_buffer_size);

    let (nr_hitsummary_tx, nr_hitsummary_rx) =
        mpsc::channel::<ParseOutput>(config.base_buffer_size);

    let nr_artifact_split_task = tokio::spawn(async move {
        let mut pairs = nr_pairs_artifact;

        while let Some(pair) = pairs.next().await {
            nr_deduped_tx
                .send(ParseOutput::Bytes(Bytes::from(pair.dedup)))
                .await
                .map_err(|_| {
                    anyhow!(
                    "rapsearch2.deduped.m8 artifact receiver dropped"
                )
                })?;

            nr_hitsummary_tx
                .send(ParseOutput::Bytes(Bytes::from(pair.summary)))
                .await
                .map_err(|_| {
                    anyhow!(
                    "rapsearch2.hitsummary.tab artifact receiver dropped"
                )
                })?;
        }

        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(nr_artifact_split_task);

    let nr_deduped_path =
        out_dir.join("rapsearch2.deduped.m8");

    let nr_hitsummary_path =
        out_dir.join("rapsearch2.hitsummary.tab");

    let nr_deduped_write_task =
        write_byte_stream_to_file(
            &nr_deduped_path,
            ReceiverStream::new(nr_deduped_rx),
            config.clone(),
            StreamDataType::JustBytes,
            "rapsearch2_deduped_m8",
            true,
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let nr_hitsummary_write_task =
        write_byte_stream_to_file(
            &nr_hitsummary_path,
            ReceiverStream::new(nr_hitsummary_rx),
            config.clone(),
            StreamDataType::JustBytes,
            "rapsearch2_hitsummary",
            true,
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(nr_deduped_write_task);
    cleanup_tasks.push(nr_hitsummary_write_task);

    let nr_hit_summary_handle = tokio::spawn(summarize_hits(
        config.clone(),
        nr_pairs_summarize,
        duplicate_clusters.clone(),
        0,
    ));

    let nr_counts_task = tokio::spawn({
        let config = config.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let nr_pairs_taxon = nr_pairs_taxon;
        async move {
            let start = Instant::now();
            info!("[run] generate_taxon_counts(NR) started");
            let res = generate_taxon_counts(
                config,
                nr_pairs_taxon,
                duplicate_clusters,
                should_keep_filter,
                "NR".to_string(),
                None,
            )
                .await;
            info!(
        "[run] generate_taxon_counts(NR) finished after {:?}",
        start.elapsed()
    );
            res
        }
    });

    let annot_concurrency =
        compute_phase_concurrency(&config, "generate_annotated_fasta", 0.4, 4.0, 128, 8);

    // Read non-host R1/R2 as proper parsed Fastq records
    let (interleaved_rx, read_stats_task) = read_fastq(
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        PairingMode::Relaxed,
        None, // technology
        u64::MAX,
        None,
        None,
        "run - Read non-host R1/R2 as proper parsed Fastq records",
        &config,
    )
        .map_err(|e| PipelineError::Other(e))?;

    // Wrap the stats task so it matches cleanup_tasks type
    let stats_wrapper = tokio::spawn(async move {
        match read_stats_task.await {
            Ok(Ok(_stats)) => Ok(()),
            Ok(Err(e)) => Err(e),
            Err(join_err) => Err(anyhow::anyhow!("read_stats_task join failed: {}", join_err)),
        }
    });
    cleanup_tasks.push(stats_wrapper);

    let interleaved_stream = ReceiverStream::new(interleaved_rx);

    // Build the small accession maps for NT/NR in parallel
    let (nt_map, nr_map) = tokio::try_join!(
    collect_hit_summary_to_accession_map_concurrent(config.clone(), nt_pairs_initial),
    collect_hit_summary_to_accession_map_concurrent(config.clone(), nr_pairs_initial)
)?;

    // Generate the annotated FASTA stream
    let (
        annotated_rx,
        unidentified_rx,
        unique_unidentified_rx,
        mut annot_tasks,
        mut annot_rxs,
    ) = generate_annotated_fasta(
        config.clone(),
        interleaved_stream,
        duplicate_clusters.clone(),
        nt_map,
        nr_map,
        annot_concurrency,
    )
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
        "Annotated FASTA from files failed: {}",
        e
    ))
        })?;

    cleanup_tasks.append(&mut annot_tasks);
    cleanup_receivers.append(&mut annot_rxs);

    let assembly_dir = out_dir.join("assembly");

    let annotated_path = out_dir.join("annotated_merged.fa");
    let unidentified_path = out_dir.join("unidentified.fa");
    let unique_unidentified_path = assembly_dir.join("unique_unidentified.fa");

    info!("[run] starting fanout_to_channels for annotated FASTA (2 private channels)");

    let (annotated_rxs, annotated_done_rx) = fanout_to_channels(
        ReceiverStream::new(annotated_rx),
        2,
        "initial_annotated",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    let mut annotated_rxs_iter = annotated_rxs.into_iter();
    let initial_annotated_file_stream = annotated_rxs_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;
    let initial_annotated_taxon_stream = annotated_rxs_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    let (unidentified_rxs, unidentified_done_rx) = fanout_to_channels(
        ReceiverStream::new(unidentified_rx),
        2,
        "initial_unidentified",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;

    let mut unidentified_rxs_iter = unidentified_rxs.into_iter();
    let initial_unidentified_file_stream = unidentified_rxs_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;
    let initial_unidentified_taxon_stream = unidentified_rxs_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    // Part of Experimental phase but starting this as early as possible to avvoid stalling the fanout
    let experimental_taxid_handle = tokio::spawn({
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let mapped = ReceiverStream::new(initial_annotated_taxon_stream);
        let unid = ReceiverStream::new(initial_unidentified_taxon_stream);
        let nt = nt_pairs_taxid_initial;
        let nr = nr_pairs_taxid_initial;
        async move {
            generate_taxid_fasta(config, mapped, unid, nt, nr, lineage_map).await
        }
    });

    // Write tasks are now eager consumers → EOF propagates instantly
    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(initial_annotated_file_stream),
        annotated_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_annotated_file_stream_to_write",
    ));

    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(initial_unidentified_file_stream),
        unidentified_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_unidentified_file_stream_to_write",
    ));

    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(unique_unidentified_rx),
        unique_unidentified_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "unique_unidentified_file_stream_to_write",
    ));

    info!("[run] annotated FASTA wiring complete (fanout_to_channels used everywhere)");
    // *******************
    // Post-processing
    // *******************

    // Assembly stats
    let assembly_out_dir = out_dir.join("assembly");

    let (
        assembly_outputs,
        assembly_bam_out_stream,
        post_assembly_cleanup_tasks,
        post_assembly_cleanup_receivers,
        post_assembly_temp_files,
        post_assembly_temp_dir,
    ) = process_assembly(
        config.clone(),
        &assembly_out_dir,
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        duplicate_clusters.clone(),
        paired,
        4, // this can become as CLI arg if needed
    )
        .await?;

    cleanup_tasks.extend(post_assembly_cleanup_tasks);
    // cleanup_receivers.extend(post_assembly_cleanup_receivers);
    temp_files.extend(post_assembly_temp_files);
    final_temp_dirs.push(post_assembly_temp_dir);

    let coverage_json_path = assembly_out_dir.join("contig_coverage.json");
    let coverage_summary_csv_path =
        assembly_out_dir.join("contig_coverage_summary.csv");

    let _generate_assembly_coverage_result =
        generate_assembly_coverage(
            config.clone(),
            assembly_bam_out_stream,
            &assembly_outputs.contigs_fasta,
            &coverage_json_path,
            &coverage_summary_csv_path,
        )
            .await
            .map_err(|e| anyhow!(
        "generate_assembly_coverage failed: {}",
        e
    ))?;

    for rx in post_assembly_cleanup_receivers {
        rx.await??;
    }

    let nt_file = config
        .args
        .nt
        .clone()
        .ok_or(PipelineError::MissingArgument(
            "NT file required for fasta extraction".into(),
        ))?;

    let nt_offset_db_file =
        config
            .args
            .nt_offset_db
            .clone()
            .ok_or(PipelineError::MissingArgument(
                "NT offset DB file required for fasta extraction".into(),
            ))?;

    let nr_file = config
        .args
        .nr
        .clone()
        .ok_or(PipelineError::MissingArgument(
            "NR file required for fasta extraction".into(),
        ))?;

    let nr_offset_db_file =
        config
            .args
            .nr_offset_db
            .clone()
            .ok_or(PipelineError::MissingArgument(
                "NR offset DB file required for fasta extraction".into(),
            ))?;

    // Blast contigs
    // NT and NR prep in parallel
    let nt_prep_handle = tokio::spawn({
        let config = config.clone();
        let nt_file = nt_file.clone();
        let nt_offset_db_file = nt_offset_db_file.clone();
        async move {
            let (nt_read_dict_map, nt_accession_dict_noarc, nt_selected_genera, _nt_total_reads) =
                nt_hit_summary_handle
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!("NT hit summary task panicked: {}", e))
                    })?
                    .map_err(|e| {
                        PipelineError::Other(anyhow!("NT hit summary parsing failed: {}", e))
                    })?;

            let (nt_ref_fasta_path, nt_ref_fasta_temp_dir) =
                build_reference_fasta_from_selected_genera(
                    config.clone(),
                    &nt_selected_genera,
                    NT_TAG,
                    &PathBuf::from(nt_file),
                    &PathBuf::from(nt_offset_db_file),
                )
                    .await
                    .map_err(|e| PipelineError::Other(e.into()))?;
            Ok::<
                (
                    AHashMap<String, Arc<ReadHit>>,
                    Arc<AHashMap<String, AccessionHit>>,
                    PathBuf,
                    TempDir,
                ),
                PipelineError,
            >((
                nt_read_dict_map,
                Arc::new(nt_accession_dict_noarc),
                nt_ref_fasta_path,
                nt_ref_fasta_temp_dir,
            ))
        }
    });

    let nr_prep_handle = tokio::spawn({
        let config = config.clone();
        let nr_file = nr_file.clone();
        let nr_offset_db_file = nr_offset_db_file.clone();
        async move {
            let (nr_read_dict_map, nr_accession_dict_noarc, nr_selected_genera, _nr_total_reads) =
                nr_hit_summary_handle
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!("NR hit summary task panicked: {}", e))
                    })?
                    .map_err(|e| {
                        PipelineError::Other(anyhow!("NR hit summary parsing failed: {}", e))
                    })?;

            let (nr_ref_fasta_path, nr_ref_fasta_temp_dir) =
                build_reference_fasta_from_selected_genera(
                    config.clone(),
                    &nr_selected_genera,
                    NR_TAG,
                    &PathBuf::from(nr_file),
                    &PathBuf::from(nr_offset_db_file),
                )
                    .await
                    .map_err(|e| PipelineError::Other(e.into()))?;
            Ok::<
                (
                    AHashMap<String, Arc<ReadHit>>,
                    Arc<AHashMap<String, AccessionHit>>,
                    PathBuf,
                    TempDir,
                ),
                PipelineError,
            >((
                nr_read_dict_map,
                Arc::new(nr_accession_dict_noarc),
                nr_ref_fasta_path,
                nr_ref_fasta_temp_dir,
            ))
        }
    });

    let (nt_prep_res, nr_prep_res) = tokio::try_join!(nt_prep_handle, nr_prep_handle)
        .map_err(|e| PipelineError::Other(anyhow!("Reference prep task panicked: {e}")))?;

    let (nt_read_dict_map, nt_accession_dict, nt_ref_fasta_path, nt_ref_fasta_temp_dir) =
        nt_prep_res?;
    let (nr_read_dict_map, nr_accession_dict, nr_ref_fasta_path, nr_ref_fasta_temp_dir) =
        nr_prep_res?;

    let nt_read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>> =
        Arc::new(Mutex::new(nt_read_dict_map));
    let nr_read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>> =
        Arc::new(Mutex::new(nr_read_dict_map));

    info!("NT ref path: {}", nt_ref_fasta_path.display());
    info!("NR ref path: {}", nr_ref_fasta_path.display());

    final_temp_dirs.push(nt_ref_fasta_temp_dir);
    final_temp_dirs.push(nr_ref_fasta_temp_dir);

    let (nt_counts, nr_counts) = tokio::try_join!(
        async {
            nt_counts_task
                .await
                .map_err(|e| PipelineError::Other(anyhow!(e.to_string())))
        },
        async {
            nr_counts_task
                .await
                .map_err(|e| PipelineError::Other(anyhow!(e.to_string())))
        },
    )?;

    let nt_counts = nt_counts?;
    let nr_counts = nr_counts?;

    // Persist the individual initial NT/NR taxon-count artifacts required
    // by the legacy non-host-alignment workflow.
    //
    // nt_counts / nr_counts are already the completed outputs of
    // generate_taxon_counts(), so no additional stream or counting pass
    // is required.
    let nt_counts_path =
        out_dir.join("gsnap_counts_with_dcr.json");

    let nr_counts_path =
        out_dir.join("rapsearch2_counts_with_dcr.json");

    let nt_counts_json = serde_json::to_vec(&nt_counts)
        .map_err(|e| {
            PipelineError::Other(anyhow!(
            "Failed to serialize gsnap_counts_with_dcr.json: {}",
            e
        ))
        })?;

    let nr_counts_json = serde_json::to_vec(&nr_counts)
        .map_err(|e| {
            PipelineError::Other(anyhow!(
            "Failed to serialize rapsearch2_counts_with_dcr.json: {}",
            e
        ))
        })?;

    let nt_counts_write_task = {
        let path = nt_counts_path.clone();

        tokio::spawn(async move {
            tokio::fs::write(&path, nt_counts_json)
                .await
                .map_err(|e| {
                    anyhow!(
                    "Failed to write {}: {}",
                    path.display(),
                    e
                )
                })?;

            info!(
            "Wrote initial NT taxon counts to {}",
            path.display()
        );

            Ok::<(), anyhow::Error>(())
        })
    };

    let nr_counts_write_task = {
        let path = nr_counts_path.clone();

        tokio::spawn(async move {
            tokio::fs::write(&path, nr_counts_json)
                .await
                .map_err(|e| {
                    anyhow!(
                    "Failed to write {}: {}",
                    path.display(),
                    e
                )
                })?;

            info!(
            "Wrote initial NR taxon counts to {}",
            path.display()
        );

            Ok::<(), anyhow::Error>(())
        })
    };

    cleanup_tasks.push(nt_counts_write_task);
    cleanup_tasks.push(nr_counts_write_task);

    let combined_path =
        out_dir.join("taxon_counts_with_dcr.json");

    let nt_counts_for_combine = nt_counts.clone();
    let nr_counts_for_combine = nr_counts.clone();

    let combine_handle = tokio::spawn(async move {
        let (_combined_path, write_json_task) = combine_taxon_counts(
            &nt_counts_for_combine,
            &nr_counts_for_combine,
            &[], // no merged yet
            combined_path,
        )
            .await
            .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts failed: {}", e)))?;

        write_json_task.await.map_err(|e| {
            PipelineError::Other(anyhow!("combine_taxon_counts write task failed: {}", e))
        })?
    });

    cleanup_tasks.push(combine_handle);



    let nt_blast_concurrency = compute_phase_concurrency(
        &config,
        "nt_blast_contigs",
        2.0,
        5.0,
        128,
        8,
    );
    info!("NT blast_contigs concurrency: {}", nt_blast_concurrency);

    let nt_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let read2contig = assembly_outputs.read2contig.clone();
        let nt_pairs_blast = nt_pairs_blast;
        let sample_base = sample_base.clone();

        async move {
            blast_contigs(
                config,
                NT_TAG,
                &sample_base,
                nt_pairs_blast,
                nt_read_dict.clone(),
                nt_accession_dict,
                nt_counts,
                &contigs_fasta_path,
                read2contig,
                &nt_ref_fasta_path,
                duplicate_clusters,
                lineage_map,
                should_keep_filter,
                4,
                nt_blast_concurrency,
            )
                .await
        }
    });

    let nr_blast_concurrency = compute_phase_concurrency(
        &config,
        "nr_blast_contigs",
        2.0,
        5.0,
        128,
        8,
    );
    info!("NR blast_contigs concurrency: {}", nr_blast_concurrency);

    let nr_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let read2contig = assembly_outputs.read2contig.clone();
        let nr_pairs_blast = nr_pairs_blast;
        let sample_base = sample_base.clone();

        async move {
            blast_contigs(
                config,
                NR_TAG,
                &sample_base,
                nr_pairs_blast,
                nr_read_dict.clone(),
                nr_accession_dict,
                nr_counts,
                &contigs_fasta_path,
                read2contig,
                &nr_ref_fasta_path,
                duplicate_clusters,
                lineage_map,
                should_keep_filter,
                4,
                nr_blast_concurrency,
            )
                .await
        }
    });

    // ── NT post: fanout refined ReducedRead pairs ─────────────────────────────
    // Branches:
    //   merge     → compute_merged_taxon_counts
    //   map       → collect_m8_accession_map (use pair.dedup)
    //   viz       → coverage viz m8 side (use pair.dedup)
    //   taxid     → generate_taxid_fasta (use pair.summary)
    //   coverage  → coverage viz hitsummary (use pair.summary)
    let nt_post_handle = tokio::spawn({
        let config = config.clone();
        async move {
            let nt_res = nt_handle.await.map_err(|e| {
                PipelineError::Other(anyhow!("NT blast_contigs task panicked: {e}"))
            })??;

            let (
                nt_read_dict,
                nt_refined_counts,
                nt_contig_summary,
                nt_refined_pairs_out,
                nt_refined_m8_top_stream_out,
                nt_cleanup_tasks,
                nt_cleanup_receivers,
                nt_temp_files,
            ) = nt_res;

            let mut cleanup_tasks = nt_cleanup_tasks;
            let mut cleanup_receivers = nt_cleanup_receivers;
            let temp_files = nt_temp_files;

            let (nt_pair_rxs, nt_pairs_done_rx) = fanout_to_channels(
                ReceiverStream::new(nt_refined_pairs_out),
                4,
                "nt_refined_pairs",
                &config,
                StreamDataType::JustBytes,
            )
                .await?;
            cleanup_receivers.push(nt_pairs_done_rx);

            let mut it = nt_pair_rxs.into_iter();
            let nt_pairs_merge = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
            let nt_pairs_map = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
            let nt_pairs_taxid_refined = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
            let nt_pairs_coverage =
                ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);

            Ok::<_, PipelineError>((
                nt_read_dict,
                nt_refined_counts,
                nt_contig_summary,
                nt_pairs_merge,
                nt_pairs_map,
                nt_pairs_taxid_refined,
                nt_pairs_coverage,
                nt_refined_m8_top_stream_out,
                cleanup_tasks,
                cleanup_receivers,
                temp_files,
            ))
        }
    });

    // ── NR post: fanout refined ReducedRead pairs ─────────────────────────────
    // Branches:
    //   merge → compute_merged_taxon_counts
    //   map   → collect_m8_accession_map (pair.dedup)
    //   taxid → generate_taxid_fasta (pair.summary)
    // top m8 is unused on NR — drain so the unbounded sender cannot stall
    let nr_post_handle = tokio::spawn({
        let config = config.clone();
        async move {
            let nr_res = nr_handle.await.map_err(|e| {
                PipelineError::Other(anyhow!("NR blast_contigs task panicked: {e}"))
            })??;

            let (
                nr_read_dict,
                nr_refined_counts,
                nr_contig_summary,
                nr_refined_pairs_out,
                nr_refined_m8_top_stream_out,
                nr_cleanup_tasks,
                nr_cleanup_receivers,
                nr_temp_files,
            ) = nr_res;

            let mut cleanup_tasks = nr_cleanup_tasks;
            let mut cleanup_receivers = nr_cleanup_receivers;
            let temp_files = nr_temp_files;

            tokio::spawn(async move {
                let mut rx = UnboundedReceiverStream::new(nr_refined_m8_top_stream_out);
                let mut n = 0u64;
                while rx.next().await.is_some() {
                    n += 1;
                }
                info!("[nr] drained {} top m8 lines (unused)", n);
            });

            let (nr_pair_rxs, nr_pairs_done_rx) = fanout_to_channels(
                ReceiverStream::new(nr_refined_pairs_out),
                3,
                "nr_refined_pairs",
                &config,
                StreamDataType::JustBytes,
            )
                .await?;
            cleanup_receivers.push(nr_pairs_done_rx);

            let mut it = nr_pair_rxs.into_iter();
            let nr_pairs_merge = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
            let nr_pairs_map = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);
            let nr_pairs_taxid_refined = ReceiverStream::new(it.next().ok_or(PipelineError::EmptyStream)?);

            Ok::<_, PipelineError>((
                nr_read_dict,
                nr_refined_counts,
                nr_contig_summary,
                nr_pairs_merge,
                nr_pairs_map,
                nr_pairs_taxid_refined,
                cleanup_tasks,
                cleanup_receivers,
                temp_files,
            ))
        }
    });

    let (nt_post_res, nr_post_res) = tokio::join!(nt_post_handle, nr_post_handle);

    let nt_post_res = nt_post_res.map_err(|e| {
        PipelineError::Other(anyhow!("NT blast_contigs postprocess task panicked: {e}"))
    })??;
    let nr_post_res = nr_post_res.map_err(|e| {
        PipelineError::Other(anyhow!("NR blast_contigs postprocess task panicked: {e}"))
    })??;


    let (
        _nt_read_dict,
        nt_refined_counts,
        nt_contig_summary,
        nt_pairs_merge,
        nt_pairs_map,
        nt_pairs_taxid_refined,
        nt_pairs_coverage,
        nt_refined_m8_top_stream_out,
        nt_cleanup_tasks,
        nt_cleanup_receivers,
        nt_blast_temp_dir,
    ) = nt_post_res;
    cleanup_tasks.extend(nt_cleanup_tasks);
    if let Some(td) = nt_blast_temp_dir {
        final_temp_dirs.push(td);
    }

    let (
        _nr_read_dict,
        nr_refined_counts,
        nr_contig_summary,
        nr_pairs_merge,
        nr_pairs_map,
        nr_pairs_taxid_refined,
        nr_cleanup_tasks,
        nr_cleanup_receivers,
        nr_blast_temp_dir,
    ) = nr_post_res;
    cleanup_tasks.extend(nr_cleanup_tasks);
    if let Some(td) = nr_blast_temp_dir {
        final_temp_dirs.push(td);
    }

    cleanup_receivers.extend(nt_cleanup_receivers);
    cleanup_receivers.extend(nr_cleanup_receivers);

    for rx in nt_call_cleanup_receivers {
        rx.await??;
    }
    for rx in nr_call_cleanup_receivers {
        rx.await??;
    }

    // Merge on refined pairs (m8 + hit summary stay locked in ReducedRead)
    let merged_fut = compute_merged_taxon_counts(
        config.clone(),
        nt_pairs_merge,
        nt_contig_summary,
        nr_pairs_merge,
        nr_contig_summary,
        lineage_map.clone(),
        should_keep_filter.clone(),
        duplicate_clusters.clone(),
        out_dir.join("refined.m8"),
        out_dir.join("refined.hitsummary.tab"),
        out_dir.join("merged_taxon_counts_with_dcr.json"),
        out_dir.join("assembly").join("combined_contig_summary.json"),
    );

    let ((merged_taxon_counts,_final_contigs, merged_cleanup_tasks), nt_refined_map, nr_refined_map) =
        tokio::try_join!(
        merged_fut,
        async {
            collect_m8_accession_map(nt_pairs_map)
                .await
                .map_err(PipelineError::Other)
        },
        async {
            collect_m8_accession_map(nr_pairs_map)
                .await
                .map_err(PipelineError::Other)
        },
    )?;
    cleanup_tasks.extend(merged_cleanup_tasks);

    let refined_combined_path =
        out_dir
            .join("assembly")
            .join("refined_taxon_counts_with_dcr.json");

    let (_path, refined_write_json_task) = combine_taxon_counts(
        &nt_refined_counts,
        &nr_refined_counts,
        &merged_taxon_counts,
        refined_combined_path,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts failed: {e}")))?;
    cleanup_tasks.push(refined_write_json_task);

    let assembly_dir = out_dir.join("assembly");
    tokio::fs::create_dir_all(&assembly_dir)
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create assembly dir: {}", e)))?;

    // Non-host READS (WDL host_filter), not contigs
    let nonhost_path = out_dir.join("nonhost_R1.fastq");
    let nonhost_file = tokio::fs::File::open(&nonhost_path)
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!("Failed to open nonhost_R1.fastq: {}", e))
        })?;

    let nonhost_rx = parse_bytes::<TokioFile>(nonhost_file, &config, StreamDataType::JustBytes)
        .await
        .map_err(|e| PipelineError::Other(anyhow!("parse_bytes nonhost failed: {}", e)))?;

    let nonhost_stream = ReceiverStream::new(nonhost_rx);

    let annot_concurrency =
        compute_phase_concurrency(&config, "refined_annotate", 0.4, 4.0, 128, 8);

    let (
        mapped_contigs_rx,
        unidentified_contigs_rx,
        unique_unidentified_rx,
        mut annot_tasks,
        mut annot_rxs,
    ) = generate_annotated_fasta(
        config.clone(),
        nonhost_stream,
        duplicate_clusters.clone(),
        nt_refined_map,
        nr_refined_map,
        annot_concurrency,
    )
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
        "Refined generate_annotated_fasta failed: {}",
        e
    ))
        })?;

    cleanup_tasks.append(&mut annot_tasks);
    cleanup_receivers.append(&mut annot_rxs);

    // Refined mapped FASTA:
    //   branch 1 -> generate_taxid_fasta
    //   branch 2 -> required WDL artifact
    let (refined_mapped_rxs, refined_mapped_done_rx) =
        fanout_to_channels(
            ReceiverStream::new(mapped_contigs_rx),
            2,
            "refined_annotated",
            &config,
            StreamDataType::JustBytes,
        )
            .await?;

    cleanup_receivers.push(refined_mapped_done_rx);

    let mut refined_mapped_rxs = refined_mapped_rxs.into_iter();

    let refined_mapped_taxid_stream =
        refined_mapped_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?;

    let refined_mapped_file_stream =
        refined_mapped_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?;


    // Refined unidentified FASTA:
    //   branch 1 -> generate_taxid_fasta
    //   branch 2 -> required WDL artifact
    let (refined_unidentified_rxs, refined_unidentified_done_rx) =
        fanout_to_channels(
            ReceiverStream::new(unidentified_contigs_rx),
            2,
            "refined_unidentified",
            &config,
            StreamDataType::JustBytes,
        )
            .await?;

    cleanup_receivers.push(refined_unidentified_done_rx);

    let mut refined_unidentified_rxs =
        refined_unidentified_rxs.into_iter();

    let refined_unidentified_taxid_stream =
        refined_unidentified_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?;

    let refined_unidentified_file_stream =
        refined_unidentified_rxs
            .next()
            .ok_or(PipelineError::EmptyStream)?;


    // Start artifact consumers immediately.
    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(refined_mapped_file_stream),
        assembly_dir.join("refined_annotated_merged.fa"),
        config.clone(),
        StreamDataType::JustBytes,
        "refined_annotated_merged",
    ));

    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(refined_unidentified_file_stream),
        assembly_dir.join("refined_unidentified.fa"),
        config.clone(),
        StreamDataType::JustBytes,
        "refined_unidentified",
    ));


    // Existing internal artifact remains unchanged.
    let unique_path =
        assembly_dir.join("unique_refined_unidentified.fa");

    let unique_write = write_fasta_stream_to_file(
        ReceiverStream::new(unique_unidentified_rx),
        unique_path,
        config.clone(),
        StreamDataType::JustBytes,
        "unique_refined_unidentified",
    );

    cleanup_tasks.push(unique_write);


    // Hit maps from refined pairs (summary side only). Already ReceiverStream<ReducedRead>.
    let (taxid_mapped_rx, taxid_combined_rx, load_nt_task, load_nr_task, taxid_main_task) =
        generate_taxid_fasta(
            config.clone(),
            ReceiverStream::new(refined_mapped_taxid_stream),
            ReceiverStream::new(refined_unidentified_taxid_stream),
            nt_pairs_taxid_refined,
            nr_pairs_taxid_refined,
            lineage_map.clone(),
        )
            .await
            .map_err(|e| PipelineError::Other(anyhow!("generate_taxid_fasta failed: {}", e)))?;

    cleanup_tasks.push(load_nt_task);
    cleanup_tasks.push(load_nr_task);
    cleanup_tasks.push(taxid_main_task);


    let (taxid_mapped_rxs, taxid_mapped_done_rx) = fanout_to_channels(
        ReceiverStream::new(taxid_mapped_rx),
        2,
        "taxid_mapped",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;
    cleanup_receivers.push(taxid_mapped_done_rx);

    let mut taxid_mapped_iter = taxid_mapped_rxs.into_iter();
    let taxid_mapped_file = taxid_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;
    let taxid_mapped_locator = taxid_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;

    let mapped_path =
        assembly_dir.join("refined_taxid_annot_mapped_only.fasta");

    let combined_path =
        assembly_dir.join("refined_taxid_annot.fasta");

    let write_mapped_handle = write_fasta_stream_to_file(
        ReceiverStream::new(taxid_mapped_file),
        mapped_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "taxid_mapped_file_stream_to_write",
    );
    cleanup_tasks.push(write_mapped_handle);

    let write_combined_handle = write_fasta_stream_to_file(
        ReceiverStream::new(taxid_combined_rx),
        combined_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "taxid_combined_file_stream_to_write",
    );
    cleanup_tasks.push(write_combined_handle);

    let assembly_dir = out_dir.join("assembly");
    info!("Starting generate_taxid_locator on combined FASTA stream");
    let (locator_outputs, mut locator_tasks, _locator_receivers) = generate_taxid_locator(
        config.clone(),
        taxid_mapped_locator, // directly pass the receiver
        assembly_dir,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("generate_taxid_locator failed: {}", e)))?;
    info!("generate_taxid_locator finished");
    cleanup_tasks.append(&mut locator_tasks);
    info!("Taxid locator files generated: {:?}", locator_outputs);

    // *******************
    // Experimental
    // *******************

    let (
        initial_taxid_mapped_rx,
        initial_taxid_combined_rx,
        initial_load_nt_task,
        initial_load_nr_task,
        initial_taxid_main_task,
    ) = experimental_taxid_handle
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
            "Experimental generate_taxid_fasta task panicked: {e}"
        ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
            "Experimental generate_taxid_fasta failed: {e}"
        ))
        })?;

    cleanup_tasks.push(initial_load_nt_task);
    cleanup_tasks.push(initial_load_nr_task);
    cleanup_tasks.push(initial_taxid_main_task);

    let initial_mapped_path = out_dir.join("taxid_annot_mapped_only.fasta");
    let initial_combined_path = out_dir.join("taxid_annot.fasta");

    info!("[run] starting fanout_to_channels for initial taxid FASTA (mapped + combined)");

    let (initial_mapped_rxs, initial_mapped_done_rx) = fanout_to_channels(
        ReceiverStream::new(initial_taxid_mapped_rx),
        3,
        "initial_taxid_mapped",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;
    cleanup_receivers.push(initial_mapped_done_rx);

    let mut initial_mapped_iter = initial_mapped_rxs.into_iter();
    let initial_taxid_mapped_file = initial_mapped_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_mapped_locator = initial_mapped_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_mapped_count = initial_mapped_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    let (initial_combined_rxs, initial_combined_done_rx) = fanout_to_channels(
        ReceiverStream::new(initial_taxid_combined_rx),
        2,
        "initial_taxid_combined",
        &config,
        StreamDataType::JustBytes,
    )
        .await?;
    cleanup_receivers.push(initial_combined_done_rx);

    let mut initial_combined_iter = initial_combined_rxs.into_iter();
    let initial_taxid_combined_file = initial_combined_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_combined_count = initial_combined_iter
        .next()
        .ok_or(PipelineError::EmptyStream)?;

    let initial_write_mapped_handle = write_fasta_stream_to_file(
        ReceiverStream::new(initial_taxid_mapped_file), // Clone rx if needed for logging
        initial_mapped_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_taxid_mapped_file_stream_to_write",
    );
    cleanup_tasks.push(initial_write_mapped_handle);

    let initial_write_combined_handle = write_fasta_stream_to_file(
        ReceiverStream::new(initial_taxid_combined_file),
        initial_combined_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_taxid_combined_file_stream_to_write",
    );
    cleanup_tasks.push(initial_write_combined_handle);

    let initial_mapped_count = stream_record_counter(initial_taxid_mapped_count, false)
        .await?
        .reads;

    let initial_combined_count = stream_record_counter(initial_taxid_combined_count, false)
        .await?
        .reads;

    info!(
    "Experimental: {} mapped, {} combined records",
    initial_mapped_count,
    initial_combined_count
);

    info!("Starting generate_taxid_locator on second");
    let (initial_locator_outputs, mut initial_locator_tasks, _initial_locator_receivers) =
        generate_taxid_locator(
            config.clone(),
            initial_taxid_mapped_locator,
            out_dir.clone(),
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!("Experimental generate_taxid_locator failed: {}", e))
            })?;
    info!("second generate_taxid_locator odone ");

    // Now safe to await the annotated and unidentified fanout tasks
    annotated_done_rx.await??;
    unidentified_done_rx.await??;

    cleanup_tasks.append(&mut initial_locator_tasks);
    info!(
        "Experimental taxid locator files generated: {:?}",
        initial_locator_outputs
    );

    let nt_info_db_path = config
        .args
        .nt_info_tab
        .clone()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("nt_info.tab"));

    let skip_coverage_viz = {
        let mut skip = false;

        if !coverage_json_path.exists() {
            warn!("Skipping coverage viz: contig_coverage_json missing (likely due to SPAdes assembly failure)");
            skip = true;
        } else if coverage_json_path
            .metadata()
            .map(|m| m.len() == 0)
            .unwrap_or(true)
        {
            warn!("Skipping coverage viz: contig_coverage_json is empty");
            skip = true;
        }

        if !assembly_outputs.contig_stats_json.exists() {
            warn!("Skipping coverage viz: contig_stats_json missing");
            skip = true;
        } else if assembly_outputs
            .contig_stats_json
            .metadata()
            .map(|m| m.len() == 0)
            .unwrap_or(true)
        {
            warn!("Skipping coverage viz: contig_stats_json is empty");
            skip = true;
        }

        skip
    };

    let out_dir_for_viz = out_dir.clone();
    let gsnap_m8_path_for_viz = out_dir.join("gsnap.m8");

    let coverage_task = if !skip_coverage_viz {
        info!("All coverage viz inputs present — spawning visualization task");

        // Coverage viz requires the completed original gsnap.m8 artifact.
        // Synchronize on its dedicated writer only here, at the point where
        // the file becomes a real computational dependency.
        nt_m8_write_task
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "gsnap.m8 writer task panicked before coverage viz: {}",
                e
            ))
            })?
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "gsnap.m8 writer failed before coverage viz: {}",
                e
            ))
            })?;

        let gsnap_m8_file = TokioFile::open(&gsnap_m8_path_for_viz)
            .await
            .map_err(|e| {
                PipelineError::IOError(format!(
                    "Failed to open completed gsnap.m8 for coverage viz: {}",
                    e
                ))
            })?;

        let gsnap_m8_rx = parse_lines(
            gsnap_m8_file,
            &config,
            StreamDataType::JustBytes,
        )
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Failed to stream gsnap.m8 for coverage viz: {}",
                e
            ))
            })?;

        // Hit summary bytes from refined pairs
        let hit_for_viz = {
            let mut pairs = nt_pairs_coverage;
            let (tx, rx) = mpsc::channel(65_536);

            tokio::spawn(async move {
                while let Some(p) = pairs.next().await {
                    if tx
                        .send(ParseOutput::Bytes(Bytes::from(p.summary)))
                        .await
                        .is_err()
                    {
                        break;
                    }
                }
            });

            rx
        };

        Some(tokio::spawn(async move {
            generate_coverage_viz(
                ReceiverStream::new(hit_for_viz),
                UnboundedReceiverStream::new(nt_refined_m8_top_stream_out), // pure contig top m8
                coverage_json_path,
                assembly_outputs.contig_stats_json,
                assembly_outputs.contigs_fasta,
                ReceiverStream::new(gsnap_m8_rx),
                nt_info_db_path,
                out_dir_for_viz,
                Some(MAX_NUM_BINS_COVERAGE),
                Some(NUM_ACCESSIONS_PER_TAXON),
                Some(MIN_CONTIG_SIZE),
                false,
            )
                .await
        }))
    } else {
        // Must still drain so fanout does not stall
        tokio::spawn(async move {
            let mut a = nt_pairs_coverage;
            let mut c = UnboundedReceiverStream::new(nt_refined_m8_top_stream_out);
            while a.next().await.is_some() {}
            while c.next().await.is_some() {}
        });
        warn!("Skipping coverage viz: required inputs missing/empty (likely assembly failure)");
        cleanup_tasks.push(nt_m8_write_task); // if this isnt needed, push it onto the cleanup pile
        None
    };

    // *******************
    // Results retrieval
    // *******************

    // Clone paths needed in multiple concurrent async blocks.
    let out_dir_host = out_dir.clone();
    let out_dir_kallisto = out_dir.clone();

    let validation_and_qc_fut = async move {

        let stats = join_with_error_handling(val_count_task).await?;
        info!(
            "Processed {} raw reads, {} raw bases, {} validated reads, {} validated bases, {} undersized, {} oversized reads",
            stats.raw_reads, stats.raw_bases, stats.validated, stats.validated_bases, stats.undersized, stats.oversized
        );

        let _qc_fastp_read_count = match qc_count_result_rx.await {
            Ok(Ok(count)) => {
                info!("Received fastp read count: {}", count);
                count
            }
            Ok(Err(e)) => {
                error!("Error getting fastp read count: {}", e);
                0
            }
            Err(e) => {
                error!("Count receiver dropped: {}", e);
                0
            }
        };

        let ercc_mapped_count = ercc_count_rx
            .await
            .map_err(|e| PipelineError::Other(anyhow!("ERCC count receiver failed: {}", e)))?;

        Ok::<_, PipelineError>((stats, ercc_mapped_count))
    };

    let host_bam_and_counts_fut = async move {
        ercc_bt2_bam_write_handle.await??;
        host_bt2_bam_write_handle.await??;

        if let Some(handle) = optional_human_bam_write_handle {
            handle.await??;
        }

        let host_bt2_counts = host_bt2_count_rx
            .await
            .map_err(|e| PipelineError::Other(anyhow!("Host bt2 counts receiver failed: {}", e)))?;
        info!("Host bt2: mapped counts: {:?}", host_bt2_counts);

        let host_hisat2_counts = host_hisat2_count_rx.await.map_err(|e| {
            PipelineError::Other(anyhow!("Host hisat2 counts receiver failed: {}", e))
        })?;
        info!("Host hisat2: mapped counts: {:?}", host_hisat2_counts);

        if paired {
            info!(
                "Computing host insert size statistics from {}",
                host_bt2_bam_path.display()
            );

            let stats = compute_insert_size_stats_from_bam(host_bt2_bam_path.clone(), None)
                .await
                .context("Failed to compute insert size stats")?;

            let metrics_path = out_dir_host.join("host_insert_metrics.txt");
            let mut f = std::fs::File::create(&metrics_path)
                .context("Failed to create insert metrics file")?;

            writeln!(f, "## METRICS")?;
            writeln!(f, "MEAN_INSERT_SIZE\t{:.2}", stats.mean)?;
            writeln!(f, "MEDIAN_INSERT_SIZE\t{:.0}", stats.median)?;
            writeln!(f, "STANDARD_DEVIATION\t{:.2}", stats.stddev)?;
            writeln!(f, "TOTAL_PROPER_PAIRS\t{}", stats.total_proper_pairs)?;
            writeln!(f, "\n## HISTOGRAM\tinsert_size\tcount")?;
            for &(size, count) in &stats.insert_sizes {
                writeln!(f, "{}\t{}", size, count)?;
            }

            let plot_path = out_dir_host.join("host_insert_size_histogram.png");
            if let Err(e) = plot_insert_sizes(&stats.insert_sizes, sample_base.as_str(), &plot_path)
            {
                warn!("Failed to plot insert size histogram: {}", e);
            }

            info!(
                "Host insert size stats written: mean {:.2}, median {:.0}, {} proper pairs",
                stats.mean, stats.median, stats.total_proper_pairs
            );
        }

        Ok::<_, PipelineError>((host_bt2_counts, host_hisat2_counts))
    };

    let kallisto_and_taxon_fut = async move {
        let kallisto_exit = join_with_error_handling(kallisto_exit_task).await;
        match kallisto_exit {
            Ok(_) => info!("Kallisto completed successfully"),
            Err(e) => {
                warn!(
                    "Kallisto failed (non-fatal, possibly no ERCC spiked in): {}",
                    e
                );
            }
        }

        let kallisto_results_task =
            kallisto_results(out_dir_kallisto.join("kallisto"), kallisto_ercc_tx)
                .await
                .map_err(|e| PipelineError::Other(e.into()))?;

        join_with_error_handling(kallisto_results_task)
            .await
            .map_err(|e| PipelineError::Other(e.into()))?;

        let kallisto_ercc_counts = kallisto_ercc_rx
            .await
            .map_err(|e| PipelineError::Other(anyhow!("ERCC counts receiver failed: {}", e)))?;

        let kallisto_dir = out_dir_kallisto.join("kallisto");
        tokio::fs::create_dir_all(&kallisto_dir).await?;

        let abundance_path = kallisto_dir.join("abundance.tsv");
        let abundance_file = TokioFile::create(&abundance_path).await?;
        let mut abundance_writer = BufWriter::new(abundance_file);

        abundance_writer
            .write_all(b"target_id\tlength\teff_length\test_counts\ttpm\n")
            .await?;

        for (target_id, est_counts) in &kallisto_ercc_counts.ercc_counts {
            let line = format!("{}\t0\t0\t{:.6}\t0.000000\n", target_id, est_counts);
            abundance_writer.write_all(line.as_bytes()).await?;
        }
        abundance_writer.flush().await?;
        info!(
            "Wrote Kallisto abundance.tsv with {} ERCC transcripts",
            kallisto_ercc_counts.ercc_counts.len()
        );

        let ercc_path = out_dir_kallisto.join("kallisto_ERCC_counts_tsv");
        let ercc_file = TokioFile::create(&ercc_path).await?;
        let mut ercc_writer = BufWriter::new(ercc_file);

        ercc_writer.write_all(b"target_id\test_counts\n").await?;

        for (target_id, est_counts) in &kallisto_ercc_counts.ercc_counts {
            let line = format!("{}\t{:.6}\n", target_id, est_counts);
            ercc_writer.write_all(line.as_bytes()).await?;
        }
        ercc_writer.flush().await?;
        info!("Wrote kallisto_ERCC_counts_tsv");

        Ok::<_, PipelineError>((kallisto_ercc_counts, ()))
    };

    let dedup_and_subsample_counts_fut = async move {
        let uniques_count = dedup_count_rx
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Dedup count receiver failed: {}",
                e
            ))
            })?;

        let unique_reads =
            uniques_count * if paired { 2 } else { 1 };

        info!(
        "Uniques count after dedup (pairs counted separately): {}",
        unique_reads
    );

        let subsample_count = subsample_count_rx
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                "Subsample count receiver failed: {}",
                e
            ))
            })?;

        let subsample_reads =
            subsample_count * if paired { 2 } else { 1 };

        info!(
        "Subsampled reads (pairs counted separately): {}",
        subsample_reads
    );

        Ok::<_, PipelineError>((
            unique_reads,
            subsample_reads,
        ))
    };

    let (
        _validation_and_qc,
        _host_bam_and_counts,
        _kallisto_and_taxon,
        _dedup_and_subsample_counts,
    ) = tokio::try_join!(
    validation_and_qc_fut,
    host_bam_and_counts_fut,
    kallisto_and_taxon_fut,
    dedup_and_subsample_counts_fut,
)?;

    // *******************
    // Cleanup
    // *******************

    // Helper to handle downcast + NotFound tolerance in one place
    fn handle_cleanup_error(e: anyhow::Error, context: &str) -> Result<(), PipelineError> {
        match e.downcast::<std::io::Error>() {
            Ok(io_err) => {
                if io_err.kind() == std::io::ErrorKind::NotFound {
                    warn!(
                        "Cleanup ignored missing file/dir (non-fatal) [{}]: {}",
                        context, io_err
                    );
                    Ok(())
                } else {
                    // Other io errors are fatal
                    Err(PipelineError::Other(anyhow::Error::new(io_err)))
                }
            }
            Err(remaining) => {
                // Not an io::Error → treat as fatal
                Err(PipelineError::Other(remaining))
            }
        }
    }

    // Wait for all cleanup tasks
    let cleanup_results = try_join_all(cleanup_tasks)
        .await
        .map_err(|join_err| PipelineError::Other(join_err.into()))?;

    // Process each task result
    for res in cleanup_results {
        if let Err(e) = res {
            handle_cleanup_error(e, "task")?;
        }
    }

    // Process cleanup receivers (oneshot channels)
    for rx in cleanup_receivers {
        match rx.await {
            Ok(Ok(_)) => {} // success

            Ok(Err(e)) => {
                handle_cleanup_error(e, "receiver")?;
            }

            Err(oneshot_err) => {
                // Channel broken (panic, drop, empty stream, etc.)
                // Usually safe to ignore if stream was empty
                warn!(
                    "Cleanup receiver channel failed (possibly empty): {}",
                    oneshot_err
                );
                // continue (non-fatal)
            }
        }
    }

    // Coverage viz — keep non-fatal as before
    if let Some(coverage_task) = coverage_task {
        match coverage_task.await {
            Ok(Ok(_)) => info!("Coverage visualization completed successfully"),
            Ok(Err(e)) => warn!("Coverage viz failed (non-fatal): {}", e),
            Err(join_err) => warn!("Coverage viz task join failed (non-fatal): {}", join_err),
        }
    }

    // Final temp resource drops — defensive, never fail the pipeline
    for tf in temp_files {
        if let Err(e) = tf.close() {
            if e.kind() == std::io::ErrorKind::NotFound {
                debug!("Temp file already gone during close: {}", e);
            } else {
                warn!("Temp file close failed (non-fatal): {}", e);
            }
        }
    }

    for td in final_temp_dirs {
        if let Err(e) = std::fs::remove_dir_all(td.path()) {
            if e.kind() == std::io::ErrorKind::NotFound {
                debug!("Temp dir already gone: {}", e);
            } else {
                warn!("Temp dir removal failed (non-fatal): {}", e);
            }
        }
    }

    drop(temp_paths);

    info!("Finished short read mNGS pipeline (cleanup completed, some non-fatal warnings may have occurred).");
    Ok(())
}
