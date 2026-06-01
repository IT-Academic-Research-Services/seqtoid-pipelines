use std::cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader};
use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, AtomicU64, Ordering as AtomicOrdering};
use std::sync::{Arc, Mutex};

use ahash::AHashMap;
use ahash::RandomState as AHashRandomState;
use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use dashmap::DashMap;
use rand::seq::SliceRandom;
use fst::Map;
use futures::future::try_join_all;
use log::{self, debug, error, info, warn};
use needletail::parse_fastx_file;
use noodles::bam::r#async::io::Reader as BamAsyncReader;
use noodles::bam::record::Record;
use noodles::sam::alignment::record::cigar::{op::Kind as OpKind};
use rand::prelude::*;
use rand_core::{OsRng, RngCore};

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use tempfile::NamedTempFile;
use tempfile::TempDir;
use tokio::fs;
use tokio::fs::{File as TokioFile};
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader as TokioBufReader, BufWriter};
use tokio::process::Command;
use tokio::sync::mpsc;
use tokio::sync::mpsc::{Sender};
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio::process::Child;
use tokio::time::{sleep, Duration, Instant};
use tokio::try_join;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::{StreamExt};
use tokio_util::io::StreamReader;
use twox_hash::XxHash64;
use rayon::prelude::*;


use crate::config::defs::{DiamondSubcommand, KallistoSubcommand, Lineage, PipelineError, ReadCountingMode, ReadStats, RunConfig, SamtoolsSubcommand, StreamDataType, Taxid, BLASTN_TAG, BLASTX_TAG, BOWTIE2_TAG,  DIAMOND_TAG, FASTP_TAG, HISAT2_TAG, KALLISTO_TAG, MAKEBLASTDB_TAG, MINIMAP2_TAG, NR_TAG, NT_TAG, READ_COUNTING_MODE, SAMTOOLS_TAG, SPADES_TAG, MMSEQS_TAG, ClusterInfo, DuplicateClusters, PairingMode, NRAlignmentBackend, SORT_TAG};
use crate::utils::blast::{ generate_taxon_count_json_from_m8, AggBucket, M8Record,
                          TaxonCount, ContigSummaryEntry, SpeciesAlignmentResults, read_id_from_m8_line, shard_for_read_id, summarize_m8_hits, process_record_pair, merge_aggregations, build_taxon_counts_list};
use crate::utils::command::blastn::BlastnConfig;
use crate::utils::command::blastx::BlastxConfig;
use crate::utils::command::bowtie2::{bowtie2_index_prep, Bowtie2Config};
use crate::utils::command::diamond::{diamond_index_prep, DiamondConfig, compute_optimal_block_size};
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::hisat2::{hisat2_index_prep, Hisat2Config};
use crate::utils::command::kallisto::KallistoConfig;
use crate::utils::command::makeblastdb::MakeblastdbConfig;
use crate::utils::command::minimap2::{minimap2_index_prep, Minimap2Config};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::mmseqs::{MmseqsBackend, MmseqsConfig, MmseqsSubcommand};
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::coverage_viz::generate_coverage_viz;
use crate::utils::fastx::{raw_read_count, read_fasta,
                          read_fastq, stream_record_counter, write_fasta_stream_to_file,
                          SequenceRecord, generate_taxid_fasta, generate_taxid_locator,
                          filter_fastq_to_bytes_stream, parse_byte_stream_to_fastq, write_combined_fastq};
use crate::utils::file::{choose_temp_dir, file_path_manipulator,
                         file_size, rename_file_path, resolve_optional_path,
                         validate_file_inputs, write_byte_stream_to_file, write_parse_output_to_file,
                         write_vecu8_to_file, materialize_to_temp};
use crate::utils::paf::{parse_paf_batch_to_m8};
use crate::utils::plotting::plot_insert_sizes;
use crate::utils::sambam::{generate_info_from_bam_stream, compute_insert_size_stats_from_bam};
use crate::utils::streams::{deinterleave_fastq_stream, join_with_error_handling,
                            parse_child_output, parse_fastq, read_child_output_to_vec, spawn_cmd,
                            stream_to_cmd, stream_to_file,
                            ChannelReader, ChildStream, ParseMode,
                            ParseOutput, ToBytes, monitor_stream, batch_rayon_process, parse_bytes, parse_lines, fanout_to_channels};
use crate::utils::streams::deinterleave_fastq_stream_to_fifos;
use crate::utils::system::{compute_phase_concurrency, compute_batch_size};
use crate::utils::taxonomy::{build_should_keep_filter, get_top_m8_nr,
                             get_top_m8_nt, load_taxid_lineages_db};

const MAX_SPADES_WORK_DIR: u64 = 500_000_000;

const MAX_NUM_BINS_COVERAGE: usize = 500;
const NUM_ACCESSIONS_PER_TAXON: usize = 10;
const MIN_CONTIG_SIZE: u64 = 500;

const MIN_REF_FASTA_SIZE: u64 = 25;
const MIN_ASSEMBLED_CONTIG_SIZE: u64 = 25;


#[derive(Debug)]
pub struct KallistoResults {
    pub ercc_counts: Vec<(String, f64)>, // (target_id, est_counts)
    pub transcript_to_gene: Vec<(String, String)>, // (transcript_id, gene_id)
}

// SampleItem for reservoir sampling
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct SampleItem {
    key: f64, // For A-Res sampling
    records: Vec<SequenceRecord>, // 1 for SE, 2 for PE
    weight: u64, // Duplicate count
    hash_key: Option<u64>, // For disk mode
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

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CombinedTaxonCount {
    pub tax_id: i32,
    pub tax_level: u8,
    pub genus_taxid: i32,
    pub family_taxid: i32,
    pub nt: Option<TaxonMetrics>,
    pub nr: Option<TaxonMetrics>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TaxonMetrics {
    pub count: u64,           // Adjusted reads assigned (with DCR)
    pub nonunique_count: u64,
    pub unique_count: u64,
    pub dcr: f64,             // Duplicate compression ratio
    pub percent_identity: f64,
    pub alignment_length: f64,
    pub e_value: f64,
    pub base_count: u64,
    pub source_count_type: Option<Vec<String>>,
}

#[derive(Debug)]
pub struct AssemblyHandle {
    pub spades_task: JoinHandle<Result<(), anyhow::Error>>,
    pub work_dir: PathBuf,
    pub out_dir: PathBuf,
}

#[derive(Debug)]
pub struct CoverageOutputs {
    pub contigs_fasta: PathBuf,
    pub contigs_all_fasta: PathBuf,
    pub scaffolds_fasta: PathBuf,
    pub sam_path: PathBuf,
    pub contig_stats_json: PathBuf,
    pub coverage_json: PathBuf,
    pub coverage_summary_csv: PathBuf,
    pub contig_stats: Arc<HashMap<String, u64>>,
    pub read2contig: Arc<HashMap<String, String>>,
}

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


    pub fn to_full_tab_line(&self, read_id: &str) -> String {
        format!("{}\t{}", read_id, self.to_tab_string())
    }

    pub fn from_tab_line(line: &str) -> Result<(String, Self)> {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 28 {
            return Err(anyhow!("Invalid hit summary line: expected 28 fields, got {}", fields.len()));
        }

        let read_id = fields[0].to_string();
        let hit = Self {
            level: fields[1].parse()?,
            taxid: fields[2].parse()?,
            accession_id: fields[3].to_string(),
            species_taxid: fields[8].parse()?,
            genus_taxid: fields[9].parse()?,
            family_taxid: fields[10].parse()?,
            contig_id: if fields[13].is_empty() { None } else { Some(fields[13].to_string()) },
            contig_accession_id: if fields[14].is_empty() { None } else { Some(fields[14].to_string()) },
            contig_species_taxid: fields[15].parse()?,
            contig_genus_taxid: fields[16].parse()?,
            contig_family_taxid: fields[17].parse()?,
            from_assembly: fields[25] == "1",
            source_count_type: if fields[27].is_empty() { None } else { Some(fields[27].to_string()) },
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


#[allow(dead_code)]
#[derive(Debug, Clone)]
struct Cluster {
    rep_id: String,
    size: u64,
    members: Vec<String>,
}


#[allow(dead_code)]
#[derive(Clone)]
struct WeightedSampleItem {
    key: f64,
    record: SequenceRecord,
    weight: u64,
}

/// Called read_fastq in single or paired FASTQ's and streams interleaved output
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `file1_path` - Path to R1 or single ended FASTQ.
/// * `file2_path` - Optional path to R2.
///  * `sample_base_buf` - PathBuf of sample basename.
/// * `out_dir` - Base dir for output files.
///
/// # Returns
/// ParseOutput validated interleaved FASTQ stream, vecs of cleaup takss and recoevers, handle for raw and validated read counts.
async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> anyhow::Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>, JoinHandle<anyhow::Result<u64, anyhow::Error>>, JoinHandle<anyhow::Result<ReadStats, anyhow::Error>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();


    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated.fq.gz"),
        "_",
    );

    let raw_count_task = raw_read_count(file1_path.clone(), file2_path.clone());


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
    let (val_streams, val_router_handle) = fanout_to_channels(
        val_rx_stream,
        2,
        "validate_input",
        &config,
        StreamDataType::IlluminaFastq
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_tasks.push(val_router_handle);

    let mut val_streams_iter = val_streams.into_iter();

    let val_ercc_bowtie2_filter_out_stream = ReceiverStream::new(
        val_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let val_file_stream = ReceiverStream::new(
        val_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );


    let val_file_write_task = write_byte_stream_to_file(
        &validated_interleaved_file_path,
        val_file_stream,
        config.clone(),
        StreamDataType::IlluminaFastq,
        "validate_input",
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(val_file_write_task);

    Ok((val_ercc_bowtie2_filter_out_stream,  cleanup_tasks, cleanup_receivers, raw_count_task, val_count_task))
}

/// bowtie2 filter function where the passing stream contains the unmapped reads
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
/// * `bt2_index_path` - Path to Bowtie2 index.
/// * `paired` - Whether the input is paired-end.
/// * `bowtie2_options` - Additional Bowtie2 options as a HashMap (e.g., HashMap::from([("--very-sensitive-local".to_string(), None)])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - unmapped FASTQ stream.
/// - Optional receiver for the total mapped count (u64) if `count_mapped` is true.
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.

async fn bowtie2_align_and_sort_stream(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    bowtie2_options: HashMap<String, Option<String>>,
    output_bam_path: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,   // BAM → for fastq extraction
    JoinHandle<Result<()>>,        // BAM write task
    oneshot::Receiver<u64>,        // count rx
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    TempDir
), PipelineError> {

    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    ).await?;

    // ─────────────────────────────
    // Bowtie2
    // ─────────────────────────────

    let bt2_config = Bowtie2Config {
        bt2_index_path,
        paired,
        r1_path: None,
        r2_path: None,
        option_fields: bowtie2_options,
    };

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config))?;

    let (bt2_child, bt2_stream_task, bt2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        BOWTIE2_TAG,
        bt2_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None
    ).await?;

    cleanup_tasks.push(bt2_stream_task);
    cleanup_tasks.push(bt2_err_task);

    let bt2_out = {
        let mut guard = bt2_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config).await?
    };

    // ─────────────────────────────
    // samtools sort -n, needs name sorting before conversion
    // ─────────────────────────────

    let sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),          // name-sorted
            ("-u".to_string(), None),          // uncompressed
            ("-O".to_string(), Some("bam".to_string())),
            ("-T".to_string(), Some(temp_dir.path().to_string_lossy().to_string())),
            ("-".to_string(), None),
        ]),
    };

    let sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&sort_config))?;

    let (sort_child, sort_task, sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out,
        SAMTOOLS_TAG,
        sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None
    ).await?;

    cleanup_tasks.push(sort_task);
    cleanup_tasks.push(sort_err_task);

    let sorted_bam_stream = {
        let mut guard = sort_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, &config).await?
    };

    // ─────────────────────────────
    // split stream
    // ─────────────────────────────

    let (streams, router_handle) = fanout_to_channels(
        ReceiverStream::new(sorted_bam_stream),
        3,
        "bt2_split",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    // track the task instead of a oneshot receiver
    cleanup_tasks.push(router_handle);

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
            ).await?;
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
        temp_dir
    ))
}

async fn bowtie2_filter_stream(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    bowtie2_options: HashMap<String, Option<String>>,
    output_bam_path: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,
    oneshot::Receiver<u64>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    JoinHandle<Result<()>>,
    PathBuf,
    Vec<TempDir>
), PipelineError> {

    let (
        fastq_input_stream,
        bam_write_task,
        count_rx,
        mut cleanup_tasks,
        cleanup_receivers,
        temp_dir
    ) = bowtie2_align_and_sort_stream(
        config.clone(),
        input_stream,
        bt2_index_path,
        paired,
        bowtie2_options,
        output_bam_path.clone(),
    ).await?;

    // samtools fastq → stdout
    let fastq_flag = if paired { "-f13" } else { "-f4" };

    let config_fastq = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (fastq_flag.to_string(), None),
            ("-".to_string(), None),
        ]),
    };

    let args = generate_cli(SAMTOOLS_TAG, &config, Some(&config_fastq))?;

    let (child, stdin_task, err_task) = stream_to_cmd(
        config.clone(),
        fastq_input_stream.into_inner(),
        SAMTOOLS_TAG,
        args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None
    ).await?;

    cleanup_tasks.push(stdin_task);
    cleanup_tasks.push(err_task);

    let out_stream = {
        let mut guard = child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Fastq, &config).await?
    };

    Ok((
        ReceiverStream::new(out_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        bam_write_task,
        output_bam_path,
        vec![temp_dir]
    ))
}


async fn bowtie2_filter_files(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    bowtie2_options: HashMap<String, Option<String>>,
    output_bam_path: PathBuf,
) -> Result<(
    PathBuf,
    Option<PathBuf>,
    oneshot::Receiver<u64>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    JoinHandle<Result<()>>,
    PathBuf,
    Vec<TempDir>
), PipelineError> {

    let (
        fastq_input_stream,
        bam_write_task,
        count_rx,
        mut cleanup_tasks,
        cleanup_receivers,
        temp_dir
    ) = bowtie2_align_and_sort_stream(
        config.clone(),
        input_stream,
        bt2_index_path,
        paired,
        bowtie2_options,
        output_bam_path.clone(),
    ).await?;

    let fastq_temp_dir = choose_temp_dir(
        config.input_size, // FASTQ size ~ input size
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        2,
        false,
    ).await?;

    let fq1 = fastq_temp_dir.path().join("bt2_R1.fastq");
    let fq2 = if paired {
        Some(fastq_temp_dir.path().join("bt2_R2.fastq"))
    } else {
        None
    };

    let mut fields = HashMap::new();
    fields.insert("-f".to_string(), Some(if paired { "13".to_string() } else { "4".to_string() }));

    if paired {
        fields.insert("-1".to_string(), Some(fq1.to_string_lossy().to_string()));
        fields.insert("-2".to_string(), Some(fq2.as_ref().unwrap().to_string_lossy().to_string()));
        fields.insert("-0".to_string(), Some("/dev/null".to_string()));
        fields.insert("-s".to_string(), Some("/dev/null".to_string()));
    } else {
        fields.insert("-o".to_string(), Some(fq1.to_string_lossy().to_string()));
    }

    let config_fastq = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: fields,
    };

    let args = generate_cli(SAMTOOLS_TAG, &config, Some(&config_fastq))?;

    let (child, stdin_task, err_task) = stream_to_cmd(
        config.clone(),
        fastq_input_stream.into_inner(),
        SAMTOOLS_TAG,
        args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None
    ).await?;

    cleanup_tasks.push(stdin_task);
    cleanup_tasks.push(err_task);

    {
        let mut guard = child.lock().await;
        guard.wait().await?;
    }

    Ok((
        fq1,
        fq2,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        bam_write_task,
        output_bam_path,
        vec![temp_dir, fastq_temp_dir]
    ))
}

/// HISAT2 filter
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream.
/// * `hisat2_index_path` - Path to HISAT2 index.
/// * `paired` - Whether the input is paired-end.
/// * `hisat2_options` - Additional HISAT2 options as a HashMap (e.g., HashMap::from([("--no-spliced-alignment".to_string(), None)])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - Unmapped FASTQ stream.
/// - Optional receiver for the total mapped count (u64) if count is needed.
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.

async fn hisat2_filter(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    hisat2_index_path: PathBuf,
    paired: bool,
    hisat2_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
    headroom: u64,
) -> Result<(
    ReceiverStream<ParseOutput>,
    oneshot::Receiver<u64>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    Vec<TempDir>
), PipelineError> {

    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size * headroom,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false
    ).await?;

    // ─────────────────────────────
    // 1. HISAT2 → SAM
    // ─────────────────────────────

    let sam_path = temp_dir.path().join("hisat2.sam");

    let hisat2_config = Hisat2Config {
        hisat2_index_path: hisat2_index_path.clone(),
        option_fields: hisat2_options,
        r1_path: r1_path.to_string_lossy().to_string(),
        r2_path: r2_path_opt.as_ref().map(|p| p.to_string_lossy().to_string()),
    };

    let mut hisat2_args = generate_cli(HISAT2_TAG, &config, Some(&hisat2_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: HISAT2_TAG.to_string(),
            error: e.to_string(),
        })?;

    hisat2_args.push("-S".to_string());
    hisat2_args.push(sam_path.to_string_lossy().to_string());

    let (mut hisat2_child, hisat2_err_task) = spawn_cmd(
        config.clone(),
        HISAT2_TAG,
        hisat2_args,
        config.args.verbose,
        None
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: HISAT2_TAG.to_string(),
        error: e.to_string(),
    })?;

    cleanup_tasks.push(hisat2_err_task);

    let status = hisat2_child.wait().await
        .map_err(|e| PipelineError::Other(e.into()))?;

    if !status.success() {
        return Err(PipelineError::ToolExecution {
            tool: HISAT2_TAG.to_string(),
            error: "HISAT2 execution failed".to_string(),
        });
    }

    // ─────────────────────────────
    // 2. sort -n → BAM
    // ─────────────────────────────

    let bam_path = temp_dir.path().join("hisat2_sorted.bam");

    let sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-o".to_string(), Some(bam_path.to_string_lossy().to_string())),
            ("-@".to_string(), Some("8".to_string())),
            ("-T".to_string(), Some(temp_dir.path().to_string_lossy().to_string())),
        ]),
    };

    let mut sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    sort_args.push(sam_path.to_string_lossy().to_string());

    let (mut sort_child, sort_err_task) = spawn_cmd(
        config.clone(),
        SAMTOOLS_TAG,
        sort_args,
        config.args.verbose,
        None
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;

    cleanup_tasks.push(sort_err_task);

    sort_child.wait().await
        .map_err(|e| PipelineError::Other(e.into()))?;

    if let Some(out) = output_bam_path {
        tokio::fs::copy(&bam_path, &out).await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
    }

    // ─────────────────────────────
    // 3. samtools fastq → R1/R2 files
    // ─────────────────────────────

    let fq1_path = temp_dir.path().join("hisat2_filtered_R1.fastq");

    let fq2_path_opt = if paired {
        Some(temp_dir.path().join("hisat2_filtered_R2.fastq"))
    } else {
        None
    };

    let mut fastq_fields = HashMap::new();
    fastq_fields.insert(
        "-f".to_string(),
        Some(if paired { "13".to_string() } else { "4".to_string() }),
    );

    if paired {
        fastq_fields.insert("-1".to_string(), Some(fq1_path.to_string_lossy().to_string()));
        fastq_fields.insert("-2".to_string(), Some(fq2_path_opt.as_ref().unwrap().to_string_lossy().to_string()));
        fastq_fields.insert("-0".to_string(), Some("/dev/null".to_string()));
        fastq_fields.insert("-s".to_string(), Some("/dev/null".to_string()));
    } else {
        fastq_fields.insert("-o".to_string(), Some(fq1_path.to_string_lossy().to_string()));
    }

    let fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: fastq_fields,
    };

    let mut fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    fastq_args.push(bam_path.to_string_lossy().to_string());

    let (mut fastq_child, fastq_err_task) = spawn_cmd(
        config.clone(),
        SAMTOOLS_TAG,
        fastq_args,
        config.args.verbose,
        None
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;

    cleanup_tasks.push(fastq_err_task);

    fastq_child.wait().await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // ─────────────────────────────
    // 4. convert output FASTQ → stream (for downstream)
    // ─────────────────────────────

    let (rx, read_task) = read_fastq(
        fq1_path.clone(),
        fq2_path_opt.clone(),
        PairingMode::Strict,
        None,
        u64::MAX,
        None,
        None,
        "hisat2_filter",
        &config,
    ).map_err(|e| PipelineError::Other(e))?;

    cleanup_tasks.push(tokio::spawn(async move {
        read_task.await??;
        Ok(())
    }));

    let out_stream = ReceiverStream::new(rx);

    // ─────────────────────────────
    // 5. mapped count
    // ─────────────────────────────

    let (count_tx, count_rx) = oneshot::channel();

    let count_task = tokio::spawn(async move {
        let flag = if paired { "-F13" } else { "-F4" };
        let output = Command::new("samtools")
            .args(["view", "-c", flag, bam_path.to_str().unwrap()])
            .output()
            .await?;

        let count = String::from_utf8_lossy(&output.stdout)
            .trim()
            .parse::<u64>()
            .unwrap_or(0);

        let _ = count_tx.send(count);
        Ok::<(), anyhow::Error>(())
    });

    cleanup_tasks.push(count_task);

    Ok((
        out_stream,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        vec![temp_dir]
    ))
}

/// QC's input stream using FASTP
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
///
/// # Returns
///
async fn fastp_qc(
    config: Arc<RunConfig>,
    paired: bool,
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<(ReceiverStream<ParseOutput>,  Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, oneshot::Receiver<Result<u64, anyhow::Error>>), PipelineError>{
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    let qc_fastp_config_view = FastpConfig {

        //These default QC thresholds are loosely based on the pre-2022 CZI pipeline using PriceSeq & LZW
        command_fields: HashMap::from([
            ("--dont_eval_duplication".to_string(), None),
            ("--length_required".to_string(), Some("35".to_string())),
            ("--qualified_quality_phred".to_string(), Some("17".to_string())),
            ("--unqualified_percent_limit".to_string(), Some("15".to_string())),
            ("--n_base_limit".to_string(), Some("15".to_string())),
            ("--complexity_threshold".to_string(), Some("60".to_string())),
            // ("--sdust_complexity_filter".to_string(), None),
        ]),

        paired: paired

    };

    let qc_fastp_args = generate_cli(FASTP_TAG, &config, Some(&qc_fastp_config_view))
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
        None
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
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Tee the byte stream for counting
    let tee_count_input = ReceiverStream::new(qc_fastp_out_stream);

    let (tee_count_streams, tee_count_router_handle) = fanout_to_channels(
        tee_count_input,
        2,
        "fastp_count",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "fanout_to_channels".to_string(),
            error: e.to_string(),
        })?;

    // track the task instead of oneshot receiver
    cleanup_tasks.push(tee_count_router_handle);

    let mut tee_count_streams_iter = tee_count_streams.into_iter();

    let qc_out_stream = ReceiverStream::new(
        tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let qc_count_stream = ReceiverStream::new(
        tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let (count_tx, count_result_rx) = oneshot::channel::<Result<u64, anyhow::Error>>();

    let count_task = tokio::spawn(async move {
        match stream_record_counter(qc_count_stream.into_inner(), false).await {
            Ok(count) => {
                debug!("Fastp output reads: {}", count);
                let _ = count_tx.send(Ok(count)); // Send the count result
                Ok(())
            }
            Err(e) => {
                let _ = count_tx.send(Err(e));
                Err(anyhow!("Failed to count fastp output reads"))
            },
        }
    });
    cleanup_tasks.push(count_task);

    Ok((qc_out_stream, cleanup_tasks, cleanup_receivers, count_result_rx))
}

/// Runs Kallisto quant
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns
async fn kallisto_quant(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    out_dir: PathBuf,
    paired: bool,
    sample_base_buf: PathBuf,
) -> Result<(
    oneshot::Sender<KallistoResults>,
    oneshot::Receiver<KallistoResults>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    JoinHandle<Result<(), anyhow::Error>>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    // Convert sample_base_buf to &str
    let sample_base = sample_base_buf
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| PipelineError::IOError(format!("Invalid sample base path: {}", sample_base_buf.display())))?
        .to_string();

    // Create Kallisto output dir (idempotent, fast on NVMe)
    let kallisto_out_dir = out_dir.join("kallisto");
    fs::create_dir_all(&kallisto_out_dir).await
        .map_err(|e| PipelineError::IOError(format!("Failed to create Kallisto output directory: {}", e)))?;
    debug!("Kallisto output dir created/confirmed: {}", kallisto_out_dir.display());

    // Convert raw byte stream to FASTQ
    let byte_rx = input_stream.into_inner();
    let byte_reader = ChannelReader::new(byte_rx);
    let fastq_rx = parse_fastq(byte_reader, &config, StreamDataType::IlluminaFastq).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "parse_fastq".to_string(),
            error: e.to_string(),
        })?;
    let fastq_stream = ReceiverStream::new(fastq_rx);

    debug!("Deinterleaving FASTQ stream for sample: {}", sample_base);
    let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        fastq_stream,
        &sample_base,
        paired,
    ).await?;

    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle_opt {
        cleanup_tasks.push(r2_handle);
    }

    let kallisto_config = if paired {
        let r2_fifo = r2_fifo_opt.as_ref().ok_or_else(|| {
            PipelineError::InvalidConfig("Paired mode requested but R2 FIFO not created".to_string())
        })?;
        KallistoConfig {
            subcommand: KallistoSubcommand::Quant,
            subcommand_fields: HashMap::from([
                ("R1".to_string(), Some(r1_fifo.to_string_lossy().to_string())),
                ("R2".to_string(), Some(r2_fifo.to_string_lossy().to_string())),
            ]),
            output_dir: kallisto_out_dir.clone(),
            reproducible: false,
        }
    } else {
        KallistoConfig {
            subcommand: KallistoSubcommand::Quant,
            subcommand_fields: HashMap::from([
                ("--single".to_string(), None),
                ("-l".to_string(), Some("200".to_string())),
                ("-s".to_string(), Some("20".to_string())),
                ("R1".to_string(), Some(r1_fifo.to_string_lossy().to_string())),
            ]),
            output_dir: kallisto_out_dir.clone(),
            reproducible: false,
        }
    };

    let kallisto_args = generate_cli(KALLISTO_TAG, &config, Some(&kallisto_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    debug!("Spawning Kallisto with args: {:?}", kallisto_args);

    let (mut kallisto_child, kallisto_err_task) = spawn_cmd(
        config.clone(),
        KALLISTO_TAG,
        kallisto_args,
        config.args.verbose,
        None
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(kallisto_err_task);

    // Create Kallisto exit task with safe FIFO cleanup
    let r1_fifo_clone = r1_fifo.clone();
    let r2_fifo_opt_clone = r2_fifo_opt.clone();
    let kallisto_exit_task = tokio::spawn(async move {
        let status = kallisto_child.wait().await
            .map_err(|e| anyhow!("Kallisto wait failed: {}", e))?;
        if !status.success() {
            return Err(anyhow!("Kallisto exited with code: {:?}", status.code()));
        }

        // Clean up R1 FIFO
        tokio::fs::remove_file(&r1_fifo_clone).await
            .map_err(|e| anyhow!("Failed to remove R1 FIFO {}: {}", r1_fifo_clone.display(), e))?;

        // Clean up R2 FIFO only if it exists
        if let Some(r2_fifo) = r2_fifo_opt_clone {
            tokio::fs::remove_file(&r2_fifo).await
                .map_err(|e| anyhow!("Failed to remove R2 FIFO {}: {}", r2_fifo.display(), e))?;
        }

        Ok(())
    });

    // Create oneshot channel for results
    let (ercc_tx, ercc_rx) = oneshot::channel::<KallistoResults>();

    Ok((ercc_tx, ercc_rx, cleanup_tasks, cleanup_receivers, kallisto_exit_task))
}

/// Parses results of kallisto
///
/// # Arguments
///
/// * 'kallisto_out_dir': PathBuf to output dir of kallisto
/// *  'ercc_tx': oneshot::Sender<KallistoResults>,
///
/// # Returns
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
                    warn!("Waiting for kallisto abundance.tsv ({} attempts left)...", attempts + 1);
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
                warn!("Malformed kallisto line ({} fields): {}", fields.len(), line);
                skipped += 1;
                continue;
            }

            let target_id = fields[0].to_string();

            let est_counts: f64 = match fields[3].parse() {
                Ok(v) => v,
                Err(e) => {
                    warn!("Invalid est_counts '{}' → {} (line: {})", fields[3], e, line);
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

/// Minimap2 filter
///
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream (paired or single-end).
/// * `paired` - Whether the input is paired-end.
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - Interleaved FASTQ stream : unmapped reads.
/// - Receiver for the total mapped read count
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
/// - Optional temporary FASTA file for host reference.
/// - Optional temporary minimap2 index file.
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
        Vec<TempDir>
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();
    let mut temp_dirs: Vec<TempDir> = Vec::new();
    let mut named_temp_files: Vec<NamedTempFile> = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false
    ).await?;


    let (_host_ref_fasta_path, host_index_path, host_ref_temp, host_index_temp, host_temp_dir,  host_ref_tasks) =
        minimap2_index_prep(
            &config,
            &config.ram_temp_dir,
            config.args.host_sequence.clone(),
            config.args.host_index.clone(),
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

    let minimap2_config = Minimap2Config {
        minimap2_index_path: host_index_path,
        r1_path: None,
        r2_path: None,
        option_fields: HashMap::from([
            ("-ax".to_string(), Some("sr".to_string())), // Short-read preset
            ("-k".to_string(), Some("15".to_string())), // Default k-mer
            ("-w".to_string(), Some("10".to_string())), // Default window
            ("-s".to_string(), Some("30".to_string())), // Match HISAT2's score threshold
            ("-A".to_string(), Some("2".to_string())), // Match score
            ("-B".to_string(), Some("4".to_string())), // Mismatch penalty
            ("-O".to_string(), Some("4,24".to_string())), // Gap open
            ("-E".to_string(), Some("2,1".to_string())), // Gap extension
        ]),
        num_threads: None,
    };
    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&minimap2_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;


    let (minimap2_child, minimap2_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::IlluminaFastq, // Assuming interleaved FASTQ
        config.args.verbose,
        None
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
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-T".to_string(), Some(temp_dir.path().to_string_lossy().to_string())),
            ("-".to_string(), None),
        ]),
    };

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
        None
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
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    use tokio_stream::wrappers::ReceiverStream;

    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_router_handle) = fanout_to_channels(
        bam_rx_stream,
        num_tees,
        "minimap2_bam_split",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    // track task instead of oneshot receiver
    cleanup_tasks.push(bam_router_handle);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = ReceiverStream::new(
            bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?
        );

        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            stream,
            config.clone(),
            StreamDataType::JustBytes,
            "mm2_filter_bam"
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;

        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = ReceiverStream::new(
        bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );
    let mapped_flag = if paired { "-F12".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None),
            (mapped_flag, None),
            ("-".to_string(), None),
        ]),
    };
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
        None
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
            let _count_lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout, &config).await?;
            Ok(())
        }
    });
    cleanup_tasks.push(count_future);

    // Extract unmapped FASTQ
    let unmapped_stream = ReceiverStream::new(
        bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None),
        ]),
    };
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
        None
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
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            &config,
        )
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
        temp_dirs
    ))
}


async fn dedup(
    config: Arc<RunConfig>,
    input_stream: mpsc::Receiver<ParseOutput>,
    paired: bool,
    prefix_len: Option<usize>,
    out_dir: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,  // uniques stream
    oneshot::Receiver<u64>,      // unique count rx
    Arc<DashMap<String, ClusterInfo>>,  // duplicate_clusters
    Vec<JoinHandle<Result<(), anyhow::Error>>>,  // cleanup tasks
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,  // cleanup rx
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    // CPU-bound hashing, moderate memory
    let num_workers = compute_phase_concurrency(
        &config,
        "dedup",
        1.0,  // ~1GB per worker
        1.0,  // 1 thread per core
        config.max_cores,
        4,    // min 4 workers
    );
    let num_shards = (num_workers * 16).max(1024);  // Good contention balance

    let shards: Vec<Arc<Mutex<ahash::AHashMap<u64, (String, u64, Vec<String>)>>>> = (0..num_shards)
        .map(|_| Arc::new(Mutex::new(ahash::AHashMap::new())))
        .collect();

    let total_count = Arc::new(AtomicU64::new(0));
    let unique_count_atomic = Arc::new(AtomicU64::new(0));

    // Channel for writing clusters.csv incrementally
    let (csv_tx, csv_rx) = mpsc::channel::<(String, String)>(config.base_buffer_size * num_workers);
    let csv_path = out_dir.join("clusters.csv");
    let csv_task = tokio::spawn(async move {
        let mut file = BufWriter::new(
            TokioFile::create(&csv_path).await
                .context("Failed to create clusters.csv")?,
        );
        file.write_all(b"representative read id,read id\n").await
            .context("CSV header write failed")?;
        let mut rx = csv_rx;
        while let Some((rep, member)) = rx.recv().await {
            file.write_all(format!("{},{}\n", rep, member).as_bytes()).await
                .context("CSV write failed")?;
        }
        file.flush().await.context("CSV flush failed")?;
        Ok(())
    });
    cleanup_tasks.push(csv_task);

    let (uniques_tx, uniques_rx) = mpsc::channel(config.base_buffer_size * num_workers);

    let mut worker_handles: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::with_capacity(num_workers);

    let worker_txs: Vec<mpsc::Sender<ParseOutput>> = (0..num_workers)
        .map(|_| {
            let (tx, rx) = mpsc::channel(config.base_buffer_size);
            let csv_tx_clone = csv_tx.clone();
            let handle = tokio::spawn(dedup_worker(
                rx,
                shards.clone(),
                uniques_tx.clone(),
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

    // Distributor: round-robin to workers
    let mut stream = ReceiverStream::new(input_stream);
    let mut i = 0usize;
    if paired {
        while let Some(r1_item) = stream.next().await {
            let r1 = match r1_item {
                ParseOutput::Fastq(rec) => rec,
                _ => continue,
            };

            if let Some(r2_item) = stream.next().await {
                let r2 = match r2_item {
                    ParseOutput::Fastq(rec) => rec,
                    _ => continue,
                };

                if r1.id() != r2.id() {
                    return Err(PipelineError::InvalidFastqFormat(format!(
                        "Mismatched pair IDs: R1={}, R2={}",
                        r1.id(), r2.id()
                    )));
                }

                let worker_tx = &worker_txs[i % num_workers];
                worker_tx.send(ParseOutput::Fastq(r1)).await
                    .map_err(|e| PipelineError::Other(anyhow!("Send R1 failed: {}", e)))?;
                worker_tx.send(ParseOutput::Fastq(r2)).await
                    .map_err(|e| PipelineError::Other(anyhow!("Send R2 failed: {}", e)))?;
                i += 1;
            } else {
                return Err(PipelineError::InvalidFastqFormat("Missing R2 in paired stream".to_string()));
            }
        }
    } else {
        while let Some(item) = stream.next().await {
            let worker_tx = &worker_txs[i % num_workers];
            worker_tx.send(item).await
                .map_err(|e| PipelineError::Other(anyhow!("Send failed: {}", e)))?;
            i += 1;
        }
    }

    drop(worker_txs);

    for handle in worker_handles {
        handle.await
            .map_err(|e| PipelineError::Other(anyhow!("Worker join failed: {}", e)))?
            .context("Worker error")?;
    }

    // Drop csv_tx to close csv_rx
    drop(csv_tx);

    // Build duplicate_clusters from shards (rep_id -> ClusterInfo)
    let duplicate_clusters = Arc::new(DashMap::with_capacity(num_shards));
    for shard in &shards {
        let guard = shard.lock()
            .map_err(|e| PipelineError::Other(anyhow!("Shard lock poisoned: {}", e)))?;
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

    // Write duplicate_cluster_sizes.tsv
    let tsv_path = out_dir.join("duplicate_cluster_sizes.tsv");
    let duplicate_clusters_clone = duplicate_clusters.clone();
    let tsv_task = tokio::spawn(async move {
        let mut file = BufWriter::new(
            TokioFile::create(&tsv_path).await
                .context("Failed to create TSV file")?,
        );
        file.write_all(b"representative read id\tcluster size\n").await
            .context("TSV header write failed")?;
        for entry in duplicate_clusters_clone.iter() {
            let (rep_id, info) = entry.pair();
            file.write_all(format!("{}\t{}\n", rep_id, info.size).as_bytes()).await
                .context("TSV write failed")?;
        }
        file.flush().await.context("TSV flush failed")?;
        Ok(())
    });
    cleanup_tasks.push(tsv_task);

    let unique_count = unique_count_atomic.load(AtomicOrdering::Relaxed);
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(unique_count)
        .map_err(|_| PipelineError::Other(anyhow!("Unique count send failed")))?;

    info!("Dedup complete: {} total records, {} unique clusters", total_count.load(AtomicOrdering::Relaxed), unique_count);

    Ok((
        ReceiverStream::new(uniques_rx),
        count_rx,
        duplicate_clusters,
        cleanup_tasks,
        cleanup_receivers,
    ))
}

async fn dedup_worker(
    mut rx: mpsc::Receiver<ParseOutput>,
    shards: Vec<Arc<Mutex<ahash::AHashMap<u64, (String, u64, Vec<String>)>>>>,
    uniques_tx: mpsc::Sender<ParseOutput>,
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
            _ => continue,
        };

        if paired {
            if orphan_r1.is_none() {
                orphan_r1 = Some(record);
                continue;
            }

            let r1 = orphan_r1.take().unwrap();
            let r2 = record;

            if r1.id() != r2.id() {
                return Err(anyhow!("Mismatched pair IDs in worker: {} != {}", r1.id(), r2.id()));
            }

            let id = r1.id().to_string();

            let r1_prefix = &r1.seq()[0..prefix_len.map_or(r1.seq().len(), |l| l.min(r1.seq().len()))];
            let r2_prefix = &r2.seq()[0..prefix_len.map_or(r2.seq().len(), |l| l.min(r2.seq().len()))];

            let mut hasher = XxHash64::default();
            hasher.write(r1_prefix);
            hasher.write_u8(0);
            hasher.write(r2_prefix);
            let hash_key = hasher.finish();

            let shard_idx = (hash_key as usize) % shards.len();

            // ─── Critical: only hold lock for the minimal time ───
            let (rep_id, is_new) = {
                let mut shard = shards[shard_idx].lock()
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
            };  // ← guard dropped here

            total_count.fetch_add(1, AtomicOrdering::Relaxed);

            // Now safe to await
            csv_tx.send((rep_id, id)).await
                .context("CSV send failed")?;

            if is_new {
                uniques_tx.send(ParseOutput::Fastq(r1)).await
                    .context("Send R1 unique failed")?;
                uniques_tx.send(ParseOutput::Fastq(r2)).await
                    .context("Send R2 unique failed")?;
            }
        } else {
            let prefix = &record.seq()[0..prefix_len.map_or(record.seq().len(), |l| l.min(record.seq().len()))];

            let mut hasher = XxHash64::default();
            hasher.write(prefix);
            let hash_key = hasher.finish();

            let shard_idx = (hash_key as usize) % shards.len();

            // Same pattern: minimal lock scope
            let (rep_id, is_new) = {
                let mut shard = shards[shard_idx].lock()
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

            csv_tx.send((rep_id, record.id().to_string())).await
                .context("CSV send failed")?;

            if is_new {
                uniques_tx.send(ParseOutput::Fastq(record)).await
                    .context("Send unique failed")?;
            }
        }
    }

    if orphan_r1.is_some() {
        return Err(anyhow!("Orphan R1 in worker at end"));
    }

    Ok(())
}

#[allow(dead_code)]
async fn subsample_weighted(
    config: Arc<RunConfig>,
    uniques_stream: mpsc::Receiver<ParseOutput>,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    seed: u64,
    max_subsample: u64,
) -> Result<(
    ReceiverStream<ParseOutput>,
    oneshot::Receiver<u64>,
    JoinHandle<Result<(), anyhow::Error>>,
), PipelineError> {
    if max_subsample == 0 {
        let (tx, rx) = mpsc::channel(1);
        drop(tx);
        let (count_tx, count_rx) = oneshot::channel();
        let _ = count_tx.send(0);
        return Ok((ReceiverStream::new(rx), count_rx, tokio::spawn(async { Ok(()) })));
    }

    let mut rng = config.rng.clone();

    let mut reservoir: Vec<WeightedSampleItem> = Vec::with_capacity(max_subsample as usize);
    let mut stream = ReceiverStream::new(uniques_stream);

    while let Some(item) = stream.next().await {
        let record = match item {
            ParseOutput::Fastq(rec) => rec,
            _ => continue,
        };

        let rep_id = record.id().to_string();
        let weight = duplicate_clusters
            .get(&rep_id)
            .map_or(1u64, |c| c.size as u64);


        let key = rng.random::<f64>().powf(1.0 / weight as f64);

        if reservoir.len() < max_subsample as usize {
            reservoir.push(WeightedSampleItem { key, record, weight });
        } else if let Some(min_item) = reservoir.iter().min_by(|a, b| a.key.partial_cmp(&b.key).unwrap()) {
            if key > min_item.key {
                // Replace the smallest key
                if let Some(idx) = reservoir.iter().position(|r| r.key == min_item.key) {
                    reservoir[idx] = WeightedSampleItem { key, record, weight };
                }
            }
        }
    }

    // Sort by key for consistent output order (same as original)
    reservoir.sort_by(|a, b| a.key.partial_cmp(&b.key).unwrap_or(Ordering::Equal));

    let count = reservoir.len() as u64;
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(count).map_err(|_| PipelineError::Other(anyhow!("Count send failed")))?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let send_task = tokio::spawn(async move {
        for item in reservoir {
            tx.send(ParseOutput::Fastq(item.record)).await?;
        }
        Ok(())
    });

    Ok((ReceiverStream::new(rx), count_rx, send_task))
}




async fn subsample_uniform(
    config: Arc<RunConfig>,
    uniques_rx: mpsc::Receiver<ParseOutput>,
    max_subsample: u64,
    paired: bool,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, JoinHandle<Result<()>>)> {
    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let (count_tx, count_rx) = oneshot::channel();

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(uniques_rx);
        let mut records = Vec::new();
        let mut total = 0u64;

        while let Some(item) = stream.next().await {
            match &item {
                ParseOutput::Fastq(_) => {
                    records.push(item);
                    total += 1;
                }
                other => {
                    return Err(anyhow!(
                        "subsample_uniform: unexpected non-FASTQ item: {:?}",
                        other
                    ));
                }
            }
        }

        // Safety check: paired streams must be even-length
        if paired && total % 2 != 0 {
            return Err(anyhow!("Odd number of FASTQ records in paired mode"));
        }

        let num_units = if paired { total / 2 } else { total };
        let subsample_units = (max_subsample).min(num_units);
        let subsample_size = subsample_units * if paired { 2 } else { 1 };

        if count_tx.send(total).is_err() {
            warn!("Subsample count send failed (receiver might have dropped)");
        }

        let mut rng = config.rng.clone();

        if paired {
            let pairs: Vec<(ParseOutput, ParseOutput)> = records
                .chunks_exact(2)
                .map(|chunk| (chunk[0].clone(), chunk[1].clone()))
                .collect();

            let mut indices: Vec<usize> = (0..pairs.len()).collect();
            indices.shuffle(&mut rng);

            let sampled_indices: HashSet<usize> =
                indices.into_iter().take(subsample_units as usize).collect();

            for (i, (r1, r2)) in pairs.into_iter().enumerate() {
                if sampled_indices.contains(&i) {
                    tx.send(r1).await.map_err(|_| anyhow!("Send R1 failed"))?;
                    tx.send(r2).await.map_err(|_| anyhow!("Send R2 failed"))?;
                }
            }
        } else {
            let mut indices: Vec<usize> = (0..records.len()).collect();
            indices.shuffle(&mut rng);

            let sampled_indices: HashSet<usize> =
                indices.into_iter().take(subsample_units as usize).collect();

            for (i, item) in records.into_iter().enumerate() {
                if sampled_indices.contains(&i) {
                    tx.send(item).await.map_err(|_| anyhow!("Send failed"))?;
                }
            }
        }

        info!(
            "Subsample distributor processed {} uniques, sampled {}",
            total, subsample_size
        );

        Ok(())
    });

    Ok((ReceiverStream::new(rx), count_rx, task))
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
pub async fn minimap2_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
) -> Result<(
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
), PipelineError> {
    let mut cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();

    let base_buffer_size = config.base_buffer_size;
    let channel_buffer = (base_buffer_size * 8).max(1024);

    // 1) Discover NT split .mmi chunks and sort largest-first.
    let nt_split_dir = PathBuf::from(config.args.nt_split_dir.clone());
    let mut chunk_paths: Vec<PathBuf> = Vec::new();
    let mut max_size_bytes: u64 = 0;

    let mut entries = fs::read_dir(&nt_split_dir)
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Cannot read NT split dir: {}", e)))?;

    while let Ok(Some(entry)) = entries.next_entry().await {
        let path = entry.path();
        if path.is_file() && path.extension().map_or(false, |ext| ext == "mmi") {
            if let Ok(meta) = path.metadata() {
                let size: u64 = meta.len();
                if size > max_size_bytes {
                    max_size_bytes = size;
                }
            }
            chunk_paths.push(path);
        }
    }

    if chunk_paths.is_empty() {
        return Err(PipelineError::Other(anyhow!(
            "No .mmi chunk files found in {}",
            nt_split_dir.display()
        )));
    }

    chunk_paths.sort_by_key(|p| {
        std::cmp::Reverse(p.metadata().map(|m| m.len()).unwrap_or(0))
    });

    let num_chunks = chunk_paths.len();
    info!("Found {} minimap2 NT index chunks (sorted largest first)", num_chunks);

    const TARGET_RAM_FRACTION: f64 = 0.55;        // was 0.35
    const HARD_MAX_CONCURRENCY: usize = 14;       // was 8
    const MIN_CONCURRENCY: usize = 4;
    const MIN_THREADS_PER_JOB: usize = 6;         // was 8
    const MAX_THREADS_PER_JOB: usize = 12;        // was 16

    let largest_gb = if max_size_bytes > 0 {
        ((max_size_bytes + (1 << 30) - 1) / (1 << 30)) as f64
    } else {
        70.0
    };

    // Tighter, more realistic per-job RAM estimate
    let est_gb_per_job = (largest_gb * 0.85)
        .max(35.0)
        .min(65.0);

    info!(
        "Largest NT chunk: {:.1} GiB → estimating {:.0} GiB RAM per minimap2 job",
        largest_gb, est_gb_per_job
    );

    let available_ram_gb = config.available_ram as f64 / 1_000_000_000.0;
    let ram_based = (available_ram_gb * TARGET_RAM_FRACTION / est_gb_per_job).floor() as usize;

    let mut concurrency = ram_based.clamp(MIN_CONCURRENCY, HARD_MAX_CONCURRENCY);

    // Extra safety valve
    let estimated_peak = concurrency as f64 * est_gb_per_job;
    if estimated_peak > available_ram_gb * 0.82 {
        let new_conc = concurrency.saturating_sub(2).max(MIN_CONCURRENCY);
        warn!(
            "Estimated peak {:.0} GiB exceeds 82% of available ({:.0} GiB) → reducing concurrency from {} to {}",
            estimated_peak, available_ram_gb, concurrency, new_conc
        );
        concurrency = new_conc;
    }

    let threads_per_job = (config.max_cores / concurrency.max(1))
        .clamp(MIN_THREADS_PER_JOB, MAX_THREADS_PER_JOB);

    info!(
        concat!(
            "NT minimap2 stage : {} concurrent jobs × {} threads/job\n",
            "  est peak RAM: ~{:.0} GiB ({:.0}% of available)\n",
            "  single-job baseline: ~81 s real time (measured)\n",
            "  largest chunk: {:.1} GiB → est/job: {:.0} GiB"
        ),
        concurrency,
        threads_per_job,
        estimated_peak,
        (estimated_peak / available_ram_gb * 100.0).round(),
        largest_gb,
        est_gb_per_job
    );
    // =========================================================

    // Global merged output stream.
    let (merged_tx, merged_rx) = mpsc::channel::<ParseOutput>(channel_buffer * 4);

    // Each worker gets its own job queue and its own output queue.
    let mut worker_job_txs: Vec<mpsc::Sender<(PathBuf, String)>> = Vec::with_capacity(concurrency);
    let mut worker_output_streams: Vec<ReceiverStream<ParseOutput>> = Vec::with_capacity(concurrency);
    let mut worker_handles: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::with_capacity(concurrency);

    for _worker_idx in 0..concurrency {
        let (job_tx, mut job_rx) = mpsc::channel::<(PathBuf, String)>(num_chunks.max(1));
        let (out_tx, out_rx) = mpsc::channel::<ParseOutput>(channel_buffer);

        worker_job_txs.push(job_tx);
        worker_output_streams.push(ReceiverStream::new(out_rx));

        let config = config.clone();
        let r1 = r1_path.clone();
        let r2 = r2_path_opt.clone();

        let handle = tokio::spawn(async move {
            while let Some((chunk_mmi, chunk_name)) = job_rx.recv().await {
                let mut options = HashMap::new();
                options.insert("-c".to_string(), None);
                options.insert("-x".to_string(), Some("sr".to_string()));
                options.insert("--secondary".to_string(), Some("yes".to_string()));

                let mm_config = Minimap2Config {
                    minimap2_index_path: chunk_mmi,
                    r1_path: Some(r1.clone()),
                    r2_path: r2.clone(),
                    option_fields: options,
                    num_threads: Some(threads_per_job),
                };

                let args = generate_cli(MINIMAP2_TAG, &config, Some(&mm_config))
                    .context("Failed to generate minimap2 args")?;

                debug!("Launching minimap2 {} → stdout (streaming PAF)", chunk_name);

                let (mut child, stderr_task) = spawn_cmd(
                    config.clone(),
                    MINIMAP2_TAG,
                    args,
                    config.args.verbose,
                    None
                )
                    .await
                    .context("Failed to spawn minimap2")?;

                if let Some(stdout) = child.stdout.take() {
                    let mut reader = tokio::io::BufReader::new(stdout);
                    let mut line = String::new();

                    while reader.read_line(&mut line).await? > 0 {
                        let trimmed = line.trim_end();
                        if !trimmed.is_empty() && !trimmed.starts_with('#') {
                            let mut bytes = trimmed.as_bytes().to_vec();
                            bytes.push(b'\n');

                            if out_tx.send(ParseOutput::Bytes(Bytes::from(bytes))).await.is_err() {
                                break;
                            }
                        }
                        line.clear();
                    }
                }

                let status = child.wait().await.context("minimap2 wait failed")?;
                if !status.success() {
                    warn!("minimap2 for {} exited with status: {}", chunk_name, status);
                } else {
                    info!("minimap2 for {} completed successfully", chunk_name);
                }

                let _ = stderr_task.await;
            }

            Ok::<(), anyhow::Error>(())
        });

        worker_handles.push(handle);
    }

    // Merge all worker output streams fairly.
    let merger_handle = tokio::spawn({
        let merged_tx = merged_tx;
        let mut merged_stream = futures::stream::select_all(worker_output_streams);
        async move {
            while let Some(item) = merged_stream.next().await {
                let _ = merged_tx.send(item).await;
            }
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(merger_handle);

    // Assign chunks round-robin to workers.
    for (idx, chunk_mmi) in chunk_paths.into_iter().enumerate() {
        let chunk_name = chunk_mmi
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| format!("chunk_{}", idx));

        let worker_idx = idx % worker_job_txs.len();
        let _ = worker_job_txs[worker_idx].send((chunk_mmi, chunk_name)).await;
    }

    drop(worker_job_txs);

    cleanup_tasks.extend(worker_handles);

    info!(
        "[minimap2_non_host_align] Launched {} minimap2 workers — PAF stream now feeding paf_to_m8",
        concurrency
    );

    Ok((
        ReceiverStream::new(merged_rx),
        cleanup_tasks,
        cleanup_receivers,
    ))
}



/// Converts a PAF stream to an m8 stream, splits it for file writing,
/// and returns the main m8 stream along with cleanup handles.
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream (paired or single-end).
///
/// # Returns
/// Tuple:
/// - m8 stream
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
pub async fn paf_to_m8(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    m8_path: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
)> {
    let genome_size = config.args.nt_db_size as f64;

    let paf_record_count = Arc::new(AtomicUsize::new(0));
    let m8_line_count = Arc::new(AtomicUsize::new(0));

    let paf_counter = paf_record_count.clone();
    let m8_counter = m8_line_count.clone();

    let batched_stream = batch_rayon_process(
        config.clone(),
        input_stream,
        move |batch: Vec<u8>| {
            let m8_lines = parse_paf_batch_to_m8(batch.clone(), genome_size);

            // Count input PAF lines in the batch.
            let paf_lines_in_batch = {
                let newline_count = batch.iter().filter(|&&b| b == b'\n').count();
                if batch.is_empty() {
                    0
                } else if *batch.last().unwrap_or(&b'\n') == b'\n' {
                    newline_count
                } else {
                    newline_count + 1
                }
            };

            paf_counter.fetch_add(paf_lines_in_batch, std::sync::atomic::Ordering::Relaxed);
            m8_counter.fetch_add(m8_lines.len(), std::sync::atomic::Ordering::Relaxed);

            m8_lines
        },
        "paf_to_m8",
        16 * 1024 * 1024, // 16 MiB batches
    );

    // Split the converted stream immediately; no extra intermediate mpsc queue.
    use tokio_stream::wrappers::ReceiverStream;

    let (streams, router_handle) = fanout_to_channels(
        batched_stream,
        2,
        "paf_to_m8_stream",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;


    let mut streams = streams;

    let main_stream = ReceiverStream::new(
        streams.remove(0)
    );

    let file_stream = ReceiverStream::new(
        streams.remove(0)
    );

    let write_task = write_byte_stream_to_file(
        &m8_path,
        file_stream,
        config.clone(),
        StreamDataType::JustBytes,
        "paf_to_m8_m8",
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let paf_total = paf_record_count.load(std::sync::atomic::Ordering::Relaxed);
    let m8_total = m8_line_count.load(std::sync::atomic::Ordering::Relaxed);

    info!(
        "[paf_to_m8] wired conversion stage — ~{} PAF lines seen, {} m8 lines emitted so far",
        paf_total, m8_total
    );

    Ok((
        main_stream,
        vec![write_task, router_handle],
        vec![],
    ))
}


/// Loads the taxid (taxid → [species, genus, family])
/// and acc2tax (accession bytes → taxid) DB's
/// async
/// # Arguments
///
/// * `config` - RunConfig struct
///
/// # Returns
/// Result:
///     Arc<AHashMap<Taxid, Lineage>>,
///     Arc<Map<Vec<u8>>>,
pub async fn load_lineage_and_acc2tax_maps(
    config: Arc<RunConfig>,
) -> Result<(
    Arc<AHashMap<Taxid, Lineage>>,
    Arc<Map<Vec<u8>>>,
)> {
    let lineage_path = PathBuf::from(&config.args.taxid_lineages_db);
    let acc2taxid_path = PathBuf::from(&config.args.acc2taxid_db);

    let lineage_future = async move {
        load_taxid_lineages_db(&lineage_path).await
            .map_err(|e| anyhow!("Failed to load lineage DB from {}: {}", lineage_path.display(), e))
    };

    let acc2taxid_future = async move {
        let bytes = tokio::fs::read(&acc2taxid_path).await
            .map_err(|e| anyhow!("Failed to read acc2taxid DB {}: {}", acc2taxid_path.display(), e))?;

        let map = Map::new(bytes)
            .map_err(|e| anyhow!("Failed to parse acc2taxid fst map: {}", e))?;

        Ok::<Arc<Map<Vec<u8>>>, anyhow::Error>(Arc::new(map))
    };

    let (lineage_map, acc2taxid_map) = try_join!(lineage_future, acc2taxid_future)?;

    info!(
        "Loaded taxonomy databases: {} lineages, acc2taxid map with {} entries",
        lineage_map.len(),
        acc2taxid_map.len()
    );

    Ok((lineage_map, acc2taxid_map))
}

/// Sorts an m8 byte stream by read ID (first column) using GNU sort.
/// RAM-first via choose_temp_dir, falls back to NVMe scratch.
/// Guarantees consecutive lines per read ID. Never drops data.
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

    // 2) Only now spawn GNU sort, so it reads a complete file.
    let sort_config = crate::utils::command::sort::SortConfig {
        key: "-k1,1".to_string(),
        parallel: None,
        buffer_size: Some("50%".to_string()),
        temp_dir: Some(temp_dir.path().to_string_lossy().into_owned()),
        output: Some(sorted_path.to_string_lossy().into_owned()),
        extra_fields: HashMap::new(),
        input: Some(unsorted_path.to_string_lossy().into_owned()),
    };

    let sort_args = crate::utils::command::generate_cli(
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

/// Calls hits with the same logic as call_hits_m8 in the CZI pipeline
/// NB:THIS NOW ASSUMES THE M8 IS READ ID SORTED
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `m8_inpout` - line-based m8 format input stream
/// * `taxid_lineages_db_path` - path to sled-backed originally derived from tp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
/// * `duplicate_cluster_sizes_path` - optinal duplicate cluster sizes for weighting
/// * `taxon_blacklist_path` - optional list of excluded taxa
/// * `deuterostome_path` - optional list of excluded taxa
/// * `taxon_whitelist_path` -optional list of included taxa
/// * `min_alignment_length` - will not attempt alignment under this threshold
///
///
/// # Returns
/// Tuple:
/// - m8 stream
///  -hit-summary TSV stream
/// - Vec of cleanup tasks
/// - Vec of cleanup receivers
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
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
), PipelineError> {
    let worker_count = concurrency.max(1);
    let mut cleanup_tasks: Vec<JoinHandle<Result<()>>> = Vec::new();
    let cleanup_receivers: Vec<oneshot::Receiver<Result<()>>> = Vec::new();

    info!(
        "[call_hits_m8:{}] STARTING STREAMING VERSION — workers={}, min_aln_len={}",
        tag, worker_count, min_aln_len
    );

    let (dedup_tx, dedup_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size * 256);
    let (summary_tx, summary_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size * 256);

    #[derive(Debug)]
    enum WorkerMsg {
        ProcessRead {
            read_id: String,
            hits: Vec<M8Record>,
        },
        Flush,
    }

    // One sender per worker; reads are sharded by read_id.
    let mut worker_txs: Vec<mpsc::Sender<WorkerMsg>> = Vec::with_capacity(worker_count);

    for worker_idx in 0..worker_count {
        let (worker_tx, mut worker_rx) = mpsc::channel::<WorkerMsg>(config.base_buffer_size * 16);
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
                    WorkerMsg::ProcessRead { read_id, hits } => {
                        // debug!(
                        //     "[call_hits_m8:{}] worker {} processing read_id={} hits={}",
                        //     worker_tag,
                        //     worker_idx,
                        //     read_id,
                        //     hits.len()
                        // );

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
                                .context("dedup output channel closed")?;

                            summary_tx
                                .send(ParseOutput::Bytes(Bytes::from(reduced.summary)))
                                .await
                                .context("summary output channel closed")?;
                        }
                    }
                    WorkerMsg::Flush => break,
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

            // if total_lines % 1_000_000 == 0 {
            //     info!(
            //         "[call_hits_m8:{}] coordinator saw {} lines so far (current read_id={})",
            //         coordinator_tag, total_lines, read_id
            //     );
            // }

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

                // debug!(
                //     "[call_hits_m8:{}] dispatching read_id={} with {} hits to worker {}",
                //     coordinator_tag,
                //     prev_read_id,
                //     current_hits.len(),
                //     shard
                // );

                worker_txs[shard]
                    .send(WorkerMsg::ProcessRead {
                        read_id: prev_read_id,
                        hits: std::mem::take(&mut current_hits),
                    })
                    .await
                    .context("failed to dispatch completed read to worker")?;

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

                // debug!(
                //     "[call_hits_m8:{}] dispatching final read_id={} with {} hits to worker {}",
                //     coordinator_tag,
                //     read_id,
                //     current_hits.len(),
                //     shard
                // );

                worker_txs[shard]
                    .send(WorkerMsg::ProcessRead {
                        read_id,
                        hits: current_hits,
                    })
                    .await
                    .context("failed to dispatch final read to worker")?;

                total_reads += 1;
            }
        }

        // Flush all workers after all reads have been dispatched.
        for (_, tx) in worker_txs.iter().enumerate() {
            tx.send(WorkerMsg::Flush)
                .await
                .context("failed to flush worker")?;
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
    let (dedup_branches, dedup_router_handle) = fanout_to_channels(
        dedup_stream,
        2,
        "call_hits_m8_dedup",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    // track router task (replaces cleanup_receivers)
    cleanup_tasks.push(dedup_router_handle);

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
    let (summary_branches, summary_router_handle) = fanout_to_channels(
        summary_stream,
        2,
        "call_hits_m8_summary",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    // track router task instead of oneshot receiver
    cleanup_tasks.push(summary_router_handle);

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


// Helper: runs Diamond on one file, waits for both parse and write to finish, returns only the m8 path
async fn run_diamond_single_file(
    config: Arc<RunConfig>,
    query_path: PathBuf,
    db_prefix: PathBuf,
    mut options: HashMap<String, Option<String>>,
    temp_dir: &TempDir,
    label: &str,
) -> Result<PathBuf, PipelineError> {
    let m8_path = temp_dir.path().join(format!("diamond_{}.m8", label));

    options.insert("--query".to_string(), Some(query_path.to_string_lossy().to_string()));

    let diamond_config = DiamondConfig {
        subcommand: DiamondSubcommand::Blastx,
        db: db_prefix,
        r1_path: Some(query_path),
        r2_path: None,
        subcommand_fields: options,
    };

    let diamond_args = generate_cli(DIAMOND_TAG, &config, Some(&diamond_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: DIAMOND_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut diamond_child, _diamond_stderr_task) = spawn_cmd(
        config.clone(),
        DIAMOND_TAG,
        diamond_args,
        config.args.verbose,
        None
    ).await?;

    let diamond_stdout = diamond_child.stdout.take()
        .ok_or_else(|| anyhow!("diamond stdout missing"))?;

    let (m8_tx, m8_rx) = mpsc::channel(config.base_buffer_size * 8);

    let parse_task: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
        let mut reader = TokioBufReader::new(diamond_stdout);
        let mut line = Vec::with_capacity(1024);

        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line).await? == 0 {
                break;
            }
            if !line.is_empty() && line[0] != b'#' {
                let bytes = Bytes::from(std::mem::take(&mut line));
                if m8_tx.send(ParseOutput::Bytes(bytes)).await.is_err() {
                    break;
                }
            }
        }
        Ok(())
    });

    let write_task = write_byte_stream_to_file(
        &m8_path,
        ReceiverStream::new(m8_rx),
        config.clone(),
        StreamDataType::JustBytes,
        &format!("diamond_{}", label),
    ).await?;

    // Wait for diamond process
    let status = diamond_child.wait().await?;
    if !status.success() {
        return Err(PipelineError::ToolExecution {
            tool: DIAMOND_TAG.to_string(),
            error: format!("exited with code {:?}", status.code()),
        });
    }

    // Wait for both tasks to complete
    parse_task.await
        .map_err(|e| PipelineError::Other(anyhow!("parse_task panicked: {}", e)))??;

    write_task.await
        .map_err(|e| PipelineError::Other(anyhow!("write_task panicked: {}", e)))??;

    Ok(m8_path)
}

async fn diamond_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
) -> Result<(
    mpsc::Receiver<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    Vec<TempDir>
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();

    let temp_dir = choose_temp_dir(
        100_000_000_000,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        true
    ).await?;

    let diamond_db = config.args.diamond_db
        .clone()
        .ok_or(PipelineError::MissingArgument("diamond_db required for NR alignment".into()))?;

    let (db_prefix, prep_tasks) = diamond_index_prep(Some(diamond_db), "non_host").await?;
    cleanup_tasks.extend(prep_tasks);

    let optimal_block_size = compute_optimal_block_size(&config)
        .await
        .unwrap_or(12.0);

    let threads = config.thread_allocation(DIAMOND_TAG, None);
    let mut index_chunks = if threads >= 192 { 32 }
    else if threads >= 128 { 24 }
    else if threads >= 96 { 16 }
    else if threads >= 64 { 12 }
    else { 4 };

    index_chunks = index_chunks.min((config.available_ram / 12_000_000_000) as usize);

    info!("Diamond non-host: threads={}, index-chunks={}, block-size={:.1} GB, tmpdir={}",
          threads, index_chunks, optimal_block_size, temp_dir.path().display());

    let diamond_options = HashMap::from([
        ("--mid-sensitive".to_string(), None),
        ("--block-size".to_string(), Some(format!("{:.1}", optimal_block_size))),
        ("-c".to_string(), Some(index_chunks.to_string())),
        ("-f".to_string(), Some("6".to_string())),
        ("--tmpdir".to_string(), Some(temp_dir.path().to_string_lossy().to_string())),
        ("--unal".to_string(), Some("0".to_string())),
    ]);

    // Run R1
    let r1_m8 = run_diamond_single_file(
        config.clone(),
        r1_path,
        db_prefix.clone(),
        diamond_options.clone(),
        &temp_dir,
        "nr_r1",
    ).await?;

    // Run R2 (if present)
    if let Some(r2_path) = r2_path_opt {
        let _r2_m8 = run_diamond_single_file(
            config.clone(),
            r2_path,
            db_prefix,
            diamond_options,
            &temp_dir,
            "nr_r2",
        ).await?;
    }

    // Merge results
    let merged_m8 = temp_dir.path().join("diamond_merged_nr.m8");
    let mut merged = TokioFile::create(&merged_m8)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let mut f1 = TokioFile::open(&r1_m8).await?;
    tokio::io::copy(&mut f1, &mut merged).await?;

    merged.flush().await?;

    let m8_file = TokioFile::open(&merged_m8)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let m8_rx = parse_lines(m8_file, &config, StreamDataType::JustBytes)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok((m8_rx, cleanup_tasks, cleanup_receivers, vec![temp_dir]))
}

/// Generates ytaxon counts from a called m8 stream
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `m8_stream` - line-based m8 format input stream
/// * `summary_stream` - stream of summary information from m8
/// * `duplicate_cluster_sizes_path` - optinal duplicate cluster sizes for weighting
/// * `taxon_blacklist_path` - optional list of excluded taxa
/// * `deuterostome_path` - optional list of excluded taxa
/// * `taxon_whitelist_path` -optional list of included taxa
/// * `count_type` - e.g. NT , type of DB being refernced
/// *  `source_count_type` - optional second source, e.g. NR
///
///
/// # Returns
/// vector of taxon counts
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

/// Concats the nt and nr and writes a json
/// expects the TaxonCount streucture prodiced in generatre_taxon_counts
///
/// # Arguments
///
/// * `nt_counts` - Vec<TaxonCount> from nt
/// * `nr_counts` - Vec<TaxonCount> from nr
/// * `output_path` - JSON pth
///
/// # Returns
///
///
pub async fn combine_taxon_counts(
    nt_counts: &[TaxonCount],
    nr_counts: &[TaxonCount],
    output_path: PathBuf,
) -> Result<(PathBuf, JoinHandle<Result<()>>)> {

    let nt_len = nt_counts.len();
    let nr_len = nr_counts.len();
    let total = nt_len + nr_len;

    let mut all_counts = Vec::with_capacity(total);
    all_counts.extend_from_slice(nt_counts);
    all_counts.extend_from_slice(nr_counts);

    let output_json = json!({
        "pipeline_output": {
            "taxon_counts_attributes": all_counts
        }
    });

    let payload = serde_json::to_vec_pretty(&output_json)
        .map_err(|e| anyhow!("JSON serialization failed: {}", e))?;

    if payload.is_empty() {
        return Err(anyhow!("Empty JSON payload for {}", output_path.display()));
    }

    let data_arc = Arc::new(payload);
    let write_task = write_vecu8_to_file(data_arc, &output_path, 1 << 20)  // 1 Mb
        .await
        .map_err(|e| anyhow!("Failed to spawn JSON write: {}", e))?;

    info!(
        "Combined {} taxon count entries ({} NT + {} NR); spawning write to {}",
        total,
        nt_len,
        nr_len,
        output_path.display()
    );

    Ok((output_path, write_task))
}


/// Conurrently generates small accession map contianing only the assembled contigs
///
/// # Arguments
///
/// * `summary_stream` -stream from hit summary
///
/// # Returns
/// Result of ahashmap of id:accession
///
/// Concurrent collector: hit-summary stream → tiny contig_id → best_accession map
/// Final map merge is offloaded to spawn_blocking (exactly as you requested).
pub async fn collect_hit_summary_to_accession_map_concurrent(
    config: Arc<RunConfig>,
    mut summary_stream: ReceiverStream<ParseOutput>,
) -> Result<AHashMap<String, String>> {
    let concurrency = compute_phase_concurrency(&config, "accession_map_final_merge", 0.2, 4.0, 32, 4);
    let batch_size = compute_batch_size(None, 1024, 512, concurrency);

    let (job_tx, job_rx) = mpsc::channel::<Vec<Vec<u8>>>(concurrency);
    let shared_rx = Arc::new(tokio::sync::Mutex::new(job_rx));

    let mut worker_handles = Vec::with_capacity(concurrency);
    for _ in 0..concurrency {
        let rx = shared_rx.clone();
        let handle = tokio::spawn(async move {
            let mut partial = AHashMap::with_capacity(batch_size);
            while let Some(batch) = {
                let mut g = rx.lock().await;
                g.recv().await
            } {
                let lines: Vec<(String, String)> = batch
                    .par_iter()
                    .filter_map(|bytes| {
                        let line = String::from_utf8_lossy(bytes).trim_end().to_string();
                        if line.is_empty() { return None; }
                        let fields: Vec<&str> = line.split('\t').collect();
                        if fields.len() < 2 { return None; }
                        Some((fields[0].to_string(), fields[1].to_string()))
                    })
                    .collect();
                for (k, v) in lines {
                    partial.insert(k, v);
                }
            }
            Ok::<AHashMap<String, String>, anyhow::Error>(partial)
        });
        worker_handles.push(handle);
    }

    // Producer
    let producer = tokio::spawn(async move {
        let mut batch = Vec::with_capacity(batch_size);
        while let Some(item) = summary_stream.next().await {
            if let ParseOutput::Bytes(b) = item {
                batch.push(b.to_vec());
                if batch.len() >= batch_size {
                    if job_tx.send(std::mem::take(&mut batch)).await.is_err() { break; }
                }
            }
        }
        if !batch.is_empty() {
            let _ = job_tx.send(batch).await;
        }
        drop(job_tx);
        Ok::<(), anyhow::Error>(())
    });

    // Await all workers FIRST (in async context)
    let mut partial_maps = Vec::with_capacity(worker_handles.len());
    for handle in worker_handles {
        let partial = handle.await
            .map_err(|e| anyhow!("Worker task panicked: {}", e))??;
        partial_maps.push(partial);
    }

    // FINAL MERGE OFFLOADED TO spawn_blocking
    let map = tokio::task::spawn_blocking(move || {
        let mut final_map = AHashMap::with_capacity(2_500_000);
        for partial in partial_maps {
            final_map.extend(partial);
        }
        final_map
    }).await
        .map_err(|e| anyhow!("spawn_blocking panicked: {}", e))?;

    producer.await??;

    info!("collect_hit_summary_to_accession_map_concurrent → {} entries (tiny contig map, final merge in spawn_blocking)", map.len());
    Ok(map)
}

/// Generates an annotated FASTA using deduped stream from  dedup_and_subsample
/// and the full duplicate cluster size stream
///
/// # Arguments
/// * `config` - RunConfig struct
/// * `dedup_stream` - deduped FASTQ stream
/// * `cluster_stream - full stream if cluster sizes
/// * `nt_map` -
///
/// # Returns
/// Result of hashmap of id:accession
///
pub fn parse_fasta_batch_to_annotated(
    batch: Vec<u8>,
    duplicate_clusters: &Arc<DashMap<String, ClusterInfo>>,
    nt_map: &HashMap<String, String>,
    nr_map: &HashMap<String, String>,
) -> Vec<Vec<u8>> {
    let mut results = Vec::new();
    let mut reader = needletail::parse_fastx_reader(std::io::Cursor::new(batch))
        .expect("Failed to parse fastx batch");

    while let Some(record) = reader.next() {
        if let Ok(rec) = record {
            let record: SequenceRecord = rec.into();
            let rep_id = record.id().to_string();

            let nr_acc = nr_map.get(&rep_id).cloned().unwrap_or_default();
            let nt_acc = nt_map.get(&rep_id).cloned().unwrap_or_default();
            let is_unidentified = nr_acc.is_empty() && nt_acc.is_empty();

            let dup_count = duplicate_clusters
                .get(&rep_id)
                .map(|r| r.size)
                .unwrap_or(1u64);

            let seq_arc = record.seq_arc();

            if !is_unidentified {
                let id = format!("NR:{}:NT:{}:{}", nr_acc, nt_acc, rep_id);
                let mapped_rec = SequenceRecord::Fasta {
                    id,
                    desc: None,
                    seq: seq_arc,
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
                    seq: seq_arc.clone(),
                };
                if let Ok(bytes) = fasta_rep.to_bytes() {
                    // Unidentified (expanded)
                    let mut tagged_unid = vec![1u8];
                    tagged_unid.extend_from_slice(&bytes);
                    results.push(tagged_unid.clone());

                    if dup_count > 1 {
                        for _ in 1..dup_count {
                            results.push(tagged_unid.clone());
                        }
                    }

                    // Unique unidentified
                    let mut tagged_uniq = vec![2u8];
                    tagged_uniq.extend_from_slice(&bytes);
                    results.push(tagged_uniq);
                }
            }
        }
    }
    results
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
pub fn parse_fasta_batch_to_annotated_streaming(
    batch: Vec<u8>,
    duplicate_clusters: &Arc<DashMap<String, ClusterInfo>>,
    nt_map: &AHashMap<String, String>,   // tiny contig_id → best NT accession
    nr_map: &AHashMap<String, String>,   // tiny contig_id → best NR accession
) -> Vec<Vec<u8>> {
    let mut results = Vec::with_capacity(batch.len() / 512);
    let mut reader = match needletail::parse_fastx_reader(std::io::Cursor::new(batch)) {
        Ok(r) => r,
        Err(e) => {
            error!("parse_fasta_batch_to_annotated_streaming: failed to create reader: {}", e);
            return results;
        }
    };

    while let Some(record) = reader.next() {
        let rec = match record {
            Ok(r) => r,
            Err(e) => {
                warn!("parse_fasta_batch_to_annotated_streaming: bad record: {}", e);
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

        let seq_arc = record.seq_arc();

        if !is_unidentified {
            let id = format!("NR:{}:NT:{}:{}", nr_acc, nt_acc, rep_id);
            let mapped_rec = SequenceRecord::Fasta {
                id,
                desc: None,
                seq: seq_arc,
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
                seq: seq_arc.clone(),
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

pub async fn generate_annotated_fasta_stream(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,           // contig FASTA
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    nt_map: AHashMap<String, String>,                    // tiny contig → NT accession
    nr_map: AHashMap<String, String>,                    // tiny contig → NR accession
    _concurrency: usize,
) -> Result<(
    mpsc::Receiver<ParseOutput>,               // annotated mapped contigs
    mpsc::Receiver<ParseOutput>,               // unidentified (expanded)
    mpsc::Receiver<ParseOutput>,               // unique unidentified (reps only)
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
)> {
    let channel_cap = (config.base_buffer_size / 4).max(10_000);
    let (mapped_tx, mapped_rx) = mpsc::channel(channel_cap);
    let (unid_tx, unid_rx) = mpsc::channel(channel_cap);
    let (uniq_tx, uniq_rx) = mpsc::channel(channel_cap);

    let (done_tx, done_rx) = oneshot::channel();
    let cleanup_receivers = vec![done_rx];

    // Wrap the tiny maps for the rayon closure
    let nt_map_arc = Arc::new(nt_map);
    let nr_map_arc = Arc::new(nr_map);

    let batched_stream = batch_rayon_process(
        config.clone(),
        input_stream,
        move |batch: Vec<u8>| {
            parse_fasta_batch_to_annotated_streaming(   // ← new streaming parser
                                                        batch,
                                                        &duplicate_clusters,
                                                        &nt_map_arc,
                                                        &nr_map_arc,
            )
        },
        "generate_annotated_fasta_streaming",
        16 * 1024 * 1024,
    );

    let demux_handle = tokio::spawn(async move {
        let mut stream = batched_stream;
        while let Some(item) = stream.next().await {
            if let ParseOutput::Bytes(b) = item {
                if b.is_empty() { continue; }
                let tag = b[0];
                let data = &b[1..];

                let mut reader = needletail::parse_fastx_reader(std::io::Cursor::new(data))
                    .map_err(|e| anyhow!("Failed to parse demuxed record: {}", e))?;

                if let Some(record) = reader.next() {
                    let rec = record.map_err(|e| anyhow!("Needletail error: {}", e))?;
                    let seq_rec: SequenceRecord = rec.into();
                    let out_item = ParseOutput::Fasta(seq_rec);

                    match tag {
                        0 => { let _ = mapped_tx.send(out_item).await; },
                        1 => { let _ = unid_tx.send(out_item).await; },
                        2 => { let _ = uniq_tx.send(out_item).await; },
                        _ => warn!("Unknown tag in generate_annotated_fasta_stream: {}", tag),
                    }
                }
            }
        }
        let _ = done_tx.send(Ok(()));
        info!("generate_annotated_fasta_stream completed");
        Ok::<(), anyhow::Error>(())
    });

    Ok((
        mapped_rx,
        unid_rx,
        uniq_rx,
        vec![demux_handle],
        cleanup_receivers,
    ))
}


#[allow(dead_code)]
async fn write_dummy_assembly_files(assembly_dir: &PathBuf) -> Result<(), anyhow::Error> {
    let contigs     = assembly_dir.join("contigs.fasta");
    let contigs_all = assembly_dir.join("contigs_all.fasta");
    let scaffolds   = assembly_dir.join("scaffolds.fasta");
    let bowtie_sam  = assembly_dir.join("read-contig.sam");  // or wherever you place it
    let stats_json  = assembly_dir.join("contig_stats.json");

    let failed_marker = b";ASSEMBLY FAILED\n";
    let no_info       = b"@NO INFO\n";
    let empty_json    = b"{}";

    tokio::fs::write(&contigs, failed_marker).await?;
    tokio::fs::write(&contigs_all, failed_marker).await?;
    tokio::fs::write(&scaffolds, failed_marker).await?;
    tokio::fs::write(&bowtie_sam, no_info).await?;
    tokio::fs::write(&stats_json, empty_json).await?;

    info!("Wrote dummy assembly files after SPAdes failure");
    Ok(())
}

/// Runs spades assembler
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns

async fn spades_assembly(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    assembly_out_dir: PathBuf,   // <out_dir>/assembly
) -> Result<JoinHandle<Result<()>>> {

    let spades_task = tokio::spawn(async move {
        let mut cleanup_tasks = Vec::new();

        // 1. Temp dir
        let est_temp_bytes = config.input_size + MAX_SPADES_WORK_DIR;
        let temp_dir = choose_temp_dir(
            est_temp_bytes,
            &config.ram_temp_dir,
            &config.args.nvme_scratch,
            4,
            true,
        )
            .await?;

        let stdout_log_path = assembly_out_dir.join("spades_stdout.log");

        info!(
            "[spades] starting: r1={}, r2={:?}, assembly_dir={}, temp_dir={},input_size={}, est_temp_bytes={}",
            r1_path.display(),
            r2_path_opt.as_ref().map(|p| p.display().to_string()),
            assembly_out_dir.display(),
            temp_dir.path().display(),
            config.input_size,
            est_temp_bytes
        );

        let mut options = HashMap::new();
        options.insert("--only-assembler".to_string(), None);

        let spades_config = SpadesConfig {
            r1_path: r1_path.clone(),
            r2_path_opt: r2_path_opt.clone(),
            outdir_path: assembly_out_dir.clone(),
            tempdir_path: Some(temp_dir.path().to_path_buf()),
            option_fields: options,
        };

        let spades_args = generate_cli(SPADES_TAG, &config, Some(&spades_config))?;

        info!("[spades] command args: {:?}", spades_args);

        let (mut spades_child, spades_stderr_task) = spawn_cmd(
            config.clone(),
            SPADES_TAG,
            spades_args,
            config.args.verbose,
            None
        )
            .await?;

        info!(
            "[spades] child spawned: stdout_present={}, stderr_present={}",
            spades_child.stdout.is_some(),
            spades_child.stderr.is_some()
        );

        let spades_out_stream = parse_child_output(
            &mut spades_child,
            ChildStream::Stdout,
            ParseMode::Lines,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SPADES_TAG.to_string(),
                error: e.to_string(),
            })?;

        let spades_write_task = tokio::spawn(stream_to_file(
            spades_out_stream,
            stdout_log_path.clone(),
        ));
        cleanup_tasks.push(spades_write_task);

        // stderr is drained and logged by spawn_cmd's internal stderr task
        spades_stderr_task.await??;

        // Check SPAdes exit
        let spades_exit = spades_child.wait().await?;

        info!(
            "[spades] exit status={:?} success={}  stdout_log={}",
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
            error!("SPAdes failed with exit: {:?}", spades_exit);

            // Write dummies so downstream code can detect the failure cleanly.
            let contigs_path = assembly_out_dir.join("contigs.fasta");
            let contigs_all_path = assembly_out_dir.join("contigs_all.fasta");
            let scaffolds_path = assembly_out_dir.join("scaffolds.fasta");
            let sam_path = assembly_out_dir.join("read-contig.sam");
            let failed_marker = b";ASSEMBLY FAILED\n";

            tokio::fs::write(&contigs_path, failed_marker).await?;
            tokio::fs::write(&contigs_all_path, failed_marker).await?;
            tokio::fs::write(&scaffolds_path, failed_marker).await?;
            tokio::fs::write(&sam_path, b"@NO INFO\n").await?;

            info!(
                "[spades] wrote dummy outputs after failure: contigs={}, contigs_all={}, scaffolds={}, sam={}",
                contigs_path.display(),
                contigs_all_path.display(),
                scaffolds_path.display(),
                sam_path.display()
            );

            return Err(anyhow!("SPAdes assembly failed"));
        }

        // Probe the expected SPAdes outputs before leaving this stage.
        let contigs_path = assembly_out_dir.join("contigs.fasta");
        let scaffolds_path = assembly_out_dir.join("scaffolds.fasta");
        let sam_path = assembly_out_dir.join("read-contig.sam");

        for p in [&contigs_path, &scaffolds_path, &sam_path] {
            match fs::metadata(p).await {
                Ok(meta) => info!(
                    "[spades] output exists: {} ({} bytes)",
                    p.display(),
                    meta.len()
                ),
                Err(e) => warn!(
                    "[spades] expected output missing: {} ({})",
                    p.display(),
                    e
                ),
            }
        }

        info!("SPAdes completed successfully");

        Ok(())
    });

    Ok(spades_task)
}

pub async fn process_assembly(
    config: Arc<RunConfig>,
    assembly_out_dir: &PathBuf,   // to make clear this is <out_dir>/assembly
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    paired: bool,
    assembly_headroom: u64,
) -> Result<(
    CoverageOutputs,
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
    Vec<NamedTempFile>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let cleanup_receivers = Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();

    let spades_task = spades_assembly(
        config.clone(),
        r1_path.clone(),
        r2_path_opt.clone(),
        assembly_out_dir.clone(),
    ).await?;

    let spades_result = spades_task.await;

    let spades_success = match spades_result {
        Ok(Ok(())) => {
            info!("SPAdes assembly task completed successfully");
            true
        }
        Ok(Err(e)) => {
            warn!("SPAdes assembly failed: {}", e);
            false
        }
        Err(join_err) => {
            error!("SPAdes task panicked or was cancelled: {}", join_err);
            false
        }
    };

    let contigs_path = assembly_out_dir.join("contigs.fasta");
    let contigs_all_path = assembly_out_dir.join("contigs_all.fasta");
    let scaffolds_path = assembly_out_dir.join("scaffolds.fasta");
    let sam_path = assembly_out_dir.join("read-contig.sam");
    let stats_path = assembly_out_dir.join("contig_stats.json");
    let coverage_path = assembly_out_dir.join("assembly_contig_coverage.json");
    let csv_path = assembly_out_dir.join("assembly_contig_coverage_summary.csv");


    if !spades_success {

        let dummy = b";ASSEMBLY FAILED\n";
        let sam_dummy = b"@NO INFO\n";
        let json_dummy = b"{}";
        let csv_dummy = b"No Contigs\n";

        tokio::fs::write(&contigs_path, dummy).await?;
        tokio::fs::write(&contigs_all_path, dummy).await?;
        tokio::fs::write(&scaffolds_path, dummy).await?;
        tokio::fs::write(&sam_path, sam_dummy).await?;
        tokio::fs::write(&stats_path, json_dummy).await?;
        tokio::fs::write(&coverage_path, json_dummy).await?;
        tokio::fs::write(&csv_path, csv_dummy).await?;

        info!(
            "[assembly] dummy outputs written to {}; no contig-based refinement will run",
            assembly_out_dir.display()
        );

        return Ok((
            CoverageOutputs {
                contigs_fasta: contigs_path,
                contigs_all_fasta: contigs_all_path,
                scaffolds_fasta: scaffolds_path,
                sam_path: sam_path,
                contig_stats_json: stats_path,
                coverage_json: coverage_path,
                coverage_summary_csv: csv_path,
                contig_stats: Arc::new(HashMap::new()),
                read2contig: Arc::new(HashMap::new()),
            },
            ReceiverStream::new({
                let (tx, rx) = mpsc::channel(1);
                drop(tx);
                rx
            }),
            cleanup_tasks,
            cleanup_receivers,
            temp_files,
        ));
    }

    // Success path, proceed with SPAdes refinement
    let contigs_size = file_size(&contigs_path).await?;
    info!(
        "[assembly] raw contigs size: {} bytes ({} GiB)",
        contigs_size,
        contigs_size as f64 / 1_073_741_824.0
    );

    let temp_dir = choose_temp_dir(
        contigs_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        assembly_headroom,
        false,
    )
        .await?;


    let index_dir = assembly_out_dir.join("bowtie_index");
    fs::create_dir_all(&index_dir).await?;
    let index_prefix = index_dir.join("contigs"); // will create contigs.1.bt2, etc.

    let num_index_cores: usize = RunConfig::thread_allocation(&*config, BOWTIE2_TAG, None);

    info!(
        "[assembly] building bowtie2 index: prefix={}, threads={}",
        index_prefix.display(),
        num_index_cores
    );

    let build_status = Command::new("bowtie2-build")
        .args([
            "--threads",
            &num_index_cores.to_string(),
            "--quiet",
            contigs_path.clone().to_str().unwrap(),
            index_prefix.to_str().unwrap(),
        ])
        .status()
        .await
        .map_err(|e| anyhow!("Failed to spawn bowtie2-build: {}", e))?;

    if !build_status.success() {
        return Err(PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: build_status.to_string(),
        });
    }

    info!("[assembly] bowtie2-build completed successfully");

    let bt2_options = HashMap::from([("--very-sensitive".to_string(), None)]);

    let bt2_config_view = Bowtie2Config {
        bt2_index_path: index_prefix,
        r1_path: Some(r1_path.clone()),
        r2_path: r2_path_opt.clone(),
        paired: paired,
        option_fields: bt2_options,
    };

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;

    info!("[assembly] bowtie2 args generated: {:?}", bt2_args);

    let (mut bt2_child, bt2_err_task) = spawn_cmd(
        config.clone(),
        BOWTIE2_TAG,
        bt2_args,
        config.args.verbose,
        None
    )
        .await?;

    cleanup_tasks.push(bt2_err_task);

    let bt2_out_stream = parse_child_output(
        &mut bt2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        &config,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-T".to_string(), Some(temp_dir.path().to_string_lossy().to_string())),
            ("-".to_string(), None),
        ]),
    };

    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    info!("[assembly] samtools sort args generated: {:?}", samtools_sort_args);

    let (samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
        None
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
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            &config,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    use tokio_stream::wrappers::ReceiverStream;

    let (non_host_streams, non_host_router_handle) = fanout_to_channels(
        ReceiverStream::new(samtools_sort_out_stream),
        3,
        "process_assembly_sam",
        &config,                          // ← added
        StreamDataType::JustBytes,        // ← correct (BAM bytes from samtools sort)
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    // track task
    cleanup_tasks.push(non_host_router_handle);

    let mut non_host_streams_iter = non_host_streams.into_iter();

    let bam_for_file = ReceiverStream::new(
        non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let bam_for_stats = ReceiverStream::new(
        non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let bam_for_output = ReceiverStream::new(
        non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let bam_path = assembly_out_dir.join("read-contig.bam");

    let write_sam_task = write_byte_stream_to_file(
        &bam_path,
        bam_for_file,
        config.clone(),
        StreamDataType::JustBytes,
        "process_assembly",
    )
        .await?;
    cleanup_tasks.push(write_sam_task);

    let bam_concurrency = compute_phase_concurrency(
        &config,
        "bam_info",
        0.5, // ~500 MB per thread (local maps + record chunk)
        1.0,
        128, // high cap for large clusters
        4,   // min for small machines like M5 Air
    );
    info!("BAM info concurrency {}", bam_concurrency);

    let (read2contig, contig_stats, _contig_uniques) = generate_info_from_bam_stream(
        bam_for_stats.into_inner(),
        &duplicate_clusters,
        config.args.min_contig_length,
        bam_concurrency,
    )
        .await?;

    info!(
        "[assembly] contig statistics generated: contigs={}, read2contig={}",
        contig_stats.len(),
        read2contig.len()
    );

    Ok((
        CoverageOutputs {
            contigs_fasta: contigs_path,
            contigs_all_fasta: contigs_all_path,
            scaffolds_fasta: scaffolds_path,
            sam_path: sam_path,
            contig_stats_json: stats_path,
            coverage_json: coverage_path,
            coverage_summary_csv: csv_path,
            contig_stats: Arc::new(contig_stats),
            read2contig: Arc::new(read2contig),
        },
        bam_for_output,
        cleanup_tasks,
        cleanup_receivers,
        temp_files,
    ))
}


/// Generates assembly coverage stats from BAM stream (using noodles) and results of process_assembly
/// # Arguments
///
/// * `config` - RunConfig struct
/// * `assembly_bam_rx` - Uncompressed BAM stream
/// * `contigs_fasta_path` - PAth to contigs fasta file
/// * `coverage_json_path` - JSON file containing coverage data
/// * `coverage_csv_path` - CSV of coverage summary
///
///
/// # Returns
/// Result
pub async fn generate_assembly_coverage(
    config: Arc<RunConfig>,
    assembly_bam_rx: ReceiverStream<ParseOutput>,
    contigs_fasta_path: &PathBuf,
    coverage_json_path: &PathBuf,
    coverage_csv_path: &PathBuf,
) -> Result<()> {
    // Early exit
    let meta = fs::metadata(contigs_fasta_path).await
        .map_err(|_| anyhow!("contigs metadata fail"))?;
    if meta.len() < 50 {
        fs::write(coverage_json_path, "{}").await?;
        fs::write(coverage_csv_path, "No Contigs\n").await?;
        return Ok(());
    }

    info!("Streaming coverage: {} → {}", contigs_fasta_path.display(), coverage_json_path.display());

    // Load contigs
    let (contig_order, contig_lengths): (Vec<String>, HashMap<String, usize>) = {
        let mut reader = parse_fastx_file(contigs_fasta_path)
            .map_err(|e| anyhow!("needletail open fail: {e}"))?;
        let mut order = Vec::with_capacity(100_000);
        let mut lengths = HashMap::with_capacity(100_000);

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow!("needletail parse fail: {e}"))?;
            let id = String::from_utf8_lossy(rec.id()).into_owned();
            let len = rec.seq().len();
            if len > 0 {
                order.push(id.clone());
                lengths.insert(id, len);
            }
        }
        if order.is_empty() {
            fs::write(coverage_json_path, "{}").await?;
            fs::write(coverage_csv_path, "No Contigs\n").await?;
            return Ok(());
        }
        (order, lengths)
    };

    // Pre-alloc depths
    let depths: Arc<DashMap<String, Vec<u32>>> = Arc::new(DashMap::with_capacity(contig_lengths.len()));
    for name in &contig_order {
        if let Some(&len) = contig_lengths.get(name) {
            depths.insert(name.clone(), vec![0u32; len + 1]);
        }
    }

    // Stream BAM
    let byte_stream = assembly_bam_rx.map(|item| match item {
        ParseOutput::Bytes(b) => Ok(b),
        _ => Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "non-bytes in BAM stream")),
    });
    let stream_reader = StreamReader::new(byte_stream);
    let mut bam_reader = BamAsyncReader::new(stream_reader);
    let header = bam_reader.read_header().await?;
    let mut record = Record::default();

    while bam_reader.read_record(&mut record).await? > 0 {
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let rid = match record.reference_sequence_id() {
            Some(Ok(rid)) => rid as usize,
            _ => continue,
        };

        let contig_name = match header.reference_sequences().get_index(rid) {
            Some((name, _)) => name.to_string(),
            None => continue,
        };

        let mut entry = depths.get_mut(&contig_name).unwrap_or_else(|| {
            panic!("Contig {} not pre-allocated", contig_name);
        });

        let start_pos = match record.alignment_start() {
            Some(Ok(p)) => p.get(),
            _ => continue,
        };
        let mut pos = start_pos.saturating_sub(1);

        for op in record.cigar().iter().filter_map(Result::ok) {
            let kind = op.kind();
            let len = op.len() as usize;

            match kind {
                OpKind::Match | OpKind::Insertion | OpKind::SequenceMatch | OpKind::SequenceMismatch => {
                    let end = (pos + len).min(entry.len());
                    for i in pos..end {
                        entry[i] = entry[i].saturating_add(1);
                    }
                    pos = end;
                }
                OpKind::Deletion | OpKind::Skip => {
                    pos = (pos + len).min(entry.len());
                }
                OpKind::SoftClip | OpKind::HardClip | OpKind::Pad => {}
            }
        }
    }

    // Parallel stats
    let stats: Vec<(String, Value)> = config.thread_pool.install(|| {
        contig_order.par_iter().filter_map(|name| {
            let entry = depths.get(name)?;
            let len = contig_lengths[name];
            let cov = &entry[..=len];

            let sum: u64 = cov.iter().map(|&d| d as u64).sum();
            let avg = sum as f64 / len as f64;

            let mut sorted = cov.to_vec();
            sorted.sort_unstable();

            let p25 = sorted[((len as f64 * 0.25) as usize).min(len - 1)];
            let p50 = sorted[((len as f64 * 0.50) as usize).min(len - 1)];
            let p75 = sorted[((len as f64 * 0.75) as usize).min(len - 1)];
            let min = *sorted.first()?;
            let max = *sorted.last()?;

            let avg2xcnt = cov.iter().filter(|&&d| d as f64 > 2.0 * avg).count() as f64 / len as f64;
            let cnt0 = cov.iter().filter(|&&d| d == 0).count() as f64 / len as f64;
            let cnt1 = cov.iter().filter(|&&d| d == 1).count() as f64 / len as f64;
            let cnt2 = cov.iter().filter(|&&d| d == 2).count() as f64 / len as f64;

            let coverage_array: Vec<u64> = cov.iter().map(|&d| d as u64).collect();

            Some((
                name.clone(),
                json!({
                    "coverage": coverage_array,
                    "avg": avg,
                    "min": min as u64,
                    "max": max as u64,
                    "p25": p25 as u64,
                    "p50": p50 as u64,
                    "p75": p75 as u64,
                    "avg2xcnt": avg2xcnt,
                    "cnt0": cnt0,
                    "cnt1": cnt1,
                    "cnt2": cnt2,
                    "contig_len": len as u64,
                }),
            ))
        }).collect()
    });

    // Async write — using AsyncWriteExt
    let mut json_file = BufWriter::new(TokioFile::create(coverage_json_path).await?);
    let mut csv_file = BufWriter::new(TokioFile::create(coverage_csv_path).await?);

    csv_file.write_all(b"contig_name,avg,min,max,p25,p50,p75,avg2xcnt,cnt0,cnt1,cnt2\n").await?;
    json_file.write_all(b"{").await?;

    for (i, (name, stat)) in stats.iter().enumerate() {
        if i > 0 {
            json_file.write_all(b",").await?;
        }
        let line = format!(r#""{}":{}"#, name, stat);
        json_file.write_all(line.as_bytes()).await?;
    }
    json_file.write_all(b"}").await?;
    json_file.flush().await?;

    for (name, stat) in &stats {
        let s = stat.as_object().unwrap();
        let line = format!(
            "{},{:.6},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            name,
            s["avg"].as_f64().unwrap_or(0.0),
            s["min"].as_u64().unwrap_or(0),
            s["max"].as_u64().unwrap_or(0),
            s["p25"].as_u64().unwrap_or(0),
            s["p50"].as_u64().unwrap_or(0),
            s["p75"].as_u64().unwrap_or(0),
            s["avg2xcnt"].as_f64().unwrap_or(0.0),
            s["cnt0"].as_f64().unwrap_or(0.0),
            s["cnt1"].as_f64().unwrap_or(0.0),
            s["cnt2"].as_f64().unwrap_or(0.0),
        );
        csv_file.write_all(line.as_bytes()).await?;
    }
    csv_file.flush().await?;

    info!("Coverage complete — {} contigs", stats.len());
    Ok(())
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
    let parts: Vec<&str> = header_line.split(|c: char| c == ',' && header_line.as_bytes().get(header_line.find(',').unwrap_or(0) + 1).map_or(false, |&b| b != b' ')).collect();

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
pub async fn build_reference_fasta_from_selected_genera(
    config: Arc<RunConfig>,
    selected_genera: &HashMap<i32, Vec<String>>,  // genus_taxid → [accessions]
    db_type: &str,
    source_fasta_path: &PathBuf,                  // full NT/NR FASTA (e.g. nt.fa)
    source_index_path: &PathBuf,                  // FST index for fast lookup
) -> Result<(PathBuf, TempDir)> {  // ← Change return to include TempDir
    // Flatten and deduplicate accessions
    let mut accessions: HashSet<String> = HashSet::new();
    for acc_list in selected_genera.values() {
        accessions.extend(acc_list.iter().cloned());
    }

    let temp_dir = choose_temp_dir(1024, &config.ram_temp_dir, &config.args.nvme_scratch, 4, false).await?;  // ← Create early

    if accessions.is_empty() {
        warn!("selected_genera is empty — writing empty reference FASTA");
        let empty_path = temp_dir.path().join(format!("{}_empty_ref.fasta", db_type));

        // Create and write minimal FASTA
        let file = TokioFile::create(&empty_path).await
            .with_context(|| format!("Failed to create empty ref FASTA: {}", empty_path.display()))?;
        let mut writer = BufWriter::new(file);
        writer.write_all(b">empty\nN\n").await
            .with_context(|| format!("Failed to write to empty ref FASTA: {}", empty_path.display()))?;
        writer.flush().await?;
        drop(writer);

        // Verify
        if !empty_path.exists() {
            return Err(anyhow!("Empty ref FASTA {} was not created after write", empty_path.display()));
        }

        info!("Wrote empty ref FASTA to {} (verified exists)", empty_path.display());
        return Ok((empty_path, temp_dir));  // ← Return dir to keep alive
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
    let file = File::open(source_fasta_path)
        .with_context(|| format!("Failed to open source FASTA: {}", source_fasta_path.display()))?;
    let mut reader = BufReader::new(file);

    // Estimate size (use existing temp_dir)

    let output_path = temp_dir.path().join(format!("{}_ref_from_selected_genera.fasta", db_type));
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
        warn!("{} accessions from selected_genera not found in source FASTA", missing);
    }
    if skipped_long > 0 {
        warn!("{} sequences skipped (too long)", skipped_long);
    }

    info!(
        "Reference FASTA built: {} → {} accessions extracted",
        output_path.display(),
        acc_vec.len() - missing
    );

    Ok((output_path, temp_dir))  // ← Keep dir alive
}


/// Summarizes an M8 stream as a
/// # Arguments
/// * `stream` - stream in m8 format
/// * `min_reads_per_genus` - Minimum threshold of reads per genus to be included
/// # Returns
/// read dicttionary of read id's -> Readhit
///  accesions -> Accesiohit
pub async fn summarize_hits(
    config: Arc<RunConfig>,
    mut stream: ReceiverStream<ParseOutput>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
    min_reads_per_genus: usize,
) -> Result<(
    AHashMap<String, Arc<ReadHit>>,
    AHashMap<String, AccessionHit>,
    HashMap<i32, Vec<String>>,
    usize,
)> {
    let concurrency = compute_phase_concurrency(
        &config,
        "summarize_hits",
        0.35,
        4.0,
        96,
        8,
    );

    let (batch_tx, batch_rx) = mpsc::channel::<Vec<ParseOutput>>(concurrency * 2);
    let batch_rx = Arc::new(tokio::sync::Mutex::new(batch_rx));

    let producer = tokio::spawn(async move {
        let mut batch = Vec::with_capacity(8192);
        while let Some(item) = stream.next().await {
            batch.push(item);
            if batch.len() >= 8192 {
                let _ = batch_tx.send(std::mem::take(&mut batch)).await;
            }
        }
        if !batch.is_empty() {
            let _ = batch_tx.send(batch).await;
        }
        Ok::<(), anyhow::Error>(())
    });

    let final_read_dict = Arc::new(DashMap::with_capacity(8_000_000));
    let final_accession_dict = Arc::new(DashMap::with_capacity(2_000_000));
    let final_genus_read = Arc::new(DashMap::<i32, usize>::new());           // usize, +1 per read
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
        let _tr = total_reads.clone();

        worker_handles.push(tokio::spawn(async move {
            let mut local_read_dict: AHashMap<String, Arc<ReadHit>> = AHashMap::with_capacity(8192);
            let mut local_accession_dict: AHashMap<String, AccessionHit> = AHashMap::with_capacity(4096);
            let mut local_genus_read: AHashMap<i32, usize> = AHashMap::with_capacity(2048);
            let mut local_genus_species: AHashMap<i32, HashSet<i32>> = AHashMap::with_capacity(1024);
            let mut local_genus_accessions: AHashMap<i32, HashSet<String>> = AHashMap::with_capacity(1024);

            loop {
                let batch = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };

                let batch = match batch {
                    Some(b) => b,
                    None => break,
                };

                let parsed: Vec<(String, String, i32, i32, i32, u8)> = tokio::task::spawn_blocking({
                    let batch = batch;
                    move || {
                        batch
                            .par_iter()
                            .filter_map(|item| {
                                let bytes = match item.to_bytes() {
                                    Ok(b) => b,
                                    Err(_) => return None,
                                };
                                let line = match std::str::from_utf8(&bytes) {
                                    Ok(s) => s.trim_end(),
                                    Err(_) => return None,
                                };
                                if line.is_empty() { return None; }

                                let parts: Vec<&str> = line.split('\t').collect();
                                if parts.len() != 6 { return None; }

                                let read_id = parts[0].to_string();
                                let accession_id = if parts[1].is_empty() || parts[1] == "-" {
                                    "-".to_string()
                                } else {
                                    parts[1].to_string()
                                };
                                let hit_taxid: i32 = parts[2].parse().unwrap_or(0);
                                let genus_taxid: i32 = parts[3].parse().unwrap_or(0);
                                let family_taxid: i32 = parts[4].parse().unwrap_or(0);
                                let level: u8 = parts[5].parse().unwrap_or(0);

                                Some((read_id, accession_id, hit_taxid, genus_taxid, family_taxid, level))
                            })
                            .collect()
                    }
                }).await.expect("summarize_hits batch parsing panicked");

                for (read_id, accession_id, hit_taxid, genus_taxid, family_taxid, level) in parsed {
                    let species_taxid = if level >= 3 && hit_taxid > 0 { hit_taxid } else { 0 };
                    let assigned_taxid = if level >= 3 {
                        species_taxid
                    } else if level == 2 || level == 1 {
                        genus_taxid
                    } else {
                        hit_taxid
                    };

                    let cluster_size: u64 = dup_clusters
                        .get(&read_id)
                        .map(|e| e.value().size)
                        .unwrap_or(1);

                    local_read_dict.insert(
                        read_id.clone(),
                        Arc::new(ReadHit {
                            level,
                            taxid: assigned_taxid,
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

                        // IMPORTANT: +1 per read (matches Python exactly)
                        *local_genus_read.entry(genus_taxid).or_insert(0) += 1;
                        local_genus_species.entry(genus_taxid).or_default().insert(species_taxid);
                        local_genus_accessions.entry(genus_taxid).or_default().insert(accession_id);
                    }
                }
            }

            // Merge into final structures
            for (k, v) in local_read_dict {
                f_read.insert(k, v);
            }
            for (k, v) in local_accession_dict {
                f_acc.entry(k)
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
        "summarize_hits: processed {} reads, {} unique accessions ({} workers, local maps)",
        total_reads_val, final_accession_dict.len(), concurrency
    );

    // Final selected_genera (exact original Python logic)
    let mut selected_genera: HashMap<i32, Vec<String>> = HashMap::new();
    for entry in final_genus_read.iter() {
        let genus_taxid = *entry.key();
        let read_count = *entry.value();

        if read_count < min_reads_per_genus {
            continue;
        }

        let species_count = final_genus_species.get(&genus_taxid).map(|s| s.len()).unwrap_or(0);
        if species_count <= 1 {
            continue;
        }

        if let Some(acc_set) = final_genus_accessions.get(&genus_taxid) {
            let accessions: Vec<String> = acc_set.iter().cloned().collect();
            selected_genera.insert(genus_taxid, accessions);
        }
    }

    let read_dict_final: AHashMap<String, Arc<ReadHit>> =
        Arc::try_unwrap(final_read_dict)
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

/// Helper for writing empty output of blast_contigs in no assemled results
///
/// # Arguments
/// # Returns
/// Result
async fn write_empty_blast_outputs(
    blast_m8: &PathBuf,
    blast_top_m8: &PathBuf,
    refined_m8: &PathBuf,
    refined_hit_summary: &PathBuf,
    refined_counts_with_dcr: &PathBuf,
    contig_summary_json: &PathBuf,

) -> Result<()> {
    tokio::fs::write(blast_m8, b" ").await?;
    tokio::fs::write(blast_top_m8, b" ").await?;
    tokio::fs::write(refined_m8, b" ").await?;
    tokio::fs::write(refined_hit_summary, b" ").await?;
    tokio::fs::write(refined_counts_with_dcr, b" ").await?;
    tokio::fs::write(contig_summary_json, b"[]").await?;
    Ok(())
}

/// Updates the read_dict in blast_contigs from streamed m8 data.
/// Updates it in place using a lock, does not return it
/// # Arguments
/// * `read2contig` - Assignment of reads to contigs
/// * `top_m8_stream` - stream in m8 format, top as found in get_top_m8_nt/nr
/// * `read_dict` - HashMap of ReadID -> ReadHits
/// * `lneage_map` - Hashmap os taxon id -> full lineage (family, genus, species)
/// * `accession_map` - accession → (species, genus, family)
/// * `should_keep` - combination of taxon lists and deuterostome filters
/// * `db_type` - denotes NT or NR DB's to keep intermdiate and output files distinguished form each other
/// * `contig2lineage_tx` - json sender
/// * `read2blastm8_tx` - refined m8 sender
/// * 'updated_tx` - updated dict sender
/// * `added_tx` - added to dict snde
/// # Returns
/// Result
pub async fn update_read_dict(
    config: Arc<RunConfig>,
    read2contig: Arc<HashMap<String, String>>,
    mut top_m8_stream: ReceiverStream<ParseOutput>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    accession_map: Arc<AHashMap<String, AccessionHit>>,
    should_keep: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    db_type: &str,
    contig2lineage_tx: Sender<ParseOutput>,
    read2blastm8_tx: Sender<ParseOutput>,
    updated_tx: Sender<ParseOutput>,
    added_tx: Sender<ParseOutput>,
) -> Result<()> {

    let update_concurrency = compute_phase_concurrency(
        &config,
        "update_read_dict",
        0.05,
        8.0,
        128,
        4,
    );
    info!("update_read_dict concurrency: {} jobs", update_concurrency);

    #[derive(Default)]
    struct FanoutBatch {
        updated: Vec<(String, Arc<ReadHit>)>,
        added: Vec<(String, Arc<ReadHit>)>,
        contig2lineage_lines: Vec<Vec<u8>>,
        read2blastm8_lines: Vec<Vec<u8>>,
    }

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

    // Snapshot the current read dict once so the fan-out stage does not keep
    // taking the mutex for every contig/read pair.
    let existing_reads: AHashMap<String, Arc<ReadHit>> = {
        let guard = read_dict
            .lock()
            .map_err(|e| anyhow!("read_dict lock poisoned: {}", e))?;
        guard.clone()
    };
    let existing_reads = Arc::new(existing_reads);

    // Invert read2contig once: contig -> reads.
    let mut contig2reads: AHashMap<String, Vec<String>> = AHashMap::with_capacity(read2contig.len());
    for (read_id, contig_id) in read2contig.iter() {
        contig2reads
            .entry(contig_id.clone())
            .or_default()
            .push(read_id.clone());
    }
    let contig2reads = Arc::new(contig2reads);

    // Pass 1: read top m8 stream and build one merged row per contig.
    // This mirrors the Python logic: contig2accession stores the merged row,
    // and duplicate contig rows are merged into that row.
    let mut contig2accession: AHashMap<String, (String, M8Record, Lineage)> = AHashMap::new();
    let db_type_upper = db_type.to_uppercase();

    while let Some(item) = top_m8_stream.next().await {
        let bytes = match item {
            ParseOutput::Bytes(b) => b,
            _ => continue,
        };

        let line = String::from_utf8_lossy(&bytes);
        let line_trimmed = line.trim_end();
        if line_trimmed.is_empty() {
            continue;
        }

        let m8 = if db_type == "nt" {
            M8Record::parse_line_nt(line_trimmed)?
        } else {
            M8Record::parse_line_nr(line_trimmed)?
        };

        let contig_id = m8.qname.clone();
        let accession_id = m8.tname.clone();

        let taxid = accession_map
            .get(&accession_id)
            .map(|hit| hit.taxid)
            .unwrap_or(0);

        let lineage = if taxid > 0 {
            lineage_map.get(&taxid).cloned().unwrap_or([0; 3])
        } else {
            [0; 3]
        };

        if !(should_keep)(&lineage) {
            continue;
        }

        if let Some((_, curr_row, _)) = contig2accession.get_mut(&contig_id) {
            merge_m8_rows(curr_row, &m8);
        } else {
            contig2accession.insert(
                contig_id.clone(),
                (accession_id.clone(), m8.clone(), lineage),
            );
        }
    }

    // Pass 2: parallel fan-out by contig.
    // Each worker only builds local vectors, then we merge once at the end.
    let fanout_batches: Vec<FanoutBatch> = {
        let contig2accession = Arc::new(contig2accession);
        let accession_map = accession_map.clone();
        let existing_reads = existing_reads.clone();
        let contig2reads = contig2reads.clone();

        tokio::task::spawn_blocking(move || -> Result<Vec<FanoutBatch>> {
            let batches = contig2reads
                .par_iter()
                .map(|(contig_id, read_ids)| {
                    let mut batch = FanoutBatch::default();

                    let Some((accession, m8_row, lineage)) = contig2accession.get(contig_id) else {
                        return batch;
                    };

                    let Some(acc_hit) = accession_map.get(accession) else {
                        return batch;
                    };

                    let species_taxid = acc_hit.taxid;
                    let genus_taxid = acc_hit.genus_taxid;
                    let family_taxid = acc_hit.family_taxid;

                    let mut lineage_line = json!({
                        "contig_name": contig_id,
                        "species_taxid": lineage[0],
                        "genus_taxid": lineage[1],
                        "family_taxid": lineage[2],
                    })
                        .to_string();
                    lineage_line.push('\n');
                    batch.contig2lineage_lines.push(lineage_line.into_bytes());

                    let mut m8_line = m8_row.to_tab_string();
                    m8_line.push('\n');
                    let m8_bytes = m8_line.into_bytes();

                    for read_id in read_ids {
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

                            batch.updated.push((read_id.clone(), Arc::new(hit)));
                        } else {
                            let hit = ReadHit {
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
                            };

                            batch.added.push((read_id.clone(), Arc::new(hit)));
                        }

                        batch.read2blastm8_lines.push(m8_bytes.clone());
                    }

                    batch
                })
                .collect();

            Ok(batches)
        })
            .await
            .map_err(|e| anyhow!("update_read_dict fan-out task panicked: {}", e))??
    };

    // Merge local results once.
    let mut updated_entries: Vec<(String, Arc<ReadHit>)> = Vec::new();
    let mut added_entries: Vec<(String, Arc<ReadHit>)> = Vec::new();
    let mut contig2lineage_lines: Vec<Vec<u8>> = Vec::new();
    let mut read2blastm8_lines: Vec<Vec<u8>> = Vec::new();

    for mut batch in fanout_batches {
        updated_entries.append(&mut batch.updated);
        added_entries.append(&mut batch.added);
        contig2lineage_lines.append(&mut batch.contig2lineage_lines);
        read2blastm8_lines.append(&mut batch.read2blastm8_lines);
    }

    {
        let mut dict = read_dict
            .lock()
            .map_err(|e| anyhow!("read_dict lock poisoned: {}", e))?;

        for (read_id, hit) in &updated_entries {
            dict.insert(read_id.clone(), hit.clone());
        }
        for (read_id, hit) in &added_entries {
            dict.insert(read_id.clone(), hit.clone());
        }
    }

    for line in contig2lineage_lines {
        contig2lineage_tx
            .send(ParseOutput::Bytes(Bytes::from(line)))
            .await
            .map_err(|_| anyhow!("contig2lineage_tx closed"))?;
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

    Ok(())
}

pub async fn generate_contig_summary_json(
    read2contig: Arc<HashMap<String, String>>,
    contig2lineage: AHashMap<String, [i32; 3]>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    db_type: &str,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,   // ← changed to DashMap
    min_contig_size: u64,
    output_tx: Sender<ParseOutput>,
) -> Result<()> {
    let mut genus_summary: AHashMap<i32, AHashMap<String, [u64; 2]>> = AHashMap::new();
    let mut species_summary: AHashMap<i32, AHashMap<String, [u64; 2]>> = AHashMap::new();

    // lock once
    {
        let dict = read_dict.lock().unwrap();

        for (read_id, hit_arc) in dict.iter() {
            let hit: &ReadHit = hit_arc.as_ref();

            let contig = read2contig
                .get(read_id)
                .cloned()
                .unwrap_or_else(|| "*".to_string());

            let (species_taxid, genus_taxid) = if contig != "*" {
                if let Some(lineage) = contig2lineage.get(&contig) {
                    (lineage[0], lineage[1])
                } else {
                    (hit.species_taxid, hit.genus_taxid)
                }
            } else {
                (hit.species_taxid, hit.genus_taxid)
            };


            let cluster_size = duplicate_clusters
                .get(read_id.as_str())
                .map(|entry| entry.value().size)
                .unwrap_or(1u64);

            // Species-level aggregation
            let sp_entry = species_summary
                .entry(species_taxid)
                .or_default()
                .entry(contig.clone())
                .or_insert([0, 0]);
            sp_entry[0] += 1;           // unique reads
            sp_entry[1] += cluster_size; // weighted by duplicates

            // Genus-level aggregation
            let gn_entry = genus_summary
                .entry(genus_taxid)
                .or_default()
                .entry(contig)
                .or_insert([0, 0]);
            gn_entry[0] += 1;
            gn_entry[1] += cluster_size;
        }
    } // ← MutexGuard dropped here — no await while holding lock

    // Now emit JSON — safe to await
    for (tax_level, summary_map) in [(1, species_summary), (2, genus_summary)] {
        for (taxid, contig_counts) in summary_map {
            let filtered: HashMap<String, u64> = contig_counts
                .into_iter()
                .filter(|(_, [unique, _])| *unique >= min_contig_size)
                .map(|(contig, [_, weighted])| (contig, weighted))
                .collect();

            if !filtered.is_empty() {
                let entry = json!({
                    "taxid": taxid,
                    "tax_level": tax_level,
                    "count_type": db_type.to_uppercase(),
                    "contig_counts": filtered
                });
                let line = format!("{entry}\n");
                output_tx
                    .send(ParseOutput::Bytes(Bytes::from(line)))
                    .await
                    .map_err(|_| anyhow!("contig_summary_tx dropped"))?;
            }
        }
    }

    Ok(())
}

// Helper: pull the first tab-delimited field from a complete line.
fn first_field(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes)
        .split('\t')
        .next()
        .unwrap_or("")
        .trim_end()
        .to_string()
}

// Helper: collect a stream into a map keyed by the first field,
// while preserving first-seen key order for deterministic appends.
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

// Helper: collect a stream into ordered rows only.
#[allow(dead_code)]
async fn collect_rows(
    mut stream: ReceiverStream<ParseOutput>,
) -> Result<Vec<Bytes>> {
    let mut rows = Vec::new();

    while let Some(item) = stream.next().await {
        let ParseOutput::Bytes(bytes) = item else {
            continue;
        };
        rows.push(bytes);
    }

    Ok(rows)
}

// collect a stream into keyed rows and preserve first-seen key order.
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

pub async fn generate_m8_and_hit_summary(
    config: Arc<RunConfig>,
    updated_reads_stream: ReceiverStream<ParseOutput>,
    added_reads_stream: ReceiverStream<ParseOutput>,
    blast_hits_stream: ReceiverStream<ParseOutput>,
    original_hit_summary_stream: ReceiverStream<ParseOutput>,
    original_deduped_m8_stream: ReceiverStream<ParseOutput>,
    refined_m8_tx: Sender<ParseOutput>,
    refined_hit_summary_tx: Sender<ParseOutput>,
) -> Result<()> {
    // Drain all inputs concurrently.
    let updated_fut = collect_line_map(updated_reads_stream);
    let added_fut = collect_line_map(added_reads_stream);
    let blast_hits_fut = collect_line_map(blast_hits_stream);
    let hit_rows_fut = collect_keyed_rows(original_hit_summary_stream);
    let m8_rows_fut = collect_keyed_rows(original_deduped_m8_stream);

    let (
        (updated_reads, _updated_order),
        (added_reads, added_order),
        (blast_hits, blast_order),
        (original_hit_rows, original_hit_order),
        (original_m8_rows, original_m8_order),
    ) = tokio::try_join!(
        updated_fut,
        added_fut,
        blast_hits_fut,
        hit_rows_fut,
        m8_rows_fut
    )?;

    // Original IDs for append checks.
    let original_hit_ids: HashSet<String> = original_hit_order.iter().cloned().collect();
    let original_m8_ids: HashSet<String> = original_m8_order.iter().cloned().collect();

    let thread_pool = config.thread_pool.clone();

    let (rewritten_m8, rewritten_hit_summary) = thread_pool.install(|| {
        rayon::join(
            || {
                original_m8_rows
                    .par_iter()
                    .map(|(qseqid, row)| {
                        blast_hits
                            .get(qseqid)
                            .cloned()
                            .unwrap_or_else(|| row.clone())
                    })
                    .collect::<Vec<Bytes>>()
            },
            || {
                original_hit_rows
                    .par_iter()
                    .map(|(read_id, row)| {
                        updated_reads
                            .get(read_id)
                            .or_else(|| added_reads.get(read_id))
                            .cloned()
                            .unwrap_or_else(|| row.clone())
                    })
                    .collect::<Vec<Bytes>>()
            },
        )
    });

    // Emit M8.
    let m8_tx = refined_m8_tx;
    let m8_send = async move {
        for row in rewritten_m8 {
            m8_tx
                .send(ParseOutput::Bytes(row))
                .await
                .map_err(|_| anyhow!("refined_m8_tx dropped"))?;
        }

        // Append blast rows that were not present in the original deduped M8.
        for read_id in blast_order {
            if !original_m8_ids.contains(&read_id) {
                if let Some(row) = blast_hits.get(&read_id) {
                    m8_tx
                        .send(ParseOutput::Bytes(row.clone()))
                        .await
                        .map_err(|_| anyhow!("refined_m8_tx dropped"))?;
                }
            }
        }

        drop(m8_tx);
        Ok::<(), anyhow::Error>(())
    };

    // Emit hit summary.
    let hit_tx = refined_hit_summary_tx;
    let hit_send = async move {
        for row in rewritten_hit_summary {
            hit_tx
                .send(ParseOutput::Bytes(row))
                .await
                .map_err(|_| anyhow!("refined_hit_summary_tx dropped"))?;
        }

        // Append added reads that were not already in the original hit summary.
        for read_id in added_order {
            if !original_hit_ids.contains(&read_id) {
                if let Some(row) = added_reads.get(&read_id) {
                    hit_tx
                        .send(ParseOutput::Bytes(row.clone()))
                        .await
                        .map_err(|_| anyhow!("refined_hit_summary_tx dropped"))?;
                }
            }
        }

        drop(hit_tx);
        Ok::<(), anyhow::Error>(())
    };

    tokio::try_join!(m8_send, hit_send)?;
    Ok(())
}


pub async fn blast_contigs(
    config: Arc<RunConfig>,
    db_type: &'static str,
    deduped_m8_stream: ReceiverStream<ParseOutput>,
    hit_summary_stream: ReceiverStream<ParseOutput>,
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
    mpsc::Receiver<ParseOutput>,
    mpsc::Receiver<ParseOutput>,
    mpsc::Receiver<ParseOutput>,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
    Vec<NamedTempFile>,
)> {
    info!("blast_contigs start!");
    use std::time::{Duration, Instant};

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
    let mut temp_files: Vec<NamedTempFile> = Vec::new();

    let out_dir = config.out_dir.join(format!("blast_{}", db_type));
    tokio::fs::create_dir_all(&out_dir).await?;
    info!(
        "[blast_contigs:{}] output directory ready: {}",
        db_type,
        out_dir.display()
    );

    let blast_m8_path = out_dir.join("blast.m8");
    let blast_top_m8_path = out_dir.join("blast_top.m8");
    let refined_m8_path = out_dir.join("refined.m8");
    let refined_hit_summary_path = out_dir.join("refined_hit_summary.tab");
    let refined_counts_path = out_dir.join("refined_counts_with_dcr.json");
    let contig_summary_path = out_dir.join("contig_summary.json");

    info!(
        "[blast_contigs:{}] output files prepared: blast_m8={}, blast_top_m8={}, refined_m8={}, refined_hit_summary={}, refined_counts={}, contig_summary={}",
        db_type,
        blast_m8_path.display(),
        blast_top_m8_path.display(),
        refined_m8_path.display(),
        refined_hit_summary_path.display(),
        refined_counts_path.display(),
        contig_summary_path.display()
    );

    let contig_size = file_size(assembled_contig_fasta).await?;
    let ref_size = file_size(reference_fasta).await?;

    info!(
        "[blast_contigs:{}] input sizes — assembled_contig={} bytes, reference={} bytes",
        db_type,
        contig_size,
        ref_size
    );

    if contig_size < MIN_REF_FASTA_SIZE || ref_size < MIN_ASSEMBLED_CONTIG_SIZE {
        warn!(
            "[blast_contigs:{}] skipping BLAST: contig_size={} ref_size={} (thresholds contig>= {}, ref>= {})",
            db_type,
            contig_size,
            ref_size,
            MIN_REF_FASTA_SIZE,
            MIN_ASSEMBLED_CONTIG_SIZE
        );

        write_empty_blast_outputs(
            &blast_m8_path,
            &blast_top_m8_path,
            &refined_m8_path,
            &refined_hit_summary_path,
            &refined_counts_path,
            &contig_summary_path,
        )
            .await?;

        let final_read_dict = read_dict.lock().unwrap().clone();

        let (_m8_tx, m8_rx) = mpsc::channel(1);
        let (_hit_tx, hit_rx) = mpsc::channel(1);
        let (_top_tx, top_rx) = mpsc::channel(1);

        info!(
            "[blast_contigs:{}] empty-output fast path complete after {:?}",
            db_type,
            fn_start.elapsed()
        );

        return Ok((
            final_read_dict,
            taxon_counts,
            vec![],
            m8_rx,
            hit_rx,
            top_rx,
            cleanup_tasks,
            cleanup_receivers,
            temp_files,
        ));
    }

    let temp_dir = choose_temp_dir(
        ref_size + contig_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        blast_headroom,
        true,
    )
        .await?;
    info!(
        "[blast_contigs:{}] temp dir chosen: {}",
        db_type,
        temp_dir.path().display()
    );

    let blastdb_suffix = format!("{}_blastindex", db_type);
    let blastdb_ram_path = NamedTempFile::with_suffix_in(blastdb_suffix, &temp_dir)
        .map_err(|e| PipelineError::Other(e.into()))?;
    let blastdb_path = blastdb_ram_path.path().to_owned();
    temp_files.push(blastdb_ram_path);

    info!(
        "[blast_contigs:{}] blastdb path: {}",
        db_type,
        blastdb_path.display()
    );


    info!("blast_contigs right before making its db");

    // === 1. Build config ===
    let makeblastdb_config = MakeblastdbConfig {
        input: reference_fasta.clone(),
        dbtype: if db_type == NT_TAG {
            "nucl".to_string()
        } else {
            "prot".to_string()
        },
        output: blastdb_path.clone(),
        option_fields: HashMap::new(),
    };


    let makeblastdb_args = generate_cli(MAKEBLASTDB_TAG, &config, Some(&makeblastdb_config))
        .map_err(|e| anyhow::anyhow!("Failed to generate makeblastdb args: {}", e))?;

    let log_path = out_dir.join(format!("{}_makeblastdb.log", db_type));

    let (mut make_child, make_err_task) = spawn_cmd(
        config.clone(),
        MAKEBLASTDB_TAG,
        makeblastdb_args,
        config.args.verbose,
        Some(log_path.clone()),
    )
        .await
        .map_err(|e| anyhow::anyhow!("Failed to spawn makeblastdb: {}", e))?;

    cleanup_tasks.push(make_err_task);

    info!("[makeblastdb:{}] waiting for database creation...", db_type);

    let status = make_child.wait().await
        .map_err(|e| anyhow::anyhow!("makeblastdb wait failed: {}", e))?;

    if !status.success() {
        return Err(anyhow::anyhow!(
        "makeblastdb failed with status: {:?}. Check log: {}",
        status.code(),
        log_path.display()
    ));
    }

    info!("[makeblastdb:{}] completed successfully", db_type);

    // verify db
    let db_prefix = blastdb_path.clone();
    let expected_extensions = if db_type == NT_TAG {
        vec!["nin", "nhr", "nog", "nsq"]
    } else {
        vec!["pin", "phr", "pog", "psq"]
    };

    let mut db_ok = false;
    for ext in expected_extensions {
        let path = db_prefix.with_extension(ext);
        if path.exists() && path.metadata().map(|m| m.len() > 0).unwrap_or(false) {
            db_ok = true;
            debug!("Found valid db file: {}", path.display());
            break;
        }
    }

    if !db_ok {
        return Err(anyhow::anyhow!(
        "makeblastdb succeeded but produced no valid index files at prefix {} (check log: {})",
        db_prefix.display(),
        log_path.display()
    ));
    }

    let blast_command = if db_type == NT_TAG {
        BLASTN_TAG
    } else {
        BLASTX_TAG
    };


    let blast_args = if db_type == NT_TAG {
        let blastn_config = BlastnConfig {
            query: assembled_contig_fasta.clone(),
            db: blastdb_path.clone(),
            outfmt: "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".to_string(),
            evalue: 1e-10,
            max_target_seqs: 5000,
            option_fields: HashMap::new(),
        };
        generate_cli(BLASTN_TAG, &config, Some(&blastn_config))?
    } else {
        let blastx_config = BlastxConfig {
            query: assembled_contig_fasta.clone(),
            db: blastdb_path.clone(),
            outfmt: 6,
            evalue: 1e-10,
            num_alignments: 5,
            option_fields: HashMap::new(),
        };
        generate_cli(BLASTX_TAG, &config, Some(&blastx_config))?
    };

    info!(
        "[blast_contigs:{}] launching {} with query={} db={}",
        db_type,
        blast_command,
        assembled_contig_fasta.display(),
        blastdb_path.display()
    );
    info!(
        "[blast_contigs:{}] blast args prepared for {}",
        db_type,
        blast_command
    );
    info!("blast_contigs right before spawning");
    let blast_spawn_start = Instant::now();
    let (mut blast_child, err_task) = spawn_cmd(
        config.clone(),
        blast_command,
        blast_args,
        config.args.verbose,
        None
    )
        .await?;
    cleanup_tasks.push(err_task);

    info!(
        "[blast_contigs:{}] {} spawned after {:?}",
        db_type,
        blast_command,
        blast_spawn_start.elapsed()
    );
    info!("blast_contigs right before output aAAAAAAAAAA");
    let blast_out_stream = parse_child_output(
        &mut blast_child,
        ChildStream::Stdout,
        ParseMode::Lines,
        &config,
    )
        .await?;
    info!(
        "[blast_contigs:{}] child stdout parser ready (buffer_size={})",
        db_type,
        config.base_buffer_size
    );

    // materialize raw blast output BEFORE heavy fanout
    // This breaks backpressure and guarantees clean EOF for the complex downstream graph
    let blast_m8_path = materialize_to_temp(
        config.clone(),
        ReceiverStream::new(blast_out_stream),
        &format!("blast_{}_raw", db_type),
        StreamDataType::JustBytes,
        4,                    // headroom
        false,                // prefer RAM for BLAST output (usually small)
        None,                 // let heuristic decide size
    ).await?;

    // Re-open as fresh stream (clean state, no backpressure from blastn)
    let file = tokio::fs::File::open(&blast_m8_path).await
        .context("Failed to reopen materialized blast m8")?;

    let blast_m8_stream = parse_bytes(file, &config, StreamDataType::JustBytes).await?;

    info!("[blast_contigs:{}] materialized blast output to {} (ready for fanout)",
          db_type, blast_m8_path.display());


    // ────────────────────────────────────────────────
    // Top hits and results processing
    // ────────────────────────────────────────────────
    let (top_internal_tx, top_internal_rx) = mpsc::channel(1024);
    let (top_for_update_tx, top_for_update_rx) = mpsc::channel(1024);
    let (top_for_caller_tx, top_for_caller_rx) = mpsc::channel(1024);

    info!(
        "[blast_contigs:{}] top-hits fanout channels created (1024 each)",
        db_type
    );

    let top_handle = if db_type == NT_TAG {
        let top_nt_concurrency = compute_phase_concurrency(
            &config,
            "top_m8_nt",
            0.25,
            1.0,
            128,
            2,
        );
        info!(
            "[blast_contigs:{}] NT top_m8 concurrency={}",
            db_type,
            top_nt_concurrency
        );

        let avg_line_bytes = 200;
        let target_batch_mb = if config.available_ram < 64 { 100 } else { 250 };
        let est_total_lines = None;
        let top_nt_batch_size = compute_batch_size(
            est_total_lines,
            avg_line_bytes,
            target_batch_mb,
            top_nt_concurrency,
        );
        info!(
            "[blast_contigs:{}] NT top_m8 batch_size={}",
            db_type,
            top_nt_batch_size
        );

        let top_start = Instant::now();
        tokio::spawn(async move {
            info!("[blast_contigs:{}] get_top_m8_nt started", db_type);
            let res = get_top_m8_nt(
                ReceiverStream::new(blast_m8_stream),
                top_internal_tx,
                top_nt_concurrency,
                top_nt_batch_size,
            )
                .await;
            info!(
                "[blast_contigs:{}] get_top_m8_nt finished after {:?}",
                db_type,
                top_start.elapsed()
            );
            res
        })
    } else {
        let top_nr_concurrency = compute_phase_concurrency(
            &config,
            "top_m8_nr",
            0.5,
            1.5,
            128,
            2,
        );
        info!(
            "[blast_contigs:{}] NR top_m8 concurrency={}",
            db_type,
            top_nr_concurrency
        );

        let avg_line_bytes = 200;
        let target_batch_mb = if config.available_ram < 64 { 100 } else { 250 };
        let est_total_lines = None;
        let top_nr_batch_size = compute_batch_size(
            est_total_lines,
            avg_line_bytes,
            target_batch_mb,
            top_nr_concurrency,
        );
        info!(
            "[blast_contigs:{}] NR top_m8 batch_size={}",
            db_type,
            top_nr_batch_size
        );

        let top_start = Instant::now();
        tokio::spawn(async move {
            info!("[blast_contigs:{}] get_top_m8_nr started", db_type);
            let res = get_top_m8_nr(
                ReceiverStream::new(blast_m8_stream),
                top_internal_tx,
                top_nr_concurrency,
                top_nr_batch_size,
            )
                .await;
            info!(
                "[blast_contigs:{}] get_top_m8_nr finished after {:?}",
                db_type,
                top_start.elapsed()
            );
            res
        })
    };

    let forward_handle = tokio::spawn({
        let mut rx = top_internal_rx;
        let tx_update = top_for_update_tx;
        let tx_caller = top_for_caller_tx;
        async move {
            let start = Instant::now();
            let mut total = 0usize;
            let mut last_log = Instant::now();

            info!("[blast_contigs:{}] top forward task started", db_type);

            while let Some(item) = rx.recv().await {
                total += 1;

                let update_send_start = Instant::now();
                if tx_update.send(item.clone()).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] top forward: update branch dropped at item {} after {:?}",
                        db_type,
                        total,
                        start.elapsed()
                    );
                    break;
                }
                let update_wait = update_send_start.elapsed();

                let caller_send_start = Instant::now();
                if tx_caller.send(item).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] top forward: caller branch dropped at item {} after {:?}",
                        db_type,
                        total,
                        start.elapsed()
                    );
                    break;
                }
                let caller_wait = caller_send_start.elapsed();

                if update_wait >= Duration::from_millis(50)
                    || caller_wait >= Duration::from_millis(50)
                    || last_log.elapsed() >= Duration::from_secs(10)
                {
                    info!(
                        "[blast_contigs:{}] top forward progress: total={} update_wait={:?} caller_wait={:?} elapsed={:?}",
                        db_type,
                        total,
                        update_wait,
                        caller_wait,
                        start.elapsed()
                    );
                    last_log = Instant::now();
                }
            }

            info!(
                "[blast_contigs:{}] top forward complete: total={} elapsed={:?}",
                db_type,
                total,
                start.elapsed()
            );

            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(forward_handle);

    let (contig2lineage_tx, contig2lineage_rx) = mpsc::channel(1024);
    let (read2blastm8_tx, read2blastm8_rx) = mpsc::channel(1024);
    let (updated_tx, updated_rx) = mpsc::channel(1024);
    let (added_tx, added_rx) = mpsc::channel(1024);

    info!(
        "[blast_contigs:{}] downstream fanout channels created (contig2lineage/read2blastm8/updated/added)",
        db_type
    );

    let update_handle = tokio::spawn({
        let start = Instant::now();
        info!("[blast_contigs:{}] update_read_dict started", db_type);
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
                start.elapsed()
            );
            res
        }
    });
    cleanup_tasks.push(update_handle);

    let contig2lineage_task = tokio::spawn(async move {
        let start = Instant::now();
        info!("[blast_contigs:{}] contig2lineage collector started", db_type);

        let mut contig2lineage: AHashMap<String, [i32; 3]> = AHashMap::new();
        let mut lineage_stream = ReceiverStream::new(contig2lineage_rx);
        let mut total = 0usize;

        while let Some(item) = lineage_stream.next().await {
            total += 1;
            let bytes = item.to_bytes()?;
            let line = String::from_utf8_lossy(&bytes);
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                let contig = fields[0].to_string();
                let species = fields[1].parse().unwrap_or(0);
                let genus = fields[2].parse().unwrap_or(0);
                let family = fields[3].parse().unwrap_or(0);
                contig2lineage.insert(contig, [species, genus, family]);
            }

            if total % 10_000 == 0 {
                info!(
                    "[blast_contigs:{}] contig2lineage collector progress: {} records after {:?}",
                    db_type,
                    total,
                    start.elapsed()
                );
            }
        }

        info!(
            "[blast_contigs:{}] contig2lineage collector complete: {} records after {:?}",
            db_type,
            total,
            start.elapsed()
        );

        Ok::<_, anyhow::Error>(contig2lineage)
    });

    let (refined_m8_local_tx, refined_m8_local_rx) = mpsc::channel(2_000_000);
    let (refined_hit_summary_local_tx, refined_hit_summary_local_rx) = mpsc::channel(2_000_000);

    info!(
        "[blast_contigs:{}] refined local channels created (2,000,000 each)",
        db_type
    );

    let generate_m8_handle = tokio::spawn({
        let start = Instant::now();
        info!("[blast_contigs:{}] generate_m8_and_hit_summary started", db_type);
        let fut = generate_m8_and_hit_summary(
            config.clone(),
            ReceiverStream::new(updated_rx),
            ReceiverStream::new(added_rx),
            ReceiverStream::new(read2blastm8_rx),
            hit_summary_stream,
            deduped_m8_stream,
            refined_m8_local_tx.clone(),
            refined_hit_summary_local_tx.clone(),
        );

        async move {
            let res = fut.await;
            info!(
                "[blast_contigs:{}] generate_m8_and_hit_summary finished after {:?}",
                db_type,
                start.elapsed()
            );
            res
        }
    });

    let (refined_m8_tx, refined_m8_rx) = mpsc::channel(2_000_000);
    let (refined_hit_summary_tx, refined_hit_summary_rx) = mpsc::channel(2_000_000);

    let (refined_m8_for_counts_tx, refined_m8_for_counts_rx) = mpsc::channel(2_000_000);
    let (refined_hit_for_counts_tx, refined_hit_for_counts_rx) = mpsc::channel(2_000_000);

    info!(
        "[blast_contigs:{}] refined export/count split channels created",
        db_type
    );

    let m8_forward_handle = tokio::spawn({
        let mut rx = refined_m8_local_rx;
        let out_tx = refined_m8_tx.clone();
        let counts_tx = refined_m8_for_counts_tx.clone();
        async move {
            let start = Instant::now();
            let mut total = 0usize;
            let mut last_log = Instant::now();

            info!("[blast_contigs:{}] refined m8 forward task started", db_type);

            while let Some(item) = rx.recv().await {
                total += 1;

                let out_send_start = Instant::now();
                if out_tx.send(item.clone()).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] refined m8 forward: caller branch dropped at item {} after {:?}",
                        db_type,
                        total,
                        start.elapsed()
                    );
                    break;
                }
                let out_wait = out_send_start.elapsed();

                let counts_send_start = Instant::now();
                if counts_tx.send(item).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] refined m8 forward: taxon-count branch dropped at item {} after {:?}",
                        db_type,
                        total,
                        start.elapsed()
                    );
                    break;
                }
                let counts_wait = counts_send_start.elapsed();

                if out_wait >= Duration::from_millis(50)
                    || counts_wait >= Duration::from_millis(50)
                    || last_log.elapsed() >= Duration::from_secs(10)
                {
                    info!(
                        "[blast_contigs:{}] refined m8 forward progress: total={} caller_wait={:?} counts_wait={:?} elapsed={:?}",
                        db_type,
                        total,
                        out_wait,
                        counts_wait,
                        start.elapsed()
                    );
                    last_log = Instant::now();
                }
            }

            info!(
                "[blast_contigs:{}] refined m8 forward complete: total={} elapsed={:?}",
                db_type,
                total,
                start.elapsed()
            );

            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(m8_forward_handle);

    let hit_forward_handle = tokio::spawn({
        let mut rx = refined_hit_summary_local_rx;
        let out_tx = refined_hit_summary_tx.clone();
        let counts_tx = refined_hit_for_counts_tx.clone();
        async move {
            let start = Instant::now();
            let mut total = 0usize;
            let mut last_log = Instant::now();

            info!("[blast_contigs:{}] refined hit forward task started", db_type);

            while let Some(item) = rx.recv().await {
                total += 1;

                let out_send_start = Instant::now();
                if out_tx.send(item.clone()).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] refined hit forward: caller branch dropped at item {} after {:?}",
                        db_type,
                        total,
                        start.elapsed()
                    );
                    break;
                }
                let out_wait = out_send_start.elapsed();

                let counts_send_start = Instant::now();
                if counts_tx.send(item).await.is_err() {
                    warn!(
                        "[blast_contigs:{}] refined hit forward: taxon-count branch dropped at item {} after {:?}",
                        db_type,
                        total,
                        start.elapsed()
                    );
                    break;
                }
                let counts_wait = counts_send_start.elapsed();

                if out_wait >= Duration::from_millis(50)
                    || counts_wait >= Duration::from_millis(50)
                    || last_log.elapsed() >= Duration::from_secs(10)
                {
                    info!(
                        "[blast_contigs:{}] refined hit forward progress: total={} caller_wait={:?} counts_wait={:?} elapsed={:?}",
                        db_type,
                        total,
                        out_wait,
                        counts_wait,
                        start.elapsed()
                    );
                    last_log = Instant::now();
                }
            }

            info!(
                "[blast_contigs:{}] refined hit forward complete: total={} elapsed={:?}",
                db_type,
                total,
                start.elapsed()
            );

            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(hit_forward_handle);

    let taxon_count_concurrency = compute_phase_concurrency(
        &config,
        "taxon_counting",
        0.4,
        4.0,
        64,
        4,
    );
    info!(
        "[blast_contigs:{}] taxon_count_concurrency={}",
        db_type,
        taxon_count_concurrency
    );

    let taxon_count_batch_size = compute_batch_size(
        None,
        220,
        150,
        taxon_count_concurrency,
    );
    info!(
        "[blast_contigs:{}] taxon_count_batch_size={}",
        db_type,
        taxon_count_batch_size
    );

    let (refined_counts_tx, refined_counts_rx) = mpsc::channel(1024);
    info!(
        "[blast_contigs:{}] refined_counts output channel created",
        db_type
    );

    let counts_handle = tokio::spawn({
        let start = Instant::now();
        info!("[blast_contigs:{}] generate_taxon_count_json_from_m8 started", db_type);
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
                start.elapsed()
            );
            res
        }
    });

    let contig2lineage = contig2lineage_task
        .await
        .map_err(|e| PipelineError::Other(anyhow!("contig2lineage task panicked: {e}")))??;

    info!(
        "[blast_contigs:{}] contig2lineage map built with {} entries",
        db_type,
        contig2lineage.len()
    );

    let (contig_summary_tx, contig_summary_rx) = mpsc::channel(1024);
    info!(
        "[blast_contigs:{}] contig_summary output channel created",
        db_type
    );

    let contig_summary_handle = tokio::spawn({
        let start = Instant::now();
        info!("[blast_contigs:{}] generate_contig_summary_json started", db_type);
        let fut = generate_contig_summary_json(
            read2contig.clone(),
            contig2lineage,
            read_dict.clone(),
            db_type,
            duplicate_clusters.clone(),
            4,
            contig_summary_tx,
        );

        async move {
            let res = fut.await;
            info!(
                "[blast_contigs:{}] generate_contig_summary_json finished after {:?}",
                db_type,
                start.elapsed()
            );
            res
        }
    });

    info!(
        "[blast_contigs:{}] awaiting processing task joins",
        db_type
    );

    let results = try_join_all(vec![
        top_handle,
        generate_m8_handle,
        counts_handle,
        contig_summary_handle,
    ])
        .await
        .map_err(|e| PipelineError::Other(anyhow!("contig processing task panicked: {}", e)))?;

    for res in results {
        res.map_err(|e| PipelineError::Other(anyhow!("contig processing task failed: {}", e)))?;
    }

    info!(
        "[blast_contigs:{}] all processing tasks joined after {:?}",
        db_type,
        fn_start.elapsed()
    );

    let mut refined_counts = Vec::new();
    let mut rx = ReceiverStream::new(refined_counts_rx);
    let mut refined_counts_seen = 0usize;
    while let Some(item) = rx.next().await {
        refined_counts_seen += 1;
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let count: TaxonCount = serde_json::from_str(&line)?;
        refined_counts.push(count);

        if refined_counts_seen % 10_000 == 0 {
            info!(
                "[blast_contigs:{}] refined_counts decoded {} records",
                db_type,
                refined_counts_seen
            );
        }
    }
    info!(
        "[blast_contigs:{}] refined_counts complete: {} records",
        db_type,
        refined_counts_seen
    );

    let mut contig_summary = Vec::new();
    let mut rx = ReceiverStream::new(contig_summary_rx);
    let mut contig_summary_seen = 0usize;
    while let Some(item) = rx.next().await {
        contig_summary_seen += 1;
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let entry: ContigSummaryEntry = serde_json::from_str(&line)?;
        contig_summary.push(entry);

        if contig_summary_seen % 10_000 == 0 {
            info!(
                "[blast_contigs:{}] contig_summary decoded {} records",
                db_type,
                contig_summary_seen
            );
        }
    }
    info!(
        "[blast_contigs:{}] contig_summary complete: {} records",
        db_type,
        contig_summary_seen
    );

    let final_read_dict = read_dict.lock().unwrap().clone();

    info!(
        "[blast_contigs:{}] DONE after {:?} — refined_counts={}, contig_summary={}, read_dict_size={}, cleanup_tasks={}, cleanup_receivers={}, temp_files={}",
        db_type,
        fn_start.elapsed(),
        refined_counts.len(),
        contig_summary.len(),
        final_read_dict.len(),
        cleanup_tasks.len(),
        cleanup_receivers.len(),
        temp_files.len()
    );

    Ok((
        final_read_dict,
        refined_counts,
        contig_summary,
        refined_m8_rx,
        refined_hit_summary_rx,
        top_for_caller_rx,
        cleanup_tasks,
        cleanup_receivers,
        temp_files,
    ))
}


/// Extract original read ID from refined taxid-annotated FASTA header
/// # Arguments
/// * `id` -header id line

/// # Returns
/// result of id and reead index
fn extract_original_id_from_taxid_fasta(id: &str) -> Result<(String, usize)> {
    let mut header = id.to_string();
    let mut read_index = 0;

    if header.ends_with("/1") {
        read_index = 0;
        header.truncate(header.len() - 2);
    } else if header.ends_with("/2") {
        read_index = 1;
        header.truncate(header.len() - 2);
    }

    let parts: Vec<&str> = header.split(':').collect();
    let nt_pos = parts
        .iter()
        .position(|&p| p == "NT")
        .ok_or_else(|| anyhow!("No 'NT' tag in taxid header: {}", id))?;

    if nt_pos + 2 >= parts.len() {
        return Err(anyhow!("Malformed header: not enough fields after NT"));
    }

    let original_id = parts[(nt_pos + 2)..].join(":");
    let full_id = if read_index == 0 {
        format!("{}/1", original_id)
    } else {
        format!("{}/2", original_id)
    };

    Ok((full_id, read_index))
}


/// Build R1 and R2 header sets from the taxid-annotated FASTA stream
/// # Arguments
/// * `taxid_mapped_stream` stream of

/// # Returns
/// result of id and reead index
async fn build_nonhost_header_sets(
    mut taxid_mapped_stream: ReceiverStream<ParseOutput>,
    duplicate_clusters: Option<Arc<DuplicateClusters>>,
) -> Result<(HashSet<String>, Option<HashSet<String>>)> {
    tokio::task::spawn_blocking(move || {
        let mut r1_set: HashSet<String> = HashSet::new();
        let mut r2_set: HashSet<String> = HashSet::new();

        while let Some(item) = futures::executor::block_on(taxid_mapped_stream.next()) {
            if let ParseOutput::Fasta(record) = item {
                if let Ok((full_id, read_index)) = extract_original_id_from_taxid_fasta(record.id()) {
                    let target_set = if read_index == 0 { &mut r1_set } else { &mut r2_set };

                    if READ_COUNTING_MODE == ReadCountingMode::CountAll {
                        if let Some(clusters) = &duplicate_clusters {
                            let mut rep_key = full_id.clone();
                            if !clusters.contains_key(&rep_key) {
                                rep_key = rep_key
                                    .trim_end_matches("/1")
                                    .trim_end_matches("/2")
                                    .to_string();
                            }

                            if let Some(cluster) = clusters.get(&rep_key) {
                                for member in &cluster.members {
                                    target_set.insert(member.clone());
                                }
                                continue;
                            }
                        }
                    }
                    target_set.insert(full_id);
                }
            }
        }

        let r2_opt = if r2_set.is_empty() { None } else { Some(r2_set) };
        Ok((r1_set, r2_opt))
    })
        .await?
}



pub async fn generate_nonhost_fastq_from_files(
    config: Arc<RunConfig>,
    original_r1_path: PathBuf,
    original_r2_path: Option<PathBuf>,
    taxid_mapped_stream: ReceiverStream<ParseOutput>,
    duplicate_clusters: Option<Arc<DuplicateClusters>>,
    output_dir: PathBuf,
) -> Result<(PathBuf, Option<PathBuf>)> {
    let out_r1 = output_dir.join("nonhost_R1.fastq");
    let out_r2 = original_r2_path.is_some().then(|| output_dir.join("nonhost_R2.fastq"));

    // Build header sets from taxid FASTA
    let (r1_headers, r2_headers_opt) = build_nonhost_header_sets(taxid_mapped_stream, duplicate_clusters).await?;

    info!(
        "Non-host header sets built: R1={} records, R2={:?}",
        r1_headers.len(),
        r2_headers_opt.as_ref().map(|s| s.len())
    );

    let r1_headers = Arc::new(r1_headers);
    let r2_headers = r2_headers_opt.map(Arc::new);

    // Stream and filter R1
    let (r1_rx, _r1_stats_task) = read_fastq(
        original_r1_path.clone(),
        None,
        PairingMode::Relaxed,
        None,
        u64::MAX,
        None,
        None,
        "generate_nonhost_fastq_from_files_r1",
        &config,

    )?;
    let r1_input_stream = ReceiverStream::new(r1_rx);
    let r1_filtered = filter_fastq_to_bytes_stream(r1_input_stream, r1_headers).await;
    let write_r1 = write_byte_stream_to_file(&out_r1, r1_filtered, config.clone(), StreamDataType::IlluminaFastq, "nonhost_r1").await?;

    // Stream and filter R2 if present
    let write_r2 = if let (Some(r2_path), Some(r2_headers)) = (original_r2_path, r2_headers) {
        let (r2_rx, _r2_stats_task) = read_fastq(
            r2_path.clone(),
            None,
            PairingMode::Relaxed,
            None,
            u64::MAX,
            None,
            None,
            "generate_nonhost_fastq_from_files_r2",
            &config,
        )?;
        let r2_input_stream = ReceiverStream::new(r2_rx);
        let r2_filtered = filter_fastq_to_bytes_stream(r2_input_stream, r2_headers).await;
        Some(write_byte_stream_to_file(&out_r2.clone().unwrap(), r2_filtered,config.clone(), StreamDataType::IlluminaFastq, "nonhost_r2").await?)
    } else {
        None
    };

    // Await writes
    write_r1.await??;
    if let Some(task) = write_r2 {
        task.await??;
    }

    info!(
        "Non-host FASTQs generated: {} {:?}",
        out_r1.display(),
        out_r2.as_ref().map(|p| p.display())
    );

    Ok((out_r1, out_r2))
}

/// Parallel preload of NR hit-summary → DashMap (exact Python semantics, >50k lines/s on EPYC)
async fn preload_nr_alignments_parallel(
    config: Arc<RunConfig>,
    preload_rx: mpsc::Receiver<ParseOutput>,
    nr_alignment_per_read: Arc<DashMap<String, SpeciesAlignmentResults, AHashRandomState>>,
) -> Result<()> {
    let concurrency = compute_phase_concurrency(
        &config,
        "nr_preload_alignments",
        0.3,   // tiny RAM per worker
        2.0,   // almost pure parse+insert
        128,   // safe even on 256-core EPYC
        8,
    );

    let batch_size = compute_batch_size(None, 180, 200, concurrency);

    let (job_tx, job_rx) = mpsc::channel::<Vec<String>>(concurrency);
    let shared_rx = Arc::new(tokio::sync::Mutex::new(job_rx));

    let mut workers = Vec::with_capacity(concurrency);
    for _ in 0..concurrency {
        let rx = Arc::clone(&shared_rx);
        let map = Arc::clone(&nr_alignment_per_read);

        let handle = tokio::spawn(async move {
            loop {
                let batch = {
                    let mut guard = rx.lock().await;
                    guard.recv().await
                };
                let Some(batch) = batch else { break };

                batch.par_iter().for_each(|line| {
                    let trimmed = line.trim_end();
                    if trimmed.is_empty() { return; }

                    let fields: Vec<&str> = trimmed.split('\t').collect();
                    if fields.len() < 10 { return; }

                    let read_id = fields[0].to_string();
                    let contig_taxid = fields[9].parse::<Taxid>().ok();
                    let read_taxid  = fields[3].parse::<Taxid>().ok();

                    map.insert(
                        read_id,
                        SpeciesAlignmentResults { contig: contig_taxid, read: read_taxid },
                    );
                });
            }
            Ok::<(), anyhow::Error>(())
        });
        workers.push(handle);
    }

    // Producer (same lockstep style as all your other high-throughput stages)
    let mut batch = Vec::with_capacity(batch_size);
    let mut stream = ReceiverStream::new(preload_rx);
    while let Some(item) = stream.next().await {
        if let ParseOutput::Bytes(b) = item {
            let line = String::from_utf8_lossy(&b).trim_end().to_string();
            if !line.is_empty() {
                batch.push(line);
            }
        }
        if batch.len() >= batch_size {
            let _ = job_tx.send(std::mem::take(&mut batch)).await;
        }
    }
    if !batch.is_empty() {
        let _ = job_tx.send(batch).await;
    }
    drop(job_tx);

    for h in workers {
        h.await??;
    }

    info!("Preloaded {} NR alignments (parallel, {} workers)", nr_alignment_per_read.len(), concurrency);
    Ok(())
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

async fn run_mmseqs_step(
    config: Arc<RunConfig>,
    args: Vec<String>,
    label: &'static str,
) -> Result<(), PipelineError> {
    let start = Instant::now();

    let mmseqs_stderr_log = config.out_dir.join("mmseqs").join(format!("{}.stderr.log", label));

    // spawn_cmd now handles numactl internally via prepend_numactl_if_beneficial
    let (mut child, stderr_task) = spawn_cmd(
        config.clone(),
        MMSEQS_TAG,
        args,
        config.args.verbose,
        Some(mmseqs_stderr_log)
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MMSEQS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let wait_task = tokio::spawn(async move {
        let status = child.wait().await?;
        if !status.success() {
            return Err(anyhow!("mmseqs {} exited with code {:?}", label, status.code()));
        }
        Ok::<(), anyhow::Error>(())
    });

    wait_task
        .await
        .map_err(|e| PipelineError::Other(anyhow!("mmseqs {} wait task panicked: {}", label, e)))??;

    stderr_task
        .await
        .map_err(|e| PipelineError::Other(anyhow!("mmseqs {} stderr task panicked: {}", label, e)))??;

    let elapsed = start.elapsed();
    info!("[mmseqs:{}] completed in {:?}", label, elapsed);

    Ok(())
}

async fn start_gpuserver(
    db: &PathBuf,
    verbose: bool,
) -> Result<Child, PipelineError> {
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

    let createdb_args = generate_cli(MMSEQS_TAG, &config, Some(&createdb_cfg))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MMSEQS_TAG.to_string(),
            error: e.to_string(),
        })?;

    info!("[mmseqs:{}] createdb args: {:?}", label, createdb_args);
    run_mmseqs_step(config.clone(), createdb_args, "createdb").await?;

    if !query_db.exists() {
        return Err(PipelineError::Other(anyhow!(
            "mmseqs createdb did not produce queryDB: {}", query_db.display()
        )));
    }
    let meta = fs::metadata(&query_db).await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    if meta.len() < 10_000 {
        return Err(PipelineError::Other(anyhow!(
            "queryDB looks too small ({} bytes) — createdb probably failed", meta.len()
        )));
    }

    info!("Query DB created successfully: {} ({} bytes)", query_db.display(), meta.len());

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
    let search_cfg = MmseqsConfig {
        subcommand: MmseqsSubcommand::Search,
        backend,
        input: Some(query_db.clone()),
        target_db: Some(target_db.clone()),
        result_db: Some(result_db.clone()),
        output: None,
        tmp_dir: Some(tmp_dir.clone()),
        threads: Some(config.mmseqs_threads()),

        sensitivity: match backend {
            MmseqsBackend::Cpu => Some("5.7".to_string()),
            MmseqsBackend::Gpu => None,
        },

        search_type: Some("3".to_string()),
        max_seqs: Some(match backend {
            MmseqsBackend::Cpu => "1000".to_string(),
            MmseqsBackend::Gpu => "3000".to_string(),
        }),
        prefilter_mode: Some("0".to_string()),
        db_load_mode: Some("2".to_string()),
        alignment_mode: Some("3".to_string()),

        index_subset: None,
        format_output: None,
        cuda_visible_devices: None,
        option_fields: HashMap::from([
            ("-e".to_string(), Some("0.001".to_string())),
            ("--min-seq-id".to_string(), Some("0.25".to_string())),
        ]),
        gpu_server: backend == MmseqsBackend::Gpu,
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

    let convert_args = generate_cli(MMSEQS_TAG, &config, Some(&convert_cfg))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MMSEQS_TAG.to_string(),
            error: e.to_string(),
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

pub async fn mmseqs_non_host_align(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    backend: MmseqsBackend,
) -> Result<(
    tokio::sync::mpsc::Receiver<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<tokio::sync::oneshot::Receiver<Result<(), anyhow::Error>>>,
    Vec<TempDir>,
), PipelineError> {
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

    write_combined_fastq(
        r1_path,
        r2_path_opt.clone(),
        &combined_fastq,
    )
        .await?;

    let merged_m8 = mmseqs_fastq_to_m8_file(
        config.clone(),
        combined_fastq,
        &temp_dir,
        "nr_merged",
        backend,
    )
        .await?;

    let m8_file = fs::File::open(&merged_m8)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let m8_rx = parse_lines(m8_file, &config, StreamDataType::JustBytes)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok((m8_rx, cleanup_tasks, cleanup_receivers, vec![temp_dir]))
}

/// Run function for Short Read mNGS pipelines
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
///
/// # Returns
/// Result<(), PipelineError>
/// Run function for Short Read mNGS pipelines
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
///
/// # Returns
/// Result<(), PipelineError>
pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();
    let mut cleanup_tasks: Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>> = Vec::new();
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

    check_versions(versions_vec, &out_dir.clone(), &config)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Check required files
    let host_bowtie2_index: String = config.args.host_bowtie2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_bowtie2_index is required".to_string()))?;

    let host_hisat2_index: String = config.args.host_hisat2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_hisat2_index is required".to_string()))?;

    let (file1_path, file2_path, sample_base_buf, sample_base) = validate_file_inputs(&config, &cwd).await?;
    let paired = file2_path.is_some();

    debug!("file1 path: {:?}  sample base buf: {:?}  sample base: {:?}", file1_path.display(), sample_base_buf.display(), sample_base);

    let taxonomy_handle = tokio::spawn(load_lineage_and_acc2tax_maps(config.clone()));

    // Input Validation
    let (val_out_stream, validate_cleanup_tasks, validate_cleanup_receivers, raw_count_task, val_count_task) = validate_input(
        config.clone(),
        file1_path.clone(),
        file2_path.clone(),
        sample_base_buf.clone(),
        &out_dir,
    ).await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // *******************
    // Host Filtering Stage
    // *******************

    // ERCC bt2 filtering and count
    let ercc_bt2_index_path = bowtie2_index_prep(&config.args.ercc_bowtie2_index, &cwd)?;
    let ercc_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
    let (ercc_bt2_out_stream, ercc_count_rx, ercc_bt2_cleanup_tasks, ercc_bt2_cleanup_receivers, ercc_bt2_bam_write_handle, _ercc_bt2_bam_path, ercc_bt2_temp_dirs) = bowtie2_filter_stream(
        config.clone(),
        val_out_stream,
        ercc_bt2_index_path,
        paired,
        ercc_bt2_options,
        out_dir.join("ercc_bt2_aligned_sorted.bam"),
    ).await?;
    cleanup_tasks.extend(ercc_bt2_cleanup_tasks);
    cleanup_receivers.extend(ercc_bt2_cleanup_receivers);
    final_temp_dirs.extend(ercc_bt2_temp_dirs);

    // QC with fastp
    let (qc_fastp_out_stream, qc_cleanup_tasks, qc_cleanup_receivers, qc_count_result_rx) = fastp_qc(
        config.clone(),
        paired,
        ercc_bt2_out_stream,
    ).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);
    cleanup_receivers.extend(qc_cleanup_receivers);

    // Split for Kallisto and Bowtie2
    let (kallisto_streams, kallisto_router_handle) = fanout_to_channels(
        qc_fastp_out_stream,
        2,
        "kallisto_split",
        &config,
        StreamDataType::JustBytes
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    // track router task
    cleanup_tasks.push(kallisto_router_handle);

    let mut kallisto_streams_iter = kallisto_streams.into_iter();

    let kallisto_bypass_stream = ReceiverStream::new(
        kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );

    let kallisto_stream = ReceiverStream::new(
        kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?
    );


    let (kallisto_ercc_tx, kallisto_ercc_rx, kallisto_cleanup_tasks, kallisto_cleanup_receivers, kallisto_exit_task) = kallisto_quant(
        config.clone(),
        kallisto_stream,
        out_dir.clone(),
        paired,
        sample_base_buf.clone(),
    ).await?;
    cleanup_tasks.extend(kallisto_cleanup_tasks);
    cleanup_receivers.extend(kallisto_cleanup_receivers);

    // Host filtering: bt2
    let host_bt2_index_path = bowtie2_index_prep(host_bowtie2_index, &cwd)?;
    let host_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
    let (host_bt2_r1, host_bt2_r2, host_bt2_count_rx, host_bt2_cleanup_tasks, host_bt2_cleanup_receivers, host_bt2_bam_write_handle, host_bt2_bam_path, host_bt2_temp_dirs) = bowtie2_filter_files(
        config.clone(),
        kallisto_bypass_stream,
        host_bt2_index_path,
        paired,
        host_bt2_options,
        out_dir.join("host_bt2_aligned_sorted.bam"),
    ).await?;
    cleanup_tasks.extend(host_bt2_cleanup_tasks);
    cleanup_receivers.extend(host_bt2_cleanup_receivers);
    host_filter_temp_dirs.extend(host_bt2_temp_dirs);


    //host filtering hisat2
    let host_hisat2_index_path = hisat2_index_prep(host_hisat2_index, &cwd)?;
    let hisat2_options = HashMap::from([]);
    let (host_hisat2_out_stream, host_hisat2_count_rx, mut host_hisat2_cleanup_tasks, mut host_hisat2_cleanup_receivers, hisat_temp_dirs) = hisat2_filter(
        config.clone(),
        host_bt2_r1,
        host_bt2_r2,
        host_hisat2_index_path,
        paired,
        hisat2_options,
        None,
        4
    )
        .await?;
    cleanup_tasks.append(&mut host_hisat2_cleanup_tasks);
    cleanup_receivers.append(&mut host_hisat2_cleanup_receivers);
    host_filter_temp_dirs.extend(hisat_temp_dirs);


    // If host is no huma, run an additional filter stage using a human reference
    let mut optional_human_bam_write_handle: Option<JoinHandle<Result<()>>> = None;
    let mut optional_human_bam_path: Option<PathBuf> = None;

    let post_filter_stream = if config.args.human_host {
        host_hisat2_out_stream
    } else {
        let human_bowtie2_index: String = config.args.human_bowtie2_index.clone();
        let human_hisat2_index: String = config.args.human_hisat2_index.clone();

        let human_bt2_index_path = bowtie2_index_prep(human_bowtie2_index, &cwd)?;
        let human_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);

        let (
            human_bt2_r1,
            human_bt2_r2,
            _human_bt2_count_rx,
            human_bt2_cleanup_tasks,
            human_bt2_cleanup_receivers,
            human_bt2_bam_write_handle,
            human_bt2_bam_path,
            human_bt2_temp_dirs
        ) = bowtie2_filter_files(
            config.clone(),
            host_hisat2_out_stream,
            human_bt2_index_path,
            paired,
            human_bt2_options,
            out_dir.join("human_bt2_aligned_sorted.bam"),
        )
            .await?;

        cleanup_tasks.extend(human_bt2_cleanup_tasks);
        cleanup_receivers.extend(human_bt2_cleanup_receivers);
        host_filter_temp_dirs.extend(human_bt2_temp_dirs);

        // Save the BAM write handle and path for later awaiting / insert stats
        optional_human_bam_write_handle = Some(human_bt2_bam_write_handle);
        optional_human_bam_path = Some(human_bt2_bam_path);

        let human_hisat2_index_path = hisat2_index_prep(human_hisat2_index, &cwd)?;
        let hisat2_options = HashMap::new();

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
            hisat2_options,
            None,
            4,
        )
            .await?;

        cleanup_tasks.append(&mut human_hisat2_cleanup_tasks);
        cleanup_receivers.append(&mut human_hisat2_cleanup_receivers);
        host_filter_temp_dirs.extend(human_hisat_temp_dirs);

        human_hisat2_out_stream
    };

    // let post_filter_monitored = monitor_stream(post_filter_stream, "Post_Hisat2", Duration::from_secs(5));


    let (pre_dedup_parsed_stream, parse_task) = parse_byte_stream_to_fastq(
        post_filter_stream.into_inner(),
        config.base_buffer_size,
        config.args.stall_threshold,
    ).await?;

    cleanup_tasks.push(parse_task);


    let (dedup_stream, dedup_count_rx, duplicate_clusters, mut dedup_cleanup_tasks, mut dedup_cleanup_receivers) = dedup(
        config.clone(),
        pre_dedup_parsed_stream,
        paired,
        Some(70), // Prefix length for deduplication. Hardcoded for now
        out_dir.clone(),
    ).await?;
    cleanup_tasks.append(&mut dedup_cleanup_tasks);
    cleanup_receivers.append(&mut dedup_cleanup_receivers);

    let uniques_count = dedup_count_rx.await?;
    let unique_reads = uniques_count * if paired { 2 } else { 1 };
    info!("Uniques count after dedup (pairs counted separately): {}", unique_reads);

    // Separate subsample (weighted by cluster sizes; correctness: full stream propagation, no silent drops via explicit send/await)
    let (subsampled_stream, subsample_count_rx, subsample_send_task) = subsample_uniform(
        config.clone(),
        dedup_stream.into_inner(),
        config.args.max_subsample as u64,
        paired,
    ).await?;
   cleanup_tasks.push(subsample_send_task);


    let subsample_count = subsample_count_rx.await?;
    let subsample_reads = subsample_count * if paired { 2 } else { 1 };
    info!("Subsampled reads (pairs counted separately): {}", subsample_reads);


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
    let (r1_rx, r2_rx_opt, deinterleave_handle) = deinterleave_fastq_stream(
        subsampled_stream,
        paired,
        config.base_buffer_size * 4,
    )
        .await
        .map_err(|e| PipelineError::Other(e))?;

    let non_host_temp_dir = choose_temp_dir(
        config.input_size * 2,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        true,
    ).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Compute paths after temp_dir but before moving rx (only borrows with as_ref())
    let non_host_r1_path = non_host_temp_dir.path().join("nonhost_R1.fastq");
    let non_host_r2_path_opt = r2_rx_opt.as_ref().map(|_| non_host_temp_dir.path().join("nonhost_R2.fastq"));

    // Now create monitored streams (this moves r1_rx and r2_rx_opt)
    // let r1_monitored = monitor_stream(
    //     ReceiverStream::new(r1_rx),
    //     "Deinterleave_R1",
    //     Duration::from_secs(5),
    // );
    // let r2_monitored_opt = r2_rx_opt.map(|rx| {
    //     monitor_stream(
    //         ReceiverStream::new(rx),
    //         "Deinterleave_R2",
    //         Duration::from_secs(5),
    //     )
    // });

    info!("Checkpoint: deinterleaving and writing non-host R1/R2 to temp (forces upstream completion). Writing to: {:?}", non_host_temp_dir);

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
    let (r1_write_task, r2_write_task_opt) = tokio::try_join!(
    r1_write_task_fut,
    async {
        match r2_write_task_fut_opt {
            Some(fut) => fut.await.map(Some),
            None => Ok(None),
        }
    }
)?;

    // then wait for both background writes and deinterleave together
    tokio::try_join!(
    async { r1_write_task.await?; Ok::<(), anyhow::Error>(()) },
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

    info!("Checkpoint complete. Files: r1:{:?}   r2:{:?}", non_host_r1_path, non_host_r2_path_opt);


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



    let (non_host_mm2_out_stream, mut non_host_mm2_cleanup_tasks, mut non_host_mm2_cleanup_receivers) = minimap2_non_host_align(
        config.clone(),
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
    ).await?;

    cleanup_tasks.append(&mut non_host_mm2_cleanup_tasks);
    cleanup_receivers.append(&mut non_host_mm2_cleanup_receivers);


    let m8_file_path = out_dir.join(rename_file_path(&sample_base_buf, None, Some("m8"), "."));

    let (m8_stream, mut m8_cleanup_tasks, mut m8_cleanup_receivers) = paf_to_m8(
        config.clone(),
        non_host_mm2_out_stream,
        m8_file_path,
    )
        .await?;
    cleanup_tasks.append(&mut m8_cleanup_tasks);
    cleanup_receivers.append(&mut m8_cleanup_receivers);

    // ────────────────────────────────────────────────────────────────
    // Sort m8 by read ID before call_hits_m8
    // This guarantees consecutive lines per read → streaming group-by
    // ────────────────────────────────────────────────────────────────
    let m8_sorted_start = Instant::now();
    let m8_sorted = sort_m8_by_read_id(
        config.clone(),
        m8_stream,
        "nt",                    // label for logging + temp files
    ).await?;
    info!(
    "[run] sort_m8_by_read_id(nt) completed after {:?}",
    m8_sorted_start.elapsed()
);

    let (lineage_map, acc2taxid_map) = taxonomy_handle.await??;

    let nt_concurrency = compute_phase_concurrency(
        &config,
        "call_hits_nt",
        1.0,
        3.5,
        64,
        16,
    );
    info!("call hits nt concurrency {}", nt_concurrency);

    let nt_call_hits_start = Instant::now();
    let (nt_call_stream, nt_call_summary_stream, mut call_cleanup_tasks, mut call_cleanup_receivers) = call_hits_m8(
        config.clone(),
        m8_sorted,                    // ←←← NOW SORTED
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        36,
        nt_concurrency,
        "nt".to_string(),
    ).await?;
    info!(
    "[run] call_hits_m8(nt) returned after {:?}",
    nt_call_hits_start.elapsed()
);

    cleanup_tasks.append(&mut call_cleanup_tasks);
    cleanup_receivers.append(&mut call_cleanup_receivers);

    info!("[run] wrapping NT outputs with monitor_stream");

    let nt_call_stream = monitor_stream(
        nt_call_stream,
        "nt_call_stream_from_call_hits_m8",
        Duration::from_secs(10),
    );

    let nt_call_summary_stream = monitor_stream(
        nt_call_summary_stream,
        "nt_call_summary_stream_from_call_hits_m8",
        Duration::from_secs(10),
    );

    let nt_split_start = Instant::now();
    info!("[run] starting fanout_to_channels for NT call m8 (2 private channels)");

    let (nt_call_rxs, nt_call_router) = fanout_to_channels(
        nt_call_stream,                    // already the monitored stream from call_hits_m8
        2, // generous buffer for m8 lines
        "nt_call",
        &config,
        StreamDataType::JustBytes
    ).await?;

    cleanup_tasks.push(nt_call_router);    // router returns the correct JoinHandle type

    info!(
        "[run] fanout_to_channels(nt_call) ready after {:?} with 2 private channels",
        nt_split_start.elapsed()
    );

    let mut nt_call_rxs_iter = nt_call_rxs.into_iter();
    let nt_call_stream = ReceiverStream::new(nt_call_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nt_blast_stream = ReceiverStream::new(nt_call_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);


    let nt_call_stream = monitor_stream(
        nt_call_stream,
        "nt_call_stream_to_generate_taxon_counts",
        Duration::from_secs(15),
    );

    let nt_blast_stream = monitor_stream(
        nt_blast_stream,
        "nt_blast_stream_to_blast_contigs",
        Duration::from_secs(15),
    );

    let nt_summary_split_start = Instant::now();
    info!("[run] starting fanout_to_channels for NT summary (6 private channels)");

    let (nt_summary_rxs, nt_summary_router) = fanout_to_channels(
        nt_call_summary_stream,
        6,
        "nt_call_summary",
        &config,
        StreamDataType::JustBytes
    ).await?;

    cleanup_tasks.push(nt_summary_router);     // now the correct type

    info!(
        "[run] fanout_to_channels(nt_call_summary) ready after {:?} with 6 private channels",
        nt_summary_split_start.elapsed()
    );

    let mut nt_summary_rxs_iter = nt_summary_rxs.into_iter();
    let nt_summary_taxon_stream       = ReceiverStream::new(nt_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nt_summary_hit_stream         = ReceiverStream::new(nt_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nt_initial_stream             = ReceiverStream::new(nt_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nt_blast_hit_stream           = ReceiverStream::new(nt_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nt_hit_summary_for_refined    = ReceiverStream::new(nt_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nt_hit_summary_for_taxid      = ReceiverStream::new(nt_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);

    info!("[run] NT split wiring complete");
    info!("[run] nt_summary_taxon_stream -> generate_taxon_counts");
    info!("[run] nt_summary_hit_stream -> summarize_hits");
    info!("[run] nt_initial_stream -> initial_taxid_fasta");
    info!("[run] nt_blast_hit_stream -> blast_contigs");

    let nt_summary_taxon_stream = monitor_stream(
        nt_summary_taxon_stream,
        "nt_summary_taxon_stream_to_generate_taxon_counts",
        Duration::from_secs(10),
    );

    let nt_summary_hit_stream = monitor_stream(
        nt_summary_hit_stream,
        "nt_summary_hit_stream_to_summarize_hits",
        Duration::from_secs(10),
    );

    let nt_initial_stream = monitor_stream(
        nt_initial_stream,
        "nt_initial_stream_to_initial_taxid_fasta",
        Duration::from_secs(15),
    );

    let nt_blast_hit_stream = monitor_stream(
        nt_blast_hit_stream,
        "nt_blast_hit_stream_to_blast_contigs",
        Duration::from_secs(15),
    );

    info!("Launching downstream consumers for NT results:");
    info!("  -> nt_call_stream (to generate_taxon_counts + nt_map_task)");
    info!("  -> nt_blast_stream (to blast_contigs later)");
    info!("  -> nt_summary_hit_stream (to summarize_hits)");
    info!("  -> nt_summary_taxon_stream (to generate_taxon_counts)");

    let nt_hit_summary_handle = tokio::spawn({
        let config = config.clone();
        let nt_summary_hit_stream = nt_summary_hit_stream;
        let duplicate_clusters = duplicate_clusters.clone();
        async move {
            let start = Instant::now();
            info!("[run] summarize_hits(nt) started");
            let res = summarize_hits(
                config.clone(),
                nt_summary_hit_stream,
                duplicate_clusters,
                0,
            ).await;
            info!("[run] summarize_hits(nt) finished after {:?}", start.elapsed());
            res
        }
    });

    let nt_counts_task = tokio::spawn({
        let config = config.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let nt_call_stream = nt_call_stream;
        let nt_summary_taxon_stream = nt_summary_taxon_stream;
        async move {
            let start = Instant::now();
            info!("[run] generate_taxon_counts(NT) started");
            info!("[run] generate_taxon_counts(NT) awaiting nt_call_stream + nt_summary_taxon_stream");
            let res = generate_taxon_counts(
                config,
                nt_call_stream,
                nt_summary_taxon_stream,
                duplicate_clusters,
                should_keep_filter,
                "NT".to_string(),
                None,
            ).await;
            info!("[run] generate_taxon_counts(NT) finished after {:?}", start.elapsed());
            res
        }
    });


    // Temporary: Skip Diamond by providing empty outputs
    let (dummy_tx, dummy_rx) = mpsc::channel::<ParseOutput>(1);
    drop(dummy_tx); // Immediately drop sender to create an empty stream
    let non_host_m8_stream = dummy_rx;
    // let mut non_host_diamond_cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    // let mut non_host_diamond_cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();

    // Diamond or MMseqs2 non_host alignment
    // let (non_host_m8_stream, mut non_host_cleanup_tasks, mut non_host_cleanup_receivers, non_host_temp_dirs) =
    //     match config.alignment_backend {
    //         NRAlignmentBackend::Diamond => {
    //             diamond_non_host_align(
    //                 config.clone(),
    //                 non_host_r1_path.clone(),
    //                 non_host_r2_path_opt.clone(),
    //             ).await?
    //         }
    //         NRAlignmentBackend::MmseqsCpu => {
    //             mmseqs_non_host_align(
    //                 config.clone(),
    //                 non_host_r1_path.clone(),
    //                 non_host_r2_path_opt.clone(),
    //                 MmseqsBackend::Cpu,
    //             ).await?
    //         }
    //         NRAlignmentBackend::MmseqsGpu => {
    //             mmseqs_non_host_align(
    //                 config.clone(),
    //                 non_host_r1_path.clone(),
    //                 non_host_r2_path_opt.clone(),
    //                 MmseqsBackend::Gpu,
    //             ).await?
    //         }
    //     };

    // cleanup_tasks.append(&mut non_host_cleanup_tasks);
    // cleanup_receivers.append(&mut non_host_cleanup_receivers);
    // final_temp_dirs.extend(non_host_temp_dirs);

    let nr_concurrency = compute_phase_concurrency(
        &config,
        "call_hits_nr",
        1.0,           // ~1 GB per thread max (mostly transient)
        3.5,
        64,            // higher cap for CPU-bound phase
        16,            // min for meaningful parallelism
    );

    info!("call hits nr concurrency {}", nr_concurrency);

    // ────────────────────────────────────────────────────────────────
    // Sort NR m8 by read ID before call_hits_m8
    // Guarantees consecutive lines per read → enables true streaming group-by
    // ────────────────────────────────────────────────────────────────
    let nr_sort_start = Instant::now();
    let nr_m8_sorted = sort_m8_by_read_id(
        config.clone(),
        ReceiverStream::new(non_host_m8_stream),
        "nr",
    ).await?;
    info!(
    "[run] sort_m8_by_read_id(nr) completed after {:?}",
    nr_sort_start.elapsed()
);

    let (nr_call_stream, nr_call_summary_stream, mut nr_call_cleanup_tasks, mut nr_call_cleanup_receivers) =
        call_hits_m8(
            config.clone(),
            nr_m8_sorted,
            sample_base_buf.clone(),
            lineage_map.clone(),
            acc2taxid_map.clone(),
            should_keep_filter.clone(),
            0,
            nr_concurrency,
            "nr".to_string(),
        ).await?;
    cleanup_tasks.append(&mut nr_call_cleanup_tasks);
    cleanup_receivers.append(&mut nr_call_cleanup_receivers);

    let nr_split_start = Instant::now();
    info!("[run] starting fanout_to_channels for NT call m8 (2 private channels)");

    let (nr_call_rxs, nr_call_router) = fanout_to_channels(
        nr_call_stream,                    // already the monitored stream from call_hits_m8
        2,
        "nr_call",
        &config,
        StreamDataType::JustBytes
    ).await?;

    cleanup_tasks.push(nr_call_router);    // router returns the correct JoinHandle type

    info!(
    "[run] fanout_to_channels(nr_call) ready after {:?} with 2 private channels",
    nr_split_start.elapsed()
);

    let mut nr_call_rxs_iter = nr_call_rxs.into_iter();
    let nr_call_stream = ReceiverStream::new(nr_call_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_blast_stream = ReceiverStream::new(nr_call_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);

    let nr_summary_split_start = Instant::now();
    info!("[run] starting fanout_to_channels for NR summary (6 private channels)");
    let (nr_summary_rxs, nr_summary_router) = fanout_to_channels(
        nr_call_summary_stream,
        6,
        "nr_call_summary",
        &config,
        StreamDataType::JustBytes
    ).await?;

    cleanup_tasks.push(nr_summary_router);     // now the correct type

    info!(
     "[run] fanout_to_channels(nr_call_summary) ready after {:?} with 6 private channels",
     nr_summary_split_start.elapsed()
 );

    let mut nr_summary_rxs_iter = nr_summary_rxs.into_iter();
    let nr_summary_taxon_stream       = ReceiverStream::new(nr_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_summary_hit_stream         = ReceiverStream::new(nr_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_initial_stream             = ReceiverStream::new(nr_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_blast_hit_stream           = ReceiverStream::new(nr_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_hit_summary_for_refined    = ReceiverStream::new(nr_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);
    let nr_hit_summary_for_taxid      = ReceiverStream::new(nr_summary_rxs_iter.next().ok_or(PipelineError::EmptyStream)?);

    let nr_hit_summary_handle = tokio::spawn(summarize_hits(
        config.clone(),
        nr_summary_hit_stream,
        duplicate_clusters.clone(),
        0,
    ));

    let nr_counts_task = tokio::spawn({
        let config = config.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        async move {
            generate_taxon_counts(
                config, // RunConfig
                nr_call_stream, //m8_stream
                nr_summary_taxon_stream, //summary_stream
                duplicate_clusters, //duplicate_clusters
                should_keep_filter, //should_keep_filter
                "NR".to_string(), // count_type
                None, //source_count type, only relevantif/when i later implement the merged path (NT + NR together) pass Some("NT".to_string()) or Some("NR".to_string()) on the respective streams.
            ).await
        }
    });


    let annot_concurrency = compute_phase_concurrency(
        &config,
        "generate_annotated_fasta",
        0.4,
        4.0,
        128,
        8,
    );

    // Read non-host R1/R2 as proper parsed Fastq records (this fixes the Non-FASTQ error)
    let (interleaved_rx, read_stats_task) = read_fastq(
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        PairingMode::Relaxed,
        None,                    // technology
        u64::MAX,
        None,
        None,
        "run - Read non-host R1/R2 as proper parsed Fastq records",
        &config
    ).map_err(|e| PipelineError::Other(e))?;

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
        collect_hit_summary_to_accession_map_concurrent(
            config.clone(),
            nt_initial_stream
        ),
        collect_hit_summary_to_accession_map_concurrent(
            config.clone(),
            nr_initial_stream
        )
    )?;

    // Generate the annotated FASTA stream
    let (annotated_rx, unidentified_rx, unique_unidentified_rx, mut annot_tasks, mut annot_rxs) =
        generate_annotated_fasta_stream(
            config.clone(),
            interleaved_stream,
            duplicate_clusters.clone(),
            nt_map,                     // small contig → NT accession map
            nr_map,                     // small contig → NR accession map
            annot_concurrency,
        )
            .await
            .map_err(|e| PipelineError::Other(anyhow!("Annotated FASTA from files failed: {}", e)))?;

    cleanup_tasks.append(&mut annot_tasks);
    cleanup_receivers.append(&mut annot_rxs);

    let assembly_dir = out_dir.join("assembly");

    let annotated_path = assembly_dir.join("annotated_merged.fa");
    let unidentified_path = assembly_dir.join("unidentified.fa");
    let unique_unidentified_path = assembly_dir.join("unique_unidentified.fa");

    info!("[run] starting fanout_to_channels for annotated FASTA (2 private channels)");

    let (annotated_rxs, annotated_router) = fanout_to_channels(
        ReceiverStream::new(annotated_rx),
        2,
        "initial_annotated",
        &config,
        StreamDataType::JustBytes
    ).await?;

    cleanup_tasks.push(annotated_router);

    let mut annotated_rxs_iter = annotated_rxs.into_iter();
    let initial_annotated_file_stream = annotated_rxs_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_annotated_taxon_stream = annotated_rxs_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (unidentified_rxs, unidentified_router) = fanout_to_channels(
        ReceiverStream::new(unidentified_rx),
        2,
        "initial_unidentified",
        &config,
        StreamDataType::JustBytes
    ).await?;

    cleanup_tasks.push(unidentified_router);

    let mut unidentified_rxs_iter = unidentified_rxs.into_iter();
    let initial_unidentified_file_stream = unidentified_rxs_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_unidentified_taxon_stream = unidentified_rxs_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Add monitors so you can actually see progress (no more 30-minute silences)
    let initial_annotated_file_stream = monitor_stream(
        ReceiverStream::new(initial_annotated_file_stream),
        "initial_annotated_file_stream_to_write",
        Duration::from_secs(15),
    );

    let initial_unidentified_file_stream = monitor_stream(
        ReceiverStream::new(initial_unidentified_file_stream),
        "initial_unidentified_file_stream_to_write",
        Duration::from_secs(15),
    );

    // Write tasks are now eager consumers → EOF propagates instantly
    cleanup_tasks.push(write_fasta_stream_to_file(
        initial_annotated_file_stream,
        annotated_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_annotated_file_stream_to_write",
    ));

    cleanup_tasks.push(write_fasta_stream_to_file(
        initial_unidentified_file_stream,
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


    let (assembly_outputs, assembly_bam_out_stream,
        post_assembly_cleanup_tasks,
        post_assembly_cleanup_receivers,
        post_assembly_temp_files) = process_assembly(
        config.clone(),
        &assembly_out_dir,
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        duplicate_clusters.clone(),
        paired,
        4 // this can become as CLI arg if needed
    ).await?;

    cleanup_tasks.extend(post_assembly_cleanup_tasks);
    cleanup_receivers.extend(post_assembly_cleanup_receivers);
    temp_files.extend(post_assembly_temp_files);


    let generate_assembly_coverage_result = generate_assembly_coverage(
        config.clone(),
        assembly_bam_out_stream,
        &assembly_outputs.contigs_fasta,
        &assembly_outputs.coverage_json,
        &assembly_outputs.coverage_summary_csv,
    )
        .await
        .map_err(|e| anyhow!("generate_assembly_coverage failed: {}", e))?;
    eprintln!("coverage result {:?}", generate_assembly_coverage_result);


    let nt_file = config
        .args
        .nt
        .clone()
        .ok_or(PipelineError::MissingArgument(
            "NT file required for fasta extraction".into(),
        ))?;

    let nt_offset_db_file = config
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

    let nr_offset_db_file = config
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
            let (
                nt_read_dict_map,
                nt_accession_dict_noarc,
                nt_selected_genera,
                _nt_total_reads,
            ) = nt_hit_summary_handle
                .await
                .map_err(|e| PipelineError::Other(anyhow!("NT hit summary task panicked: {}", e)))?
                .map_err(|e| PipelineError::Other(anyhow!("NT hit summary parsing failed: {}", e)))?;

            let (nt_ref_fasta_path, nt_ref_fasta_temp_dir) = build_reference_fasta_from_selected_genera(
                config.clone(),
                &nt_selected_genera,
                NT_TAG,
                &PathBuf::from(nt_file),
                &PathBuf::from(nt_offset_db_file),
            )
                .await
                .map_err(|e| PipelineError::Other(e.into()))?;
            Ok::<(AHashMap<String, Arc<ReadHit>>, Arc<AHashMap<String, AccessionHit>>, PathBuf, TempDir), PipelineError>((
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
            let (
                nr_read_dict_map,
                nr_accession_dict_noarc,
                nr_selected_genera,
                _nr_total_reads,
            ) = nr_hit_summary_handle
                .await
                .map_err(|e| PipelineError::Other(anyhow!("NR hit summary task panicked: {}", e)))?
                .map_err(|e| PipelineError::Other(anyhow!("NR hit summary parsing failed: {}", e)))?;

            let (nr_ref_fasta_path, nr_ref_fasta_temp_dir) = build_reference_fasta_from_selected_genera(
                config.clone(),
                &nr_selected_genera,
                NR_TAG,
                &PathBuf::from(nr_file),
                &PathBuf::from(nr_offset_db_file),
            )
                .await
                .map_err(|e| PipelineError::Other(e.into()))?;
            Ok::<(AHashMap<String, Arc<ReadHit>>, Arc<AHashMap<String, AccessionHit>>, PathBuf, TempDir), PipelineError>((
                nr_read_dict_map,
                Arc::new(nr_accession_dict_noarc),
                nr_ref_fasta_path,
                nr_ref_fasta_temp_dir,
            ))
        }
    });

    let (nt_prep_res, nr_prep_res) = tokio::try_join!(nt_prep_handle, nr_prep_handle)
        .map_err(|e| PipelineError::Other(anyhow!("Reference prep task panicked: {e}")))?;

    let (nt_read_dict_map, nt_accession_dict, nt_ref_fasta_path, nt_ref_fasta_temp_dir) = nt_prep_res?;
    let (nr_read_dict_map, nr_accession_dict, nr_ref_fasta_path, nr_ref_fasta_temp_dir) = nr_prep_res?;

    let nt_read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>> = Arc::new(Mutex::new(nt_read_dict_map));
    let nr_read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>> = Arc::new(Mutex::new(nr_read_dict_map));

    info!("NT ref path: {}", nt_ref_fasta_path.display());
    info!("NR ref path: {}", nr_ref_fasta_path.display());

    final_temp_dirs.push(nt_ref_fasta_temp_dir);
    final_temp_dirs.push(nr_ref_fasta_temp_dir);


    let (nt_counts, nr_counts) = tokio::try_join!(
        async { nt_counts_task.await.map_err(|e| PipelineError::Other(anyhow!(e.to_string()))) },
        async { nr_counts_task.await.map_err(|e| PipelineError::Other(anyhow!(e.to_string()))) },
    )?;
    let nt_counts = nt_counts?;
    let nr_counts = nr_counts?;


    let combined_path = out_dir.join(rename_file_path(
        &sample_base_buf,
        None,
        Some("taxon_counts_with_dcr.json"),
        "_",
    ));

    let nt_counts_for_combine = nt_counts.clone();
    let nr_counts_for_combine = nr_counts.clone();

    let combine_handle = tokio::spawn(async move {
        let (_combined_path, write_json_task) = combine_taxon_counts(
            &nt_counts_for_combine,
            &nr_counts_for_combine,
            combined_path,
        )
            .await
            .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts failed: {}", e)))?;

        write_json_task
            .await
            .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts write task failed: {}", e)))?
    });

    cleanup_tasks.push(combine_handle);


    let nt_blast_concurrency = compute_phase_concurrency(
        &config,
        "nt_blast_contigs",
        2.0,           // ~2 GB per thread — BLAST can be very memory-hungry (index + large queries)
        5.0,           // strong CPU-bound phase (alignment + scoring)
        128,           // high cap — BLAST scales decently to 128 on EPYC
        8,             // min — still want parallelism on MacBook Air
    );
    info!("NT blast_contigs concurrency: {}", nt_blast_concurrency);

    let nt_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let read2contig = assembly_outputs.read2contig.clone();

        async move {
            blast_contigs(
                config,
                NT_TAG,
                nt_blast_stream.into(),
                nt_blast_hit_stream,
                nt_read_dict.clone(),
                nt_accession_dict,
                nt_counts,
                &contigs_fasta_path.clone(),
                read2contig,
                &nt_ref_fasta_path,
                duplicate_clusters,
                lineage_map,
                should_keep_filter,
                4,
                nt_blast_concurrency
            ).await
        }
    });

    let nr_blast_concurrency = compute_phase_concurrency(
        &config,
        "nr_blast_contigs",
        2.0,           // ~2 GB per thread — BLAST can be very memory-hungry (index + large queries)
        5.0,           // strong CPU-bound phase (alignment + scoring)
        128,           // high cap — BLAST scales decently to 128 on EPYC
        8,             // min — still want parallelism on MacBook Air
    );
    info!("NR blast_contigs concurrency: {}", nr_blast_concurrency);

    let nr_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let read2contig = assembly_outputs.read2contig.clone();

        async move {
            blast_contigs(
                config,
                NR_TAG,
                nr_blast_stream,
                nr_blast_hit_stream,
                nr_read_dict.clone(),
                nr_accession_dict,
                nr_counts,
                &contigs_fasta_path.clone(),
                read2contig,
                &nr_ref_fasta_path,
                duplicate_clusters,
                lineage_map,
                should_keep_filter,
                4,
                nr_blast_concurrency
            ).await
        }
    });

    let nt_post_handle = tokio::spawn({
        let config = config.clone();
        async move {
            let nt_res = nt_handle
                .await
                .map_err(|e| PipelineError::Other(anyhow!("NT blast_contigs task panicked: {e}")))??;

            let (
                nt_read_dict,
                nt_refined_counts,
                nt_contig_summary,
                nt_refined_m8_stream_out,
                nt_refined_hit_summary_stream_out,
                nt_refined_m8_top_stream_out,
                nt_cleanup_tasks,
                nt_cleanup_receivers,
                nt_temp_files,
            ) = nt_res;

            let mut cleanup_tasks = nt_cleanup_tasks;
            let cleanup_receivers = nt_cleanup_receivers;
            let temp_files = nt_temp_files;

            let (nt_m8_rxs, nt_m8_router) = fanout_to_channels(
                ReceiverStream::new(nt_refined_m8_stream_out),
                3,
                "nt_m8",
                &config,
                StreamDataType::JustBytes
            ).await?;
            cleanup_tasks.push(nt_m8_router);

            let mut nt_m8_iter = nt_m8_rxs.into_iter();
            let nt_m8_merge = nt_m8_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nt_m8_map = nt_m8_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nt_m8_viz = nt_m8_iter.next().ok_or(PipelineError::EmptyStream)?;

            let (nt_hitsummary_rxs, nt_hitsummary_router) = fanout_to_channels(
                ReceiverStream::new(nt_refined_hit_summary_stream_out),
                3,
                "nt_hitsummary",
                &config,
                StreamDataType::JustBytes
            ).await?;
            cleanup_tasks.push(nt_hitsummary_router);

            let mut nt_hitsummary_iter = nt_hitsummary_rxs.into_iter();
            let nt_hit_summary_merge = nt_hitsummary_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nt_hit_summary_taxid = nt_hitsummary_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nt_hit_summary_coverage = nt_hitsummary_iter.next().ok_or(PipelineError::EmptyStream)?;

            Ok::<_, PipelineError>((
                nt_read_dict,
                nt_refined_counts,
                nt_contig_summary,
                nt_m8_merge,
                nt_m8_map,
                nt_m8_viz,
                nt_hit_summary_merge,
                nt_hit_summary_taxid,
                nt_hit_summary_coverage,
                nt_refined_m8_top_stream_out,
                cleanup_tasks,
                cleanup_receivers,
                temp_files,
            ))
        }
    });

    let nr_post_handle = tokio::spawn({
        let config = config.clone();
        async move {
            let nr_res = nr_handle
                .await
                .map_err(|e| PipelineError::Other(anyhow!("NR blast_contigs task panicked: {e}")))??;

            let (
                nr_read_dict,
                nr_refined_counts,
                nr_contig_summary,
                nr_refined_m8_stream_out,
                nr_refined_hit_summary_stream_out,
                _nr_refined_m8_top_stream_out,
                nr_cleanup_tasks,
                nr_cleanup_receivers,
                nr_temp_files,
            ) = nr_res;

            let mut cleanup_tasks = nr_cleanup_tasks;
            let cleanup_receivers = nr_cleanup_receivers;
            let temp_files = nr_temp_files;

            let (nr_m8_rxs, nr_m8_router) = fanout_to_channels(
                ReceiverStream::new(nr_refined_m8_stream_out),
                2,
                "nr_m8",
                &config,
                StreamDataType::JustBytes
            ).await?;
            cleanup_tasks.push(nr_m8_router);

            let mut nr_m8_iter = nr_m8_rxs.into_iter();
            let nr_m8_merge = nr_m8_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nr_m8_map = nr_m8_iter.next().ok_or(PipelineError::EmptyStream)?;

            let (nr_hitsummary_rxs, nr_hitsummary_router) = fanout_to_channels(
                ReceiverStream::new(nr_refined_hit_summary_stream_out),
                3,
                "nr_hitsummary",
                &config,
                StreamDataType::JustBytes
            ).await?;
            cleanup_tasks.push(nr_hitsummary_router);

            let mut nr_hitsummary_iter = nr_hitsummary_rxs.into_iter();
            let nr_hit_summary_merge = nr_hitsummary_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nr_hit_summary_preload = nr_hitsummary_iter.next().ok_or(PipelineError::EmptyStream)?;
            let nr_hit_summary_taxid = nr_hitsummary_iter.next().ok_or(PipelineError::EmptyStream)?;

            Ok::<_, PipelineError>((
                nr_read_dict,
                nr_refined_counts,
                nr_contig_summary,
                nr_m8_merge,
                nr_m8_map,
                nr_hit_summary_merge,
                nr_hit_summary_preload,
                nr_hit_summary_taxid,
                cleanup_tasks,
                cleanup_receivers,
                temp_files,
            ))
        }
    });

    let (nt_post_res, nr_post_res) = tokio::join!(nt_post_handle, nr_post_handle);

    let nt_post_res = nt_post_res
        .map_err(|e| PipelineError::Other(anyhow!("NT blast_contigs postprocess task panicked: {e}")))??;
    let nr_post_res = nr_post_res
        .map_err(|e| PipelineError::Other(anyhow!("NR blast_contigs postprocess task panicked: {e}")))??;

    let (
        _nt_read_dict,
        nt_refined_counts,
        nt_contig_summary,
        nt_m8_merge,
        _nt_m8_map,
        nt_m8_viz,
        nt_hit_summary_merge,
        nt_hit_summary_taxid,
        nt_hit_summary_coverage,
        nt_refined_m8_top_stream_out,
        nt_cleanup_tasks,
        nt_cleanup_receivers,
        nt_temp_files,
    ) = nt_post_res;
    cleanup_tasks.extend(nt_cleanup_tasks);
    cleanup_receivers.extend(nt_cleanup_receivers);
    temp_files.extend(nt_temp_files);

    let (
        _nr_read_dict,
        nr_refined_counts,
        nr_contig_summary,
        nr_m8_merge,
        _nr_m8_map,
        nr_hit_summary_merge,
        nr_hit_summary_preload,
        nr_hit_summary_taxid,
        nr_cleanup_tasks,
        nr_cleanup_receivers,
        nr_temp_files,
    ) = nr_post_res;
    cleanup_tasks.extend(nr_cleanup_tasks);
    cleanup_receivers.extend(nr_cleanup_receivers);
    temp_files.extend(nr_temp_files);

    // let spades_task = spades_assembly(
    //     config.clone(),
    //     non_host_r1_path.clone(),
    //     non_host_r2_path_opt.clone(),
    //     &out_dir,
    // ).await?;
    //
    // let spades_completion = spades_task.await;
    //
    // match spades_completion {
    //     Ok(Ok(())) => {
    //         info!("SPAdes assembly task completed successfully");
    //     }
    //     Ok(Err(e)) => {
    //         warn!("SPAdes assembly failed: {}", e);
    //     }
    //     Err(join_err) => {
    //         error!("SPAdes task panicked or was cancelled: {}", join_err);
    //     }
    // }

    // ────────────────────────────────────────────────────────────────
    // PRELOAD NR alignments (parallel version)
    // ────────────────────────────────────────────────────────────────
    let nr_alignment_per_read: Arc<DashMap<String, SpeciesAlignmentResults, AHashRandomState>> =
        Arc::new(DashMap::with_capacity_and_hasher(
            80_000_000,
            AHashRandomState::new()
        ));

    let preload_handle = tokio::spawn(preload_nr_alignments_parallel(
        config.clone(),
        nr_hit_summary_preload,
        Arc::clone(&nr_alignment_per_read),
    ));

    // ────────────────────────────────────────────────────────────────
    // Spawn merged taxon counts (unchanged except we await preload first)
    // ────────────────────────────────────────────────────────────────
    let merged_cleanup_tasks = compute_merged_taxon_counts(
        config.clone(),
        nt_m8_merge,
        nt_hit_summary_merge,
        nt_contig_summary,
        nr_m8_merge,
        nr_hit_summary_merge,
        nr_contig_summary,
        lineage_map.clone(),
        should_keep_filter.clone(),
        duplicate_clusters.clone(),
        out_dir.join("refined.m8"),
        out_dir.join("refined.hitsummary.tab"),
        out_dir.join("refined_taxon_counts_with_dcr.json"),
        out_dir.join("assembly_combined_contig_summary.json"),
        nr_alignment_per_read,
    ).await?;
    cleanup_tasks.extend(merged_cleanup_tasks);

    // Wait for preload before proceeding (keeps exact ordering semantics from Python)
    let _ = preload_handle.await
        .map_err(|e| PipelineError::Other(anyhow!("NR preload task panicked: {}", e)))??;


    let refined_combined_path = out_dir.join(rename_file_path(&sample_base_buf, None, Some("refined_taxon_counts_with_dcr.json"), "_"));
    let (_refined_combined_path, refined_write_json_task) = combine_taxon_counts(
        &nt_refined_counts,
        &nr_refined_counts,
        refined_combined_path,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts failed: {}", e)))?;
    cleanup_tasks.push(refined_write_json_task);


    // ────────────────────────────────────────────────────────────────
    // NEW: Build tiny contig-level maps for the REFINED path
    // ────────────────────────────────────────────────────────────────
    let (nt_refined_map, nr_refined_map) = tokio::try_join!(
        collect_hit_summary_to_accession_map_concurrent(
            config.clone(),
            nt_hit_summary_for_refined
        ),
        collect_hit_summary_to_accession_map_concurrent(
            config.clone(),
            nr_hit_summary_for_refined
        )
    )?;

    // Contig FASTA stream (from assembly)
    // ───────────────────────────────────────────────────────────────
    let contigs_fasta = assembly_outputs.contigs_fasta.clone();
    let contigs_file = tokio::fs::File::open(&contigs_fasta).await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to open contigs.fasta: {}", e)))?;

    let contigs_rx = parse_bytes::<TokioFile>(contigs_file, &config, StreamDataType::JustBytes).await
        .map_err(|e| PipelineError::Other(anyhow!("parse_bytes failed: {}", e)))?;

    let contigs_stream = ReceiverStream::new(contigs_rx);

    // Empty cluster stream (contigs have no duplicates)
    let (_cluster_tx, _cluster_rx) = mpsc::channel::<ParseOutput>(1);
    drop(_cluster_tx);

    let assembly_dir = out_dir.join("assembly");
    tokio::fs::create_dir_all(&assembly_dir).await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create assembly dir: {}", e)))?;


    let nr_annot_concurrency = compute_phase_concurrency(
        &config,
        "nr_annot_concurrency",
        0.4,
        4.0,
        128,
        8,
    );

    let (
        mapped_contigs_rx,
        unidentified_contigs_rx,
        _unique_unidentified_rx,
        mut annot_tasks,
        mut annot_rxs,
    ) = generate_annotated_fasta_stream(
        config.clone(),
        contigs_stream,
        duplicate_clusters.clone(),
        nt_refined_map,                     // tiny contig → NT accession
        nr_refined_map,                     // tiny contig → NR accession
        nr_annot_concurrency,
    ).await
        .map_err(|e| PipelineError::Other(anyhow!("Refined generate_annotated_fasta failed: {}", e)))?;
    cleanup_tasks.append(&mut annot_tasks);
    cleanup_receivers.append(&mut annot_rxs);


    let (
        taxid_mapped_rx,
        taxid_combined_rx,
        load_nt_task,
        load_nr_task,
        taxid_main_task,
    ) = generate_taxid_fasta(
        config.clone(),
        ReceiverStream::new(mapped_contigs_rx),
        ReceiverStream::new(unidentified_contigs_rx),
        ReceiverStream::new(nt_hit_summary_taxid),
        ReceiverStream::new(nr_hit_summary_taxid),
        lineage_map.clone(),
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("generate_taxid_fasta failed: {}", e)))?;


    cleanup_tasks.push(load_nt_task);
    cleanup_tasks.push(load_nr_task);
    cleanup_tasks.push(taxid_main_task);


    let (taxid_mapped_rxs, taxid_mapped_router) = fanout_to_channels(
        ReceiverStream::new(taxid_mapped_rx),
        2,
        "taxid_mapped",
        &config,
        StreamDataType::JustBytes
    ).await?;
    cleanup_tasks.push(taxid_mapped_router);

    let mut taxid_mapped_iter = taxid_mapped_rxs.into_iter();
    let taxid_mapped_file   = taxid_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;
    let taxid_mapped_locator = taxid_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;


    let mapped_path = assembly_dir.join("refined_taxid_annot_mapped_only_fasta");
    let combined_path = assembly_dir.join("refined_taxid_annot_fasta");


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
    ) = generate_taxid_fasta(
        config.clone(),
        ReceiverStream::new(initial_annotated_taxon_stream),
        ReceiverStream::new(initial_unidentified_taxon_stream),
        nt_hit_summary_for_taxid,
        nr_hit_summary_for_taxid,
        lineage_map.clone(),
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Experimental generate_taxid_fasta failed: {}", e)))?;

    cleanup_tasks.push(initial_load_nt_task);
    cleanup_tasks.push(initial_load_nr_task);
    cleanup_tasks.push(initial_taxid_main_task);

    let initial_mapped_path = out_dir.join("taxid_annot_mapped_only.fasta");
    let initial_combined_path = out_dir.join("taxid_annot.fasta");


    info!("[run] starting fanout_to_channels for initial taxid FASTA (mapped + combined)");

    let (initial_mapped_rxs, initial_mapped_router) = fanout_to_channels(
        ReceiverStream::new(initial_taxid_mapped_rx),
        3,
        "initial_taxid_mapped",
        &config,
        StreamDataType::JustBytes
    ).await?;
    cleanup_tasks.push(initial_mapped_router);

    let mut initial_mapped_iter = initial_mapped_rxs.into_iter();
    let initial_taxid_mapped_file = initial_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_mapped_locator = initial_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_mapped_count = initial_mapped_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (initial_combined_rxs, initial_combined_router) = fanout_to_channels(
        ReceiverStream::new(initial_taxid_combined_rx),
        2,
        "initial_taxid_combined",
        &config,
        StreamDataType::JustBytes
    ).await?;
    cleanup_tasks.push(initial_combined_router);

    let mut initial_combined_iter = initial_combined_rxs.into_iter();
    let initial_taxid_combined_file = initial_combined_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_combined_count = initial_combined_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Add monitors so we finally see progress logs again
    let initial_taxid_mapped_file = monitor_stream(
        ReceiverStream::new(initial_taxid_mapped_file),
        "initial_taxid_mapped_file_to_write",
        Duration::from_secs(15),
    );
    let initial_taxid_combined_file = monitor_stream(
        ReceiverStream::new(initial_taxid_combined_file),
        "initial_taxid_combined_file_to_write",
        Duration::from_secs(15),
    );



    let initial_write_mapped_handle = write_fasta_stream_to_file(
        initial_taxid_mapped_file, // Clone rx if needed for logging
        initial_mapped_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_taxid_mapped_file_stream_to_write",
    );
    cleanup_tasks.push(initial_write_mapped_handle);

    let initial_write_combined_handle = write_fasta_stream_to_file(
        initial_taxid_combined_file,
        initial_combined_path.clone(),
        config.clone(),
        StreamDataType::JustBytes,
        "initial_taxid_combined_file_stream_to_write",
    );
    cleanup_tasks.push(initial_write_combined_handle);


    let initial_mapped_count = stream_record_counter(initial_taxid_mapped_count, false).await?;
    let initial_combined_count = stream_record_counter(initial_taxid_combined_count, false).await?;
    info!("Experimental: {} mapped, {} combined records", initial_mapped_count, initial_combined_count);

    info!("Starting generate_taxid_locator on second");
    let (initial_locator_outputs, mut initial_locator_tasks, _initial_locator_receivers) = generate_taxid_locator(
        config.clone(),
        initial_taxid_mapped_locator,
        out_dir.clone(),
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Experimental generate_taxid_locator failed: {}", e)))?;
    info!("second generate_taxid_locator odone ");
    cleanup_tasks.append(&mut initial_locator_tasks);
    info!("Experimental taxid locator files generated: {:?}", initial_locator_outputs);

    let nt_info_db_path = config.args.nt_info_tab
        .clone()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("nt_info.tab"));


    let skip_coverage_viz = {
        let mut skip = false;

        if !assembly_outputs.coverage_json.exists() {
            warn!("Skipping coverage viz: contig_coverage_json missing (likely due to SPAdes assembly failure)");
            skip = true;
        } else if assembly_outputs.coverage_json.metadata().map(|m| m.len() == 0).unwrap_or(true) {
            warn!("Skipping coverage viz: contig_coverage_json is empty");
            skip = true;
        }

        if !assembly_outputs.contig_stats_json.exists() {
            warn!("Skipping coverage viz: contig_stats_json missing");
            skip = true;
        } else if assembly_outputs.contig_stats_json.metadata().map(|m| m.len() == 0).unwrap_or(true) {
            warn!("Skipping coverage viz: contig_stats_json is empty");
            skip = true;
        }

        skip
    };

    let out_dir_for_viz = out_dir.clone();
    let coverage_task = if !skip_coverage_viz {
        info!("All coverage viz inputs present — spawning visualization task");
        Some(tokio::spawn(async move {
            generate_coverage_viz(
                ReceiverStream::new(nt_hit_summary_coverage),
                ReceiverStream::new(nt_refined_m8_top_stream_out),
                assembly_outputs.coverage_json,
                assembly_outputs.contig_stats_json,
                ReceiverStream::new(nt_m8_viz),
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
        warn!("Skipping coverage viz: required inputs missing/empty (likely assembly failure)");
        None
    };





    // *******************
    // Results retrieval
    // *******************

    // Clone paths needed in multiple concurrent async blocks.
    let out_dir_host = out_dir.clone();
    let out_dir_kallisto = out_dir.clone();

    let validation_and_qc_fut = async move {
        let raw_count = join_with_error_handling(raw_count_task).await?;
        info!("Processed {} raw reads (additive from R1 and R2 if paired)", raw_count);

        let stats = join_with_error_handling(val_count_task).await?;
        info!(
        "Processed {} validated, {} undersized, {} oversized reads",
        stats.validated, stats.undersized, stats.oversized
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

        Ok::<_, PipelineError>((raw_count, stats, ercc_mapped_count))
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

        let host_hisat2_counts = host_hisat2_count_rx
            .await
            .map_err(|e| PipelineError::Other(anyhow!("Host hisat2 counts receiver failed: {}", e)))?;
        info!("Host hisat2: mapped counts: {:?}", host_hisat2_counts);

        if paired {
            info!(
            "Computing host insert size statistics from {}",
            host_bt2_bam_path.display()
        );

            let stats = compute_insert_size_stats_from_bam(
                host_bt2_bam_path.clone(),
                None,
            )
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
            if let Err(e) = plot_insert_sizes(&stats.insert_sizes, sample_base.as_str(), &plot_path) {
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
                warn!("Kallisto failed (non-fatal, possibly no ERCC spiked in): {}", e);
            }
        }

        let kallisto_results_task = kallisto_results(out_dir_kallisto.join("kallisto"), kallisto_ercc_tx)
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

        ercc_writer
            .write_all(b"target_id\test_counts\n")
            .await?;

        for (target_id, est_counts) in &kallisto_ercc_counts.ercc_counts {
            let line = format!("{}\t{:.6}\n", target_id, est_counts);
            ercc_writer.write_all(line.as_bytes()).await?;
        }
        ercc_writer.flush().await?;
        info!("Wrote kallisto_ERCC_counts_tsv");

        Ok::<_, PipelineError>((kallisto_ercc_counts, ()))
    };

    let (_validation_and_qc, _host_bam_and_counts, _kallisto_and_taxon) =
        tokio::try_join!(
        validation_and_qc_fut,
        host_bam_and_counts_fut,
        kallisto_and_taxon_fut
    )?;

    // *******************
    // Cleanup
    // *******************

    // Helper to handle downcast + NotFound tolerance in one place
    fn handle_cleanup_error(e: anyhow::Error, context: &str) -> Result<(), PipelineError> {
        match e.downcast::<std::io::Error>() {
            Ok(io_err) => {
                if io_err.kind() == std::io::ErrorKind::NotFound {
                    warn!("Cleanup ignored missing file/dir (non-fatal) [{}]: {}", context, io_err);
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
            Ok(Ok(_)) => {}  // success

            Ok(Err(e)) => {
                handle_cleanup_error(e, "receiver")?;
            }

            Err(oneshot_err) => {
                // Channel broken (panic, drop, empty stream, etc.)
                // Usually safe to ignore if stream was empty
                warn!("Cleanup receiver channel failed (possibly empty): {}", oneshot_err);
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
