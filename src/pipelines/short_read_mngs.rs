use std::cmp::Reverse;
use std::cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd};
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader, ErrorKind};
use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::collections::hash_map::DefaultHasher;
use std::process::{Command as StdCommand, Stdio};

use ahash::AHashMap;
use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use dashmap::DashMap;
use rand::seq::SliceRandom;
use fst::Map;
use futures::future::try_join_all;
use log::{self, debug, error, info, warn, LevelFilter};
use needletail::parse_fastx_file;
use noodles::bam::r#async::io::Reader as BamAsyncReader;
use noodles::bam::record::Record;
use noodles::sam::alignment::record::cigar::{op::Kind as OpKind, Op};
use once_cell::sync::Lazy;
use rand::prelude::*;
use rand_core::{OsRng, RngCore};
use rayon::prelude::*;
use regex::Regex;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use tempfile::NamedTempFile;
use tempfile::TempDir;
use tokio::fs;
use tokio::fs::{File as TokioFile, OpenOptions as TokioOpenOptions};
use tokio::io::{AsyncBufReadExt, AsyncRead, AsyncReadExt, AsyncWriteExt, BufReader as TokioBufReader, BufWriter};
use tokio::process::Command;
use tokio::sync::mpsc;
use tokio::sync::mpsc::{channel, Sender};
use tokio::sync::oneshot;
use tokio::sync::Semaphore;
use tokio::task::JoinHandle;
use tokio::process::Child;
use tokio::time::{sleep, Duration, Instant};
use tokio::try_join;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::{Stream, StreamExt};
use tokio_util::io::StreamReader;
use tokio::runtime::Builder as RuntimeBuilder;
use twox_hash::XxHash64;
use rayon::prelude::*;
use crate::cli::Technology;
use crate::config::defs::{DiamondSubcommand, KallistoSubcommand, Lineage, PipelineError, ReadCountingMode,
                          ReadStats, RunConfig, SamtoolsStats, SamtoolsSubcommand, StreamDataType,
                          Taxid, BCFTOOLS_TAG, BLASTN_TAG, BLASTX_TAG, BOWTIE2_TAG,
                          CZID_DEDUP_TAG, DIAMOND_TAG, FASTP_TAG, HISAT2_TAG, KALLISTO_TAG,
                          KRAKEN2_TAG, LOG_NORMAL_POSITIVE_DOUBLE, MAFFT_TAG, MAKEBLASTDB_TAG,
                          MINIMAP2_TAG, MIN_NORMAL_POSITIVE_DOUBLE, NR_TAG,
                          NT_TAG, PIGZ_TAG, QUAST_TAG, READ_COUNTING_MODE,
                          SAMTOOLS_TAG, SEQKIT_TAG, SPADES_TAG, STAR_TAG, ClusterInfo, DuplicateClusters};
use crate::utils::blast::{consensus_level, generate_taxon_count_json_from_m8, AggBucket, M8Record,
                          TaxonCount, ContigSummaryEntry, compute_merged_taxon_counts, SpeciesAlignmentResults};
use crate::utils::command::blastn::{BlastnArgGenerator, BlastnConfig};
use crate::utils::command::blastx::{BlastxArgGenerator, BlastxConfig};
use crate::utils::command::bowtie2::{bowtie2_index_prep, Bowtie2Config};
use crate::utils::command::czid_dedup::CzidDedupConfig;
use crate::utils::command::diamond::{diamond_index_prep, DiamondArgGenerator, DiamondConfig, compute_optimal_block_size};
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::hisat2::{hisat2_index_prep, Hisat2Config};
use crate::utils::command::kallisto::KallistoConfig;
use crate::utils::command::makeblastdb::{MakeblastdbArgGenerator, MakeblastdbConfig};
use crate::utils::command::minimap2::{minimap2_index_prep, Minimap2ArgGenerator, Minimap2Config};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::coverage_viz::generate_coverage_viz;
use crate::utils::fastx::{compare_read_ids, parse_header, raw_read_count, read_fasta,
                          read_fastq, stream_record_counter, write_fasta_stream_to_file,
                          SequenceRecord, generate_taxid_fasta, generate_taxid_locator,
                          filter_fastq_to_bytes_stream, parse_byte_stream_to_fastq};
use crate::utils::file::{available_space_for_path, choose_temp_dir, file_path_manipulator,
                         file_size, rename_file_path, resolve_optional_path,
                         validate_file_inputs, write_byte_stream_to_file, write_parse_output_to_file,
                         write_vecu8_to_file};
use crate::utils::paf::PafRecord;
use crate::utils::plotting::plot_insert_sizes;
use crate::utils::sambam::{generate_info_from_bam_stream, compute_insert_size_stats_from_bam};
use crate::utils::stats::parse_samtools_stats;
use crate::utils::streams::{create_fifo, deinterleave_fastq_stream, interleave_fastq_streams, join_with_error_handling,
                            parse_child_output, parse_fasta, parse_fastq, read_child_output_to_vec, spawn_cmd,
                            stream_to_cmd, stream_to_file, t_junction, write_to_fifo,
                            ChannelReader, ChildStream, ParseMode,
                            ParseOutput, ToBytes};
use crate::utils::streams::deinterleave_fastq_stream_to_fifos;
use crate::utils::system::{detect_ram, compute_phase_concurrency};
use crate::utils::taxonomy::{build_should_keep_filter, get_top_m8_nr,
                             get_top_m8_nt, load_taxid_lineages_db, validate_taxid_lineage};

const UNMAPPED_HEADER_PREFIX: &str = ">NR::NT::";
const MAX_ACCESSION_SEQUENCE_LEN: u64 = 100_000_000;
const EST_BYTES_PER_ACCESSION: u64 = 20_000; // ~10k seq + header
const MAX_PARSE_ERRORS: usize = 100; // Threshold before failing
const MAX_SPADES_WORK_DIR: u64 = 500_000_000;

const MAX_NUM_BINS_COVERAGE: usize = 500;
const NUM_ACCESSIONS_PER_TAXON: usize = 10;
const MIN_CONTIG_SIZE: u64 = 500;

const MIN_REF_FASTA_SIZE: u64 = 25;
const MIN_ASSEMBLED_CONTIG_SIZE: u64 = 25;
const MIN_ASSEMBLY_READ_COUNT: usize = 56;

const DIAMOND_TEMP:u64 = 20 * 1024 * 1024 * 1024; // 20 GB

const MM2_TARGET_THREADS_PER : usize = 8;

const EST_GB_PER_MM2_JOB_STD: f64 = 80.0;   // conservative for 1 TB, gives ~6 jobs
const EST_GB_PER_MM2_JOB_EXTRA: f64 = 100.0;  // for 1.5 TB, gives ~10–12 jobs
const MM2_JOB_RAM_THRESHOLD: usize = 1400;


static FIX_COMMA_REGEXP: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r",(?=[^\s])").unwrap()
});

#[derive(Debug)]
pub struct KallistoResults {
    pub ercc_counts: Vec<(String, f64)>, // (target_id, est_counts)
    pub transcript_to_gene: Vec<(String, String)>, // (transcript_id, gene_id)
}

// SampleItem for reservoir sampling
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
    pub contigs_ram_fasta: PathBuf,
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
    let mut cleanup_receivers = Vec::new();


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
        Some(config.args.technology.clone()),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        config.base_buffer_size,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);

    // Split the byte stream for fastp and write
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(val_done_rx);

    let mut val_streams_iter = val_streams.into_iter();
    let val_ercc_bowtie2_filter_out_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let val_file_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let val_file_write_task = write_byte_stream_to_file(
        &validated_interleaved_file_path,
        ReceiverStream::new(val_file_stream),
        Some(config.base_buffer_size),
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(val_file_write_task);

    Ok((ReceiverStream::new(val_ercc_bowtie2_filter_out_stream),  cleanup_tasks, cleanup_receivers, raw_count_task, val_count_task))
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

async fn bowtie2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    bowtie2_options: HashMap<String, Option<String>>,
    output_bam_path: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,                    // unmapped FASTQ stream
    oneshot::Receiver<u64>,                         // mapped read count
    Vec<JoinHandle<Result<(), anyhow::Error>>>,     // cleanup tasks
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, // cleanup receivers
    JoinHandle<Result<()>>,                         // BAM write task handle
    PathBuf,                                        // written BAM path
    Vec<TempDir>
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    ).await?;

    // BT2
    let bt2_config_view = Bowtie2Config {
        bt2_index_path: bt2_index_path.clone(),
        paired,
        r1_path: None,
        r2_path: None,
        option_fields: bowtie2_options,
    };

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bt2_child, bt2_stream_task, bt2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        BOWTIE2_TAG,
        bt2_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bt2_stream_task);
    cleanup_tasks.push(bt2_err_task);

    let bt2_out_stream = {
        let mut guard = bt2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BOWTIE2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),          // name-sorted
            ("-u".to_string(), None),          // uncompressed
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

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
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
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // ────────────────────────────────────────────────────────────────
    // T-junction: split for unmapped FASTQ, mapped count, and BAM write
    // ────────────────────────────────────────────────────────────────

    let num_tees = 3; // unmapped → fastq, count, bam file write
    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = t_junction(
        bam_rx_stream,
        num_tees,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "bowtie2_bam_split".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let bam_file_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let config_for_write = config.clone();
    let output_bam_path_for_write = output_bam_path.clone();

    // ────────────────────────────────────────────────────────────────
    // Write full BAM to disk
    // ────────────────────────────────────────────────────────────────

    let bam_write_task: JoinHandle<Result<()>> = tokio::spawn(async move {
        write_byte_stream_to_file(
            &output_bam_path_for_write,
            ReceiverStream::new(bam_file_stream),
            Some(config_for_write.base_buffer_size),
        )
            .await
            .map_err(|e| anyhow!("Failed to write BAM file {}: {}", output_bam_path_for_write.display(), e))?;

        Ok(())
    });

    // ────────────────────────────────────────────────────────────────
    // Count mapped reads (rough approximation via line count)
    // ────────────────────────────────────────────────────────────────

    let (count_tx, count_rx) = oneshot::channel();
    let count_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(count_stream);
        let mut line_count = 0u64;
        while let Some(item) = stream.next().await {
            if let ParseOutput::Bytes(b) = item {
                if !b.is_empty() {
                    line_count += 1;
                }
            }
        }
        let count = line_count / if paired { 8 } else { 4 }; // rough pair estimate
        let _ = count_tx.send(count);
        Ok(())
    });
    cleanup_tasks.push(count_task);

    // ────────────────────────────────────────────────────────────────
    // Samtools fastq: extract unmapped reads
    // ────────────────────────────────────────────────────────────────

    let unmapped_flag = if paired { "-f13" } else { "-f4" };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag.to_string(), None),
            ("-".to_string(), None),
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
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
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let mut temp_dirs = vec![temp_dir];

    Ok((
        ReceiverStream::new(unmapped_fastq_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        bam_write_task,
        output_bam_path,
        temp_dirs
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
    input_stream: ReceiverStream<ParseOutput>,
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
    let mut cleanup_receivers = Vec::new();

    let required_space = config.input_size * headroom;
    let ram_available = available_space_for_path(&config.ram_temp_dir).await.unwrap_or(0);
    let use_ram = required_space <= ram_available;

    let temp_dir = choose_temp_dir(
        config.input_size * headroom,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false
    ).await?;

    // 1. Deinterleave input to files for HISAT2 input
    let (r1_rx, r2_rx_opt, deint_handle) = deinterleave_fastq_stream(
        input_stream,
        paired,
        config.base_buffer_size,
    ).await.map_err(|e| PipelineError::Other(e))?;

    let r1_path = temp_dir.path().join("hisat2_input_r1.fastq");
    let r1_write_task = write_byte_stream_to_file(
        &r1_path,
        ReceiverStream::new(r1_rx),
        Some(config.base_buffer_size),
    ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
    cleanup_tasks.push(r1_write_task);

    let r2_path_opt = if paired {
        let r2_path = temp_dir.path().join("hisat2_input_r2.fastq");
        let r2_write_task = write_byte_stream_to_file(
            &r2_path,
            ReceiverStream::new(r2_rx_opt.unwrap()),
            Some(config.base_buffer_size),
        ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(r2_write_task);
        Some(r2_path)
    } else {
        None
    };

    deint_handle.await.map_err(|e| PipelineError::Other(e.into()))??;

    // 2. HISAT2 alignment to temp SAM
    let sam_path = temp_dir.path().join("hisat2.sam");

    let hisat2_config = Hisat2Config {
        hisat2_index_path: hisat2_index_path.clone(),
        option_fields: hisat2_options,
        r1_path: r1_path.to_string_lossy().to_string(),
        r2_path: r2_path_opt.as_ref().map(|p| p.to_string_lossy().to_string()),
    };

    let mut hisat2_args = generate_cli(HISAT2_TAG, &config, Some(&hisat2_config))
        .map_err(|e| PipelineError::ToolExecution { tool: HISAT2_TAG.to_string(), error: e.to_string() })?;

    hisat2_args.push("-S".to_string());
    hisat2_args.push(sam_path.to_string_lossy().to_string());
    info!("Hisat2 args: {:?}", hisat2_args);

    let (mut hisat2_child, hisat2_err_task) = spawn_cmd(
        config.clone(),
        HISAT2_TAG,
        hisat2_args.clone(),
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: HISAT2_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(hisat2_err_task);

    let hisat2_status = hisat2_child.wait().await.map_err(|e| PipelineError::Other(e.into()))?;
    if !hisat2_status.success() {
        error!("HISAT2 failed with exit code {}", hisat2_status.code().unwrap_or(-1));
        return Err(PipelineError::ToolExecution { tool: HISAT2_TAG.to_string(), error: "HISAT2 execution failed".to_string() });
    } else {
        info!("HISAT2 completed successfully");
    }

    let sam_size = tokio::fs::metadata(&sam_path).await?.len();
    info!("HISAT2 SAM file size: {} bytes (should be >0)", sam_size);
    if sam_size == 0 {
        error!("HISAT2 did not write to SAM file - check args: {:?}", hisat2_args.clone());
        return Err(PipelineError::Other(anyhow!("Empty SAM output from HISAT2")));
    }

    // 3. samtools sort -n SAM → BAM
    let bam_path = temp_dir.path().join("hisat2_sorted.bam");

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None),
            ("-o".to_string(), Some(bam_path.to_string_lossy().to_string())),
            ("-@".to_string(), Some("8".to_string())),
            ("-l".to_string(), Some("1".to_string())),
            ("-T".to_string(), Some(temp_dir.path().to_string_lossy().to_string())),
        ]),
    };

    let mut samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution { tool: SAMTOOLS_TAG.to_string(), error: e.to_string() })?;

    samtools_sort_args.push(sam_path.to_string_lossy().to_string());

    let (mut sort_child, sort_err_task) = spawn_cmd(
        config.clone(),
        SAMTOOLS_TAG,
        samtools_sort_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(sort_err_task);

    sort_child.wait().await.map_err(|e| PipelineError::Other(e.into()))?;

    if let Some(bam_out) = output_bam_path {
        tokio::fs::copy(&bam_path, &bam_out).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
    }

    // 4. Pipe: sorted BAM → fixmate -m → fastq -f 13/-f 4 (to files)
    let fq1_path = temp_dir.path().join("hisat2_host_filtered1.fastq");


    let fq2_path_opt = if paired {
        let fq2 = temp_dir.path().join("hisat2_host_filtered2.fastq");
        Some(fq2)
    } else {
        None
    };

    // Fixmate config and args
    let fixmate_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fixmate,
        subcommand_fields: HashMap::from([
            ("-m".to_string(), None),
        ]),
    };

    let mut fixmate_args = generate_cli(SAMTOOLS_TAG, &config, Some(&fixmate_config))
        .map_err(|e| PipelineError::ToolExecution { tool: SAMTOOLS_TAG.to_string(), error: e.to_string() })?;

    fixmate_args.push(bam_path.to_string_lossy().to_string());  // input
    fixmate_args.push("-".to_string());  // stdout

    let (mut fixmate_child, fixmate_err_task) = spawn_cmd(
        config.clone(),
        SAMTOOLS_TAG,
        fixmate_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution { tool: SAMTOOLS_TAG.to_string(), error: e.to_string() })?;
    cleanup_tasks.push(fixmate_err_task);

    // Parse fixmate stdout (BAM) as stream
    let fixmate_out_stream = parse_child_output(
        &mut fixmate_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await.map_err(|e| PipelineError::ToolExecution { tool: SAMTOOLS_TAG.to_string(), error: e.to_string() })?;

    // Fastq config and args
    let mut fastq_subcommand_fields: HashMap<String, Option<String>> = HashMap::new();
    fastq_subcommand_fields.insert("-f".to_string(), Some(if paired { "13".to_string() } else { "4".to_string() }));
    if paired {
        fastq_subcommand_fields.insert("-1".to_string(), Some(fq1_path.to_string_lossy().to_string()));
        fastq_subcommand_fields.insert("-2".to_string(), Some(fq2_path_opt.as_ref().unwrap().to_string_lossy().to_string()));
        fastq_subcommand_fields.insert("-0".to_string(), Some("/dev/null".to_string()));
        fastq_subcommand_fields.insert("-s".to_string(), Some("/dev/null".to_string()));
    } else {
        fastq_subcommand_fields.insert("-o".to_string(), Some(fq1_path.to_string_lossy().to_string()));
    }

    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: fastq_subcommand_fields,
    };

    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution { tool: SAMTOOLS_TAG.to_string(), error: e.to_string() })?;

    let (mut fastq_child, fastq_stdin_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        fixmate_out_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution { tool: SAMTOOLS_TAG.to_string(), error: e.to_string() })?;
    cleanup_tasks.push(fastq_stdin_task);
    cleanup_tasks.push(fastq_err_task);

    {
        let mut fastq_guard = fastq_child.lock().await;
        fastq_guard.wait().await.map_err(|e: std::io::Error| PipelineError::Other(e.into()))?;
    }
    // Debug: Count lines in FASTQ to check for uneven R1/R2 (4 lines per read)
    let r1_lines = Command::new("wc").arg("-l").arg(&fq1_path).output().await?;
    let r1_line_str = String::from_utf8_lossy(&r1_lines.stdout).trim().to_string();
    let r1_records = r1_line_str.parse::<u64>().unwrap_or(0) / 4;

    let (r2_records, r2_line_str) = if let Some(ref fq2_path) = fq2_path_opt {
        let r2_lines = Command::new("wc").arg("-l").arg(fq2_path).output().await?;
        let r2_str = String::from_utf8_lossy(&r2_lines.stdout).trim().to_string();
        (r2_str.parse::<u64>().unwrap_or(0) / 4, r2_str)
    } else {
        (0, "0".to_string())
    };

    info!("FASTQ debug: R1 lines={}, est. records={}; R2 lines={}, est. records={}", r1_line_str, r1_records, r2_line_str, r2_records);

    if paired && r1_records != r2_records {
        warn!("Uneven FASTQ records: R1={}, R2={} — possible truncation in samtools fastq", r1_records, r2_records);
    }

    // Debug FASTQ sizes
    let fq1_size = tokio::fs::metadata(&fq1_path).await?.len();
    info!("Extracted FASTQ1 size: {} bytes (~{} reads est.)", fq1_size, fq1_size / 600);
    if let Some(ref fq2_path) = fq2_path_opt {
        let fq2_size = tokio::fs::metadata(fq2_path).await?.len();
        info!("Extracted FASTQ2 size: {} bytes (~{} reads est.)", fq2_size, fq2_size / 600);
    }

    // 5. Read FASTQ files into stream
    let (unmapped_receiver, read_task) = read_fastq(
        fq1_path,
        fq2_path_opt,
        None,
        u64::MAX,
        None,
        None,
        config.base_buffer_size,
    ).map_err(|e| PipelineError::Other(e))?;

    let unmapped_stream = ReceiverStream::new(unmapped_receiver);

    let read_task_wrapped = tokio::spawn(async move {
        match read_task.await {
            Ok(Ok(stats)) => {
                info!("read_fastq completed with stats: {:?}", stats);
                Ok(())
            }
            Ok(Err(e)) => Err(e),
            Err(join_err) => Err(anyhow::anyhow!("read_task join failed: {}", join_err)),
        }
    });
    cleanup_tasks.push(read_task_wrapped);

    // 6. Count mapped reads from BAM
    let (count_tx, count_rx) = oneshot::channel::<u64>();

    let bam_count_task = tokio::spawn(async move {
        let flag = if paired { "-F13" } else { "-F4" };
        let count_args = vec!["view".to_string(), "-c".to_string(), flag.to_string(), bam_path.to_string_lossy().to_string()];
        let output = Command::new("samtools")
            .args(&count_args)
            .output()
            .await
            .map_err(|e| anyhow!("samtools count failed: {}", e))?;
        let count_str = String::from_utf8_lossy(&output.stdout).trim().to_string();
        let count: u64 = count_str.parse().unwrap_or(0);
        let _ = count_tx.send(count);
        Ok(())
    });
    cleanup_tasks.push(bam_count_task);
    let mut temp_dirs = vec![temp_dir];

    Ok((
        unmapped_stream,
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        temp_dirs
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
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<(ReceiverStream<ParseOutput>,  Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, oneshot::Receiver<Result<u64, anyhow::Error>>), PipelineError>{
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let qc_fastp_config_view = FastpConfig {

        //These default QC thresholds are loosely based on the pre-2022 CZI pipeline using PriceSeq & LZW
        command_fields: HashMap::from([
            ("--dont_eval_duplication".to_string(), None),
            ("--length_required".to_string(), Some("35".to_string())),
            ("--qualified_quality_phred".to_string(), Some("17".to_string())),
            ("--unqualified_percent_limit".to_string(), Some("15".to_string())),
            ("--n_base_limit".to_string(), Some("15".to_string())),
            ("--complexity_threshold".to_string(), Some("60".to_string())),
        ]),
    };

    let qc_fastp_args = generate_cli(FASTP_TAG, &config, Some(&qc_fastp_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut qc_fastp_child, qc_fastp_stream_task, qc_fastp_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        FASTP_TAG,
        qc_fastp_args,
        StreamDataType::JustBytes,
        config.args.verbose,
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
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Tee the byte stream for counting
    let tee_count_input = ReceiverStream::new(qc_fastp_out_stream);
    let (mut tee_count_streams, tee_count_done_rx) = t_junction(
        tee_count_input,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "fastp_count".to_string(),
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "t_junction".to_string(),
            error: e.to_string(),
        })?;
    cleanup_receivers.push(tee_count_done_rx);

    let mut tee_count_streams_iter = tee_count_streams.into_iter();
    let qc_out_stream = tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let qc_count_stream = tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (count_tx, count_result_rx) = oneshot::channel::<Result<u64, anyhow::Error>>();

    let count_task = tokio::spawn(async move {
        match stream_record_counter(qc_count_stream, false).await {
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

    Ok((ReceiverStream::new(qc_out_stream), cleanup_tasks, cleanup_receivers, count_result_rx))
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
    let mut cleanup_receivers = Vec::new();

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

    let kallisto_index = config.args.kallisto_index.clone()
        .ok_or(PipelineError::MissingArgument("kallisto_index required".to_string()))?;

    // Convert raw byte stream to FASTQ
    let byte_rx = input_stream.into_inner();
    let byte_reader = ChannelReader::new(byte_rx);
    let fastq_rx = parse_fastq(byte_reader, config.base_buffer_size).await
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
    let mut cleanup_receivers = Vec::new();
    let mut temp_dirs: Vec<TempDir> = Vec::new();
    let mut named_temp_files: Vec<NamedTempFile> = Vec::new();

    let temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false
    ).await?;


    let (_host_ref_fasta_path, host_index_path, host_ref_temp, host_index_temp, host_temp_dir,  mut host_ref_tasks) =
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


    let (mut minimap2_child, minimap2_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::IlluminaFastq, // Assuming interleaved FASTQ
        config.args.verbose,
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
            config.base_buffer_size,
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

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
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
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };
    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);
    let (bam_streams, bam_done_rx) = t_junction(
        bam_rx_stream,
        num_tees,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "minimap2_bam_split".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
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

    let (mut count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        bam_count_stream,
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (count_tx, count_rx) = oneshot::channel::<u64>();
    let count_future = tokio::spawn(async move {
        let mut guard = count_child_arc.lock().await;
        let count_lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        count_tx.send(mapped_count).map_err(|_| anyhow!("Failed to send mapped count"))?;
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Extract unmapped FASTQ
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
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

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
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
            config.base_buffer_size,
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


/// Combined de-duplication and subsambler.
/// Uses and A-Res algorithm and hashing to creat a min heap BinaryHeap
/// Runs in memory when possible, but has a disk option when
/// file size and memory constraints demand.
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Strucutred FASTQ stream
/// * 'paired` - Bool controlling single or paired end data
/// * `seed` - Seed value for rng allowing possibility of deterministic output
/// * `max_subsample` - Limit of reads for subsampling
/// * `prefix_len` - Length of prefix subsample of read
/// * `out_dir` - Output directory for intermediate files
///
/// # Returns
/// Result of:
/// Tuple of:
/// Output stream
/// Count stream
/// Vec of tkio tasks
/// Vec of tokio receivers
/// with Pipeline Error
async fn dedup_and_subsample(
    config: Arc<RunConfig>,
    input_stream: mpsc::Receiver<ParseOutput>,
    paired: bool,
    seed: u64,
    max_subsample: u64,
    prefix_len: Option<usize>,
    out_dir: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,
    oneshot::Receiver<u64>,
    ReceiverStream<ParseOutput>,
    Arc<HashMap<String, ClusterInfo>>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let mut count_map: HashMap<u64, u64> = HashMap::new();
    let mut dedup_map: HashMap<u64, (u64, Vec<SequenceRecord>)> = HashMap::new();
    let mut total_count = 0u64;

    let mut stream = ReceiverStream::new(input_stream);
    let mut first_chunk = true;
    let mut orphan_r1: Option<SequenceRecord> = None;
    let mut r1_count: u64 = 0;
    let mut r2_count: u64 = 0;

    while let Some(item) = stream.next().await {
        let record = match item {
            ParseOutput::Fastq(rec) => rec,
            ParseOutput::Bytes(bytes) => {
                if first_chunk {
                    if !bytes.starts_with(b"@") {
                        return Err(PipelineError::InvalidFastqFormat("Raw byte stream does not start with FASTQ header '@'".to_string()));
                    }
                    first_chunk = false;
                }

                let cursor = std::io::Cursor::new(&*bytes);
                let mut parser = needletail::parse_fastx_reader(cursor)
                    .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

                if let Some(rec_res) = parser.next() {
                    let rec = rec_res.map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;
                    SequenceRecord::from(rec)
                } else {
                    continue;
                }
            }
            _ => {
                warn!("Skipping unexpected ParseOutput variant in dedup stream");
                continue;
            }
        };

        if paired {
            if let Some(r1) = orphan_r1.take() {
                r1_count += 1;
                r2_count += 1;

                let mut records = vec![r1.clone(), record.clone()];
                let mut seq_bytes = r1.seq().to_vec();
                seq_bytes.extend_from_slice(record.seq());

                let mut hasher = XxHash64::with_seed(seed);
                let write_len = prefix_len.map_or(seq_bytes.len(), |l| l.min(seq_bytes.len()));
                hasher.write(&seq_bytes[0..write_len]);
                let hash_key = hasher.finish();

                let entry = dedup_map.entry(hash_key).or_insert((0, records.clone()));
                entry.0 += 1;
                entry.1 = records;

                total_count += 1;
            } else {
                // Store as pending R1
                orphan_r1 = Some(record);
            }
        } else {
            let mut records = vec![record.clone()];
            let seq_bytes = record.seq().to_vec();

            let mut hasher = XxHash64::with_seed(seed);
            let write_len = prefix_len.map_or(seq_bytes.len(), |l| l.min(seq_bytes.len()));
            hasher.write(&seq_bytes[0..write_len]);
            let hash_key = hasher.finish();

            let entry = dedup_map.entry(hash_key).or_insert((0, records.clone()));
            entry.0 += 1;
            entry.1 = records;

            total_count += 1;
        }
    }

    // Handle any final orphan R1 (e.g., uneven small file)
    if let Some(orphan) = orphan_r1 {
        warn!("Orphan R1 at end of paired stream: {} (skipped)", orphan.id());
    }

    // Fail on severe imbalance (real corruption)
    let imbalance = (r1_count as i64 - r2_count as i64).abs();
    if imbalance > 1 {
        error!("Severe paired stream imbalance: R1={}, R2={} (diff={})", r1_count, r2_count, imbalance);
        return Err(PipelineError::InvalidFastqFormat(format!(
            "Paired stream imbalance: R1={}, R2={}, diff={}",
            r1_count, r2_count, imbalance
        )));
    }

    info!(
        "Dedup complete: {} valid pairs processed (R1={}, R2={}), {} unique hashes",
        r1_count.min(r2_count),
        r1_count,
        r2_count,
        dedup_map.len()
    );

    eprintln!("total count for input {}", total_count);

    if total_count == 0 {
        info!("No reads in dedup input — returning empty outputs");
        let (empty_tx, empty_rx) = mpsc::channel(1);
        drop(empty_tx);
        let (cluster_tx, cluster_rx) = mpsc::channel(1);
        drop(cluster_tx);
        let (count_tx, count_rx) = oneshot::channel();
        let _ = count_tx.send(0);
        return Ok((ReceiverStream::new(empty_rx), count_rx, ReceiverStream::new(cluster_rx), Arc::new(HashMap::new()), cleanup_tasks, cleanup_receivers));
    }

    count_map = dedup_map.iter().map(|(&k, &(w, _))| (k, w)).collect();

    let subsample_size = max_subsample.min(count_map.len() as u64);
    info!("Subsampling {} items (max_subsample={}, unique hashes={})", subsample_size, max_subsample, count_map.len());

    let mut heap: BinaryHeap<Reverse<SampleItem>> = BinaryHeap::new();
    let mut rng = config.rng.clone();

    if subsample_size < count_map.len() as u64 {
        for (&hash_key, &weight) in &count_map {
            let key = rng.random::<f64>().powf(1.0 / weight as f64);
            if heap.len() < subsample_size as usize {
                heap.push(Reverse(SampleItem { key, records: vec![], weight, hash_key: Some(hash_key) }));
            } else if key > heap.peek().unwrap().0.key {
                heap.pop();
                heap.push(Reverse(SampleItem { key, records: vec![], weight, hash_key: Some(hash_key) }));
            }
        }
    } else {
        for (&hash_key, &weight) in &count_map {
            heap.push(Reverse(SampleItem { key: 0.0, records: vec![], weight, hash_key: Some(hash_key) }));
        }
    }

    let mut sampled: Vec<SampleItem> = heap.into_iter().map(|r| r.0).collect();
    sampled.sort_by(|a, b| a.key.partial_cmp(&b.key).unwrap_or(Ordering::Equal));

    for item in sampled.iter_mut() {
        if let Some(hash_key) = item.hash_key {
            if let Some((_, records)) = dedup_map.get(&hash_key) {
                item.records = records.clone();
            }
        }
    }

    info!("Sampled {} items; first has {} records", sampled.len(), sampled.get(0).map_or(0, |i| i.records.len()));

    let mut duplicate_clusters: HashMap<String, ClusterInfo> = HashMap::with_capacity(dedup_map.len());
    for (&_hash_key, &(size, ref records)) in &dedup_map {
        if records.is_empty() { continue; }
        let rep_id = records[0].id().to_string();
        let members: Vec<String> = records.iter().map(|r| r.id().to_string()).collect();

        duplicate_clusters.insert(rep_id, ClusterInfo { size, members });
    }

    let (cluster_tx, cluster_rx) = mpsc::channel(config.base_buffer_size);

    tokio::spawn(async move {
        for (&_hash_key, &(weight, ref records)) in &dedup_map {
            if records.is_empty() { continue; }
            let rep_id = records[0].id();
            let rep_line = format!("{},{}\n", rep_id, rep_id);
            let _ = cluster_tx.send(ParseOutput::Bytes(Arc::new(rep_line.into_bytes()))).await;
            for rec in records.iter().skip(1) {
                let line = format!("{},{}\n", rec.id(), rep_id);
                let _ = cluster_tx.send(ParseOutput::Bytes(Arc::new(line.into_bytes()))).await;
            }
        }
    });

    let tsv_path = out_dir.join("duplicate_cluster_sizes.tsv");
    let duplicate_clusters_for_tsv = duplicate_clusters.clone();

    let tsv_write_task = tokio::spawn(async move {
        let mut file = TokioOpenOptions::new()
            .write(true)
            .create(true)
            .open(&tsv_path)
            .await
            .map_err(|e| anyhow!("TSV open error: {}", e))?;

        for (rep_id, cluster) in duplicate_clusters_for_tsv.iter() {
            file.write_all(format!("{}\t{}\n", rep_id, cluster.size).as_bytes())
                .await
                .map_err(|e| anyhow!("TSV write error: {}", e))?;
        }

        Ok(())
    });
    cleanup_tasks.push(tsv_write_task);

    let unique_count = sampled.len() as u64;
    eprintln!("unique count = {}", unique_count);
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(unique_count).map_err(|_| PipelineError::Other(anyhow!("Failed to send subsample count")))?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let send_task = tokio::spawn(async move {
        for item in sampled {
            for record in item.records {
                tx.send(ParseOutput::Fastq(record)).await.map_err(|_| anyhow!("Send failed"))?;
            }
        }
        Ok(())
    });
    cleanup_tasks.push(send_task);

    Ok((
        ReceiverStream::new(rx),
        count_rx,
        ReceiverStream::new(cluster_rx),
        Arc::new(duplicate_clusters),
        cleanup_tasks,
        cleanup_receivers,
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
    ReceiverStream<ParseOutput>,  // cluster stream (for CSV-like)
    Arc<HashMap<String, ClusterInfo>>,  // duplicate_clusters
    Vec<JoinHandle<Result<(), anyhow::Error>>>,  // cleanup tasks
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,  // cleanup rx
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let mut dedup_map: HashMap<u64, (u64, Vec<SequenceRecord>)> = HashMap::new();
    let mut total_count = 0u64;

    // ───────────────────────────────────────────────────────────────
    // Debug: counters + hash key log
    // ───────────────────────────────────────────────────────────────
    let mut fastq_count: u64 = 0;           // Total Fastq records seen
    let mut pair_count: u64 = 0;            // Successful pairs formed
    let mut hash_keys: Vec<u64> = Vec::with_capacity(1_000_000);  // Pre-alloc

    let mut stream = ReceiverStream::new(input_stream);
    let mut orphan_r1: Option<SequenceRecord> = None;
    let mut r1_count: u64 = 0;
    let mut r2_count: u64 = 0;

    while let Some(item) = stream.next().await {
        let record = match item {
            ParseOutput::Fastq(rec) => rec,
            _ => continue,
        };

        fastq_count += 1;  // 1. Count every Fastq record

        if paired {
            if let Some(r1) = orphan_r1.take() {

                if r1.id() != record.id() {
                    return Err(PipelineError::InvalidFastqFormat(format!(
                        "Mismatched pair IDs: R1={}, R2={}",
                        r1.id(), record.id()
                    )));
                }

                pair_count += 1;  // 2. Count successful pairing

                let r1_prefix = &r1.seq()[..prefix_len.map_or(r1.seq().len(), |l| l.min(r1.seq().len()))]; //prefix length limitation
                let r2_prefix = &record.seq()[..prefix_len.map_or(record.seq().len(), |l| l.min(record.seq().len()))];

                let mut hasher = DefaultHasher::new();
                hasher.write(r1_prefix);
                hasher.write(&[0u8]);
                hasher.write(r2_prefix);
                let hash_key = hasher.finish();

                hash_keys.push(hash_key);  // 3. Collect every hash key

                let entry = dedup_map.entry(hash_key).or_insert((0, vec![]));
                entry.0 += 1;
                entry.1.push(r1.clone());
                entry.1.push(record.clone());

                total_count += 1;
            } else {
                orphan_r1 = Some(record);
            }
        } else {
            let seq_bytes = record.seq().to_vec();

            let mut hasher = DefaultHasher::new();
            let write_len = prefix_len.map_or(seq_bytes.len(), |l| l.min(seq_bytes.len()));
            hasher.write(&seq_bytes[0..write_len]);
            let hash_key = hasher.finish();

            hash_keys.push(hash_key);  // 3. Collect every hash key

            let entry = dedup_map.entry(hash_key).or_insert((0, vec![]));
            entry.0 += 1;
            entry.1.push(record.clone());

            total_count += 1;
        }
    }

    // Log counts to console immediately (before any failure)
    eprintln!("DEBUG: Total Fastq records seen: {}", fastq_count);
    eprintln!("DEBUG: Successful pairs formed: {}", pair_count);
    eprintln!("DEBUG: Total hash keys collected: {}", hash_keys.len());

    // Write hash keys to file (async, no blocking)
    let hash_path = out_dir.join("dedup_hash_keys.txt");
    let hash_keys_clone = hash_keys.clone();  // Clone for move
    let hash_write_task = tokio::spawn(async move {
        let mut file = BufWriter::new(TokioFile::create(&hash_path).await?);
        for key in hash_keys_clone {
            file.write_all(format!("{}\n", key).as_bytes()).await?;
        }
        file.flush().await?;
        Ok(())
    });
    cleanup_tasks.push(hash_write_task);

    // Align orphan/imbalance: Error on >1 diff (no skip, like czid-dedup TryFrom)
    if orphan_r1.is_some() || (r1_count as i64 - r2_count as i64).abs() > 1 {
        return Err(PipelineError::InvalidFastqFormat(format!(
            "Paired imbalance: R1={}, R2={}, orphan={}",
            r1_count, r2_count, orphan_r1.is_some() as u8
        )));
    }

    info!("Dedup complete: {} valid pairs, {} unique hashes", total_count, dedup_map.len());

    let mut duplicate_clusters: HashMap<String, ClusterInfo> = HashMap::with_capacity(dedup_map.len());
    for (&_hash, &(size, ref records)) in &dedup_map {
        if records.is_empty() { continue; }
        let rep_id = records[0].id().to_string();
        let members: Vec<String> = records.iter().map(|r| r.id().to_string()).collect();
        duplicate_clusters.insert(rep_id, ClusterInfo { size, members });
    }

    // Write cluster sizes TSV
    let tsv_path = out_dir.join("duplicate_cluster_sizes.tsv");
    let duplicate_clusters_for_tsv = duplicate_clusters.clone();
    let tsv_task = tokio::spawn(async move {
        let mut file = BufWriter::new(TokioFile::create(&tsv_path).await?);
        for (rep_id, cluster) in duplicate_clusters_for_tsv.iter() {
            file.write_all(format!("{}\t{}\n", rep_id, cluster.size).as_bytes()).await?;
        }
        file.flush().await?;
        Ok(())
    });
    cleanup_tasks.push(tsv_task);

    // Uniques stream: First rep per cluster
    let (uniques_tx, uniques_rx) = mpsc::channel(config.base_buffer_size);
    let value = dedup_map.clone();
    let uniques_task = tokio::spawn(async move {
        for (&_hash, &(_, ref records)) in &value {
            if !records.is_empty() {
                uniques_tx.send(ParseOutput::Fastq(records[0].clone())).await?;
            }
        }
        Ok(())
    });
    cleanup_tasks.push(uniques_task);

    // Cluster stream
    let (cluster_tx, cluster_rx) = mpsc::channel(config.base_buffer_size);
    let dedup_map_for_cluster = dedup_map.clone();
    let cluster_task = tokio::spawn(async move {
        for (&_hash, &(_, ref records)) in &dedup_map_for_cluster {
            if records.is_empty() { continue; }
            let rep_id = records[0].id();
            let rep_line = format!("{},{}\n", rep_id, rep_id);
            cluster_tx.send(ParseOutput::Bytes(Arc::new(rep_line.into_bytes()))).await?;
            for rec in records.iter().skip(1) {
                let line = format!("{},{}\n", rec.id(), rep_id);
                cluster_tx.send(ParseOutput::Bytes(Arc::new(line.into_bytes()))).await?;
            }
        }
        Ok(())
    });
    cleanup_tasks.push(cluster_task);

    let unique_count = dedup_map.len() as u64;
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(unique_count).map_err(|_| PipelineError::Other(anyhow!("Count send failed")))?;

    Ok((
        ReceiverStream::new(uniques_rx),
        count_rx,
        ReceiverStream::new(cluster_rx),
        Arc::new(duplicate_clusters),
        cleanup_tasks,
        cleanup_receivers,
    ))
}

async fn subsample_weighted(
    config: Arc<RunConfig>,
    uniques_stream: mpsc::Receiver<ParseOutput>,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    seed: u64,
    max_subsample: u64,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, JoinHandle<Result<(), anyhow::Error>>), PipelineError> {
    let mut heap: BinaryHeap<Reverse<SampleItem>> = BinaryHeap::new();
    let mut rng = config.rng.clone();

    let mut stream = ReceiverStream::new(uniques_stream);
    let mut sampled_items = Vec::new();

    while let Some(item) = stream.next().await {
        let record = match item {
            ParseOutput::Fastq(rec) => rec,
            _ => continue,
        };

        let rep_id = record.id().to_string();
        let weight = duplicate_clusters.get(&rep_id).map_or(1, |c| c.size);

        let key = rng.random::<f64>().powf(1.0 / weight as f64);
        sampled_items.push(SampleItem { key, records: vec![record], weight, hash_key: None });
    }

    let subsample_size = max_subsample.min(sampled_items.len() as u64);
    for mut item in sampled_items.into_iter() {
        if heap.len() < subsample_size as usize {
            heap.push(Reverse(item));
        } else if item.key > heap.peek().unwrap().0.key {
            heap.pop();
            heap.push(Reverse(item));
        }
    }

    let mut sampled = heap.into_iter().map(|r| r.0).collect::<Vec<_>>();
    sampled.sort_by(|a, b| a.key.partial_cmp(&b.key).unwrap_or(Ordering::Equal));

    let count = sampled.len() as u64;
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(count).map_err(|_| PipelineError::Other(anyhow!("Count send failed")))?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let send_task: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
        for item in sampled {
            for rec in item.records {
                tx.send(ParseOutput::Fastq(rec)).await?;
            }
        }
        Ok(())
    });

    Ok((ReceiverStream::new(rx), count_rx, send_task))
}

async fn subsample_uniform(
    config: Arc<RunConfig>,
    uniques_stream: mpsc::Receiver<ParseOutput>,
    max_subsample: u64,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, JoinHandle<Result<(), anyhow::Error>>), PipelineError> {
    let mut records = Vec::new();
    let mut stream = ReceiverStream::new(uniques_stream);
    while let Some(item) = stream.next().await {
        if let ParseOutput::Fastq(rec) = item {
            records.push(rec);
        }
    }

    let subsample_size = max_subsample.min(records.len() as u64) as usize;

    // Uniform random sampling: shuffle indices and take first N
    let mut rng = config.rng.clone();
    let mut indices: Vec<usize> = (0..records.len()).collect();
    indices.shuffle(&mut rng);  // Shuffle all indices randomly
    let sampled_indices: Vec<usize> = indices.into_iter().take(subsample_size).collect();

    let count = sampled_indices.len() as u64;
    let (count_tx, count_rx) = oneshot::channel();
    count_tx.send(count).map_err(|_| PipelineError::Other(anyhow!("Count send failed")))?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let send_task: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
        for idx in sampled_indices {
            tx.send(ParseOutput::Fastq(records[idx].clone())).await?;
        }
        Ok(())
    });

    Ok((ReceiverStream::new(rx), count_rx, send_task))
}

/// Alignment of deduped/subsampled stream to a non-host DB (canoncially NT)
/// for output in PAF format, allowing downstream metagenomic calling.
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
async fn minimap2_non_host_align(
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
    let mut cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();

    let base_buffer_size = config.base_buffer_size;
    let channel_buffer = base_buffer_size * 4;

    // 2. Discover NT split .mmi chunks + sort largest-first (unchanged)
    let nt_split_dir = PathBuf::from(config.args.nt_split_dir.clone());
    let mut chunk_paths: Vec<PathBuf> = Vec::new();
    let mut max_size_bytes: u64 = 0;

    let mut entries = tokio::fs::read_dir(&nt_split_dir)
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Cannot read NT split dir: {}", e)))?;

    while let Some(entry) = entries.next_entry().await.transpose() {
        if let Ok(entry) = entry {
            let path = entry.path();
            if path.is_file() && path.extension().map_or(false, |ext| ext == "mmi") {
                if let Ok(meta) = path.metadata() {
                    let size = meta.len();
                    if size > max_size_bytes {
                        max_size_bytes = size;
                    }
                }
                chunk_paths.push(path);
            }
        }
    }

    if chunk_paths.is_empty() {
        return Err(PipelineError::Other(anyhow!("No .mmi chunk files found in {}", nt_split_dir.display())));
    }

    chunk_paths.sort_by_key(|p| std::cmp::Reverse(p.metadata().map(|m| m.len()).unwrap_or(0)));

    let num_chunks = chunk_paths.len();
    info!("Found {} minimap2 NT index chunks (sorted largest first)", num_chunks);

    // 3. Estimate memory per job (unchanged)
    const SAFETY_MULTIPLIER: f64 = 0.95;
    const MIN_ESTIMATE_GB: u64 = 20;

    let largest_gb = if max_size_bytes > 0 {
        ((max_size_bytes + (1 << 30) - 1) / (1 << 30)) as u64
    } else {
        70
    };

    let est_gb_per_job = ((largest_gb as f64 * SAFETY_MULTIPLIER).ceil() as u64).max(MIN_ESTIMATE_GB);

    info!("Largest NT chunk: {} GiB → estimating {:.0} GiB RAM per minimap2 job", largest_gb, est_gb_per_job);

    // 4. Compute concurrency & threads-per-job (unchanged)
    const MIN_THREADS_PER_JOB: usize = 4;
    const MAX_THREADS_PER_JOB: usize = 12;


    let concurrency = compute_phase_concurrency(
        &config,
        "minimap2_nt",
        est_gb_per_job as f64,
        6.0,
        20,
        4,
    );

    let threads_per_job = (config.max_cores / concurrency)
        .clamp(MIN_THREADS_PER_JOB, MAX_THREADS_PER_JOB);

    info!(
    "NT minimap2 stage: using {} concurrent jobs × {} threads (est peak ~{} GiB)",
    concurrency, threads_per_job, concurrency as u64 * est_gb_per_job as u64
);


    // ────────────────────────────────────────────────────────────────
    // 5. Rayon thread pool — enforces exact concurrency
    // ────────────────────────────────────────────────────────────────
    let rayon_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(concurrency)
        .thread_name(|i| format!("mm2-rayon-{}", i))
        .build()
        .map_err(|e| PipelineError::Other(anyhow!("Rayon pool failed: {}", e)))?;

    info!("Rayon pool ready with {} threads for minimap2 chunks", concurrency);

    // Prepare work items
    let work_items: Vec<_> = chunk_paths.into_iter().enumerate().collect();

    // Parallel execution — each rayon thread runs one minimap2 job
    let results: Vec<Result<(tokio::sync::mpsc::Receiver<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>), anyhow::Error>> =
        rayon_pool.install(|| {
            work_items.par_iter().map(|(idx, chunk_mmi)| -> Result<_, anyhow::Error> {
                let chunk_name = chunk_mmi.file_name()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| format!("chunk_{}", idx));

                let mut options = std::collections::HashMap::new();
                options.insert("-c".to_string(), None);
                options.insert("-x".to_string(), Some("sr".to_string()));
                options.insert("--secondary".to_string(), Some("yes".to_string()));

                let mm_config = Minimap2Config {
                    minimap2_index_path: chunk_mmi.clone(),
                    r1_path: Some(r1_path.clone()),
                    r2_path: r2_path_opt.clone(),
                    option_fields: options,
                    num_threads: Some(threads_per_job),
                };

                let args = generate_cli(MINIMAP2_TAG, &config, Some(&mm_config))
                    .context("Failed to generate minimap2 args")?;

                info!("Rayon launching minimap2 for {}: minimap2 {}", chunk_name, args.join(" "));

                // Create per-thread Tokio runtime to allow .await inside sync rayon worker
                let rt = tokio::runtime::Builder::new_current_thread()
                    .enable_all()
                    .build()
                    .context("Failed to create per-thread tokio runtime")?;

                // Run async spawn_cmd and parse inside the runtime
                let (mut child, stderr_task) = rt.block_on(spawn_cmd(
                    config.clone(),
                    MINIMAP2_TAG,
                    args,
                    config.args.verbose,
                ))?;

                if let Some(pid) = child.id() {
                    info!("minimap2 PID for {}: {}", chunk_name, pid);
                }

                let paf_receiver = rt.block_on(parse_child_output(
                    &mut child,
                    ChildStream::Stdout,
                    ParseMode::Lines,
                    config.base_buffer_size,
                ))?;

                info!("Chunk {} finished", chunk_name);

                Ok((paf_receiver, stderr_task))
            }).collect()
        });

    // 6. Convert results to ReceiverStream
    let mut partial_paf_receivers = Vec::new();

    for res in results {
        let (paf_rx, stderr_task) = res?;
        partial_paf_receivers.push(ReceiverStream::new(paf_rx));
        cleanup_tasks.push(stderr_task);
    }

    // 7. Merge all PAF streams
    let (merged_tx, merged_rx) = mpsc::channel(channel_buffer);

    let gather_handle = tokio::spawn(PafRecord::merge_paf_streams(
        partial_paf_receivers,
        merged_tx,
        channel_buffer,
    ));

    cleanup_tasks.push(gather_handle);

    Ok((
        ReceiverStream::new(merged_rx),
        cleanup_tasks,
        cleanup_receivers,
    ))
}


    // ────────────────────────────────────────────────────────────────
    //   TEMP BYPASS: SKIP REAL MINIMAP2 SCATTER-GATHER TO TEST OOM
    //   (comment out or delete this whole block when done debugging)
    // ────────────────────────────────────────────────────────────────

    // ────────────────────────────────────────────────────────────────
    //   TEMP BYPASS: SKIP REAL MINIMAP2 SCATTER-GATHER TO TEST OOM
    //   Comment out or delete this whole block when finished debugging
    // ────────────────────────────────────────────────────────────────

//     info!("=== MINIMAP2 SCATTER-GATHER BYPASSED FOR OOM DEBUG ===");
//     info!("Returning fake single-line ParseOutput stream instead of real PAF");
//
//     use tokio::sync::mpsc;
//
//     let (fake_tx, fake_rx) = mpsc::channel::<ParseOutput>(4);
//
//     // Create a fake ParseOutput::Bytes containing one valid-looking PAF line
//     let fake_paf_content = b"read_000001\t150\t0\t150\t+\tNT_dummy_accession\t200\t10\t160\t140\t93\t0\n".to_vec();
//
//     fake_tx.send(ParseOutput::Bytes(fake_paf_content.into()))
//         .await
//         .expect("Failed to send fake ParseOutput in debug bypass");
//
//     drop(fake_tx);  // Close the channel → downstream sees EOF
//
//     let fake_paf_stream = ReceiverStream::new(fake_rx);
//
//     info!("Fake ParseOutput::Bytes PAF stream created — pipeline continues without real minimap2");
//
//     Ok((
//         fake_paf_stream,
//         cleanup_tasks,
//         cleanup_receivers,
//         Some(temp_dir),
//     ))
// }

// ────────────────────────────────────────────────────────────────
    //   End of temporary OOM-test bypass
    // ────────────────────────────────────────────────────────────────

/// Converts a PAF stream to an m8 stream, splits it with t_junction for file writing,
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
    mut input_stream: ReceiverStream<ParseOutput>,
    m8_path: PathBuf,
) -> Result<(
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
)> {
    let (m8_tx, m8_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);

    let config_clone = config.clone();

    let conversion_handle = tokio::spawn(async move {
        while let Some(item) = input_stream.next().await {
            let line_bytes = match item {
                ParseOutput::Bytes(b) => b,
                _ => {
                    debug!("Skipping non-bytes item in PAF stream");
                    continue;
                }
            };

            let line = match String::from_utf8((*line_bytes).clone()) {
                Ok(s) => s,
                Err(e) => {
                    debug!("Invalid UTF-8 in PAF line: {}", e);
                    continue;
                }
            };

            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            // eprintln!("paf: {}", line);
            let record = match PafRecord::parse_line(&line) {
                Ok(r) => r,
                Err(e) => {
                    debug!("Failed to parse PAF line: {} — {}", e, line);
                    continue;
                }
            };

            let genome_size = config_clone.args.nt_db_size as f64;
            let m8_line = record.to_m8_line(genome_size);
            // eprintln!("m8: {}", m8_line);
            // eprintln!("");
            let m8_bytes = (m8_line + "\n").into_bytes();

            if m8_tx.send(ParseOutput::Bytes(Arc::new(m8_bytes))).await.is_err() {
                return Err(anyhow!("m8 channel send failed"));
            }
        }
        Ok(())
    });

    let m8_stream = ReceiverStream::new(m8_rx);

    let (mut streams, done_rx) = t_junction(
        m8_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "paf_to_m8_stream".to_string(),
        None,
    )
        .await?;

    let main_stream = streams.remove(0);
    let file_stream = streams.remove(0);

    let write_task = write_byte_stream_to_file(
        &m8_path,
        ReceiverStream::new(file_stream),
        Some(config.base_buffer_size),
    )
        .await?;

    Ok((
        ReceiverStream::new(main_stream),
        vec![conversion_handle, write_task],
        vec![done_rx],
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

/// Calls hits with the same logic as call_hits_m8 in the CZI pipeline
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
pub async fn call_hits_m8_stream(
    config: Arc<RunConfig>,
    mut m8_input: ReceiverStream<ParseOutput>,
    sample_base_buf: PathBuf,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    acc2taxid_map: Arc<Map<Vec<u8>>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    min_aln_len: u64,
    concurrency: usize,
) -> Result<(
    ReceiverStream<ParseOutput>,
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let (dedup_tx, dedup_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);
    let (summary_tx, summary_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);

    // 1. Batching task — collect m8 lines into batches
    const BATCH_SIZE: usize = 50_000;

    let (batch_tx, mut batch_rx) = mpsc::unbounded_channel::<Vec<Vec<u8>>>();

    let batching_task = tokio::spawn(async move {
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(BATCH_SIZE);
        let mut line_count = 0u64;

        while let Some(item) = m8_input.next().await {
            line_count += 1;
            let line_bytes = match item {
                ParseOutput::Bytes(b) => b.to_vec(),
                _ => continue,
            };

            if line_bytes.is_empty() || line_bytes[0] == b'#' {
                continue;
            }

            batch.push(line_bytes);

            if batch.len() >= BATCH_SIZE {
                let _ = batch_tx.send(std::mem::take(&mut batch));
                batch = Vec::with_capacity(BATCH_SIZE);
            }
        }

        if !batch.is_empty() {
            let _ = batch_tx.send(batch);
        }

        info!("Batching complete: {} lines", line_count);
        Ok(())
    });
    cleanup_tasks.push(batching_task);

    // 2. Drain all batches into a Vec (blocking drain in tokio task)
    let all_batches = {
        let mut batches = Vec::new();
        while let Some(batch) = batch_rx.recv().await {
            batches.push(batch);
        }
        batches
    };

    // 3. Rayon parallel processing of batches
    let rayon_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(concurrency)
        .thread_name(|i| format!("m8-worker-{}", i))
        .build()
        .map_err(|e| PipelineError::Other(anyhow!("Rayon pool failed: {}", e)))?;

    info!("Rayon pool ready with {} threads for m8 processing", concurrency);

    let processed_batches = rayon_pool.install(|| {
        all_batches.par_iter().map(|batch| {
            let mut local_dedup: Vec<String> = Vec::new();
            let mut local_summary: Vec<String> = Vec::new();

            let lineage_map_clone = lineage_map.clone();
            let acc2taxid_map_clone = acc2taxid_map.clone();
            let should_keep_clone = should_keep_filter.clone();

            let mut read_groups: HashMap<String, Vec<M8Record>> = HashMap::with_capacity(batch.len() / 10);

            for line_bytes in batch {
                let line = match String::from_utf8(line_bytes.clone()) {
                    Ok(s) => s.trim_end().to_string(),
                    Err(_) => continue,
                };

                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }

                let rec = match M8Record::parse_line_nr(&line) {
                    Ok(r) => r,
                    Err(e) => {
                        warn!("Failed to parse m8 line: {} — {}", e, line);
                        continue;
                    }
                };

                if rec.alen < min_aln_len {
                    continue;
                }

                read_groups.entry(rec.qname.clone()).or_default().push(rec);
            }

            for (read_id, mut hits) in read_groups {
                let mut valid_hits = Vec::new();
                for hit in &hits {
                    if let Some(taxid_u64) = acc2taxid_map_clone.get(hit.tname.as_bytes()) {
                        let taxid = taxid_u64 as i32;
                        if taxid == 0 {
                            continue;
                        }
                        let lineage = lineage_map_clone.get(&taxid).cloned().unwrap_or([-1; 3]);
                        if !(should_keep_clone)(&lineage) {
                            continue;
                        }
                        valid_hits.push(hit.clone());
                    }
                }

                if valid_hits.is_empty() {
                    continue;
                }

                valid_hits.sort_by(|a, b| b.bitscore.partial_cmp(&a.bitscore).unwrap_or(Ordering::Equal));
                let best = valid_hits[0].clone();

                let (tax_level, consensus_taxid, consensus_hits) = consensus_level(
                    &hits,
                    &*lineage_map_clone,
                    &*acc2taxid_map_clone,
                    &should_keep_clone,
                ).unwrap_or_else(|e| {
                    warn!("Consensus failed for {}: {}", read_id, e);
                    (0, 0, Vec::new())
                });

                let dedup_line = format!(
                    "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}",
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
                );
                local_dedup.push(dedup_line);

                let first_lineage = consensus_hits
                    .first()
                    .and_then(|h| acc2taxid_map_clone.get(h.tname.as_bytes()))
                    .and_then(|taxid_u64| {
                        let taxid = taxid_u64 as i32;
                        lineage_map_clone.get(&taxid)
                    })
                    .cloned()
                    .unwrap_or([-1; 3]);

                let species = first_lineage[0] as i64;
                let genus = first_lineage[1] as i64;
                let family = first_lineage[2] as i64;

                let summary_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    read_id, best.tname, species, genus, family, tax_level
                );
                local_summary.push(summary_line);
            }

            (local_dedup, local_summary)
        }).collect::<Vec<_>>()
    });

    // Merge batches sequentially into output channels (preserves order)
    let merge_task = tokio::spawn(async move {
        for (dedup_batch, summary_batch) in processed_batches {
            for line in dedup_batch {
                let _ = dedup_tx.send(ParseOutput::Bytes(Arc::new((line + "\n").into_bytes()))).await;
            }
            for line in summary_batch {
                let _ = summary_tx.send(ParseOutput::Bytes(Arc::new((line + "\n").into_bytes()))).await;
            }
        }
        Ok(())
    });
    cleanup_tasks.push(merge_task);

    // File writes (unchanged)
    let m8_dedup_file_path = config.out_dir.join(rename_file_path(&sample_base_buf, None, Some("dedup.m8"), "_"));
    let summary_file_path = config.out_dir.join(rename_file_path(&sample_base_buf, None, Some("summary.txt"), "_"));

    let dedup_stream = ReceiverStream::new(dedup_rx);
    let (mut dedup_branches, dedup_done) = t_junction(
        dedup_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "call_hits_m8_dedup".to_string(),
        None,
    ).await?;
    cleanup_receivers.push(dedup_done);
    let dedup_main = dedup_branches.remove(0);
    let dedup_file_stream = dedup_branches.remove(0);

    let call_file_write_task = write_byte_stream_to_file(
        &m8_dedup_file_path,
        ReceiverStream::new(dedup_file_stream),
        Some(config.base_buffer_size),
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    cleanup_tasks.push(call_file_write_task);

    let summary_stream = ReceiverStream::new(summary_rx);
    let (mut summary_branches, summary_done) = t_junction(
        summary_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "call_hits_m8_summary".to_string(),
        None,
    ).await?;
    cleanup_receivers.push(summary_done);
    let summary_main = summary_branches.remove(0);
    let summary_file_stream = summary_branches.remove(0);

    let summary_file_write_task = write_byte_stream_to_file(
        &summary_file_path,
        ReceiverStream::new(summary_file_stream),
        Some(config.base_buffer_size),
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    cleanup_tasks.push(summary_file_write_task);

    Ok((
        ReceiverStream::new(dedup_main),
        ReceiverStream::new(summary_main),
        cleanup_tasks,
        cleanup_receivers,
    ))
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
    let mut cleanup_receivers = Vec::new();

    // Temp dir
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

    // Safety cap for Epyc nodes
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
        ("--verbose".to_string(), None),
    ]);

    let diamond_config = DiamondConfig {
        subcommand: DiamondSubcommand::Blastx,
        db: db_prefix,
        r1_path: Some(r1_path.clone()),
        r2_path: r2_path_opt.clone(),
        subcommand_fields: diamond_options,
    };

    let diamond_args = generate_cli(DIAMOND_TAG, &config, Some(&diamond_config))
        .map_err(|e| PipelineError::ToolExecution { tool: DIAMOND_TAG.to_string(), error: e.to_string() })?;


    let (mut diamond_child, diamond_stderr_task) = spawn_cmd(
        config.clone(),
        DIAMOND_TAG,
        diamond_args,
        config.args.verbose,
    ).await?;

    // Take pipes upfront (before any move)
    let diamond_stdout = diamond_child.stdout.take().ok_or_else(|| anyhow!("diamond stdout missing"))?;
    let diamond_stderr = diamond_child.stderr.take().ok_or_else(|| anyhow!("diamond stderr missing"))?;

    // Parse stdout to m8 stream
    let (m8_tx, m8_rx) = mpsc::channel(config.base_buffer_size * 8);

    let parse_task = tokio::spawn(async move {
        let mut reader = TokioBufReader::new(diamond_stdout);
        let mut line = Vec::with_capacity(1024);

        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line).await? == 0 {
                break;
            }
            if !line.is_empty() && line[0] != b'#' {
                let arc_line = Arc::new(line.clone());
                if m8_tx.send(ParseOutput::Bytes(arc_line)).await.is_err() {
                    break;
                }
            }
        }
        Ok(())
    });

    // Wait for diamond exit
    let wait_task = tokio::spawn(async move {
        let status = diamond_child.wait().await?;
        if !status.success() {
            return Err(anyhow!("diamond exited with code {:?}", status.code()));
        }
        Ok(())
    });

    // Cleanup
    cleanup_tasks.push(parse_task);
    cleanup_tasks.push(wait_task);
    cleanup_tasks.push(diamond_stderr_task);  // your stderr logging

    let mut temp_dirs = vec![temp_dir];

    Ok((m8_rx, cleanup_tasks, cleanup_receivers, temp_dirs))
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
    mut m8_stream: ReceiverStream<ParseOutput>,
    mut summary_stream: ReceiverStream<ParseOutput>,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    count_type: String,               // e.g. "NT"
    source_count_type: Option<String>,// e.g. "NR" for the other DB
) -> Result<Vec<TaxonCount>, PipelineError> {

    // Expected format (tab-separated):
    // read_id   accession   species   genus   family   level
    let mut read_summaries: HashMap<String, (u8, String, Vec<i32>)> = HashMap::new();

    while let Some(item) = summary_stream.next().await {
        let line = match item {
            ParseOutput::Bytes(b) => {
                let s = String::from_utf8(b.to_vec())
                    .map_err(|e| PipelineError::Other(anyhow!("UTF-8 error: {}", e)))?;
                let s = s.trim_end();
                if s.is_empty() || s.starts_with('#') {
                    continue;
                }
                s.to_string()
            }
            _ => continue,
        };

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 6 {
            warn!("Skipping malformed summary line (expected 6 cols): {}", line);
            continue;
        }

        let read_id = cols[0].to_string();
        let accession = cols[1].to_string();
        let species: i32 = cols[2].parse().unwrap_or(0);
        let genus: i32 = cols[3].parse().unwrap_or(0);
        let family: i32 = cols[4].parse().unwrap_or(0);
        let level: u8 = cols[5].parse().unwrap_or(0);

        let consensus_taxid: i32 = match level {
            1 => species,
            2 => genus,
            3 => family,
            _ => 0,
        };

        let raw_lineage = vec![species, genus, family];
        let cleaned = validate_taxid_lineage(&raw_lineage, consensus_taxid, level);
        read_summaries.insert(read_id, (level, accession, cleaned));
    }
    info!("Loaded {} read summaries", read_summaries.len());


    // consume the **deduped m8** stream (one line per read)
    // We only need percent-identity, alignment-length and e-value.
    let mut read_metrics: HashMap<String, (f64, u64, f64)> = HashMap::new();

    while let Some(item) = m8_stream.next().await {
        let line = match item {
            ParseOutput::Bytes(b) => String::from_utf8(b.to_vec())
                .map_err(|e| anyhow!("m8 line not UTF-8: {}", e))?,
            _ => continue,
        };

        let rec = M8Record::parse_line_nr(&line)
            .map_err(|e| anyhow!("failed to parse m8 line: {}", e))?;
        read_metrics.insert(rec.qname, (rec.pident, rec.alen, rec.evalue));
    }
    info!("Loaded {} read metrics", read_metrics.len());


    let mut aggregation: HashMap<Vec<i32>, AggBucket> = HashMap::new();
    let mut base_count_per_read: HashMap<String, u64> = HashMap::new();

    for (read_id, (level, _accession, lineage)) in read_summaries {
        let cluster_size = duplicate_clusters
            .get(&read_id)
            .map(|cluster| cluster.size)
            .unwrap_or(1u64);

        // Metrics are guaranteed to exist because the m8 line came from the same read
        if let Some((pident, alen, mut raw_evalue)) = read_metrics.get(&read_id).cloned() {
            if !should_keep_filter(&lineage) {
                continue;
            }

            let evalue = if raw_evalue <= 0.0 || !raw_evalue.is_finite() {
                MIN_NORMAL_POSITIVE_DOUBLE
            } else if raw_evalue <= MIN_NORMAL_POSITIVE_DOUBLE {
                MIN_NORMAL_POSITIVE_DOUBLE
            } else {
                raw_evalue
            };
            let evalue_log10 = evalue.log10();


            if level == 1 {
                base_count_per_read.insert(read_id.clone(), alen * cluster_size);
            }

            // walk the lineagerootward
            let mut agg_key = lineage.clone();
            while !agg_key.is_empty() {
                let bucket = aggregation.entry(agg_key.clone()).or_default();

                bucket.nonunique_count += cluster_size;
                bucket.unique_count    += 1;
                bucket.base_count      += *base_count_per_read.get(&read_id).unwrap_or(&0);
                bucket.sum_percent_identity += pident;
                bucket.sum_alignment_length += alen as f64;
                bucket.sum_e_value += evalue_log10;

                if let Some(src) = &source_count_type {
                    bucket.source_count_type.insert(src.clone());
                }

                // chop the leafward rank
                agg_key.remove(0);
            }
        }
    }

    let mut taxon_counts = Vec::new();

    for (mut agg_key, bucket) in aggregation {
        if bucket.unique_count == 0 { continue; }

        // Number of ranks we still have in the key (3 = species, 2 = genus, 1 = family)
        let remaining = agg_key.len() as u8;
        let tax_level = 4 - remaining;               // 1 = species, 2 = genus, 3 = family
        let tax_id    = agg_key[0];

        // Helper to pull a rank with a negative fallback
        macro_rules! rank_or_neg {
            ($idx:expr, $neg:expr) => {
                if tax_level <= $idx { agg_key.get($idx as usize - tax_level as usize).copied().unwrap_or($neg) }
                else { $neg }
            };
        }
        let genus_taxid  = rank_or_neg!(2, -200);
        let family_taxid = rank_or_neg!(3, -300);

        let count = if READ_COUNTING_MODE == ReadCountingMode::CountAll {
            bucket.nonunique_count
        } else {
            bucket.unique_count
        };

        let dcr = bucket.nonunique_count as f64 / bucket.unique_count as f64;
        let percent_identity = bucket.sum_percent_identity / bucket.unique_count as f64;
        let alignment_length = bucket.sum_alignment_length / bucket.unique_count as f64;
        let e_value          = bucket.sum_e_value / bucket.unique_count as f64;

        let source_vec = if bucket.source_count_type.is_empty() {
            None
        } else {
            let mut v: Vec<String> = bucket.source_count_type.into_iter().collect();
            v.sort_unstable();
            Some(v)
        };

        taxon_counts.push(TaxonCount {
            tax_id,
            tax_level,
            genus_taxid,
            family_taxid,
            count,
            nonunique_count: bucket.nonunique_count,
            unique_count:    bucket.unique_count,
            dcr,
            percent_identity,
            alignment_length,
            e_value,
            count_type: count_type.clone(),
            base_count: bucket.base_count,
            source_count_type: source_vec,
        });
    }

    Ok(taxon_counts)
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


/// Generates an accession map from streamed m8 data
///
/// # Arguments
///
/// * `m8_stream` -stream in m8 format as defined in blast.rs
///
/// # Returns
/// Result of hashmap of id:accession
///
async fn collect_m8_to_accession_map(
    mut m8_stream: ReceiverStream<ParseOutput>,
) -> Result<HashMap<String, String>> {
    let mut map: HashMap<String, String> = HashMap::new();
    let mut item_count = 0;
    let start_time = Instant::now();

    while let Some(item) = m8_stream.next().await {
        item_count += 1;
        if let ParseOutput::Bytes(bytes) = item {
            let line = String::from_utf8_lossy(&*bytes).to_string();
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() >= 2 {
                let read_id = fields[0].to_string();
                let accession = fields[1].to_string();
                if map.insert(read_id.clone(), accession).is_some() {
                    warn!("Duplicate read_id in deduped m8: {}", read_id);
                }
            } else {
                warn!("Invalid m8 line: {}", line);
            }
        } else {
            return Err(anyhow!("Non-Bytes in m8 stream"));
        }

        if item_count % 100_000 == 0 {
            debug!("Processed {} m8 lines in {:?}", item_count, start_time.elapsed());
        }
    }

    info!("Collected {} accessions from m8 stream", map.len());
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
pub async fn generate_annotated_fasta(
    config: Arc<RunConfig>,
    mut dedup_stream: ReceiverStream<ParseOutput>,
    mut cluster_stream: ReceiverStream<ParseOutput>,
    nt_map: HashMap<String, String>,
    nr_map: HashMap<String, String>,
) -> Result<(
    mpsc::Receiver<ParseOutput>,               // mapped (annotated) contigs
    mpsc::Receiver<ParseOutput>,               // unidentified (all, including expanded duplicates)
    mpsc::Receiver<ParseOutput>,               // unique unidentified (representatives only)
    Vec<JoinHandle<Result<()>>>,               // cleanup tasks
    Vec<oneshot::Receiver<Result<()>>>,        // optional receivers
)> {
    // Channels for the three output streams
    let (mapped_tx, mapped_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size / 100); // ~10k records
    let (unidentified_tx, unidentified_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size / 100);
    let (unique_unidentified_tx, unique_unidentified_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size / 100);

    let (done_tx, done_rx) = tokio::sync::oneshot::channel();
    let mut cleanup_receivers = vec![done_rx];

    // Buffer of representative sequences: rep_id → (seq_vec, is_unidentified)
    let mut rep_buffer: HashMap<String, (Vec<u8>, bool)> = HashMap::with_capacity(1_000_000);
    let mut csv_line = String::with_capacity(256);

    let mut processed_reps = 0u64;
    let mut unidentified_reps = 0u64;
    let mut expanded_duplicates = 0u64;
    let start = tokio::time::Instant::now();

    let process_task = tokio::spawn(async move {
        while let Some(item) = dedup_stream.next().await {
            let record = match item {
                ParseOutput::Fastq(rec) | ParseOutput::Fasta(rec) => rec,
                _ => continue,
            };

            processed_reps += 1;
            if processed_reps % 100_000 == 0 {
                log::debug!("generate_annotated_fasta: processed {} reps", processed_reps);
            }

            let rep_id = record.id().to_string();
            let seq = record.seq().to_vec(); // Clone once
            let nr_acc = nr_map.get(&rep_id).cloned().unwrap_or_default();
            let nt_acc = nt_map.get(&rep_id).cloned().unwrap_or_default();
            let is_unidentified = nr_acc.is_empty() && nt_acc.is_empty();

            // Always buffer the representative
            rep_buffer.insert(rep_id.clone(), (seq.clone(), is_unidentified));

            // -----------------------------------------------------------------
            // Write representative to appropriate stream(s)
            // -----------------------------------------------------------------
            if !is_unidentified {
                // Annotated contig
                let header = format!(">NR:{}:NT:{}:{}\n", nr_acc, nt_acc, rep_id);
                let fasta = SequenceRecord::Fasta {
                    id: format!("NR:{}:NT:{}:{}", nr_acc, nt_acc, rep_id),
                    desc: None,
                    seq: Arc::new(seq),
                };
                mapped_tx
                    .send(ParseOutput::Fasta(fasta))
                    .await
                    .map_err(|_| anyhow!("mapped_tx dropped"))?;
            } else {
                // Unidentified representative → both unidentified streams
                unidentified_reps += 1;
                let header = format!("{}{}\n", UNMAPPED_HEADER_PREFIX, rep_id);
                let fasta = SequenceRecord::Fasta {
                    id: rep_id.clone(),
                    desc: None,
                    seq: Arc::new(seq.clone()),
                };

                // To full unidentified (all members)
                unidentified_tx
                    .send(ParseOutput::Fasta(fasta.clone()))
                    .await
                    .map_err(|_| anyhow!("unidentified_tx dropped"))?;

                // To unique unidentified (only reps)
                unique_unidentified_tx
                    .send(ParseOutput::Fasta(fasta))
                    .await
                    .map_err(|_| anyhow!("unique_unidentified_tx dropped"))?;
            }

            // -----------------------------------------------------------------
            // Expand duplicates for unidentified clusters
            // -----------------------------------------------------------------
            loop {
                let csv_opt = cluster_stream.next().await;
                let csv_item = match csv_opt {
                    Some(item) => item,
                    None => break,
                };

                if let ParseOutput::Bytes(bytes) = csv_item {
                    csv_line.clear();
                    csv_line.push_str(&String::from_utf8_lossy(&bytes));

                    if !csv_line.ends_with('\n') {
                        continue; // partial line
                    }

                    let line = std::mem::take(&mut csv_line);
                    let parts: Vec<&str> = line.trim().split(',').collect();
                    if parts.len() != 2 {
                        continue;
                    }

                    let member_id = parts[0];
                    let csv_rep_id = parts[1];

                    if csv_rep_id == rep_id && is_unidentified {
                        if let Some((seq, _)) = rep_buffer.get(&rep_id) {
                            // Preserve /1 /2 suffix if present
                            let (base, suffix) = if member_id.ends_with("/1") || member_id.ends_with("/2") {
                                let len = member_id.len();
                                (&member_id[..len - 2], &member_id[len - 2..])
                            } else {
                                (member_id, "")
                            };
                            let header = format!("{}{}{}\n", UNMAPPED_HEADER_PREFIX, base, suffix);
                            let fasta = SequenceRecord::Fasta {
                                id: format!("{}{}", base, suffix),
                                desc: None,
                                seq: Arc::new(seq.clone()),
                            };

                            unidentified_tx
                                .send(ParseOutput::Fasta(fasta))
                                .await
                                .map_err(|_| anyhow!("unidentified_tx dropped (duplicate)"))?;

                            expanded_duplicates += 1;
                        }
                    }
                }
            }
        }

        // Close all senders
        drop(mapped_tx);
        drop(unidentified_tx);
        drop(unique_unidentified_tx);

        log::info!(
            "generate_annotated_fasta: {} reps, {} unidentified ({} expanded duplicates) in {:?}",
            processed_reps,
            unidentified_reps,
            expanded_duplicates,
            start.elapsed()
        );

        let _ = done_tx.send(Ok(()));
        Ok(())
    });

    let mut cleanup_tasks = vec![process_task];

    Ok((
        mapped_rx,
        unidentified_rx,
        unique_unidentified_rx,
        cleanup_tasks,
        cleanup_receivers,
    ))
}

async fn generate_annotated_fasta_from_files(
    config: Arc<RunConfig>,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    cluster_stream: ReceiverStream<ParseOutput>,
    nt_map: HashMap<String, String>,
    nr_map: HashMap<String, String>,
) -> Result<(
    mpsc::Receiver<ParseOutput>,
    mpsc::Receiver<ParseOutput>,
    mpsc::Receiver<ParseOutput>,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
)> {
    let mut cleanup_tasks = Vec::new();

    // Parse R1 to raw receiver
    let r1_file = TokioFile::open(&r1_path).await
        .map_err(|e| anyhow!("Failed to open R1 {}: {}", r1_path.display(), e))?;
    let r1_receiver = parse_fastq(r1_file, config.base_buffer_size).await?;

    let interleaved_receiver = if let Some(r2_path) = r2_path_opt {
        let r2_file = TokioFile::open(&r2_path).await
            .map_err(|e| anyhow!("Failed to open R2 {}: {}", r2_path.display(), e))?;
        let r2_receiver = parse_fastq(r2_file, config.base_buffer_size).await?;

        // Interleave raw receivers
        let (inter_rx, inter_task) = interleave_fastq_streams(
            r1_receiver,
            r2_receiver,
            config.base_buffer_size,
        ).await?;

        cleanup_tasks.push(inter_task);
        inter_rx  // raw receiver
    } else {
        r1_receiver  // single-end
    };

    // Wrap the final receiver in ReceiverStream for generate_annotated_fasta
    let interleaved_stream = ReceiverStream::new(interleaved_receiver);

    // Call existing stream-based function
    generate_annotated_fasta(
        config,
        interleaved_stream,
        cluster_stream,
        nt_map,
        nr_map,
    ).await
}

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
    out_dir: &PathBuf,
) -> Result<JoinHandle<Result<()>>> {
    let assembly_dir = out_dir.join("assembly");
    fs::create_dir_all(&assembly_dir).await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create assembly dir: {}", e)))?;

    let spades_task = tokio::spawn(async move {
        let mut cleanup_tasks = Vec::new();

        // 1. Temp dir
        let est_temp_bytes = config.input_size + MAX_SPADES_WORK_DIR;
        let temp_dir = choose_temp_dir(
            est_temp_bytes,
            &config.ram_temp_dir,
            &config.args.nvme_scratch,
            4,
            true
        ).await?;


        // 3. Run SPAdes
        let spades_work_dir = TempDir::new_in(&temp_dir)?;
        let spades_work_path = spades_work_dir.path().to_path_buf();
        let stdout_log_path = assembly_dir.join("spades_stdout.log");

        let mut options = HashMap::new();
        options.insert("--only-assembler".to_string(), None);

        let spades_config = SpadesConfig {
            r1_path: r1_path.clone(),
            r2_path_opt: r2_path_opt.clone(),
            outdir_path: spades_work_path.clone(),
            option_fields: options,
        };

        let spades_args = generate_cli(SPADES_TAG, &config, Some(&spades_config))?;
        let (mut spades_child, spades_stderr_task) = spawn_cmd(
            config.clone(),
            SPADES_TAG,
            spades_args,
            config.args.verbose,
        ).await?;

        let spades_out_stream = parse_child_output(
            &mut spades_child,
            ChildStream::Stdout,
            ParseMode::Lines,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SPADES_TAG.to_string(),
                error: e.to_string(),
            })?;

        let spades_write_task = tokio::spawn(stream_to_file(spades_out_stream, stdout_log_path.clone()));
        cleanup_tasks.push(spades_write_task);

        // Await stderr task (drains and logs stderr)
        spades_stderr_task.await??;

        // Check SPAdes exit
        let spades_exit = spades_child.wait().await?;
        if !spades_exit.success() {
            error!("SPAdes failed with exit: {:?}", spades_exit);
            // Write dummies
            let contigs_path = assembly_dir.join("contigs.fasta");
            let contigs_all_path = assembly_dir.join("contigs_all.fasta");
            let scaffolds_path = assembly_dir.join("scaffolds.fasta");
            let sam_path = assembly_dir.join("read-contig.sam");
            let failed_marker = b";ASSEMBLY FAILED\n";
            tokio::fs::write(&contigs_path, failed_marker).await?;
            tokio::fs::write(&contigs_all_path, failed_marker).await?;
            tokio::fs::write(&scaffolds_path, failed_marker).await?;
            tokio::fs::write(&sam_path, b"@NO INFO\n").await?;
            return Err(anyhow!("SPAdes assembly failed"));
        }

        info!("SPAdes completed successfully");

        Ok(())
    });

    Ok(spades_task)
}

pub async fn process_assembly(
    config: Arc<RunConfig>,
    out_dir: &PathBuf,
    work_dir: &PathBuf,
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
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
    let mut cleanup_receivers = Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();

    let raw_contigs_path = work_dir.join("contigs.fasta");
    let raw_scaffolds_path = work_dir.join("scaffolds.fasta");
    let (empty_tx, empty_rx) = mpsc::channel::<ParseOutput>(1);

    let spades_success = async {
        match fs::metadata(&raw_contigs_path).await {
            Ok(meta) => meta.len() > 0 && meta.is_file(),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => false,
            Err(e) => {
                warn!("Failed to read contigs.fasta metadata: {}", e);
                false
            }
        }
    }
        .await;

    if !spades_success {
        let reason = if !raw_contigs_path.exists() {
            "contigs.fasta does not exist"
        } else {
            "contigs.fasta is empty"
        };

        warn!(
            "SPAdes assembly failed: {} — writing dummy outputs and skipping contig-based refinement",
            reason
        );

        let dummy = ";ASSEMBLY FAILED";
        let empty_sam = "@NO INFO\n";
        let empty_json = "{}";

        fs::write(out_dir.join("contigs.fasta"), dummy).await?;
        fs::write(out_dir.join("contigs.ram.fasta"), dummy).await?;
        fs::write(out_dir.join("contigs_all.fasta"), dummy).await?;
        fs::write(out_dir.join("scaffolds.fasta"), dummy).await?;
        fs::write(out_dir.join("read-contig.sam"), empty_sam).await?;
        fs::write(out_dir.join("contig_stats.json"), empty_json).await?;

        return Ok((
            CoverageOutputs {
                contigs_fasta: out_dir.join("contigs.fasta"),
                contigs_ram_fasta: out_dir.join("contigs.ram.fasta"),
                contigs_all_fasta: out_dir.join("contigs_all.fasta"),
                scaffolds_fasta: out_dir.join("scaffolds.fasta"),
                sam_path: out_dir.join("read-contig.sam"),
                contig_stats_json: out_dir.join("contig_stats.json"),
                coverage_json: out_dir.join("assembly_contig_coverage.json"),
                coverage_summary_csv: out_dir.join("assembly_contig_coverage_summary.csv"),
                contig_stats: Arc::new(HashMap::new()),
                read2contig: Arc::new(HashMap::new()),
            },
            ReceiverStream::new(empty_rx),
            cleanup_tasks,
            cleanup_receivers,
            temp_files,
        ));
    }

    let contigs_all_out = out_dir.join("contigs_all.fasta");
    fs::copy(&raw_contigs_path, &contigs_all_out).await?;

    let contigs_out = out_dir.join("contigs.fasta");
    let contigs_size = file_size(&raw_contigs_path).await?;

    let temp_dir = choose_temp_dir(
        contigs_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        assembly_headroom,
        false
    )
        .await?;

    let ram_fasta_path = NamedTempFile::with_suffix_in(".fa", &temp_dir)
        .map_err(|e| PipelineError::Other(e.into()))?;
    let ram_path = ram_fasta_path.path().to_owned();

    // Contig length filtering
    let rx = read_fasta(
        raw_contigs_path.clone(),
        u64::MAX,
        Some(config.args.min_contig_length),
        None,
        8192,
    )
        .map_err(|e| PipelineError::InvalidFastaFormat(e.to_string()))?;

    let write_handle = write_fasta_stream_to_file(
        ReceiverStream::new(rx),
        ram_path.clone(),
        config.base_buffer_size,
    );

    write_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Writer task panicked: {}", e)))?
        .map_err(|e| PipelineError::Other(anyhow!("Writer task failed: {}", e)))?;

    fs::copy(&ram_path, &contigs_out).await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to copy contigs to output: {}", e)))?;

    temp_files.push(ram_fasta_path);

    let scaffolds_out = out_dir.join("scaffolds.fasta");
    if raw_scaffolds_path.exists() {
        fs::copy(&raw_scaffolds_path, &scaffolds_out).await?;
    } else {
        fs::write(&scaffolds_out, ";NO SCAFFOLDS").await?;
    }

    let index_dir = out_dir.join("bowtie_index");
    fs::create_dir_all(&index_dir).await?;
    let index_prefix = index_dir.join("contigs");  //  will create contigs.1.bt2, contigs.2.bt2, ...

    let num_index_cores: usize = RunConfig::thread_allocation(&*config, BOWTIE2_TAG, None);

    let build_status = Command::new("bowtie2-build")
        .args([
            "--threads",
            &num_index_cores.to_string(),
            "--quiet",
            ram_path.clone().to_str().unwrap(),
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

    let (mut bt2_child, bt2_err_task) = spawn_cmd(
        config.clone(),
        BOWTIE2_TAG,
        bt2_args,
        config.args.verbose,
    ).await?;

    cleanup_tasks.push(bt2_err_task);

    let bt2_out_stream = parse_child_output(
        &mut bt2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
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

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
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
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let (non_host_streams, non_host_done_rx) = t_junction(
        ReceiverStream::new(samtools_sort_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "process_assembly_sam".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(non_host_done_rx);

    let mut non_host_streams_iter = non_host_streams.into_iter();
    let bam_for_file = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let bam_for_stats = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let bam_for_output = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let sam_path = out_dir.join("read-contig.bam");

    let write_sam_task = write_byte_stream_to_file(
        &sam_path,
        ReceiverStream::new(bam_for_file),
        Some(config.base_buffer_size),
    )
        .await?;
    cleanup_tasks.push(write_sam_task);

    let (read2contig, contig_stats, _contig_uniques) = generate_info_from_bam_stream(
        bam_for_stats,
        &duplicate_clusters,
        config.args.min_contig_length,
        &config.thread_pool
    )
        .await?;

    Ok((
        CoverageOutputs {
            contigs_fasta: out_dir.join("contigs.fasta"),
            contigs_ram_fasta: ram_path.clone(),
            contigs_all_fasta: out_dir.join("contigs_all.fasta"),
            scaffolds_fasta: out_dir.join("scaffolds.fasta"),
            sam_path: out_dir.join("read-contig.sam"),
            contig_stats_json: out_dir.join("contig_stats.json"),
            coverage_json: out_dir.join("assembly_contig_coverage.json"),
            coverage_summary_csv: out_dir.join("assembly_contig_coverage_summary.csv"),
            contig_stats: Arc::new(contig_stats),
            read2contig: Arc::new(read2contig),
        },
        ReceiverStream::new(bam_for_output),
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
        ParseOutput::Bytes(b) => Ok(Bytes::from(b.to_vec())),
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
    let estimated_size = (acc_vec.len() as u64) * 2000; // ~2KB per accession avg

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
    mut stream: ReceiverStream<ParseOutput>,
    min_reads_per_genus: usize,
) -> Result<(
    AHashMap<String, Arc<ReadHit>>,
    AHashMap<String, AccessionHit>,
    HashMap<i32, Vec<String>>,
    usize,
)> {
    let read_dict: DashMap<String, Arc<ReadHit>> = DashMap::with_capacity(8_000_000);
    let accession_dict: DashMap<String, AccessionHit> = DashMap::with_capacity(2_000_000);

    let genus_read_counts: DashMap<i32, usize> = DashMap::new();
    let genus_species: DashMap<i32, HashSet<i32>> = DashMap::new();
    let genus_accessions: DashMap<i32, HashSet<String>> = DashMap::new();

    let mut total_reads: usize = 0;

    while let Some(item) = stream.next().await {
        let bytes = match item {
            ParseOutput::Bytes(b) => b,
            _ => return Err(anyhow!("Expected Bytes in summarize_hits stream")),
        };

        let line = match std::str::from_utf8(&bytes) {
            Ok(s) => s.trim_end(),
            Err(_) => {
                warn!("Non-UTF8 line in hit summary stream");
                continue;
            }
        };

        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() != 6 {
            warn!(
        "Malformed hit_summary line (expected exactly 6 fields, got {}): {}",
        parts.len(),
        line
    );
            continue;
        }

        let read_id       = parts[0].to_string();
        let accession_id  = if parts[1].is_empty() || parts[1] == "-" {
            "-".to_string()
        } else {
            parts[1].to_string()
        };
        let hit_taxid: i32 = parts[2].parse()
            .map_err(|_| warn!("Invalid hit_taxid: {}", parts[2])).ok().unwrap_or(0);
        let genus_taxid: i32 = parts[3].parse()
            .map_err(|_| warn!("Invalid genus_taxid: {}", parts[3])).ok().unwrap_or(0);
        let family_taxid: i32 = parts[4].parse()
            .map_err(|_| warn!("Invalid family_taxid: {}", parts[4])).ok().unwrap_or(0);
        let level: u8 = parts[5].parse()
            .map_err(|_| warn!("Invalid level: {}", parts[5])).ok().unwrap_or(0);

        let species_taxid = if level >= 3 && hit_taxid > 0 {
            hit_taxid                     // level 3+ means species consensus already achieved
        } else {
            0                             // only genus/family → species unknown
        };

        // assign highest resolved level
        let assigned_taxid = if level >= 3 {
            species_taxid
        } else if level == 2 {
            genus_taxid
        } else if level == 1 {
            genus_taxid
        } else {
            hit_taxid                     // fallback (should be rare)
        };

        read_dict.insert(
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
                source_count_type: None
            }),
        );

        total_reads += 1;

        if accession_id != "-" && genus_taxid > 0 {
            accession_dict
                .entry(accession_id.clone())
                .or_insert_with(|| AccessionHit {
                    taxid: species_taxid,
                    genus_taxid,
                    family_taxid,
                    count: 0,
                })
                .count += 1;

            *genus_read_counts.entry(genus_taxid).or_insert(0) += 1;
            genus_species.entry(genus_taxid).or_default().insert(species_taxid);
            genus_accessions.entry(genus_taxid).or_default().insert(accession_id);
        }
    }

    info!(
        "summarize_hits: processed {} reads, {} unique accessions",
        total_reads,
        accession_dict.len()
    );

    let mut selected_genera: HashMap<i32, Vec<String>> = HashMap::new();
    for entry in genus_read_counts.iter() {
        let genus_taxid = *entry.key();
        let read_count = *entry.value();

        if read_count < min_reads_per_genus {
            continue;
        }

        let species_set = genus_species.get(&genus_taxid);
        let species_count = species_set.as_ref().map(|s| s.len()).unwrap_or(0);

        if species_count <= 1 {
            continue;
        }

        if let Some(acc_set) = genus_accessions.get(&genus_taxid) {
            let accessions: Vec<String> = acc_set.iter().cloned().collect();
            selected_genera.insert(genus_taxid, accessions);
        }
    }

    let read_dict_final: AHashMap<String, Arc<ReadHit>> = read_dict.into_iter().collect();
    let accession_dict_final: AHashMap<String, AccessionHit> = accession_dict.into_iter().collect();

    Ok((
        read_dict_final,
        accession_dict_final,
        selected_genera,
        total_reads,
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
    while let Some(item) = top_m8_stream.next().await {
        let bytes = match item {
            ParseOutput::Bytes(b) => b,
            _ => continue,
        };

        let line = String::from_utf8_lossy(&**bytes);
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
        let accession = m8.tname.clone();

        let taxid = accession_map
            .get(&accession)
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

        let json = json!({
            "contig_name": contig_id,
            "species_taxid": lineage[0],
            "genus_taxid": lineage[1],
            "family_taxid": lineage[2],
        });
        let mut json_line = serde_json::to_string(&json)?;
        json_line.push('\n');
        contig2lineage_tx
            .send(ParseOutput::Bytes(Arc::new(json_line.into_bytes())))
            .await?;

        // Send refined M8 line
        let mut m8_line = m8.to_tab_string();
        m8_line.push('\n');
        read2blastm8_tx
            .send(ParseOutput::Bytes(Arc::new(m8_line.into_bytes())))
            .await?;

        // Get all reads in this contig
        let reads_in_contig: Vec<String> = read2contig
            .iter()
            .filter(|&(_, c)| c == &contig_id)
            .map(|(r, _)| r.clone())
            .collect();

        let shared_hit: Arc<ReadHit> = Arc::new(ReadHit {
            level: 1,
            taxid: lineage[0],
            accession_id: accession.clone(),
            species_taxid: lineage[0],
            genus_taxid: lineage[1],
            family_taxid: lineage[2],
            contig_id: Some(contig_id.clone()),
            contig_accession_id: Some(accession.clone()),
            contig_species_taxid: lineage[0],
            contig_genus_taxid: lineage[1],
            contig_family_taxid: lineage[2],
            from_assembly: true,
            source_count_type: Some(db_type.to_uppercase()),
        });


        for read_id in reads_in_contig {
            let final_line = shared_hit.to_full_tab_line(&read_id);

            let was_present = {
                let mut dict = read_dict.lock().unwrap();
                let was_present = dict.contains_key(&read_id);
                dict.insert(read_id.clone(), shared_hit.clone());
                was_present
            };

            let bytes = Arc::new(final_line.into_bytes());
            if was_present {
                let _ = updated_tx.send(ParseOutput::Bytes(bytes)).await;
            } else {
                let _ = added_tx.send(ParseOutput::Bytes(bytes)).await;
            }
        }

    }

    Ok(())
}



pub async fn generate_contig_summary_json(
    read2contig: Arc<HashMap<String, String>>,
    contig2lineage: AHashMap<String, [i32; 3]>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    db_type: &str,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    min_contig_size: u64,
    mut output_tx: Sender<ParseOutput>,
) -> Result<()> {
    let mut genus_summary: AHashMap<i32, AHashMap<String, [u64; 2]>> = AHashMap::new();
    let mut species_summary: AHashMap<i32, AHashMap<String, [u64; 2]>> = AHashMap::new();

    //lock once
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
                .map(|cluster| cluster.size)
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
                    .send(ParseOutput::Bytes(Arc::new(line.into_bytes())))
                    .await
                    .map_err(|_| anyhow!("contig_summary_tx dropped"))?;
            }
        }
    }

    Ok(())
}


pub async fn generate_m8_and_hit_summary(
    mut updated_reads_stream: ReceiverStream<ParseOutput>,
    mut added_reads_stream: ReceiverStream<ParseOutput>,
    mut blast_hits_stream: ReceiverStream<ParseOutput>,
    mut original_hit_summary_stream: ReceiverStream<ParseOutput>,
    mut original_deduped_m8_stream: ReceiverStream<ParseOutput>,
    refined_m8_tx: Sender<ParseOutput>,
    refined_hit_summary_tx: Sender<ParseOutput>,
) -> Result<()> {


    // copy everything from original deduped M8
    while let Some(item) = original_deduped_m8_stream.next().await {
        refined_m8_tx.send(item).await
            .map_err(|_| anyhow!("refined_m8_tx dropped"))?;
    }

    // append new blast hits
    while let Some(item) = blast_hits_stream.next().await {
        refined_m8_tx.send(item).await
            .map_err(|_| anyhow!("refined_m8_tx dropped"))?;
    }
    drop(refined_m8_tx); // signals EOF


    // stream through the original hit summary and replace rows that were updated
    let mut updated_reads = AHashMap::new();
    while let Some(item) = updated_reads_stream.next().await {
        if let ParseOutput::Bytes(bytes) = item {
            let line = String::from_utf8_lossy(&*bytes);
            let read_id = line.split('\t').next().unwrap_or("").to_string();
            updated_reads.insert(read_id, bytes);
        }
    }

    let mut added_reads = AHashMap::new();
    while let Some(item) = added_reads_stream.next().await {
        if let ParseOutput::Bytes(bytes) = item {
            let line = String::from_utf8_lossy(&*bytes);
            let read_id = line.split('\t').next().unwrap_or("").to_string();
            added_reads.insert(read_id, bytes);
        }
    }

    //walk the original hit summary once
    while let Some(item) = original_hit_summary_stream.next().await {
        if let ParseOutput::Bytes(bytes) = item {
            let line = String::from_utf8_lossy(&*bytes);
            let read_id = line.split('\t').next().unwrap_or("").to_string();

            // pefer updated version, then added, otherwise keep original
            let to_send = updated_reads.get(&read_id)
                .or_else(|| added_reads.get(&read_id))
                .map_or(bytes.clone(), |v| v.clone());

            refined_hit_summary_tx.send(ParseOutput::Bytes(to_send)).await
                .map_err(|_| anyhow!("refined_hit_summary_tx dropped"))?;
        }
    }

    // append  added reads not in the original hit summary
    for (read_id, bytes) in added_reads {
        if !updated_reads.contains_key(&read_id) {
            refined_hit_summary_tx.send(ParseOutput::Bytes(bytes)).await
                .map_err(|_| anyhow!("refined_hit_summary_tx dropped"))?;
        }
    }

    drop(refined_hit_summary_tx);
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
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    blast_headroom: u64,
) -> Result<(
    AHashMap<String, Arc<ReadHit>>,
    Vec<TaxonCount>,
    Vec<ContigSummaryEntry>,
    mpsc::Receiver<ParseOutput>,   // refined reassigned M8
    mpsc::Receiver<ParseOutput>,   // refined hit summary
    mpsc::Receiver<ParseOutput>,   // top hit per contig M8 (blast_top_m8)
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
    Vec<NamedTempFile>,
)> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();

    let out_dir = config.out_dir.join(format!("blast_{}", db_type));
    tokio::fs::create_dir_all(&out_dir).await?;

    let blast_m8_path = out_dir.join("blast.m8");
    let blast_top_m8_path = out_dir.join("blast_top.m8");
    let refined_m8_path = out_dir.join("refined.m8");
    let refined_hit_summary_path = out_dir.join("refined_hit_summary.tab");
    let refined_counts_path = out_dir.join("refined_counts_with_dcr.json");
    let contig_summary_path = out_dir.join("contig_summary.json");

    let contig_size = file_size(assembled_contig_fasta).await?;
    let ref_size = file_size(reference_fasta).await?;
    if contig_size < MIN_REF_FASTA_SIZE || ref_size < MIN_ASSEMBLED_CONTIG_SIZE {
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
        true
    )
        .await?;

    let blastdb_suffix = format!("{}_blastindex", db_type);
    let blastdb_ram_path = NamedTempFile::with_suffix_in(blastdb_suffix, &temp_dir)
        .map_err(|e| PipelineError::Other(e.into()))?;
    let blastdb_path = blastdb_ram_path.path().to_owned();
    temp_files.push(blastdb_ram_path);

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
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAKEBLASTDB_TAG.to_string(),
            error: e.to_string(),
        })?;
    debug!("Spawning makeblastdb with args: {:?}", makeblastdb_args);

    let (_, makeblastdb_err_task) = spawn_cmd(
        config.clone(),
        MAKEBLASTDB_TAG,
        makeblastdb_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAKEBLASTDB_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(makeblastdb_err_task);

    let blast_command = if db_type == NT_TAG { BLASTN_TAG } else { BLASTX_TAG };

    let blast_args = if db_type == NT_TAG {
        let blastn_config = BlastnConfig {
            query: assembled_contig_fasta.clone(),
            db: blastdb_path.clone(),
            outfmt: "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
                .to_string(),
            evalue: 1e-10,
            max_target_seqs: 5000,
            option_fields: HashMap::new(),
        };

        generate_cli(BLASTN_TAG, &config, Some(&blastn_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: BLASTN_TAG.to_string(),
                error: e.to_string(),
            }
        })?
    } else {
        let blastx_config = BlastxConfig {
            query: assembled_contig_fasta.clone(),
            db: blastdb_path.clone(),
            outfmt: 6,
            evalue: 1e-10,
            num_alignments: 5,
            option_fields: HashMap::new(),
        };

        generate_cli(BLASTX_TAG, &config, Some(&blastx_config)).map_err(|e| {
            PipelineError::ToolExecution {
                tool: BLASTX_TAG.to_string(),
                error: e.to_string(),
            }
        })?
    };

    let (mut blast_child, err_task) = spawn_cmd(
        config.clone(),
        blast_command,
        blast_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: blast_command.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(err_task);

    let blast_out_stream = parse_child_output(
        &mut blast_child,
        ChildStream::Stdout,
        ParseMode::Lines,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "blastn/x".to_string(),
            error: e.to_string(),
        })?;

    // Internal channel: get_top_m8_* writes here
    let (top_internal_tx, top_internal_rx) = mpsc::channel(1024);

    // Public channels: one for update_read_dict, one for caller (coverage viz)
    let (top_for_update_tx, top_for_update_rx) = mpsc::channel(1024);
    let (top_for_caller_tx, top_for_caller_rx) = mpsc::channel(1024);

    let top_handle = if db_type == NT_TAG {
        tokio::spawn(get_top_m8_nt(ReceiverStream::new(blast_out_stream), top_internal_tx))
    } else {
        tokio::spawn(get_top_m8_nr(ReceiverStream::new(blast_out_stream), top_internal_tx))
    };

    // Forward top hits to both consumers (update_read_dict and caller)
    let forward_handle = tokio::spawn({
        let mut rx = top_internal_rx;
        let tx_update = top_for_update_tx;
        let tx_caller = top_for_caller_tx;
        async move {
            while let Some(item) = rx.recv().await {
                let _ = tx_update.send(item.clone()).await;
                let _ = tx_caller.send(item).await;
            }
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(forward_handle);

    let (contig2lineage_tx, contig2lineage_rx) = mpsc::channel(1024);
    let (read2blastm8_tx, read2blastm8_rx) = mpsc::channel(1024);
    let (updated_tx, updated_rx) = mpsc::channel(1024);
    let (added_tx, added_rx) = mpsc::channel(1024);

    let update_handle = tokio::spawn(update_read_dict(
        read2contig.clone(),
        ReceiverStream::new(top_for_update_rx),  // top hits only — matches Python
        read_dict.clone(),
        lineage_map.clone(),
        accession_map.clone(),
        should_keep_filter.clone(),
        db_type,
        contig2lineage_tx,
        read2blastm8_tx,
        updated_tx,
        added_tx,
    ));

    // Local mpsc channels used by generate_m8_and_hit_summary
    let (refined_m8_local_tx, refined_m8_local_rx) = mpsc::channel(2_000_000);
    let (refined_hit_summary_local_tx, refined_hit_summary_local_rx) = mpsc::channel(2_000_000);

    let generate_m8_handle = tokio::spawn(generate_m8_and_hit_summary(
        ReceiverStream::new(updated_rx),
        ReceiverStream::new(added_rx),
        ReceiverStream::new(read2blastm8_rx),
        hit_summary_stream,
        deduped_m8_stream,
        refined_m8_local_tx.clone(),
        refined_hit_summary_local_tx.clone(),
    ));

    // Output mpsc channels (returned to caller)
    let (refined_m8_tx, refined_m8_rx) = mpsc::channel(2_000_000);
    let (refined_hit_summary_tx, refined_hit_summary_rx) = mpsc::channel(2_000_000);

    // Internal channels for taxon counting
    let (refined_m8_for_counts_tx, refined_m8_for_counts_rx) = mpsc::channel(2_000_000);
    let (refined_hit_for_counts_tx, refined_hit_for_counts_rx) = mpsc::channel(2_000_000);

    // Forwarding: local → output + counts
    let m8_forward_handle = tokio::spawn({
        let mut rx = refined_m8_local_rx;
        let out_tx = refined_m8_tx.clone();
        let counts_tx = refined_m8_for_counts_tx.clone();
        async move {
            while let Some(item) = rx.recv().await {
                let _ = out_tx.send(item.clone()).await;
                let _ = counts_tx.send(item).await;
            }
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(m8_forward_handle);

    let hit_forward_handle = tokio::spawn({
        let mut rx = refined_hit_summary_local_rx;
        let out_tx = refined_hit_summary_tx.clone();
        let counts_tx = refined_hit_for_counts_tx.clone();
        async move {
            while let Some(item) = rx.recv().await {
                let _ = out_tx.send(item.clone()).await;
                let _ = counts_tx.send(item).await;
            }
            Ok::<(), anyhow::Error>(())
        }
    });
    cleanup_tasks.push(hit_forward_handle);

    // Internal taxon count aggregation
    let (refined_counts_tx, refined_counts_rx) = mpsc::channel(1024);
    let counts_handle = tokio::spawn(generate_taxon_count_json_from_m8(
        ReceiverStream::new(refined_m8_for_counts_rx),
        ReceiverStream::new(refined_hit_for_counts_rx),
        db_type,
        lineage_map.clone(),
        should_keep_filter.clone(),
        duplicate_clusters.clone(),
        refined_counts_tx,
    ));

    update_handle.await??;

    let mut contig2lineage: AHashMap<String, [i32; 3]> = AHashMap::new();
    let mut lineage_stream = ReceiverStream::new(contig2lineage_rx);
    while let Some(item) = lineage_stream.next().await {
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
    }

    let (contig_summary_tx, contig_summary_rx) = mpsc::channel(1024);
    let contig_summary_handle = tokio::spawn(generate_contig_summary_json(
        read2contig.clone(),
        contig2lineage,
        read_dict.clone(),
        db_type,
        duplicate_clusters.clone(),
        4,
        contig_summary_tx,
    ));

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

    let mut refined_counts = Vec::new();
    let mut rx = ReceiverStream::new(refined_counts_rx);
    while let Some(item) = rx.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let count: TaxonCount = serde_json::from_str(&line)?;
        refined_counts.push(count);
    }

    let mut contig_summary = Vec::new();
    let mut rx = ReceiverStream::new(contig_summary_rx);
    while let Some(item) = rx.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let entry: ContigSummaryEntry = serde_json::from_str(&line)?;
        contig_summary.push(entry);
    }

    let final_read_dict = read_dict.lock().unwrap().clone();

    Ok((
        final_read_dict,
        refined_counts,
        contig_summary,
        refined_m8_rx,
        refined_hit_summary_rx,
        top_for_caller_rx,  // blast_top_m8 — ready for generate_coverage_viz
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
        None,
        u64::MAX,
        None,
        None,
        config.base_buffer_size * 4, // Large chunks for speed
    )?;
    let r1_input_stream = ReceiverStream::new(r1_rx);
    let r1_filtered = filter_fastq_to_bytes_stream(r1_input_stream, r1_headers).await;
    let write_r1 = write_byte_stream_to_file(&out_r1, r1_filtered, Some(config.base_buffer_size)).await?;

    // Stream and filter R2 if present
    let write_r2 = if let (Some(r2_path), Some(r2_headers)) = (original_r2_path, r2_headers) {
        let (r2_rx, _r2_stats_task) = read_fastq(
            r2_path.clone(),
            None,
            None,
            u64::MAX,
            None,
            None,
            config.base_buffer_size * 4,
        )?;
        let r2_input_stream = ReceiverStream::new(r2_rx);
        let r2_filtered = filter_fastq_to_bytes_stream(r2_input_stream, r2_headers).await;
        Some(write_byte_stream_to_file(&out_r2.clone().unwrap(), r2_filtered, Some(config.base_buffer_size)).await?)
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
    let mut temp_paths: Vec<PathBuf> = Vec::new();

    info!("Starting short read mNGS pipeline.");

    // *******************
    // Setup and Validation
    // *******************

    // External tools check
    check_versions(vec![BOWTIE2_TAG, MINIMAP2_TAG, KALLISTO_TAG, SPADES_TAG, MAKEBLASTDB_TAG,
                        BLASTN_TAG, BLASTX_TAG, DIAMOND_TAG], &out_dir.clone())
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Check required files
    let host_bowtie2_index: String = config.args.host_bowtie2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_bowtie2_index is required".to_string()))?;

    let host_hisat2_index: String = config.args.host_hisat2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_hisat2_index is required".to_string()))?;

    let (file1_path, file2_path, sample_base_buf, sample_base) = validate_file_inputs(&config, &cwd).await?;
    let paired = file2_path.is_some();

    let seed = config.args.seed.unwrap_or_else(|| {
        let mut bytes = [0u8; 8];
        OsRng.fill_bytes(&mut bytes);
        let random_seed = u64::from_le_bytes(bytes);
        random_seed
    });

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
    let (ercc_bt2_out_stream, ercc_count_rx, ercc_bt2_cleanup_tasks, ercc_bt2_cleanup_receivers, ercc_bt2_bam_write_handle, ercc_bt2_bam_path, ercc_bt2_temp_dirs) = bowtie2_filter(
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
        ercc_bt2_out_stream,
    ).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);
    cleanup_receivers.extend(qc_cleanup_receivers);

    // Split for Kallisto and Bowtie2
    let (mut kallisto_streams, kallisto_split_done_rx) = t_junction(
        qc_fastp_out_stream,
        2,
        config.base_buffer_size * 10,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::JustBytes,
        "kallisto_split".to_string(),
        None,
    ).await?;
    cleanup_receivers.push(kallisto_split_done_rx);

    let mut kallisto_streams_iter = kallisto_streams.into_iter();
    let kallisto_bypass_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let kallisto_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let kallisto_stream = ReceiverStream::new(kallisto_stream);

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
    let (host_bt2_out_stream, host_bt2_count_rx, host_bt2_cleanup_tasks, host_bt2_cleanup_receivers, host_bt2_bam_write_handle, host_bt2_bam_path, host_bt2_temp_dirs) = bowtie2_filter(
        config.clone(),
        ReceiverStream::new(kallisto_bypass_stream),
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
        host_bt2_out_stream,
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
            human_bt2_out_stream,
            _human_bt2_count_rx,
            human_bt2_cleanup_tasks,
            human_bt2_cleanup_receivers,
            human_bt2_bam_write_handle,
            human_bt2_bam_path,
            human_bt2_temp_dirs
        ) = bowtie2_filter(
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
            human_bt2_out_stream,
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

    let (pre_dedup_parsed_stream, parse_task) = parse_byte_stream_to_fastq(
        post_filter_stream.into_inner(),
        config.base_buffer_size,
        config.args.stall_threshold,
    ).await?;

    cleanup_tasks.push(parse_task);

    let (dedup_stream, dedup_count_rx, cluster_stream, duplicate_clusters, mut dedup_cleanup_tasks, mut dedup_cleanup_receivers) = dedup(
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
    let non_host_deinterleave_task = tokio::spawn(deinterleave_fastq_stream(
        subsampled_stream,
        paired,
        config.base_buffer_size * 4,
    ));

    let (r1_rx, r2_rx_opt, deinterleave_handle) = non_host_deinterleave_task.await??;

    let non_host_temp_dir = choose_temp_dir(
        config.input_size * 2,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        true,
    ).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let non_host_r1_path = non_host_temp_dir.path().join("nonhost_R1.fastq");
    let non_host_r2_path_opt = r2_rx_opt.as_ref().map(|_| non_host_temp_dir.path().join("nonhost_R2.fastq"));

    info!("Checkpoint: deinterleaving and writing non-host R1/R2 to temp (forces upstream completion)");

    let r1_write_task = write_parse_output_to_file(
        &non_host_r1_path,
        ReceiverStream::new(r1_rx),
        Some(config.base_buffer_size * 4),
    ).await?;

    r1_write_task.await??;

    let r2_write_task_opt = r2_rx_opt.map(|r2_rx| {
        let r2_path = non_host_r2_path_opt.as_ref().expect("R2 path should exist when r2_rx is Some");
        write_parse_output_to_file(
            r2_path,
            ReceiverStream::new(r2_rx),
            Some(config.base_buffer_size * 4),
        )
    });

    if let Some(task_future) = r2_write_task_opt {
        let task_handle = task_future.await?;
        task_handle.await?;
    }

    deinterleave_handle.await??;

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

    // Now spades/mm2/diamond/etc. from files (use original paths — they weren't moved)
    let spades_task = spades_assembly(
        config.clone(),
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        &out_dir,
    ).await?;

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

    let (lineage_map, acc2taxid_map) = taxonomy_handle.await??;

    let nt_concurrency = compute_phase_concurrency(
        &config,
        "call_hits_nt",
        1.0,           // ~1 GB per thread max (mostly transient)
        3.5,
        64,            // higher cap for CPU-bound phase
        16,            // min for meaningful parallelism
    );
    info!("call hits nt concurrency {}", nt_concurrency);

    let (call_stream, call_summary_stream, mut call_cleanup_tasks, mut call_cleanup_receivers) = call_hits_m8_stream(
        config.clone(),
        m8_stream,
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        36,
        nt_concurrency,
    ).await?;
    cleanup_tasks.append(&mut call_cleanup_tasks);
    cleanup_receivers.append(&mut call_cleanup_receivers);

    let (nt_streams, nt_done_rx) = t_junction(
        call_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nt_call".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nt_done_rx);

    let mut nt_streams_iter = nt_streams.into_iter();
    let nt_call_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_m8_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_blast_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (nt_summary_streams, nt_summary_done_rx) = t_junction(
        call_summary_stream,
        4,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nt_call_summary".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nt_summary_done_rx);
    let mut nt_summary_streams_iter = nt_summary_streams.into_iter();
    let nt_summary_taxon_stream = nt_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_summary_hit_stream = nt_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_initial_stream = nt_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_blast_hit_stream = nt_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let nt_hit_summary_handle = tokio::spawn(summarize_hits(ReceiverStream::new(nt_summary_hit_stream), 0));

    let nt_counts = generate_taxon_counts(
        config.clone(),
        ReceiverStream::new(nt_call_stream),
        ReceiverStream::new(nt_summary_taxon_stream),
        duplicate_clusters.clone(),
        should_keep_filter.clone(),
        "NT".to_string(),
        None,
    ).await?;
    info!("NT taxon counts: {} entries", nt_counts.len());

    let nt_map = collect_m8_to_accession_map(ReceiverStream::new(nt_m8_stream)).await?;

    // Temporary: Skip Diamond by providing empty outputs
    let (dummy_tx, dummy_rx) = mpsc::channel::<ParseOutput>(1);
    drop(dummy_tx); // Immediately drop sender to create an empty stream
    let non_host_diamond_m8_stream = dummy_rx;
    let mut non_host_diamond_cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut non_host_diamond_cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();


    // Diamond non_host alignment
    // let (non_host_diamond_m8_stream, mut non_host_diamond_cleanup_tasks, mut non_host_diamond_cleanup_receivers, non_host_diamond_temp_dirs) = diamond_non_host_align(
    //     config.clone(),
    //     ReceiverStream::new(non_host_dmnd_stream),
    //     paired,
    //     sample_base.clone()
    // ).await?;
    // cleanup_tasks.append(&mut non_host_diamond_cleanup_tasks);
    // cleanup_receivers.append(&mut non_host_diamond_cleanup_receivers);
    // temp_dirs.extend(non_host_diamond_temp_dirs);

    let nr_concurrency = compute_phase_concurrency(
        &config,
        "call_hits_nr",
        1.0,           // ~1 GB per thread max (mostly transient)
        3.5,
        64,            // higher cap for CPU-bound phase
        16,            // min for meaningful parallelism
    );

    info!("call hits nr concurrency {}", nr_concurrency);


    let (nr_call_stream, nr_call_summary_stream, mut nr_call_cleanup_tasks, mut nr_call_cleanup_receivers) = call_hits_m8_stream(
        config.clone(),
        ReceiverStream::new(non_host_diamond_m8_stream),
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        0,
        nr_concurrency,
    ).await?;
    cleanup_tasks.append(&mut nr_call_cleanup_tasks);
    cleanup_receivers.append(&mut nr_call_cleanup_receivers);

    let (nr_streams, nr_done_rx) = t_junction(
        nr_call_stream,
        4,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nr_call".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nr_done_rx);

    let mut nr_streams_iter = nr_streams.into_iter();
    let nr_call_stream = nr_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_m8_stream = nr_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_initial_stream = nr_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_blast_stream = nr_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let (nr_summary_streams, nr_summary_done_rx) = t_junction(
        nr_call_summary_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nr_call_summary".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nr_summary_done_rx);
    let mut nr_summary_streams_iter = nr_summary_streams.into_iter();
    let nr_summary_taxon_stream = nr_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_summary_hit_stream = nr_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_blast_hit_stream = nr_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let nr_hit_summary_handle = tokio::spawn(summarize_hits(ReceiverStream::new(nr_summary_hit_stream), 0));

    let nr_counts = generate_taxon_counts(
        config.clone(),
        ReceiverStream::new(nr_call_stream),
        ReceiverStream::new(nr_summary_taxon_stream),
        duplicate_clusters.clone(),
        should_keep_filter.clone(),
        "NR".to_string(),
        None,
    ).await?;
    info!("NR taxon counts: {} entries", nr_counts.len());

    let nr_map = collect_m8_to_accession_map(ReceiverStream::new(nr_m8_stream)).await?;

    let combined_path = out_dir.join(rename_file_path(&sample_base_buf, None, Some("taxon_counts_with_dcr.json"), "_"));
    let (_combined_path, write_json_task) = combine_taxon_counts(
        &nt_counts,
        &nr_counts,
        combined_path,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts failed: {}", e)))?;
    cleanup_tasks.push(write_json_task);


    let (annotated_rx, unidentified_rx, unique_unidentified_rx, mut annot_tasks, mut annot_rxs) = generate_annotated_fasta_from_files(
        config.clone(),
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        cluster_stream,
        nt_map,
        nr_map,
    ).await
        .map_err(|e| PipelineError::Other(anyhow!("Annotated FASTA from files failed: {}", e)))?;

    cleanup_tasks.append(&mut annot_tasks);
    cleanup_receivers.append(&mut annot_rxs);

    let assembly_dir = out_dir.join("assembly");

    let annotated_path = assembly_dir.join("annotated_merged.fa");
    let unidentified_path = assembly_dir.join("unidentified.fa");
    let unique_unidentified_path = assembly_dir.join("unique_unidentified.fa");


    let (initial_annotated_streams, initial_annotated_done_rx) = t_junction(
        ReceiverStream::new(annotated_rx),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nr_call_summary".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(initial_annotated_done_rx);
    let mut initial_annotated_streams_iter = initial_annotated_streams.into_iter();
    let initial_annotated_file_stream = initial_annotated_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_annotated_taxon_stream = initial_annotated_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let (initial_unidentified_streams, initial_unidentified_done_rx) = t_junction(
        ReceiverStream::new(unidentified_rx),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nr_call_summary".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(initial_unidentified_done_rx);
    let mut initial_unidentified_streams_iter = initial_unidentified_streams.into_iter();
    let initial_unidentified_file_stream = initial_unidentified_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_unidentified_taxon_stream = initial_unidentified_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(initial_annotated_file_stream),
        annotated_path.clone(),
        config.base_buffer_size,
    ));

    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(initial_unidentified_file_stream),
        unidentified_path.clone(),
        config.base_buffer_size,
    ));

    cleanup_tasks.push(write_fasta_stream_to_file(
        ReceiverStream::new(unique_unidentified_rx),
        unique_unidentified_path.clone(),
        config.base_buffer_size,
    ));


    // *******************
    // Post-processing
    // *******************

    // Assembly stats
    let assembly_out_dir = out_dir.join("assembly");  // Assuming assembly_handle.out_dir was this
    let assembly_work_dir = assembly_out_dir.join("spades");  // Match Python's subdir

    let spades_completion = spades_task.await;

    match spades_completion {
        Ok(Ok(())) => {
            info!("SPAdes assembly task completed successfully");
        }
        Ok(Err(e)) => {
            warn!("SPAdes assembly failed: {}", e);
        }
        Err(join_err) => {
            error!("SPAdes task panicked or was cancelled: {}", join_err);
        }
    }

    let (assembly_outputs, assembly_bam_out_stream,
        mut post_assembly_cleanup_tasks,
        mut post_assembly_cleanup_receivers,
        mut post_assembly_temp_files) = process_assembly(
        config.clone(),
        &assembly_out_dir,
        &assembly_work_dir,
        non_host_r1_path.clone(),
        non_host_r2_path_opt.clone(),
        duplicate_clusters.clone(),
        paired,
        4 // this can become as CLI arg if needed
    ).await?;

    cleanup_tasks.extend(post_assembly_cleanup_tasks);
    cleanup_receivers.extend(post_assembly_cleanup_receivers);
    temp_files.extend(post_assembly_temp_files);

    let _ = fs::remove_dir_all(assembly_work_dir).await;

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
    //NT
    let (
        nt_read_dict_noarc,
        nt_accession_dict_noarc,
        nt_selected_genera,
        nt_total_reads,
    ) = nt_hit_summary_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("NT hit summary task panicked: {}", e)))?
        .map_err(|e| PipelineError::Other(anyhow!("NT hit summary parsing failed: {}", e)))?;

    let nt_read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>> = Arc::new(Mutex::new(nt_read_dict_noarc));
    let nt_accession_dict = Arc::new(nt_accession_dict_noarc);

    //nt accessions
    let (nt_ref_fasta_path, nt_ref_fasta_temp_dir) = build_reference_fasta_from_selected_genera(
        config.clone(),
        &nt_selected_genera,
        NT_TAG,
        &PathBuf::from(nt_file),
        &PathBuf::from(nt_offset_db_file),
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    eprintln!("nt path {}", nt_ref_fasta_path.display());
    final_temp_dirs.push(nt_ref_fasta_temp_dir);

    // //NR
    let (
        nr_read_dict_noarc,
        nr_accession_dict_noarc,
        nr_selected_genera,
        nr_total_reads,
    ) = nr_hit_summary_handle
        .await
        .map_err(|e| PipelineError::Other(anyhow!("NT hit summary task panicked: {}", e)))?
        .map_err(|e| PipelineError::Other(anyhow!("NT hit summary parsing failed: {}", e)))?;

    let nr_read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>> = Arc::new(Mutex::new(nr_read_dict_noarc));
    let nr_accession_dict = Arc::new(nr_accession_dict_noarc);


    let (nr_ref_fasta_path, nr_ref_fasta_temp_dir) = build_reference_fasta_from_selected_genera(
        config.clone(),
        &nr_selected_genera,
        NR_TAG,
        &PathBuf::from(nr_file),
        &PathBuf::from(nr_offset_db_file),
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    eprintln!("nr path {}", nr_ref_fasta_path.display());
    final_temp_dirs.push(nr_ref_fasta_temp_dir);

    let nt_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_ram_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let read2contig = assembly_outputs.read2contig.clone();

        async move {
            blast_contigs(
                config,
                NT_TAG,
                ReceiverStream::new(nt_blast_stream),
                ReceiverStream::new(nt_blast_hit_stream),
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
            ).await
        }
    });

    let nr_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_ram_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_clusters = duplicate_clusters.clone();
        let read2contig = assembly_outputs.read2contig.clone();

        async move {
            blast_contigs(
                config,
                NR_TAG,
                ReceiverStream::new(nr_blast_stream),
                ReceiverStream::new(nr_blast_hit_stream),
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
            ).await
        }
    });


    let (nt_res, nr_res) = tokio::try_join!(nt_handle, nr_handle)
        .map_err(|e| PipelineError::Other(anyhow!("blast_contigs task panicked: {e}")))?;

    let (nt_read_dict, nt_refined_counts,
        nt_contig_summary, nt_refined_m8_stream_out,
        nt_refined_hit_summary_stream_out, nt_refined_m8_top_stream_out,
        nt_cleanup_tasks, nt_cleanup_receivers,
        nt_temp_files) = nt_res?;

    cleanup_tasks.extend(nt_cleanup_tasks);
    cleanup_receivers.extend(nt_cleanup_receivers);
    temp_files.extend(nt_temp_files);


    let (nt_m8_streams, nt_m8_rx) = t_junction(
        ReceiverStream::new(nt_refined_m8_stream_out),
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nt_m8".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nt_m8_rx);

    let mut nt_m8_streams_iter = nt_m8_streams.into_iter();
    let nt_m8_merge = nt_m8_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_m8_map = nt_m8_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_m8_viz = nt_m8_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (nt_hitsummary_streams, nt_hitsummary_rx) = t_junction(
        ReceiverStream::new(nt_refined_hit_summary_stream_out),
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nt_hitsummary_rx);

    let mut nt_hitsummary_streams_iter = nt_hitsummary_streams.into_iter();
    let nt_hit_summary_merge = nt_hitsummary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_hit_summary_taxid = nt_hitsummary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_hit_summary_coverage = nt_hitsummary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let (nr_read_dict, nr_refined_counts,
        nr_contig_summary, nr_refined_m8_stream_out,
        nr_refined_hit_summary_stream_out, nr_refined_m8_top_stream_out,
        nr_cleanup_tasks, nr_cleanup_receivers,
        nr_temp_files) = nr_res?;
    cleanup_tasks.extend(nr_cleanup_tasks);
    cleanup_receivers.extend(nr_cleanup_receivers);
    temp_files.extend(nr_temp_files);


    let (nr_m8_streams, nr_m8_rx) = t_junction(
        ReceiverStream::new(nr_refined_m8_stream_out),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "nr_m8".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nr_m8_rx);

    let mut nr_m8_streams_iter = nr_m8_streams.into_iter();
    let nr_m8_merge = nr_m8_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_m8_map = nr_m8_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let (nr_hitsummary_streams, nr_hitsummary_rx) = t_junction(
        ReceiverStream::new(nr_refined_hit_summary_stream_out),
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nr_hitsummary_rx);

    let mut nr_hitsummary_streams_iter = nr_hitsummary_streams.into_iter();
    let nr_hit_summary_merge = nr_hitsummary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_hit_summary_preload = nr_hitsummary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_hit_summary_taxid = nr_hitsummary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    //prelod
    let mut nr_alignment_per_read: HashMap<String, SpeciesAlignmentResults> = HashMap::with_capacity(80_000_000);

    let mut preload_stream = ReceiverStream::new(nr_hit_summary_preload);
    while let Some(item) = preload_stream.next().await {
        let bytes = item.to_bytes()?;
        let line = String::from_utf8_lossy(&bytes);
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 10 {
            warn!("Malformed NR hit summary line during preload: {}", trimmed);
            continue;
        }

        let read_id = fields[0].to_string();
        let contig_taxid = fields[9].parse::<Taxid>().ok();
        let read_taxid = fields[3].parse::<Taxid>().ok();

        nr_alignment_per_read.insert(
            read_id,
            SpeciesAlignmentResults {
                contig: contig_taxid,
                read: read_taxid,
            },
        );
    }
    info!("Preloaded {} NR alignments into memory for merge", nr_alignment_per_read.len());


    let compute_merged_taxon_handle = tokio::spawn(compute_merged_taxon_counts(
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
        nr_alignment_per_read
    ));


    let refined_combined_path = out_dir.join(rename_file_path(&sample_base_buf, None, Some("refined_taxon_counts_with_dcr.json"), "_"));
    let (_refined_combined_path, refined_write_json_task) = combine_taxon_counts(
        &nt_refined_counts,
        &nr_refined_counts,
        refined_combined_path,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("combine_taxon_counts failed: {}", e)))?;
    cleanup_tasks.push(refined_write_json_task);


    // Build refined accession maps
    let nt_refined_map = collect_m8_to_accession_map(ReceiverStream::new(nt_m8_map)).await
        .map_err(|e| PipelineError::Other(anyhow!("NT refined accession map failed: {}", e)))?;

    let nr_refined_map = collect_m8_to_accession_map(ReceiverStream::new(nr_m8_map)).await
        .map_err(|e| PipelineError::Other(anyhow!("NR refined accession map failed: {}", e)))?;


    // Contig FASTA stream (from assembly)
    // ───────────────────────────────────────────────────────────────
    let contigs_fasta = assembly_outputs.contigs_fasta.clone();
    let contigs_file = tokio::fs::File::open(&contigs_fasta).await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to open contigs.fasta: {}", e)))?;

    let contigs_rx = parse_fasta(contigs_file, 32768).await
        .map_err(|e| PipelineError::Other(anyhow!("parse_fasta failed: {}", e)))?;

    let mut contigs_stream = ReceiverStream::new(contigs_rx);

    // Empty cluster stream (contigs have no duplicates)
    let (_cluster_tx, cluster_rx) = mpsc::channel::<ParseOutput>(1);
    drop(_cluster_tx);
    let contigs_cluster_stream = ReceiverStream::new(cluster_rx);

    let assembly_dir = out_dir.join("assembly");
    tokio::fs::create_dir_all(&assembly_dir).await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create assembly dir: {}", e)))?;


    let (
        mapped_contigs_rx,
        unidentified_contigs_rx,
        unique_unidentified_rx,
        mut annot_tasks,
        mut annot_rxs,
    ) = generate_annotated_fasta(
        config.clone(),
        contigs_stream,
        contigs_cluster_stream,
        nt_refined_map.into_iter().collect(),
        nr_refined_map.into_iter().collect(),
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


    let (taxid_mapped_streams, taxid_mapped_rx) = t_junction(
        ReceiverStream::new(taxid_mapped_rx),
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "taxid_mapped".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(taxid_mapped_rx);

    let mut taxid_mapped_streams_iter = taxid_mapped_streams.into_iter();
    let taxid_mapped_file = taxid_mapped_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let taxid_mapped_locator = taxid_mapped_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let taxid_mapped_nonhost = taxid_mapped_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let mapped_path = assembly_dir.join("refined_taxid_annot_mapped_only_fasta");
    let combined_path = assembly_dir.join("refined_taxid_annot_fasta");


    let write_mapped_handle = write_fasta_stream_to_file(
        ReceiverStream::new(taxid_mapped_file),
        mapped_path.clone(),
        config.base_buffer_size,
    );
    cleanup_tasks.push(write_mapped_handle);

    let write_combined_handle = write_fasta_stream_to_file(
        ReceiverStream::new(taxid_combined_rx),
        combined_path.clone(),
        config.base_buffer_size,
    );
    cleanup_tasks.push(write_combined_handle);


    let assembly_dir = out_dir.join("assembly");

    let (locator_outputs, mut locator_tasks, _locator_receivers) = generate_taxid_locator(
        taxid_mapped_locator, // directly pass the receiver
        assembly_dir,
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("generate_taxid_locator failed: {}", e)))?;

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
        ReceiverStream::new(initial_annotated_taxon_stream),
        ReceiverStream::new(initial_unidentified_taxon_stream),
        ReceiverStream::new(nt_initial_stream),
        ReceiverStream::new(nr_initial_stream),
        lineage_map.clone(),
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Experimental generate_taxid_fasta failed: {}", e)))?;

    cleanup_tasks.push(initial_load_nt_task);
    cleanup_tasks.push(initial_load_nr_task);
    cleanup_tasks.push(initial_taxid_main_task);

    let initial_mapped_path = out_dir.join("taxid_annot_mapped_only.fasta");
    let initial_combined_path = out_dir.join("taxid_annot.fasta");

    let (initial_taxid_mapped_streams, initial_taxid_mapped_rx) = t_junction(
        ReceiverStream::new(initial_taxid_mapped_rx),
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "initial_taxid_mapped".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(initial_taxid_mapped_rx);

    let mut initial_taxid_mapped_streams_iter = initial_taxid_mapped_streams.into_iter();
    let initial_taxid_mapped_file = initial_taxid_mapped_streams_iter.next().ok_or(PipelineError::EmptyStream)?; // Optional write
    let initial_taxid_mapped_locator = initial_taxid_mapped_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let initial_taxid_mapped_count = initial_taxid_mapped_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let (initial_taxid_combined_streams, initial_taxid_combined_rx) = t_junction(
        ReceiverStream::new(initial_taxid_combined_rx),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        config.base_backpressure_pause,
        StreamDataType::IlluminaFastq,
        "initial_taxid_combined".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(initial_taxid_combined_rx);

    let mut initial_taxid_combined_streams_iter = initial_taxid_combined_streams.into_iter();
    let initial_taxid_combined_file = initial_taxid_combined_streams_iter.next().ok_or(PipelineError::EmptyStream)?; // Optional write
    let initial_taxid_combined_count = initial_taxid_combined_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let initial_write_mapped_handle = write_fasta_stream_to_file(
        ReceiverStream::new(initial_taxid_mapped_file), // Clone rx if needed for logging
        initial_mapped_path.clone(),
        config.base_buffer_size,
    );
    cleanup_tasks.push(initial_write_mapped_handle);

    let initial_write_combined_handle = write_fasta_stream_to_file(
        ReceiverStream::new(initial_taxid_combined_file),
        initial_combined_path.clone(),
        config.base_buffer_size,
    );
    cleanup_tasks.push(initial_write_combined_handle);


    let initial_mapped_count = stream_record_counter(initial_taxid_mapped_count, false).await?;
    let initial_combined_count = stream_record_counter(initial_taxid_combined_count, false).await?;
    info!("Experimental: {} mapped, {} combined records", initial_mapped_count, initial_combined_count);


    let (initial_locator_outputs, mut initial_locator_tasks, _initial_locator_receivers) = generate_taxid_locator(
        initial_taxid_mapped_locator,
        out_dir.clone(),
    )
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Experimental generate_taxid_locator failed: {}", e)))?;

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


    let clusters_opt = if READ_COUNTING_MODE == ReadCountingMode::CountAll {
        Some(duplicate_clusters.clone())
    } else {
        None
    };

    // shutting this off as: 1. needletail is more picky than seqkit and finds mismatcxhes too easily
    // 2. It is unlikely we'll be allowed to keep the non host FASTQ's
    // let (nonhost_r1, nonhost_r2_opt) = generate_nonhost_fastq_from_files(
    //     config.clone(),
    //     file1_path.clone(),
    //     file2_path.clone(),
    //     ReceiverStream::new(taxid_mapped_nonhost),
    //     clusters_opt,
    //     out_dir.clone(),
    // ).await?;


    // *******************
    // Results retrieval
    // *******************

    let raw_count = join_with_error_handling(raw_count_task).await?;
    info!("Processed {} raw reads (additive from R1 and R2 if paired)", raw_count);

    let stats = join_with_error_handling(val_count_task).await?;
    info!("Processed {} validated, {} undersized, {} oversized reads",
             stats.validated, stats.undersized, stats.oversized);

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


    let ercc_mapped_count = ercc_count_rx.await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC count receiver failed: {}", e)))?;

    ercc_bt2_bam_write_handle.await??;
    host_bt2_bam_write_handle.await??;
    if let Some(handle) = optional_human_bam_write_handle {
        handle.await??;
    }

    let host_bt2_counts = host_bt2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host bt2 counts receiver failed: {}", e)))?;
    info!("Host bt2: mapped counts: {:?}", host_bt2_counts);


    // insert size, now host BAM written out
    if paired {
        info!("Computing host insert size statistics from {}", host_bt2_bam_path.display());

        let stats = compute_insert_size_stats_from_bam(
            host_bt2_bam_path.clone(),
            None,
            &config.thread_pool,
        )
            .await
            .context("Failed to compute insert size stats")?;

        // Write Picard-style metrics file
        let metrics_path = out_dir.join("host_insert_metrics.txt");
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

        let plot_path = out_dir.join("host_insert_size_histogram.png");
        if let Err(e) = plot_insert_sizes(&stats.insert_sizes, sample_base.as_str(), &plot_path) {
            warn!("Failed to plot insert size histogram: {}", e);
        }

        info!("Host insert size stats written: mean {:.2}, median {:.0}, {} proper pairs",
          stats.mean, stats.median, stats.total_proper_pairs);
    }

    let host_hisat2_counts = host_hisat2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host hisat2 counts receiver failed: {}", e)))?;
    info!("Host hisat2: mapped counts: {:?}", host_hisat2_counts);


    // Await Kallisto exit and process results. Allow graceful exit even if kallisto finds nothing.
    // Cannot asusme ERCC's spiked in.
    let kallisto_exit = join_with_error_handling(kallisto_exit_task).await;
    match kallisto_exit {
        Ok(_) => info!("Kallisto completed successfully"),
        Err(e) => {
            warn!("Kallisto failed (non-fatal, possibly no ERCC spiked in): {}", e);
        }
    }

    let kallisto_results_task = kallisto_results(out_dir.join("kallisto"), kallisto_ercc_tx)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    join_with_error_handling(kallisto_results_task).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Collect Kallisto results
    let kallisto_ercc_counts = kallisto_ercc_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC counts receiver failed: {}", e)))?;

    let kallisto_dir = out_dir.join("kallisto");
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
    info!("Wrote Kallisto abundance.tsv with {} ERCC transcripts", kallisto_ercc_counts.ercc_counts.len());

    let ercc_path = out_dir.join("kallisto_ERCC_counts_tsv");
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

    let compute_merged_taxon_result = compute_merged_taxon_handle.await??;

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
