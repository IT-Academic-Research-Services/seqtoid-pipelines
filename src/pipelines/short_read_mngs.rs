use std::cmp::Reverse;
use std::cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd};
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader, ErrorKind};
use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use ahash::AHashMap;
use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use dashmap::DashMap;
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
use tokio::io::{AsyncBufReadExt, AsyncReadExt, AsyncWriteExt, BufReader as TokioBufReader, BufWriter};
use tokio::process::Command;
use tokio::sync::mpsc;
use tokio::sync::mpsc::{channel, Sender};
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio::time::{sleep, Duration, Instant};
use tokio::try_join;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;
use tokio_util::io::StreamReader;
use twox_hash::XxHash64;


use crate::config::defs::{DiamondSubcommand, KallistoSubcommand, Lineage, PipelineError, ReadCountingMode,
                          ReadStats, RunConfig, SamtoolsStats, SamtoolsSubcommand, StreamDataType,
                          Taxid, BCFTOOLS_TAG, BLASTN_TAG, BLASTX_TAG, BOWTIE2_TAG,
                          CZID_DEDUP_TAG, DIAMOND_TAG, FASTP_TAG, HISAT2_TAG, KALLISTO_TAG,
                          KRAKEN2_TAG, LOG_NORMAL_POSITIVE_DOUBLE, MAFFT_TAG, MAKEBLASTDB_TAG,
                          MINIMAP2_TAG, MIN_NORMAL_POSITIVE_DOUBLE, NR_TAG,
                          NT_TAG, PIGZ_TAG, QUAST_TAG, READ_COUNTING_MODE,
                          SAMTOOLS_TAG, SEQKIT_TAG, SPADES_TAG, STAR_TAG};
use crate::utils::blast::{consensus_level, generate_taxon_count_json_from_m8, AggBucket, M8Record,
                          TaxonCount};
use crate::utils::command::blastn::{BlastnArgGenerator, BlastnConfig};
use crate::utils::command::blastx::{BlastxArgGenerator, BlastxConfig};
use crate::utils::command::bowtie2::{bowtie2_index_prep, Bowtie2Config};
use crate::utils::command::czid_dedup::CzidDedupConfig;
use crate::utils::command::diamond::{diamond_index_prep, DiamondArgGenerator, DiamondConfig};
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::hisat2::{hisat2_index_prep, Hisat2Config};
use crate::utils::command::kallisto::KallistoConfig;
use crate::utils::command::makeblastdb::{MakeblastdbArgGenerator, MakeblastdbConfig};
use crate::utils::command::minimap2::{minimap2_index_prep, Minimap2ArgGenerator, Minimap2Config};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::fastx::{compare_read_ids, parse_header, raw_read_count, read_fasta,
                          read_fastq, stream_record_counter, write_fasta_stream_to_file, SequenceRecord};
use crate::utils::file::{available_space_for_path, choose_temp_dir, file_path_manipulator,
                         file_size, rename_file_path, resolve_optional_path,
                         validate_file_inputs, write_byte_stream_to_file, write_parse_output_to_temp_file,
                         write_vecu8_to_file};
use crate::utils::paf::PafRecord;
use crate::utils::plotting::plot_insert_sizes;
use crate::utils::sambam::generate_info_from_bam_stream;
use crate::utils::stats::parse_samtools_stats;
use crate::utils::streams::{create_fifo, deinterleave_fastq_stream, interleave_fastq_streams, join_with_error_handling,
                            parse_child_output, parse_fastq, read_child_output_to_vec, spawn_cmd,
                            stream_to_cmd, stream_to_file, t_junction, write_to_fifo,
                            ChannelReader, ChildStream, ParseMode,
                            ParseOutput, ToBytes};
use crate::utils::streams::deinterleave_fastq_stream_to_fifos;
use crate::utils::taxonomy::{build_should_keep_filter, get_top_m8_nr,
                             get_top_m8_nt, load_taxid_lineages_db, validate_taxid_lineage};

const UNMAPPED_HEADER_PREFIX: &str = ">NR::NT::";
const MAX_ACCESSION_SEQUENCE_LEN: u64 = 100_000_000;
const EST_BYTES_PER_ACCESSION: u64 = 20_000; // ~10k seq + header
const MAX_PARSE_ERRORS: usize = 100; // Threshold before failing
const MAX_SPADES_WORK_DIR: u64 = 500_000_000;

const MIN_REF_FASTA_SIZE: u64 = 25;
const MIN_ASSEMBLED_CONTIG_SIZE: u64 = 25;

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
}

impl ReadHit {
    /// Convert a ReadHit into the exact tab-separated line format used by the original CZID pipeline
    /// This matches what HitSummaryWriter writes in the Python code
    pub fn to_tab_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            // 1. read_id — not stored in ReadHit, will be added by caller
            "", // placeholder — we fill it in update_read_dict
            self.level,
            self.taxid,
            &self.accession_id,
            "", // alignment length — not used in hit summary
            "", // percent identity — not used
            "", // bitscore — not used
            "", // evalue — not used
            self.species_taxid,
            self.genus_taxid,
            self.family_taxid,
            "", // name — not used
            "", // common_name — not used
            self.contig_id.as_deref().unwrap_or(""),
            self.contig_accession_id.as_deref().unwrap_or(""),
            self.contig_species_taxid,
            self.contig_genus_taxid,
            self.contig_family_taxid,
            "", // contig_length — not used
            "", // contig_name — not used
            "", // contig_accession — duplicate
            "", // contig_taxid — not used
            "", // contig_species_taxid — duplicate
            "", // contig_genus_taxid — duplicate
            "", // contig_family_taxid — duplicate
            if self.from_assembly { "1" } else { "0" },
            "1", // count — always 1 per read
            "1"  // count_type — not used downstream
        )
    }
}

#[derive(Default, Clone, Debug)]
pub struct AccessionHit {
    pub taxid: i32, // i.e. species_taxid
    pub genus_taxid: i32,
    pub family_taxid: i32,
    pub count: u64,
}

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
        100,
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
    output_bam_path: Option<PathBuf>,
    compute_insert_stats: bool,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, Option<oneshot::Receiver<SamtoolsStats>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let paired_compute_insert_stats = paired && compute_insert_stats;

    // BT2
    let bt2_config_view = Bowtie2Config {
        bt2_index_path: bt2_index_path.clone(),
        paired: paired,
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

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted for fastq extraction
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
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

    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 } + if paired_compute_insert_stats { 1 } else { 0 };
    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = t_junction(
        bam_rx_stream,
        num_tees,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
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
    let bam_file_stream = if output_bam_path.is_some() {
        Some(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?)
    } else {
        None
    };
    let insert_stats_stream = if paired_compute_insert_stats {
        Some(bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?)
    } else {
        None
    };

    // Handle BAM output if requested
    if let (Some(bam_path), Some(stream)) = (output_bam_path, bam_file_stream) {
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count mapped reads
    let (count_tx, count_rx) = oneshot::channel();
    let count_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(count_stream);
        let mut line_count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(line) if !line.is_empty() => line_count += 1,
                _ => continue,
            }
        }
        let count = line_count as u64 / if paired { 8 } else { 4 };
        count_tx.send(count).map_err(|_| anyhow!("Failed to send count"))?;
        Ok(())
    });
    cleanup_tasks.push(count_task);

    // Samtools fastq to extract unmapped reads
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
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

    // Compute insert size metrics if paired and compute_insert_stats is true
    let insert_stats_rx = if paired_compute_insert_stats {
        let insert_stats_stream = insert_stats_stream.unwrap();

        // Split insert_stats_stream for counting and processing
        let (mut stats_streams, stats_done_rx) = t_junction(
            ReceiverStream::new(insert_stats_stream),
            2,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes,
            "insert_stats_split".to_string(),
            None,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?;
        cleanup_receivers.push(stats_done_rx);

        let mut stats_streams_iter = stats_streams.into_iter();
        let count_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let process_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

        // Count mapped reads in BAM
        let (count_tx, count_rx) = oneshot::channel();
        let count_task = tokio::spawn(async move {
            let mut stream = ReceiverStream::new(count_stream);
            let mut line_count = 0;
            while let Some(item) = stream.next().await {
                match item {
                    ParseOutput::Bytes(line) if !line.is_empty() => line_count += 1,
                    _ => continue,
                }
            }
            count_tx.send(line_count).map_err(|_| anyhow!("Failed to send BAM line count"))?;
            Ok(())
        });
        cleanup_tasks.push(count_task);

        let bam_line_count = count_rx
            .await
            .map_err(|e| PipelineError::Other(anyhow!("BAM count receiver failed: {}", e)))?;

        if bam_line_count == 0 {
            warn!("No mapped reads in BAM for insert size stats; skipping samtools stats");
            let (stats_tx, stats_rx) = oneshot::channel();
            stats_tx
                .send(SamtoolsStats {
                    summary: HashMap::new(),
                    insert_sizes: Vec::new(),
                })
                .map_err(|_| PipelineError::Other(anyhow!("Failed to send empty stats")))?;
            Some(stats_rx)
        } else {
            // Coordinate-sort BAM for samtools stats
            let samtools_coord_sort_config = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Sort,
                subcommand_fields: HashMap::from([
                    ("-u".to_string(), None),
                    ("-O".to_string(), Some("bam".to_string())),
                    ("-".to_string(), None),
                ]),
            };
            let samtools_coord_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_coord_sort_config))
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;

            let (mut coord_sort_child, coord_sort_task, coord_sort_err_task) = stream_to_cmd(
                config.clone(),
                process_stream,
                SAMTOOLS_TAG,
                samtools_coord_sort_args,
                StreamDataType::JustBytes,
                config.args.verbose,
            )
                .await
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;
            cleanup_tasks.push(coord_sort_task);
            cleanup_tasks.push(coord_sort_err_task);

            let coord_sort_out_stream = {
                let mut guard = coord_sort_child.lock().await;
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

            // Run samtools stats
            let samtools_stats_config = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Stats,
                subcommand_fields: HashMap::new(),
            };
            let samtools_stats_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_stats_config))
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;

            let (mut stats_child, stats_stream_task, stats_err_task) = stream_to_cmd(
                config.clone(),
                coord_sort_out_stream,
                SAMTOOLS_TAG,
                samtools_stats_args,
                StreamDataType::JustBytes,
                config.args.verbose,
            )
                .await
                .map_err(|e| PipelineError::ToolExecution {
                    tool: SAMTOOLS_TAG.to_string(),
                    error: e.to_string(),
                })?;
            cleanup_tasks.push(stats_stream_task);
            cleanup_tasks.push(stats_err_task);

            let stats_out_stream = {
                let mut guard = stats_child.lock().await;
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

            // Parse samtools stats
            let (stats_tx, stats_rx) = oneshot::channel();
            let stats_task = tokio::spawn(async move {
                let stats = parse_samtools_stats(stats_out_stream)
                    .await
                    .map_err(|e| anyhow!("Failed to parse samtools stats: {}", e))?;
                stats_tx.send(stats).map_err(|_| anyhow!("Failed to send samtools stats"))?;
                Ok(())
            });
            cleanup_tasks.push(stats_task);

            Some(stats_rx)
        }
    } else {
        None
    };

    Ok((
        ReceiverStream::new(unmapped_fastq_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        insert_stats_rx,
    ))
}


/// QC's input stream using FASTP
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
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
        100,
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
        // Retry reading abundance.tsv (handles NVMe flush delays)
        let mut retries = 5;
        loop {
            match fs::metadata(&abundance_path).await {
                Ok(meta) if meta.len() > 0 => {
                    debug!("Found abundance.tsv: {} bytes", meta.len());
                    break;
                }
                _ => {
                    if retries == 0 {
                        return Err(anyhow!("No abundance.tsv found at {}", abundance_path.display()));
                    }
                    warn!("abundance.tsv not ready (retry {}); sleeping 2s", retries);
                    sleep(Duration::from_secs(2)).await;
                    retries -= 1;
                }
            }
        }

        let abundance_file = TokioFile::open(&abundance_path).await
            .map_err(|e| anyhow!("Failed to open abundance.tsv: {}", e))?;
        let reader = tokio::io::BufReader::new(abundance_file);
        let mut lines = reader.lines();
        let mut ercc_counts = Vec::new();
        let mut transcript_to_gene = Vec::new();

        debug!("Starting to read abundance.tsv from: {}", abundance_path.display());

        if lines.next_line().await?.is_none() {
            return Err(anyhow!("Empty abundance.tsv"));
        }

        // Parse each line: target_id, length, eff_length, est_counts, tpm
        while let Some(line) = lines.next_line().await? {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                return Err(anyhow!("Invalid abundance.tsv line: {}", line));
            }
            let target_id = fields[0].to_string();
            let est_counts = fields[3].parse::<f64>()
                .map_err(|e| anyhow!("Failed to parse est_counts in line {}: {}", line, e))?;
            ercc_counts.push((target_id.clone(), est_counts));
            transcript_to_gene.push((target_id, format!("gene_{}", fields[0])));
        }

        debug!("Finished reading kallisto's abundance.tsv: {} transcripts", ercc_counts.len());
        let results = KallistoResults {
            ercc_counts,
            transcript_to_gene,
        };
        ercc_tx.send(results)
            .map_err(|_| anyhow!("Failed to send Kallisto results"))?;
        debug!("Sent Kallisto results");
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
        Option<NamedTempFile>,
        Option<NamedTempFile>,
        Option<TempDir>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

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

    let minimap2_config = Minimap2Config {
        minimap2_index_path: host_index_path,
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
        100,
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

    Ok((
        ReceiverStream::new(unmapped_fastq_stream),
        count_rx,
        cleanup_tasks,
        cleanup_receivers,
        host_ref_temp,
        host_index_temp,
        host_temp_dir
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
    seed: u64, // Used for hashing, rng from config
    max_subsample: u64,
    prefix_len: Option<usize>,
    out_dir: PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, ReceiverStream<ParseOutput>, HashMap<String, u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let mut count_map: HashMap<u64, u64> = HashMap::new();
    let mut dedup_map: HashMap<u64, (u64, Vec<SequenceRecord>)> = HashMap::new();
    let mut total_count = 0u64;

    // RAM-based deduplication
    let mut stream = ReceiverStream::new(input_stream);
    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Fastq(record) => {
                let mut records = vec![record.clone()];
                let mut seq_bytes = (*record.seq()).to_vec(); // Clone Vec<u8> to owned Vec<u8>
                if paired {
                    if let Some(next_item) = stream.next().await {
                        match next_item {
                            ParseOutput::Fastq(next_record) => {
                                records.push(next_record.clone());
                                seq_bytes.extend_from_slice(&next_record.seq()); // Extend with &[u8]
                            }
                            _ => {
                                return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in paired stream".to_string()));
                            }
                        }
                    } else {
                        return Err(PipelineError::InvalidFastqFormat("Unpaired read in paired stream".to_string()));
                    }
                }
                let mut hasher = XxHash64::with_seed(seed);
                let write_len = prefix_len.map(|len| len.min(seq_bytes.len())).unwrap_or(seq_bytes.len());
                hasher.write(&seq_bytes[0..write_len]);
                let hash_key = hasher.finish();
                let entry = dedup_map.entry(hash_key).or_insert((0, records.clone()));
                entry.0 += 1;
                entry.1 = records;
                total_count += 1;
            }
            _ => {
                return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in dedup stream".to_string()));
            }
        }
    }
    count_map = dedup_map.iter().map(|(&k, &(w, _))| (k, w)).collect();
    debug!("RAM-based count complete: {} unique hashes from {} total {}", count_map.len(), total_count, if paired { "pairs" } else { "reads" });

    let subsample_size = max_subsample.min(count_map.len() as u64);

    let mut heap: BinaryHeap<Reverse<SampleItem>> = BinaryHeap::new();
    let mut rng = config.rng.clone(); // Clone the seeded RNG from config

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
        // All uniques
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

    let mut duplicate_cluster_sizes: HashMap<String, u64> = HashMap::new();
    for (&hash_key, &(weight, ref records)) in &dedup_map {
        if !records.is_empty() {
            duplicate_cluster_sizes.insert(records[0].id().to_string(), weight);
        }
    }

    let (cluster_tx, cluster_rx) = mpsc::channel(config.base_buffer_size);

    tokio::spawn(async move {
        for (&hash_key, &(weight, ref records)) in &dedup_map {
            if records.is_empty() { continue; }
            let rep_id = records[0].id();
            let rep_line = format!("{},{}\n", rep_id, rep_id);
            if cluster_tx.send(ParseOutput::Bytes(Arc::new(rep_line.into_bytes()))).await.is_err() {
                break;
            }
            for rec in records.iter().skip(1) {
                let line = format!("{},{}\n", rec.id(), rep_id);
                if cluster_tx.send(ParseOutput::Bytes(Arc::new(line.into_bytes()))).await.is_err() {
                    break;
                }
            }
        }
    });

    let tsv_path = out_dir.join("duplicate_cluster_sizes.tsv");
    let duplicate_cluster_sizes_clone = duplicate_cluster_sizes.clone(); // move into task
    let tsv_write_task = tokio::spawn(async move {
        let mut file = TokioOpenOptions::new()
            .write(true)
            .create(true)
            .open(&tsv_path)
            .await
            .map_err(|e| anyhow!("TSV open error: {}", e))?;

        for (read_id, weight) in duplicate_cluster_sizes_clone {
            file.write_all(format!("{}\t{}\n", read_id, weight).as_bytes())
                .await
                .map_err(|e| anyhow!("TSV write error: {}", e))?;
        }

        Ok(())
    });
    cleanup_tasks.push(tsv_write_task);

    let unique_count = sampled.len() as u64;
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

    Ok((ReceiverStream::new(rx), count_rx, ReceiverStream::new(cluster_rx), duplicate_cluster_sizes, cleanup_tasks, cleanup_receivers))
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
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
        Option<NamedTempFile>,
        Option<NamedTempFile>,
        Option<TempDir>
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let (_non_host_ref_fasta_path, non_host_index_path, non_host_ref_temp, non_host_index_temp, non_host_temp_dir,mut non_host_ref_tasks) =
        minimap2_index_prep(
            &config,
            &config.ram_temp_dir,
            config.args.target_sequence.clone(),
            config.args.target_index.clone(),
            "non-host",
        )
            .await?;
    try_join_all(non_host_ref_tasks).await?;

    debug!("Non-host index: {}", non_host_index_path.display());

    let minimap2_config = Minimap2Config {
        minimap2_index_path: non_host_index_path,
        option_fields: HashMap::from([
            ("-cx".to_string(), Some("sr".to_string())), // Short-read preset, PAF out
            ("--secondary".to_string(), Some("yes".to_string())), // allow secondary alignments
        ]),
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
        StreamDataType::IlluminaFastq,
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
            ParseMode::Lines,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    Ok((
        ReceiverStream::new(minimap2_out_stream),
        cleanup_tasks,
        cleanup_receivers,
        non_host_ref_temp,
        non_host_index_temp,
        non_host_temp_dir
    ))
}

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
        100,
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
    min_aln_len: u64
) -> Result<(
    ReceiverStream<ParseOutput>,
    ReceiverStream<ParseOutput>,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
)> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let (dedup_tx, dedup_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);
    let (summary_tx, summary_rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);

    let processing_task = tokio::spawn(async move {

        let lineage_map_clone = lineage_map.clone();
        let acc2taxid_map_clone = acc2taxid_map.clone();
        let should_keep_clone = should_keep_filter.clone();
        let mut read_groups: HashMap<String, Vec<M8Record>> = HashMap::new();

        while let Some(item) = m8_input.next().await {
            let line = match item {
                ParseOutput::Bytes(b) => {
                    let bytes = b.to_vec();
                    match String::from_utf8(bytes) {
                        Ok(s) => s.trim_end().to_string(),
                        Err(e) => {
                            debug!("m8 line not UTF-8: {}", e);
                            continue;
                        }
                    }
                }
                _ => continue,
            };

            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            let rec = match M8Record::parse_line_nr(&line) {  // ←short-read alignment processing,12 column, thus  nr, not nt
                Ok(r) => r,
                Err(e) => {
                    warn!("Failed to parse GSNAP m8 line: {} — {}", e, line);
                    continue;
                }
            };

            if rec.alen < min_aln_len {
                continue;
            }

            read_groups.entry(rec.qname.clone()).or_default().push(rec);
        }

        for (read_id, hits) in read_groups {
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
                } else {
                    // debug!("Missing taxid for accession: {}", hit.tname);
                }
            }

            if valid_hits.is_empty() {
                continue;
            }

            valid_hits.sort_by(|a, b| b.bitscore.partial_cmp(&a.bitscore).unwrap_or(Ordering::Equal));
            let best = valid_hits[0].clone();

            let (tax_level, consensus_taxid, consensus_hits) = consensus_level(
                &hits,
                &*lineage_map_clone,  // ← &AHashMap
                &*acc2taxid_map_clone,
                &should_keep_clone,
            )?;

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
            dedup_tx
                .send(ParseOutput::Bytes(Arc::new((dedup_line + "\n").into_bytes())))
                .await?;

            // Extract S/G/F from first hit's lineage
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
            summary_tx
                .send(ParseOutput::Bytes(Arc::new(summary_line.into_bytes())))
                .await?;
        }

        Ok(())
    });
    cleanup_tasks.push(processing_task);

    let m8_dedup_file_path = config.out_dir.join(rename_file_path(&sample_base_buf, None, Some("dedup.m8"), "_"));
    let summary_file_path = config.out_dir.join(rename_file_path(&sample_base_buf, None, Some("summary.txt"), "_"));

    let dedup_stream = ReceiverStream::new(dedup_rx);
    let (mut dedup_branches, dedup_done) = t_junction(
        dedup_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
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
        100,
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

/// Non-host DIAMOND alignment – write to a temp file, then stream it.
async fn diamond_non_host_align(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    sample_base: String,
) -> Result<(
    mpsc::Receiver<ParseOutput>,
    Vec<JoinHandle<Result<(), anyhow::Error>>>,
    Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    Option<NamedTempFile>,
    Option<NamedTempFile>,
    Option<TempDir>,
), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let diamond_db = config
        .args
        .diamond_db
        .clone()
        .ok_or(PipelineError::MissingArgument(
            "diamond_db required for NR alignment".into(),
        ))?;
    let (db_prefix, prep_tasks) = diamond_index_prep(Some(diamond_db), "non_host").await?;
    cleanup_tasks.extend(prep_tasks);


    let temp_output = NamedTempFile::new_in(&config.ram_temp_dir)
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create temp file: {}", e)))?;
    let temp_path = temp_output.path().to_path_buf();


    let diamond_options = HashMap::from([
        ("--mid-sensitive".to_string(), None),
        ("-f".to_string(), Some("6".to_string())),
        ("-o".to_string(), Some(temp_path.to_string_lossy().to_string())),
    ]);

    let diamond_config = DiamondConfig {
        subcommand: DiamondSubcommand::Blastx,
        db: db_prefix,
        subcommand_fields: diamond_options,
    };

    let diamond_args = generate_cli(DIAMOND_TAG, &config, Some(&diamond_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: DIAMOND_TAG.to_string(),
            error: e.to_string(),
        })?;

    debug!("DIAMOND args: {:?}", diamond_args);

    let (diamond_child_arc, diamond_stream_task, diamond_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        DIAMOND_TAG,
        diamond_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await?;

    cleanup_tasks.push(diamond_err_task);

    diamond_stream_task.await??;

    let diamond_child_clone = diamond_child_arc.clone();
    let wait_task = tokio::spawn(async move {
        let mut guard = diamond_child_clone.lock().await;
        let status = guard.wait().await.map_err(|e| anyhow!("Diamond wait error: {}", e))?;
        if !status.success() {
            return Err(anyhow!("DIAMOND exited with code {:?}", status.code()));
        }
        Ok(())
    });
    wait_task.await??;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);

    let read_task = tokio::spawn(async move {
        // Open the file we just created
        let file = TokioFile::open(&temp_path)
            .await
            .map_err(|e| anyhow!("Failed to open DIAMOND output {}: {}", temp_path.display(), e))?;

        let mut reader = TokioBufReader::new(file);
        let mut line = Vec::<u8>::new();

        loop {
            line.clear();
            match reader.read_until(b'\n', &mut line).await {
                Ok(0) => {
                    info!("DIAMOND temp-file EOF");
                    break;
                }
                Ok(_) => {
                    // Clone the line into an Arc so the rest of the pipeline can share it
                    let arc_line = Arc::new(line.clone());
                    if tx.send(ParseOutput::Bytes(arc_line)).await.is_err() {
                        warn!("DIAMOND m8 receiver dropped");
                        break;
                    }
                }
                Err(e) => {
                    error!("DIAMOND temp-file read error: {}", e);
                    return Err(anyhow!("Read error: {}", e));
                }
            }
        }

        drop(temp_output);

        Ok(())
    });
    cleanup_tasks.push(read_task);


    Ok((rx, cleanup_tasks, cleanup_receivers, None, None, None))
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
    duplicate_cluster_sizes: Arc<HashMap<String, u64>>,
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
        let cluster_size = *duplicate_cluster_sizes.get(&read_id).unwrap_or(&1u64);

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
    out_dir: &PathBuf,
    sample_base: PathBuf,
) -> Result<(
    PathBuf,                     // annotated_merged.fa
    PathBuf,                     // unidentified.fa
    Option<PathBuf>,             // unique_unidentified.fa
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
)> {

    let annotated_path = out_dir.join(rename_file_path(&sample_base, None, Some("annotated_merged.fa"), "_"));
    let unidentified_path = out_dir.join(rename_file_path(&sample_base, None, Some("unidentified.fa"), "_"));
    let unique_unid_path = out_dir.join(rename_file_path(&sample_base, None, Some("unique_unidentified.fa"), "_"));

    let annotated_file = TokioFile::create(&annotated_path).await?;
    let unid_file = TokioFile::create(&unidentified_path).await?;
    let unique_file = TokioFile::create(&unique_unid_path).await?;

    let mut annotated_writer = BufWriter::with_capacity(config.base_buffer_size, annotated_file);
    let mut unid_writer = BufWriter::with_capacity(config.base_buffer_size, unid_file);
    let mut unique_writer = BufWriter::with_capacity(config.base_buffer_size, unique_file);


    let (done_tx, done_rx) = oneshot::channel();
    let mut cleanup_receivers = vec![done_rx];


    // Helper macro – writes a FASTA entry (header + seq + newline)
    macro_rules! write_fasta {
        ($writer:expr, $header:expr, $seq:expr) => {{
            $writer.write_all($header.as_bytes()).await?;
            $writer.write_all($seq).await?;   // $seq is &[u8]
            $writer.write_all(b"\n").await?;
        }};
    }

    let mut rep_buffer: HashMap<String, (Vec<u8>, bool)> = HashMap::with_capacity(1_000_000);
    let mut csv_line = String::with_capacity(256);

    let mut processed_reps = 0u64;
    let mut unidentified_reps = 0u64;
    let mut expanded_duplicates = 0u64;
    let start = Instant::now();

    let process_task = tokio::spawn(async move {
        while let Some(item) = dedup_stream.next().await {
            let record = match item {
                ParseOutput::Fastq(rec) => rec,
                _ => continue,
            };

            processed_reps += 1;
            if processed_reps % 100_000 == 0 {
                debug!("generate_annotated_fasta: processed {} reps", processed_reps);
            }

            let rep_id = record.id().to_string();
            let seq = record.seq().to_vec();
            let nr_acc = nr_map.get(&rep_id).cloned().unwrap_or_default();
            let nt_acc = nt_map.get(&rep_id).cloned().unwrap_or_default();
            let is_unidentified = nr_acc.is_empty() && nt_acc.is_empty();

            rep_buffer.insert(rep_id.clone(), (seq.clone(), is_unidentified));

            if !is_unidentified {
                let header = format!(">NR:{}:NT:{}:{}\n", nr_acc, nt_acc, rep_id);
                write_fasta!(annotated_writer, header, &seq);
            } else {
                unidentified_reps += 1;
                let header = format!("{}{}\n", UNMAPPED_HEADER_PREFIX, rep_id);
                write_fasta!(unid_writer, header, &seq);
                write_fasta!(unique_writer, header, &seq);
            }

            // Drain CSV lines for this rep
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
                        continue;
                    }

                    let line = std::mem::take(&mut csv_line);
                    let parts: Vec<&str> = line.trim().split(',').collect();
                    if parts.len() != 2 { continue; }

                    let member_id = parts[0];
                    let csv_rep_id = parts[1];

                    if csv_rep_id == rep_id && is_unidentified {
                        if let Some((seq, _)) = rep_buffer.get(&rep_id) {
                            let (base, suffix) = if member_id.ends_with("/1") || member_id.ends_with("/2") {
                                let len = member_id.len();
                                (&member_id[..len - 2], &member_id[len - 2..])
                            } else {
                                (member_id, "")
                            };
                            let header = format!("{}{}{}\n", UNMAPPED_HEADER_PREFIX, base, suffix);
                            write_fasta!(unid_writer, header, seq);
                            expanded_duplicates += 1;
                        }
                    }
                }
            }
        }

        annotated_writer.flush().await?;
        unid_writer.flush().await?;
        unique_writer.flush().await?;

        info!(
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
        annotated_path,
        unidentified_path,
        Some(unique_unid_path),
        cleanup_tasks,
        cleanup_receivers,
    ))
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
    input_stream: ReceiverStream<ParseOutput>,
    out_dir: PathBuf,
    paired: bool,
    sample_base_buf: PathBuf,
)  -> Result<(
    AssemblyHandle,
    Vec<JoinHandle<Result<()>>>,
    Vec<oneshot::Receiver<Result<()>>>,
),  PipelineError>

{
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let spades_out_dir = out_dir.join("assembly");
    fs::create_dir_all(&spades_out_dir).await?;
    let spades_work_dir =  config.ram_temp_dir.join("spades");
    fs::create_dir_all(&spades_work_dir).await?;

    let input_file_path = spades_work_dir.join(rename_file_path(&sample_base_buf, None, Some("assembly_input.fq"), "_"));

    let mut total_input_bytes : u64 = 0;

    let write_input_task = tokio::spawn({
        let input_file_path = input_file_path.clone();

        async move {
            let file = TokioFile::create(&input_file_path).await
                .map_err(|e| anyhow!("Cannot create SPAdes input file {}: {}", input_file_path.display(), e))?;
            let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, file);

            let mut stream = input_stream;
            let mut total_bytes = 0u64;

            while let Some(item) = stream.next().await {

                let bytes = item.to_bytes()
                    .map_err(|e| anyhow!("Failed to convert ParseOutput to FASTQ bytes: {}", e))?;

                writer.write_all(&bytes).await
                    .map_err(|e| anyhow!("Write error to {}: {}", input_file_path.display(), e))?;

                total_bytes += bytes.len() as u64;
            }

            writer.flush().await
                .map_err(|e| anyhow!("Flush error to {}: {}", input_file_path.display(), e))?;

            info!(
            "SPAdes input written: {} ({} bytes)",
            input_file_path.display(),
            total_bytes
        );

            if total_bytes == 0 {
                warn!("SPAdes input file is empty – no non-host reads");
            }
            total_input_bytes = total_bytes;
            Ok::<PathBuf, anyhow::Error>(input_file_path)
        }

    });

    // Decide SPAdes work dir based on actual input size
    let spades_work_dir = if total_input_bytes <= MAX_SPADES_WORK_DIR {
        config.ram_temp_dir.join("spades")
    } else {
        warn!("SPAdes input {} bytes ({:.1} MB) → too large for RAM, using disk temp",
          total_input_bytes, total_input_bytes as f64 / 1e6);
        std::env::temp_dir().join("spades")
    };

    fs::create_dir_all(&spades_work_dir).await?;

    let input_file_path = write_input_task.await??;

    let spades_config = SpadesConfig {
        input_path: input_file_path.clone(),
        outdir_path: spades_out_dir.clone(),
        paired: paired,

        option_fields: HashMap::from([
            ("--only-assembler".to_string(), None),
            ("--isolate".to_string(), None),
        ]),
    };
    let spades_args = generate_cli(SPADES_TAG, &config, Some(&spades_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SPADES_TAG.to_string(),
            error: e.to_string(),
        })?;


    let (mut spades_child, spades_err_task) = spawn_cmd(
        config.clone(),
        SPADES_TAG,
        spades_args,
        config.args.verbose,
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SPADES_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(spades_err_task);


    let spades_task = tokio::spawn(async move {
        let output = spades_child.wait_with_output().await?;
        eprintln!("SPAdes exit: {:?}", output.status);
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            eprintln!("SPAdes FAILED:\n{}", stderr);
            return Err(anyhow!("SPAdes failed: {}", stderr));
        }
        Ok(())
    });

    Ok((AssemblyHandle {
        spades_task,
        work_dir: spades_work_dir,
        out_dir: spades_out_dir}, cleanup_tasks, cleanup_receivers))
}


pub async fn process_assembly(
    config: Arc<RunConfig>,
    out_dir: &PathBuf,
    work_dir: &PathBuf,
    bowtie_stream: ReceiverStream<ParseOutput>,
    duplicate_cluster_sizes: Arc<HashMap<String, u64>>,
    paired: bool,
    assembly_headroom: u64,
) -> Result<(
           CoverageOutputs,
           ReceiverStream<ParseOutput>,
           Vec<JoinHandle<Result<()>>>,
           Vec<oneshot::Receiver<Result<()>>>,
           Vec<NamedTempFile>
),  PipelineError>

 {
    let mut cleanup_tasks = Vec::new();let mut cleanup_receivers = Vec::new();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();

    let raw_contigs_path    = work_dir.join("contigs.fasta");
    let raw_scaffolds_path  = work_dir.join("scaffolds.fasta");

     let (empty_tx, empty_rx) = mpsc::channel::<ParseOutput>(1);

    let spades_success = if raw_contigs_path.exists() {
        fs::metadata(&raw_contigs_path).await?.len() > 0
    } else {
        false
    };

    // if it failed write empty data
    if !spades_success {

        let dummy = ";ASSEMBLY FAILED";
        let empty_sam = "@NO INFO\n";
        let empty_json = "{}";

        fs::write(out_dir.join("contigs.fasta"), dummy).await?;
        fs::write(out_dir.join("contigs.ram.fasta"), dummy).await?;
        fs::write(out_dir.join("contigs_all.fasta"), dummy).await?;
        fs::write(out_dir.join("scaffolds.fasta"), dummy).await?;
        fs::write(out_dir.join("read-contig.sam"), empty_sam).await?;
        fs::write(out_dir.join("contig_stats.json"), empty_json).await?;

        return Ok((CoverageOutputs {
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
        }, ReceiverStream::new(empty_rx), cleanup_tasks,  cleanup_receivers, temp_files));
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
     ).await?;

     let ram_fasta_path = NamedTempFile::with_suffix_in(".fa", &temp_dir)
         .map_err(|e| PipelineError::Other(e.into()))?;
     let ram_path = ram_fasta_path.path().to_owned(); // PathBuf, so its owned

    //contig length filtering
    let rx = read_fasta(
        raw_contigs_path.clone(),
        u64::MAX,           // no record limit
        Some(config.args.min_contig_length),
        None,
        8192,
    ).map_err(|e| PipelineError::InvalidFastaFormat(e.to_string()))?;

     let (write_handle, _ram_path) = write_fasta_stream_to_file(
         ReceiverStream::new(rx),
         ram_path.clone(),
         config.base_buffer_size,
     ).await?;

     write_handle.await
         .map_err(|e| PipelineError::Other(anyhow!("FASTA write task panicked: {}", e)))?
         .map_err(|e| PipelineError::Other(anyhow!("FASTA write failed: {}", e)))?;

     fs::copy(&ram_path, contigs_out).await
         .map_err(|e| PipelineError::Other(anyhow!("Failed to copy contigs to output: {}", e)))?;

     temp_files.push(ram_fasta_path);

    let scaffolds_out = out_dir.join("scaffolds.fasta");
    let spades_scaffolds = work_dir.join("scaffolds.fasta");
    if spades_scaffolds.exists() {
        fs::copy(&spades_scaffolds, &scaffolds_out).await?;
    } else {
        fs::write(&scaffolds_out, ";NO SCAFFOLDS").await?;
    }

    // bt2 index
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
        return Err(PipelineError::ToolExecution {tool: BOWTIE2_TAG.to_string(), error: build_status.to_string()});
    }

    let bt2_options = HashMap::from([
        ("--very-sensitive".to_string(), None),
    ]);

    // BT2
    let bt2_config_view = Bowtie2Config {
        bt2_index_path: index_prefix,
        paired: paired,
        option_fields: bt2_options,
    };

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;
    eprintln!("Bowtie 2 args: {:?}", bt2_args);
    let (mut bt2_child, bt2_stream_task, bt2_err_task) = stream_to_cmd(
        config.clone(),
        bowtie_stream.into_inner(),
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

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted for fastq extraction
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())), // uncompressed BAM for noodles to consume
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
         100,
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

     let (read2contig, contig_stats) = generate_info_from_bam_stream(
         bam_for_stats,
         &duplicate_cluster_sizes,
         config.args.min_contig_length,
     ).await?;

    Ok((CoverageOutputs {
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
    }, ReceiverStream::new(bam_for_output), cleanup_tasks,  cleanup_receivers, temp_files))



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
    let parts: Vec<String> = header_line[1..]
        .split('\x01')
        .map(|part| {
            let trimmed = part.trim_start();
            FIX_COMMA_REGEXP.replace_all(trimmed, " ").to_string()
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
pub async fn build_reference_fasta_from_selected_genera(
    config: Arc<RunConfig>,
    selected_genera: &HashMap<i32, Vec<String>>,  // genus_taxid → [accessions]
    db_type: &str,
    source_fasta_path: &PathBuf,                  // full NT/NR FASTA (e.g. nt.fa)
    source_index_path: &PathBuf,                  // FST index for fast lookup
) -> Result<PathBuf> {
    // Flatten and deduplicate accessions
    let mut accessions: HashSet<String> = HashSet::new();
    for acc_list in selected_genera.values() {
        accessions.extend(acc_list.iter().cloned());
    }

    if accessions.is_empty() {
        warn!("selected_genera is empty — writing empty reference FASTA");
        let temp_dir = choose_temp_dir(1024, &config.ram_temp_dir, &config.args.nvme_scratch, 4).await?;
        let empty_path = temp_dir.join(format!("{}_empty_ref.fasta", db_type));
        fs::write(&empty_path, b">empty\nN\n").await?;
        return Ok(empty_path);
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

    // Estimate size and pick best temp location
    let estimated_size = (acc_vec.len() as u64) * 2000; // ~2KB per accession avg
    let temp_dir = choose_temp_dir(
        estimated_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
    ).await?;

    let output_path = temp_dir.join(format!("{}_ref_from_selected_genera.fasta", db_type));
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
    drop(writer); // ensure written

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

    Ok(output_path)
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
        if parts.len() < 7 {
            warn!("Malformed hit_summary line (<7 fields): {}", line);
            continue;
        }

        let read_id = parts[0].to_string();

        let level = match parts[1].parse::<u8>() {
            Ok(l) => l,
            Err(_) => {
                warn!("Invalid level field: {}", parts[1]);
                continue;
            }
        };

        let taxid = match parts[2].parse::<i32>() {
            Ok(t) => t,
            Err(_) => {
                warn!("Invalid taxid field: {}", parts[2]);
                continue;
            }
        };

        let accession_id = if parts[3] == "-" || parts[3].is_empty() {
            "-".to_string()
        } else {
            parts[3].to_string() // keep full versioned accession
        };

        let species_taxid = parts[4].parse::<i32>().unwrap_or(0);
        let genus_taxid = parts[5].parse::<i32>().unwrap_or(0);
        let family_taxid = parts[6].parse::<i32>().unwrap_or(0);

        read_dict.insert(
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
///
/// # Arguments
/// * `read2contig` - Assignment of reads to contigs
/// * `top_m8_stream` - stream in m8 format, top as found in get_top_m8_nt/nr
/// * `read_dict` - HashMap of ReadID -> ReadHits
/// * `accession_dict` - accession → (species, genus, family)
/// * `updated_read_tx` - sendiner

/// # Returns
/// Result of hashset of the accessions
pub async fn update_read_dict(
    read2contig: Arc<HashMap<String, String>>,
    mut top_m8_stream: ReceiverStream<ParseOutput>,
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    accession_map: Arc<AHashMap<String, AccessionHit>>,
    should_keep: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    db_type: &str,
    mut contig2lineage_tx: Sender<ParseOutput>,
    mut read2blastm8_tx: Sender<ParseOutput>,
    mut updated_tx: Sender<ParseOutput>,
    mut added_tx: Sender<ParseOutput>,
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
        });

        let hit_tab_suffix = shared_hit.to_tab_string();
        // Precompute the base hit summary line


        for read_id in reads_in_contig {
            let final_line = format!("{}\t{}", read_id, &hit_tab_suffix);

            let was_present = {
                let mut dict = read_dict.lock().unwrap();
                let was_present = dict.contains_key(&read_id);
                // This is now just an atomic ref-count increment — essentially free
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
    duplicate_cluster_sizes: Arc<HashMap<String, u64>>,
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

            let cluster_size = duplicate_cluster_sizes
                .get(read_id)
                .copied()
                .unwrap_or(1);

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


pub async fn blast_contigs (
    config: Arc<RunConfig>,
    db_type: &'static str,
    deduped_m8_stream: ReceiverStream<ParseOutput>, // deduped_m8
    hit_summary_stream: ReceiverStream<ParseOutput>, // deduped_m8
    read_dict: Arc<Mutex<AHashMap<String, Arc<ReadHit>>>>,
    accession_map: Arc<AHashMap<String, AccessionHit>> ,  //hit_summary, but already done
    taxon_counts: Vec<TaxonCount>, // orig_counts_with_dcr
    assembled_contig_fasta: &PathBuf, //assembled_contig, CoverageOutputs -> contigs_ram_fasta
    read2contig: Arc<HashMap<String, String>>,
    reference_fasta: &PathBuf, // reference_fasta
    duplicate_cluster_sizes: Arc<HashMap<String, u64>>, //duplicate_cluster_sizes_path
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
    should_keep_filter: Arc<impl Fn(&[i32]) -> bool + Send + Sync + 'static>,
    blast_headroom: u64


) -> Result<( AHashMap<String, Arc<ReadHit>>, Vec<TaxonCount>, Vec<ContigSummaryEntry>, Vec<JoinHandle<Result<()>>>,
             Vec<oneshot::Receiver<Result<()>>>, Vec<NamedTempFile> )>{
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
        ).await?;

        let final_read_dict = read_dict.lock().unwrap().clone();
        return Ok(( final_read_dict, taxon_counts, vec![], cleanup_tasks, cleanup_receivers, temp_files ));
    }

    let temp_dir = choose_temp_dir(
        ref_size + contig_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        blast_headroom,
    ).await?;

    let blastdb_suffix = format!("{}_blastindex", db_type);
    let blastdb_ram_path = NamedTempFile::with_suffix_in(blastdb_suffix, &temp_dir)
        .map_err(|e| PipelineError::Other(e.into()))?;
    let blastdb_path = blastdb_ram_path.path().to_owned();
    temp_files.push(blastdb_ram_path);

    let makeblastdb_config = MakeblastdbConfig {
        input: reference_fasta.clone(),
        dbtype: if db_type == NT_TAG { "nucl".to_string() } else { "prot".to_string() },
        output: blastdb_path.clone(),
        option_fields: HashMap::new()
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
    ).await
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
            outfmt: "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".to_string(),
            evalue: 1e-10,
            max_target_seqs: 5000,
            option_fields: HashMap::new(),
        };

        generate_cli(BLASTN_TAG, &config, Some(&blastn_config))
            .map_err(|e| PipelineError::ToolExecution {
                tool: BLASTN_TAG.to_string(),
                error: e.to_string(),
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

        generate_cli(BLASTX_TAG, &config, Some(&blastx_config))
            .map_err(|e| PipelineError::ToolExecution {
                tool: BLASTX_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // let envs = if db_type == NT_TAG {
    //     vec![("BATCH_SIZE".to_string(), "10000".to_string())]
    // } else {
    //     vec![]
    // };

    let (mut blast_child, err_task) = spawn_cmd(
        config.clone(),
        blast_command,
        blast_args,
        config.args.verbose,
    ).await
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
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "blastn/x".to_string(),
            error: e.to_string(),
        })?;

    // Now feed blast_out_stream directly into top_m8
    let (top_tx, top_rx) = channel(1024);
    let top_handle = if db_type == NT_TAG {
        tokio::spawn(get_top_m8_nt(ReceiverStream::new(blast_out_stream), top_tx))
    } else {
        tokio::spawn(get_top_m8_nr(ReceiverStream::new(blast_out_stream), top_tx))
    };

    let (contig2lineage_tx, contig2lineage_rx) = channel(1024);
    let (read2blastm8_tx, read2blastm8_rx) = channel(1024);

    let (updated_tx, updated_rx) = channel(1024);
    let (added_tx, added_rx) = channel(1024);

    let update_handle = tokio::spawn(update_read_dict(
        read2contig.clone(),
        ReceiverStream::new(top_rx),
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

    let (refined_m8_tx, refined_m8_rx) = channel(2_000_000);
    let (refined_hit_summary_tx, refined_hit_summary_rx) = channel(2_000_000);

    let (updated_tx, updated_rx) = channel(1024);
    let (added_tx, added_rx) = channel(1024);

    let generate_m8_handle = tokio::spawn(generate_m8_and_hit_summary(
        ReceiverStream::new(updated_rx),        // not used any more – we’ll ignore it
        ReceiverStream::new(added_rx),         // not used any more
        ReceiverStream::new(read2blastm8_rx),   // new BLAST hits (read → best hit)
        hit_summary_stream,                   // original hit summary (teed)
        deduped_m8_stream,                   // original deduped M8
        refined_m8_tx.clone(),
        refined_hit_summary_tx.clone(),
    ));

    let (refined_counts_tx, refined_counts_rx) = channel(1024);
    let counts_handle = tokio::spawn(generate_taxon_count_json_from_m8(
        ReceiverStream::new(refined_m8_rx),
        ReceiverStream::new(refined_hit_summary_rx),
        db_type,
        lineage_map.clone(),
        should_keep_filter.clone(),
        duplicate_cluster_sizes.clone(),
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
            let contig  = fields[0].to_string();
            let species  = fields[1].parse().unwrap_or(0);
            let genus    = fields[2].parse().unwrap_or(0);
            let family   = fields[3].parse().unwrap_or(0);
            contig2lineage.insert(contig, [species, genus, family]);
        }
    }


    let (contig_summary_tx, contig_summary_rx) = channel(1024);
    let contig_summary_handle = tokio::spawn(generate_contig_summary_json(
        read2contig.clone(),
        contig2lineage,
        read_dict.clone(),
        db_type,
        duplicate_cluster_sizes.clone(),
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
        cleanup_tasks,
        cleanup_receivers,
        temp_files
    ))
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
    let mut temp_dirs: Vec<TempDir> = Vec::new();

    info!("Starting short read mNGS pipeline.");

    // *******************
    // Setup and Validation
    // *******************

    // External tools check
    check_versions(vec![BOWTIE2_TAG, MINIMAP2_TAG, KALLISTO_TAG, SPADES_TAG, MAKEBLASTDB_TAG,
                        BLASTN_TAG, BLASTX_TAG])
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Check required files
    let host_bowtie2_index: String = config.args.host_bowtie2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_bowtie2_index is required".to_string()))?;

    let (file1_path, file2_path, sample_base_buf, sample_base) = validate_file_inputs(&config, &cwd)?;
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
        file1_path,
        file2_path,
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
    let (ercc_bt2_out_stream, ercc_count_rx, ercc_bt2_cleanup_tasks, ercc_bt2_cleanup_receivers, _ercc_bt2_insert_stats_rx) = bowtie2_filter(
        config.clone(),
        val_out_stream,
        ercc_bt2_index_path,
        paired,
        ercc_bt2_options,
        None,
        false
    ).await?;
    cleanup_tasks.extend(ercc_bt2_cleanup_tasks);
    cleanup_receivers.extend(ercc_bt2_cleanup_receivers);

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
        100,
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
    let (host_bt2_out_stream, host_bt2_count_rx, host_bt2_cleanup_tasks, host_bt2_cleanup_receivers, host_bt2_insert_stats_rx) = bowtie2_filter(
        config.clone(),
        ReceiverStream::new(kallisto_bypass_stream),
        host_bt2_index_path,
        paired,
        host_bt2_options,
        None,
        true
    ).await?;
    cleanup_tasks.extend(host_bt2_cleanup_tasks);
    cleanup_receivers.extend(host_bt2_cleanup_receivers);

    // Host filtering: mm2
    let (host_mm2_out_stream, host_mm2_count_rx, mut host_mm2_cleanup_tasks, mut host_mm2_cleanup_receivers, host_ref_temp, host_index_temp, host_temp_dir) = minimap2_filter(
        config.clone(),
        host_bt2_out_stream,
        paired,
        None, // No BAM output; set to Some(PathBuf::from("host_mm2.bam")) if needed
    )
        .await?;
    cleanup_tasks.append(&mut host_mm2_cleanup_tasks);
    cleanup_receivers.append(&mut host_mm2_cleanup_receivers);

    if let Some(temp) = host_ref_temp {
        temp_files.push(temp);
    }
    if let Some(temp) = host_index_temp {
        temp_files.push(temp);
    }
    if let Some(temp) = host_temp_dir {
        temp_dirs.push(temp);
    }

    // If host is no huma, run an additional filter stage using a human reference
    let post_filter_stream  = if config.args.human_host {
        host_mm2_out_stream
    }
    else {
        let human_bowtie2_index: String = config.args.human_bowtie2_index.clone();

        // human filtering: bt2
        let human_bt2_index_path = bowtie2_index_prep(human_bowtie2_index, &cwd)?;
        let human_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
        let (human_bt2_out_stream, human_bt2_count_rx, human_bt2_cleanup_tasks, human_bt2_cleanup_receivers, _human_bt2_insert_stats_rx) = bowtie2_filter(
            config.clone(),
            host_mm2_out_stream,
            human_bt2_index_path,
            paired,
            human_bt2_options,
            None,
            false
        ).await?;
        cleanup_receivers.extend(human_bt2_cleanup_receivers);

        // human filtering: mm2
        let (human_mm2_out_stream, human_mm2_count_rx, mut human_mm2_cleanup_tasks, mut human_mm2_cleanup_receivers, human_ref_temp, human_index_temp, human_index_dir) = minimap2_filter(
            config.clone(),
            human_bt2_out_stream,
            paired,
            None, // No BAM output; set to Some(PathBuf::from("human_mm2.bam")) if needed
        )
            .await?;
        cleanup_tasks.append(&mut human_mm2_cleanup_tasks);
        cleanup_receivers.append(&mut human_mm2_cleanup_receivers);

        if let Some(temp) = human_ref_temp {
            temp_files.push(temp);
        }
        if let Some(temp) = human_index_temp {
            temp_files.push(temp);
        }
        if let Some(temp) = human_index_dir {
            temp_dirs.push(temp);
        }

        human_mm2_out_stream

    };

    let (dedup_stream, dedup_count_rx, cluster_stream,duplicate_cluster_sizes_noarc, dedup_cleanup_tasks, dedup_cleanup_receivers) = dedup_and_subsample(
        config.clone(),
        post_filter_stream.into_inner(),
        paired,
        seed,
        config.args.max_subsample as u64,
        Some(75), // Prefix length for deduplication. Hardcoded fior now
        out_dir.clone(),
    ).await?;
    cleanup_tasks.extend(dedup_cleanup_tasks);
    cleanup_receivers.extend(dedup_cleanup_receivers);

    let duplicate_cluster_sizes = Arc::new(duplicate_cluster_sizes_noarc);

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


    // split inpout stream for mm2 and dmnd
    let (non_host_streams, non_host_done_rx) = t_junction(
        dedup_stream,
        5,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,
        "non_host_output".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(non_host_done_rx);

    let mut non_host_streams_iter = non_host_streams.into_iter();
    let non_host_mm2_stream = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let non_host_dmnd_stream = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let non_host_annot_stream = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let non_host_assembly_stream = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let non_host_coverage_bt2_stream = non_host_streams_iter.next().ok_or(PipelineError::EmptyStream)?;



    // This is part of post-process starting here for concurrency

    let (assembly_handle, mut assembly_cleanup_tasks, mut assembly_cleanup_receivers) = spades_assembly(
        config.clone(),
        ReceiverStream::new(non_host_assembly_stream),
        out_dir.clone(),
        paired,
        sample_base_buf.clone(),
    ).await?;

    cleanup_tasks.append(&mut assembly_cleanup_tasks);
    cleanup_receivers.append(&mut assembly_cleanup_receivers);


    // Minimap2 non_host alignment
    let (non_host_mm2_out_stream, mut non_host_mm2_cleanup_tasks, mut non_host_mm2_cleanup_receivers, non_host_ref_temp, non_host_index_temp, non_host_index_dir) = minimap2_non_host_align(
        config.clone(),
        ReceiverStream::new(non_host_mm2_stream),
    )
        .await?;
    cleanup_tasks.append(&mut non_host_mm2_cleanup_tasks);
    cleanup_receivers.append(&mut non_host_mm2_cleanup_receivers);

    if let Some(temp) = non_host_ref_temp {
        temp_files.push(temp);
    }
    if let Some(index_temp) = non_host_index_temp {
        temp_files.push(index_temp);
    }
    if let Some(index_temp) = non_host_index_dir {
        temp_dirs.push(index_temp);
    }

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

    let (call_stream, call_summary_stream, mut call_cleanup_tasks, mut call_cleanup_receivers) = call_hits_m8_stream(
        config.clone(),
        m8_stream,
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        36
    ).await?;
    cleanup_tasks.append(&mut call_cleanup_tasks);
    cleanup_receivers.append(&mut call_cleanup_receivers);

    let (nt_streams, nt_done_rx) = t_junction(
        call_stream,
        4,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,
        "nt_call".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(nt_done_rx);

    let mut nt_streams_iter =nt_streams.into_iter();
    let nt_call_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_m8_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_acc_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nt_blast_stream = nt_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (nt_summary_streams, nt_summary_done_rx) = t_junction(
        call_summary_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
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
    let nt_blast_hit_stream = nt_summary_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let nt_hit_summary_handle = tokio::spawn(summarize_hits(ReceiverStream::new(nt_summary_hit_stream), 0));

    let nt_counts = generate_taxon_counts(
        config.clone(),
        ReceiverStream::new(nt_call_stream),
        ReceiverStream::new(nt_summary_taxon_stream),
        duplicate_cluster_sizes.clone(),
        should_keep_filter.clone(),
        "NT".to_string(),
        None,
    ).await?;
    info!("NT taxon counts: {} entries", nt_counts.len());

    let nt_map = collect_m8_to_accession_map(ReceiverStream::new(nt_m8_stream)).await?;


    // Diamond non_host alignment
    let (non_host_diamond_m8_stream, mut non_host_diamond_cleanup_tasks, mut non_host_diamond_cleanup_receivers, diamond_ref_temp, diamond_index_temp, diamond_index_dir) = diamond_non_host_align(
        config.clone(),
        ReceiverStream::new(non_host_dmnd_stream),
        paired,
        sample_base
    ).await?;
    cleanup_tasks.append(&mut non_host_diamond_cleanup_tasks);
    cleanup_receivers.append(&mut non_host_diamond_cleanup_receivers);


    let (nr_call_stream, nr_call_summary_stream, mut nr_call_cleanup_tasks, mut nr_call_cleanup_receivers) = call_hits_m8_stream(
        config.clone(),
        ReceiverStream::new(non_host_diamond_m8_stream),
        sample_base_buf.clone(),
        lineage_map.clone(),
        acc2taxid_map.clone(),
        should_keep_filter.clone(),
        0,
    ).await?;
    cleanup_tasks.append(&mut nr_call_cleanup_tasks);
    cleanup_receivers.append(&mut nr_call_cleanup_receivers);

    let (nr_streams, nr_done_rx) = t_junction(
        nr_call_stream,
        4,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
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
    let nr_acc_stream = nr_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let nr_blast_stream = nr_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let (nr_summary_streams, nr_summary_done_rx) = t_junction(
        nr_call_summary_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
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
        duplicate_cluster_sizes.clone(),
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

    let (_annotated_path, _unidentified_path, _unique_path_opt, mut annot_tasks, mut annot_rxs) =
        generate_annotated_fasta(
            config.clone(),
            ReceiverStream::new(non_host_annot_stream),
            cluster_stream,
            nt_map, nr_map,
            &out_dir, sample_base_buf.clone(),
        ).await?;
    cleanup_tasks.append(&mut annot_tasks);
    cleanup_receivers.append(&mut annot_rxs);

    // *******************
    // Post-processing
    // *******************

    // Assembly stats

    let assembly_out_dir = &assembly_handle.out_dir;
    let assembly_work_dir =  &assembly_handle.work_dir;

    assembly_handle.spades_task.await??;

    let (assembly_outputs, assembly_bam_out_stream,
        mut post_assembly_cleanup_tasks,
        mut post_assembly_cleanup_receivers, mut post_assemly_temp_files) = process_assembly(
        config.clone(),
        assembly_out_dir,
        assembly_work_dir,
        ReceiverStream::new(non_host_coverage_bt2_stream),
        duplicate_cluster_sizes.clone(),
        paired,
        4 // this can become as CLI arg if needed
    ).await?;

    cleanup_tasks.extend(post_assembly_cleanup_tasks);
    cleanup_receivers.extend(post_assembly_cleanup_receivers);
    temp_files.extend(post_assemly_temp_files);

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
    let nt_ref_fasta_path = build_reference_fasta_from_selected_genera(
        config.clone(),
        &nt_selected_genera,
        NT_TAG,
        &PathBuf::from(nt_file),
        &PathBuf::from(nt_offset_db_file),
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    eprintln!("nt path {}", nt_ref_fasta_path.display());


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


    let nr_ref_fasta_path = build_reference_fasta_from_selected_genera(
        config.clone(),
        &nr_selected_genera,
        NR_TAG,
        &PathBuf::from(nr_file),
        &PathBuf::from(nr_offset_db_file),
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    eprintln!("nr path {}", nr_ref_fasta_path.display());

    let nt_handle = tokio::spawn({
        let contigs_fasta_path = assembly_outputs.contigs_ram_fasta.clone();
        let config = config.clone();
        let lineage_map = lineage_map.clone();
        let should_keep_filter = should_keep_filter.clone();
        let duplicate_cluster_sizes = duplicate_cluster_sizes.clone();
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
                duplicate_cluster_sizes,
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
        let duplicate_cluster_sizes = duplicate_cluster_sizes.clone();
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
                duplicate_cluster_sizes,
                lineage_map,
                should_keep_filter,
                4,
            ).await
        }
    });



    //  Do everything else that does NOT depend on the results
    //    (e.g. writing intermediate JSONs, logging, cleanup prep…)



    let (nt_res, nr_res) = tokio::try_join!(nt_handle, nr_handle)
        .map_err(|e| PipelineError::Other(anyhow!("blast_contigs task panicked: {e}")))?;

    let (nt_read_dict, nt_refined_counts, nt_contig_summary,
        nt_cleanup_tasks, nt_cleanup_receivers,
        nt_temp_files) = nt_res?;

    cleanup_tasks.extend(nt_cleanup_tasks);
    cleanup_receivers.extend(nt_cleanup_receivers);
    temp_files.extend(nt_temp_files);

    let (nr_read_dict, nr_refined_counts, nr_contig_summary,
        nr_cleanup_tasks, nr_cleanup_receivers,
        nr_temp_files) = nr_res?;
    cleanup_tasks.extend(nr_cleanup_tasks);
    cleanup_receivers.extend(nr_cleanup_receivers);
    temp_files.extend(nr_temp_files);


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

    let host_bt2_counts = host_bt2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host bt2 counts receiver failed: {}", e)))?;
    info!("Host bt2: mapped counts: {:?}", host_bt2_counts);

    if let Some(insert_stats_rx) = host_bt2_insert_stats_rx {
        match insert_stats_rx.await {
            Ok(insert_stats) if !insert_stats.insert_sizes.is_empty() => {
                info!("Host bt2 insert size stats: {:?}", insert_stats.insert_sizes);
                let output_path = out_dir.join("insert_size_histogram.png");
                let sample_name = sample_base_buf
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .to_string();
                let plot_task = tokio::spawn(async move {
                    plot_insert_sizes(&insert_stats.insert_sizes, &sample_name, &output_path)
                        .map_err(|e| anyhow!("Failed to plot insert sizes: {}", e))?;
                    Ok(())
                });
                cleanup_tasks.push(plot_task);
            }
            Ok(_) => warn!("No insert size stats available; skipping plotting"),
            Err(e) => error!("Host bt2 insert stats receiver failed: {}", e),
        }
    }

    let host_mm2_counts = host_mm2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host minimap2 counts receiver failed: {}", e)))?;
    info!("Host minimap2: mapped counts: {:?}", host_mm2_counts);

    let dedup_count = dedup_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("dedup count receiver failed: {}", e)))?;
    info!("Dedup: unique reads/pairs: {}", dedup_count);

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
    // info!("Kallisto ERCC counts: {:?}", kallisto_ercc_counts);

    // *******************
    // Cleanup
    // *******************
    let results = try_join_all(cleanup_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }
    for receiver in cleanup_receivers {
        receiver
            .await
            .map_err(|e| PipelineError::Other(e.into()))?
            .map_err(|e| PipelineError::Other(e))?;
    }

    drop(temp_files);
    drop(temp_dirs);

    info!("Finished short read mNGS pipeline.");
    Ok(())
}
