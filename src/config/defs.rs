// src/config/defs.rs
use std::path::PathBuf;
use std::collections::HashMap;
use std::cmp::min;
use std::sync::Arc;

use anyhow::anyhow;
use lazy_static::lazy_static;
use rayon::ThreadPool;
use tokio::sync::Semaphore;
use rand::rngs::StdRng;
use tokio::task::JoinError;
use serde_json::Error as SerdeJsonError;
use serde::{Deserialize, Serialize};
use log::LevelFilter;
use once_cell::sync::Lazy;

use crate::cli::Arguments;

// External software
pub const GZIP_EXT: &str = "gz";
pub const FASTP_TAG: &str = "fastp";
pub const PIGZ_TAG: &str = "pigz";
pub const H5DUMP_TAG: &str = "h5dump";
pub const MINIMAP2_TAG: &str = "minimap2";
pub const SAMTOOLS_TAG: &str = "samtools";
pub const BCFTOOLS_TAG: &str = "bcftools";
pub const KRAKEN2_TAG: &str = "kraken2";
pub const IVAR_TAG: &str = "ivar";
pub const MUSCLE_TAG: &str = "muscle";
pub const MAFFT_TAG: &str = "mafft";
pub const QUAST_TAG: &str = "quast.py";
pub const NUCMER_TAG: &str = "nucmer";
pub const SHOW_COORDS_TAG: &str = "show-coords";
pub const SEQKIT_TAG: &str = "seqkit";
pub const BOWTIE2_TAG: &str = "bowtie2";
pub const HISAT2_TAG: &str = "hisat2";
pub const KALLISTO_TAG: &str = "kallisto";
pub const STAR_TAG: &str = "STAR";
pub const CZID_DEDUP_TAG: &str = "czid-dedup";
pub const DIAMOND_TAG: &str = "diamond";
pub const SPADES_TAG: &str = "spades.py";
pub const BLASTN_TAG: &str = "blastn";
pub const BLASTX_TAG: &str = "blastx";
pub const MAKEBLASTDB_TAG: &str = "makeblastdb";
pub const SORT_TAG: &str = "sort";
pub const MMSEQS_TAG: &str = "mmseqs";

pub const NT_TAG: &str = "nt";
pub const NR_TAG: &str = "nr";

// Taxonomy defs
pub type Taxid = i32;
pub type Lineage = [i32; 3];
pub const INVALID_CALL_BASE_ID: i32 = -100;

pub const LOG_NORMAL_POSITIVE_DOUBLE: f64 = 1e-200;
pub const MIN_NORMAL_POSITIVE_DOUBLE: f64 = f64::MIN_POSITIVE;

pub const CONFORMING_PREAMBLE: &str = ">family_nr:-300:family_nt:-300:genus_nr:-200:genus_nt:-200:species_nr:-100:species_nt:-100:";

pub const SMT_MODEST_MULTIPLIER: f32 = 1.3_f32;

pub static SIMD_LEVEL: Lazy<SimdLevel> = Lazy::new(|| {
    crate::utils::system::detect_simd_level()
});

lazy_static! {
    pub static ref TOOL_VERSIONS: HashMap<&'static str, f32> = {
        let mut m = HashMap::new();
        m.insert(SAMTOOLS_TAG, 1.19);
        m.insert(BCFTOOLS_TAG, 1.19);
        m.insert(MINIMAP2_TAG, 2.24);
        m.insert(KRAKEN2_TAG, 2.1);
        m.insert(PIGZ_TAG, 2.8);
        m.insert(FASTP_TAG, 1.0);
        m.insert(MAFFT_TAG, 7.5);
        m.insert(QUAST_TAG, 5.20);
        m.insert(SEQKIT_TAG, 2.10);
        m.insert(BOWTIE2_TAG, 2.50);
        m.insert(KALLISTO_TAG, 0.5);
        m.insert(HISAT2_TAG, 2.20);
        m.insert(STAR_TAG, 2.7);
        m.insert(CZID_DEDUP_TAG, 0.1);
        m.insert(DIAMOND_TAG, 2.1);
        m.insert(SPADES_TAG, 4.2);
        m.insert(BLASTN_TAG, 2.12);
        m.insert(BLASTX_TAG, 2.12);
        m.insert(MAKEBLASTDB_TAG, 2.12);
        m.insert(MMSEQS_TAG, 0.0);
        m
    };
}

lazy_static! {
    static ref TOOL_THREAD_CAPS: HashMap<&'static str, usize> = {
        let mut m = HashMap::new();
        m.insert("bowtie2", 64);
        m.insert("hisat2", 64);
        m.insert("minimap2", 64);
        m.insert("samtools", 32);     // Sort is I/O-bound
        m.insert("spades", 128);      // Compute-heavy
        m.insert("diamond", 256);  // warp factor 10
        m.insert("fastp", 32);        // I/O-bound
        m.insert("pigz", 16);         // Compression scales poorly >16
        m.insert("kraken2", 64);      // Memory/I/O heavy
        m
    };
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SamtoolsSubcommand {
    View,
    Fastq,
    Stats,
    Sort,
    Index,
    Mpileup,
    Consensus,
    Depth,
    Ampliconclip,
    Fixmate
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BcftoolsSubcommand {
    Call,
    Consensus,
    Mpileup,
    View
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IvarSubcommand {
    Consensus
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeqkitSubcommand {
    Stats,
    Rmdup
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KallistoSubcommand {
    Index,
    Quant
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DiamondSubcommand {
    Blastx
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoreAllocation {
    Maximal,
    High,
    Low,
    Minimal
}

// Enum to specify the data type for tuning batch sizes
#[derive(Clone, Copy, Debug)]
pub enum StreamDataType {
    JustBytes,        // Streamed Vec<u8> from samtools, minimap2, BWA, etc.
    IlluminaFastq, // SequenceRecord for Illumina FASTQ or FASTA
    OntFastq,      // SequenceRecord for ONT FASTQ or FASTA
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ReadCountingMode {
    CountAll,
    CountUnique,
}
pub const READ_COUNTING_MODE: ReadCountingMode = ReadCountingMode::CountUnique;


// For specifying which read (or read pairs are validated fore size.
#[derive(Debug)]
pub struct ReadStats {
    pub undersized: u64,
    pub validated: u64,
    pub oversized: u64,
    pub unpaired_r1: u64,
    pub unpaired_r2: u64,
}

/// Samtools stats output
#[derive(Debug)]
pub struct SamtoolsStats {
    pub summary: HashMap<String, String>,
    pub insert_sizes: Vec<(u32, u64)>,
    pub total_pairs: u64
}

// Static Filenames
pub const NUCMER_DELTA: &str = "alignment.delta";

// Static Parameters
pub const IVAR_QUAL_THRESHOLD: usize = 20;
pub const IVAR_FREQ_THRESHOLD: f64 = 0.75;

pub const FASTA_TAG: &str = "fasta";
pub const FASTQ_TAG: &str = "fastq";
pub const FASTA_EXTS: &[&'static str] = &["fasta", "fa", "fna", "faa", "ffn", "frn"];
pub const FASTQ_EXTS: &[&'static str] = &["fastq", "fq"];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdLevel {
    Scalar,
    Avx2,
    Avx512,
}

#[derive(Debug, Clone)]
pub struct GpuInfo {
    pub index: usize,              // 0-based
    pub name: String,              // e.g. "NVIDIA H100 80GB HBM3" or "Apple M2"
    pub memory_mib: Option<u64>,   // total VRAM in MiB (None if unknown)
    pub is_discrete: bool,         // true for dedicated card, false for integrated
    pub driver: Option<String>,    // e.g. "550.90.07" or None
}

#[derive(Debug, Clone)]
pub struct GpuDetection {
    pub count: usize,
    pub gpus: Vec<GpuInfo>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NRAlignmentBackend {
    Diamond,          // always CPU
    MmseqsCpu,
    MmseqsGpu,
}


#[derive(Clone, Debug)]
pub struct RunConfig {
    pub cwd: PathBuf,
    pub ram_temp_dir: PathBuf,
    pub out_dir: PathBuf,
    pub args: Arguments,
    pub thread_pool: Arc<ThreadPool>,
    pub maximal_semaphore: Arc<Semaphore>,
    pub base_buffer_size: usize,
    pub input_size: u64,
    pub physical_cores: usize,
    pub max_cores: usize,
    pub available_ram: u64,
    pub rng: StdRng,
    pub log_level: LevelFilter,
    pub base_backpressure_pause: u64,
    pub simd: SimdLevel,
    pub gpu_info: GpuDetection,
    pub has_gpu: bool,
    pub alignment_backend: NRAlignmentBackend,
    pub run_id: String,
    pub efs_base_dir: PathBuf,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TaxonSeqLocation {
    pub taxid: i32,
    pub first_byte: u64,
    pub last_byte: u64,
    pub hit_type: String,
}

#[derive(Clone)]
pub struct ClusterInfo {
    pub size: u64,                  // Cluster size (used by generate_taxon_counts, process_assembly)
    pub members: Vec<String>,       // All member IDs (rep at [0]) — used by non-host in CountAll
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PairingMode {
    Strict,
    Relaxed,
}

pub type DuplicateClusters = HashMap<String, ClusterInfo>;

impl RunConfig {
    pub fn get_core_allocation(&self, tag: &str, subcommand: Option<&str>) -> CoreAllocation {
        match (tag, subcommand) {
            (MINIMAP2_TAG, _) | (KRAKEN2_TAG, _) | (MAFFT_TAG, _) | (NUCMER_TAG, _) | (FASTP_TAG, _)
            | (PIGZ_TAG, _) | (BOWTIE2_TAG, _) | (KALLISTO_TAG, _) | (DIAMOND_TAG, _) |
            (SPADES_TAG, _) | (BLASTN_TAG, _) | (BLASTX_TAG, _) | (HISAT2_TAG, _)
            | (MAKEBLASTDB_TAG, _) | (MMSEQS_TAG, _) => CoreAllocation::Maximal,
            (SAMTOOLS_TAG, Some("sort")) | (BCFTOOLS_TAG, Some("mpileup")) |
            (BCFTOOLS_TAG, Some("call")) | (QUAST_TAG, _) | (MUSCLE_TAG, _) => CoreAllocation::High,
            (SAMTOOLS_TAG, Some("view")) | (SAMTOOLS_TAG, Some("stats")) |
            (SAMTOOLS_TAG, Some("depth")) | (BCFTOOLS_TAG, Some("view")) | (SEQKIT_TAG, _) => CoreAllocation::Low,
            (IVAR_TAG, _) | (SHOW_COORDS_TAG, _) | (H5DUMP_TAG, _) => CoreAllocation::Minimal,
            _ => CoreAllocation::Minimal,
        }
    }

    pub fn mmseqs_threads(&self) -> usize {
        let base = self.thread_allocation(MMSEQS_TAG, Some("search"));

        // Scale based on physical scale of the machine
        let threads = match self.max_cores {
            0..=64   => base.min(56),                    // small machines / MacBook
            65..=128 => base.min(112),                   // r6id.32xlarge
            129..=192 => base.min(160),                  // r8id.48xlarge class
            _        => base.min(224),                   // 256+ core EPYC permanent nodes
        };

        // Final safety cap (never use 100% of cores)
        threads.min(self.max_cores * 9 / 10)
    }

    pub fn thread_allocation(&self, tag: &str, subcommand: Option<&str>) -> usize {
        let max_cores = min(self.max_cores, self.args.threads.max(1));

        let mut allocation = match self.get_core_allocation(tag, subcommand) {
            CoreAllocation::Maximal => max_cores,
            CoreAllocation::High => ((max_cores as f32 * 0.75) as usize).max(1),
            CoreAllocation::Low => (max_cores / 3).max(1),
            CoreAllocation::Minimal => 1,
        };

        let prefer_physical = match tag {
            DIAMOND_TAG | MINIMAP2_TAG | KRAKEN2_TAG | KALLISTO_TAG | MAFFT_TAG | MMSEQS_TAG => false, // allow modest SMT
            BOWTIE2_TAG | HISAT2_TAG | SPADES_TAG | FASTP_TAG | PIGZ_TAG => true,      // physical only
            _ => false,  // default allow
        };

        if prefer_physical && self.args.use_smt { // If a program does poorly with SMT, don't allow it even if use_smt true
            allocation = allocation.min(self.physical_cores);
        } else if !prefer_physical && self.args.use_smt {
            // Modest SMT: allow up to ~1.3× physical cores

            let modest_max = (self.physical_cores as f32 * SMT_MODEST_MULTIPLIER).ceil() as usize;
            allocation = allocation.min(modest_max);
        }

        // Apply tool-specific cap
        if let Some(cap) = TOOL_THREAD_CAPS.get(tag) {
            allocation = allocation.min(*cap);
        }

        // Per-tool overrides
        match (tag, subcommand) {
            (PIGZ_TAG, _) => allocation.min(16),
            (FASTP_TAG, _) => allocation.min(32),
            (BCFTOOLS_TAG, Some("mpileup")) => allocation.min(16),
            (HISAT2_TAG, _) => allocation.min(64), //scales a bit worse than bt2
            _ => allocation,
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum PipelineError {

    // File-related errors
    #[error("File not found: {0}")]
    FileNotFound(PathBuf),

    #[error("Should be a directory: {0}")]
    NotDirectory(PathBuf),

    #[error("I/O error: {0}")]
    IOError(String),

    #[error("Wrong extension: {0}")]
    WrongExtension(String),

    // Format-related errors
    #[error("Invalid FASTQ format in {0}")]
    InvalidFastqFormat(String),

    #[error("Invalid FASTA format in {0}")]
    InvalidFastaFormat(String),

    #[error("Invalid configuration: {0}")]
    InvalidConfig(String),

    // Stream-related errors
    #[error("Stream data dropped unexpectedly")]
    StreamDataDropped,

    #[error("Empty stream encountered")]
    EmptyStream,

    // Tool and execution errors
    #[error("Tool execution failed: {tool} with error: {error}")]
    ToolExecution { tool: String, error: String },

    #[error("Argument missing: {0}")]
    MissingArgument(String),

    // Bioinformatics-specific errors
    #[error("Reference sequence retrieval failed: {0}")]
    ReferenceRetrievalFailed(String),

    #[error("No sequences matched target taxonomy ID: {0}")]
    NoTargetSequences(String),

    #[error("No alignments.")]
    NoAlignments,

    // ¯\_(ツ)_/¯
    #[error("{0}")]
    Other(#[from] anyhow::Error),
}

impl From<std::io::Error> for PipelineError {
    fn from(err: std::io::Error) -> Self {
        PipelineError::IOError(err.to_string())
    }
}

impl From<JoinError> for PipelineError {
    fn from(err: JoinError) -> Self {
        PipelineError::Other(anyhow::Error::new(err))
    }
}

impl From<SerdeJsonError> for PipelineError {
    fn from(err: SerdeJsonError) -> Self {
        PipelineError::Other(anyhow::Error::new(err))
    }
}

impl From<tokio::sync::oneshot::error::RecvError> for PipelineError {
    fn from(err: tokio::sync::oneshot::error::RecvError) -> Self {
        PipelineError::Other(anyhow!("Oneshot receive failed: {}", err))
    }
}