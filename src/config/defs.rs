// src/config/defs.rs
use std::path::PathBuf;
use crate::cli::Arguments;
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::cmp::min;
use std::sync::Arc;
use rayon::ThreadPool;
use tokio::sync::Semaphore;
use std::io;
use tokio::task::JoinError; // Added for JoinError
use serde_json::Error as SerdeJsonError; // Added for serde_json::Error

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
    Grep
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KallistoSubcommand {
    Index,
    Quant
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


// For specifying which read (or read pairs are validated fore size.
#[derive(Debug)]
pub struct ReadStats {
    pub undersized: u64,
    pub validated: u64,
    pub oversized: u64,
}

/// Samtools stats output
#[derive(Debug)]
pub struct SamtoolsStats {
    pub summary: HashMap<String, String>,
    pub insert_sizes: Vec<(u32, u64)>,
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

#[derive(Clone, Debug)]
pub struct RunConfig {
    pub cwd: PathBuf,
    pub ram_temp_dir: PathBuf,
    pub out_dir: PathBuf,
    pub args: Arguments,
    pub thread_pool: Arc<ThreadPool>,
    pub maximal_semaphore: Arc<Semaphore>,
    pub base_buffer_size: usize,
    pub input_size_mb: u64,
}

impl RunConfig {
    pub fn get_core_allocation(&self, tag: &str, subcommand: Option<&str>) -> CoreAllocation {
        match (tag, subcommand) {
            (MINIMAP2_TAG, _) | (KRAKEN2_TAG, _) | (MAFFT_TAG, _) | (NUCMER_TAG, _) | (FASTP_TAG, _) | (PIGZ_TAG, _) | (BOWTIE2_TAG, _) | (KALLISTO_TAG, _) => CoreAllocation::Maximal,  // Keep as-is for full potential
            (SAMTOOLS_TAG, Some("sort")) | (BCFTOOLS_TAG, Some("mpileup")) |
            (BCFTOOLS_TAG, Some("call")) | (QUAST_TAG, _) | (MUSCLE_TAG, _)  => CoreAllocation::High,
            (SAMTOOLS_TAG, Some("view")) | (SAMTOOLS_TAG, Some("stats")) |
            (SAMTOOLS_TAG, Some("depth")) | (BCFTOOLS_TAG, Some("view")) | (SEQKIT_TAG, _) => CoreAllocation::Low,
            (IVAR_TAG, _) | (SHOW_COORDS_TAG, _) | (H5DUMP_TAG, _) => CoreAllocation::Minimal,
            _ => CoreAllocation::Minimal,
        }
    }

    pub fn thread_allocation(&self, tag: &str, subcommand: Option<&str>) -> usize {
        let max_cores = min(num_cpus::get(), self.args.threads);
        let mut allocation = match self.get_core_allocation(tag, subcommand) {
            CoreAllocation::Maximal => max_cores,
            CoreAllocation::High => ((max_cores as f32 * 0.75) as usize).max(1),
            CoreAllocation::Low => (max_cores / 3).max(1),
            CoreAllocation::Minimal => 1,
        };

        match (tag, subcommand) {
            (PIGZ_TAG, _) => allocation.min(16), // Cap at 16: Compression scales poorly >16
            (FASTP_TAG, _) => allocation.min(32), // Cap at 32: QC/filtering I/O-bound
            (BCFTOOLS_TAG, Some("mpileup")) => allocation.min(16), // Cap at 16: Diminishing returns for pileup
            _ => allocation,
        }
    }

    pub fn get_buffer_size(&self, file_size_mb: u64) -> usize {
        if file_size_mb > 10_000 { // >10GB
            (self.base_buffer_size / 10).max(5_000) // ~5k-50k records (~5-50MB for Illumina)
        } else {
            self.base_buffer_size // ~100k-1M records (~100MB-1GB)
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum PipelineError {
    #[error("File not found: {0}")]
    FileNotFound(PathBuf),
    #[error("Invalid FASTQ format in {0}")]
    InvalidFastqFormat(String),
    #[error("I/O error: {0}")]
    IOError(String),
    #[error("Tool execution failed: {tool} with error: {error}")]
    ToolExecution { tool: String, error: String },
    #[error("Stream data dropped unexpectedly")]
    StreamDataDropped,
    #[error("Invalid configuration: {0}")]
    InvalidConfig(String),
    #[error("Reference sequence retrieval failed: {0}")]
    ReferenceRetrievalFailed(String),
    #[error("Empty stream encountered")]
    EmptyStream,
    #[error("No sequences matched target taxonomy ID: {0}")]
    NoTargetSequences(String),
    #[error("No alignments.")]
    NoAlignments,
    #[error("Argument missing: {0}")]
    MissingArgument(String),
    #[error("Wrong extension: {0}")]
    WrongExtension(String),
    #[error("Invalid extension: {0}")]
    Other(#[from] anyhow::Error), // Wraps external errors
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