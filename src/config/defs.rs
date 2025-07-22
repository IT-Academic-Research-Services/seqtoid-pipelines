use std::path::PathBuf;
use crate::cli::Arguments;
use lazy_static::lazy_static;
use std::collections::HashMap;

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


lazy_static! {
    pub static ref TOOL_VERSIONS: HashMap<&'static str, f32> = {
        let mut m = HashMap::new();
        m.insert(SAMTOOLS_TAG, 1.20);
        m.insert(BCFTOOLS_TAG, 1.20);
        m.insert(MINIMAP2_TAG, 2.24);
        m.insert(KRAKEN2_TAG, 2.1);
        m.insert(PIGZ_TAG, 2.8);
        m.insert(FASTP_TAG, 1.0);
        m.insert(MAFFT_TAG, 7.5);
        m.insert(QUAST_TAG, 5.20);
        m.insert(SEQKIT_TAG, 2.10);

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
    Depth
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

// Static Filenames
pub const NUCMER_DELTA: &str = "alignment.delta";


// Static Parameters

pub const IVAR_QUAL_THRESHOLD: usize = 20;
pub const IVAR_FREQ_THRESHOLD: f64 = 0.75;

pub const FASTA_TAG : &str = "fasta";
pub const FASTQ_TAG : &str = "fastq";
pub const FASTA_EXTS: &[&'static str] = &["fasta", "fa", "fna", "faa", "ffn", "frn"];
pub const FASTQ_EXTS: &[&'static str] = &["fastq", "fq"];


pub struct RunConfig  {
    pub cwd: PathBuf,
    pub ram_temp_dir: PathBuf,
    pub args: Arguments
}