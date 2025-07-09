use std::path::PathBuf;
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


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SamtoolsSubcommand {
    View,
    Fastq,
    Stats,
    Sort,
    Index,
    Mpileup,
    Consensus
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