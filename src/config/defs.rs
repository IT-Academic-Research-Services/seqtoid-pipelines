pub const GZIP_EXT: &str = "gz";
pub const FASTP_TAG: &str = "fastp";
pub const PIGZ_TAG: &str = "pigz";
pub const H5DUMP_TAG: &str = "h5dump";
pub const MINIMAP2_TAG: &str = "minimap2";
pub const SAMTOOLS_TAG: &str = "samtools";

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SamtoolsSubcommand {
    View,
    Fastq
}
pub const BCFTOOLS_TAG: &str = "bcftools";

pub const FASTA_TAG : &str = "fasta";
pub const FASTQ_TAG : &str = "fastq";
pub const FASTA_EXTS: &[&'static str] = &["fasta", "fa", "fna", "faa", "ffn", "frn"];
pub const FASTQ_EXTS: &[&'static str] = &["fastq", "fq"];