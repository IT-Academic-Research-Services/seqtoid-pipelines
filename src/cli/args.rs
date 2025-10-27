use clap::{Parser, ValueEnum};

#[derive(Debug, Clone, ValueEnum, Default, PartialEq)]
pub enum Technology {
    #[default]
    Illumina,
    ONT,
}

#[derive(Debug, Clone, ValueEnum, Default, PartialEq)]
pub enum TargetType {
    #[default]
    Viral,
    Bacterial,
}

#[derive(Parser, Debug, Clone, Default)]
#[command(name = "myapp", version = "1.0")]
pub struct Arguments {

    #[arg(short, long)]
    pub module: String,

    #[arg(short = 'v', long = "verbose", action)]
    pub verbose: bool,

    #[arg(short = 'i', long = "file1")]
    pub file1: Option<String>,

    #[arg(short = 'I', long = "file2")]
    pub file2: Option<String>,

    #[arg(short = 'o', long = "out", help = "Output directory for all generated files. If not specified, a directory named '<sample_base>_YYYYMMDD' will be created in the current working directory.")]
    pub out_dir: Option<String>,

    #[arg(long, default_value_t = 50000000)]
    pub max_reads: usize,

    #[clap(long, help = "Optional fixed seed for reproducibility; defaults to OS entropy")]
    pub seed: Option<u64>,

    #[arg(long = "technology", default_value = "illumina", value_enum)]
    pub technology: Technology,

    #[arg(long)]
    pub min_read_len: Option<usize>,

    #[arg(long)]
    pub med_read_len: Option<usize>,

    #[arg(long)]
    pub max_read_len: Option<usize>,

    #[arg(long, default_value_t = 64)]
    pub threads: usize,

    #[arg(short = 'q', long = "quality", default_value_t = 30)]
    pub quality: u8,

    #[arg(long, default_value_t = 10_000_000)]
    pub stall_threshold: u64,

    #[arg(long, default_value_t = 0)]
    pub stream_sleep_ms: u64,

    #[arg(short = 'a', long)]  // For host removal
    pub host_sequence : Option<String>,

    #[arg(long, help = "Optional path to pre-built minimap2 index for host reference (e.g., hg38.mmi)")]
    pub host_index: Option<String>,

    #[arg(short = 't', long)] // For target aligning
    pub target_sequence : Option<String>,

    #[arg(long, help = "Optional path to pre-built minimap2 index for target reference (e.g., covid.mmi)")]
    pub target_index: Option<String>,

    #[arg(long)]
    pub target_taxid : Option<String>,

    #[arg(long = "index")]
    pub ref_index : Option<String>,

    #[arg(long, default_value_t = false)]
    pub limit_align_threads: bool,

    #[arg(long, default_value_t = false)]
    pub dont_filter_reads: bool,

    #[arg(short = 'k', long = "kdb")]
    pub kraken_db : Option<String>,

    #[arg(long)]
    pub adapter_fasta : Option<String>,

    #[arg(long, default_value = "ercc_sequences.fasta")]
    pub ercc_sequence : Option<String>,

    #[arg(long, help = "Optional path to pre-built minimap2 index for ERCC reference (e.g., ercc.mmi)")]
    pub ercc_index: Option<String>,

    #[arg( long = "target_type", default_value = "viral", value_enum)]
    pub target_type: TargetType,

    #[arg(long, default_value_t = 10)]
    pub min_depth: usize,

    #[arg(long, default_value_t = 0.0006)] // assumes about 20 mutations between 2 random samples
    pub bcftools_call_theta: f64,  // (this is an overestimate to increase sensitivity)

    #[arg(long, default_value = "artic_v3_primers.bed")]
    pub primer_bed_path : Option<String>,

    #[arg(long, default_value = "ercc.bowtie2.tar")]
    pub ercc_bowtie2_index : String,

    #[arg(long)]
    pub kallisto_index : Option<String>,

    #[arg(long)]
    pub host_bowtie2_index : Option<String>,

    #[arg(long, default_value = "human.bowtie2.tar")]
    pub human_bowtie2_index : String,

    #[arg(long)]
    pub host_hisat2_index : Option<String>,

    #[arg(long, default_value = "human.hisat2.tar")]
    pub human_hisat2_index : String,

    #[arg(long)]
    pub host_star_index : Option<String>,

    #[arg(long, default_value = "human.star.tar.gz")]
    pub human_star_index : String,

    #[arg(long, default_value_t = true)]
    pub human_host: bool,

    #[arg(long, default_value_t = 1000000)]
    pub max_subsample: usize,

    #[arg(long, default_value_t = 1_528_186_360_278)]
    pub nt_db_size: usize,

    #[arg(long, default_value = "taxid-lineages.db")]
    pub taxid_lineages_db: String,

    #[arg(long, default_value = "accession2_taxid_sled.db")]
    pub acc2taxid_db: String,

    #[arg(long)]
    pub taxon_whitelist: Option<String>,

    #[arg(long)]
    pub taxon_blacklist: Option<String>,

    #[arg(long)]
    pub deuterostome_list: Option<String>,

    #[clap(long)]
    pub duplicate_cluster_sizes: Option<String>,

    #[arg(long, default_value_t = 36)]
    pub min_alignment_length: u64,

    #[clap(
        long,
        value_delimiter = ',',
        value_parser = clap::value_parser!(i32),
        default_value = "23,760",
        help = "Comma-separated list of test taxids (e.g., 562,287)"
    )]
    pub test_taxids: Vec<i32>,
}