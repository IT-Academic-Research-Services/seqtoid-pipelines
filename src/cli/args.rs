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

    #[arg(short = 't', long = "technology", default_value = "illumina", value_enum)]
    pub technology: Technology,

    #[arg(long)]
    pub min_read_len: Option<usize>,

    #[arg(long)]
    pub max_read_len: Option<usize>,

    #[arg(long, default_value_t = 64)]
    pub threads: usize,

    #[arg(short = 'q', long = "quality", default_value_t = 30)]
    pub quality: u8,

    #[arg(long, default_value_t = 10000)]
    pub stall_threshold: u64,

    #[arg(long, default_value_t = 0)]
    pub stream_sleep_ms: u64,

    #[arg(short = 'a', long = "host_accession")]  // For host removal
    pub host_accession : Option<String>,

    #[arg(long)]  // For host removal
    pub host_sequence : Option<String>,

    #[arg(short = 'r', long = "ref_accession")]  // For target aligning
    pub ref_accession : Option<String>,

    #[arg(long)] // For target aligning
    pub ref_sequence : Option<String>,

    #[arg(short = 'd', long = "db")]
    pub ref_db : Option<String>,

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
    pub ercc_sequences : String,

    #[arg( long = "target_type", default_value = "viral", value_enum)]
    pub target_type: TargetType,

    #[arg(long, default_value_t = 10)]
    pub min_depth: usize,

    #[arg(long, default_value_t = 0.0006)] // assumes about 20 mutations between 2 random samples
    pub bcftools_call_theta: f64,  // (this is an overestimate to increase sensitivity)
    
}