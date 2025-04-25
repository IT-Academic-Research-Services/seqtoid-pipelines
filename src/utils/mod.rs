pub mod fastx;
pub mod file;
pub mod streams;
pub mod command;

use clap::{Parser, ValueEnum};

#[derive(Parser, Default, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Arguments {
    #[arg(short, long)]
    pub module: String,

    #[arg(short = 'i', long = "file1", required = true)]
    pub file1: String,

    #[arg(short = 'I', long = "file2")]
    pub file2: Option<String>,

    #[arg(long, default_value_t = 50000000)]
    pub max_reads: usize,
    
    #[arg(short = 't', long = "technology", default_value = "illumina", value_enum)]
    pub technology: Technology,
    
    #[arg(long)]
    pub min_read_len: Option<usize>,
    
    #[arg(long)]
    pub max_read_len: Option<usize>,
    
    #[arg(long, default_value_t = 4)]
    pub threads: usize,
    
    #[arg(short = 'q', long = "quality", default_value_t = 30)]
    pub quality: u8,

    #[arg(long, default_value_t = 15)]
    pub stall_threshold: u64,

    #[arg(long, default_value_t = 0)]
    pub stream_sleep_ms: u64,
}



#[derive(Debug, Clone, ValueEnum, Default, PartialEq)]
pub enum Technology {
    #[default]
    Illumina,
    ONT
}
