pub mod fastq;
mod file;
mod streams;

use std::io::{self, Read};
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

}

#[derive(Debug, Clone, ValueEnum, Default, PartialEq)]
pub enum Technology {
    #[default]
    illumina,
    ONT
}
