pub mod fastq;
mod file;


use std::io::{self, Read};
use clap::Parser;

#[derive(Parser, Default, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Arguments {
    #[arg(short, long)]
    pub module: String,

    #[arg(short = 'i', long = "file1", required = true)]
    pub file1: String,

    #[arg(short = 'I', long = "file2")]
    pub file2: Option<String>,

}

