
use crate::utils::{Arguments, fastq};
use seq_io::fastq::{Reader, Record};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;
use flate2::read::GzDecoder;

pub fn run(args: &Arguments) {
    eprintln!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);
    println!("File 1: {}", args.file1);
    match &args.file2 {
        Some(file) => println!("File 2: {}", file),
        None => println!("File2 not given"),
    }
}