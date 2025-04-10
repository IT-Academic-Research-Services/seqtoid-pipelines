
use crate::utils::{Arguments, fastq};
use seq_io::fastq::{Reader, Record};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;
use flate2::read::GzDecoder;
use crate::utils::fastq::read_and_interleave_fastq;

pub async fn run(args: &Arguments) -> anyhow::Result<()> {
    eprintln!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);
    println!("File 1: {}", args.file1);
    match &args.file2 {
        Some(file) => println!("File 2: {}", file),
        None => println!("File2 not given"),
    }

    let mut rx = read_and_interleave_fastq(&args.file1, args.file2.as_deref())?;


    while let Some(record) = rx.recv().await {

        match record.id() {
            Ok(id) => println!("ID: {}", id),
            Err(e) => println!("Error getting ID: {}", e),
        }

        println!("Sequence: {}", String::from_utf8_lossy(record.seq()));
        println!("Quality: {}", String::from_utf8_lossy(record.qual()));
        println!("---"); // Separator between records
    }

    Ok(())

}
