
use crate::utils::{Arguments, fastx};
use seq_io::fastq::{Reader, Record};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;
use flate2::read::GzDecoder;
use crate::utils::fastx::{record_counter, read_and_interleave_sequences};
use std::time::Instant;
use crate::utils::fastx::SequenceRecord;


pub async fn run(args: &Arguments) -> anyhow::Result<()> {
    eprintln!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);
    println!("File 1: {}", args.file1);
    match &args.file2 {
        Some(file) => println!("File 2: {}", file),
        None => println!("File2 not given"),
    }
    let technology = Some(args.technology.clone());

    // let start = Instant::now();
    // let fq1count = record_counter(&args.file1);
    // let duration = start.elapsed();
    // println!("Num records: {:?} Time in ms: {}", fq1count,  duration.as_millis());
    
    let mut rx = read_and_interleave_sequences(&args.file1, args.file2.as_deref(), technology)?;
    // 
    // 
    while let Some(record) = rx.recv().await {

        println!("ID: {}", record.id());

        println!("---"); // Separator between records
    }

    Ok(())

}
