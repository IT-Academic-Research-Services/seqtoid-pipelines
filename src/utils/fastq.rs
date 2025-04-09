
use seq_io::fastq::Reader;
use std::fs::File;
use std::io::{self, BufReader, Read};
use flate2::read::GzDecoder;
use crate::utils::file::is_gzipped;

/// Returns a seq_io-based FASTQ reader
pub fn fastq_reader(path: &str) -> io::Result<Reader<Box<dyn Read>>> {
    let file = File::open(path)?;
    let reader: Reader<Box<dyn Read>> = if is_gzipped(path)? {
        let file = File::open(path)?; // Reopen due to cursor movement
        let decoder = GzDecoder::new(file);
        Reader::new(Box::new(decoder) as Box<dyn Read>)
    } else {
        let buffered = BufReader::new(file);
        Reader::new(Box::new(buffered) as Box<dyn Read>)
    };

    Ok(reader)
}

pub fn read_and_interleave_fastq(fastq1_path: &str, fastq2_path: Option<&str>, max_reads: usize) -> Result<bool, std::io::Error> {
    let mut paired_end = false;

    let fastq1_reader = fastq_reader(fastq1_path)?;

    match fastq2_path {
        Some(p) => {
            let fastq2_reader = fastq_reader(p)?;
            paired_end = true;
        }
        None => println!("No file path provided"),
    }
    Ok(true)
}
