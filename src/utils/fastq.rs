use seq_io::fastq::{Reader, Record, OwnedRecord};
use std::fs::File;
use std::io::{self, BufReader, Read};
use flate2::read::GzDecoder;
use tokio::sync::mpsc;
use crate::utils::file::is_gzipped;
use crate::utils::Technology;

pub enum FastqReader {
    Uncompressed(BufReader<File>),
    Gzipped(GzDecoder<File>),
}

impl Read for FastqReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            FastqReader::Uncompressed(r) => r.read(buf),
            FastqReader::Gzipped(r) => r.read(buf),
        }
    }
}

pub fn fastq_reader(path: &str) -> io::Result<Reader<FastqReader>> {
    let file = File::open(path)?;
    let reader = if is_gzipped(path)? {
        Reader::new(FastqReader::Gzipped(GzDecoder::new(file)))
    } else {
        Reader::new(FastqReader::Uncompressed(BufReader::new(file)))
    };
    Ok(reader)
}

pub fn read_and_interleave_fastq(
    fastq1_path: &str,
    fastq2_path: Option<&str>,
    technology: Technology,
) -> anyhow::Result<mpsc::Receiver<OwnedRecord>> {
    let fastq1_path = fastq1_path.to_string();
    let fastq2_path = fastq2_path.map(String::from);

    let (tx, rx) = mpsc::channel(100);

    if let Some(fastq2_path) = fastq2_path { // paired-ended
        tokio::spawn(async move {
            let reader1 = fastq_reader(&fastq1_path).unwrap();
            let reader2 = fastq_reader(&fastq2_path).unwrap();
            let mut records1 = reader1.into_records();
            let mut records2 = reader2.into_records();

            while let (Some(Ok(r1)), Some(Ok(r2))) = (records1.next(), records2.next()) {
                let r1_owned = r1.to_owned();
                let r2_owned = r2.to_owned();


                if technology == Technology::illumina {
                    if !compare_read_ids(r1_owned.id(), r2_owned.id()) {
                        eprintln!("Mismatch!");
                    }
                    
                }
                
                if tx.send(r1_owned).await.is_err() {
                    eprintln!("Failed to send R1 record");
                    break;
                }
                if tx.send(r2_owned).await.is_err() {
                    eprintln!("Failed to send R2 record");
                    break;
                }

            }
        });
    } else { // single-ended
        let tx = tx.clone();
        tokio::spawn(async move {
            if let Ok(reader) = fastq_reader(&fastq1_path) {
                for result in reader.into_records() {
                    if let Ok(record) = result {
                        if tx.send(record.to_owned()).await.is_err() {
                            eprintln!("Failed to send record");
                            break;
                        }
                    }
                }
            }
        });
    }

    Ok(rx)
}

fn compare_read_ids(
    id1_result: Result<&str, impl std::error::Error>,
    id2_result: Result<&str, impl std::error::Error>,
) -> bool {
    match (id1_result, id2_result) {
        (Ok(id1), Ok(id2)) => {
            // Try Casava 1.8+ format first (space-separated)
            let id1_parts: Vec<&str> = id1.splitn(2, ' ').collect();
            let id2_parts: Vec<&str> = id2.splitn(2, ' ').collect();

            if id1_parts.len() == 2 && id2_parts.len() == 2 {
                let read_id1 = id1_parts[0];
                let read_id2 = id2_parts[0];
                if read_id1 != read_id2 {
                    return false;
                }
                let attr1_parts: Vec<&str> = id1_parts[1].split(':').collect();
                let attr2_parts: Vec<&str> = id2_parts[1].split(':').collect();
                if attr1_parts.is_empty() || attr2_parts.is_empty() {
                    eprintln!("Invalid attributes format in R1 or R2");
                    return false;
                }
                let read_num1 = attr1_parts[0];
                let read_num2 = attr2_parts[0];
                return (read_num1 == "1" && read_num2 == "2") || (read_num1 == "2" && read_num2 == "1");
            }

            // Fallback to /1 and /2 format
            if id1.ends_with("/1") && id2.ends_with("/2") {
                let base_id1 = id1.trim_end_matches("/1");
                let base_id2 = id2.trim_end_matches("/2");
                return base_id1 == base_id2;
            } else if id1.ends_with("/2") && id2.ends_with("/1") {
                let base_id1 = id1.trim_end_matches("/2");
                let base_id2 = id2.trim_end_matches("/1");
                return base_id1 == base_id2;
            }

            eprintln!("Unsupported read ID format or mismatched pairs");
            false
        }
        (Err(e1), _) => {
            eprintln!("Error getting R1 ID: {}", e1);
            false
        }
        (_, Err(e2)) => {
            eprintln!("Error getting R2 ID: {}", e2);
            false
        }
    }
}