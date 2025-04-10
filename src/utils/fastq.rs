use seq_io::fastq::{Reader, Record, OwnedRecord};
use std::fs::File;
use std::io::{self, BufReader, Read};
use flate2::read::GzDecoder;
use tokio::sync::mpsc;
use crate::utils::file::is_gzipped;

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