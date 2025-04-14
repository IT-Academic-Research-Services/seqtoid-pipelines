use std::cmp::max;
use seq_io::fasta::{Reader as FastaReader, OwnedRecord as FastaOwnedRecord};
use seq_io::fastq::{Reader as FastqReader, OwnedRecord as FastqOwnedRecord};
use std::fs::File;
use std::io::{self, BufReader, Read};
use flate2::read::GzDecoder;
use crate::utils::file::is_gzipped;
use crate::utils::Technology;


/// Defines FASTA and FASTQ as part of a unified FASTX structure.
#[derive(Clone)]
pub enum SequenceRecord {
    Fasta {
        id: String,
        desc: Option<String>,
        seq: Vec<u8>,
    },
    Fastq {
        id: String,
        desc: Option<String>,
        seq: Vec<u8>,
        qual: Vec<u8>,
    },
}

/// Maps id and seq to the correct file type.
impl SequenceRecord {
    pub fn id(&self) -> &str {
        match self {
            SequenceRecord::Fasta { id, .. } => id,
            SequenceRecord::Fastq { id, .. } => id,
        }
    }

    pub fn seq(&self) -> &[u8] {
        match self {
            SequenceRecord::Fasta { seq, .. } => seq,
            SequenceRecord::Fastq { seq, .. } => seq,
        }
    }
}

impl From<FastaOwnedRecord> for SequenceRecord {
    fn from(record: FastaOwnedRecord) -> Self {
        let (id, desc) = parse_header(&record.head, '>');
        SequenceRecord::Fasta {
            id,
            desc,
            seq: record.seq,
        }
    }
}

impl From<FastqOwnedRecord> for SequenceRecord {
    fn from(record: FastqOwnedRecord) -> Self {
        let (id, desc) = parse_header(&record.head, '@');
        SequenceRecord::Fastq {
            id,
            desc,
            seq: record.seq,
            qual: record.qual,
        }
    }
}

/// Parses a FASTX header.
///
///
/// # Arguments
///
/// * `head` - Header line of a FASTX record.
/// * 'prefix' - Leading, defining character of the header. > for FASTA, @ for FASTQ.
///
/// # Returns
/// Tuple: (id, desc) split of header on whitespace.
///
fn parse_header(head: &[u8], prefix: char) -> (String, Option<String>) {
    let head_str = String::from_utf8_lossy(head).into_owned();
    let parts: Vec<&str> = head_str.splitn(2, |c: char| c.is_whitespace()).collect();
    let id = parts[0].trim_start_matches(prefix).to_string();
    let desc = parts.get(1).map(|s| s.to_string()).filter(|s| !s.is_empty());
    (id, desc)
}

/// Custom reader enum for handling compressed/uncompressed files
pub enum FileReader {
    Uncompressed(BufReader<File>),
    Gzipped(GzDecoder<File>),
}

/// Trait implementation of reading from either a compressed or uncompressed file.
impl Read for FileReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            FileReader::Uncompressed(r) => r.read(buf),
            FileReader::Gzipped(r) => r.read(buf),
        }
    }
}

// Enum to hold either FASTA or FASTQ reader
pub enum SequenceReader {
    Fasta(FastaReader<FileReader>),
    Fastq(FastqReader<FileReader>),
}

// Function to create a SequenceReader based on file type
pub fn sequence_reader(path: &str) -> io::Result<SequenceReader> {
    let file = File::open(path)?;
    let is_gz = is_gzipped(path)?;
    let reader = if is_gz {
        FileReader::Gzipped(GzDecoder::new(file))
    } else {
        FileReader::Uncompressed(BufReader::new(file))
    };

    // Determine format by extension (or add explicit format detection if needed)
    let is_fasta = path.to_lowercase().ends_with(".fasta") || path.to_lowercase().ends_with(".fa");
    Ok(if is_fasta {
        SequenceReader::Fasta(FastaReader::new(reader))
    } else {
        SequenceReader::Fastq(FastqReader::new(reader))
    })
}

/// Counts the number of records in a FASTQ.
///
///
/// # Arguments
///
/// * `fastq_path` - Valid path to a fastq file.
///
/// # Returns
/// u64: Number of records in the FASTQ.
///
pub fn record_counter(path: &str) -> io::Result<u64> {
    let mut counter = 0;
    match sequence_reader(path)? {
        SequenceReader::Fasta(reader) => {
            for _ in reader.into_records() {
                counter += 1;
            }
        }
        SequenceReader::Fastq(reader) => {
            for _ in reader.into_records() {
                counter += 1;
            }
        }
    }
    Ok(counter)
}


/// Asynchronously outputs a stream from one or two FASTQ files.
/// If two FASTQ, the stream is interleaved and a header check is performed 
/// to ensure each pair of reads is really an R1/R2 pair.
/// # Arguments
///
/// * `fastq_path` - Valid path to a fastq file.
///
/// # Returns
/// Receiver<Owned Record> stream on an async stream.
///
pub fn read_and_interleave_sequences(
    path1: &str,
    path2: Option<&str>,
    technology: Option<Technology>,
    max_reads: usize,
) -> anyhow::Result<tokio::sync::mpsc::Receiver<SequenceRecord>> {
    let path1 = path1.to_string();
    let path2 = path2.map(String::from);
    let (tx, rx) = tokio::sync::mpsc::channel(100);
    let mut read_counter = 0;
    match (path2, sequence_reader(&path1)?) {
        (Some(path2), SequenceReader::Fastq(_)) => {
            tokio::spawn(async move {
                let reader1 = sequence_reader(&path1).unwrap();
                let reader2 = sequence_reader(&path2).unwrap();
                let (mut records1, mut records2) = match (reader1, reader2) {
                    (SequenceReader::Fastq(r1), SequenceReader::Fastq(r2)) => {
                        (r1.into_records(), r2.into_records())
                    }
                    _ => {
                        eprintln!("Paired-end mode requires FASTQ files");
                        return;
                    }
                };
                while let (Some(Ok(r1)), Some(Ok(r2))) = (records1.next(), records2.next()) {
                    let r1_owned: SequenceRecord = r1.to_owned().into();
                    let r2_owned: SequenceRecord = r2.to_owned().into();
                    if let Some(Technology::illumina) = technology {
                        if !compare_read_ids(r1_owned.id(), r2_owned.id()) {
                            eprintln!("Read ID mismatch in paired-end FASTQ");
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
                    read_counter += 1;
                    if read_counter >= max_reads {
                        break;
                    }
                }
            });
        }
        (None, SequenceReader::Fastq(reader)) => {
            tokio::spawn(async move {
                for result in reader.into_records() {
                    if let Ok(record) = result {
                        if tx.send(record.to_owned().into()).await.is_err() {
                            eprintln!("Failed to send FASTQ record");
                            break;
                        }
                    }
                }
            });
        }
        (None, SequenceReader::Fasta(reader)) => {
            tokio::spawn(async move {
                for result in reader.into_records() {
                    if let Ok(record) = result {
                        if tx.send(record.to_owned().into()).await.is_err() {
                            eprintln!("Failed to send FASTA record");
                            break;
                        }
                    }
                }
            });
        }
        (Some(_), SequenceReader::Fasta(_)) => {
            return Err(anyhow::anyhow!("Paired-end mode not supported for FASTA"));
        }
    }

    Ok(rx)
}


/// Compares the headers ot two FASTQ reads.
/// If two FASTQ, the stream is interleaved and a header check is performed 
/// to ensure each pair of reads is really an R1/R2 pair.
/// # Arguments
///
/// * `id1_result`: &str - ID string for read 1
/// * `id2_result`: &str - ID string for read 2
///
/// # Returns
/// bool: true if reads are a matched pair.
///
fn compare_read_ids(
    id1: &str,
    id2: &str,
) -> bool {
    
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_sequence_reader_fasta() -> io::Result<()> {
        let mut tmp = NamedTempFile::new_in(std::env::temp_dir())?;
        let path = tmp.path().with_extension("fasta");
        std::fs::rename(tmp.path(), &path)?;
        writeln!(tmp, ">seq1 testFASTA\nATCG")?;
        tmp.flush()?;
        
        let reader = sequence_reader(path.to_str().ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Invalid path"))?)?;
        match reader {
            SequenceReader::Fasta(_) => Ok(()),
            _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fasta reader")),
        }
    }

    #[test]
    #[test]
    fn test_sequence_reader_fastq() -> io::Result<()> {
        let mut tmp = NamedTempFile::new_in(std::env::temp_dir())?;
        let path = tmp.path().with_extension("fastq");
        std::fs::rename(tmp.path(), &path)?;
        writeln!(tmp, "@seq1\nATCG\n+\nIIII")?;
        tmp.flush()?;

        let reader = sequence_reader(path.to_str().ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Invalid path"))?)?;
        match reader {
            SequenceReader::Fastq(_) => Ok(()),
            _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fastq reader")),
        }
    }
}