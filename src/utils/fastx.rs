use seq_io::fasta::{Reader as FastaReader, OwnedRecord as FastaOwnedRecord};
use seq_io::fastq::{Reader as FastqReader, OwnedRecord as FastqOwnedRecord};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use crate::utils::file::{extension_remover, is_gzipped, FileReader, WriteToFile};
use crate::utils::Technology;
use std::collections::HashMap;
use lazy_static::lazy_static;

const FASTA_TAG : &str = "fasta";
const FASTQ_TAG : &str = "fastq";
const FASTA_EXTS: &[&'static str] = &["fasta", "fa", "fna", "faa", "ffn", "frn"];
const FASTQ_EXTS: &[&'static str] = &["fastq", "fq"];

lazy_static! {
    static ref R1_R2_TAGS: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        m.insert("R1", "R2");
        m.insert("r1", "r2");
        m.insert("1", "2");
        m.insert("F", "R");
        m.insert("f", "r");
        m.insert("FWD", "REV");
        m.insert("fwd", "rev");
        m.insert("PE1", "PE2");
        m.insert("pe1", "pe2");
        m.insert("READ1", "READ2");
        m.insert("read1", "read2");
        m
    };
}

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


/// Enum to hold either FASTA or FASTQ reader
pub enum SequenceReader {
    Fasta(FastaReader<FileReader>),
    Fastq(FastqReader<FileReader>),
}

/// Creates a SequenceReader for either FASTA or FASTQ files.
///
///
/// # Arguments
///
/// * `path`: &PATHBuf - Valid path to a fastx file.
///
/// # Returns
/// io::Result<SequenceReader>: Result bearing the correct SequenceReader.
///
pub fn sequence_reader(path: &PathBuf) -> io::Result<SequenceReader> {
    let file = File::open(path)?;
    let is_gz = is_gzipped(path)?;
    let reader = if is_gz {
        FileReader::Gzipped(GzDecoder::new(file))
    } else {
        FileReader::Uncompressed(BufReader::new(file))
    };

    let is_fasta = fastx_filetype(path)?;
    match is_fasta.as_str() {
        FASTA_TAG => Ok(SequenceReader::Fasta(FastaReader::new(reader))),
        FASTQ_TAG => Ok(SequenceReader::Fastq(FastqReader::new(reader))),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unsupported file type for path: {:?}", path),
        )),
    }
}


/// Implementation for SequenceRecord
impl WriteToFile for SequenceRecord {
    fn write_to_file<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        match self {
            SequenceRecord::Fastq { id, desc, seq, qual } => {
                write_fastq_record(writer, id, desc.as_deref(), seq, qual)
            }
            SequenceRecord::Fasta { id, desc, seq } => {
                write_fasta_record(writer, id, desc.as_deref(), seq)
            }
        }
    }
}

pub fn write_fasta_record<W: Write>(
    writer: &mut W,
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
) -> io::Result<()> {
    // Write header
    writer.write_all(b">")?;
    writer.write_all(id.as_bytes())?;
    if let Some(desc) = desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    writer.write_all(b"\n")?;

    // Write sequence (with line wrapping, e.g., 80 chars per line)
    for chunk in seq.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

pub fn write_fastq_record<W: Write>(
    writer: &mut W,
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    // Write header
    writer.write_all(b"@")?;
    writer.write_all(id.as_bytes())?;
    if let Some(desc) = desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    writer.write_all(b"\n")?;

    // Write sequence
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;

    // Write separator
    writer.write_all(b"+")?;
    writer.write_all(b"\n")?;

    // Write quality scores
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}


/// Determines if a file path is a FASTA, FASTQ, or neither.
/// Checks extensions, not the body.
///
/// # Arguments
///
/// * `path` - Header line of a FASTX record.
///
/// # Returns
/// Result<String>. Ok fastq or fasta, or err.
///
fn fastx_filetype(path: &PathBuf) -> io::Result<String> {
    let (_, extensions) = extension_remover(path);

    for ext in &extensions {
        if FASTA_EXTS.iter().any(|&e| e.eq_ignore_ascii_case(ext)) {
            return Ok(FASTA_TAG.to_string());
        }

        if FASTQ_EXTS.iter().any(|&e| e.eq_ignore_ascii_case(ext)) {
            return Ok(FASTQ_TAG.to_string());
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "File '{}' has invalid extension(s) '{:?}'. Expected FASTA ({:?}) or FASTQ ({:?}).",
            path.display(),
            extensions,
            FASTA_EXTS,
            FASTQ_EXTS
        ),
    ))
}


#[derive(Debug, PartialEq)]
pub struct R1R2Result {
    pub delimiter: Option<char>,
    pub r1_tag: Option<String>,
    pub file_name: Option<String>,
    pub index: Option<usize>,
    
}

/// Tries to locate
/// Checks extensions, not the body.
///
/// # Arguments
///
/// * `path` - Header line of a FASTX record.
///
/// # Returns
/// Result<String>. Ok fastq or fasta, or err.
///
pub fn r1r2_base(path: &PathBuf)  -> R1R2Result {

    let delimiters = ['_', '.', '-'];
    let (stem, extensions) = extension_remover(&path);
    
    match stem.file_name() {
        Some(name) => {
            match name.to_str() {
                Some(filename) => {
                    for &delimiter in delimiters.iter() {
                        let parts: Vec<&str> = filename.split(delimiter).collect();
                        for (index, part) in parts.iter().enumerate() {
                            if R1_R2_TAGS.contains_key(part) {
                                let r1_tag = part.to_string();
                                let prefix_parts = &parts[..index];
                                let prefix = if prefix_parts.is_empty() {
                                    String::new()
                                } else {
                                    prefix_parts.join(&delimiter.to_string())
                                };

                                let new_file = format!("{}.{}", prefix, extensions.join("."));
                                
                                return R1R2Result {
                                    delimiter: Some(delimiter),
                                    r1_tag: Some(r1_tag),
                                    file_name: Some(new_file),
                                    index: Some(index),
                                };
                            }
                        }
                    }
                }
                None => {
                    return R1R2Result {
                        delimiter: None,
                        r1_tag: None,
                        file_name: None,
                        index: None,
                    };
                }
            }
        },
        None => {
            return R1R2Result {
                delimiter: None,
                r1_tag: None,
                file_name: None,
                index: None,
            };
        }
    }
    
    return R1R2Result {
        delimiter: None,
        r1_tag: None,
        file_name: None,
        index: None,
    };
    
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
pub fn record_counter(path: &PathBuf) -> io::Result<u64> {
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

/// Enum to hold either an uncompressed or gzipped file writer
pub enum FileWriter {
    Uncompressed(BufWriter<File>),
    Gzipped(GzEncoder<BufWriter<File>>),
}

impl Write for FileWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            FileWriter::Uncompressed(w) => w.write(buf),
            FileWriter::Gzipped(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            FileWriter::Uncompressed(w) => w.flush(),
            FileWriter::Gzipped(w) => w.flush(),
        }
    }
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
    path1: PathBuf,
    path2: Option<PathBuf>,
    technology: Option<Technology>,
    max_reads: usize,
    min_read_len: Option<usize>,
    max_read_len: Option<usize>,
) -> anyhow::Result<tokio::sync::mpsc::Receiver<SequenceRecord>> {

    
    let (tx, rx) = tokio::sync::mpsc::channel(100);
    let mut read_counter = 0;
    match (path2, sequence_reader(&path1)?) {
        (Some(path2), SequenceReader::Fastq(_)) => {

            if let Some(Technology::ONT) = technology {
                return Err(anyhow::anyhow!("Paired-end mode not supported for ONT technology!"));
            }

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

                    if let Some(min_len) = min_read_len {
                        if r1_owned.seq().len() < min_len || r2_owned.seq().len() < min_len {
                            eprintln!("Read length below minimum: {}", min_len);
                            return; // Drop tx to close channel
                        }
                    }
                    if let Some(max_len) = max_read_len {
                        if r1_owned.seq().len() > max_len || r2_owned.seq().len() > max_len {
                            eprintln!("Read length above maximum: {}", max_len);
                            return;
                        }
                    }
                    
                    if let Some(Technology::Illumina) = technology {
                        if !compare_read_ids(r1_owned.id(), r2_owned.id()) {
                            eprintln!("Read ID mismatch in paired-end FASTQ");
                            return;
                        }
                    }
                    if tx.send(r1_owned).await.is_err() {
                        eprintln!("Failed to send R1 record");
                        return;
                    }
                    if tx.send(r2_owned).await.is_err() {
                        eprintln!("Failed to send R2 record");
                        return;
                    }
                    read_counter += 1;
                    if read_counter >= max_reads {
                        return;
                    }
                }
            });
        }
        (None, SequenceReader::Fastq(reader)) => {
            tokio::spawn(async move {
                for result in reader.into_records() {
                    if let Ok(record) = result {
                        let r1_owned: SequenceRecord = record.to_owned().into();

                        if let Some(min_len) = min_read_len {
                            if r1_owned.seq().len() < min_len {
                                eprintln!("Read length below minimum: {}", min_len);
                                return;
                            }
                        }
                        if let Some(max_len) = max_read_len {
                            if r1_owned.seq().len() > max_len {
                                eprintln!("Read length above maximum: {}", max_len);
                                return; 
                            }
                        }
                        
                        if tx.send(r1_owned.into()).await.is_err() {
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
                        let r1_owned: SequenceRecord = record.to_owned().into();

                        if let Some(min_len) = min_read_len {
                            if r1_owned.seq().len() < min_len {
                                eprintln!("Read length below minimum: {}", min_len);
                                return;
                            }
                        }
                        if let Some(max_len) = max_read_len {
                            if r1_owned.seq().len() > max_len {
                                eprintln!("Read length above maximum: {}", max_len);
                                return;
                            }
                        }
                        if tx.send(r1_owned.to_owned().into()).await.is_err() {
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
        
        let reader_result = sequence_reader(&path);
        match reader_result {
            Ok(reader) => {
                match reader {
                    SequenceReader::Fasta(reader) => { Ok(()) },
                    _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fasta reader")),
                }
            }
            Err(e) => Err(e),
        }
        
    }
    
    #[test]
    fn test_sequence_reader_fastq() -> io::Result<()> {
        let mut tmp = NamedTempFile::new_in(std::env::temp_dir())?;
        let path = tmp.path().with_extension("fastq");
        std::fs::rename(tmp.path(), &path)?;
        writeln!(tmp, "@seq1\nATCG\n+\nIIII")?;
        tmp.flush()?;

        let reader_result = sequence_reader(&path);
        match reader_result {
            Ok(reader) => {
                match reader {
                    SequenceReader::Fastq(reader) => { Ok(()) },
                    _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fastq reader")),
                }
            }
            Err(e) => Err(e),
        }
    }
}