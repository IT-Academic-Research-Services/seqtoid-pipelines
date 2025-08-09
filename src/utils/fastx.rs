use std::fs::File;
use std::io::Write;
use std::io::{self, BufReader};
use std::path::PathBuf;
use std::sync::Arc;
use std::collections::HashSet;
use anyhow::{Result, anyhow};
use flate2::read::GzDecoder;
use crate::utils::file::{extension_remover, is_gzipped, FileReader};
use crate::cli::Technology;
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::utils::sequence::{DNA, normal_phred_qual_string}; 
use futures::Stream;
use tokio_stream::{self as stream};
use crate::config::defs::{FASTA_TAG, FASTQ_TAG, FASTA_EXTS, FASTQ_EXTS};
use seq_io::fasta::{Reader as FastaReader, OwnedRecord as FastaOwnedRecord};
use seq_io::fastq::{Reader as FastqReader, OwnedRecord as FastqOwnedRecord};
use futures::stream::StreamExt;
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::streams::{ParseOutput};
use tokio::time::Duration;
use futures::future::try_join_all;
use tokio::task::JoinHandle;

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

/// Define size of contig thresholds to mimic what Quast does
lazy_static! {
    pub static ref CONTIG_THRESHOLDS: Vec<usize> = vec![0, 1000, 5000, 10000, 25000, 50000];
}



/// Defines FASTA and FASTQ as part of a unified FASTX structure.
#[derive(Clone, Debug)]
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
    #[allow(dead_code)]
    pub fn qual(&self) -> &[u8] {
        match self {
            SequenceRecord::Fasta { .. } => &[],
            SequenceRecord::Fastq { qual, .. } => qual,
        }
    }

    #[allow(dead_code)]
    pub fn desc(&self) -> Option<&str> {
        match self {
            SequenceRecord::Fasta { desc, .. } => desc.as_deref(),
            SequenceRecord::Fastq { desc, .. } => desc.as_deref(),
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


/// HashSet based sequence validator
///
///
/// # Arguments
///
/// * `seq`: &[u8] ref to vec of bytes
/// * 'valid_bases: vec of allowed bases
///
/// # Returns
/// Result<(), String> for error if any
pub fn validate_sequence(seq: &Arc<Vec<u8>>, valid_bases: &[u8]) -> Result<()> {
    let valid_bases_set: HashSet<u8> = valid_bases.iter().copied().collect();
    if let Some(&invalid_base) = seq.iter().find(|&&b| !valid_bases_set.contains(&b)) {
        Err(anyhow::anyhow!("Invalid nucleotide '{}' found in sequence", invalid_base as char))
    } else {
        Ok(())
    }
}


/// PArallel HashSet based sequence validator
/// for large sequences (probably over 1 billion bases)
///
/// # Arguments
///
/// * `seq`: &[u8] ref to vec of bytes
/// * 'valid_bases: vec of allowed bases
/// * num threads
///
/// # Returns
/// Result<(), String> for error if any
pub async fn validate_sequence_parallel(seq: Arc<Vec<u8>>, valid_bases: &[u8], num_threads: usize) -> Result<()> {
    let valid_bases_set: HashSet<u8> = valid_bases.iter().copied().collect();
    let chunk_size = (seq.len() + num_threads - 1) / num_threads;
    let mut handles: Vec<JoinHandle<Result<(), String>>> = Vec::new();

    for i in 0..num_threads {
        let start = i * chunk_size;
        let end = (start + chunk_size).min(seq.len());
        if start >= seq.len() {
            break;
        }
        let seq_arc = Arc::clone(&seq);
        let valid_bases_set = valid_bases_set.clone();
        let handle = tokio::spawn(async move {
            let chunk = &seq_arc[start..end];
            if let Some(&invalid_base) = chunk.iter().find(|&&b| !valid_bases_set.contains(&b)) {
                Err(format!("Invalid nucleotide '{}' found in sequence", invalid_base as char))
            } else {
                Ok(())
            }
        });
        handles.push(handle);
    }

    let results = try_join_all(handles).await?;
    for result in results {
        result.map_err(|e| anyhow::anyhow!(e))?;
    }
    Ok(())
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
                            if (*R1_R2_TAGS).contains_key(part) {
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
pub fn parse_header(head: &[u8], prefix: char) -> (String, Option<String>) {
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
#[allow(dead_code)]
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


/// Generates FASTQ (not yet FASTA) records
/// # Arguments
///
/// * `num_records` - Number of SequenceRecords to make.
/// * 'seq_len' = Length of each seq.
/// * 'mean' = Mean quality of bases
/// * 'stdev' = St Dev of quality of bases. 
///
/// # Returns
/// Stream<Item = SequenceRecord>
pub fn fastx_generator(num_records: usize, seq_len: usize, mean: f32, stdev: f32) -> impl Stream<Item = SequenceRecord> {
    let records: Vec<SequenceRecord> = if seq_len == 0 {
        Vec::new() // Empty vector for zero read size
    } else {
        (0..num_records)
            .map(|i| {
                let seq = DNA::random_sequence(seq_len);
                let qual = normal_phred_qual_string(seq_len, mean, stdev);
                SequenceRecord::Fastq {
                    id: format!("read{}", i + 1),
                    desc: None,
                    seq: seq.into_bytes(),
                    qual: qual.into_bytes(),
                }
            })
            .collect()
    };
    stream::iter(records)
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
    let (tx, rx) = tokio::sync::mpsc::channel(10000);
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

                let mut r1_count = 0;
                let mut r2_count = 0;
                let mut last_progress = tokio::time::Instant::now();

                while let (Some(r1_result), Some(r2_result)) = (records1.next(), records2.next()) {
                    if last_progress.elapsed() > tokio::time::Duration::from_secs(15) {
                        eprintln!("Stall detected at {} read pairs", read_counter);
                        last_progress = tokio::time::Instant::now();
                    }

                    match (r1_result, r2_result) {
                        (Ok(r1), Ok(r2)) => {
                            let r1_owned: SequenceRecord = r1.to_owned().into();
                            let r2_owned: SequenceRecord = r2.to_owned().into();

                            if let Some(min_len) = min_read_len {
                                if r1_owned.seq().len() < min_len || r2_owned.seq().len() < min_len {
                                    eprintln!("Read length below minimum: {} at pair {}", min_len, read_counter + 1);
                                    return;
                                }
                            }
                            if let Some(max_len) = max_read_len {
                                if r1_owned.seq().len() > max_len || r2_owned.seq().len() > max_len {
                                    eprintln!("Read length above maximum: {} at pair {}", max_len, read_counter + 1);
                                    return;
                                }
                            }

                            if let Some(Technology::Illumina) = technology {
                                if !compare_read_ids(r1_owned.id(), r2_owned.id()) {
                                    eprintln!("Read ID mismatch at R1 count {}: {} vs {}", r1_count + 1, r1_owned.id(), r2_owned.id());
                                    return;
                                }
                            }

                            if tx.send(r1_owned).await.is_err() {
                                eprintln!("Failed to send R1 record at count {}", r1_count + 1);
                                return;
                            }
                            r1_count += 1;

                            if tx.send(r2_owned).await.is_err() {
                                eprintln!("Failed to send R2 record at count {}", r2_count + 1);
                                return;
                            }
                            r2_count += 1;

                            read_counter += 1;
                            if read_counter >= max_reads {
                                eprintln!("Reached max reads: {}", max_reads);
                                return;
                            }
                        }
                        (Err(e), _) => {
                            eprintln!("Error reading R1 at count {}: {}", r1_count + 1, e);
                            return;
                        }
                        (_, Err(e)) => {
                            eprintln!("Error reading R2 at count {}: {}", r2_count + 1, e);
                            return;
                        }
                    }
                }

            });
        }
        (None, SequenceReader::Fastq(reader)) => {
            tokio::spawn(async move {
                let mut last_progress = tokio::time::Instant::now();
                for result in reader.into_records() {
                    if last_progress.elapsed() > tokio::time::Duration::from_secs(15) {
                        eprintln!("Stall detected at {} FASTQ records", read_counter);
                        last_progress = tokio::time::Instant::now();
                    }
                    match result {
                        Ok(record) => {
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
                            if tx.send(r1_owned).await.is_err() {
                                eprintln!("Failed to send FASTQ record at count {}", read_counter + 1);
                                return;
                            }
                            read_counter += 1;
                            if read_counter >= max_reads {
                                eprintln!("Reached max reads: {}", max_reads);
                                return;
                            }
                        }
                        Err(e) => {
                            eprintln!("Error reading FASTQ record at count {}: {}", read_counter + 1, e);
                            return;
                        }
                    }
                }
            });
        }
        (None, SequenceReader::Fasta(reader)) => {
            tokio::spawn(async move {
                let mut last_progress = tokio::time::Instant::now();
                for result in reader.into_records() {
                    if last_progress.elapsed() > tokio::time::Duration::from_secs(15) {
                        eprintln!("Stall detected at {} FASTA records", read_counter);
                        last_progress = tokio::time::Instant::now();
                    }
                    match result {
                        Ok(record) => {
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
                            if tx.send(r1_owned).await.is_err() {
                                eprintln!("Failed to send FASTA record at count {}", read_counter + 1);
                                return;
                            }
                            read_counter += 1;
                            if read_counter >= max_reads {
                                eprintln!("Reached max reads: {}", max_reads);
                                return;
                            }
                        }
                        Err(e) => {
                            eprintln!("Error reading FASTA record at count {}: {}", read_counter + 1, e);
                            return;
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
fn compare_read_ids(id1: &str, id2: &str) -> bool {
    // Extract ID parts without @
    let id_part1 = id1.trim_start_matches('@').splitn(2, ' ').next().unwrap_or("");
    let id_part2 = id2.trim_start_matches('@').splitn(2, ' ').next().unwrap_or("");

    // Check for identical IDs (SRA and some Casava 1.8+ cases)
    if id_part1 == id_part2 {
        return true;
    }

    // Check for /1 and /2 format (pre-Casava 1.8 Illumina and custom formats)
    let full_id_part1 = id1.splitn(2, ' ').next().unwrap_or("");
    let full_id_part2 = id2.splitn(2, ' ').next().unwrap_or("");


    if (full_id_part1.ends_with("/1") && full_id_part2.ends_with("/2")) ||
        (full_id_part1.ends_with("/2") && full_id_part2.ends_with("/1")) {
        let base1 = &full_id_part1[..full_id_part1.len() - 2];
        let base2 = &full_id_part2[..full_id_part2.len() - 2];
        if base1 == base2 {
            return true;
        }
    }

    false
}

/// Writes out a FASTA file to a FIFO pipe.
///
/// # Arguments
///
/// * `fasta_path` - Valid path to a FASTA file.
/// * 'fifo_path` - APath used by named FIFO pipe.
///
/// # Returns
/// Result()
///
pub fn write_fasta_to_fifo(fasta_path: &PathBuf, fifo_path: &PathBuf) -> Result<()> {
    let mut reader = match sequence_reader(fasta_path)? {
        SequenceReader::Fasta(reader) => reader,
        _ => return Err(anyhow!("Input file {} is not a FASTA file", fasta_path.display())),
    };

    let mut fifo_file = std::fs::File::create(fifo_path)?;
    for record_result in reader.into_records() {
        let record = record_result.map_err(|e| anyhow!("Error reading FASTA: {}", e))?;
        let seq_record: SequenceRecord = record.to_owned().into();
        let fasta_line = format!(">{}\n{}\n", seq_record.id(), String::from_utf8_lossy(seq_record.seq()));
        fifo_file.write_all(fasta_line.as_bytes())?;
    }
    fifo_file.flush()?;
    Ok(())
}


/// Filters a stream of FASTQ records based on a predicate applied to the ID line.
///
/// # Arguments
/// * `input_rx` - Receiver of parsed FASTQ records (from `parse_fastq` or similar).
/// * `buffer_size` - Size of the output channel buffer.
/// * `filter_fn` - Function to filter records based on their ID.
///
/// # Returns
/// Tuple of (mpsc of sequence records, task join handle result)
pub fn parse_and_filter_fastq_id(
    input_rx: mpsc::Receiver<ParseOutput>,
    buffer_size: usize,
    filter_fn: impl Fn(&str) -> bool + Send + 'static,
) -> (mpsc::Receiver<SequenceRecord>, tokio::task::JoinHandle<Result<(), anyhow::Error>>) {
    let (filtered_tx, filtered_rx) = mpsc::channel(buffer_size);

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(input_rx);
        let mut count = 0;

        while let Some(item) = stream.next().await {
            if let ParseOutput::Fastq(record) = item {
                if filter_fn(record.id()) {
                    if filtered_tx.send(record).await.is_err() {
                        eprintln!("Failed to send filtered FASTQ record at count {}", count + 1);
                        return Err(anyhow::anyhow!(
                            "Failed to send filtered FASTQ record at count {}",
                            count + 1
                        ));
                    }
                    count += 1;
                }
            } else {
                eprintln!("Unexpected ParseOutput::Bytes at count {}", count + 1);
                continue; // Skip non-FASTQ items
            }

            // Periodic sleep to prevent tight loops
            if count % 1000 == 0 {
                tokio::time::sleep(Duration::from_millis(1)).await;
            }
        }

        eprintln!("Filtered {} FASTQ records from Kraken2 classified stream", count);
        Ok(())
    });

    (filtered_rx, task)
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
                    SequenceReader::Fasta(_reader) => { Ok(()) },
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
                    SequenceReader::Fastq(_reader) => { Ok(()) },
                    _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fastq reader")),
                }
            }
            Err(e) => Err(e),
        }
    }
}