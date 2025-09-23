use std::fs::File;
use std::io::Write;
use std::io::Cursor;
use std::io::{self, BufReader};
use std::path::PathBuf;
use std::sync::Arc;
use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use anyhow::{Result, anyhow};
use flate2::read::GzDecoder;
use crate::utils::file::{extension_remover, is_gzipped, FileReader};
use crate::cli::Technology;
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::utils::sequence::{DNA, normal_phred_qual_string}; 
use futures::Stream;
use tokio_stream::{self as stream};
use crate::config::defs::{FASTA_TAG, FASTQ_TAG, FASTA_EXTS, FASTQ_EXTS, ReadStats};
use needletail::{parse_fastx_file, FastxReader, parser::{SequenceRecord as NeedletailSequenceRecord, FastqReader}};
use futures::stream::StreamExt;
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::streams::ParseOutput;
use tokio::fs::File as TokioFile;
use tokio::io::{AsyncWriteExt, BufWriter};
use tokio::process::Command as TokioCommand;
use tokio::sync::mpsc::Receiver;
use tokio::time::{Duration, Instant};
use futures::future::try_join_all;
use tokio::task::JoinHandle;
use memchr::memmem;
use memchr::memchr;

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

lazy_static! {
    pub static ref CONTIG_THRESHOLDS: Vec<usize> = vec![0, 1000, 5000, 10000, 25000, 50000];
}



/// Defines FASTA and FASTQ as part of a unified FASTX structure.
#[derive(Clone, Debug)]
pub enum SequenceRecord {
    Fasta {
        id: String,
        desc: Option<String>,
        seq: Arc<Vec<u8>>,
    },
    Fastq {
        id: String,
        desc: Option<String>,
        seq: Arc<Vec<u8>>,
        qual: Arc<Vec<u8>>,
    },
}

impl SequenceRecord {
    pub fn id(&self) -> &str {
        match self {
            SequenceRecord::Fasta { id, .. } => id,
            SequenceRecord::Fastq { id, .. } => id,
        }
    }

    pub fn seq(&self) -> &[u8] {
        match self {
            SequenceRecord::Fasta { seq, .. } => &**seq,
            SequenceRecord::Fastq { seq, .. } => &**seq,
        }
    }

    #[allow(dead_code)]
    pub fn qual(&self) -> &[u8] {
        match self {
            SequenceRecord::Fasta { .. } => &[],
            SequenceRecord::Fastq { qual, .. } => &**qual,
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


/// Enum to hold either FASTA or FASTQ reader
pub enum SequenceReader {
    Fasta(Box<dyn needletail::FastxReader>),
    Fastq(Box<dyn needletail::FastxReader>),
}

impl From<needletail::parser::SequenceRecord<'_>> for SequenceRecord {
    fn from(record: needletail::parser::SequenceRecord<'_>) -> Self {
        let (id, desc) = parse_header(record.id(), if record.qual().is_some() { '@' } else { '>' });
        if let Some(qual) = record.qual() {
            SequenceRecord::Fastq {
                id,
                desc,
                seq: Arc::new(record.seq().to_vec()),
                qual: Arc::new(qual.to_vec()),
            }
        } else {
            SequenceRecord::Fasta {
                id,
                desc,
                seq: Arc::new(record.seq().to_vec()),
            }
        }
    }
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
    let reader = needletail::parse_fastx_file(path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Failed to open FASTX file {}: {}", path.display(), e)))?;

    let is_fasta = fastx_filetype(path)?;
    match is_fasta.as_str() {
        FASTA_TAG => Ok(SequenceReader::Fasta(reader)),
        FASTQ_TAG => Ok(SequenceReader::Fastq(reader)),
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
    pub prefix: Option<String>,
    
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
                                    prefix: Some(prefix),
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
                        prefix: None,
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
                prefix: None,
            };
        }
    }
    
    return R1R2Result {
        delimiter: None,
        r1_tag: None,
        file_name: None,
        index: None,
        prefix: None,
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


/// Reads a FASTQ file.
///
/// # Arguments
///
/// * `path1` - Valid path to a FASTQ file.
/// * `path2` - Optional, for R2 file of paired end reads..
/// * 'technology' - Illumina or Nanopore/
/// * 'max_reads' - Cutoff to allow truncated file reading.
/// * 'min_seq_length' - Reads under this value are discarded.
/// * 'max_seq_length' - Reads over this value are discarded.
/// * 'chunk_size' - Size of segments of input data to be processed.
/// # Returns
///  Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<u64, anyhow::Error>>)> : Byte stream of FASTQ input and task handle with validated read count.
///
pub fn read_fastq(
    path1: PathBuf,
    path2: Option<PathBuf>,
    technology: Option<Technology>,
    max_reads: u64,
    min_read_len: Option<usize>,
    max_read_len: Option<usize>,
    chunk_size: usize,
) -> Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<ReadStats, anyhow::Error>>)> {
    let (tx, rx) = mpsc::channel(chunk_size / 1024); // ~1KB per item

    let task = tokio::spawn(async move {
        eprintln!("read_fastq: Starting with path1={:?}, path2={:?}, max_reads={}, min_read_len={:?}, max_read_len={:?}, chunk_size={}",
                  path1, path2, max_reads, min_read_len, max_read_len, chunk_size);

        let mut undersized_reads: u64 = 0;
        let mut oversized_reads: u64 = 0;
        let mut validated_reads: u64 = 0;

        match path2 {
            None => {
                // Single-end: Stream chunks from one file
                eprintln!("read_fastq: Opening single-end FASTQ: {:?}", path1);
                let mut reader = parse_fastx_file(&path1)
                    .map_err(|e| {
                        anyhow!("Failed to open FASTQ {}: {}", path1.display(), e)
                    })?;
                let mut buffer = Vec::with_capacity(chunk_size);

                while let Some(result) = reader.next() {
                    if validated_reads + undersized_reads + oversized_reads >= max_reads {
                        eprintln!("read_fastq: Reached max_reads={}", max_reads);
                        break;
                    }
                    let record = result.map_err(|e| {
                        anyhow!("Parse error at read {}: {}", validated_reads + undersized_reads + oversized_reads + 1, e)
                    })?;

                    // Ensure FASTQ (not FASTA)
                    if record.qual().is_none() {
                        return Err(anyhow!("Expected FASTQ, got FASTA at read {}", validated_reads + undersized_reads + oversized_reads + 1));
                    }

                    // Log read details
                    let seq_len = record.seq().len();
                    let _id_str = std::str::from_utf8(record.id()).unwrap_or("<invalid utf8>");

                    // Optional length filters
                    if let Some(min) = min_read_len {
                        if seq_len < min {
                            eprintln!("read_fastq: Skipping read {}: seq_len={} < min={}", validated_reads + undersized_reads + oversized_reads + 1, seq_len, min);
                            undersized_reads += 1;
                            continue;
                        }
                    }
                    if let Some(max) = max_read_len {
                        if seq_len > max {
                            eprintln!("read_fastq: Skipping read {}: seq_len={} > max={}", validated_reads + undersized_reads + oversized_reads + 1, seq_len, max);
                            oversized_reads += 1;
                            continue;
                        }
                    }

                    // Construct FASTQ record
                    let mut record_bytes = Vec::new();
                    record_bytes.push(b'@');
                    record_bytes.extend_from_slice(record.id());
                    record_bytes.push(b'\n');
                    record_bytes.extend_from_slice(&record.seq());
                    record_bytes.push(b'\n');
                    record_bytes.push(b'+');
                    record_bytes.push(b'\n');
                    record_bytes.extend_from_slice(record.qual().unwrap());
                    record_bytes.push(b'\n');

                    buffer.extend_from_slice(&record_bytes);
                    validated_reads += 1;

                    // Send if chunk full
                    if buffer.len() >= chunk_size {
                        if tx.send(ParseOutput::Bytes(Arc::new(std::mem::take(&mut buffer)))).await.is_err() {
                            return Err(anyhow!("Failed to send byte chunk at read {}", validated_reads));
                        }
                    }
                }
                // Send remaining
                if !buffer.is_empty() {
                    if tx.send(ParseOutput::Bytes(Arc::new(buffer))).await.is_err() {
                        return Err(anyhow!("Failed to send final byte chunk"));
                    }
                }
                eprintln!("read_fastq: Processed {} single-end reads (undersized: {}, validated: {}, oversized: {})",
                          validated_reads + undersized_reads + oversized_reads, undersized_reads, validated_reads, oversized_reads);
                if validated_reads == 0 {
                    eprintln!("read_fastq: Warning: No reads processed from {:?}", path1);
                }
                Ok(ReadStats {
                    undersized: undersized_reads,
                    validated: validated_reads,
                    oversized: oversized_reads,
                })
            }
            Some(path2) => {
                if let Some(Technology::ONT) = technology {
                    return Err(anyhow!("Paired-end not supported for ONT"));
                }
                // Paired-end: Interleave raw record bytes
                eprintln!("read_fastq: Opening paired-end FASTQ: R1={:?}, R2={:?}", path1, path2);
                let mut reader1 = parse_fastx_file(&path1)
                    .map_err(|e| {
                        anyhow!("Failed to open R1 FASTQ {}: {}", path1.display(), e)
                    })?;
                let mut reader2 = parse_fastx_file(&path2)
                    .map_err(|e| {
                        anyhow!("Failed to open R2 FASTQ {}: {}", path2.display(), e)
                    })?;
                let mut buffer = Vec::with_capacity(chunk_size);

                loop {
                    if validated_reads + undersized_reads + oversized_reads >= max_reads {
                        eprintln!("read_fastq: Reached max_reads={}", max_reads);
                        break;
                    }
                    // Read R1
                    let r1_opt = reader1.next();
                    let r1_record = match r1_opt {
                        Some(Ok(rec)) => rec,
                        Some(Err(e)) => {
                            return Err(anyhow!("R1 parse error at read {}: {}", validated_reads + undersized_reads + oversized_reads + 1, e));
                        }
                        None => {
                            eprintln!("read_fastq: R1 stream ended at read {}", validated_reads + undersized_reads + oversized_reads + 1);
                            break;
                        }
                    };
                    if r1_record.qual().is_none() {
                        return Err(anyhow!("Expected FASTQ for R1 at read {}", validated_reads + undersized_reads + oversized_reads + 1));
                    }

                    // Read R2
                    let r2_opt = reader2.next();
                    let r2_record = match r2_opt {
                        Some(Ok(rec)) => rec,
                        Some(Err(e)) => {
                            return Err(anyhow!("R2 parse error at read {}: {}", validated_reads + undersized_reads + oversized_reads + 1, e));
                        }
                        None => {
                            return Err(anyhow!("R2 ended before R1 at read {}", validated_reads + undersized_reads + oversized_reads + 1));
                        }
                    };
                    if r2_record.qual().is_none() {
                        return Err(anyhow!("Expected FASTQ for R2 at read {}", validated_reads + undersized_reads + oversized_reads + 1));
                    }

                    // Log read details
                    let r1_id = std::str::from_utf8(r1_record.id()).unwrap_or("<invalid utf8>");
                    let r2_id = std::str::from_utf8(r2_record.id()).unwrap_or("<invalid utf8>");
                    let r1_len = r1_record.seq().len();
                    let r2_len = r2_record.seq().len();

                    // Relaxed ID check: Strip common suffixes
                    let r1_head = r1_record.id();
                    let r2_head = r2_record.id();
                    if !compare_read_ids_bytes(r1_head, r2_head) {
                        return Err(anyhow!("ID mismatch at pair {}: {:?} vs {:?}", validated_reads + undersized_reads + oversized_reads + 1, r1_id, r2_id));
                    }

                    // Length filters
                    if let Some(min) = min_read_len {
                        if r1_len < min || r2_len < min {
                            undersized_reads += 1;
                            continue;
                        }
                    }
                    if let Some(max) = max_read_len {
                        if r1_len > max || r2_len > max {
                            oversized_reads += 1;
                            continue;
                        }
                    }

                    // Construct R1 FASTQ
                    let mut r1_bytes = Vec::new();
                    r1_bytes.push(b'@');
                    r1_bytes.extend_from_slice(r1_record.id());
                    r1_bytes.push(b'\n');
                    r1_bytes.extend_from_slice(&r1_record.seq());
                    r1_bytes.push(b'\n');
                    r1_bytes.push(b'+');
                    r1_bytes.push(b'\n');
                    r1_bytes.extend_from_slice(r1_record.qual().unwrap());
                    r1_bytes.push(b'\n');

                    // Construct R2 FASTQ
                    let mut r2_bytes = Vec::new();
                    r2_bytes.push(b'@');
                    r2_bytes.extend_from_slice(r2_record.id());
                    r2_bytes.push(b'\n');
                    r2_bytes.extend_from_slice(&r2_record.seq());
                    r2_bytes.push(b'\n');
                    r2_bytes.push(b'+');
                    r2_bytes.push(b'\n');
                    r2_bytes.extend_from_slice(r2_record.qual().unwrap());
                    r2_bytes.push(b'\n');

                    // Interleave: R1 then R2
                    buffer.extend_from_slice(&r1_bytes);
                    buffer.extend_from_slice(&r2_bytes);
                    validated_reads += 1;

                    // Send if chunk full
                    if buffer.len() >= chunk_size {
                        if tx.send(ParseOutput::Bytes(Arc::new(std::mem::take(&mut buffer)))).await.is_err() {
                            return Err(anyhow!("Failed to send byte chunk at read {}", validated_reads));
                        }
                    }
                }
                // Send remaining
                if !buffer.is_empty() {
                    if tx.send(ParseOutput::Bytes(Arc::new(buffer))).await.is_err() {
                        return Err(anyhow!("Failed to send final byte chunk"));
                    }
                }
                // Check for extra R2
                if reader2.next().is_some() {
                    return Err(anyhow!("R2 has extra reads after R1 ended"));
                }
                eprintln!("read_fastq: Processed {} paired-end reads (undersized: {}, validated: {}, oversized: {})",
                          validated_reads + undersized_reads + oversized_reads, undersized_reads, validated_reads, oversized_reads);
                if validated_reads == 0 {
                    eprintln!("read_fastq: Warning: No reads processed from R1={:?}, R2={:?}", path1, path2);
                }
                Ok(ReadStats {
                    undersized: undersized_reads,
                    validated: validated_reads,
                    oversized: oversized_reads,
                })
            }
        }
    });

    Ok((rx, task))
}


/// Reads a FASTA file.
///
///
/// # Arguments
///
/// * `path` - Valid path to a FASTA file.
/// * 'max_records' - Cutoff to allow truncated file reading.
/// * 'min_seq_length' - Reads under this value are discarded.
/// * 'max_seq_length' - Reads over this value are discarded.
/// * 'chunk_size' - Size of segments of input data to be processed.
/// # Returns
///  Result<mpsc::Receiver<ParseOutput>>: Byte stream of FASTA input.
pub fn read_fasta(
    path: PathBuf,
    max_records: u64,
    min_seq_len: Option<usize>,
    max_seq_len: Option<usize>,
    chunk_size: usize,
) -> Result<mpsc::Receiver<ParseOutput>> {
    let (tx, rx) = mpsc::channel(chunk_size / 1024); // ~1KB per item

    tokio::spawn(async move {
        eprintln!("read_fasta: Starting with path={:?}, max_records={}, min_seq_len={:?}, max_seq_len={:?}, chunk_size={}",
                  path, max_records, min_seq_len, max_seq_len, chunk_size);

        let mut reader = parse_fastx_file(&path)
            .map_err(|e| {
                anyhow!("Failed to open FASTA {}: {}", path.display(), e)
            })?;

        let mut buffer = Vec::with_capacity(chunk_size);
        let mut record_count: u64 = 0;

        while let Some(result) = reader.next() {
            if record_count >= max_records {
                eprintln!("read_fasta: Reached max_records={}", max_records);
                break;
            }
            let record = result.map_err(|e| {
                anyhow!("Parse error at record {}: {}", record_count + 1, e)
            })?;

            // Ensure FASTA (not FASTQ)
            if record.qual().is_some() {
                return Err(anyhow!("Expected FASTA, got FASTQ at record {}", record_count + 1));
            }

            // Log record details
            let seq_len = record.seq().len();
            let _id_str = std::str::from_utf8(record.id()).unwrap_or("<invalid utf8>");

            // Optional length filters
            if let Some(min) = min_seq_len {
                if seq_len < min {
                    eprintln!("read_fasta: Skipping record {}: seq_len={} < min={}", record_count + 1, seq_len, min);
                    continue;
                }
            }
            if let Some(max) = max_seq_len {
                if seq_len > max {
                    eprintln!("read_fasta: Skipping record {}: seq_len={} > max={}", record_count + 1, seq_len, max);
                    continue;
                }
            }

            // Construct FASTA record
            let mut record_bytes = Vec::new();
            record_bytes.push(b'>');
            record_bytes.extend_from_slice(record.id());
            record_bytes.push(b'\n');
            record_bytes.extend_from_slice(&record.seq());
            record_bytes.push(b'\n');

            buffer.extend_from_slice(&record_bytes);
            record_count += 1;

            // Send if chunk full
            if buffer.len() >= chunk_size {
                if tx.send(ParseOutput::Bytes(Arc::new(std::mem::take(&mut buffer)))).await.is_err() {
                    return Err(anyhow!("Failed to send byte chunk at record {}", record_count));
                }
            }
        }

        // Send remaining
        if !buffer.is_empty() {
            if tx.send(ParseOutput::Bytes(Arc::new(buffer))).await.is_err() {
                return Err(anyhow!("Failed to send final byte chunk"));
            }
        }

        // eprintln!("read_fasta: Processed {} FASTA records", record_count);
        if record_count == 0 {
            eprintln!("read_fasta: Warning: No records processed from {:?}", path);
        }
        Ok(())
    });

    Ok(rx)
}

fn compare_read_ids_bytes(head1: &[u8], head2: &[u8]) -> bool {
    let id1 = extract_id(head1);
    let id2 = extract_id(head2);
    let id1_str = std::str::from_utf8(id1).unwrap_or("");
    let id2_str = std::str::from_utf8(id2).unwrap_or("");
    let id1_base = id1_str.trim_end_matches(|c| c == '/' || c == ' ' || c == '1' || c == '2');
    let id2_base = id2_str.trim_end_matches(|c| c == '/' || c == ' ' || c == '1' || c == '2');
    // eprintln!("read_fastq: Comparing IDs: {} vs {}", id1_base, id2_base);
    id1_base == id2_base
}

fn extract_id(head: &[u8]) -> &[u8] {
    let start = head.iter().position(|&b| b != b'@').unwrap_or(0);
    let end = head.iter().position(|&b| b == b' ' || b == b'\t').unwrap_or(head.len());
    &head[start..end]
}

/// Counts the number of records in a FASTQ.
///
///
/// # Arguments
///
/// * `path1` - R1 or single ended
/// * `path2` - R2 optional
///
/// # Returns
/// Result<u64>: Number of records in the FASTQ.
///
pub fn raw_read_count(
    path1: PathBuf,
    path2: Option<PathBuf>,
) -> JoinHandle<Result<u64, anyhow::Error>> {
    tokio::spawn(async move {
        let count1 = count_fastx_sequences(&path1).await?;

        if let Some(p2) = path2 {
            let count2 = count_fastx_sequences(&p2).await?;
            if count1 != count2 {
                return Err(anyhow!(
                    "Mismatched read counts in paired files: {} (R1) vs {} (R2)",
                    count1, count2
                ));
            }
            Ok(count1 * 2)
        } else {
            Ok(count1)
        }
    })
}

/// Helper function to count sequences in a single FASTX file using needletail
///
///
/// # Arguments
///
/// * `path1` - R1 or single ended
/// * `path2` - R2 optional
///
/// # Returns
/// Result<u64>: Number of records in the FASTQ.
///
async fn count_fastx_sequences(path: &PathBuf) -> Result<u64> {
    let path_clone = path.clone();
    let handle: JoinHandle<Result<u64>> = tokio::spawn(async move {
        let mut reader = parse_fastx_file(&path_clone)
            .map_err(|e| anyhow!("Failed to open FASTX file {}: {}", path_clone.display(), e))?;
        let mut count: u64 = 0;
        while let Some(record) = reader.next() {
            record.map_err(|e| anyhow!("Parse error at record {}: {}", count + 1, e))?;
            count += 1;
        }
        Ok(count)
    });
    handle.await?
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
                    seq: Arc::new(seq.into_bytes()),
                    qual: Arc::new(qual.into_bytes()),
                }
            })
            .collect()
    };
    stream::iter(records)
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

    let mut fifo_file = File::create(fifo_path)?;
    while let Some(record_result) = reader.next() {
        let record = record_result.map_err(|e| anyhow!("Error reading FASTA: {}", e))?;
        let seq_record: SequenceRecord = record.to_owned().into();
        let fasta_line = format!(">{}\n{}\n", seq_record.id(), String::from_utf8_lossy(seq_record.seq()));
        fifo_file.write_all(fasta_line.as_bytes())?;
    }
    fifo_file.flush()?;
    Ok(())
}

#[allow(dead_code)]
pub async fn stream_record_counter(
    rx: mpsc::Receiver<ParseOutput>,
    early_exit: bool,
) -> Result<u64, anyhow::Error> {
    let mut stream = ReceiverStream::new(rx);
    let mut counter: u64 = 0;

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(chunk) => {
                let chunk_str = String::from_utf8_lossy(&chunk);
                for line in chunk_str.lines() {
                    if line.starts_with('@') || line.starts_with('>') {
                        counter += 1;
                        if early_exit {
                            return Ok(1); // Early exit: Not empty
                        }
                    }
                }
            }
            ParseOutput::Fastq(_) => {
                counter += 1;
                if early_exit {
                    return Ok(1);
                }
            }
            ParseOutput::Fasta(_) => {
                counter += 1;
                if early_exit {
                    return Ok(1);
                }
            }
        }
    }

    Ok(counter)
}

/// Filters a stream of FASTQ records based on a predicate applied to the ID line.
///
/// # Arguments
/// * `input_rx` - Receiver of parsed FASTQ records (from `parse_fastq` or similar).
/// * `buffer_size` - Size of the output channel buffer.
/// * `pattern` - // e.g., "kraken:taxid|123"
///
/// # Returns
/// Tuple of (mpsc of sequence records, task join handle result)
pub fn parse_and_filter_fastq_id(
    input_rx: mpsc::Receiver<ParseOutput>,
    buffer_size: usize,
    pattern: String,
) -> (mpsc::Receiver<SequenceRecord>, JoinHandle<Result<(), anyhow::Error>>) {
    let pattern_bytes = pattern.into_bytes();
    let (filtered_tx, filtered_rx) = mpsc::channel(buffer_size);
    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(input_rx);
        let mut count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Fastq(record) => {
                    if memmem::find(record.id().as_bytes(), &pattern_bytes).is_some() {
                        if filtered_tx.send(record).await.is_err() {
                            eprintln!("Failed to send filtered FASTQ record at count {}", count + 1);
                            return Err(anyhow!(
                                "Failed to send filtered FASTQ record at count {}",
                                count + 1
                            ));
                        }
                        count += 1;
                    }
                }
                _ => {
                    eprintln!("Unexpected ParseOutput at count {}", count + 1);
                    continue; // Skip non-FASTQ items
                }
            }
        }
        eprintln!("Filtered {} FASTQ records from Kraken2 classified stream", count);
        Ok(())
    });
    (filtered_rx, task)
}


pub async fn parse_byte_stream_to_fastq(
    input_rx: mpsc::Receiver<ParseOutput>,
    buffer_size: usize,
    stall_threshold_secs: u64,
) -> Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>)> {
    let (tx, rx) = mpsc::channel(buffer_size);

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(input_rx);
        let mut full_buffer: Vec<u8> = Vec::new(); // Accumulate all bytes (RAM-efficient for filtered streams)
        let mut record_count = 0;
        let mut last_progress = tokio::time::Instant::now();

        while let Some(item) = stream.next().await {
            if last_progress.elapsed() > Duration::from_secs(stall_threshold_secs) {
                eprintln!("parse_byte_stream_to_fastq: Stall detected at {} records", record_count);
                last_progress = tokio::time::Instant::now();
            }

            match item {
                ParseOutput::Bytes(bytes) => {
                    // eprintln!("parse_byte_stream_to_fastq: Received {} bytes chunk", bytes.len());
                    full_buffer.extend_from_slice(&*bytes);
                }
                _ => {
                    eprintln!("parse_byte_stream_to_fastq: Unexpected non-Bytes ParseOutput at record {}", record_count + 1);
                    continue;
                }
            }
        }

        if full_buffer.is_empty() {
            eprintln!("parse_byte_stream_to_fastq: Empty byte buffer received");
            return Err(anyhow!("Empty byte buffer for FASTQ parsing"));
        }

        // eprintln!("parse_byte_stream_to_fastq: Accumulated {} total bytes for parsing", full_buffer.len());

        // Parse accumulated bytes into FASTQ records
        let cursor = Cursor::new(full_buffer);
        let mut reader = FastqReader::new(cursor);
        while let Some(result) = reader.next() {
            match result {
                Ok(record) => {
                    let owned_record: SequenceRecord = record.to_owned().into();
                    if tx.send(ParseOutput::Fastq(owned_record)).await.is_err() {
                        eprintln!("parse_byte_stream_to_fastq: Failed to send FASTQ record at count {}", record_count + 1);
                        return Err(anyhow!("Failed to send FASTQ record at count {}", record_count + 1));
                    }
                    record_count += 1;
                }
                Err(e) => {
                    eprintln!("parse_byte_stream_to_fastq: Error parsing FASTQ at count {}: {}", record_count + 1, e);
                    return Err(anyhow!("FASTQ parsing error at count {}: {}", record_count + 1, e));
                }
            }
        }

        eprintln!("parse_byte_stream_to_fastq: Parsed {} FASTQ records", record_count);
        Ok(())
    });

    Ok((rx, task))
}

/// Concatenates paired-end FASTQ records by joining R1 and R2 sequences with an 'N' and
/// corresponding quality scores with a '!' (low-quality score).
///
/// # Arguments
/// * `input_stream` - ReceiverStream of ParseOutput containing interleaved FASTQ records (R1, R2, ...).
/// * `buffer_size` - Size of the output channel buffer.
/// * `stall_threshold_secs` - Seconds before logging a stall warning.
///
/// # Returns
/// A tuple of:
/// - `ReceiverStream<ParseOutput>`: Stream of concatenated FASTQ records.
/// - `JoinHandle<Result<(), anyhow::Error>>`: Task handle for the concatenation process.
pub async fn concatenate_paired_reads(
    input_stream: ReceiverStream<ParseOutput>,
    buffer_size: usize,
    stall_threshold_secs: u64,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        JoinHandle<Result<(), anyhow::Error>>,
    ),
    anyhow::Error,
> {
    let (tx, rx) = mpsc::channel(buffer_size);
    let mut stream = input_stream;
    let mut last_progress = tokio::time::Instant::now();
    let mut pair_count = 0;

    let task = tokio::spawn(async move {
        let mut r1: Option<SequenceRecord> = None;

        while let Some(item) = stream.next().await {
            if last_progress.elapsed() > Duration::from_secs(stall_threshold_secs) {
                eprintln!("concatenate_paired_reads: Stall detected at {} read pairs", pair_count);
                last_progress = tokio::time::Instant::now();
            }

            match item {
                ParseOutput::Fastq(record) => {
                    if let Some(prev_r1) = r1.take() {
                        // This is R2; process the pair
                        if !compare_read_ids(prev_r1.id(), record.id()) {
                            eprintln!(
                                "Read ID mismatch at pair {}: {} vs {}",
                                pair_count + 1,
                                prev_r1.id(),
                                record.id()
                            );
                            return Err(anyhow!(
                                "Read ID mismatch at pair {}: {} vs {}",
                                pair_count + 1,
                                prev_r1.id(),
                                record.id()
                            ));
                        }

                        if let (
                            SequenceRecord::Fastq {
                                id: r1_id,
                                desc: r1_desc,
                                seq: r1_seq,
                                qual: r1_qual,
                            },
                            SequenceRecord::Fastq {
                                seq: r2_seq,
                                qual: r2_qual,
                                ..
                            },
                        ) = (&prev_r1, &record)
                        {
                            // Concatenate sequences: R1 + N + R2
                            let mut new_seq = Vec::with_capacity(r1_seq.len() + r2_seq.len() + 1);
                            new_seq.extend_from_slice(&r1_seq);
                            new_seq.push(b'N');
                            new_seq.extend_from_slice(&r2_seq);

                            // Concatenate quality scores: R1 + ! + R2
                            let mut new_qual = Vec::with_capacity(r1_qual.len() + r2_qual.len() + 1);
                            new_qual.extend_from_slice(&r1_qual);
                            new_qual.push(b'!'); // Low quality for 'N'
                            new_qual.extend_from_slice(&r2_qual);


                            let new_record = SequenceRecord::Fastq {
                                id: r1_id.clone(),
                                desc: r1_desc.clone(),
                                seq: Arc::new(new_seq),
                                qual: Arc::new(new_qual),
                            };

                            if tx.send(ParseOutput::Fastq(new_record)).await.is_err() {
                                eprintln!("Failed to send concatenated FASTQ record at pair {}", pair_count + 1);
                                return Err(anyhow!(
                                    "Failed to send concatenated FASTQ record at pair {}",
                                    pair_count + 1
                                ));
                            }
                            pair_count += 1;
                        } else {
                            return Err(anyhow!("Non-FASTQ record in pair at count {}", pair_count + 1));
                        }
                    } else {
                        // This is R1; store and wait for R2
                        r1 = Some(record);
                    }
                }
                _ => {
                    eprintln!("Unexpected ParseOutput at pair count {}", pair_count + 1);
                    continue;
                }
            }
        }

        // Check for unpaired R1
        if r1.is_some() {
            eprintln!("Unpaired R1 record at end of stream");
            return Err(anyhow!("Unpaired R1 record at end of stream"));
        }

        eprintln!("concatenate_paired_reads: Processed {} read pairs", pair_count);
        Ok(())
    });

    Ok((ReceiverStream::new(rx), task))
}


/// De-interleaves an R1/R2 byte stream
///
///
/// # Arguments
/// * `input_rx` - ReceiverStream of ParseOutput containing interleaved FASTQ records (R1, R2, ...).
/// * `ram_temp_dir` - /dev/shm on linux, /var usually elsewhere
/// * `sample_base` - For naming
/// * `buffer_size`
/// * `stall_threshold_secs` - Seconds before logging a stall warning.
///
/// # Returns
/// Result(r1 fifo, r2 fifo, task handle)
/// -
pub async fn deinterleave_fastq_to_fifos(
    mut input_rx: Receiver<ParseOutput>,
    ram_temp_dir: PathBuf,
    sample_base: &str,
    buffer_size: usize,
    stall_threshold_secs: u64,
) -> Result<(PathBuf, PathBuf, JoinHandle<Result<()>>)> {
    // Create FIFO paths in ram_temp_dir
    let r1_fifo_path = ram_temp_dir.join(format!("{}_kallisto_r1.fifo", sample_base));
    let r2_fifo_path = ram_temp_dir.join(format!("{}_kallisto_r2.fifo", sample_base));

    // Run mkfifo for each (async via TokioCommand)
    for path in [&r1_fifo_path, &r2_fifo_path] {
        TokioCommand::new("mkfifo")
            .arg(path)
            .status()
            .await
            .map_err(|e| anyhow!("Failed to create named pipe {:?}: {}", path, e))?;
    }

    // Clone paths for use in closure
    let r1_fifo_path_clone = r1_fifo_path.clone();
    let r2_fifo_path_clone = r2_fifo_path.clone();

    // Spawn writer task: Open FIFOs async, parse & write streaming
    let writer_task = tokio::spawn(async move {
        // Open FIFOs for writing (blocks until Kallisto or reader opens for read)
        let mut r1_file = BufWriter::new(TokioFile::create(&r1_fifo_path_clone).await?);
        let mut r2_file = BufWriter::new(TokioFile::create(&r2_fifo_path_clone).await?);

        let mut buffer: Vec<u8> = Vec::with_capacity(buffer_size * 2);
        let mut lines: Vec<Vec<u8>> = Vec::with_capacity(8); // For one pair (8 lines max)
        let mut prev_r1_id: Option<String> = None;
        let mut is_r1 = true;
        let mut pair_count = 0;
        let mut last_progress = Instant::now();

        while let Some(item) = input_rx.recv().await {
            if last_progress.elapsed() > Duration::from_secs(stall_threshold_secs) {
                eprintln!("deinterleave_fastq_to_fifos: Stall detected at {} read pairs", pair_count);
                last_progress = Instant::now();
            }

            if let ParseOutput::Bytes(chunk) = item {
                buffer.extend_from_slice(&*chunk);

                // Process complete lines from buffer
                loop {
                    if let Some(pos) = memchr(b'\n', &buffer) {
                        let line: Vec<u8> = buffer.drain(..=pos).collect();
                        lines.push(line);

                        // Got 4 lines? Parse & write record
                        if lines.len() == 4 {
                            // Extract ID (trim @ and \n)
                            let id_bytes = &lines[0][1..lines[0].len() - 1];
                            let id = String::from_utf8_lossy(id_bytes).to_string();

                            // Seq & qual (trim \n)
                            let seq = &lines[1][..lines[1].len() - 1];
                            let qual = &lines[3][..lines[3].len() - 1];

                            if seq.len() != qual.len() {
                                return Err(anyhow!("Mismatched seq/qual lengths at record {}", pair_count * 2 + if is_r1 { 1 } else { 2 }));
                            }

                            // Write FASTQ bytes to record_bytes (synchronous, Vec)
                            let mut record_bytes = Vec::with_capacity(id_bytes.len() + seq.len() + qual.len() + 6); // @ + \n + + + \n + 2\n
                            std::io::Write::write_all(&mut record_bytes, b"@")?;
                            std::io::Write::write_all(&mut record_bytes, id_bytes)?;
                            std::io::Write::write_all(&mut record_bytes, b"\n")?;
                            std::io::Write::write_all(&mut record_bytes, seq)?;
                            std::io::Write::write_all(&mut record_bytes, b"\n+\n")?;
                            std::io::Write::write_all(&mut record_bytes, qual)?;
                            std::io::Write::write_all(&mut record_bytes, b"\n")?;

                            if is_r1 {
                                r1_file.write_all(&record_bytes).await?;
                                prev_r1_id = Some(id);
                            } else {
                                // Check ID match (no drops)
                                if let Some(prev_id) = prev_r1_id.take() {
                                    if !compare_read_ids(&prev_id, &id) {
                                        return Err(anyhow!("ID mismatch at pair {}: {} vs {}", pair_count + 1, prev_id, id));
                                    }
                                } else {
                                    return Err(anyhow!("Missing prev R1 ID at pair {}", pair_count + 1));
                                }
                                r2_file.write_all(&record_bytes).await?;
                                pair_count += 1;
                            }

                            is_r1 = !is_r1;
                            lines.clear();
                        }
                    } else {
                        // Incomplete line; wait for next chunk
                        break;
                    }
                }
            } else {
                return Err(anyhow!("Unexpected non-Bytes in FASTQ stream at pair {}", pair_count + 1));
            }
        }

        // Check leftovers
        if !lines.is_empty() {
            return Err(anyhow!("Incomplete FASTQ record at end ({} lines left)", lines.len()));
        }
        if !is_r1 {
            return Err(anyhow!("Unpaired R1 at end"));
        }

        // Flush & close (signals EOF to Kallisto)
        r1_file.flush().await?;
        r2_file.flush().await?;

        eprintln!("deinterleave_fastq_to_fifos: Processed {} read pairs", pair_count);
        Ok(())
    });

    Ok((r1_fifo_path, r2_fifo_path, writer_task))
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
    #[tokio::test]
    async fn test_stream_record_counter_fastq() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let fastq_data = vec![
            ParseOutput::Bytes(b"@read1\nATCG\n+\nIIII\n".to_vec().into()),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
        ];

        tokio::spawn(async move {
            for item in fastq_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_record_counter(rx, false).await?;
        assert_eq!(count, 2, "Should count 2 FASTQ records");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_record_counter_fasta_early_exit() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let fasta_data = vec![
            ParseOutput::Bytes(b">seq1\nATCG\n".to_vec().into()),
            ParseOutput::Fasta(SequenceRecord::Fasta {
                id: "seq2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
            }),
        ];

        tokio::spawn(async move {
            for item in fasta_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_record_counter(rx, true).await?;
        assert_eq!(count, 1, "Should stop at first FASTA record with early_exit");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_record_counter_empty() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        drop(tx); // Close immediately
        let count = stream_record_counter(rx, false).await?;
        assert_eq!(count, 0, "Should count 0 for empty stream");
        Ok(())
    }

    #[tokio::test]
    async fn test_concatenate_paired_reads() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/1".to_string(),
                desc: None,
                seq: Arc::new(b"ATCG".to_vec()),
                qual: Arc::new(b"IIII".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
        ];
        tokio::spawn(async move {
            for record in records {
                tx.send(record).await.unwrap();
            }
        });
        let (concat_rx, task) = concatenate_paired_reads(ReceiverStream::new(rx), 10, 10).await?;
        let mut concat_records = Vec::new();
        let mut stream = concat_rx;  // Directly use concat_rx (already a ReceiverStream)
        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
            concat_records.push(record);
        }
        task.await??;
        assert_eq!(concat_records.len(), 1);
        assert_eq!(concat_records[0].id(), "read1/1");
        assert_eq!(concat_records[0].seq(), b"ATCGNGCTA");
        assert_eq!(concat_records[0].qual(), b"IIII!HHHH");
        Ok(())
    }

}