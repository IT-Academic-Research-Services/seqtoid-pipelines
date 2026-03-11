use std::fs::File;
use std::io::Cursor;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::io::{BufWriter as StdBufWriter, Write};

use std::path::PathBuf;
use std::sync::Arc;
use std::collections::HashSet;
use anyhow::{Result, anyhow, Context};
use log::{self, debug, info, error, warn};

use std::collections::HashMap;
use ahash::AHashMap;
use lazy_static::lazy_static;
use crate::utils::sequence::{DNA, normal_phred_qual_string}; 
use futures::Stream;
use tokio_stream::{self as stream};

use needletail::{parse_fastx_file, FastxReader, parser::{FastqReader}};
use futures::stream::StreamExt;
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tokio::time::{Duration};
use tokio::fs::File as TokioFile;
use futures::future::try_join_all;
use tokio::io::AsyncWriteExt;
use tokio::task::JoinHandle;
use tokio::sync::mpsc::Receiver;
use memchr::memmem;
use fst::{MapBuilder};

use crate::utils::streams::{ParseOutput, ToBytes};
use crate::utils::file::{extension_remover};
use crate::cli::Technology;
use crate::config::defs::{FASTA_TAG, FASTQ_TAG, FASTA_EXTS, FASTQ_EXTS, ReadStats, Taxid, Lineage, CONFORMING_PREAMBLE, PipelineError, RunConfig};
use crate::utils::taxonomy::{get_valid_lineage, generate_locator_work, combine_taxon_loc_json};
use crate::utils::system::compute_phase_concurrency;

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

    pub fn seq_arc(&self) -> Arc<Vec<u8>> {
        match self {
            SequenceRecord::Fasta { seq, .. } => seq.clone(),
            SequenceRecord::Fastq { seq, .. } => seq.clone(),
        }
    }

    pub fn seq_owned(&self) -> Vec<u8> {
        match self {
            SequenceRecord::Fasta { seq, .. } => (**seq).clone(),
            SequenceRecord::Fastq { seq, .. } => (**seq).clone(),
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
        debug!("read_fastq: Starting with path1={:?}, path2={:?}, max_reads={}, min_read_len={:?}, max_read_len={:?}, chunk_size={}",
                  path1, path2, max_reads, min_read_len, max_read_len, chunk_size);

        let mut undersized_reads: u64 = 0;
        let mut oversized_reads: u64 = 0;
        let mut validated_reads: u64 = 0;
        let mut unpaired_r1: u64 = 0;
        let mut unpaired_r2: u64 = 0;

        match path2 {
            None => {
                // Single-end: unchanged
                debug!("read_fastq: Opening single-end FASTQ: {:?}", path1);
                let mut reader = parse_fastx_file(&path1)
                    .map_err(|e| anyhow!("Failed to open FASTQ {}: {}", path1.display(), e))?;
                let mut buffer = Vec::with_capacity(chunk_size);

                while let Some(result) = reader.next() {
                    if validated_reads + undersized_reads + oversized_reads >= max_reads {
                        debug!("read_fastq: Reached max_reads={}", max_reads);
                        break;
                    }
                    let record = result.map_err(|e| {
                        anyhow!("Parse error at read {}: {}", validated_reads + undersized_reads + oversized_reads + 1, e)
                    })?;

                    if record.qual().is_none() {
                        return Err(anyhow!("Expected FASTQ, got FASTA at read {}", validated_reads + undersized_reads + oversized_reads + 1));
                    }

                    let seq_len = record.seq().len();

                    if let Some(min) = min_read_len {
                        if seq_len < min {
                            undersized_reads += 1;
                            continue;
                        }
                    }
                    if let Some(max) = max_read_len {
                        if seq_len > max {
                            oversized_reads += 1;
                            continue;
                        }
                    }

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

                    if buffer.len() >= chunk_size {
                        if tx.send(ParseOutput::Bytes(Arc::new(std::mem::take(&mut buffer)))).await.is_err() {
                            return Err(anyhow!("Failed to send byte chunk at read {}", validated_reads));
                        }
                    }
                }
                if !buffer.is_empty() {
                    if tx.send(ParseOutput::Bytes(Arc::new(buffer))).await.is_err() {
                        return Err(anyhow!("Failed to send final byte chunk"));
                    }
                }
                debug!("read_fastq: Processed {} single-end reads (undersized: {}, validated: {}, oversized: {})",
                          validated_reads + undersized_reads + oversized_reads, undersized_reads, validated_reads, oversized_reads);
                if validated_reads == 0 {
                    warn!("read_fastq: Warning: No reads processed from {:?}", path1);
                }
                Ok(ReadStats {
                    undersized: undersized_reads,
                    validated: validated_reads,
                    oversized: oversized_reads,
                    unpaired_r1: 0,
                    unpaired_r2: 0,
                })
            }
            Some(path2) => {
                if let Some(Technology::ONT) = technology {
                    return Err(anyhow!("Paired-end not supported for ONT"));
                }
                debug!("read_fastq: Opening paired-end FASTQ: R1={:?}, R2={:?}", path1, path2);
                let mut reader1 = parse_fastx_file(&path1)
                    .map_err(|e| anyhow!("Failed to open R1 FASTQ {}: {}", path1.display(), e))?;
                let mut reader2 = parse_fastx_file(&path2)
                    .map_err(|e| anyhow!("Failed to open R2 FASTQ {}: {}", path2.display(), e))?;
                let mut buffer = Vec::with_capacity(chunk_size);

                loop {
                    if validated_reads + undersized_reads + oversized_reads + unpaired_r1 + unpaired_r2 >= max_reads {
                        debug!("read_fastq: Reached max_reads={}", max_reads);
                        break;
                    }

                    // Try to read R1
                    let r1_opt = reader1.next();
                    let r1_record = match r1_opt {
                        Some(Ok(rec)) => Some(rec),
                        Some(Err(e)) => {
                            warn!("R1 parse error at read {}: {} — discarding unpaired",
                                  validated_reads + undersized_reads + oversized_reads + unpaired_r1 + unpaired_r2 + 1, e);
                            unpaired_r1 += 1;
                            continue;
                        }
                        None => None,
                    };

                    // Try to read R2
                    let r2_opt = reader2.next();
                    let r2_record = match r2_opt {
                        Some(Ok(rec)) => Some(rec),
                        Some(Err(e)) => {
                            warn!("R2 parse error at read {}: {} — discarding unpaired",
                                  validated_reads + undersized_reads + oversized_reads + unpaired_r1 + unpaired_r2 + 1, e);
                            unpaired_r2 += 1;
                            continue;
                        }
                        None => None,
                    };

                    match (r1_record, r2_record) {
                        (Some(r1), Some(r2)) => {
                            // Both valid → process pair
                            if r1.qual().is_none() || r2.qual().is_none() {
                                return Err(anyhow!("Expected FASTQ for pair at read {}",
                                                    validated_reads + undersized_reads + oversized_reads + 1));
                            }

                            if !compare_read_ids_bytes(r1.id(), r2.id()) {
                                warn!("ID mismatch at pair {} — discarding unpaired",
                                      validated_reads + undersized_reads + oversized_reads + 1);
                                unpaired_r1 += 1;
                                unpaired_r2 += 1;
                                continue;
                            }

                            let r1_len = r1.seq().len();
                            let r2_len = r2.seq().len();

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

                            // Build R1 bytes
                            let mut r1_bytes = Vec::new();
                            r1_bytes.push(b'@');
                            r1_bytes.extend_from_slice(r1.id());
                            r1_bytes.push(b'\n');
                            r1_bytes.extend_from_slice(&r1.seq());
                            r1_bytes.push(b'\n');
                            r1_bytes.push(b'+');
                            r1_bytes.push(b'\n');
                            r1_bytes.extend_from_slice(r1.qual().unwrap());
                            r1_bytes.push(b'\n');

                            // Build R2 bytes
                            let mut r2_bytes = Vec::new();
                            r2_bytes.push(b'@');
                            r2_bytes.extend_from_slice(r2.id());
                            r2_bytes.push(b'\n');
                            r2_bytes.extend_from_slice(&r2.seq());
                            r2_bytes.push(b'\n');
                            r2_bytes.push(b'+');
                            r2_bytes.push(b'\n');
                            r2_bytes.extend_from_slice(r2.qual().unwrap());
                            r2_bytes.push(b'\n');

                            // Interleave: R1 then R2
                            buffer.extend_from_slice(&r1_bytes);
                            buffer.extend_from_slice(&r2_bytes);
                            validated_reads += 1;

                            if buffer.len() >= chunk_size {
                                if tx.send(ParseOutput::Bytes(Arc::new(std::mem::take(&mut buffer)))).await.is_err() {
                                    return Err(anyhow!("Failed to send byte chunk at read {}", validated_reads));
                                }
                            }
                        }
                        (Some(_), None) => {
                            unpaired_r1 += 1;
                            if unpaired_r1 <= 10 {
                                warn!("Unpaired R1 read after R2 ended — discarding");
                            }
                        }
                        (None, Some(_)) => {
                            unpaired_r2 += 1;
                            if unpaired_r2 <= 10 {
                                warn!("Unpaired R2 read after R1 ended — discarding");
                            }
                        }
                        (None, None) => break,
                    }
                }

                // Send remaining buffer
                if !buffer.is_empty() {
                    if tx.send(ParseOutput::Bytes(Arc::new(buffer))).await.is_err() {
                        return Err(anyhow!("Failed to send final byte chunk"));
                    }
                }

                if unpaired_r1 > 0 || unpaired_r2 > 0 {
                    warn!("Discarded {} unpaired R1 reads and {} unpaired R2 reads from raw input", unpaired_r1, unpaired_r2);
                }

                debug!("read_fastq: Processed {} paired-end reads (undersized: {}, validated: {}, oversized: {}, unpaired R1: {}, unpaired R2: {})",
                          validated_reads + undersized_reads + oversized_reads + unpaired_r1 + unpaired_r2,
                          undersized_reads, validated_reads, oversized_reads, unpaired_r1, unpaired_r2);

                if validated_reads == 0 {
                    warn!("read_fastq: Warning: No valid paired reads processed from R1={:?}, R2={:?}", path1, path2);
                }

                Ok(ReadStats {
                    undersized: undersized_reads,
                    validated: validated_reads,
                    oversized: oversized_reads,
                    unpaired_r1,
                    unpaired_r2,
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
        debug!("read_fasta: Starting with path={:?}, max_records={}, min_seq_len={:?}, max_seq_len={:?}, chunk_size={}",
                  path, max_records, min_seq_len, max_seq_len, chunk_size);

        let mut reader = parse_fastx_file(&path)
            .map_err(|e| {
                anyhow!("Failed to open FASTA {}: {}", path.display(), e)
            })?;

        let mut buffer = Vec::with_capacity(chunk_size);
        let mut record_count: u64 = 0;

        while let Some(result) = reader.next() {
            if record_count >= max_records {
                debug!("read_fasta: Reached max_records={}", max_records);
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
                    debug!("read_fasta: Skipping record {}: seq_len={} < min={}", record_count + 1, seq_len, min);
                    continue;
                }
            }
            if let Some(max) = max_seq_len {
                if seq_len > max {
                    debug!("read_fasta: Skipping record {}: seq_len={} > max={}", record_count + 1, seq_len, max);
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

        debug!("read_fasta: Processed {} FASTA records", record_count);
        if record_count == 0 {
            warn!("read_fasta: Warning: No records processed from {:?}", path);
        }
        Ok(())
    });

    Ok(rx)
}


/// Write any ParseOutput stream that contains FASTA records to a file.
/// The caller decides the destination path (RAM temp, scratch, final output, …).
///
/// # Arguments
///
/// * `stream` - ParseOutput stream, must eb fASTA foirmat
/// * 'dest_path' - location to write to
/// * 'mbuffer_capacity'
/// # Returns      Vec<JoinHandle<Result<()>>>,
pub fn write_fasta_stream_to_file(
    stream: ReceiverStream<ParseOutput>,
    dest_path: PathBuf,
    buffer_capacity: usize,
) -> JoinHandle<Result<()>> {
    let dest_path_clone = dest_path.clone();

    tokio::spawn(async move {
        let file = TokioFile::create(&dest_path_clone)
            .await
            .map_err(|e| anyhow!("Cannot create FASTA file {}: {}", dest_path_clone.display(), e))?;

        let mut writer = tokio::io::BufWriter::with_capacity(buffer_capacity, file);

        let mut stream = stream;
        while let Some(item) = stream.next().await {
            let bytes = item
                .to_bytes()
                .map_err(|e| anyhow!("Failed to convert FASTA record to bytes: {}", e))?;

            writer
                .write_all(&bytes)
                .await
                .map_err(|e| anyhow!("Write error to {}: {}", dest_path_clone.display(), e))?;
        }

        writer
            .flush()
            .await
            .map_err(|e| anyhow!("Flush error to {}: {}", dest_path_clone.display(), e))?;

        info!("FASTA written to {}", dest_path_clone.display());
        Ok(())
    })
}

fn compare_read_ids_bytes(id1: &[u8], id2: &[u8]) -> bool {
    // Fast path: exact byte match (very common when IDs are identical)
    if id1 == id2 {
        return true;
    }

    // Convert to str with fallback (safe, no panic)
    let s1 = std::str::from_utf8(id1).unwrap_or("");
    let s2 = std::str::from_utf8(id2).unwrap_or("");

    // Expanded suffix list — ordered roughly by frequency / importance
    // Each group is labeled with the platforms / formats it covers
    let suffixes: [&str; 24] = [
        // Classic Illumina / most common overall
        "/1", "/2",

        // Older Illumina, some custom pipelines
        ".1", ".2",

        // Illumina with space + colon (HiSeq, NovaSeq, NextSeq)
        " 1:", " 2:",

        // Full Illumina flowcell/lane info suffix (very common on NovaSeq)
        " 1:N:0:", " 2:N:0:",

        // Variant with slash instead of space (some demux tools)
        "/1:N:0:", "/2:N:0:",

        // BGISEQ, DNBSEQ, MGISEQ, many Chinese platforms
        "_1", "_2",

        // Very old Illumina CASAVA format (pre-2011, rare but still seen)
        " 1#0/1", " 2#0/2",

        // Nanopore, Oxford Nanopore, some custom / long-read paired
        "|1", "|2",

        // Rare 10x Genomics / Parse / custom variants
        " 1#", " 2#",

        // 10x Genomics, Parse Biosciences, some single-cell demux
        "/R1", "/R2",

        // Sometimes appears in headers instead of file names
        "_R1", "_R2",

        // Rare variant without N (some older or custom)
        " 1:0", " 2:0",
    ];

    let mut base1 = s1;
    let mut base2 = s2;

    // Strip the longest matching suffix from each ID
    // We break after the first match to avoid over-trimming
    for &suffix in &suffixes {
        if base1.ends_with(suffix) {
            base1 = &base1[..base1.len() - suffix.len()];
            break;
        }
    }

    for &suffix in &suffixes {
        if base2.ends_with(suffix) {
            base2 = &base2[..base2.len() - suffix.len()];
            break;
        }
    }

    // Clean up any trailing whitespace, colon, or hash that might remain
    base1 = base1.trim_end_matches(|c: char| c.is_whitespace() || c == ':' || c == '#');
    base2 = base2.trim_end_matches(|c: char| c.is_whitespace() || c == ':' || c == '#');

    // Final comparison
    base1 == base2
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
pub fn compare_read_ids(id1: &str, id2: &str) -> bool {
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
                            return Err(anyhow!(
                                "Failed to send filtered FASTQ record at count {}",
                                count + 1
                            ));
                        }
                        count += 1;
                    }
                }
                _ => {
                    warn!("Unexpected ParseOutput at count {}", count + 1);
                    continue; // Skip non-FASTQ items
                }
            }
        }
        debug!("Filtered {} FASTQ records from Kraken2 classified stream", count);
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
                warn!("parse_byte_stream_to_fastq: Stall detected at {} records", record_count);
                last_progress = tokio::time::Instant::now();
            }

            match item {
                ParseOutput::Bytes(bytes) => {
                    full_buffer.extend_from_slice(&*bytes);
                }
                _ => {
                    warn!("parse_byte_stream_to_fastq: Unexpected non-Bytes ParseOutput at record {}", record_count + 1);
                    continue;
                }
            }
        }

        if full_buffer.is_empty() {
            return Err(anyhow!("Empty byte buffer for FASTQ parsing"));
        }

        // Parse accumulated bytes into FASTQ records
        let cursor = Cursor::new(full_buffer);
        let mut reader = FastqReader::new(cursor);
        while let Some(result) = reader.next() {
            match result {
                Ok(record) => {
                    let owned_record: SequenceRecord = record.to_owned().into();
                    if tx.send(ParseOutput::Fastq(owned_record)).await.is_err() {
                        return Err(anyhow!("Failed to send FASTQ record at count {}", record_count + 1));
                    }
                    record_count += 1;
                }
                Err(e) => {
                    return Err(anyhow!("FASTQ parsing error at count {}: {}", record_count + 1, e));
                }
            }
        }

        debug!("parse_byte_stream_to_fastq: Parsed {} FASTQ records", record_count);
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
                warn!("concatenate_paired_reads: Stall detected at {} read pairs", pair_count);
                last_progress = tokio::time::Instant::now();
            }

            match item {
                ParseOutput::Fastq(record) => {
                    if let Some(prev_r1) = r1.take() {
                        // This is R2; process the pair
                        if !compare_read_ids(prev_r1.id(), record.id()) {
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
                    warn!("Unexpected ParseOutput at pair count {}", pair_count + 1);
                    continue;
                }
            }
        }

        // Check for unpaired R1
        if r1.is_some() {
            return Err(anyhow!("Unpaired R1 record at end of stream"));
        }

        debug!("concatenate_paired_reads: Processed {} read pairs", pair_count);
        Ok(())
    });

    Ok((ReceiverStream::new(rx), task))
}

/// Reads a FASTA and generates an index file for that FASTA for only the offsets of the
/// accessions in the file. This allows faster lookups of the accessions. Especially
/// useful for huge files like NT or NR.
///
/// # Arguments
/// * `fasta_path` - Path to a fasta file
/// * `index_ath` - Name of the resulting index file
///
/// # Returns
/// Result \ , Box error
pub fn build_fasta_index(
    fasta_path: &PathBuf,
    index_path: &PathBuf,
    verify: bool,
) -> Result<()> {
    let file = File::open(fasta_path)
        .with_context(|| format!("Failed to open FASTA file: {}", fasta_path.display()))?;

    let mut reader = BufReader::new(file);
    let mut line = String::with_capacity(1024);
    let mut entries = Vec::with_capacity(8_000_000);  // trying to size for NR

    let mut pos: u64 = 0;

    while reader.read_line(&mut line)? > 0 {
        if line.starts_with('>') {
            let header = line[1..].trim_start();
            let acc = header
                .split_whitespace()
                .next()
                .ok_or_else(|| anyhow!("Empty FASTA header"))?
                .split('.')
                .next()
                .unwrap()
                .to_owned();

            entries.push((acc, pos));

            // Optional verification
            if verify {
                let mut byte = [0u8; 1];
                reader.seek(SeekFrom::Start(pos))?;
                reader.read_exact(&mut byte)?;
                if byte[0] != b'>' {
                    return Err(anyhow!(
                        "VERIFICATION FAILED at offset {pos}: expected '>', got 0x{:02x}",
                        byte[0]
                    ));
                }
                // Return to where we were (after the line we just read)
                reader.seek(SeekFrom::Start(pos + line.len() as u64))?;
            }
        }
        pos += line.len() as u64;
        line.clear();
    }

    // Dedup + sort (keeps first occurrence)
    entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    entries.dedup_by(|a, b| a.0 == b.0);

    info!("Writing FST index with {} unique accessions → {}", entries.len(), index_path.display());

    let mut builder = MapBuilder::new(StdBufWriter::new(
        File::create(index_path)
            .with_context(|| format!("Cannot create index file: {}", index_path.display()))?,
    ))?;

    for (acc, offset) in entries {
        builder.insert(acc, offset)?;
    }

    builder.finish()?;
    info!("FST index built successfully");

    Ok(())
}


pub async fn generate_taxid_fasta(
    config: Arc<RunConfig>,
    mapped_contigs_stream: ReceiverStream<ParseOutput>,
    unidentified_contigs_stream: ReceiverStream<ParseOutput>,
    nt_hit_summary_stream: ReceiverStream<ParseOutput>,
    nr_hit_summary_stream: ReceiverStream<ParseOutput>,
    lineage_map: Arc<AHashMap<Taxid, Lineage>>,
) -> Result<(
    mpsc::Receiver<ParseOutput>,           // mapped contigs with taxid headers
    mpsc::Receiver<ParseOutput>,           // combined (mapped + conformed unmapped)
    JoinHandle<Result<()>>,                // load_nt_task
    JoinHandle<Result<()>>,                // load_nr_task
    JoinHandle<Result<()>>,                // main_task
)> {
    // Output channels
    let (mapped_tx, mapped_rx) = mpsc::channel(1024);
    let (combined_tx, combined_rx) = mpsc::channel(2048);

    // Channels for completed hit maps
    let (nt_hits_tx, mut nt_hits_rx) = mpsc::channel::<Arc<AHashMap<String, (Taxid, u8)>>>(1);
    let (nr_hits_tx, mut nr_hits_rx) = mpsc::channel::<Arc<AHashMap<String, (Taxid, u8)>>>(1);

    // Load NT hit summary → build map
    let load_nt_task = tokio::spawn({
        let nt_stream = nt_hit_summary_stream;
        async move {
            let mut hits = AHashMap::with_capacity(100_000);
            let mut stream = nt_stream;
            while let Some(item) = stream.next().await {
                if let ParseOutput::Bytes(bytes) = item {
                    let line = String::from_utf8_lossy(&bytes).trim_end().to_string();
                    if line.is_empty() { continue; }
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 4 { continue; }
                    let read_id = fields[0].to_string();
                    let level: u8 = fields[1].parse().unwrap_or(255);
                    let taxid: Taxid = fields[3].parse().unwrap_or(-1);
                    hits.insert(read_id, (taxid, level));
                }
            }
            nt_hits_tx.send(Arc::new(hits)).await.map_err(|_| anyhow!("NT hits tx dropped"))?;
            Ok(())
        }
    });

    // Load NR hit summary → build map
    let load_nr_task = tokio::spawn({
        let nr_stream = nr_hit_summary_stream;
        async move {
            let mut hits = AHashMap::with_capacity(100_000);
            let mut stream = nr_stream;
            while let Some(item) = stream.next().await {
                if let ParseOutput::Bytes(bytes) = item {
                    let line = String::from_utf8_lossy(&bytes).trim_end().to_string();
                    if line.is_empty() { continue; }
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 4 { continue; }
                    let read_id = fields[0].to_string();
                    let level: u8 = fields[1].parse().unwrap_or(255);
                    let taxid: Taxid = fields[3].parse().unwrap_or(-1);
                    hits.insert(read_id, (taxid, level));
                }
            }
            nr_hits_tx.send(Arc::new(hits)).await.map_err(|_| anyhow!("NR hits tx dropped"))?;
            Ok(())
        }
    });

    // Main processing task
    let main_task = tokio::spawn(async move {
        // Wait for both maps
        let nt_hits = nt_hits_rx.recv().await.ok_or(anyhow!("NT hits map not received"))?;
        let nr_hits = nr_hits_rx.recv().await.ok_or(anyhow!("NR hits map not received"))?;

        let concurrency = compute_phase_concurrency(
            &config,
            "generate_taxid_fasta",
            0.5,           // ~0.5 GB per thread max
            4.0,           // 4 threads per core
            64,            // cap
            8,             // min for meaningful parallelism
        );

        const BATCH_SIZE: usize = 1000;

        // Process annotated (mapped) contigs in parallel
        {
            let (batch_tx, batch_rx) = mpsc::channel::<Vec<SequenceRecord>>(concurrency * 2);
            let batch_rx = Arc::new(tokio::sync::Mutex::new(batch_rx));
            
            let mut worker_handles = Vec::with_capacity(concurrency);
            for _ in 0..concurrency {
                let rx = batch_rx.clone();
                let m_tx = mapped_tx.clone();
                let c_tx = combined_tx.clone();
                let nt_hits = nt_hits.clone();
                let nr_hits = nr_hits.clone();
                let lineage_map = lineage_map.clone();

                worker_handles.push(tokio::spawn(async move {
                    loop {
                        let batch = {
                            let mut lock = rx.lock().await;
                            lock.recv().await
                        };

                        let batch = match batch {
                            Some(b) => b,
                            None => break,
                        };

                        for rec in batch {
                            let annotated_id = rec.id();
                            let parts: Vec<&str> = annotated_id.split(':').collect();
                            let contig_id = parts.last().unwrap_or(&"").to_string();

                            let nr_lineage = get_valid_lineage(&nr_hits, &lineage_map, &contig_id);
                            let nt_lineage = get_valid_lineage(&nt_hits, &lineage_map, &contig_id);

                            let new_header = format!(
                                ">family_nr:{}:family_nt:{}:genus_nr:{}:genus_nt:{}:species_nr:{}:species_nt:{}:{}\n",
                                nr_lineage[2], nt_lineage[2],
                                nr_lineage[1], nt_lineage[1],
                                nr_lineage[0], nt_lineage[0],
                                annotated_id
                            );

                            let (new_id, new_desc) = parse_header(new_header.trim_start_matches('>').as_bytes(), '>');
                            let seq_arc = rec.seq_arc();
                            let new_rec = SequenceRecord::Fasta {
                                id: new_id,
                                desc: new_desc,
                                seq: seq_arc,
                            };

                            m_tx.send(ParseOutput::Fasta(new_rec.clone())).await.map_err(|_| anyhow!("mapped_tx dropped"))?;
                            c_tx.send(ParseOutput::Fasta(new_rec)).await.map_err(|_| anyhow!("combined_tx dropped"))?;
                        }
                    }
                    Ok::<(), anyhow::Error>(())
                }));
            }

            let mut mapped_stream = mapped_contigs_stream;
            let mut current_batch = Vec::with_capacity(BATCH_SIZE);
            while let Some(item) = mapped_stream.next().await {
                match item {
                    ParseOutput::Fasta(rec) => {
                        current_batch.push(rec);
                        if current_batch.len() >= BATCH_SIZE {
                            batch_tx.send(std::mem::take(&mut current_batch))
                                .await
                                .map_err(|_| anyhow!("batch_tx dropped unexpectedly during mapped_stream processing"))?;
                        }
                    }
                    _ => return Err(anyhow!("Unexpected item type in mapped_contigs_stream: expected Fasta")),
                }
            }
            if !current_batch.is_empty() {
                batch_tx.send(current_batch).await.map_err(|_| anyhow!("batch_tx dropped"))?;
            }
            drop(batch_tx);

            for handle in worker_handles {
                handle.await??;
            }
        }

        // Process unidentified contigs → conform headers
        // Also parallelize this since it can be many contigs
        {
            let (unid_batch_tx, unid_batch_rx) = mpsc::channel::<Vec<SequenceRecord>>(concurrency * 2);
            let unid_batch_rx = Arc::new(tokio::sync::Mutex::new(unid_batch_rx));
            let mut unid_worker_handles = Vec::with_capacity(concurrency);

            for _ in 0..concurrency {
                let rx = unid_batch_rx.clone();
                let c_tx = combined_tx.clone();
                unid_worker_handles.push(tokio::spawn(async move {
                    loop {
                        let batch = {
                            let mut lock = rx.lock().await;
                            lock.recv().await
                        };
                        let batch = match batch {
                            Some(b) => b,
                            None => break,
                        };
                        for rec in batch {
                            let conformed_header = format!("{}{}", CONFORMING_PREAMBLE, rec.id());
                            let (new_id, new_desc) = parse_header(conformed_header.as_bytes(), '>');
                            let seq_arc = rec.seq_arc();
                            let new_rec = SequenceRecord::Fasta {
                                id: new_id,
                                desc: new_desc,
                                seq: seq_arc,
                            };
                            c_tx.send(ParseOutput::Fasta(new_rec)).await.map_err(|_| anyhow!("combined_tx dropped"))?;
                        }
                    }
                    Ok::<(), anyhow::Error>(())
                }));
            }

            let mut unid_stream = unidentified_contigs_stream;
            let mut current_unid_batch = Vec::with_capacity(BATCH_SIZE);
            while let Some(item) = unid_stream.next().await {
                match item {
                    ParseOutput::Fasta(rec) => {
                        current_unid_batch.push(rec);
                        if current_unid_batch.len() >= BATCH_SIZE {
                            unid_batch_tx.send(std::mem::take(&mut current_unid_batch))
                                .await
                                .map_err(|_| anyhow!("unid_batch_tx dropped unexpectedly during unidentified_contigs_stream processing"))?;
                        }
                    }
                    _ => return Err(anyhow!("Unexpected item type in unidentified_contigs_stream: expected Fasta")),
                }
            }
            if !current_unid_batch.is_empty() {
                unid_batch_tx.send(current_unid_batch).await.map_err(|_| anyhow!("unid_batch_tx dropped"))?;
            }
            drop(unid_batch_tx);

            for handle in unid_worker_handles {
                handle.await??;
            }
        }

        drop(mapped_tx);
        drop(combined_tx);
        Ok(())
    });

    Ok((
        mapped_rx,
        combined_rx,
        load_nt_task,
        load_nr_task,
        main_task,
    ))
}



pub async fn generate_taxid_locator(
    fasta_stream_rx: tokio::sync::mpsc::Receiver<ParseOutput>,
    assembly_dir: PathBuf,
) -> Result<
    (
        Vec<PathBuf>,
        Vec<JoinHandle<Result<()>>>,
        Vec<Receiver<Result<()>>>,
    ),
    PipelineError,
> {

    let mut fasta_stream = ReceiverStream::new(fasta_stream_rx);
    let mut records: Vec<(String, String)> = Vec::with_capacity(1_000_000);

    while let Some(item) = fasta_stream.next().await {
        let bytes = item
            .to_bytes()
            .map_err(|e| PipelineError::Other(anyhow!("ParseOutput to_bytes failed: {}", e)))?;

        if bytes.is_empty() {
            continue;
        }

        let record_str = String::from_utf8(bytes.to_vec())
            .map_err(|e| PipelineError::Other(anyhow!("Invalid UTF-8 in FASTA record: {}", e)))?;

        let mut lines = record_str.lines();
        let header_line = match lines.next() {
            Some(l) => l,
            None => continue,
        };
        let seq_line = match lines.next() {
            Some(l) => l,
            None => {
                warn!("FASTA record missing sequence line");
                continue;
            }
        };

        if !header_line.starts_with('>') {
            warn!("Invalid FASTA header (missing '>'): {}", header_line);
            continue;
        }

        let header = header_line[1..].trim().to_string(); // strip '>'
        let seq = seq_line.trim().to_string();

        records.push((header, seq));
    }

    info!("Collected {} mapped contigs for taxid locator generation", records.len());

    if records.is_empty() {
        // Write empty outputs
        let empty_paths = vec![
            assembly_dir.join("refined_taxid_annot_sorted_nt.fasta"),
            assembly_dir.join("refined_taxid_annot_sorted_nr.fasta"),
        ];
        return Ok((empty_paths, vec![], vec![]));
    }

    let records_arc = Arc::new(records);

    let levels = ["species", "genus", "family"];
    let hit_types = ["NT", "NR"];

    let mut output_paths: Vec<PathBuf> = Vec::new();
    let mut locator_tasks: Vec<JoinHandle<Result<()>>> = Vec::new();

    for &level in &levels {
        for &hit_type in &hit_types {
            let taxid_field = format!("{}_{}", level, hit_type.to_lowercase());

            let fa_name = if level == "species" {
                format!("refined_taxid_annot_sorted_{}.fasta", hit_type.to_lowercase())
            } else {
                format!("refined_taxid_annot_sorted_{}_{}.fasta", level, hit_type.to_lowercase())
            };
            let output_fa = assembly_dir.join(&fa_name);
            let output_json = assembly_dir.join(format!(
                "refined_taxid_locations_{}_{}.json",
                level,
                hit_type.to_lowercase()
            ));

            let records_clone = records_arc.clone();
            let taxid_field_clone = taxid_field.clone();
            let hit_type_clone = hit_type.to_string();
            let output_fa_clone = output_fa.clone();
            let output_json_clone = output_json.clone();

            let handle = tokio::task::spawn_blocking(move || {
                generate_locator_work(
                    records_clone,
                    taxid_field_clone,
                    hit_type_clone,
                    output_fa_clone,
                    output_json_clone,
                )
            });

            locator_tasks.push(handle);

            output_paths.push(output_fa);
            output_paths.push(output_json);
        }
    }

    let json_paths: Vec<PathBuf> = output_paths
        .iter()
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("json"))
        .cloned()
        .collect();

    let combined_path = assembly_dir.join("refined_taxid_locations_combined.json");
    let combined_path_for_worker = combined_path.clone(); // ← explicit clone

    let combine_handle = tokio::task::spawn_blocking(move || {
        combine_taxon_loc_json(json_paths, combined_path_for_worker)
    });

    locator_tasks.push(combine_handle);
    output_paths.push(combined_path);

    Ok((output_paths, locator_tasks, vec![]))
}

pub async fn filter_fastq_to_bytes_stream(
    mut input_stream: ReceiverStream<ParseOutput>,
    headers: Arc<HashSet<String>>,
) -> ReceiverStream<ParseOutput> {
    let (tx, rx) = tokio::sync::mpsc::channel(2048);

    tokio::spawn(async move {
        while let Some(item) = input_stream.next().await {
            if let ParseOutput::Fastq(record) = item {
                if headers.contains(&record.id()[..]) {
                    match record.to_bytes() {
                        Ok(vec) => {
                            let arc_data = Arc::new(vec);
                            let _ = tx.send(ParseOutput::Bytes(arc_data)).await;
                        }
                        Err(e) => error!("Failed to serialize FASTQ record: {}", e),
                    }
                }
                // else: explicitly skipped
            }
        }
    });

    ReceiverStream::new(rx)
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