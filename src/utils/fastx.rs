use crate::utils::taxonomy::get_taxid;
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
use tokio::io::AsyncSeekExt;
use memchr::{memmem, memchr3};
use fst::{MapBuilder};
use once_cell::sync::Lazy;
use bytes::Bytes;

use crate::utils::streams::{ParseOutput, ToBytes};
use crate::utils::file::{extension_remover};
use crate::cli::Technology;
use crate::config::defs::{StreamDataType, FASTA_TAG, FASTQ_TAG, FASTA_EXTS, FASTQ_EXTS, ReadStats, Taxid, Lineage, CONFORMING_PREAMBLE, PipelineError, RunConfig, TaxonSeqLocation, PairingMode, SIMD_LEVEL, SimdLevel};
use crate::utils::taxonomy::{get_valid_lineage, combine_taxon_loc_json};
use crate::utils::system::{compute_phase_concurrency};

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
        seq: Bytes,
    },
    Fastq {
        id: String,
        desc: Option<String>,
        seq: Bytes,
        qual: Bytes,
    },
}

impl SequenceRecord {
    pub fn id(&self) -> &str {
        match self {
            SequenceRecord::Fasta { id, .. } => id,
            SequenceRecord::Fastq { id, .. } => id,
        }
    }

    #[allow(dead_code)]
    pub fn desc(&self) -> Option<&str> {
        match self {
            SequenceRecord::Fasta { desc, .. } => desc.as_deref(),
            SequenceRecord::Fastq { desc, .. } => desc.as_deref(),
        }
    }


    /// Cheap clone of the sequence bytes (refcount bump only).
    pub fn seq(&self) -> Bytes {
        match self {
            SequenceRecord::Fasta { seq, .. } => seq.clone(),
            SequenceRecord::Fastq { seq, .. } => seq.clone(),
        }
    }



    #[allow(dead_code)]
    pub fn qual(&self) -> Bytes {
        match self {
            SequenceRecord::Fasta { .. } => Bytes::new(),
            SequenceRecord::Fastq { qual, .. } => qual.clone(),
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
                seq: Bytes::from(record.seq().to_vec()),
                qual: Bytes::from(qual.to_vec()),
            }
        } else {
            SequenceRecord::Fasta {
                id,
                desc,
                seq: Bytes::from(record.seq().to_vec()),
            }
        }
    }
}



/// Validates that a sequence only contains allowed nucleotide bases.
/// Accepts `&[u8]`, `Bytes`, `Vec<u8>`, `&Bytes`, etc.
pub fn validate_sequence(seq: impl AsRef<[u8]>, valid_bases: &[u8]) -> Result<()> {
    let seq = seq.as_ref();
    let valid: HashSet<u8> = valid_bases.iter().copied().collect();

    if let Some(&invalid) = seq.iter().find(|b| !valid.contains(b)) {
        return Err(anyhow!(
            "Invalid nucleotide '{}' found in sequence",
            invalid as char
        ));
    }
    Ok(())
}

/// Parallel validation for very large sequences (e.g. long contigs).
/// Uses cheap `Bytes` clones to distribute chunks across threads.
pub async fn validate_sequence_parallel(
    seq: Bytes,
    valid_bases: &[u8],
    num_threads: usize,
) -> Result<()> {
    if seq.is_empty() || num_threads == 0 {
        return Ok(());
    }

    let valid: HashSet<u8> = valid_bases.iter().copied().collect();
    let len = seq.len();
    let chunk_size = (len + num_threads - 1) / num_threads;

    let mut handles = Vec::with_capacity(num_threads);

    for i in 0..num_threads {
        let start = i * chunk_size;
        if start >= len {
            break;
        }
        let end = (start + chunk_size).min(len);

        // Cheap clone (refcount bump) + zero-copy slice
        let chunk = seq.slice(start..end);
        let valid = valid.clone();

        let handle = tokio::spawn(async move {
            if let Some(&invalid) = chunk.iter().find(|b| !valid.contains(b)) {
                Err(format!("Invalid nucleotide '{}' found in sequence", invalid as char))
            } else {
                Ok(())
            }
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.await?.map_err(|e| anyhow!(e))?;
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
/// # Arguments
/// * `head`   - Header bytes (needletail has already stripped the leading `>`/`@`).
/// * `prefix` - The leading character (`>` or `@`) to strip from the id field.
/// Note: needletail already stripped the leading > or @.
/// # Returns
/// `(id, desc)` — desc is `None` when the header has no whitespace.
pub fn parse_header(head: &[u8], prefix: char) -> (String, Option<String>) {
    PARSE_HEADER_FN(head, prefix)
}

/// Scalar path — reference implementation.
fn parse_header_scalar(head: &[u8], prefix: char) -> (String, Option<String>) {
    match memchr::memchr3(b' ', b'\t', b'\n', head) {
        Some(pos) => parse_from_pos(head, pos, prefix),
        None => {
            // id-only header
            let id = String::from_utf8_lossy(head)
                .trim_start_matches(prefix)
                .to_string();
            (id, None)
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn parse_header_avx512_inner(head: &[u8], prefix: char) -> (String, Option<String>) {
    use std::arch::x86_64::*;

    let len = head.len();
    if len == 0 {
        return (String::new(), None);
    }

    let space_splat = _mm512_set1_epi8(b' ' as i8);
    let tab_splat   = _mm512_set1_epi8(b'\t' as i8);
    let nl_splat    = _mm512_set1_epi8(b'\n' as i8);

    let mut i = 0usize;
    while i + 64 <= len {
        let chunk = _mm512_loadu_si512(head.as_ptr().add(i).cast::<__m512i>());
        let mask: u64 =
            _mm512_cmpeq_epi8_mask(chunk, space_splat) |
                _mm512_cmpeq_epi8_mask(chunk, tab_splat)   |
                _mm512_cmpeq_epi8_mask(chunk, nl_splat);

        if mask != 0 {
            let pos = i + mask.trailing_zeros() as usize;
            return parse_from_pos(head, pos, prefix);
        }
        i += 64;
    }

    // Fallback: use scalar on the full header for correctness
    // (handles long headers and headers with no whitespace)
    parse_header_scalar(head, prefix)
}

/// Core logic: given a known position of the first whitespace, extract id + desc.
/// This is shared between scalar and AVX-512 paths.
fn parse_from_pos(head: &[u8], pos: usize, prefix: char) -> (String, Option<String>) {
    let id_bytes = &head[..pos];
    let desc_bytes = &head[pos + 1..];

    let id = String::from_utf8_lossy(id_bytes)
        .trim_start_matches(prefix)
        .to_string();

    let desc = if desc_bytes.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(desc_bytes).into_owned())
    };

    (id, desc)
}

#[cfg(target_arch = "x86_64")]
fn parse_header_avx512(head: &[u8], prefix: char) -> (String, Option<String>) {
    unsafe { parse_header_avx512_inner(head, prefix) }
}

static PARSE_HEADER_FN: Lazy<fn(&[u8], char) -> (String, Option<String>)> = Lazy::new(|| {
    #[cfg(target_arch = "x86_64")]
    if matches!(*SIMD_LEVEL, SimdLevel::Avx512) {
        debug!("parse_header: using AVX-512 path");
        return parse_header_avx512;
    }
    debug!("parse_header: using scalar path");
    parse_header_scalar
});




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
    pairing_mode: PairingMode,
    technology: Option<Technology>,
    max_reads: u64,
    min_read_len: Option<usize>,
    max_read_len: Option<usize>,
    tag: &str,
    config: &RunConfig,         // ← added
) -> Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<ReadStats, anyhow::Error>>)> {

    let effective_chunk = crate::utils::system::compute_buffer_size(
        config,
        "read_fastq",
        match technology {
            Some(Technology::ONT) => StreamDataType::OntFastq,
            _ => StreamDataType::IlluminaFastq,
        },
        1.1,                        // light boost for reader
    );

    let (tx, rx) = mpsc::channel(effective_chunk / 1024); // ~1KB per item

    let tag_owned = tag.to_string();
    let task = tokio::spawn(async move {
        let tag = &tag_owned;
        debug!("read_fastq: Starting with path1={:?}, path2={:?}, max_reads={}, min_read_len={:?}, max_read_len={:?}, effective_chunk={}",
                  path1, path2, max_reads, min_read_len, max_read_len, effective_chunk);

        let mut undersized_reads: u64 = 0;
        let mut oversized_reads: u64 = 0;
        let mut validated_reads: u64 = 0;
        let mut unpaired_r1: u64 = 0;
        let mut unpaired_r2: u64 = 0;

        match path2 {
            None => {
                // Single-end
                debug!("read_fastq: Opening single-end FASTQ: {:?}", path1);
                let mut reader = parse_fastx_file(&path1)
                    .map_err(|e| anyhow!("Failed to open FASTQ {}: {}", path1.display(), e))?;
                let mut buffer = Vec::with_capacity(effective_chunk);

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

                    if buffer.len() >= effective_chunk {
                        if tx.send(ParseOutput::Bytes(Bytes::from(std::mem::take(&mut buffer)))).await.is_err() {
                            return Err(anyhow!("Failed to send byte chunk at read {}", validated_reads));
                        }
                    }
                }

                if !buffer.is_empty() {
                    if tx.send(ParseOutput::Bytes(Bytes::from(buffer))).await.is_err() {
                        return Err(anyhow!("Failed to send final byte chunk"));
                    }
                }

                debug!("read_fastq: Processed {} single-end reads (undersized: {}, validated: {}, oversized: {})",
                          validated_reads + undersized_reads + oversized_reads, undersized_reads, validated_reads, oversized_reads);

                if validated_reads == 0 {
                    warn!("{}: No reads processed from {:?}", tag, path1);
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

                let mut buffer = Vec::with_capacity(effective_chunk);

                loop {
                    if validated_reads + undersized_reads + oversized_reads + unpaired_r1 + unpaired_r2 >= max_reads {
                        debug!("read_fastq: Reached max_reads={}", max_reads);
                        break;
                    }

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
                            if r1.qual().is_none() || r2.qual().is_none() {
                                return Err(anyhow!("Expected FASTQ for pair at read {}", validated_reads + undersized_reads + oversized_reads + 1));
                            }

                            if matches!(pairing_mode, PairingMode::Strict) && !compare_read_ids_bytes(r1.id(), r2.id()) {
                                warn!("{}: ID mismatch at pair {} — discarding unpaired: {} vs {}",
                                      tag,
                                      validated_reads + undersized_reads + oversized_reads + 1,
                                      String::from_utf8_lossy(r1.id()),
                                      String::from_utf8_lossy(r2.id()));
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

                            // Build R1 + R2 interleaved
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

                            buffer.extend_from_slice(&r1_bytes);
                            buffer.extend_from_slice(&r2_bytes);
                            validated_reads += 1;

                            if buffer.len() >= effective_chunk {
                                if tx.send(ParseOutput::Bytes(Bytes::from(std::mem::take(&mut buffer)))).await.is_err() {
                                    return Err(anyhow!("Failed to send byte chunk at read {}", validated_reads));
                                }
                            }
                        }
                        (Some(_), None) => { unpaired_r1 += 1; }
                        (None, Some(_)) => { unpaired_r2 += 1; }
                        (None, None) => break,
                    }
                }

                if !buffer.is_empty() {
                    if tx.send(ParseOutput::Bytes(Bytes::from(buffer))).await.is_err() {
                        return Err(anyhow!("Failed to send final byte chunk"));
                    }
                }

                if unpaired_r1 > 0 || unpaired_r2 > 0 {
                    warn!("{}: Discarded {} unpaired R1 and {} unpaired R2 reads", tag, unpaired_r1, unpaired_r2);
                }

                debug!("read_fastq: Processed {} paired-end reads (undersized: {}, validated: {}, oversized: {}, unpaired R1: {}, unpaired R2: {})",
                          validated_reads + undersized_reads + oversized_reads + unpaired_r1 + unpaired_r2,
                          undersized_reads, validated_reads, oversized_reads, unpaired_r1, unpaired_r2);

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
    config: &RunConfig,
) -> Result<mpsc::Receiver<ParseOutput>> {

    let effective_chunk = crate::utils::system::compute_buffer_size(
        config,
        "read_fasta",
        StreamDataType::JustBytes,
        1.2,                         // moderate boost for FASTA reader
    );

    let (tx, rx) = mpsc::channel(effective_chunk / 1024);

    tokio::spawn(async move {
        debug!("read_fasta: Starting with path={:?}, max_records={}, min_seq_len={:?}, max_seq_len={:?}, effective_chunk={}",
                  path, max_records, min_seq_len, max_seq_len, effective_chunk);

        let mut reader = parse_fastx_file(&path)
            .map_err(|e| anyhow!("Failed to open FASTA {}: {}", path.display(), e))?;

        let mut buffer = Vec::with_capacity(effective_chunk);
        let mut record_count: u64 = 0;

        while let Some(result) = reader.next() {
            if record_count >= max_records {
                debug!("read_fasta: Reached max_records={}", max_records);
                break;
            }
            let record = result.map_err(|e| {
                anyhow!("Parse error at record {}: {}", record_count + 1, e)
            })?;

            if record.qual().is_some() {
                return Err(anyhow!("Expected FASTA, got FASTQ at record {}", record_count + 1));
            }

            let seq_len = record.seq().len();

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

            let mut record_bytes = Vec::new();
            record_bytes.push(b'>');
            record_bytes.extend_from_slice(record.id());
            record_bytes.push(b'\n');
            record_bytes.extend_from_slice(&record.seq());
            record_bytes.push(b'\n');

            buffer.extend_from_slice(&record_bytes);
            record_count += 1;

            if buffer.len() >= effective_chunk {
                if tx.send(ParseOutput::Bytes(Bytes::from(std::mem::take(&mut buffer)))).await.is_err() {
                    return Err(anyhow!("Failed to send byte chunk at record {}", record_count));
                }
            }
        }

        if !buffer.is_empty() {
            if tx.send(ParseOutput::Bytes(Bytes::from(buffer))).await.is_err() {
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
/// Writes any FASTA stream to file with aggressive batching + large buffer.
/// Used in 9 places — keeps the exact same API so nothing else changes.
pub fn write_fasta_stream_to_file(
    stream: ReceiverStream<ParseOutput>,
    dest_path: PathBuf,
    config: Arc<RunConfig>,           // ← Arc for 'static spawn
    data_type: StreamDataType,        // ← added
    label: &str,                      // optional label for logging
) -> JoinHandle<Result<()>> {

    let dest_path_clone = dest_path.clone();
    let label_owned = label.to_string();

    tokio::spawn(async move {
        let file = TokioFile::create(&dest_path_clone)
            .await
            .map_err(|e| anyhow!("Cannot create FASTA file {}: {}", dest_path_clone.display(), e))?;

        // Dynamic hot buffer
        let effective_buffer = crate::utils::system::compute_buffer_size(
            &config,
            "write_fasta_stream_to_file",
            data_type,
            2.0,                          // very hot write path
        );

        let mut writer = tokio::io::BufWriter::with_capacity(effective_buffer, file);

        let mut stream = stream;
        let mut batch: Vec<u8> = Vec::with_capacity(effective_buffer / 4); // 4:1 batch ratio

        while let Some(item) = stream.next().await {
            let bytes = item
                .to_bytes()
                .map_err(|e| anyhow!("Failed to convert record to bytes: {}", e))?;

            batch.extend_from_slice(&bytes);

            if batch.len() >= effective_buffer / 4 {
                writer.write_all(&batch).await
                    .map_err(|e| anyhow!("Write error to {}: {}", dest_path_clone.display(), e))?;
                batch.clear();
            }
        }

        if !batch.is_empty() {
            writer.write_all(&batch).await
                .map_err(|e| anyhow!("Final write error to {}: {}", dest_path_clone.display(), e))?;
        }

        writer.flush().await
            .map_err(|e| anyhow!("Flush error to {}: {}", dest_path_clone.display(), e))?;

        let final_size = writer.stream_position().await.unwrap_or(0);
        info!("{} written to {} ({} bytes, buffer {} MiB)",
              label_owned, dest_path_clone.display(), final_size, effective_buffer / (1024*1024));

        Ok(())
    })
}

fn compare_read_ids_bytes(id1: &[u8], id2: &[u8]) -> bool {
    // Fast path: exact byte match (very common when IDs are identical)
    if id1 == id2 {
        return true;
    }

    let mut base1 = id1;
    let mut base2 = id2;

    // Strip common suffixes (byte-level)
    let suffixes: [&[u8]; 24] = [
        b"/1", b"/2",
        b".1", b".2",
        b" 1:", b" 2:",
        b" 1:N:0:", b" 2:N:0:",
        b"/1:N:0:", b"/2:N:0:",
        b"_1", b"_2",
        b" 1#0/1", b" 2#0/2",
        b"|1", b"|2",
        b" 1#", b" 2#",
        b"/R1", b"/R2",
        b"_R1", b"_R2",
        b" 1:0", b" 2:0",
    ];

    let mut stripped1 = false;
    for &suffix in &suffixes {
        if base1.ends_with(suffix) {
            base1 = &base1[..base1.len() - suffix.len()];
            stripped1 = true;
            break;
        }
    }

    let mut stripped2 = false;
    for &suffix in &suffixes {
        if base2.ends_with(suffix) {
            base2 = &base2[..base2.len() - suffix.len()];
            stripped2 = true;
            break;
        }
    }

    if (stripped1 || stripped2) && base1 == base2 {
        return true;
    }

    // Try a more aggressive approach: strip everything after the first space
    if let Some(space_idx) = base1.iter().position(|&b| b == b' ') {
        base1 = &base1[..space_idx];
    }
    if let Some(space_idx) = base2.iter().position(|&b| b == b' ') {
        base2 = &base2[..space_idx];
    }

    // Clean up any trailing characters
    while !base1.is_empty() {
        let last = base1[base1.len() - 1];
        if last == b' ' || last == b'\t' || last == b':' || last == b'#' || last == b'.' || last == b'/' || last == b'_' {
            base1 = &base1[..base1.len() - 1];
        } else {
            break;
        }
    }
    while !base2.is_empty() {
        let last = base2[base2.len() - 1];
        if last == b' ' || last == b'\t' || last == b':' || last == b'#' || last == b'.' || last == b'/' || last == b'_' {
            base2 = &base2[..base2.len() - 1];
        } else {
            break;
        }
    }

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
                    seq: Bytes::from(seq.into_bytes()),
                    qual: Bytes::from(qual.into_bytes()),
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
    compare_read_ids_bytes(id1.as_bytes(), id2.as_bytes())
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
        let fasta_line = format!(">{}\n{}\n", seq_record.id(), String::from_utf8_lossy(&seq_record.seq()));
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
                                seq: Bytes::from(new_seq),
                                qual: Bytes::from(new_qual),
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
                            let new_rec = SequenceRecord::Fasta {
                                id: new_id,
                                desc: new_desc,
                                seq: rec.seq(),   // ← cheap Bytes clone (just refcount bump)
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

        // Process unidentified contigs → conform headers (parallel + buffered)
        {
            let unidentified_concurrency = compute_phase_concurrency(
                &config,
                "unidentified_taxid_conform",
                0.5,
                4.0,
                128,
                16,
            );
            let unid_batch_size = 32_000;


            let (unid_batch_tx, unid_batch_rx) = mpsc::channel::<Vec<SequenceRecord>>(unidentified_concurrency * 3);
            let unid_batch_rx = Arc::new(tokio::sync::Mutex::new(unid_batch_rx));

            let mut unid_workers = Vec::with_capacity(unidentified_concurrency);
            for _ in 0..unidentified_concurrency {
                let rx = unid_batch_rx.clone();
                let c_tx = combined_tx.clone();

                unid_workers.push(tokio::spawn(async move {
                    loop {
                        let batch_opt = {
                            let mut guard = rx.lock().await;
                            guard.recv().await
                        };

                        let batch = match batch_opt {
                            Some(b) => b,
                            None => break,
                        };

                        for item in batch.iter() {
                            let rec: &SequenceRecord = item;

                            let conformed = format!("{}{}", CONFORMING_PREAMBLE, rec.id());
                            let (new_id, new_desc) = parse_header(conformed.as_bytes(), '>');

                            let new_rec = SequenceRecord::Fasta {
                                id: new_id,
                                desc: new_desc,
                                seq: rec.seq(),
                            };


                            let _ = c_tx.send(ParseOutput::Fasta(new_rec)).await;
                        }
                    }

                    Ok::<(), anyhow::Error>(())
                }));
            }

            // Producer: batch incoming unidentified FASTA records
            let mut unid_stream = unidentified_contigs_stream;
            let mut current_batch: Vec<SequenceRecord> = Vec::with_capacity(unid_batch_size);

            while let Some(item) = unid_stream.next().await {
                match item {
                    ParseOutput::Fasta(rec) => {
                        current_batch.push(rec);
                        if current_batch.len() >= unid_batch_size {
                            let batch_to_send = std::mem::take(&mut current_batch);
                            let _ = unid_batch_tx.send(batch_to_send).await;
                        }
                    }
                    _ => {
                        warn!("Unexpected non-Fasta item in unidentified_contigs_stream");
                    }
                }
            }

            if !current_batch.is_empty() {
                let _ = unid_batch_tx.send(current_batch).await;
            }

            drop(unid_batch_tx);

            // Wait for all workers
            for worker in unid_workers {
                worker.await??;
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
    config: Arc<RunConfig>,
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
    let concurrency = compute_phase_concurrency(
        &config,
        "generate_taxid_locator",
        0.5,
        4.0,
        128,
        16,
    );
    info!("generate_taxid_locator starting — {} output workers", concurrency);

    let levels = vec!["species", "genus", "family"];
    let hit_types = vec!["NT", "NR"];

    let mut output_paths = Vec::new();
    let mut tasks: Vec<JoinHandle<Result<()>>> = Vec::new();
    let mut record_channels = Vec::new(); // (taxid_field, sender)

    // Create one channel + one task per output type (8 total)
    for level in levels.iter() {
        for hit_type in hit_types.iter() {
            let level_owned = level.to_string();
            let hit_type_owned = hit_type.to_string();
            let taxid_field = format!("{}_{}", level_owned, hit_type_owned.to_lowercase());

            let fa_name = if level_owned == "species" {
                format!("refined_taxid_annot_sorted_{}.fasta", hit_type_owned.to_lowercase())
            } else {
                format!(
                    "refined_taxid_annot_sorted_{}_{}.fasta",
                    level_owned,
                    hit_type_owned.to_lowercase()
                )
            };
            let fa_path = assembly_dir.join(&fa_name);
            let json_path = assembly_dir.join(format!(
                "refined_taxid_locations_{}_{}.json",
                level_owned,
                hit_type_owned.to_lowercase()
            ));

            output_paths.push(fa_path.clone());
            output_paths.push(json_path.clone());

            let (tx, rx) = mpsc::channel::<(String, String)>(32768);

            // Clone what the worker needs before moving into closure
            let taxid_field_for_worker = taxid_field.clone();
            let hit_type_for_worker = hit_type_owned.clone();
            let fa_path_for_worker = fa_path.clone();
            let json_path_for_worker = json_path.clone();

            let task = tokio::spawn(async move {
                let mut records_by_taxid: AHashMap<i32, Vec<(String, String)>> = AHashMap::new();

                let mut rx = rx;
                while let Some((header, seq)) = rx.recv().await {
                    let taxid = get_taxid(&header, &taxid_field_for_worker).unwrap_or(-1);
                    records_by_taxid.entry(taxid).or_default().push((header, seq));
                }

                if records_by_taxid.is_empty() {
                    std::fs::write(&fa_path_for_worker, b"")?;
                    std::fs::write(&json_path_for_worker, b"[]")?;
                    debug!("Empty output written: {}", fa_path_for_worker.display());
                    return Ok(());
                }

                let mut sorted_taxids: Vec<i32> = records_by_taxid.keys().copied().collect();
                sorted_taxids.sort();

                let file = std::fs::File::create(&fa_path_for_worker)
                    .map_err(|e| anyhow!("create {}: {}", fa_path_for_worker.display(), e))?;
                let mut writer = std::io::BufWriter::with_capacity(32 * 1024 * 1024, file);

                let mut pos: u64 = 0;
                let mut locations = Vec::with_capacity(records_by_taxid.len());

                for taxid in sorted_taxids {
                    let recs = records_by_taxid.get(&taxid).unwrap();
                    let start = pos;

                    for (header, seq) in recs {
                        writer.write_all(b">")?;
                        writer.write_all(header.as_bytes())?;
                        writer.write_all(b"\n")?;
                        pos += 1 + header.len() as u64 + 1;

                        writer.write_all(seq.as_bytes())?;
                        writer.write_all(b"\n")?;
                        pos += seq.len() as u64 + 1;
                    }

                    locations.push(TaxonSeqLocation {
                        taxid,
                        first_byte: start,
                        last_byte: pos.saturating_sub(1),
                        hit_type: hit_type_for_worker.clone(),
                    });
                }

                writer.flush()
                    .map_err(|e| anyhow!("flush {}: {}", fa_path_for_worker.display(), e))?;

                let json_bytes = serde_json::to_vec(&locations)?;
                std::fs::write(&json_path_for_worker, json_bytes)?;

                info!(
                    "Wrote {} ({} taxids, {} records)",
                    fa_path_for_worker.display(),
                    locations.len(),
                    records_by_taxid.values().map(|v| v.len()).sum::<usize>()
                );
                Ok(())
            });

            tasks.push(task);
            record_channels.push((taxid_field, tx));
        }
    }

    // Stream parsing + selective broadcast
    let mut fasta_stream = ReceiverStream::new(fasta_stream_rx);
    while let Some(item) = fasta_stream.next().await {
        if let ParseOutput::Fasta(rec) = item {
            let header = rec.id().to_string();
            let seq = String::from_utf8(rec.seq().to_vec())
                .map_err(|e| anyhow!("Invalid UTF-8 in sequence: {}", e))?;

            // Send only to channels where this header has a valid taxid for that field
            for (field, tx) in &record_channels {
                if get_taxid(&header, field).unwrap_or(-1) != -1 {
                    let _ = tx.send((header.clone(), seq.clone())).await;
                }
            }
        }
    }

    // Close all channels
    for (_, tx) in record_channels {
        drop(tx);
    }

    // Wait for all workers
    for task in tasks {
        task.await??;
    }

    // Combine JSON locators
    let json_paths: Vec<PathBuf> = output_paths
        .iter()
        .filter(|p| p.extension().map_or(false, |e| e == "json"))
        .cloned()
        .collect();

    let combined_path = assembly_dir.join("refined_taxid_locations_combined.json");
    combine_taxon_loc_json(json_paths, combined_path.clone())?;

    output_paths.push(combined_path);

    info!("generate_taxid_locator complete");
    Ok((output_paths, vec![], vec![]))
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
                            let _ = tx.send(ParseOutput::Bytes(vec)).await;
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

pub async fn write_combined_fastq(
    r1_path: PathBuf,
    r2_path_opt: Option<PathBuf>,
    combined_path: &PathBuf,
) -> Result<(), PipelineError> {
    use tokio::io::AsyncWriteExt;

    let mut out = TokioFile::create(combined_path)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    let mut r1 = TokioFile::open(&r1_path)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    tokio::io::copy(&mut r1, &mut out)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    if let Some(r2_path) = r2_path_opt {
        let mut r2 = TokioFile::open(&r2_path)
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        tokio::io::copy(&mut r2, &mut out)
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
    }

    out.flush()
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use futures::StreamExt;
    use std::fs;
    use std::io;
    use std::path::PathBuf;
    use std::sync::Arc;
    use tempfile::tempdir;

    fn compare_header(head: &[u8], prefix: char) {
        let scalar = parse_header_scalar(head, prefix);
        let dispatched = parse_header(head, prefix);
        assert_eq!(scalar.0, dispatched.0, "id mismatch for {:?}", std::str::from_utf8(head));
        assert_eq!(scalar.1, dispatched.1, "desc mismatch for {:?}", std::str::from_utf8(head));
    }

    #[test]
    fn test_header_id_only_fasta() {
        let head = b"NC_045512.2";
        let (id, desc) = parse_header_scalar(head, '>');
        assert_eq!(id, "NC_045512.2");
        assert!(desc.is_none());
        compare_header(head, '>');
    }

    #[test]
    fn test_header_with_space_desc() {
        let head = b"NC_045512.2 Severe acute respiratory syndrome coronavirus 2";
        let (id, desc) = parse_header_scalar(head, '>');
        assert_eq!(id, "NC_045512.2");
        assert_eq!(desc.as_deref(), Some("Severe acute respiratory syndrome coronavirus 2"));
        compare_header(head, '>');
    }

    #[test]
    fn test_header_with_tab_desc() {
        let head = b"read1\tsome description";
        let (id, desc) = parse_header_scalar(head, '@');
        assert_eq!(id, "read1");
        assert_eq!(desc.as_deref(), Some("some description"));
        compare_header(head, '@');
    }

    #[test]
    fn test_header_fastq_prefix_stripped() {
        let head = b"@SRR12345.1 length=150";
        let (id, desc) = parse_header_scalar(head, '@');
        assert_eq!(id, "SRR12345.1");
        assert_eq!(desc.as_deref(), Some("length=150"));
        compare_header(head, '@');
    }

    #[test]
    fn test_header_empty_desc_after_space() {
        let head = b"contig_00001 ";
        let (id, desc) = parse_header_scalar(head, '>');
        assert_eq!(id, "contig_00001");
        assert!(desc.is_none(), "trailing space with empty desc should give None");
        compare_header(head, '>');
    }

    #[test]
    fn test_header_long_id_spans_two_simd_chunks() {
        let long_id = "a".repeat(70);
        let head = format!("{} some description", long_id);
        let (id, desc) = parse_header_scalar(head.as_bytes(), '>');
        assert_eq!(id, long_id);
        assert_eq!(desc.as_deref(), Some("some description"));
        compare_header(head.as_bytes(), '>');
    }

    #[test]
    fn test_header_long_id_only_no_whitespace() {
        let long_id = "a".repeat(100);
        let (id, desc) = parse_header_scalar(long_id.as_bytes(), '>');
        assert_eq!(id, long_id);
        assert!(desc.is_none());
        compare_header(long_id.as_bytes(), '>');
    }

    #[test]
    fn test_header_sra_style() {
        let head = b"SRR12345.1.1 SRR12345.1 length=150";
        let (id, desc) = parse_header_scalar(head, '@');
        assert_eq!(id, "SRR12345.1.1");
        assert_eq!(desc.as_deref(), Some("SRR12345.1 length=150"));
        compare_header(head, '@');
    }

    #[test]
    fn test_compare_read_ids_exact_and_suffix_variants() {
        assert!(compare_read_ids("read1", "read1"));
        assert!(compare_read_ids("read1/1", "read1/2"));
        assert!(compare_read_ids("read1_R1", "read1_R2"));
        assert!(compare_read_ids("read1.1", "read1.2"));
        assert!(compare_read_ids("read1 1:N:0:", "read1 2:N:0:"));
        assert!(!compare_read_ids("read1", "read2"));
        assert!(!compare_read_ids("read1/1", "read2/2"));
    }

    #[test]
    fn test_sequence_record_accessors() {
        let fasta = SequenceRecord::Fasta {
            id: "contig1".to_string(),
            desc: Some("some desc".to_string()),
            seq: Bytes::from_static(b"ACGT"),
        };

        assert_eq!(fasta.id(), "contig1");
        assert_eq!(fasta.seq().as_ref(), b"ACGT");
        assert_eq!(fasta.qual().as_ref(), b"");
        assert_eq!(fasta.desc(), Some("some desc"));

        let fasta_clone = fasta.clone();
        assert_eq!(fasta_clone.seq().as_ref(), b"ACGT");

        let fastq = SequenceRecord::Fastq {
            id: "read1".to_string(),
            desc: None,
            seq: Bytes::from_static(b"TGCA"),
            qual: Bytes::from_static(b"IIII"),
        };

        assert_eq!(fastq.id(), "read1");
        assert_eq!(fastq.seq().as_ref(), b"TGCA");
        assert_eq!(fastq.qual().as_ref(), b"IIII");
        assert_eq!(fastq.desc(), None);

        let fastq_clone = fastq.clone();
        assert_eq!(fastq_clone.seq().as_ref(), b"TGCA");
    }

    #[test]
    fn test_validate_sequence() {
        assert!(validate_sequence(b"ACGTACGT", b"ACGT").is_ok());

        let err = validate_sequence(b"ACNT", b"ACGT").unwrap_err();
        assert!(err.to_string().contains("Invalid nucleotide"));
    }

    #[tokio::test]
    async fn test_validate_sequence_parallel() -> Result<()> {
        validate_sequence_parallel(Bytes::from_static(b"ACGTACGTACGT"), b"ACGT", 3).await?;

        let err = validate_sequence_parallel(Bytes::from_static(b"ACGTXCGT"), b"ACGT", 3)
            .await
            .unwrap_err();
        assert!(err.to_string().contains("Invalid nucleotide"));
        Ok(())
    }

    #[tokio::test]
    async fn test_fastx_generator_yields_expected_records() -> Result<()> {
        let mut stream = fastx_generator(4, 12, 30.0, 5.0);
        let mut seen = 0usize;

        while let Some(record) = stream.next().await {
            seen += 1;
            assert_eq!(record.id(), format!("read{}", seen));
            assert_eq!(record.seq().len(), 12);
            assert_eq!(record.qual().len(), 12);
            assert!(matches!(record, SequenceRecord::Fastq { .. }));
        }

        assert_eq!(seen, 4);
        Ok(())
    }

    #[tokio::test]
    async fn test_fastx_generator_zero_length_is_empty() -> Result<()> {
        let mut stream = fastx_generator(10, 0, 30.0, 5.0);
        assert!(stream.next().await.is_none());
        Ok(())
    }

    #[test]
    fn test_r1r2_base_detects_r1_and_builds_filename() {
        let path = PathBuf::from("/tmp/sample_R1.fastq");
        let result = r1r2_base(&path);

        assert_eq!(result.delimiter, Some('_'));
        assert_eq!(result.r1_tag.as_deref(), Some("R1"));
        assert_eq!(result.prefix.as_deref(), Some("sample"));
        assert_eq!(result.index, Some(1));
        assert_eq!(result.file_name.as_deref(), Some("sample.fastq"));
    }

    #[test]
    fn test_r1r2_base_non_matching_path_returns_empty_result() {
        let path = PathBuf::from("/tmp/sample.fastq");
        let result = r1r2_base(&path);

        assert_eq!(
            result,
            R1R2Result {
                delimiter: None,
                r1_tag: None,
                file_name: None,
                index: None,
                prefix: None,
            }
        );
    }

    #[test]
    fn test_sequence_reader_fasta() -> io::Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("sample.fasta");
        fs::write(&path, b">seq1 testFASTA\nATCG\n")?;

        match sequence_reader(&path)? {
            SequenceReader::Fasta(_) => Ok(()),
            _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fasta reader")),
        }
    }

    #[test]
    fn test_sequence_reader_fastq() -> io::Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("sample.fastq");
        fs::write(&path, b"@seq1\nATCG\n+\nIIII\n")?;

        match sequence_reader(&path)? {
            SequenceReader::Fastq(_) => Ok(()),
            _ => Err(io::Error::new(io::ErrorKind::Other, "Expected Fastq reader")),
        }
    }

    #[test]
    fn test_sequence_reader_invalid_extension() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("sample.txt");
        fs::write(&path, b"@seq1\nATCG\n+\nIIII\n").unwrap();

        match sequence_reader(&path) {
            Ok(_) => panic!("expected invalid extension error"),
            Err(err) => assert_eq!(err.kind(), io::ErrorKind::InvalidData),
        }
    }

    #[tokio::test]
    async fn test_raw_read_count_single_end() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("single.fastq");
        fs::write(&path, b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nHHHH\n")?;

        let handle = raw_read_count(path, None);
        let count = handle.await??;
        assert_eq!(count, 2);
        Ok(())
    }

    #[tokio::test]
    async fn test_raw_read_count_paired_end() -> Result<()> {
        let dir = tempdir()?;
        let r1 = dir.path().join("r1.fastq");
        let r2 = dir.path().join("r2.fastq");

        fs::write(
            &r1,
            b"@read1/1\nACGT\n+\nIIII\n@read2/1\nTGCA\n+\nHHHH\n",
        )?;
        fs::write(
            &r2,
            b"@read1/2\nACGT\n+\nIIII\n@read2/2\nTGCA\n+\nHHHH\n",
        )?;

        let handle = raw_read_count(r1, Some(r2));
        let count = handle.await??;
        assert_eq!(count, 4);
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_record_counter_fastq() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let fastq_data = vec![
            ParseOutput::Bytes(b"@read1\nATCG\n+\nIIII\n".to_vec().into()),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: Bytes::from_static(b"GCTA"),
                qual: Bytes::from_static(b"HHHH"),
            }),
        ];

        tokio::spawn(async move {
            for item in fastq_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_record_counter(rx, false).await?;
        assert_eq!(count, 2);
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
                seq: Bytes::from_static(b"GCTA"),
            }),
        ];

        tokio::spawn(async move {
            for item in fasta_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_record_counter(rx, true).await?;
        assert_eq!(count, 1);
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_record_counter_empty() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        drop(tx);
        let count = stream_record_counter(rx, false).await?;
        assert_eq!(count, 0);
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_and_filter_fastq_id() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);

        tokio::spawn(async move {
            tx.send(ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "keep_me_01".to_string(),
                desc: None,
                seq: Bytes::from_static(b"ATCG"),
                qual: Bytes::from_static(b"IIII"),
            }))
                .await
                .unwrap();

            tx.send(ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "ignore_me_02".to_string(),
                desc: None,
                seq: Bytes::from_static(b"TGCA"),
                qual: Bytes::from_static(b"HHHH"),
            }))
                .await
                .unwrap();
        });

        let (filtered_rx, task) = parse_and_filter_fastq_id(rx, 10, "keep_me".to_string());
        let mut filtered_rx = filtered_rx;
        let mut results = Vec::new();

        while let Some(record) = filtered_rx.recv().await {
            results.push(record);
        }

        task.await??;

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].id(), "keep_me_01");
        assert_eq!(results[0].seq().as_ref(), b"ATCG");
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_byte_stream_to_fastq_round_trip() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        tx.send(ParseOutput::Bytes(
            b"@read1\nACGT\n+\nIIII\n".to_vec().into(),
        ))
            .await
            .unwrap();
        drop(tx);

        let (out_rx, task) = parse_byte_stream_to_fastq(rx, 10, 1).await?;
        let mut out_stream = ReceiverStream::new(out_rx);
        let mut records = Vec::new();

        while let Some(item) = out_stream.next().await {
            match item {
                ParseOutput::Fastq(record) => records.push(record),
                other => panic!("unexpected output item: {:?}", other),
            }
        }

        task.await??;

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id(), "read1");
        assert_eq!(records[0].seq().as_ref(), b"ACGT");
        assert_eq!(records[0].qual().as_ref(), b"IIII");
        Ok(())
    }

    #[tokio::test]
    async fn test_concatenate_paired_reads() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/1".to_string(),
                desc: None,
                seq: Bytes::from_static(b"ATCG"),
                qual: Bytes::from_static(b"IIII"),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/2".to_string(),
                desc: None,
                seq: Bytes::from_static(b"GCTA"),
                qual: Bytes::from_static(b"HHHH"),
            }),
        ];

        tokio::spawn(async move {
            for record in records {
                tx.send(record).await.unwrap();
            }
        });

        let (concat_rx, task) = concatenate_paired_reads(ReceiverStream::new(rx), 10, 10).await?;
        let mut concat_records = Vec::new();
        let mut stream = concat_rx;

        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
            concat_records.push(record);
        }

        task.await??;

        assert_eq!(concat_records.len(), 1);
        assert_eq!(concat_records[0].id(), "read1/1");
        assert_eq!(concat_records[0].seq().as_ref(), b"ATCGNGCTA");
        assert_eq!(concat_records[0].qual().as_ref(), b"IIII!HHHH");
        Ok(())
    }

    #[test]
    fn test_parse_header_avx_vs_scalar_equivalence() {
        use proptest::prelude::*;

        proptest!(|(s in "\\PC*", prefix in "[>@]")| {
            let bytes = s.as_bytes();
            let p = prefix.chars().next().unwrap();
            let scalar = parse_header_scalar(bytes, p);
            let dispatched = parse_header(bytes, p);

            assert_eq!(scalar, dispatched, "Mismatch on: {:?} with prefix {}", s, p);
        });
    }

    #[test]
    fn test_parse_header_equivalence() {
        let long_header = b"very_long_id_that_spans_multiple_chunks_".repeat(4);

        let cases: Vec<&[u8]> = vec![
            b"seq1 first record",
            b"seq1\tfirst record",
            b"seq1\nfirst record",
            b"seq1 first record with spaces",
            b"seq1 ",
            b"seq1",
            b"very_long_id_that_spans_multiple_chunks_",
            &long_header[..],
        ];

        for &header in &cases {
            let scalar = parse_header_scalar(header, '>');
            let actual = parse_header(header, '>');

            assert_eq!(
                scalar, actual,
                "Header mismatch on: {:?}\n  scalar: {:?}\n  actual: {:?}",
                std::str::from_utf8(header),
                scalar,
                actual
            );
        }
    }

    #[test]
    fn test_parse_header_long_and_edge_cases() {
        let long_id = "a".repeat(250);

        let test_cases: Vec<(&[u8], char)> = vec![
            (long_id.as_bytes(), '>'),
            (b"contig_00001 ", '>'),
            (b"SRR12345.1 length=150", '@'),
            (b"read1\tsome description with\ttabs", '@'),
            (b"header_with_trailing_space ", '>'),
        ];

        for (header, prefix) in test_cases {
            let scalar = parse_header_scalar(header, prefix);
            let actual = parse_header(header, prefix);

            assert_eq!(
                scalar, actual,
                "Failed on header: {:?} (prefix '{}')",
                std::str::from_utf8(header),
                prefix
            );
        }
    }
}