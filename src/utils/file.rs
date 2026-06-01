use std::fs::File;
use std::io;
use std::sync::Arc;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use log::{self, debug, info, warn};
use anyhow::{Result, anyhow};
use tempfile::TempDir;
use tokio::io::{AsyncWriteExt, BufWriter};
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::streams::ParseOutput;
use tokio::task::JoinHandle;
use tokio::fs::File as TokioFile;
use tokio::process::Command;
use tokio::io::AsyncSeekExt;
use tokio_stream::StreamExt;
use tempfile::Builder as TempfileBuilder;
use crate::utils::streams::ToBytes;
use crate::utils::fastx::SequenceRecord;
use crate::config::defs::{RunConfig, PipelineError, StreamDataType};
use sysinfo::Disks;
use nix::sys::statvfs;
use regex::Regex;





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

pub fn is_gzipped(path: &PathBuf) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut buffer = [0u8; 2];
    file.read_exact(&mut buffer)?;
    Ok(buffer == [0x1F, 0x8B]) // Gzip magic bytes
}


/// Absolut path resolver.
/// # Arguments
///
/// * `path`: &str path
/// * 'basedir': canoncially th cwd, but any absolute path dir
///
/// # Returns
/// PathBuf:absolute path
pub fn resolve_to_absolute(path: &str, base_dir: &Path) -> PathBuf {
    let p = PathBuf::from(path);
    let abs = if p.is_relative() {
        base_dir.join(p)
    } else {
        p
    };
    abs.canonicalize().unwrap_or(abs)
}


/// Uses resolve_to_absolute to resolve optionmal file paths.
/// # Arguments
///
/// * `cli`optiuonal cli arg string path
/// * 'basedir': passes to resolve_to_absolute
///
/// # Returns
/// Result of optional resolved PAthBuf
pub fn resolve_optional_path(
    cli: &Option<String>,
    base_dir: &Path,
) -> Result<Option<PathBuf>> {
    match cli {
        Some(s) => {
            let abs = resolve_to_absolute(s, base_dir);
            if !abs.exists() {
                return Err(anyhow!("File not found: {}", abs.display()));
            }
            if !abs.is_file() {
                return Err(anyhow!("Not a file: {}", abs.display()));
            }
            Ok(Some(abs))
        }
        None => Ok(None),
    }
}

/// Allows pre and post-fixes to be appended to a base file name whle preserving its dir.
/// # Arguments
///
/// * `path`: &PathBuf
/// * 'prefix': Option<&str> = added ahead of base
/// * 'postfix': Option<&str> = added after of base
///  * 'delimiter': &str = added between prefix, postfix and base.
///
/// # Returns
/// PathBuf
pub fn rename_file_path(
    path: &Path,
    prefix: Option<&str>,
    postfix: Option<&str>,
    delimiter: &str,
) -> PathBuf {
    let (stem, extensions) = extension_remover(path);
    let base = stem.file_name().and_then(|s| s.to_str()).unwrap_or("");

    let new_base = match (prefix, postfix) {
        (Some(p), Some(q)) => format!("{p}{delimiter}{base}{delimiter}{q}"),
        (Some(p), None) => format!("{p}{delimiter}{base}"),
        (None, Some(q)) => format!("{base}{delimiter}{q}"),
        (None, None) => base.to_string(),
    };

    let new_name = if extensions.is_empty() {
        new_base
    } else {
        format!("{}.{}", new_base, extensions.join("."))
    };

    PathBuf::from(new_name)
}

/// Calls file_name_manipulator to make alterations to the file name.
/// Then returns absolute path.
/// # Arguments
///
/// * `path`: PathBuf - File path
/// * 'prefix': Option<&str> = added ahead of base
/// * 'postfix': Option<&str> = added after of base
/// /// * 'postfix': Option<&str> = added after of base
// /// * 'delimiter': &str = added between prefix, postfix and base.
///
/// # Returns
/// PathBuf: modified, absolute path
pub fn file_path_manipulator(path: &PathBuf, parent_dir: Option<&PathBuf>, prefix: Option<&str>, postfix: Option<&str>, delimiter: &str) -> PathBuf {
    let resolved_path = if path.is_absolute() {
        path
    } else {
        let parent_path = parent_dir.expect("parent_dir must not be None when path is relative").canonicalize().ok().expect("Parent directory not found. {parent_dir:?}");
        &parent_path.join(path)
    };
    let absolute_file_name = file_name_manipulator(resolved_path, prefix, postfix, delimiter);
    let absolute_file_path = PathBuf::from(absolute_file_name);

    absolute_file_path
}


/// Calls extension_remover to retrieve base and extensions.
/// Then prepends prefix and appends postfix, if any.

/// # Arguments
///
/// * `path`: PathBuf - File path
/// * 'prefix': Option<&str> = added ahead of base
/// * 'postfix': Option<&str> = added after of base
/// * 'delimiter': &str = added between prefix, postfix and base.
///
/// # Returns
/// String: new file name
pub fn file_name_manipulator(path: &PathBuf, prefix: Option<&str>, postfix:Option<&str>, delimiter: &str) -> String {

    let (stem, extensions) = extension_remover(&path);
    let base = stem.to_str().unwrap_or("");
    let new_base = match (prefix, postfix) {
        (Some(p), Some(q)) => format!("{}{}{}{}{}", p, delimiter, base, delimiter, q),
        (Some(p), None) => format!("{}{}{}", p, delimiter, base),
        (None, Some(q)) => format!("{}{}{}", base, delimiter, q),
        (None, None) => base.to_string(),
    };

    let new_file_name = if extensions.is_empty() {
        new_base
    } else {
        format!("{}.{}", new_base, extensions.join("."))
    };

    new_file_name
}


/// Strips either one extension, or two if the last one is gz
/// # Arguments
///
/// * `path`: PathBuf - File path
///
/// # Returns
/// PathBuf of stripped file, extensions.
pub fn extension_remover(path: &Path) -> (PathBuf, Vec<String>) {
    let mut extensions = Vec::new();
    let mut current = path;

    while let Some(ext_os) = current.extension() {
        let ext = ext_os.to_str().expect("non-UTF8 extension");
        extensions.push(ext.to_string());

        if extensions.len() > 2 { break; }
        if extensions.len() == 2 && extensions[1] != "gz" {
            extensions.pop();
            break;
        }

        current = current.file_stem().map(Path::new).unwrap_or(Path::new(""));
    }

    extensions.reverse();

    let stem_path = if let Some(parent) = path.parent() {
        parent.join(current)
    } else {
        current.to_path_buf()
    };

    (stem_path, extensions)
}

/// Checks if a file has one of a set of extensions
/// # Arguments
///
/// * `extensions`: from extension_remover
/// * 'valid_extensions': &[&str] = extensions (from lazy_static)
///
/// # Returns
/// Bool
pub fn has_any_extension_from_path(extensions: &[String], valid_extensions: &[&str]) -> bool {
    extensions
        .iter()
        .any(|ext| valid_extensions.iter().any(|&valid| ext.eq_ignore_ascii_case(valid)))
}

/// Scan a directory for files with extensions matching
/// # Arguments
///
/// * `dir`: PathBuf - Directory path
/// * 'extensions': &[&str] = extensions (from lazy_static)
///
/// # Returns
/// Result<Vev<PathBuf>>
pub fn scan_files_with_extensions(dir: &PathBuf, valid_extensions: &[&str]) -> Result<Vec<PathBuf>> {
    if !dir.is_dir() {
        return Err(anyhow!("Provided path is not a directory: {}", dir.display()));
    }

    let mut matching_files = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        let (_, extensions) = extension_remover(&path);

        if path.is_file() && has_any_extension_from_path(&extensions, valid_extensions) {
            let full_path = file_path_manipulator(&path, Some(dir), None, None, "");
            matching_files.push(full_path);
        }
    }

    if matching_files.is_empty() {
        return Err(anyhow!("No files with extensions {:?} found in directory: {}", valid_extensions, dir.display()));
    }

    Ok(matching_files)
}


/// Creates a named FIFO pipe and writes data from a ParseOutput stream to it asynchronously.
///
/// # Arguments
///
/// - `fifo_path`: Path to the FIFO pipe to create.
/// - `input_stream`: A `ReceiverStream` yielding `ParseOutput` items (expects `Bytes` variant).
/// - `buffer_size`: Optional buffer size for the writer (defaults to 4MB if not provided).
///
/// # Returns
/// A `Result` containing a `JoinHandle` that resolves to `Result<(), anyhow::Error>` upon completion.
pub async fn write_parse_output_to_fifo(
    fifo_path: &PathBuf,
    mut input_stream: ReceiverStream<ParseOutput>,
    buffer_size: Option<usize>,
) -> Result<JoinHandle<Result<(), anyhow::Error>>> {
    if fifo_path.exists() {
        tokio::fs::remove_file(fifo_path).await // Async remove
            .map_err(|e| anyhow!("Failed to remove existing FIFO at {}: {}", fifo_path.display(), e))?;
    }
    let status = Command::new("mkfifo")
        .arg(fifo_path)
        .status()
        .await
        .map_err(|e| anyhow!("Failed to execute mkfifo for {}: {}", fifo_path.display(), e))?;
    if !status.success() {
        return Err(anyhow!("mkfifo failed with status: {}", status));
    }
    let buffer_capacity = buffer_size.unwrap_or(16 * 1024 * 1024); // Default to 16MB for large files
    let fifo_path = fifo_path.clone();
    let task = tokio::spawn(async move {
        let writer_file = TokioFile::create(&fifo_path)
            .await
            .map_err(|e| anyhow!("Failed to open FIFO at {}: {}", fifo_path.display(), e))?;
        let mut writer = BufWriter::with_capacity(buffer_capacity, writer_file);
        let mut byte_count = 0;
        while let Some(item) = input_stream.next().await {
            match item {
                ParseOutput::Bytes(data) => {
                    writer
                        .write_all(&data)
                        .await
                        .map_err(|e| anyhow!("Failed to write to FIFO at {}: {}", fifo_path.display(), e))?;
                    byte_count += data.len();
                }
                ParseOutput::Fasta(record) => {
                    if let SequenceRecord::Fasta { id, desc, seq } = &record {
                        if let Some(d) = desc {
                            writer.write_all(format!(">{} {}\n", id, d).as_bytes()).await?;
                        } else {
                            writer.write_all(format!(">{}\n", id).as_bytes()).await?;
                        }
                        writer.write_all(&**seq).await?;
                        writer.write_all(b"\n").await?;
                        byte_count += seq.len() + id.len() + 2 + desc.as_ref().map_or(0, |d| d.len() + 1);
                    }
                }
                ParseOutput::Fastq(record) => {
                    if let SequenceRecord::Fastq { id, desc, seq, qual } = &record {
                        if let Some(d) = desc {
                            writer.write_all(format!("@{} {}\n", id, d).as_bytes()).await?;
                        } else {
                            writer.write_all(format!("@{}\n", id).as_bytes()).await?;
                        }
                        writer.write_all(&**seq).await?;
                        writer.write_all(b"\n+\n").await?;
                        writer.write_all(&**qual).await?;
                        writer.write_all(b"\n").await?;
                        byte_count += seq.len() + qual.len() + id.len() + 4 + desc.as_ref().map_or(0, |d| d.len() + 1);
                    }
                }
            }
        }
        writer
            .flush()
            .await
            .map_err(|e| anyhow!("Failed to flush FIFO at {}: {}", fifo_path.display(), e))?;
        writer
            .shutdown()
            .await
            .map_err(|e| anyhow!("Failed to shutdown FIFO at {}: {}", fifo_path.display(), e))?;
        if byte_count == 0 {
            return Err(anyhow!("No data written to FIFO at {}", fifo_path.display()));
        }
        Ok(())
    });
    Ok(task)
}


/// Creates a temporary named FIFO pipe and writes data from a ParseOutput stream to it asynchronously.
///
/// # Arguments
///
/// - `input_stream`: A `ReceiverStream` yielding `ParseOutput` items (expects `Bytes` variant).
/// - `buffer_size`: Optional buffer size for the writer (defaults to 4MB if not provided).
/// - `suffix`: Optional suffix for the temporary FIFO file name.
/// - `temp_dir`: Optional directory for the temporary FIFO (if provided, uses this directory; otherwise, uses system temp dir).
///
/// # Returns
///
/// A `Result` containing a `JoinHandle` that resolves to `Result<(), anyhow::Error>` upon completion.
pub async fn write_parse_output_to_temp_fifo(
    input_stream: ReceiverStream<ParseOutput>,
    buffer_size: Option<usize>,
    suffix: Option<&str>,
    temp_dir: Option<impl AsRef<Path>>,
) -> Result<(JoinHandle<Result<(), anyhow::Error>>, PathBuf)> {
    let mut builder = TempfileBuilder::new();
    if let Some(suf) = suffix {
        builder.suffix(suf);
    }
    let temp_name = if let Some(dir) = temp_dir {
        builder.tempfile_in(dir.as_ref())?
    } else {
        builder.tempfile()?
    };
    let temp_path = temp_name.path().to_path_buf();
    tokio::fs::remove_file(&temp_path).await?; // Async remove to replace with FIFO
    drop(temp_name); // Explicitly drop to avoid holding file

    let task = write_parse_output_to_fifo(&temp_path, input_stream, buffer_size.or(Some(16 * 1024 * 1024))).await?;

    Ok((task, temp_path))
}


/// Writes a Vec<u8> to a file asynchronously.
///
/// # Arguments
///
/// - `data`: The data as a byte vector (e.g., genomic sequence).
/// - `temp_path`: Path to the file (e.g., in /dev/shm).
/// - `buffer_size`: Buffer size for the writer (e.g., 4 MB).
///
/// # Returns
///
/// A `Result` containing a `JoinHandle` that resolves to `Result<(), anyhow::Error>` upon completion.
pub async fn write_vecu8_to_file<P: AsRef<Path>>(
    data: Arc<Vec<u8>>, // Accept Arc to avoid cloning
    temp_path: P,
    buffer_size: usize,
) -> Result<JoinHandle<Result<(), anyhow::Error>>> {
    let temp_path = temp_path.as_ref().to_path_buf();
    let buffer_capacity = buffer_size.max(16 * 1024 * 1024);
    let task = tokio::spawn(async move {
        let file = TokioFile::create(&temp_path)
            .await
            .map_err(|e| anyhow!("Failed to create file at {}: {}", temp_path.display(), e))?;
        let mut writer = BufWriter::with_capacity(buffer_capacity, file);

        let byte_count = if data.len() <= 4 * 1024 * 1024 { // Write small data in one go
            writer.write_all(&**data).await?;
            data.len()
        } else {
            const CHUNK_SIZE: usize = 16 * 1024 * 1024; // Larger chunks for big data
            let mut byte_count = 0;
            for chunk in data.chunks(CHUNK_SIZE) {
                writer.write_all(chunk).await?;
                byte_count += chunk.len();
            }
            byte_count
        };

        writer.flush().await?;
        writer.shutdown().await?;
        if byte_count == 0 {
            return Err(anyhow!("No data written to file at {}", temp_path.display()));
        }
        Ok(())
    });
    Ok(task)
}


/// Creates a tempfile in the specified directory with an optional suffix and writes data from a ParseOutput stream to it asynchronously.
///
/// # Arguments
///
/// - `input_steam`:`ReceiverStream` `ParseOutput` items.
/// - `buffer_size`:
/// - `suffix`: Optional
/// - `ram_temp_dir`: Directory for the temporary file (e.g., RAM disk like /dev/shm).
///
/// # Returns
/// A `Result` containing:
/// - A `JoinHandle` that resolves to `Result<(), anyhow::Error>` upon completion.
/// - The temporary file path as `PathBuf`.
/// - The `NamedTempFile` handle (hold this to prevent deletion; auto-deletes on drop unless persisted).
pub async fn write_parse_output_to_file(
    output_path: &PathBuf,
    mut input_stream: ReceiverStream<ParseOutput>,
    buffer_size: Option<usize>,
) -> Result<JoinHandle<Result<(), anyhow::Error>>> {
    let buffer_capacity = buffer_size.unwrap_or(16 * 1024 * 1024);
    let output_path_clone = output_path.clone();

    let task = tokio::spawn(async move {
        let file = TokioFile::create(&output_path_clone)
            .await
            .map_err(|e| anyhow!("Failed to create file at {}: {}", output_path_clone.display(), e))?;
        let mut writer = BufWriter::with_capacity(buffer_capacity, file);
        let mut total_bytes = 0u64;

        while let Some(item) = input_stream.next().await {
            let bytes = item.to_bytes()?;  // Convert ANY ParseOutput to bytes (uses ToBytes impl)
            writer.write_all(&bytes)
                .await
                .map_err(|e| anyhow!("Failed to write to {}: {}", output_path_clone.display(), e))?;
            total_bytes += bytes.len() as u64;
        }

        writer.flush()
            .await
            .map_err(|e| anyhow!("Failed to flush {}: {}", output_path_clone.display(), e))?;
        writer.shutdown()
            .await
            .map_err(|e| anyhow!("Failed to shutdown {}: {}", output_path_clone.display(), e))?;

        if total_bytes == 0 {
            warn!("No data written to file at {}", output_path_clone.display());
        } else {
            debug!("Wrote {} bytes to {}", total_bytes, output_path_clone.display());
        }

        Ok(())
    });

    Ok(task)
}

pub fn resolve_existing_input_path(input: &str, cwd: &Path) -> Result<PathBuf, PipelineError> {
    let resolved = resolve_to_absolute(input, cwd);

    if !resolved.exists() {
        return Err(PipelineError::FileNotFound(resolved));
    }
    if !resolved.is_file() {
        return Err(PipelineError::InvalidConfig(format!(
            "Not a file: {}",
            resolved.display()
        )));
    }

    Ok(resolved)
}

fn basename_without_extensions(path: &Path) -> Result<String, PipelineError> {
    let mut current = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| {
            PipelineError::InvalidConfig(format!(
                "Could not derive basename from {}",
                path.display()
            ))
        })?
        .to_string();

    loop {
        let p = Path::new(&current);
        match p.file_stem().and_then(|s| s.to_str()) {
            Some(stem) if stem != current => current = stem.to_string(),
            _ => break,
        }
    }

    Ok(current)
}

fn strip_common_read_suffixes(base: &str) -> String {
    // Ordered from most specific to least specific.
    // This covers common Illumina / SRA-like layouts:
    //   sample_S1_L001_R1_001
    //   sample_L001_R1_001
    //   sample_R1_001
    //   sample_R1
    //   sample_1
    //   sample.1
    //   sample-1
    let patterns = [
        r"(?i)^(?P<base>.+?)[._-]S\d+[._-]L\d{3}[._-]R[12][._-]\d{3}$",
        r"(?i)^(?P<base>.+?)[._-]L\d{3}[._-]R[12][._-]\d{3}$",
        r"(?i)^(?P<base>.+?)[._-]R[12][._-]\d{3}$",
        r"(?i)^(?P<base>.+?)[._-]R?[12]$",
    ];

    for pat in patterns {
        let re = Regex::new(pat).expect("valid mate-suffix regex");
        if let Some(caps) = re.captures(base) {
            return caps
                .name("base")
                .map(|m| m.as_str().trim_end_matches(&['_', '-', '.'][..]).to_string())
                .unwrap_or_else(|| base.to_string());
        }
    }

    base.to_string()
}

pub fn derive_sample_base_from_file1(file1_path: &Path) -> Result<PathBuf, PipelineError> {
    let basename = basename_without_extensions(file1_path)?;
    let stripped = strip_common_read_suffixes(&basename);
    Ok(PathBuf::from(stripped))
}

pub async fn validate_file_inputs(
    config: &RunConfig,
    cwd: &PathBuf,
) -> Result<(PathBuf, Option<PathBuf>, PathBuf, String), PipelineError> {
    let file1_path = match &config.args.file1 {
        Some(file) => resolve_existing_input_path(file, cwd)?,
        None => {
            return Err(PipelineError::InvalidConfig(
                "File1 path required".to_string(),
            ))
        }
    };

    let file2_path = match &config.args.file2 {
        Some(file) => {
            let resolved = resolve_existing_input_path(file, cwd)?;
            Some(resolved)
        }
        None => None,
    };

    let sample_base_buf = derive_sample_base_from_file1(&file1_path)?;
    let sample_base = sample_base_buf.to_string_lossy().into_owned();

    if let Some(ref file2) = file2_path {
        let file2_base = derive_sample_base_from_file1(file2)?;
        if file2_base != sample_base_buf {
            warn!(
                "file1/file2 sample bases do not match: '{}' vs '{}'",
                sample_base_buf.display(),
                file2_base.display()
            );
        }
    }

    Ok((file1_path, file2_path, sample_base_buf, sample_base))
}


/// Writes a byte-based ParseOutput stream to a regular file asynchronously.
///
/// # Arguments
/// - `output_path`: Path to the output file.
/// - `input_stream`: A `ReceiverStream` yielding `ParseOutput` items (expects `Bytes` variant).
/// - `buffer_size`: Optional buffer size for the writer (defaults to 16MB if not provided).
///
/// # Returns
/// A `JoinHandle` that resolves to `Result<(), anyhow::Error>` upon completion.
pub async fn write_byte_stream_to_file(
    dest_path: &PathBuf,
    stream: ReceiverStream<ParseOutput>,
    config: Arc<RunConfig>,
    data_type: StreamDataType,
    label: &str,                    // still &str here
) -> Result<JoinHandle<Result<()>>> {

    let dest_path_clone = dest_path.clone();
    let label = label.to_string();  // ← clone to owned String

    let write_handle = tokio::spawn(async move {
        let file = TokioFile::create(&dest_path_clone)
            .await
            .map_err(|e| anyhow!("Cannot create file {}: {}", dest_path_clone.display(), e))?;

        let effective_buffer = crate::utils::system::compute_buffer_size(
            &config,
            "write_byte_stream_to_file",
            data_type,
            1.8,
        );

        let mut writer = tokio::io::BufWriter::with_capacity(effective_buffer, file);

        let mut stream = stream;
        let mut batch: Vec<u8> = Vec::with_capacity(effective_buffer / 4);

        while let Some(item) = stream.next().await {
            let bytes = item
                .to_bytes()
                .map_err(|e| anyhow!("Failed to convert to bytes: {}", e))?;

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
              label, dest_path_clone.display(), final_size, effective_buffer / (1024*1024));

        Ok(())
    });

    Ok(write_handle)
}

pub async fn file_size(path: &PathBuf) -> Result<u64> {
    let metadata = tokio::fs::metadata(path).await
        .map_err(|e| anyhow!("Failed to read file metadata {}: {}", path.display(), e))?;
    Ok(metadata.len())
}

pub async fn available_space_for_path(path: &PathBuf) -> Result<u64> {
    let target = path.canonicalize().unwrap_or_else(|_| path.clone());
    info!("Querying space for canonical path: {}", target.display());

    // Primary: nix statvfs (reliable for tmpfs/NVMe)
    match statvfs::statvfs(&target) {
        Ok(stat) => {
            let block_size = stat.block_size() as u64;               // safe cast
            let blocks_avail = stat.blocks_available() as u64;       // safe cast
            let avail = block_size * blocks_avail;
            info!("statvfs success: block_size={} blocks_avail={} → {} bytes available", block_size, blocks_avail, avail);
            return Ok(avail);
        }
        Err(e) => warn!("statvfs failed: {}", e),
    }

    // Fallback: sysinfo
    let disks = Disks::new_with_refreshed_list();
    for disk in disks.list() {
        let mp = disk.mount_point();
        if target.starts_with(mp) {
            let avail = disk.available_space();
            info!("sysinfo match on {}: {} bytes", mp.display(), avail);
            return Ok(avail);
        }
    }

    Err(anyhow!("No matching filesystem found for {}", target.display()))
}


/// pick RAM-backed temp dir or fallback to disk temp dir
/// based on avialble space qand headroom
///
/// # Arguments
/// * `estimated_bytes` –  file size, buffer size in bytes
/// * `ram_dir` – e.g., `/dev/shm` or `config.ram_temp_dir`
/// * `headroom_factor` – how much of available RAM you’re willing to use (e.g., 4 = 25%)
///
/// # Returns
/// The chosen temp directory (`ram_dir` if safe, otherwise `std::env::temp_dir()`)
pub async fn choose_temp_dir(
    estimated_bytes: u64,
    ram_dir: &PathBuf,
    nvme_scratch: &Option<String>,
    headroom_factor: u64,
    prefer_nvme: bool,
) -> Result<TempDir, PipelineError> {
    async fn check_space(path: &PathBuf, required: u64, factor: u64) -> Result<bool, PipelineError> {
        let avail = available_space_for_path(path).await?;
        info!("choose temp dir avilable bytes {}", avail);
        Ok(required <= avail / factor)
    }

    async fn try_create_temp(
        dir_path: &PathBuf,
        estimated: u64,
        factor: u64,
        label: &str,
    ) -> Option<Result<TempDir, PipelineError>> {
        if !check_space(dir_path, estimated, factor).await.ok()? {
            return None;
        }

        debug!(
            "Trying {} dir {}: {} bytes fits in {} available (/{})",
            label,
            dir_path.display(),
            estimated,
            available_space_for_path(dir_path).await.ok()?,
            factor
        );

        let temp_dir = match TempDir::new_in(dir_path) {
            Ok(td) => td,
            Err(e) => return Some(Err(PipelineError::Other(anyhow!(
                "Failed to create TempDir in {} dir {}: {}", label, dir_path.display(), e
            )))),
        };

        // ────────────────────────────────────────────────────────────────
        //  force directory visibility before returning
        // This closes the race window where DIAMOND sees ENOENT
        // ────────────────────────────────────────────────────────────────
        if let Err(e) = tokio::fs::create_dir(temp_dir.path().join("probe")).await {
            warn!("Filesystem probe failed in {}: {}. Continuing anyway.", dir_path.display(), e);
            // Not fatal — we still return the dir, just log
        } else {
            let _ = tokio::fs::remove_dir(temp_dir.path().join("probe")).await;
        }

        // Optional ultra-paranoid version (usually overkill but zero risk):
        // tokio::fs::sync_all().await.ok();   // fsync whole filesystem — slow!

        debug!("Created and probed temp dir: {}", temp_dir.path().display());
        Some(Ok(temp_dir))
    }


    let ram_attempt = try_create_temp(ram_dir, estimated_bytes, headroom_factor, "RAM").await;
    let nvme_path = nvme_scratch.as_ref().map(PathBuf::from);
    let nvme_attempt = if let Some(ref nvme) = nvme_path {
        try_create_temp(nvme, estimated_bytes, headroom_factor, "NVMe").await
    } else {
        debug!("No NVMe scratch configured — skipping");
        None
    };

    let chosen = if prefer_nvme {
        info!("choosing nvme");
        nvme_attempt.or(ram_attempt)
    } else {
        info!("choosing ram");
        ram_attempt.or(nvme_attempt)
    };

    if let Some(res) = chosen {
        return res;
    }

    debug!("Falling back to system temp dir");
    let fallback = TempDir::new()
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create system TempDir: {}", e)))?;

    // Probe fallback too
    let _ = tokio::fs::create_dir(fallback.path().join("probe")).await;
    let _ = tokio::fs::remove_dir(fallback.path().join("probe")).await;

    Ok(fallback)
}

// Helper function for writing byte streams to files and force and EOF
// Meant for small files writing to /dev/shm or nvme scratch
//
pub async fn materialize_stream_to_file(
    config: Arc<RunConfig>,
    stream: ReceiverStream<ParseOutput>,
    output_path: PathBuf,
    label: &str,
    data_type: StreamDataType,
) -> Result<PathBuf> {
    let write_handle = write_byte_stream_to_file(
        &output_path,
        stream,
        config.clone(),
        data_type,
        label,
    ).await?;

    write_handle.await??; // Wait for completion + propagate errors

    debug!("Materialized {} → {} (size: {} bytes)",
          label, output_path.display(), file_size(&output_path).await?);

    Ok(output_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[tokio::test]
    async fn test_write_vecu8_to_file() -> std::io::Result<()> {
        let test_data: Arc<Vec<u8>> = Arc::new(vec![1, 2, 3, 4]);
        let temp_name = NamedTempFile::new()?;
        let temp_path = temp_name.into_temp_path();
        let write_task = write_vecu8_to_file(test_data.clone(), &temp_path, 10000)
            .await
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        write_task
            .await
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        let mut read_file = std::fs::File::open(&temp_path)?;
        let mut chk_buffer = Vec::new();
        read_file.read_to_end(&mut chk_buffer)?;
        assert_eq!(chk_buffer, *test_data);
        std::fs::remove_file(&temp_path)?;
        Ok(())
    }

    #[tokio::test]
    async fn test_write_vecu8_to_file_empty() -> std::io::Result<()> {
        let test_data: Arc<Vec<u8>> = Arc::new(vec![]);
        let temp_name = NamedTempFile::new()?;
        let temp_path = temp_name.into_temp_path();
        let write_task = write_vecu8_to_file(test_data.clone(), &temp_path, 10000)
            .await
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        let result = write_task
            .await
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        assert!(result.is_err(), "Expected error for empty data");
        assert_eq!(result.unwrap_err().to_string(), format!("No data written to file at {}", temp_path.display()));
        Ok(())
    }


}