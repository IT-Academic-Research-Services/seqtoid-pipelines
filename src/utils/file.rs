use std::fs::File;
use std::io;
use std::sync::Arc;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use crate::config::defs::GZIP_EXT;
use anyhow::{Result, anyhow};
use tempfile::NamedTempFile;
use tokio::io::{AsyncWriteExt, BufWriter};
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::streams::ParseOutput;
use tokio::task::JoinHandle;
use tokio::fs::File as TokioFile;
use tokio::process::Command;
use tokio_stream::StreamExt;
use tempfile::Builder as TempfileBuilder;
use crate::utils::streams::ToBytes;
use crate::utils::fastx::{SequenceRecord, r1r2_base};
use crate::config::defs::{RunConfig, PipelineError};




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
pub fn extension_remover(path: &PathBuf) -> (PathBuf, Vec<String>) {
    let path = Path::new(path);
    let parent = path.parent().unwrap_or(Path::new(""));
    let mut current = path;
    let mut extensions = Vec::new();

    while let Some(ext) = current.extension().and_then(|e| e.to_str()) {

        match extensions.len() {
            1 => {
                if extensions[0] == "gz" {
                    extensions.push(ext.to_string());
                } else {
                    break;
                }
            }
            0 => {
                extensions.push(ext.to_string());
            }
            _ => {
                break;
            }
        }

        current = current.file_stem().map(Path::new).unwrap_or(Path::new(""));
    }

    let stem = current.to_string_lossy().into_owned();
    let out_path_buf = if parent == Path::new("") {
        PathBuf::from(stem)
    } else {
        parent.join(stem)
    };
    extensions.reverse();
    (out_path_buf, extensions)
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
        let mut writer_file = TokioFile::create(&fifo_path)
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
pub async fn write_parse_output_to_temp_file<P: AsRef<Path>>(
    input_stream: ReceiverStream<ParseOutput>,
    buffer_size: Option<usize>,
    suffix: Option<&str>,
    ram_temp_dir: P,
) -> Result<(JoinHandle<Result<(), anyhow::Error>>, PathBuf, tempfile::NamedTempFile)> {
    let mut temp_name = match suffix {
        Some(s) => TempfileBuilder::new()
            .suffix(s)
            .tempfile_in(ram_temp_dir.as_ref()),
        None => TempfileBuilder::new().tempfile_in(ram_temp_dir.as_ref()),
    }
        .map_err(|e| anyhow!("Failed to create temp file in {}: {}", ram_temp_dir.as_ref().display(), e))?;
    let temp_path = temp_name.path().to_path_buf();

    let writer_file = TokioFile::create(&temp_path).await?;
    let buffer_capacity = buffer_size.unwrap_or(16 * 1024 * 1024);
    let temp_path_clone = temp_path.clone();

    let task = tokio::spawn(async move {
        let mut writer = BufWriter::with_capacity(buffer_capacity, writer_file);
        let mut input_stream = input_stream;
        let mut byte_count = 0;

        while let Some(item) = input_stream.next().await {
            match item {
                ParseOutput::Bytes(data) => {
                    writer.write_all(&data).await?;
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

        writer.flush().await?;
        writer.shutdown().await?;
        if byte_count == 0 {
            return Err(anyhow!("No data written to temp file at {}", temp_path_clone.display()));
        }
        Ok(())
    });

    Ok((task, temp_path, temp_name))
}

pub fn validate_file_inputs(
    config: &RunConfig,
    cwd: &PathBuf,
) -> Result<(PathBuf, Option<PathBuf>, PathBuf, String), PipelineError> {
    let file1_path: PathBuf = match &config.args.file1 {
        Some(file) => {
            let file1_full_path = file_path_manipulator(&PathBuf::from(file), Some(cwd), None, None, "");
            if file1_full_path.exists() {
                file1_full_path
            } else {
                return Err(PipelineError::FileNotFound(file1_full_path));
            }
        }
        None => return Err(PipelineError::InvalidConfig("File1 path required".to_string())),
    };

    let sample_base: String;
    let file1_r1r2 = r1r2_base(&file1_path);
    sample_base = match file1_r1r2.prefix {
        Some(prefix) => prefix,
        None => {
            eprintln!("No R1 tag found. Using bare file 1 stem as sample_base.");
            file1_path.to_string_lossy().into_owned()
        }
    };

    if !file1_path.exists() {
        return Err(PipelineError::FileNotFound(file1_path));
    }

    let sample_base_buf: PathBuf = PathBuf::from(&sample_base);
    let (no_ext_sample_base_buf, _) = extension_remover(&sample_base_buf);
    let no_ext_sample_base = no_ext_sample_base_buf.to_string_lossy().into_owned();

    let file2_path: Option<PathBuf> = match &config.args.file2 {
        Some(file) => {
            let file2_full_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd.clone()), None, None, "");
            if file2_full_path.exists() {
                Some(file2_full_path)
            } else {
                eprintln!("File2 path does not exist: {}", file2_full_path.display());
                None
            }
        }
        None => None,
    };

    Ok((file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;
    use tokio::sync::mpsc;

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

    #[tokio::test]
    async fn test_write_parse_output_to_temp_file_fastq() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let fastq = SequenceRecord::Fastq {
            id: "read1".to_string(),
            desc: None,
            seq: Arc::new(b"ATCG".to_vec()),
            qual: Arc::new(b"IIII".to_vec()),
        };
        tokio::spawn(async move { tx.send(ParseOutput::Fastq(fastq)).await.unwrap(); });
        let (task, path, _temp) = write_parse_output_to_temp_file(ReceiverStream::new(rx), None, Some(".fq"), std::env::temp_dir()).await?;
        task.await??;
        let file = std::fs::File::open(&path)?;
        let mut contents = String::new();
        std::io::BufReader::new(file).read_to_string(&mut contents)?;
        assert_eq!(contents, "@read1\nATCG\n+\nIIII\n");
        Ok(())
    }
}