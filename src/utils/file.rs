use std::fs::File;
use std::io;
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
/// * `path`: PAthBuf - File path
/// * 'prefix': Option<&str> = added ahead of base
/// * 'postfix': Option<&str> = added after of base
/// /// * 'postfix': Option<&str> = added after of base
// /// * 'delimiter': &str = added between prefix, postfix and base.
///
/// # Returns
/// PathBuf: modified, absolute path
pub fn file_path_manipulator(path: &PathBuf, parent_dir: &PathBuf, prefix: Option<&str>, postfix: Option<&str>, delimiter: &str) -> PathBuf {
    
    let absolute_path = parent_dir.canonicalize().ok().expect("Parent directory not found. {parent_dir:?}");
    let new_file_name = file_name_manipulator(path, prefix, postfix, delimiter);
    let new_file_path = PathBuf::from(new_file_name);
    let absolute_file_path = absolute_path.join(&new_file_path);

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

pub fn extension_remover(path: &PathBuf) -> (PathBuf, Vec<String>) {
    let path = Path::new(&path);
    let mut current = path;
    let mut extensions = Vec::new();

    while let Some(ext) = current.extension().and_then(|e| e.to_str()) {

        match extensions.len() {
            1 => {
                if extensions[0] == GZIP_EXT {
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
    
    let out_path_buf = PathBuf::from(current.to_string_lossy().into_owned());
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
            let full_path = file_path_manipulator(&path, dir, None, None, "");
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
        std::fs::remove_file(fifo_path)?;
    }
    
    let status = Command::new("mkfifo")
        .arg(fifo_path)
        .status()
        .await
        .map_err(|e| anyhow!("Failed to execute mkfifo for {}: {}", fifo_path.display(), e))?;
    if !status.success() {
        return Err(anyhow!("mkfifo failed with status: {}", status));
    }
    
    let buffer_capacity = buffer_size.unwrap_or(4 * 1024 * 1024); // 4MB
    
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
                _ => {
                    return Err(anyhow!("Expected Bytes in stream, got unexpected data at {}", fifo_path.display()));
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

        eprintln!("Wrote {} bytes to FIFO at {}", byte_count, fifo_path.display());
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
///
/// # Returns
///
/// A `Result` containing a `JoinHandle` that resolves to `Result<(), anyhow::Error>` upon completion.
pub async fn write_parse_output_to_temp(
    input_stream: ReceiverStream<ParseOutput>,
    buffer_size: Option<usize>,
) -> Result<(JoinHandle<Result<(), anyhow::Error>>, PathBuf)> {

    let temp_name = NamedTempFile::new()?;
    let temp_path = temp_name.path().to_path_buf();
    temp_name.close()?;     // Immediately close the regular file to avoid conflicts with mkfifo
    
    let task = write_parse_output_to_fifo(&temp_path, input_stream, buffer_size).await?;

    Ok((task, temp_path))
}


