use std::fs::File;
use std::io;
use std::io::{BufReader, Read, BufWriter, Write};
use std::path::{Path, PathBuf};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use crate::GZIP_EXT;


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

pub enum FileWriter {
    Uncompressed(BufWriter<File>),
    Gzipped(GzEncoder<BufWriter<File>>),
}

/// Trait to abstract writing items to a file
pub trait WriteToFile {
    fn write_to_file<W: Write>(&self, writer: &mut W) -> io::Result<()>;
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

/// Implementation for generic byte-like types
impl<T> WriteToFile for T
where
    T: AsRef<[u8]>,
{
    fn write_to_file<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(self.as_ref())?;
        writer.write_all(b"\n")?;
        Ok(())
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
    
    let mut out_path_buf = PathBuf::from(current.to_string_lossy().into_owned());
    extensions.reverse();
    (out_path_buf, extensions)
}