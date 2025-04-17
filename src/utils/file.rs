use std::fs;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use crate::GZIP_EXT;

pub fn is_gzipped(path: &PathBuf) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut buffer = [0u8; 2];
    file.read_exact(&mut buffer)?;
    Ok(buffer == [0x1F, 0x8B]) // Gzip magic bytes
}

pub fn write_fastq_record(
    file: &mut File,
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    if let Some(desc) = desc {
        writeln!(file, "@{} {}", id, desc)?;
    } else {
        writeln!(file, "@{}", id)?;
    }
    file.write_all(seq)?;
    writeln!(file)?;
    writeln!(file, "+")?;
    file.write_all(qual)?;
    writeln!(file)?;
    Ok(())
}

pub fn write_fasta_record(
    file: &mut File,
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
) -> io::Result<()> {
    if let Some(desc) = desc {
        writeln!(file, ">{} {}", id, desc)?;
    } else {
        writeln!(file, ">{}", id)?;
    }
    file.write_all(seq)?;
    writeln!(file)?;
    Ok(())
}


/// Calls extension_remover to retrieve base and extensions.
/// Then prepends prefix and appends postfix, if any.
/// Then returns absolute path.
/// # Arguments
///
/// * `path`: PAthBuf - File path
/// * 'prefix': Option<&str> = added ahead of base
/// * 'postfix': Option<&str> = added after of base
///
/// # Returns
/// Option<PathBuf>: modified, absolute path
///
pub fn file_path_manipulator(path: PathBuf, prefix: Option<&str>, postfix: Option<&str>) -> Option<PathBuf> {
    
    let absolute_path = match fs::canonicalize(path) {
        Ok(path) => path,
        Err(_) => return None,
    };

    
    let (stem, extensions) = extension_remover(&absolute_path);
    
    let base = stem.to_str().unwrap_or("");
    let new_base = format!(
        "{}{}{}",
        prefix.unwrap_or(""),
        base,
        postfix.unwrap_or("")
    );

    let new_file_name = if extensions.is_empty() {
        new_base
    } else {
        format!("{}.{}", new_base, extensions.join("."))
    };
    
    let parent = absolute_path
        .parent()
        .unwrap_or_else(|| Path::new(""));
    let new_path = parent.join(new_file_name);

    match fs::canonicalize(&new_path) {
        Ok(canonical_path) => Some(canonical_path),
        Err(_) => None,
    }
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