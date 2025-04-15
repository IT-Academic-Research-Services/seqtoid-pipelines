use std::fs;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::PathBuf;
use crate::GZIP_EXT;

pub fn is_gzipped(path: &str) -> io::Result<bool> {
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
/// * `path`: &str - File path
/// * 'prefix': Option<&str> = added ahead of base
/// * 'postfix': Option<&str> = added after of base
///
/// # Returns
/// Option<PathBuf>: modified, absolute path
///
pub fn file_name_manipulator(path: &str, prefix: Option<&str>, postfix: Option<&str>) -> Option<PathBuf> {
    
    let absolute_path = match fs::canonicalize(path) {
        Ok(path) => path,
        Err(_) => return None, // Return None if path doesn't exist or canonicalization fails
    };

    let base_file_name = match absolute_path.file_name().and_then(|name| name.to_str()) {
        Some(name) => name,
        None => return None, // Return None if no file name or non-UTF-8
    };

    let (base_parts, extensions) = extension_remover(base_file_name);
    let base = base_parts.into_iter().collect::<Vec<_>>().join(".");
    
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
        .unwrap_or_else(|| std::path::Path::new(""));
    let new_path = parent.join(new_file_name);

    match fs::canonicalize(&new_path) {
        Ok(canonical_path) => Some(canonical_path),
        Err(_) => Some(new_path), 
    }
}


/// Separates a file name, removing its extension and returning as a vec.
/// If the last ext is .gz, will include the part before it (e.g. fastq.gz)
/// # Arguments
///
/// * `path`: &str - File path
///
/// # Returns
/// (Vec<String>, Vec<String>): base, extension(s)
///
fn extension_remover(path: &str) -> (Vec<String>, Vec<String>) {


    let path_parts: Vec<&str> = path.split('.').collect();
    let mut extensions: Vec<String> = Vec::new();
    let mut base: Vec<String> = Vec::new();
    let mut count = 0;
    let mut gzip_present = false;

    for part in path_parts.iter().rev() {

        if count == 0 {
            if *part == GZIP_EXT {
                extensions.push(part.to_string());
                gzip_present = true;
            } else {
                extensions.push(part.to_string());
            }
        } else if count == 1  {
            if gzip_present {
                extensions.push(part.to_string());
            }
            else {
                base.push(part.to_string());
            }

        }
        else {
            base.push(part.to_string());
        }
        count += 1; // Increment count to process next part
    }
    
    let revbase : Vec<String> = base.iter().rev().cloned().collect();
    let revext : Vec<String> = extensions.iter().rev().cloned().collect();
    (revbase, revext)
}
