use std::fs;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

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

    let file_name = match absolute_path.file_name().and_then(|name| name.to_str()) {
        Some(name) => name,
        None => return None, 
    };

    let path = Path::new(file_name);
    let mut extensions = Vec::new();
    let mut current = path;

    while let Some(ext) = current.extension().and_then(|e| e.to_str()) {
        extensions.push(ext);
        current = Path::new(current.file_stem().unwrap_or_default());
    }
    extensions.reverse();
    
    let base = current.to_str().unwrap_or("");
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

