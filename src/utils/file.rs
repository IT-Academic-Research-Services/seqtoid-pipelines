use std::fs::File;
use std::io;
use std::io::{Read, Write};

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