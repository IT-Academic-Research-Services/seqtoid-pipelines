use std::fs::File;
use std::io;
use std::io::Read;

pub fn is_gzipped(path: &str) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut buffer = [0u8; 2];
    file.read_exact(&mut buffer)?;
    Ok(buffer == [0x1F, 0x8B])
}