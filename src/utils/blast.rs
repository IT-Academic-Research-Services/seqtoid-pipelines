// BLAST-related file functions and structures
use anyhow::{anyhow, Result};


/// Single BLAST m8 line
#[derive(Debug, Clone)]
struct M8Record {
    qname: String,
    tname: String,
    pident: f64,
    alen: u64,
    mismatch: u64,
    gapopen: u64,
    qstart: u64,
    qend: u64,
    tstart: u64,
    tend: u64,
    evalue: f64,
    bitscore: f64,
}

impl M8Record {
    fn parse_line(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        macro_rules! next {
            () => {
                fields.next().ok_or_else(|| anyhow!("missing field"))?
            };
        }
        Ok(Self {
            qname: next!().to_string(),
            tname: next!().to_string(),
            pident: next!().parse()?,
            alen: next!().parse()?,
            mismatch: next!().parse()?,
            gapopen: next!().parse()?,
            qstart: next!().parse()?,
            qend: next!().parse()?,
            tstart: next!().parse()?,
            tend: next!().parse()?,
            evalue: next!().parse()?,
            bitscore: next!().parse()?,
        })
    }
}