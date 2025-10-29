// Functions and definitions for the minimap2-associated PAF file format
use anyhow::{Result, anyhow};

const LAMBDA: f64 = 1.58;
const K: f64 = 0.1;

#[derive(Debug, Clone)]
pub struct PafRecord {
    qname: String,
    qlen: u64,
    qstart: u64,
    qend: u64,
    strand: char,
    tname: String,
    tlen: u64,
    tstart: u64,
    tend: u64,
    nmatch: u64,
    alen: u64,
    mapq: u64,
    tags: std::collections::HashMap<String, String>,
}

impl PafRecord {
    pub fn parse_line(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        let qname = fields.next().ok_or_else(|| anyhow!("Missing qname"))?.to_string();
        let qlen = fields.next().ok_or_else(|| anyhow!("Missing qlen"))?.parse()?;
        let qstart = fields.next().ok_or_else(|| anyhow!("Missing qstart"))?.parse()?;
        let qend = fields.next().ok_or_else(|| anyhow!("Missing qend"))?.parse()?;
        let strand = fields.next().ok_or_else(|| anyhow!("Missing strand"))?.chars().next().ok_or_else(|| anyhow!("Invalid strand"))?;
        let tname = fields.next().ok_or_else(|| anyhow!("Missing tname"))?.to_string();
        let tlen = fields.next().ok_or_else(|| anyhow!("Missing tlen"))?.parse()?;
        let tstart = fields.next().ok_or_else(|| anyhow!("Missing tstart"))?.parse()?;
        let tend = fields.next().ok_or_else(|| anyhow!("Missing tend"))?.parse()?;
        let nmatch = fields.next().ok_or_else(|| anyhow!("Missing nmatch"))?.parse()?;
        let alen = fields.next().ok_or_else(|| anyhow!("Missing alen"))?.parse()?;
        let mapq = fields.next().ok_or_else(|| anyhow!("Missing mapq"))?.parse()?;

        let mut tags = std::collections::HashMap::new();
        for tag_str in fields {
            let parts: Vec<&str> = tag_str.splitn(3, ':').collect();
            if parts.len() == 3 {
                let key = parts[0].to_string();
                let _typ = parts[1]; // We can check if needed
                let val = parts[2].to_string();
                tags.insert(key, val);
            }
        }

        Ok(Self {
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, tags,
        })
    }

    fn calc_bitscore(&self) -> f64 {
        let nonmatch = self.tags.get("NM").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let alen = self.alen as f64;
        let score = alen - 2.0 * nonmatch;
        (score * LAMBDA - K.ln()) / 2.0f64.ln()
    }

    fn calc_evalue(&self, genome_size: f64) -> f64 {
        let nonmatch = self.tags.get("NM").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let alen = self.alen as f64;
        let score = (alen - 2.0 * nonmatch).max(0.0);
        K * alen * genome_size * (-LAMBDA * score).exp()
    }

    fn calc_gap_openings(&self) -> u64 {
        let cigar = self.tags.get("cg").cloned().unwrap_or_default();
        let mut go = 0;
        for c in cigar.chars() {
            if c == 'I' || c == 'D' {
                go += 1;
            }
        }
        go
    }

    fn percent_identity(&self) -> f64 {
        if self.alen == 0 {
            0.0
        } else {
            (self.nmatch as f64 / self.alen as f64) * 100.0
        }
    }

    pub fn to_m8_line(&self, genome_size: f64) -> String {
        let mut tstart = self.tstart;
        let mut tend = self.tend;
        if self.strand == '-' {
            std::mem::swap(&mut tstart, &mut tend);
        }
        let qstart_1 = self.qstart + 1;
        let mut tstart_adj = tstart;
        let mut tend_adj = tend;
        if self.strand == '+' {
            tstart_adj += 1;
        } else {
            tend_adj += 1;
        }
        let nonmatch = self.tags.get("NM").and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
        let gap_openings = self.calc_gap_openings();
        let percent_ident = self.percent_identity();
        let bitscore = self.calc_bitscore();
        let evalue = self.calc_evalue(genome_size);

        format!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.3}\n",
            self.qname, self.tname, percent_ident, self.alen,
            nonmatch, gap_openings, qstart_1, self.qend,
            tstart_adj, tend_adj, evalue, bitscore
        )
    }
}