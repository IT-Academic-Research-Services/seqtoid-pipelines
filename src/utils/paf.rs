// Functions and definitions for the minimap2-associated PAF file format
use anyhow::{anyhow, Result};
use log::{info, debug, warn};
use std::cmp::Reverse;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::sync::mpsc::{self, Sender};
use tokio_stream::{StreamExt, wrappers::ReceiverStream, StreamMap};  
use dashmap::DashMap;
use rayon::prelude::*;

use crate::utils::streams::ParseOutput;

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
        let line = line.trim_end();
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

    pub fn alignment_score(&self) -> i32 {
        self.tags
            .get("AS")
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0)
    }

    pub fn to_paf_line(&self) -> String {
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname, self.qlen, self.qstart, self.qend, self.strand,
            self.tname, self.tlen, self.tstart, self.tend, self.nmatch, self.alen, self.mapq,
        );

        for (key, value) in &self.tags {
            line.push_str(&format!("\t{}:Z:{}", key, value));
        }

        line
    }

    pub async fn merge_paf_files(
        paf_paths: &[PathBuf],
        tx: Sender<ParseOutput>,
    ) -> Result<()> {
        info!("Merging {} raw PAF files → streaming PAF output (no conversion)", paf_paths.len());

        let start = std::time::Instant::now();
        let mut total_lines = 0u64;
        let mut skipped = 0u64;

        let mut handles = Vec::new();

        for path in paf_paths {
            let path = path.clone();
            let tx_clone = tx.clone();

            handles.push(tokio::spawn(async move {
                let file = File::open(&path).await
                    .map_err(|e| anyhow!("Failed to open PAF {}: {}", path.display(), e))?;

                let mut reader = BufReader::new(file);
                let mut line = String::new();
                let mut local_lines = 0u64;
                let mut local_skipped = 0u64;

                while reader.read_line(&mut line).await? > 0 {
                    local_lines += 1;
                    let trimmed = line.trim();

                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        local_skipped += 1;
                        line.clear();
                        continue;
                    }

                    // Forward raw PAF line unchanged
                    let bytes = (trimmed.to_string() + "\n").into_bytes();
                    if tx_clone.send(ParseOutput::Bytes(Arc::new(bytes))).await.is_err() {
                        return Ok(());
                    }

                    line.clear();
                }

                info!("Merged {}: {} lines ({} skipped)", path.display(), local_lines, local_skipped);
                Ok::<(), anyhow::Error>(())
            }));
        }

        for h in handles {
            h.await??;
        }

        info!("Raw PAF merge complete: {} total lines in {:.2?}", total_lines, start.elapsed());
        Ok(())
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
        let score = alen - 2.0 * nonmatch;
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

    fn extract_accession(&self) -> String {
        self.tname
            .split(':')
            .find(|part| part.starts_with("NT:") || part.starts_with("NR:"))
            .and_then(|s| s.split('|').next())
            .map(|s| s.trim().to_string())
            .or_else(|| {
                self.tname
                    .split(|c| c == ' ' || c == '\t')
                    .next()
                    .map(|s| s.strip_prefix('>').unwrap_or(s).trim().to_string())
            })
            .unwrap_or_default()
    }

    pub fn to_m8_line(&self, genome_size: f64) -> String {
        let accession = self.extract_accession();
        let base_accession = accession.split('.').next().unwrap_or(&accession).to_string();

        if accession.is_empty() {
            return String::new();
        }

        let (mut tstart, mut tend) = (self.tstart, self.tend);
        if self.strand == '-' {
            std::mem::swap(&mut tstart, &mut tend);
        }
        let qstart_1 = self.qstart + 1;
        let tstart_adj = if self.strand == '+' { tstart + 1 } else { tstart };
        let tend_adj   = if self.strand == '+' { tend } else { tend + 1 };

        let nonmatch = self.tags.get("NM")
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0);
        let gap_openings = self.calc_gap_openings();
        let percent_ident = self.percent_identity();
        let bitscore = self.calc_bitscore();
        let evalue = self.calc_evalue(genome_size);

        format!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.3}\n",
            self.qname, base_accession, percent_ident, self.alen,
            nonmatch, gap_openings, qstart_1, self.qend,
            tstart_adj, tend_adj, evalue, bitscore
        )
    }
}

pub fn parse_paf_batch_to_m8(batch: Vec<u8>, genome_size: f64) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && !line.starts_with(b"#"))
        .flat_map(|line_bytes: &[u8]| {
            if let Ok(line) = std::str::from_utf8(line_bytes) {
                if let Ok(record) = PafRecord::parse_line(line) {
                    let m8_line = record.to_m8_line(genome_size);
                    return vec![(m8_line + "\n").into_bytes()];
                }
            }
            vec![]
        })
        .collect()
}