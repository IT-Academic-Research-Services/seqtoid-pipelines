// Functions and definitions for the minimap2-associated PAF file format
use anyhow::{anyhow, Result};
use log::{info, debug};
use std::cmp::Reverse;
use std::collections::HashMap;
use std::sync::Arc;

use tokio::sync::mpsc::Sender;
use tokio_stream::{StreamExt, wrappers::ReceiverStream, StreamMap};  

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

    pub async fn merge_paf_streams(
        streams: Vec<ReceiverStream<ParseOutput>>,
        tx: Sender<ParseOutput>,
        _buffer_size: usize,
    ) -> Result<()> {
        info!("Starting merge_paf_streams with {} streams", streams.len());
        let mut query_hits: HashMap<String, Vec<PafRecord>> = HashMap::new();
        let mut stream_map = StreamMap::new();

        for (idx, stream) in streams.into_iter().enumerate() {
            stream_map.insert(idx, stream);
        }

        let mut record_count = 0;
        let start_time = tokio::time::Instant::now();

        // Collect all PAF records from all partial streams concurrently
        while let Some((_idx, item)) = stream_map.next().await {
            match item {
                ParseOutput::Bytes(bytes) => {
                    let line = String::from_utf8_lossy(&*bytes).trim().to_string();
                    for line in line.lines() {
                        let trimmed = line.trim();
                        if trimmed.is_empty() {
                            continue;
                        }
                        let record = PafRecord::parse_line(trimmed)?;
                        let unique_queries = query_hits.len();
                        let hits = query_hits
                            .entry(record.qname.clone())
                            .or_default();
                        
                        hits.push(record);
                        record_count += 1;
                        if record_count % 100_000 == 0 {
                            info!("merge_paf_streams: processed {} records, {} unique queries, elapsed {:?}", 
                                record_count, unique_queries, start_time.elapsed());
                        }

                        if hits.len() > 20 {
                            hits.sort_by_key(|h| Reverse(h.alignment_score()));
                            hits.truncate(6);
                        }
                    }
                }
                _ => {
                    return Err(anyhow!("Unexpected non-Bytes item in PAF stream"));
                }
            }
        }

        info!("merge_paf_streams: finished collecting {} records ({} unique queries) in {:?}", 
            record_count, query_hits.len(), start_time.elapsed());

        // Sort queries alphabetically for deterministic output
        let mut queries: Vec<String> = query_hits.keys().cloned().collect();
        queries.sort();

        info!("merge_paf_streams: emitting merged records...");
        let emit_start = tokio::time::Instant::now();
        let mut emitted_count = 0;

        // For each query: sort by descending AS score, keep top 6 hits
        for query in queries {
            if let Some(mut hits) = query_hits.remove(&query) {
                hits.sort_by_key(|h| Reverse(h.alignment_score()));
                hits.truncate(6);  // 1 primary + up to 5 secondaries

                for hit in hits {
                    let paf_line = hit.to_paf_line();
                    tx.send(ParseOutput::Bytes(Arc::new(paf_line.into_bytes())))
                        .await
                        .map_err(|e| anyhow!("Failed to send merged PAF line: {}", e))?;
                    emitted_count += 1;
                }
            }
        }

        info!("merge_paf_streams: finished emitting {} records in {:?}", emitted_count, emit_start.elapsed());

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