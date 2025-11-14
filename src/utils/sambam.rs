use std::collections::{HashMap, HashSet};

use futures::StreamExt;
use anyhow::{Result, anyhow};
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;

use crate::utils::streams::ParseOutput;


pub async fn stream_sam_alignment_counter(
    rx: mpsc::Receiver<ParseOutput>,
    early_exit: bool,
) -> Result<u64, anyhow::Error> {
    let mut stream = ReceiverStream::new(rx);
    let mut counter: u64 = 0;

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(chunk) => {
                let chunk_str = String::from_utf8_lossy(&chunk);
                for line in chunk_str.lines() {
                    if !line.starts_with('@') {
                        counter += 1;
                        if early_exit {
                            return Ok(1); // Early exit: Not empty
                        }
                    }
                }
            }
            _ => return Err(anyhow!("Unexpected ParseOutput in SAM stream; expected Bytes")),
        }
    }

    Ok(counter)
}


pub async fn generate_info_from_sam_stream(
    rx: mpsc::Receiver<ParseOutput>,
    duplicate_cluster_sizes: &HashMap<String, u64>,
    min_contig_size: usize,
) -> Result<HashMap<String, u64>> {
    let mut stream = ReceiverStream::new(rx);
    let mut contig_stats = HashMap::new();
    let mut contig_unique_counts = HashMap::new();
    let mut seen = HashSet::new();

    while let Some(item) = stream.next().await {
        let bytes = match item {
            ParseOutput::Bytes(b) => b,
            _ => return Err(anyhow!("Expected Bytes in SAM stream")),
        };

        let line = std::str::from_utf8(&bytes)
            .map_err(|e| anyhow!("SAM line not UTF-8: {}", e))?;

        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let read = fields[0];
        if seen.contains(read) {
            continue;
        }
        seen.insert(read.to_string());

        let contig = fields[2];
        if contig == "*" {
            continue;
        }

        let cluster_size = duplicate_cluster_sizes.get(read).copied().unwrap_or(1);
        *contig_stats.entry(contig.to_string()).or_insert(0) += cluster_size;
        *contig_unique_counts.entry(contig.to_string()).or_insert(0) += 1;
    }

    // Filter by min_contig_size
    contig_stats.retain(|contig, _| {
        contig_unique_counts.get(contig).copied().unwrap_or(0) >= min_contig_size
    });

    Ok(contig_stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tokio::sync::mpsc;
    use tokio_stream::wrappers::ReceiverStream;

    #[tokio::test]
    async fn test_stream_sam_alignment_counter() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let sam_data = vec![
            ParseOutput::Bytes(b"@HD\tVN:1.6\n@SQ\tSN:ERCC-00002\tLN:1061\nread1\t77\t*\t0\t0\t*\t*\t0\t0\tATCG\tIIII\n".to_vec().into()),
        ];

        tokio::spawn(async move {
            for item in sam_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_sam_alignment_counter(rx, true).await?;
        assert_eq!(count, 1, "Should count 1 SAM alignment with early_exit");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_sam_alignment_counter_empty() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let sam_data = vec![ParseOutput::Bytes(b"@HD\tVN:1.6\n@SQ\tSN:ERCC-00002\tLN:1061\n".to_vec().into())];

        tokio::spawn(async move {
            for item in sam_data {
                tx.send(item).await.unwrap();
            }
        });

        let count = stream_sam_alignment_counter(rx, false).await?;
        assert_eq!(count, 0, "Should count 0 for header-only SAM");
        Ok(())
    }
}
