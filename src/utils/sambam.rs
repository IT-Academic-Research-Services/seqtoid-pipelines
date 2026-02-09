use anyhow::{anyhow, Result};
use bytes::Bytes;
use futures::StreamExt;
use std::collections::{HashMap, HashSet};
use std::io;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tokio_util::io::StreamReader;
use noodles::bam::r#async::io::Reader as BamAsyncReader;
use noodles::bam::record::Record;
use noodles::bam::record::Record as BamRecord;
use noodles::sam::Header;

use crate::utils::streams::ParseOutput;
use crate::config::defs::{ClusterInfo, DuplicateClusters};

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


fn anyhow_to_io(e: anyhow::Error) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e)
}



pub async fn generate_info_from_bam_stream(
    rx: mpsc::Receiver<ParseOutput>,
    duplicate_clusters: &Arc<HashMap<String, ClusterInfo>>,
    min_contig_size: usize,
    thread_pool: &Arc<rayon::ThreadPool>,
) -> Result<(
    HashMap<String, String>,           // read2contig
    HashMap<String, u64>,              // contig_stats (cluster-adjusted)
    HashMap<String, usize>,            // contig_unique_counts
)> {
    let byte_stream = ReceiverStream::new(rx).map(|item| match item {
        ParseOutput::Bytes(arc) => Ok(Bytes::from((*arc).clone())),
        _ => Err(anyhow_to_io(anyhow!("BAM stream received non-Bytes variant"))),
    });

    let stream_reader = StreamReader::new(byte_stream);
    let mut bam_reader = BamAsyncReader::new(stream_reader);

    let header = bam_reader.read_header().await.map_err(|e| anyhow!("BAM header error: {}", e))?;

    // Shared structures
    let read2contig = Arc::new(Mutex::new(HashMap::with_capacity(20_000_000)));
    let contig_stats = Arc::new(Mutex::new(HashMap::with_capacity(1024)));
    let contig_unique_counts = Arc::new(Mutex::new(HashMap::with_capacity(1024)));
    let seen_reads = Arc::new(Mutex::new(HashSet::with_capacity(20_000_000)));

    // Collect all records
    let mut records: Vec<BamRecord> = Vec::new();
    let mut record = BamRecord::default();
    loop {
        match bam_reader.read_record(&mut record).await {
            Ok(0) => break,
            Ok(_) => records.push(record.clone()),
            Err(e) => return Err(anyhow!("BAM record read error: {}", e)),
        }
    }

    // Chunk and process in parallel
    let chunk_size = (records.len() / 64).max(1); // Adjust divisor based on cores
    let chunks: Vec<_> = records.chunks(chunk_size).collect();

    thread_pool.install(|| {
        chunks.par_iter().for_each(|chunk| {
            let mut local_read2contig = HashMap::new();
            let mut local_contig_stats = HashMap::new();
            let mut local_contig_unique_counts = HashMap::new();
            let mut local_seen = HashSet::new();

            for record in chunk.iter() {  // ← .iter() gives &BamRecord
                let flags = record.flags();
                if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                    continue;
                }

                let read_name = record
                    .name()
                    .and_then(|b| std::str::from_utf8(b.as_ref()).ok())
                    .unwrap_or("*")
                    .to_string();

                if !local_seen.insert(read_name.clone()) {
                    continue;
                }

                let rid = match record.reference_sequence_id() {
                    Some(Ok(rid)) => rid,
                    _ => continue,
                };

                let contig_name = match header.reference_sequences().get_index(rid) {
                    Some((name, _)) => name.to_string(),
                    None => continue,
                };

                local_read2contig.insert(read_name.clone(), contig_name.clone());

                let cluster_size = duplicate_clusters
                    .get(&read_name)
                    .map(|cluster| cluster.size)
                    .unwrap_or(1u64);

                *local_contig_stats.entry(contig_name.clone()).or_insert(0) += cluster_size;
                *local_contig_unique_counts.entry(contig_name).or_insert(0) += 1;
            }

            // Merge local into global (lock once per chunk)
            {
                let mut global_read2contig = read2contig.lock().unwrap();
                let mut global_contig_stats = contig_stats.lock().unwrap();
                let mut global_contig_unique_counts = contig_unique_counts.lock().unwrap();
                let mut global_seen = seen_reads.lock().unwrap();

                global_read2contig.extend(local_read2contig);
                for (k, v) in local_contig_stats {
                    *global_contig_stats.entry(k).or_insert(0) += v;
                }
                for (k, v) in local_contig_unique_counts {
                    *global_contig_unique_counts.entry(k).or_insert(0) += v;
                }
                global_seen.extend(local_seen);
            }
        });
    });

    // Extract final results
    let read2contig = Arc::try_unwrap(read2contig).unwrap().into_inner().unwrap();
    let mut contig_stats = Arc::try_unwrap(contig_stats).unwrap().into_inner().unwrap();
    let mut contig_unique_counts = Arc::try_unwrap(contig_unique_counts).unwrap().into_inner().unwrap();

    // Apply min unique read filter
    contig_stats.retain(|contig_name, _| {
        contig_unique_counts.get(contig_name).copied().unwrap_or(0) >= min_contig_size
    });

    contig_unique_counts.retain(|contig_name, _| contig_stats.contains_key(contig_name));

    Ok((read2contig, contig_stats, contig_unique_counts))
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
