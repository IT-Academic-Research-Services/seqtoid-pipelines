use anyhow::{anyhow, Result};
use bytes::Bytes;
use futures::StreamExt;
// use noodles::bam::r#async::io::Reader as BamAsyncReader;
// use noodles::bam::record::Record;
// use noodles::sam::Header;
use std::collections::{HashMap, HashSet};
use std::io;
use std::sync::Arc;
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tokio_util::io::StreamReader;
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


// fn anyhow_to_io(e: anyhow::Error) -> io::Error {
//     io::Error::new(io::ErrorKind::InvalidData, e)
// }

// pub async fn generate_info_from_bam_stream(
//     rx: mpsc::Receiver<ParseOutput>,
//     duplicate_cluster_sizes: &HashMap<String, u64>,
//     min_contig_size: usize,
// ) -> Result<HashMap<String, u64>> {
// 
//     let byte_stream = ReceiverStream::new(rx).map(|item| match item {
//         ParseOutput::Bytes(arc) => Ok(Bytes::from((*arc).clone())), // Clone Vec<u8>, then to Bytes
//         _ => Err(anyhow_to_io(anyhow!(
//             "BAM stream received non-Bytes variant — data loss prevented"
//         ))),
//     });
// 
//     let stream_reader = StreamReader::new(byte_stream);
//     let mut bam_reader = BamAsyncReader::new(stream_reader);
// 
//     let header = bam_reader
//         .read_header()
//         .await
//         .map_err(|e| anyhow!("BAM header error: {e}"))?;
// 
//     let mut contig_stats = HashMap::with_capacity(1024);
//     let mut contig_unique_counts = HashMap::with_capacity(1024);
//     let mut seen_reads = HashSet::with_capacity(20_000_000); // Scales to 100M+ reads, fits 1.5TB RAM
//     let mut record = Record::default(); // Reused: zero per-record allocs
// 
//     loop {
//         match bam_reader.read_record(&mut record).await {
//             Ok(0) => break, // EOF
//             Ok(_) => {
//                 let flags = record.flags();
//                 if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
//                     continue;
//                 }
// 
//                 let read_name = record
//                     .name()
//                     .and_then(|b| std::str::from_utf8(b.as_ref()).ok())
//                     .unwrap_or("*")
//                     .to_string();
// 
//                 if !seen_reads.insert(read_name.clone()) {
//                     continue;
//                 }
// 
//                 let rid = record
//                     .reference_sequence_id()
//                     .ok_or_else(|| anyhow!("Missing reference_sequence_id"))??;
// 
//                 let contig_name = header
//                     .reference_sequences()
//                     .get_index(rid)
//                     .ok_or_else(|| anyhow!("Invalid reference ID {rid:?}"))?
//                     .0
//                     .to_string();
// 
//                 let cluster_size = duplicate_cluster_sizes.get(&read_name).copied().unwrap_or(1);
//                 *contig_stats.entry(contig_name.clone()).or_insert(0u64) += cluster_size;
//                 *contig_unique_counts.entry(contig_name).or_insert(0usize) += 1;
//             }
//             Err(e) => return Err(anyhow!("BAM record error: {e}")), // Propagate: no silent drops
//         }
//     }
// 
//     // Retain contigs with >= min_contig_size unique reads
//     contig_stats.retain(|contig_name, _| {
//         contig_unique_counts.get(contig_name).copied().unwrap_or(0) >= min_contig_size
//     });
// 
//     Ok(contig_stats)
// }

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
