use std::io::{self, Write};
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tokio::sync::{broadcast, oneshot};
use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncBufReadExt, AsyncReadExt, AsyncWriteExt, BufReader, BufWriter};
use tokio::process::{Child, ChildStdout, Command};
use tokio::time::{sleep, Duration, Instant};
use tokio_stream::{Stream, StreamExt};
use tokio_stream::wrappers::BroadcastStream;
use crate::utils::fastx::SequenceRecord;
use tokio::sync::broadcast::{Receiver, Sender};
use std::pin::Pin;





pub trait ToBytes {
    fn to_bytes(&self) -> std::io::Result<Vec<u8>>;
}

impl<T> ToBytes for T
where
    T: AsRef<[u8]>,
{
    fn to_bytes(&self) -> io::Result<Vec<u8>> {
        let mut bytes = self.as_ref().to_vec();
        bytes.push(b'\n');
        Ok(bytes)
    }
}

impl ToBytes for SequenceRecord {
    fn to_bytes(&self) -> io::Result<Vec<u8>> {
        let mut buffer = Vec::new();
        match self {
            SequenceRecord::Fastq { id, desc, seq, qual } => {
                if let Some(desc) = desc {
                    writeln!(buffer, "@{} {}", id, desc)?;
                } else {
                    writeln!(buffer, "@{}", id)?;
                }
                buffer.extend_from_slice(seq);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(qual);
                buffer.push(b'\n');
            }
            SequenceRecord::Fasta { id, desc, seq } => {
                if let Some(desc) = desc {
                    writeln!(buffer, ">{} {}", id, desc)?;
                } else {
                    writeln!(buffer, ">{}", id)?;
                }
                buffer.extend_from_slice(seq);
                buffer.push(b'\n');
            }
        }
        Ok(buffer)
    }
}


pub async fn t_junction<S>(
    input: S,
    n_outputs: usize,
    stall_threshold: u64,
    stream_sleep_ms: Option<u64>,
) -> Result<(Vec<BroadcastStream<SequenceRecord>>, oneshot::Receiver<()>)>
where
    S: Stream<Item = SequenceRecord> + Unpin + Send + 'static,
{
    let (tx, _) = tokio::sync::broadcast::channel(10000);
    let output_rxs: Vec<_> = (0..n_outputs)
        .map(|_| BroadcastStream::new(tx.subscribe()))
        .collect();
    let (done_tx, done_rx) = oneshot::channel();

    // Pin the input stream to ensure it is Unpin and Send
    let mut input = Box::pin(input);

    tokio::spawn(async move {
        let mut count = 0;
        while let Some(record) = input.next().await {
            if let Err(e) = tx.send(record) {
                eprintln!("Broadcast error at record {}: {}", count, e);
                // Continue to ensure all records are sent to remaining subscribers
            }
            count += 1;

            // Log stalls
            if count % stall_threshold == 0 {
                eprintln!("Broadcast stall detected at {} records", count);
                if let Some(sleep_ms) = stream_sleep_ms {
                    tokio::time::sleep(std::time::Duration::from_millis(sleep_ms)).await;
                }
            }
        }
        eprintln!("Finished broadcasting {} records", count);
        let _ = done_tx.send(());
    });

    Ok((output_rxs, done_rx))
}

pub async fn stream_to_cmd(
    mut stream: BroadcastStream<SequenceRecord>,
    cmd_tag: &str,
    cmd_args: Vec<&str>,
) -> Result<Child> {
    eprintln!("Spawning {} with args: {:?}", cmd_tag, cmd_args);
    let mut cmd = Command::new(cmd_tag);
    cmd.args(&cmd_args)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    let mut child = cmd.spawn().map_err(|e| anyhow!("Failed to spawn {}: {}", cmd_tag, e))?;

    eprintln!("Spawned {}: stdin={:?}, stdout={:?}, stderr={:?}",
              cmd_tag,
              child.stdin.is_some(),
              child.stdout.is_some(),
              child.stderr.is_some()
    );

    let stdin = child.stdin.take().ok_or_else(|| anyhow!("Failed to open stdin for {}", cmd_tag))?;
    let cmd_tag_owned = cmd_tag.to_string(); // Clone to String for 'static lifetime
    tokio::spawn(async move {
        let mut stdin = tokio::io::BufWriter::new(stdin);
        let mut count = 0;
        while let Some(Ok(record)) = stream.next().await {
            if let SequenceRecord::Fastq { id, desc, seq, qual } = record {
                if let Err(e) = stdin.write_all(b"@").await {
                    eprintln!("Error writing ID to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing ID to {} stdin: {}", cmd_tag_owned, e));
                }
                if let Err(e) = stdin.write_all(id.as_bytes()).await {
                    eprintln!("Error writing ID to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing ID to {} stdin: {}", cmd_tag_owned, e));
                }
                if let Some(desc) = desc {
                    if let Err(e) = stdin.write_all(b" ").await {
                        eprintln!("Error writing desc to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                        return Err(anyhow!("Error writing desc to {} stdin: {}", cmd_tag_owned, e));
                    }
                    if let Err(e) = stdin.write_all(desc.as_bytes()).await {
                        eprintln!("Error writing desc to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                        return Err(anyhow!("Error writing desc to {} stdin: {}", cmd_tag_owned, e));
                    }
                }
                if let Err(e) = stdin.write_all(b"\n").await {
                    eprintln!("Error writing newline to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing newline to {} stdin: {}", cmd_tag_owned, e));
                }
                if let Err(e) = stdin.write_all(&seq).await {
                    eprintln!("Error writing seq to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing seq to {} stdin: {}", cmd_tag_owned, e));
                }
                if let Err(e) = stdin.write_all(b"\n+\n").await {
                    eprintln!("Error writing plus to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing plus to {} stdin: {}", cmd_tag_owned, e));
                }
                if let Err(e) = stdin.write_all(&qual).await {
                    eprintln!("Error writing qual to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing qual to {} stdin: {}", cmd_tag_owned, e));
                }
                if let Err(e) = stdin.write_all(b"\n").await {
                    eprintln!("Error writing newline to {} stdin at record {}: {}", cmd_tag_owned, count, e);
                    return Err(anyhow!("Error writing newline to {} stdin: {}", cmd_tag_owned, e));
                }
                count += 1;
            }
        }
        if let Err(e) = stdin.flush().await {
            eprintln!("Error flushing {} stdin: {}", cmd_tag_owned, e);
            return Err(anyhow!("Error flushing {} stdin: {}", cmd_tag_owned, e));
        }
        eprintln!("Finished writing {} records to {} stdin", count, cmd_tag_owned);
        Ok::<(), anyhow::Error>(())
    });

    Ok(child)
}
/// Takes output from stream_to_cmd and outputs it as seq_io records.
///
/// # Arguments
///
/// * `sydout' - Child process stdout.
///
/// # Returns
/// Result<BroadcastStream<SequenceRecord>>
// pub async fn parse_child_stdout_to_fastq(
//     stdout: ChildStdout,
// ) -> Result<BroadcastStream<SequenceRecord>> {
//     let (tx, rx) = broadcast::channel(100);
//     let mut lines = BufReader::with_capacity(1024 * 1024, stdout).lines();
//
//     tokio::spawn(async move {
//         while let Some(line) = lines.next_line().await.unwrap_or(None) {
//             if line.starts_with('@') {
//                 let id = line[1..].to_string();
//                 let seq = match lines.next_line().await {
//                     Ok(Some(s)) => s,
//                     _ => {
//                         eprintln!("Missing sequence line for FASTQ record");
//                         break;
//                     }
//                 };
//                 let plus = match lines.next_line().await {
//                     Ok(Some(p)) if p == "+" => p,
//                     _ => {
//                         eprintln!("Invalid or missing '+' line for FASTQ record");
//                         break;
//                     }
//                 };
//                 let qual = match lines.next_line().await {
//                     Ok(Some(q)) => q,
//                     _ => {
//                         eprintln!("Missing quality line for FASTQ record");
//                         break;
//                     }
//                 };
//                 if seq.len() != qual.len() {
//                     eprintln!("Sequence and quality lengths mismatch");
//                     break;
//                 }
//                 let record = SequenceRecord::Fastq {
//                     id,
//                     desc: None,
//                     seq: seq.into_bytes(),
//                     qual: qual.into_bytes(),
//                 };
//                 if tx.send(record).is_err() {
//                     eprintln!("Failed to send SequenceRecord");
//                     break;
//                 }
//             }
//         }
//         // Lines dropped, channel remains open
//     });
//
//     Ok(BroadcastStream::new(rx))
// }

pub async fn parse_child_stdout_to_fastq<R: AsyncRead + Unpin>(
    reader: R,
    sender: tokio::sync::mpsc::Sender<SequenceRecord>,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(1024 * 1024, reader); // 1MB buffer
    let mut count = 0;
    let mut buffer = String::new();

    loop {
        // Read ID line
        buffer.clear();
        let bytes_read = match reader.read_line(&mut buffer).await {
            Ok(bytes) => bytes,
            Err(e) => {
                eprintln!("Error reading ID line at record {}: {}", count, e);
                return Err(anyhow!("Failed to read ID line: {}", e));
            }
        };
        if bytes_read == 0 {
            println!("Reached EOF at record {}", count);
            break; // EOF
        }
        let id_line = buffer.trim_end();
        if !id_line.starts_with('@') {
            eprintln!("Invalid FASTQ format at record {}: expected '@', got '{}'", count, id_line);
            return Err(anyhow!("Invalid FASTQ format: expected '@', got '{}'", id_line));
        }

        // Parse ID and desc
        let (id, desc) = match id_line[1..].split_once(' ') {
            Some((id, desc)) => (id.to_string(), Some(desc.to_string())),
            None => (id_line[1..].to_string(), None),
        };

        // Read sequence
        buffer.clear();
        if let Err(e) = reader.read_line(&mut buffer).await {
            eprintln!("Error reading sequence at record {}: {}", count, e);
            return Err(anyhow!("Failed to read sequence: {}", e));
        }
        let seq = buffer.trim_end().as_bytes().to_vec();
        if seq.is_empty() {
            eprintln!("Missing sequence at record {}", count);
            return Err(anyhow!("Missing sequence for record {}", count));
        }

        // Read plus line
        buffer.clear();
        if let Err(e) = reader.read_line(&mut buffer).await {
            eprintln!("Error reading plus line at record {}: {}", count, e);
            return Err(anyhow!("Failed to read plus line: {}", e));
        }
        let plus = buffer.trim_end();
        if plus != "+" {
            eprintln!("Invalid FASTQ format at record {}: expected '+', got '{}'", count, plus);
            return Err(anyhow!("Invalid FASTQ format: expected '+', got '{}'", plus));
        }

        // Read quality
        buffer.clear();
        if let Err(e) = reader.read_line(&mut buffer).await {
            eprintln!("Error reading quality at record {}: {}", count, e);
            return Err(anyhow!("Failed to read quality: {}", e));
        }
        let qual = buffer.trim_end().as_bytes().to_vec();
        if qual.is_empty() {
            eprintln!("Missing quality at record {}", count);
            return Err(anyhow!("Missing quality for record {}", count));
        }

        // Validate sequence and quality lengths
        if seq.len() != qual.len() {
            eprintln!("Sequence length ({}) != quality length ({}) at record {}", seq.len(), qual.len(), count);
            return Err(anyhow!("Sequence length ({}) != quality length ({}) for record {}", seq.len(), qual.len(), count));
        }

        // Create SequenceRecord::Fastq
        let record = SequenceRecord::Fastq {
            id,
            desc,
            seq,
            qual,
        };

        // Log first 16 records
        if count < 16 {
            println!(
                "Parsed record {}: id={}, seq_len={}, qual_len={}, desc={:?}",
                count,
                record.id(),
                record.seq().len(),
                record.qual().len(),
                record.desc()
            );
        }

        // Send record
        if let Err(e) = sender.send(record).await {
            eprintln!("Failed to send record {}: {}", count, e);
            return Err(anyhow!("Failed to send record: {}", e));
        }
        count += 1;

        // Throttle to prevent channel backpressure
        if count % 1000 == 0 {
            sleep(Duration::from_millis(1)).await;
        }
    }

    println!("Total records parsed: {}", count);
    Ok(())
}


pub async fn parse_child_stdout_to_bytes(
    stdout: ChildStdout,
) -> Result<BroadcastStream<Vec<u8>>> {
    let (tx, rx) = broadcast::channel(100);
    let mut reader = BufReader::with_capacity(1024 * 1024, stdout);

    tokio::spawn(async move {
        let mut buffer = vec![0u8; 256 * 1024];
        loop {
            match reader.read(&mut buffer).await {
                Ok(0) => break, // EOF
                Ok(n) => {
                    let chunk = buffer[..n].to_vec();
                    if tx.send(chunk).is_err() {
                        eprintln!("Failed to send byte chunk to broadcast channel");
                        break;
                    }
                }
                Err(e) => {
                    eprintln!("Error reading stdout: {}", e);
                    break;
                }
            }
        }
    });

    Ok(BroadcastStream::new(rx))
}

/// Writes a stream to a file.
/// Capture stdout and return from function.
///
/// # Arguments
///
/// * `input_stream' - BroadcastStream<T>
/// * 'output_opatj' - PathBuf for output file

///
/// # Returns
/// Result<()>
// pub async fn stream_to_file<T>(
//     mut input_stream: BroadcastStream<T>,
//     output_path: PathBuf,
// ) -> anyhow::Result<()>
// where
//     T: ToBytes + Clone + Send + 'static,
// {
//     let file = File::create(&output_path)
//         .await
//         .map_err(|e| anyhow!("Failed to create file {}: {}", output_path.display(), e))?;
//     let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
//
//     let mut total_bytes = 0;
//     let mut last_progress = Instant::now();
//
//     while let Some(result) = input_stream.next().await {
//         if last_progress.elapsed() > Duration::from_secs(15) {
//             eprintln!("File writer stall detected at {} bytes", total_bytes);
//             last_progress = Instant::now();
//         }
//
//         let item = result.map_err(|e| anyhow!("Broadcast stream error: {}", e))?;
//         let bytes = item.to_bytes().map_err(|e| anyhow!("Failed to convert to bytes: {}", e))?;
//         writer
//             .write_all(&bytes)
//             .await
//             .map_err(|e| anyhow!("Failed to write to file {}: {}", output_path.display(), e))?;
//         total_bytes += bytes.len();
//     }
//
//     writer
//         .flush()
//         .await
//         .map_err(|e| anyhow!("Failed to flush file {}: {}", output_path.display(), e))?;
//
//     eprintln!(
//         "File stream for {} completed, wrote {} bytes",
//         output_path.display(),
//         total_bytes
//     );
//
//     Ok(())
// }
pub async fn stream_to_file(
    mut rx: tokio::sync::mpsc::Receiver<SequenceRecord>,
    path: PathBuf,
) -> Result<()> {
    let mut file = File::create(&path).await?;
    let mut count = 0;
    while let Some(record) = rx.recv().await {
        if count < 16 {
            println!("Writing record {}: id={}, seq_len={}", count, record.id(), record.seq().len());
        }
        if let SequenceRecord::Fastq { id, desc, seq, qual } = record {
            file.write_all(b"@").await?;
            file.write_all(id.as_bytes()).await?;
            if let Some(desc) = desc {
                file.write_all(b" ").await?;
                file.write_all(desc.as_bytes()).await?;
            }
            file.write_all(b"\n").await?;
            file.write_all(&seq).await?;
            file.write_all(b"\n+\n").await?;
            file.write_all(&qual).await?;
            file.write_all(b"\n").await?;
            count += 1;
        }
    }
    if count < 47108 {
        eprintln!("stream_to_file stopped early at {} sequences, expected ~47108", count);
    }
    println!("Wrote {} sequences to {}", count, path.display());
    Ok(())
}

pub async fn read_child_stdout(mut child: Child) -> io::Result<()> {
    let mut stdout = child
        .stdout
        .take()
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to open stdout"))?;

    let mut buffer = [0; 1024];
    loop {
        match stdout.read(&mut buffer).await {
            Ok(0) => {
                println!("Child process stdout closed");
                break;
            }
            Ok(n) => {
                let output = String::from_utf8_lossy(&buffer[..n]);
                print!("Child stdout: {}", output);
            }
            Err(e) => {
                eprintln!("Failed to read from stdout: {}", e);
                break;
            }
        }
    }

    let status = child.wait().await?;
    println!("Child exited with: {}", status);

    Ok(())
}

pub async fn stream_bytes_to_file(mut rx: tokio_stream::wrappers::BroadcastStream<Vec<u8>>, path: PathBuf) -> Result<()> {
    let mut file = File::create(&path).await?;
    let mut total_bytes = 0;
    while let Some(Ok(bytes)) = rx.next().await {
        file.write_all(&bytes).await?;
        total_bytes += bytes.len();
    }
    println!("Wrote {} bytes to {}", total_bytes, path.display());
    Ok(())
}

pub async fn debug_fastp_stdout<R: AsyncRead + Unpin>(reader: R, path: PathBuf) -> Result<()> {
    let mut file = File::create(&path).await?;
    let mut reader = BufReader::new(reader);
    tokio::io::copy(&mut reader, &mut file).await?;
    println!("Wrote fastp raw stdout to {}", path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tokio::process::Command;

    #[tokio::test]
    async fn test_parse_child_stdout_to_fastq() -> Result<()> {
        let mut cmd = Command::new("echo");
        cmd.arg("@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nHHHH\n");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let stdout = child.stdout.take().unwrap();
        let stream = parse_child_stdout_to_fastq(stdout).await?;
        let mut stream = stream;
        let mut records = Vec::new();
        while let Some(Ok(record)) = stream.next().await {
            records.push(record);
        }
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "read1");
        assert_eq!(records[0].seq(), b"ATCG");
        assert_eq!(records[1].id(), "read2");
        assert_eq!(records[1].seq(), b"GCTA");
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_child_stdout_to_bytes() -> Result<()> {
        let mut cmd = Command::new("echo");
        cmd.arg("test data");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let stdout = child.stdout.take().unwrap();
        let stream = parse_child_stdout_to_bytes(stdout).await?;
        let mut stream = stream;
        let mut chunks = Vec::new();
        while let Some(Ok(chunk)) = stream.next().await {
            chunks.push(chunk);
        }
        let data = chunks.concat();
        assert!(String::from_utf8_lossy(&data).contains("test data"));
        Ok(())
    }
}