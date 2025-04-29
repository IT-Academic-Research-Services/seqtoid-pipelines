use std::fs;
use std::path::Path;
use anyhow::anyhow;
use anyhow::Result;
use std::path::PathBuf;
use tokio::fs::File;
use tokio::task;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt, AsyncBufReadExt, BufReader};
use tokio::process::{Child, ChildStdout, Command};
use tokio::sync::{broadcast, oneshot};
use tokio::time::{self, Duration};
use tokio_stream::{Stream, StreamExt};
use tokio_stream::wrappers::BroadcastStream;
use tokio::sync::mpsc;
use crate::utils::fastx::{SequenceRecord, fastx_generator};
use tokio::task::JoinHandle;


pub trait ToBytes {
    fn to_bytes(&self) -> Result<Vec<u8>>;
}

impl ToBytes for Vec<u8> {
    fn to_bytes(&self) -> Result<Vec<u8>> {
        Ok(self.clone())
    }
}

impl ToBytes for SequenceRecord {
    fn to_bytes(&self) -> Result<Vec<u8>> {
        let mut buffer = Vec::new();
        match self {
            SequenceRecord::Fastq { id, desc, seq, qual } => {
                if let Some(desc) = desc {
                    buffer.extend_from_slice(format!("@{} {}\n", id, desc).as_bytes());
                } else {
                    buffer.extend_from_slice(format!("@{}\n", id).as_bytes());
                }
                buffer.extend_from_slice(seq);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(qual);
                buffer.push(b'\n');
            }
            SequenceRecord::Fasta { id, desc, seq } => {
                if let Some(desc) = desc {
                    buffer.extend_from_slice(format!(">{} {}\n", id, desc).as_bytes());
                } else {
                    buffer.extend_from_slice(format!(">{}\n", id).as_bytes());
                }
                buffer.extend_from_slice(seq);
                buffer.push(b'\n');
            }
        }
        Ok(buffer)
    }
}

/// Generates any number of output streams from a single input stream.
/// The streams are all asynchronous.
///
/// # Arguments
///
/// * `input_stream`: An asynchronous stream yielding items of type `T`.
/// * `num_streams`: Number of output streams to generate.
/// * `stall_threshold_secs`: Seconds before logging a stall.
/// * `sleep_duration_ms`: Optional milliseconds to sleep per item (None for no sleep).
///
/// # Returns
/// A `Result` containing a tuple of:
/// - A vector of `BroadcastStream<T>` for downstream processing.
/// - A `oneshot::Receiver<()>` to await task completion.
pub async fn t_junction<S, T>(
    input: S,
    n_outputs: usize,
    buffer_size: usize,
    stall_threshold: u64,
    stream_sleep_ms: Option<u64>,
) -> Result<(Vec<BroadcastStream<T>>, oneshot::Receiver<Result<(), anyhow::Error>>)>
where
    S: Stream<Item = T> + Unpin + Send + 'static,
    T: ToBytes + Clone + Send + Sync + 'static,
{
    let (done_tx, done_rx) = oneshot::channel::<Result<(), anyhow::Error>>();

    let (tx, _) = broadcast::channel(buffer_size);
    let output_rxs: Vec<_> = (0..n_outputs)
        .map(|_| BroadcastStream::new(tx.subscribe()))
        .collect();

    let mut input = Box::pin(input);

    tokio::spawn(async move {
        // Check for zero subscribers
        if n_outputs == 0 {
            let _ = done_tx.send(Err(anyhow!("No subscribers: cannot process stream")));
            return;
        }

        let mut count = 0;
        while let Some(item) = input.next().await {
            // Check for lagging subscribers
            if tx.receiver_count() == 0 {
                panic!("All subscribers dropped before stream completion");
            }
            // Send item and handle errors
            match tx.send(item) {
                Ok(_) => (),
                Err(broadcast::error::SendError(_)) => {
                    panic!("Broadcast channel lagged: failed to send item {}", count + 1);
                }
            }
            count += 1;
            if count % stall_threshold == 0 {
                if let Some(sleep_ms) = stream_sleep_ms {
                    tokio::time::sleep(Duration::from_millis(sleep_ms)).await;
                }
            }
        }
        // Verify all items were sent to active subscribers
        if tx.receiver_count() > 0 {
            let _ = done_tx.send(Ok(()));
        } else {
            panic!("No active subscribers at stream completion");
        }
    });

    Ok((output_rxs, done_rx))
}

/// Asynchronously spawn an external process and feed it a stream as stdin.
/// Capture stdout and return from function.
///
/// # Arguments
///
/// * `stream' - Receiver stream: tokio::mpsc
/// * 'command' - command to shell out
/// * 'args' = args for shelled out command
///
/// # Returns
/// tokio::process::Command containing stdout and stderr
pub async fn stream_to_cmd<T: ToBytes + Clone + Send + Sync + 'static>(
    mut rx: BroadcastStream<T>,
    cmd_tag: &str,
    args: Vec<&str>,
) -> Result<(Child, JoinHandle<Result<(), anyhow::Error>>)> {
    let cmd_tag_owned = cmd_tag.to_string();
    let mut child = Command::new(&cmd_tag_owned)
        .args(&args)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow!("Failed to spawn {}: {}", cmd_tag_owned, e))?;

    let stdin = child
        .stdin
        .take()
        .ok_or_else(|| anyhow!("Failed to open stdin for {}", cmd_tag_owned))?;

    let task = tokio::spawn(async move {
        let mut writer = tokio::io::BufWriter::new(stdin);
        while let Some(result) = rx.next().await {
            match result {
                Ok(item) => {
                    let bytes = item.to_bytes()?;
                    writer.write_all(&bytes).await?;
                    writer.flush().await?;
                }
                Err(e) => return Err(anyhow!("Broadcast stream error: {}", e)),
            }
        }
        writer.flush().await?;
        writer.shutdown().await?;
        Ok(())
    });

    Ok((child, task))
}

/// Takes output from stream_to_cmd and outputs it as seq_io records.
///
/// # Arguments
///
/// * `sydout' - Child process stdout.
///
/// # Returns
/// Result<BroadcastStream<SequenceRecord>>
pub async fn parse_child_stdout_to_fastq<R: AsyncRead + Unpin>(
    reader: R,
    sender: tokio::sync::mpsc::Sender<SequenceRecord>,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(1024 * 1024, reader);
    let mut buffer = String::new();
    let mut count = 0;

    loop {
        buffer.clear();
        let bytes_read = reader.read_line(&mut buffer).await?;
        if bytes_read == 0 {
            break;
        }
        let id_line = buffer.trim_end();
        if !id_line.starts_with('@') {
            return Err(anyhow::anyhow!("Invalid FASTQ format: expected '@', got '{}'", id_line));
        }

        let (id, desc) = match id_line[1..].split_once(' ') {
            Some((id, desc)) => (id.to_string(), Some(desc.to_string())),
            None => (id_line[1..].to_string(), None),
        };

        buffer.clear();
        reader.read_line(&mut buffer).await?;
        let seq = buffer.trim_end().as_bytes().to_vec();
        if seq.is_empty() {
            return Err(anyhow::anyhow!("Missing sequence"));
        }

        buffer.clear();
        reader.read_line(&mut buffer).await?;
        let plus = buffer.trim_end();
        if plus != "+" {
            return Err(anyhow::anyhow!("Invalid FASTQ format: expected '+', got '{}'", plus));
        }

        buffer.clear();
        reader.read_line(&mut buffer).await?;
        let qual = buffer.trim_end().as_bytes().to_vec();
        if qual.is_empty() {
            return Err(anyhow::anyhow!("Missing quality"));
        }

        if seq.len() != qual.len() {
            return Err(anyhow::anyhow!("Sequence length ({}) != quality length ({})", seq.len(), qual.len()));
        }

        let record = SequenceRecord::Fastq {
            id,
            desc,
            seq,
            qual,
        };

        sender.send(record).await?;
        count += 1;

        if count % 1000 == 0 {
            tokio::time::sleep(Duration::from_millis(1)).await;
        }
    }

    Ok(())
}

/// Takes output from stream_to_cmd and outputs it as a byte stream.
///
/// # Arguments
///
/// * `stdout' - Child process stdout.
///
/// # Returns
/// Result<BroadcastStream<Vec<u8>>>
pub async fn parse_child_stdout_to_bytes(stdout: ChildStdout) -> Result<BroadcastStream<Vec<u8>>> {
    let (tx, rx) = broadcast::channel(100);
    let mut reader = BufReader::with_capacity(1024 * 1024, stdout);

    tokio::spawn(async move {
        let mut buffer = vec![0u8; 256 * 1024];
        loop {
            match reader.read(&mut buffer).await {
                Ok(0) => break,
                Ok(n) => {
                    let chunk = buffer[..n].to_vec();
                    let _ = tx.send(chunk);
                }
                Err(_) => break,
            }
        }
    });

    Ok(BroadcastStream::new(rx))
}

/// Writes a stream of SequenceRecords to a file.
/// Capture stdout and return from function.
///
/// # Arguments
///
/// * `rx' - tokio::sync::mpsc::Receiver<SequenceRecord>
/// * 'path' - PathBuf for output file

///
/// # Returns
/// Result<()>
pub async fn stream_sequence_records_to_file(mut rx: tokio::sync::mpsc::Receiver<SequenceRecord>, path: PathBuf) -> Result<()> {
    let mut file = File::create(&path).await?;
    while let Some(record) = rx.recv().await {
        let bytes = record.to_bytes()?;
        file.write_all(&bytes).await?;
    }
    file.flush().await?;
    Ok(())
}

/// Writes a stream of bytes to a file.
/// Capture stdout and return from function.
///
/// # Arguments
///
/// * `rx' - BBroadcastStream<Vec<u8>>
/// * 'path' - PathBuf for output file

///
/// # Returns
/// Result<()>
pub async fn stream_bytes_to_file(mut rx: BroadcastStream<Vec<u8>>, path: PathBuf) -> Result<()> {
    let mut file = File::create(&path).await?;
    while let Some(Ok(bytes)) = rx.next().await {
        file.write_all(&bytes).await?;
    }
    file.flush().await?;
    Ok(())
}

/// A sink for testing output streams.
/// Reads child stdout to screen.
///
/// # Arguments
///
/// * `child' - Child process from stream_to_cmd.
///
/// # Returns
/// io::Result<()>
pub async fn read_child_stdout(mut child: Child) -> Result<()> {
    let mut stdout = child
        .stdout
        .take()
        .ok_or_else(|| anyhow::anyhow!("Failed to open stdout"))?;

    let mut buffer = [0; 1024];
    loop {
        match stdout.read(&mut buffer).await {
            Ok(0) => break,
            Ok(_) => {
                print!("{}", String::from_utf8_lossy(&buffer));
            },
            Err(_) => break,
        }
    }

    child.wait().await?;
    Ok(())
}


/// Reads child stdout to a Vec<String>.
///
/// # Arguments
///
/// * `child' - Child process from stream_to_cmd.
///
/// # Returns
/// Result<Vec<String>>
pub async fn read_child_stdout_to_vec(mut child: Child) -> Result<Vec<String>> {
    let mut stdout = child
        .stdout
        .take()
        .ok_or_else(|| anyhow::anyhow!("Failed to open stdout"))?;

    let mut reader = BufReader::new(&mut stdout);
    let mut lines = Vec::new();
    let mut buffer = String::new();

    loop {
        buffer.clear();
        match reader.read_line(&mut buffer).await {
            Ok(0) => break, // EOF reached
            Ok(_) => {
                // Remove trailing newline and store the line
                let line = buffer.trim_end().to_string();
                if !line.is_empty() {
                    lines.push(line);
                }
            }
            Err(e) => return Err(anyhow::anyhow!("Failed to read stdout: {}", e)),
        }
    }

    // Wait for the child process to finish
    let status = child.wait().await?;
    if !status.success() {
        return Err(anyhow::anyhow!("Child process exited with non-zero status: {}", status));
    }

    Ok(lines)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tokio::process::Command;

    #[tokio::test]
    async fn test_t_junction_zero_streams() -> Result<()> {
        let stream = fastx_generator(10, 143, 35.0, 3.0);
        let (outputs, done_rx) = t_junction(stream, 0, 50000, 10000, Some(1)).await?;
        assert_eq!(outputs.len(), 0);
        let result = done_rx.await?;
        assert!(result.is_err(), "Expected error for zero subscribers");
        assert_eq!(
            result.unwrap_err().to_string(),
            "No subscribers: cannot process stream"
        );
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_two_records() -> Result<()> {
        let records = vec![
            SequenceRecord::Fastq {
                id: "read1".to_string(),
                desc: None,
                seq: b"ATCG".to_vec(),
                qual: b"IIII".to_vec(),
            },
            SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: b"GCTA".to_vec(),
                qual: b"HHHH".to_vec(),
            },
        ];
        let stream = tokio_stream::iter(records);
        let (mut outputs, done_rx) = t_junction(stream, 2, 50000, 10000, Some(1)).await?;
        let mut output1 = outputs.pop().unwrap();
        let mut output2 = outputs.pop().unwrap();
        let mut records1 = Vec::new();
        let mut records2 = Vec::new();
        while let Some(Ok(record)) = output1.next().await {
            records1.push(record);
        }
        while let Some(Ok(record)) = output2.next().await {
            records2.push(record);
        }
        assert_eq!(records1.len(), 2);
        assert_eq!(records2.len(), 2);
        assert_eq!(records1[0].id(), "read1");
        assert_eq!(records2[0].id(), "read1");
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_long_stream() -> Result<()> {
        let stream = fastx_generator(10000, 143, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 2, 50000, 10000, Some(1)).await?;
        let mut output1 = outputs.pop().unwrap();
        let mut output2 = outputs.pop().unwrap();
        let mut records1 = Vec::new();
        let mut records2 = Vec::new();
        while let Some(Ok(record)) = output1.next().await {
            records1.push(record);
        }
        while let Some(Ok(record)) = output2.next().await {
            records2.push(record);
        }
        assert_eq!(records1.len(), 10000);
        assert_eq!(records2.len(), 10000);

        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_ten_thousand_records_ten_streams() -> Result<()> {
        let stream = fastx_generator(10000, 143, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 100, 50000,10000, Some(1)).await?;
        let mut records = Vec::new();
        for _output in &outputs {
            let record :Vec<SequenceRecord> = Vec::new();
            records.push(record);
        }

        for i in 0..records.len() {
            while let Some(Ok(record)) = outputs[i].next().await {
                records[i].push(record);
            }
        }

        for i in 0..records.len() {
            assert_eq!(records[i].len(), 10000)
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_empty_stream() -> Result<()> {
        let stream = fastx_generator(0, 50, 35.0, 3.0);
        let (outputs, done_rx) = t_junction(stream, 2, 50000, 10000, Some(1)).await?;
        for mut output in outputs {
            assert!(output.next().await.is_none(), "Empty stream should yield no items");
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_single_record() -> Result<()> {
        let stream = fastx_generator(1, 50, 35.0, 3.0);
        let (outputs, done_rx) = t_junction(stream, 2, 50000, 10000, Some(1)).await?;
        let mut handles = Vec::new();
        for output in outputs {
            handles.push(task::spawn(async move {
                let mut records = Vec::new();
                let mut stream = output;
                while let Some(Ok(record)) = stream.next().await {
                    records.push(record);
                }
                Ok::<_, anyhow::Error>(records)
            }));
        }
        let all_records = time::timeout(Duration::from_secs(10), async {
            let mut all_records = Vec::new();
            for handle in handles {
                let records = handle.await??;
                all_records.push(records);
            }
            Ok::<_, anyhow::Error>(all_records)
        })
            .await??;
        for records in &all_records {
            assert_eq!(records.len(), 1, "Should have one record");
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_slow_consumer() -> Result<()> {
        let stream = fastx_generator(1000, 50, 35.0, 3.0);
        let (outputs, done_rx) = t_junction(stream, 2, 50, 10, Some(100)).await?;
        let mut handles = Vec::new();
        for (i, output) in outputs.into_iter().enumerate() {
            let handle = task::spawn(async move {
                let mut records = Vec::new();
                let mut stream = output;
                while let Some(Ok(record)) = stream.next().await {
                    records.push(record);
                    if i == 1 {
                        time::sleep(Duration::from_millis(5)).await; // Simulate slow consumer
                    }
                }
                Ok::<_, anyhow::Error>(records)
            });
            handles.push(handle);
        }
        let all_records = time::timeout(Duration::from_secs(30), async {
            let mut all_records = Vec::new();
            for handle in handles {
                let records = handle.await??;
                all_records.push(records);
            }
            Ok::<_, anyhow::Error>(all_records)
        })
            .await??;
        for records in &all_records {
            assert_eq!(records.len(), 1000, "Should have all records");
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_million_records_ten_streams() -> Result<()> {
        let num_records = 1000000;
        // 1M records, 2 streams, stall 1000, No sleep
        let stream = fastx_generator(num_records, 143, 35.0, 3.0);
        let (outputs, done_rx) = t_junction(stream, 2, 50000, 10000, Some(1)).await?;


        let mut handles = Vec::new();
        for output in outputs {
            let handle = task::spawn(async move {
                let mut records = Vec::new();
                let mut stream = output;
                while let Some(Ok(record)) = stream.next().await {
                    records.push(record);
                }
                Ok::<_, anyhow::Error>(records)
            });
            handles.push(handle);
        }
        let all_records = time::timeout(Duration::from_secs(60), async {
            let mut all_records = Vec::new();
            for handle in handles {
                let records = handle.await??;
                all_records.push(records);
            }
            Ok::<_, anyhow::Error>(all_records)
        })
            .await??;

        for (i, records) in all_records.iter().enumerate() {
            assert_eq!(
                records.len(),
                num_records,
                "Output {} should have {} records",
                i,
                num_records
            );
        }

        // Verify records match across outputs
        for i in 0..num_records {
            assert_eq!(
                all_records[0][i].id(),
                all_records[1][i].id(),
                "Record {} IDs should match",
                i
            );
            assert_eq!(
                all_records[0][i].seq(),
                all_records[1][i].seq(),
                "Record {} sequences should match",
                i
            );
            assert_eq!(
                all_records[0][i].qual(),
                all_records[1][i].qual(),
                "Record {} quality scores should match",
                i
            );
        }


        done_rx.await??;
        Ok(())
    }


    #[tokio::test]
    async fn test_stream_to_cmd_valid() -> Result<()> {
        let stream = fastx_generator(100, 50, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "cat", vec![]).await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??; // Check for stream errors
        done_rx.await??; // Check for t_junction errors
        child.wait().await?;
        assert!(!output.is_empty(), "Output should contain data");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_valid_cat() -> Result<()> {
        let stream = fastx_generator(2, 10, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "cat", vec![]).await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");
        assert!(!output.is_empty(), "Output should contain FASTQ data");
        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1"), "Output should contain first read ID");
        assert!(output_str.contains("@read2"), "Output should contain second read ID");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_invalid_command() -> Result<()> {
        let stream = fastx_generator(2, 10, 35.0, 3.0);
        let (mut outputs, _done_rx) = t_junction(stream, 1,  50000, 10000, Some(1)).await?;
        let result = stream_to_cmd(outputs.pop().unwrap(), "nonexistent_cmd", vec![]).await;
        assert!(result.is_err(), "Should fail for invalid command");
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Failed to spawn nonexistent_cmd"), "Error should mention command name");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_empty_stream() -> Result<()> {
        let stream = fastx_generator(0, 10, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "cat", vec![]).await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");
        assert!(output.is_empty(), "Output should be empty for empty stream");
        Ok(())
    }


    #[tokio::test]
    async fn test_stream_to_cmd_large_stream() -> Result<()> {
        let num_records = 10000;
        let stream = fastx_generator(num_records, 50, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "cat", vec![]).await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");

        // Parse output as FASTQ records
        let output_str = String::from_utf8_lossy(&output);
        let mut lines = output_str.lines().peekable();
        let mut record_count = 0;

        while let Some(line) = lines.next() {
            if line.starts_with('@') {
                // Expect sequence, +, and quality lines
                let seq_line = lines.next();
                let plus_line = lines.next();
                let qual_line = lines.next();
                if seq_line.is_none() || plus_line != Some("+") || qual_line.is_none() {
                    return Err(anyhow!("Invalid FASTQ format at record {}", record_count + 1));
                }
                record_count += 1;
            } else {
                return Err(anyhow!("Unexpected line in FASTQ output: {}", line));
            }
        }

        assert_eq!(
            record_count, num_records,
            "Output should contain {} records, found {}",
            num_records, record_count
        );
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_no_output() -> Result<()> {
        let stream = fastx_generator(2, 10, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "true", vec![]).await?;
        task.await??;
        done_rx.await??;
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");
        let stdout = child.stdout.take();
        assert!(stdout.is_none() || tokio::io::copy(&mut stdout.unwrap(), &mut Vec::new()).await? == 0, "No output expected");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_premature_exit() -> Result<()> {
        let stream = fastx_generator(10, 10, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "head", vec!["-n", "1"]).await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        let task_result = task.await?;
        let status = child.wait().await?;
        done_rx.await??;
        assert!(status.success(), "Child process should exit successfully");
        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1"), "Output should contain first read ID");
        assert!(!output_str.contains("@read2"), "Output should not contain second read ID");
        assert!(task_result.is_ok() || task_result.is_err(), "Task may error due to broken pipe");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_resource_cleanup() -> Result<()> {
        let stream = fastx_generator(5, 10, 35.0, 3.0);
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (mut child, task) = stream_to_cmd(outputs.pop().unwrap(), "cat", vec![]).await?;
        task.await??;
        done_rx.await??;
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");
        assert!(child.stdin.is_none(), "Stdin should be closed");
        Ok(())
    }
    #[tokio::test]
    async fn test_parse_child_stdout_to_fastq() -> Result<()> {
        let mut cmd = Command::new("echo");
        cmd.arg("@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nHHHH\n");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let stdout = child.stdout.take().unwrap();
        let (tx, mut rx) = tokio::sync::mpsc::channel(100);
        tokio::spawn(parse_child_stdout_to_fastq(stdout, tx));
        let mut records = Vec::new();
        while let Some(record) = rx.recv().await {
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
        while let Some(Ok(chunk)) = stream.next().await {
            let data = chunk;
            assert!(String::from_utf8_lossy(&data).contains("test data"));
            break;
        }
        Ok(())
    }
    
    #[tokio::test]
    async fn test_stream_sequence_records_to_file() -> Result<()> {
        let _ = fs::remove_file("stream_to_file_test.fq");
        let (tx, rx) = mpsc::channel(100);
        let records = fastx_generator(2, 10, 35.0, 3.0);
        let mut records_stream = Box::pin(records);
        
        let send_task = tokio::spawn(async move {
            while let Some(record) = records_stream.next().await {
                if tx.send(record).await.is_err() {
                    return Err(anyhow::anyhow!("Failed to send record"));
                }
            }
            Ok::<(), anyhow::Error>(())
        });
        
        let write_task = tokio::spawn(async move {
            stream_sequence_records_to_file(rx, Path::new("stream_to_file_test.fq").to_path_buf()).await
        });
        
        send_task.await??;
        write_task.await??;
        
        assert!(Path::new("stream_to_file_test.fq").exists(), "Output file was not created");
        
        let content = fs::read_to_string("stream_to_file_test.fq")?;
        let num_records = content.lines().filter(|line| line.starts_with('@')).count();
        assert_eq!(num_records, 2, "Expected 2 FASTQ records, found {}", num_records);
        
        fs::remove_file("stream_to_file_test.fq")?;

        Ok(())
    }

    #[tokio::test]
    async fn test_stream_sequence_records_to_file_no_records() -> Result<()> {
        let _ = fs::remove_file("stream_to_file_norecord_test.fq");
        let (tx, rx) = mpsc::channel(100);
        let records = fastx_generator(0, 10, 35.0, 3.0);
        let mut records_stream = Box::pin(records);

        let send_task = tokio::spawn(async move {
            while let Some(record) = records_stream.next().await {
                if tx.send(record).await.is_err() {
                    return Err(anyhow::anyhow!("Failed to send record"));
                }
            }
            Ok::<(), anyhow::Error>(())
        });

        let write_task = tokio::spawn(async move {
            stream_sequence_records_to_file(rx, Path::new("stream_to_file_norecord_test.fq").to_path_buf()).await
        });

        send_task.await??;
        write_task.await??;

        assert!(Path::new("stream_to_file_norecord_test.fq").exists(), "Output file was not created");

        let content = fs::read_to_string("stream_to_file_norecord_test.fq")?;
        let num_records = content.lines().filter(|line| line.starts_with('@')).count();
        assert_eq!(num_records, 0, "Expected 2 FASTQ records, found {}", num_records);

        fs::remove_file("stream_to_file_norecord_test.fq")?;

        Ok(())
    }

    #[tokio::test]
    async fn test_read_child_to_stdout() -> Result<()> {

        let records = vec![
            SequenceRecord::Fastq {
                id: "read1".to_string(),
                desc: None,
                seq: b"ATCG".to_vec(),
                qual: b"IIII".to_vec(),
            },
            SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: b"GCTA".to_vec(),
                qual: b"HHHH".to_vec(),
            },
        ];
        let stream = tokio_stream::iter(records);
        
 
        let (mut outputs, done_rx) = t_junction(stream, 1, 50000, 10000, Some(1)).await?;
        let (child, task) = stream_to_cmd(outputs.pop().unwrap(), "cat", vec![]).await?;

        task.await??;
        if let Err(e) = read_child_stdout(child).await {
            return Err(e);
        }
        
        Ok(())
    }
    
    
}


#[cfg(stress_test)]
mod stress_tests {
    use super::*;
    use tokio::time::Duration;
    use crate::utils::fastx::fastx_generator;

    #[tokio::test]
    async fn test_t_junction_stress() -> Result<()> {
        let buffer_sizes = vec![10_000, 50_000, 100_000, 1_000_000];
        let stall_thresholds = vec![100, 1_000, 10_000, 50_000];
        let sleep_ms = vec![0, 1, 2, 5, 10];
        let stream_nums = vec![2, 5, 10];
        let num_reads = vec![10000, 100000, 1000000];
        let read_sizes  = vec![100, 1000, 10000];

        let mut log = std::fs::File::create("stress_test.log")?;
        writeln!("Buffer size\tStall\tSleep\tStreams\tReads\tSize\tTime\tMemory", buffer_size, stall, sleep, stream_num, num_read, read_size, elapsed, memory_used);
        let sys = System::new_all();

        for buffer_size in &buffer_sizes {
            for stall in &stall_thresholds {
                for sleep in &sleep_ms {
                    for stream_num in &stream_nums {
                        for num_read in &num_reads {
                            for read_size in &read_sizes {
                                let start = Instant::now();
                                let memory_before = sys.used_memory();
                                eprintln!("Buffer size: {}  Stall: {}  Sleep: {}  Streams: {} Reads: {}  Size: {}", buffer_size, stall, sleep, stream_num, num_read, read_size);
                                let stream = fastx_generator(num_read, read_size, 35.0, 3.0); 
                                let (mut outputs, done_rx) = t_junction(
                                    stream,
                                    stream_num, 
                                    *stall,
                                    if *sleep == 0 { None } else { Some(*sleep) },
                                ).await?;

                                let mut records = vec![Vec::new(), Vec::new()];
                                for i in 0..*stream_num {
                                    while let Some(Ok(record)) = outputs[i].next().await {
                                        records[i].push(record);
                                    }
                                }
                                let elapsed = start.elapsed();
                                let memory_after = sys.used_memory();
                                let memory_used = (memory_after - memory_before) / 1024 / 1024; // MB

                                for i in 0..*stream_num {
                                    assert_eq!(num_read, records[i].len());
                                }
                                done_rx.await??;

                                writeln!("{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}", buffer_size, stall, sleep, stream_num, num_read, read_size, elapsed, memory_used);
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}