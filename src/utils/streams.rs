use sysinfo::System;
use anyhow::anyhow;
use anyhow::Result;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
use std::io;
use std::os::fd::AsRawFd;


use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt, AsyncBufReadExt, BufReader, BufWriter};
use tokio::process::{Child, Command};
use tokio::sync::{mpsc, oneshot};
use tokio::time::{Duration, sleep};
use tokio_stream::{Stream, StreamExt};
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::Semaphore;
use tokio::task::JoinHandle;
use tokio::sync::Notify;
use tokio::fs::File as TokioFile;
use tokio::fs::{metadata, remove_file};
use tokio::fs::OpenOptions as TokioOpenOptions;
use std::os::unix::fs::FileTypeExt;
use futures::TryFutureExt;
use uuid::Uuid;

use crate::utils::fastx::{SequenceRecord, parse_header};
use crate::config::defs::{PipelineError, StreamDataType};
use crate::config::defs::{CoreAllocation, RunConfig};



pub trait ToBytes {
    fn to_bytes(&self) -> Result<Vec<u8>>;
}

impl ToBytes for SequenceRecord {
    fn to_bytes(&self) -> Result<Vec<u8>> {
        let mut buffer = Vec::new();  // Alloc only when needed (e.g., for stdin pipes)
        match self {
            SequenceRecord::Fastq { id, desc, seq, qual } => {
                if let Some(desc) = desc {
                    buffer.extend_from_slice(format!("@{} {}\n", id, desc).as_bytes());
                } else {
                    buffer.extend_from_slice(format!("@{}\n", id).as_bytes());
                }
                buffer.extend_from_slice(&**seq);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(&**qual);
                buffer.push(b'\n');
            }
            SequenceRecord::Fasta { id, desc, seq } => {
                if let Some(desc) = desc {
                    buffer.extend_from_slice(format!(">{} {}\n", id, desc).as_bytes());
                } else {
                    buffer.extend_from_slice(format!(">{}\n", id).as_bytes());
                }
                buffer.extend_from_slice(&**seq);
                buffer.push(b'\n');
            }
        }
        Ok(buffer)
    }
}

impl ToBytes for ParseOutput {
    fn to_bytes(&self) -> Result<Vec<u8>> {
        match self {
            ParseOutput::Fastq(record) => record.to_bytes(),
            ParseOutput::Fasta(record) => record.to_bytes(),
            ParseOutput::Bytes(bytes) => {
                // Avoid cloning if possible; return Vec only if needed
                Arc::try_unwrap(bytes.clone()).map_err(|_| anyhow!("Cannot unwrap Arc with multiple references"))
                    .or_else(|_| Ok((**bytes).clone())) // Fallback to clone if Arc is shared
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum ChildStream {
    Stdout,
    Stderr,
}

#[derive(Clone, Copy, Debug)]
pub enum ParseMode {
    Fastq,
    Fasta,
    Bytes,
    Lines
}

#[derive(Clone, Debug)]
pub enum ParseOutput {
    Fastq(SequenceRecord),
    Fasta(SequenceRecord),
    Bytes(Arc<Vec<u8>>),
}

/// Converts mpsc::Receiver<ParseOutput> type output to AsyncRead suitable to parse_fasta. parse_fastq.
pub struct ChannelReader {
    rx: mpsc::Receiver<ParseOutput>,
    buffer: Vec<u8>,
    closed: bool,
}

impl ChannelReader {
    pub fn new(rx: mpsc::Receiver<ParseOutput>) -> Self {
        Self {
            rx,
            buffer: Vec::new(),
            closed: false,
        }
    }
}

impl AsyncRead for ChannelReader {
    fn poll_read(
        mut self: std::pin::Pin<&mut Self>,
        cx: &mut std::task::Context<'_>,
        buf: &mut tokio::io::ReadBuf<'_>,
    ) -> std::task::Poll<io::Result<()>> {
        if self.closed {
            return std::task::Poll::Ready(Ok(()));
        }

        // If buffer has data, copy to output
        if !self.buffer.is_empty() {
            let to_copy = std::cmp::min(self.buffer.len(), buf.remaining());
            buf.put_slice(&self.buffer[0..to_copy]);
            self.buffer.drain(0..to_copy);
            return std::task::Poll::Ready(Ok(()));
        }

        // Fetch next chunk from channel
        match self.rx.poll_recv(cx) {
            std::task::Poll::Ready(Some(item)) => {
                match item {
                    ParseOutput::Bytes(bytes) => {
                        self.buffer = Arc::try_unwrap(bytes).unwrap_or_else(|arc| (*arc).clone());
                        let to_copy = std::cmp::min(self.buffer.len(), buf.remaining());
                        buf.put_slice(&self.buffer[0..to_copy]);
                        self.buffer.drain(0..to_copy);
                        std::task::Poll::Ready(Ok(()))
                    }
                    _ => std::task::Poll::Ready(Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Expected ParseOutput::Bytes, got another variant",
                    ))),
                }
            }
            std::task::Poll::Ready(None) => {
                self.closed = true;
                std::task::Poll::Ready(Ok(()))
            }
            std::task::Poll::Pending => std::task::Poll::Pending,
        }
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
/// - A vector of `mpsc::Receiver<T>` for downstream processing.
/// - A `oneshot::Receiver<()>` to await task completion.
pub async fn t_junction<S>(
    input: S,
    n_outputs: usize,
    base_buffer_size: usize,
    stall_threshold: u64,
    stream_sleep_ms: Option<u64>,
    backpressure_pause_ms: u64,
    data_type: StreamDataType,
    label: String,
    notify: Option<Arc<Notify>>,
) -> Result<(Vec<mpsc::Receiver<ParseOutput>>, oneshot::Receiver<Result<(), anyhow::Error>>)>
where
    S: Stream<Item = ParseOutput> + Unpin + Send + 'static,
{
    if n_outputs == 0 {
        return Err(anyhow!("No subscribers: cannot process stream in {}", label));
    }

    let mut system = System::new_all();
    system.refresh_memory();
    let mut available_ram = system.available_memory();
    if available_ram == 0 {
        available_ram = system.total_memory().max(8_000_000) / 4; // Assume 8GB min if total=0
    }

    const MAX_PROCESSES: usize = 4;
    const RAM_FRACTION: f64 = 0.5;
    const MIN_BUFFER_PER_STREAM: usize = 5_000;
    const MAX_BUFFER_PER_STREAM: usize = 1_000_000; // Cap at ~1M records to limit RAM

    let record_size = match data_type {
        StreamDataType::IlluminaFastq => 1_000, // ~1KB per FASTQ record
        StreamDataType::OntFastq => 10_000,     // ~10KB per ONT read
        StreamDataType::JustBytes => 500,       // ~500B for SAM/BAM/VCF
    };

    let min_buffer_size = MIN_BUFFER_PER_STREAM * n_outputs.max(1);
    let max_buffer_size = if available_ram > 0 {
        let calculated = ((available_ram as f64 * RAM_FRACTION / MAX_PROCESSES as f64) / record_size as f64) as usize;
        calculated.max(min_buffer_size)
    } else {
        eprintln!("Warning: Failed to detect available RAM in {}, using fallback buffer size", label);
        (base_buffer_size * n_outputs.max(1) * 2).max(min_buffer_size)
    };
    let buffer_size = (base_buffer_size * n_outputs.max(1)).clamp(min_buffer_size, max_buffer_size.min(MAX_BUFFER_PER_STREAM));
    let (done_tx, done_rx) = oneshot::channel::<Result<(), anyhow::Error>>();
    let mut output_txs: Vec<(usize, mpsc::Sender<ParseOutput>)> = Vec::with_capacity(n_outputs);
    let mut output_rxs = Vec::with_capacity(n_outputs);

    // Use bounded channels with adaptive capacity
    for i in 0..n_outputs {
        let (tx, rx) = mpsc::channel(buffer_size);
        output_txs.push((i, tx));
        output_rxs.push(rx);
    }

    tokio::spawn(async move {
        let mut input = Box::pin(input);
        let mut item_count = 0;
        let mut last_progress_time = Instant::now();
        let mut dropped_receivers = Vec::new();

        while let Some(item) = input.next().await {
            item_count += 1;
            let mut active_txs = Vec::new();
            let mut backpressure_detected = false;

            for (i, tx) in output_txs.into_iter() {
                // Check buffer fullness for backpressure
                if tx.capacity() < buffer_size / 10 { // <10% capacity left
                    backpressure_detected = true;
                    // eprintln!("{}: Backpressure detected on receiver {} at item {} (capacity {}/{}", label, i, item_count, tx.capacity(), buffer_size);
                }
                match tx.try_send(item.clone()) {
                    Ok(()) => active_txs.push((i, tx)),
                    Err(mpsc::error::TrySendError::Full(_)) => {
                        eprintln!("{}: Receiver {} full at item {}, retrying after pause", label, i, item_count);
                        sleep(Duration::from_millis(backpressure_pause_ms)).await;
                        if tx.send(item.clone()).await.is_err() {
                            eprintln!("{}: Receiver {} dropped at item {}", label, i, item_count);
                            dropped_receivers.push(i);
                        } else {
                            active_txs.push((i, tx));
                        }
                    }
                    Err(mpsc::error::TrySendError::Closed(_)) => {
                        eprintln!("{}: Receiver {} dropped at item {}", label, i, item_count);
                        dropped_receivers.push(i);
                    }
                }
            }

            output_txs = active_txs;
            if output_txs.is_empty() {
                eprintln!("{}: All receivers dropped at item {}. Receivers dropped: {:?}", label, item_count, dropped_receivers);
                let _ = done_tx.send(Err(anyhow!("{}: All receivers dropped at item {}, data loss occurred", label, item_count)));
                return;
            }

            // Adaptive throttling on backpressure
            if backpressure_detected {
                let pause_ms = backpressure_pause_ms.max(stream_sleep_ms.unwrap_or(0) * 2);
                // eprintln!("{}: Pausing for {}ms due to backpressure at item {}", label, pause_ms, item_count);
                sleep(Duration::from_millis(pause_ms)).await;
            } else if item_count % stall_threshold == 0 {
                eprintln!("{}: Processed {} items, checking for stalls", label, item_count);
                if let Some(sleep_ms) = stream_sleep_ms {
                    sleep(Duration::from_millis(sleep_ms)).await;
                }
            }

            // Periodic stall check
            if last_progress_time.elapsed() > Duration::from_secs(stall_threshold) {
                eprintln!("{}: Stall detected at item {}", label, item_count);
                last_progress_time = Instant::now();
            }
        }

        // Ensure all receivers process remaining data
        for (i, tx) in output_txs.into_iter() {
            let _ = tx; // Move Sender to drop it, signaling EOF
            // eprintln!("{}: Closed sender for receiver {} after {} items", label, i, item_count);
        }

        if let Some(n) = notify {
            n.notify_waiters();
        }

        if !dropped_receivers.is_empty() {
            let _ = done_tx.send(Err(anyhow!("{}: {} receivers dropped mid-stream: {:?}", label, dropped_receivers.len(), dropped_receivers)));
        } else {
            let _ = done_tx.send(Ok(()));
        }
    });

    Ok((output_rxs, done_rx))
}

/// Asynchronously spawn an external process and feed it a stream as stdin.
/// Capture stdout and return from function.
///
/// # Arguments
///
/// * `rx` - Receiver stream: tokio::mpsc
/// * `command` - command to shell out
/// * `args` - args for shelled out command
/// * `data_type` - Type of data being streamed (SAM/BAM, Illumina FASTQ, ONT FASTQ)
///
/// # Returns
/// tokio::process::Command containing stdout and stderr

pub async fn stream_to_cmd(
    config: Arc<RunConfig>,
    rx: mpsc::Receiver<ParseOutput>,
    cmd_tag: &str,
    args: Vec<String>,
    data_type: StreamDataType,
    verbose: bool,
) -> Result<(Arc<tokio::sync::Mutex<Child>>, JoinHandle<Result<(), anyhow::Error>>, JoinHandle<Result<(), anyhow::Error>>)> {
    let core_allocation = config.get_core_allocation(cmd_tag, None);
    let _permit = if core_allocation == CoreAllocation::Maximal {
        Some(config.maximal_semaphore.clone().acquire_owned().await?)
    } else {
        None
    };

    let (batch_size_bytes, writer_capacity) = match data_type {
        StreamDataType::JustBytes => (65_536, 65_536),
        StreamDataType::IlluminaFastq => (8_388_608, 8_388_608),
        StreamDataType::OntFastq => (1_048_576, 1_048_576),
    };

    let cmd_tag_owned = cmd_tag.to_string();
    let cmd_tag_err_owned = cmd_tag.to_string();

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
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| anyhow!("Failed to capture stderr for {}", cmd_tag_owned))?;

    let child = Arc::new(tokio::sync::Mutex::new(child));
    let child_clone = child.clone();

    let stdin_task = tokio::spawn(async move {
        let mut writer = BufWriter::with_capacity(writer_capacity, stdin);
        let mut batch = Vec::with_capacity(batch_size_bytes);
        let mut total_written = 0;

        let mut stream = ReceiverStream::new(rx);
        while let Some(item) = stream.next().await {
            match &item {
                ParseOutput::Bytes(arc_bytes) => batch.extend_from_slice(&**arc_bytes),
                ParseOutput::Fastq(record) => {  // Direct match on inner enum
                    if let SequenceRecord::Fastq { id, desc, seq, qual } = record {
                        if let Some(d) = desc {
                            batch.extend_from_slice(format!("@{} {}\n", id, d).as_bytes());
                        } else {
                            batch.extend_from_slice(format!("@{}\n", id).as_bytes());
                        }
                        batch.extend_from_slice(&**seq);
                        batch.extend_from_slice(b"\n+\n");
                        batch.extend_from_slice(&**qual);
                        batch.push(b'\n');
                    }
                }
                ParseOutput::Fasta(record) => {
                    if let SequenceRecord::Fasta { id, desc, seq } = record {
                        if let Some(d) = desc {
                            batch.extend_from_slice(format!(">{} {}\n", id, d).as_bytes());
                        } else {
                            batch.extend_from_slice(format!(">{}\n", id).as_bytes());
                        }
                        batch.extend_from_slice(&**seq);
                        batch.push(b'\n');
                    }
                }
            }

            if batch.len() >= batch_size_bytes {
                if let Err(e) = writer.write_all(&batch).await {
                    if e.kind() == std::io::ErrorKind::BrokenPipe {
                        // Check status before ignoring (non-blocking)
                        let mut guard = child_clone.lock().await;
                        if let Ok(Some(status)) = guard.try_wait() {
                            if status.success() {
                                eprintln!("Ignoring BrokenPipe in {} after successful exit (code {:?})", cmd_tag_owned, status.code());
                                break; // Safe: No loss, child done
                            } else {
                                return Err(anyhow!("BrokenPipe in {} with failed child exit: {:?}", cmd_tag_owned, status));
                            }
                        } else {
                            // Child not exited yet; treat as error (potential loss)
                            return Err(anyhow!("BrokenPipe in {} before child exit: {}", cmd_tag_owned, e));
                        }
                    } else {
                        return Err(anyhow!("Write error in {}: {}", cmd_tag_owned, e));
                    }
                }
                writer.flush().await?;
                total_written += batch.len();
                batch.clear();
            }
        }

        if !batch.is_empty() {
            if let Err(e) = writer.write_all(&batch).await {
                if e.kind() == std::io::ErrorKind::BrokenPipe {
                    let mut guard = child_clone.lock().await;
                    if let Ok(Some(status)) = guard.try_wait() {
                        if status.success() {
                            eprintln!("Ignoring final BrokenPipe in {} after successful exit (code {:?})", cmd_tag_owned, status.code());
                        } else {
                            return Err(anyhow!("Final BrokenPipe in {} with failed child exit: {:?}", cmd_tag_owned, status));
                        }
                    } else {
                        return Err(anyhow!("Final BrokenPipe in {} before child exit: {}", cmd_tag_owned, e));
                    }
                } else {
                    return Err(anyhow!("Final write error in {}: {}", cmd_tag_owned, e));
                }
            } else {
                writer.flush().await?;
                total_written += batch.len();
            }
        }

        writer.flush().await?;
        drop(writer); // Close stdin, signal EOF

        // Final check: Await child completion and verify success
        let mut guard = child_clone.lock().await;
        let status = guard.wait().await?;
        if !status.success() {
            return Err(anyhow!("Child {} failed after stream: exit {:?}", cmd_tag_owned, status));
        }
        Ok(())
    });

    let stderr_task = tokio::spawn(async move {
        let mut reader = BufReader::with_capacity(1_048_576, stderr);
        let mut buffer = vec![0u8; 8192];
        if verbose {
            loop {
                match reader.read(&mut buffer).await {
                    Ok(0) => break,
                    Ok(n) => {
                        let stderr_chunk = String::from_utf8_lossy(&buffer[..n]);
                        eprintln!("[{} stderr]: {}", cmd_tag_err_owned, stderr_chunk);
                    }
                    Err(e) => {
                        eprintln!("Error reading stderr for {}: {}", cmd_tag_err_owned, e);
                        return Err(anyhow!("Failed to read stderr: {}", e));
                    }
                }
            }
        } else {
            while reader.read(&mut buffer).await? != 0 {}
        }
        Ok(())
    });

    Ok((child, stdin_task, stderr_task))
}

pub async fn spawn_cmd(
    config: Arc<RunConfig>,
    cmd_tag: &str,
    args: Vec<String>,
    verbose: bool,
) -> Result<(Child, JoinHandle<Result<(), anyhow::Error>>)> {
    let core_allocation = config.get_core_allocation(cmd_tag, None);
    let _permit = if core_allocation == CoreAllocation::Maximal {
        Some(config.maximal_semaphore.clone().acquire_owned().await?)
    } else {
        None
    };

    let cmd_tag_owned = cmd_tag.to_string();
    let mut child = Command::new(&cmd_tag_owned)
        .args(&args)
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow!("Failed to spawn {}: {}", cmd_tag_owned, e))?;

    let stderr_task = {
        let stderr = child
            .stderr
            .take()
            .ok_or_else(|| anyhow!("Failed to capture stderr for {}", cmd_tag_owned))?;
        let cmd_tag_clone = cmd_tag_owned.clone();
        tokio::spawn(async move {
            let mut reader = BufReader::with_capacity(1024 * 1024, stderr);
            let mut buffer = vec![0u8; 8192];
            if verbose {
                loop {
                    match reader.read(&mut buffer).await {
                        Ok(0) => break,
                        Ok(n) => {
                            let stderr_chunk = String::from_utf8_lossy(&buffer[..n]);
                            eprintln!("[{} stderr]: {}", cmd_tag_clone, stderr_chunk);
                        }
                        Err(e) => {
                            eprintln!("Error reading stderr for {}: {}", cmd_tag_clone, e);
                            return Err(anyhow!("Failed to read stderr: {}", e));
                        }
                    }
                }
            } else {
                while reader.read(&mut buffer).await? != 0 {}
            }
            Ok(())
        })
    };

    // Permit is automatically dropped when the function exits, releasing it
    Ok((child, stderr_task))
}

/// Parse the output of a stream, either Fastq or simple bytes.
///
/// # Arguments
///
/// * `child` - tokio::process::Child
/// * `stream` - stream from child
/// * `mode` - ParseMode enum
/// * `buffer_size` - stream buffer size to pass to helpers
///
/// # Returns
/// Result<mpsc::Receiver<ParseOutput>>
pub async fn parse_child_output(
    child: &mut Child,
    stream: ChildStream,
    mode: ParseMode,
    buffer_size: usize,
) -> Result<mpsc::Receiver<ParseOutput>> {
    match (stream, mode) {
        (ChildStream::Stdout, ParseMode::Fastq) => {
            let stdout = child.stdout.take().ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_fastq(stdout, buffer_size).await
        }
        (ChildStream::Stderr, ParseMode::Fastq) => {
            let stderr = child.stderr.take().ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_fastq(stderr, buffer_size).await
        }
        (ChildStream::Stdout, ParseMode::Fasta) => {
            let stdout = child.stdout.take().ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_fasta(stdout, buffer_size).await
        }
        (ChildStream::Stderr, ParseMode::Fasta) => {
            let stderr = child.stderr.take().ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_fasta(stderr, buffer_size).await
        }
        (ChildStream::Stdout, ParseMode::Bytes) => {
            let stdout = child.stdout.take().ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_bytes(stdout, buffer_size).await
        }
        (ChildStream::Stderr, ParseMode::Bytes) => {
            let stderr = child.stderr.take().ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_bytes(stderr, buffer_size).await
        }
        (ChildStream::Stdout, ParseMode::Lines) => {
            let stdout = child.stdout.take().ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_lines(stdout, buffer_size).await
        }
        (ChildStream::Stderr, ParseMode::Lines) => {
            let stderr = child.stderr.take().ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_lines(stderr, buffer_size).await
        }
    }
}

/// Helper function to parse FASTQ data
///
/// # Arguments
///
/// * `reader` - reading stream
/// * `buffer_size` - stream buffer size
///
/// # Returns
/// Result<mpsc::Receiver<ParseOutput>>
pub async fn parse_fastq<R: AsyncRead + Unpin + Send + 'static>(
    reader: R,
    buffer_size: usize,
) -> Result<mpsc::Receiver<ParseOutput>> {
    let (tx, rx) = mpsc::channel(buffer_size);
    let mut reader = BufReader::with_capacity(1024 * 1024, reader);
    let mut buffer = String::new();
    let mut count = 0;

    tokio::spawn(async move {
        loop {
            buffer.clear();
            let bytes_read = match reader.read_line(&mut buffer).await {
                Ok(n) => n,
                Err(e) => {
                    eprintln!("parse_fastq: Error reading FASTQ line: {}", e);
                    return;
                }
            };
            if bytes_read == 0 {
                eprintln!("parse_fastq: Reached end of input, processed {} records", count);
                break;
            }
            let id_line = buffer.trim_end();
            if !id_line.starts_with('@') {
                eprintln!("parse_fastq: Invalid FASTQ format: expected '@', got '{}'", id_line);
                return;
            }

            let id = id_line[1..].to_string();
            let desc = None;

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                eprintln!("parse_fastq: Error reading sequence line");
                return;
            }
            let seq = buffer.trim_end().as_bytes().to_vec();
            if seq.is_empty() {
                eprintln!("parse_fastq: Missing sequence");
                return;
            }

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                eprintln!("parse_fastq: Error reading plus line");
                return;
            }
            let plus = buffer.trim_end();
            if plus != "+" {
                eprintln!("parse_fastq: Invalid FASTQ format: expected '+', got '{}'", plus);
                return;
            }

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                eprintln!("parse_fastq: Error reading quality line");
                return;
            }
            let qual = buffer.trim_end().as_bytes().to_vec();
            if qual.len() != seq.len() {
                eprintln!("parse_fastq: Sequence and quality lengths do not match: seq_len={}, qual_len={}", seq.len(), qual.len());
                return;
            }

            let record = SequenceRecord::Fastq {
                id,
                desc,
                seq: Arc::new(seq),
                qual: Arc::new(qual),
            };

            if tx.send(ParseOutput::Fastq(record)).await.is_err() {
                eprintln!("parse_fastq: Receiver dropped after {} records", count);
                break;
            }
            count += 1;
        }
        eprintln!("parse_fastq: Completed parsing, sent {} FASTQ records", count);
    });

    Ok(rx)
}


/// Helper function to parse FASTA data
///
/// # Arguments
///
/// * `reader` - reading stream
/// * `buffer_size` - stream buffer size
///
/// # Returns
/// Result<mpsc::Receiver<ParseOutput>>
pub async fn parse_fasta<R: AsyncRead + Unpin + Send + 'static>(
    reader: R,
    buffer_size: usize,
) -> Result<mpsc::Receiver<ParseOutput>> {
    let (tx, rx) = mpsc::channel(buffer_size);
    let mut reader = BufReader::with_capacity(1024 * 1024, reader);
    let mut buffer = String::new();
    let mut current_id = None;
    let mut current_desc = None;
    let mut current_seq = Vec::new();

    tokio::spawn(async move {
        while reader.read_line(&mut buffer).await.is_ok() {
            let line = buffer.trim_end();
            if line.is_empty() {
                break;
            }
            if line.starts_with('>') {
                if let Some(id) = current_id.take() {
                    let record = SequenceRecord::Fasta {
                        id,
                        desc: current_desc.take(),
                        seq: Arc::new(current_seq),
                    };
                    if tx.send(ParseOutput::Fasta(record)).await.is_err() {
                        return Err(anyhow!("Receiver dropped during FASTA parsing, data loss detected"));
                    }
                    current_seq = Vec::new();
                }
                let (id, desc) = parse_header(line.as_bytes(), '>');
                current_id = Some(id);
                current_desc = desc;
            } else if current_id.is_some() {
                current_seq.extend_from_slice(line.as_bytes());
            }
            buffer.clear();
        }
        if let Some(id) = current_id.take() {
            let record = SequenceRecord::Fasta {
                id,
                desc: current_desc.take(),
                seq: Arc::new(current_seq), // Move final current_seq
            };
            if tx.send(ParseOutput::Fasta(record)).await.is_err() {
                return Err(anyhow!("Receiver dropped during final FASTA parsing, data loss detected"));
            }
        }
        Ok::<(), anyhow::Error>(())
    });
    Ok(rx)
}

/// Helper function to parse byte data
///
/// # Arguments
///
/// * `reader` - reading stream
/// * `buffer_size` - stream buffer size
///
/// # Returns
/// Result<mpsc::Receiver<ParseOutput>>
pub async fn parse_bytes<R: AsyncRead + Unpin + Send + 'static>(
    reader: R,
    buffer_size: usize,
) -> Result<mpsc::Receiver<ParseOutput>> {
    let (tx, rx) = mpsc::channel(buffer_size);

    tokio::spawn(async move {
        let mut reader = BufReader::with_capacity(1024 * 1024, reader);
        let mut buffer = vec![0u8; 8192];
        loop {
            let bytes_read = match reader.read(&mut buffer).await {
                Ok(n) => n,
                Err(e) => {
                    eprintln!("Error reading bytes: {}", e);
                    return;
                }
            };
            if bytes_read == 0 {
                break;
            }
            if tx.send(ParseOutput::Bytes(Arc::new(buffer[..bytes_read].to_vec()))).await.is_err() {
                eprintln!("Receiver dropped while sending bytes");
                break;
            }
        }
    });

    Ok(rx)
}

/// Helper function to parse line data
///
/// # Arguments
///
/// * `reader` - reading stream
/// * `buffer_size` - stream buffer size
///
/// # Returns
/// Result<mpsc::Receiver<ParseOutput>>
pub async fn parse_lines<R: AsyncRead + Unpin + Send + 'static>(
    reader: R,
    buffer_size: usize,
) -> Result<mpsc::Receiver<ParseOutput>> {
    let (tx, rx) = mpsc::channel(buffer_size);
    let mut reader = BufReader::with_capacity(1024 * 1024, reader);

    tokio::spawn(async move {
        let mut line = String::new();
        loop {
            line.clear();
            let bytes_read = match reader.read_line(&mut line).await {
                Ok(n) => n,
                Err(e) => {
                    eprintln!("Error reading line: {}", e);
                    return;
                }
            };
            if bytes_read == 0 {
                break;
            }
            let trimmed = line.trim_end().to_string();
            if !trimmed.is_empty() {
                if tx.send(ParseOutput::Bytes(Arc::new(trimmed.into_bytes()))).await.is_err() {
                    eprintln!("Receiver dropped while sending line");
                    break;
                }
            }
        }
    });

    Ok(rx)
}

/// Takes a ReceiverStream and writes it to a file
///
/// # Arguments
///
/// * `rx` - ReceiverStream<ParseOutput>, parsed by the parsing functions above.
/// * `path` - PathBuf to file.
///
/// # Returns
/// Result<()>
pub async fn stream_to_file(
    rx: mpsc::Receiver<ParseOutput>,
    path: PathBuf,
) -> Result<(), anyhow::Error> {
    let mut file = BufWriter::with_capacity(4 * 1024 * 1024, tokio::fs::File::create(path).await?);
    let mut stream = ReceiverStream::new(rx);

    while let Some(item) = stream.next().await {
        match &item {
            ParseOutput::Fastq(record) => {
                let bytes = record.to_bytes()?;
                file.write_all(&bytes).await?;
            }
            ParseOutput::Fasta(record) => {
                let bytes = record.to_bytes()?;
                file.write_all(&bytes).await?;
            }
            ParseOutput::Bytes(arc_bytes) => {
                file.write_all(&*arc_bytes).await?;
            }
        }
    }

    file.flush().await?;
    Ok(())
}

/// Reads a child process's stdout or stderr into a Vec<String>.
///
/// # Arguments
///
/// * `child` - Child process with stdout/stderr.
/// * `stream` - Selects stdout or stderr (default: stdout).
///
/// # Returns
/// Result<Vec<String>> containing the lines from the selected stream.
pub async fn read_child_output_to_vec(child: &mut Child, stream: ChildStream) -> Result<Vec<String>> {
    let rx = parse_child_output(child, stream, ParseMode::Bytes, 1000).await?;
    let mut stream = ReceiverStream::new(rx);
    let mut lines = Vec::new();

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(chunk) => {
                let text = String::from_utf8_lossy(&chunk);
                for line in text.lines() {
                    let trimmed = line.trim_end();
                    if !trimmed.is_empty() {
                        lines.push(trimmed.to_string());
                    }
                }
            }
            ParseOutput::Fastq(_) => {
                return Err(anyhow!("Unexpected Fastq record when parsing bytes"));
            }
            ParseOutput::Fasta(_) => {
                return Err(anyhow!("Unexpected Fastq record when parsing bytes"));
            }
        }
    }

    let status = child.wait().await?;
    if !status.success() {
        return Err(anyhow!("Child process exited with non-zero status: {}", status));
    }

    Ok(lines)
}


/// Creates an async FIFO pipe using Tokio command.
///
/// # Arguments
///
/// * `path` - PathBuf to to the new FIFO
///
/// # Returns
/// Result(): whether the fifo was created successfully.
pub async fn create_fifo(path: &PathBuf) -> Result<()> {
    let status = Command::new("mkfifo")
        .arg(path)
        .status()
        .await
        .map_err(|e| anyhow!("Failed to execute mkfifo: {}", e))?;
    if !status.success() {
        return Err(anyhow!("mkfifo failed with exit code: {:?}", status.code()));
    }
    Ok(())
}

/// Writes a byte stream to a FIFO
///
/// # Arguments
///
/// * `rx` - Steram of byte data
/// * `fifo_path` - PathBuf to to the new FIFO
///
/// # Returns
/// Result(): whether the fifo wrote successfully
pub async fn write_to_fifo(mut rx: mpsc::Receiver<ParseOutput>, fifo_path: PathBuf) -> Result<()> {
    // Open for write-only, without create/truncate to avoid issues with existing FIFOs
    let file = TokioOpenOptions::new()
        .write(true)
        .open(&fifo_path)
        .await
        .map_err(|e| anyhow!("Failed to open FIFO {} for write: {}", fifo_path.display(), e))?;

    // Use 16MB buffer for high-throughput runs on 1.5TB RAM nodes
    let mut writer = BufWriter::with_capacity(16_777_216, file);
    eprintln!("Opened FIFO for writing: {}", fifo_path.display());

    let mut bytes_written = 0;
    while let Some(item) = rx.recv().await {
        let bytes = item.to_bytes()
            .map_err(|e| anyhow!("Failed to convert item to bytes: {}", e))?;
        writer.write_all(&bytes).await
            .map_err(|e| {
                if e.kind() == std::io::ErrorKind::BrokenPipe {
                    anyhow!("Broken pipe writing to FIFO {}: possible early reader exit", fifo_path.display())
                } else {
                    anyhow!("Write to FIFO {} failed: {}", fifo_path.display(), e)
                }
            })?;
        bytes_written += bytes.len();

        // Log every 10MB for monitoring high-throughput runs
        if bytes_written % 10_485_760 == 0 {
            eprintln!("Wrote {} bytes to FIFO: {}", bytes_written, fifo_path.display());
        }
    }

    writer.flush().await
        .map_err(|e| anyhow!("Flush to FIFO {} failed: {}", fifo_path.display(), e))?;
    eprintln!("Finished writing {} bytes to FIFO: {}", bytes_written, fifo_path.display());
    Ok(())
}

pub async fn deinterleave_fastq_stream(
    input_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    buffer_size: usize,
) -> Result<(mpsc::Receiver<ParseOutput>, Option<mpsc::Receiver<ParseOutput>>, JoinHandle<Result<(), anyhow::Error>>)> {
    let (r1_tx, r1_rx) = mpsc::channel::<ParseOutput>(buffer_size);
    let r2_tx_rx_opt = if paired {
        Some(mpsc::channel::<ParseOutput>(buffer_size))
    } else {
        None
    };
    let r2_tx_opt = r2_tx_rx_opt.as_ref().map(|(tx, _rx)| tx.clone());
    let r2_rx_opt = r2_tx_rx_opt.map(|(_tx, rx)| rx);

    let deinterleave_handle = tokio::spawn(async move {
        let mut stream = input_stream;
        let mut is_r1 = true; // Track expected read pair
        let mut record_count = 0;
        while let Some(item) = stream.next().await {
            record_count += 1;
            if let ParseOutput::Fastq(record) = item {
                let bytes = record.to_bytes().map_err(|e| anyhow!("Failed to convert FASTQ to bytes: {}", e))?;
                let bytes_output = ParseOutput::Bytes(Arc::new(bytes));
                if paired {
                    if is_r1 {
                        r1_tx.send(bytes_output.clone()).await
                            .map_err(|_| anyhow!("R1 send failed"))?;
                    } else {
                        r2_tx_opt.as_ref().unwrap().send(bytes_output).await
                            .map_err(|_| anyhow!("R2 send failed"))?;
                    }
                    is_r1 = !is_r1; // Toggle for next read
                } else {
                    r1_tx.send(bytes_output).await
                        .map_err(|_| anyhow!("Single-end send failed"))?;
                }
            } else {
                return Err(anyhow!("Non-FASTQ item in deinterleave stream"));
            }
        }

        if paired && !is_r1 {
            return Err(anyhow!("Incomplete paired-end stream: {} records, expected even number", record_count));
        }
        eprintln!("deinterleave_fastq_stream: Completed deinterleaving {} records", record_count);
        Ok(())
    });

    Ok((r1_rx, r2_rx_opt, deinterleave_handle))
}

/// Deinterleaves an interleaved FASTQ stream into two separate streams, writing each to a FIFO in RAM temp dir.
/// assumes input is interleaved if paired. for single-end, all goes to R1 FIFO.
///
/// # Arguments
/// * `config` - RunConfig struct
/// * `input_stream` - Interleaved FASTQ stream
/// * `sample_base` - basename for FIFOs
/// * `paired` - If false, routes all to R1 FIFO (R2 FIFO created but unused).
///
/// # Returns
/// Tuple:
/// (r1_fifo_path, r2_fifo_path, deinterleave_handle, r1_write_handle, r2_write_handle)
pub async fn deinterleave_fastq_stream_to_fifos(
    config: Arc<RunConfig>,
    mut input_stream: ReceiverStream<ParseOutput>,
    sample_base: &str,
    paired: bool,
) -> Result<(PathBuf, PathBuf, JoinHandle<Result<(), anyhow::Error>>, JoinHandle<Result<(), anyhow::Error>>, Option<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let ram_temp_dir = config.ram_temp_dir.clone();
    let r1_fifo = ram_temp_dir.join(format!("{}_R1.fq", sample_base));
    let r2_fifo = ram_temp_dir.join(format!("{}_R2.fq", sample_base));

    // Proactively remove existing FIFOs with retry
    eprintln!("Cleaning up existing FIFOs: R1={}, R2={}", r1_fifo.display(), r2_fifo.display());
    for fifo in [&r1_fifo, &r2_fifo] {
        for attempt in 1..=3 {
            if metadata(fifo).await.is_ok() {
                eprintln!("Attempt {} to remove existing FIFO: {}", attempt, fifo.display());
                match remove_file(fifo).await {
                    Ok(_) => eprintln!("Successfully removed FIFO: {}", fifo.display()),
                    Err(e) => {
                        eprintln!("Failed to remove FIFO {} on attempt {}: {}", fifo.display(), attempt, e);
                        if attempt < 3 {
                            sleep(Duration::from_millis(100)).await; // Brief delay before retry
                            continue;
                        }
                        return Err(PipelineError::IOError(format!("Failed to remove existing FIFO {} after {} attempts: {}", fifo.display(), attempt, e)));
                    }
                }
            } else {
                eprintln!("No existing FIFO found: {}", fifo.display());
                break;
            }
        }
    }

    // Create FIFOs
    eprintln!("Creating FIFOs: R1={}, R2={}", r1_fifo.display(), r2_fifo.display());
    create_fifo(&r1_fifo).await
        .map_err(|e| PipelineError::IOError(format!("Failed to create R1 FIFO {}: {}", r1_fifo.display(), e)))?;
    if !metadata(&r1_fifo).await.map_err(|e| PipelineError::IOError(format!("R1 FIFO check failed: {}", e)))?.file_type().is_fifo() {
        return Err(PipelineError::IOError(format!("R1 path {} is not a FIFO", r1_fifo.display())));
    }

    if paired {
        create_fifo(&r2_fifo).await
            .map_err(|e| PipelineError::IOError(format!("Failed to create R2 FIFO {}: {}", r2_fifo.display(), e)))?;
        if !metadata(&r2_fifo).await.map_err(|e| PipelineError::IOError(format!("R2 FIFO check failed: {}", e)))?.file_type().is_fifo() {
            return Err(PipelineError::IOError(format!("R2 path {} is not a FIFO", r2_fifo.display())));
        }
    }


    let buffer_size = (config.base_buffer_size * 10).min(100_000);
    let (r1_tx, r1_rx) = mpsc::channel::<ParseOutput>(buffer_size);
    let (r2_tx, r2_rx) = mpsc::channel::<ParseOutput>(buffer_size);

    // Deinterleave task
    let deinterleave_handle = tokio::spawn(async move {
        let mut is_r1 = true;
        let mut record_count = 0;
        while let Some(item) = input_stream.next().await {
            record_count += 1;
            match item {
                ParseOutput::Fastq(record) => {
                    if paired {
                        if is_r1 {
                            r1_tx.send(ParseOutput::Fastq(record)).await
                                .map_err(|_| PipelineError::StreamDataDropped)
                                .map_err(Into::<anyhow::Error>::into)?;
                        } else {
                            r2_tx.send(ParseOutput::Fastq(record)).await
                                .map_err(|_| PipelineError::StreamDataDropped)
                                .map_err(Into::<anyhow::Error>::into)?;
                        }
                        is_r1 = !is_r1;
                    } else {
                        r1_tx.send(ParseOutput::Fastq(record)).await
                            .map_err(|_| PipelineError::StreamDataDropped)
                            .map_err(Into::<anyhow::Error>::into)?;
                    }
                }
                _ => return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in deinterleave stream".to_string()))
                    .map_err(Into::<anyhow::Error>::into),
            }
        }
        if paired && !is_r1 {
            return Err(PipelineError::InvalidFastqFormat(format!(
                "Incomplete paired-end stream: {} records, expected even number", record_count
            )))
                .map_err(Into::<anyhow::Error>::into);
        }
        eprintln!("deinterleave_fastq_stream_to_fifos: Completed deinterleaving {} records", record_count);
        Ok(())
    });

    // Writer tasks
    let r1_write_handle = tokio::spawn(write_to_fifo(r1_rx, r1_fifo.clone()));
    let r2_write_handle = if paired {
        Some(tokio::spawn(write_to_fifo(r2_rx, r2_fifo.clone())))
    }
    else {
        None
    };

    // Defer cleanup to pipeline completion (handled in run cleanup_tasks)
    Ok((r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle))
}




pub async fn deinterleave_fastq_stream_to_process_sub<S>(
    config: Arc<RunConfig>,
    input_stream: S,
    sample_id: &str,
    paired: bool,
) -> Result<(
    String,                                // r1_fd_str
    Option<String>,                       // r2_fd_str
    JoinHandle<Result<(), anyhow::Error>>, // combined_handle
    Option<JoinHandle<Result<(), anyhow::Error>>>, // unused (kept for compatibility)
    Arc<Notify>,                          // writers_ready
    PathBuf,                              // r1_fifo
    Option<PathBuf>,                      // r2_fifo
)>
where
    S: Stream<Item = ParseOutput> + Unpin + Send + 'static,
{
    // Split stream into two (R1 and R2 for paired, R1 only for single-end)
    let num_outputs = if paired { 2 } else { 1 }; // Fix: Only 1 or 2 outputs
    let (output_rxs, done_rx) = t_junction(
        input_stream,
        num_outputs,
        config.base_buffer_size * 10, // Large buffer for 1.5TB RAM
        config.args.stall_threshold,
        Some(10), // Throttle for large inputs
        100,      // Backpressure pause
        StreamDataType::IlluminaFastq,
        format!("deinterleave_process_sub_{}", sample_id),
        None,
    ).await?;

    let mut rxs = output_rxs.into_iter();
    let r1_write_rx = rxs.next().ok_or(PipelineError::EmptyStream)?;
    let r2_write_rx = if paired { rxs.next() } else { None };

    let writers_ready = Arc::new(Notify::new());

    // Create temp FIFO for R1 (use RAM-based temp dir for speed)
    let r1_fifo = config.ram_temp_dir.join(format!("{}_r1_{}.fifo", sample_id, Uuid::new_v4()));
    create_fifo(&r1_fifo).await.map_err(|e| anyhow!("Failed to create R1 FIFO: {}", e))?;

    // Spawn cat process for R1, reading from FIFO
    let mut r1_child = Command::new("cat")
        .arg(&r1_fifo)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
    let r1_child_stdout = r1_child.stdout.take().ok_or(anyhow!("No stdout for R1 cat"))?;
    let r1_fd = r1_child_stdout.as_raw_fd();
    let r1_fd_str = format!("/dev/fd/{}", r1_fd);
    eprintln!("Created R1 file descriptor: {}", r1_fd_str);

    // Hold stdout to prevent premature closure
    let r1_child_stdout = Arc::new(r1_child_stdout); // Keep alive
    let r1_cat_handle: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
        r1_child.wait().await.map_err(|e| anyhow!("R1 cat process failed: {}", e))?;
        eprintln!("R1 cat process completed");
        Ok(())
    });

    // Spawn R1 writer
    let r1_fifo_clone = r1_fifo.clone(); // Clone for use in closure
    let writers_ready_clone = writers_ready.clone();
    let r1_child_stdout_clone = r1_child_stdout.clone(); // Keep stdout alive
    let r1_handle = tokio::spawn(async move {
        eprintln!("Opened R1 process substitution for writing");
        writers_ready_clone.notify_one(); // Signal ready
        let mut expect_r1 = true;
        let mut input = ReceiverStream::new(r1_write_rx);
        let mut r1_writer = BufWriter::with_capacity(16_777_216, TokioFile::create(&r1_fifo_clone).await?);
        let mut bytes_written = 0;

        while let Some(item) = input.next().await {
            match item {
                ParseOutput::Fastq(_) => {
                    if paired {
                        if expect_r1 {
                            let bytes = item.to_bytes()?;
                            r1_writer.write_all(&bytes).await?;
                            bytes_written += bytes.len();
                            expect_r1 = false;
                        } else {
                            // Skip R2 records for R1 stream
                            expect_r1 = true;
                            continue;
                        }
                    } else {
                        let bytes = item.to_bytes()?;
                        r1_writer.write_all(&bytes).await?;
                        bytes_written += bytes.len();
                    }
                }
                _ => return Err(anyhow!("Non-FASTQ in R1 stream: {:?}", item)),
            }
        }
        r1_writer.flush().await?;
        if paired && !expect_r1 {
            return Err(anyhow!("Incomplete paired-end FASTQ: missing R2"));
        }
        eprintln!("Finished writing {} bytes to R1 process substitution", bytes_written);
        drop(r1_child_stdout_clone); // Explicitly drop stdout after writing
        done_rx.await??; // Check for stream drops
        Ok(())
    });

    // R2 handling for paired-end
    let (r2_fd_str, r2_handle, r2_cat_handle, r2_fifo) = if paired {
        let r2_fifo = config.ram_temp_dir.join(format!("{}_r2_{}.fifo", sample_id, Uuid::new_v4()));
        create_fifo(&r2_fifo).await.map_err(|e| anyhow!("Failed to create R2 FIFO: {}", e))?;

        let r2_fifo_clone = r2_fifo.clone(); // Clone for use in closure
        let mut r2_child = Command::new("cat")
            .arg(&r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;
        let r2_child_stdout = r2_child.stdout.take().ok_or(anyhow!("No stdout for R2 cat"))?;
        let r2_fd = r2_child_stdout.as_raw_fd();
        let r2_fd_str = format!("/dev/fd/{}", r2_fd);
        eprintln!("Created R2 file descriptor: {}", r2_fd_str);

        let r2_child_stdout = Arc::new(r2_child_stdout); // Keep alive
        let r2_cat_handle: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
            r2_child.wait().await.map_err(|e| anyhow!("R2 cat process failed: {}", e))?;
            eprintln!("R2 cat process completed");
            Ok(())
        });

        let writers_ready_clone = writers_ready.clone();
        let r2_child_stdout_clone = r2_child_stdout.clone(); // Keep stdout alive
        let r2_handle = tokio::spawn(async move {
            eprintln!("Opened R2 process substitution for writing");
            writers_ready_clone.notify_one(); // Signal ready
            let mut input = ReceiverStream::new(r2_write_rx.unwrap());
            let mut r2_writer = BufWriter::with_capacity(16_777_216, TokioFile::create(&r2_fifo_clone).await?);
            let mut bytes_written = 0;
            let mut expect_r2 = false;

            while let Some(item) = input.next().await {
                match item {
                    ParseOutput::Fastq(_) => {
                        if expect_r2 {
                            let bytes = item.to_bytes()?;
                            r2_writer.write_all(&bytes).await?;
                            bytes_written += bytes.len();
                            expect_r2 = false;
                        } else {
                            // Skip R1 records for R2 stream
                            expect_r2 = true;
                            continue;
                        }
                    }
                    _ => return Err(anyhow!("Non-FASTQ in R2 stream: {:?}", item)),
                }
            }
            r2_writer.flush().await?;
            if expect_r2 {
                return Err(anyhow!("Incomplete paired-end FASTQ: missing R1"));
            }
            eprintln!("Finished writing {} bytes to R2 process substitution", bytes_written);
            drop(r2_child_stdout_clone); // Explicitly drop stdout after writing
            Ok(())
        });

        (Some(r2_fd_str), Some(r2_handle), Some(r2_cat_handle), Some(r2_fifo))
    } else {
        (None, None, None, None)
    };

    // Combine handles for R1 and R2 (if paired)
    let combined_handle = tokio::spawn(async move {
        r1_handle.await??;
        if let Some(r2) = r2_handle {
            r2.await??;
        }
        r1_cat_handle.await??;
        if let Some(r2_cat) = r2_cat_handle {
            r2_cat.await??;
        }
        Ok(())
    });

    Ok((r1_fd_str, r2_fd_str, combined_handle, None, writers_ready, r1_fifo, r2_fifo))
}


/// Combines multiple streams into one
///
/// # Arguments
///
/// * `streams` - Vec of Receiver<ParseOutput>
/// * `buffer_size` - Size of buffer in bytes.
///
/// # Returns
/// Receiver stream and result
pub async fn y_junction(
    streams: Vec<mpsc::Receiver<ParseOutput>>,
    buffer_size: usize,
) -> Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>)> {

    let (tx, rx) = mpsc::channel(buffer_size);

    let combined_task = tokio::spawn(async move {
        for mut rx in streams {
            while let Some(item) = rx.recv().await {
                if tx.send(item).await.is_err() {
                    return Err(anyhow!("Failed to send item to combined stream"));
                }
            }
        }
        Ok::<(), anyhow::Error>(())
    });

    Ok((rx, combined_task))
}


pub async fn convert_fasta_stream_to_sequence_record(
    rx: mpsc::Receiver<ParseOutput>,
    buffer_size: usize,
) -> Result<(mpsc::Receiver<SequenceRecord>, JoinHandle<Result<(), anyhow::Error>>)> {
    let (stats_tx, stats_rx) = mpsc::channel(buffer_size);
    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(rx);
        while let Some(item) = stream.next().await {
            if let ParseOutput::Fasta(record) = item {
                match record {
                    SequenceRecord::Fasta { .. } => {
                        if stats_tx.send(record).await.is_err() {
                            eprintln!("Failed to send FASTA record to stats");
                            return Err(anyhow::anyhow!("Failed to send FASTA record"));
                        }
                    }
                    SequenceRecord::Fastq { .. } => {
                        eprintln!("Unexpected FASTQ record in consensus_stats_stream");
                    }
                }
            } else {
                eprintln!("Unexpected non-FASTA item in stream");
            }
        }
        Ok(())
    });
    Ok((stats_rx, task))
}


pub async fn bytes_to_lines(
    rx: mpsc::Receiver<ParseOutput>,
    buffer_size: usize,
) -> Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>)> {
    let (tx, output_rx) = mpsc::channel(buffer_size);
    let mut leftover = Vec::new();

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(rx);
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(arc_bytes) => {
                    let mut bytes = arc_bytes.as_ref().to_vec();
                    if !leftover.is_empty() {
                        leftover.append(&mut bytes);
                        bytes = leftover;
                        leftover = Vec::new();
                    }

                    let mut start = 0;
                    for i in 0..bytes.len() {
                        if bytes[i] == b'\n' {
                            let line = Arc::new(bytes[start..=i].to_vec());
                            if tx.send(ParseOutput::Bytes(line)).await.is_err() {
                                eprintln!("Warning: Receiver dropped in bytes_to_lines, stopping line processing");
                                break;
                            }
                            start = i + 1;
                        }
                    }

                    if start < bytes.len() {
                        leftover = bytes[start..].to_vec();
                    }
                }
                _ => {
                    return Err(anyhow!("Expected ParseOutput::Bytes, got unexpected variant"));
                }
            }
        }

        if !leftover.is_empty() {
            if tx.send(ParseOutput::Bytes(Arc::new(leftover))).await.is_err() {
                eprintln!("Warning: Receiver dropped when sending final leftover line");
            }
        }

        Ok(())
    });

    Ok((output_rx, task))
}

/// Helper function to allow error handling from the counting functions.
pub async fn join_with_error_handling<T>(task: JoinHandle<anyhow::Result<T, anyhow::Error>>) -> anyhow::Result<T, PipelineError> {
    match task.await {
        Ok(Ok(result)) => Ok(result),
        Ok(Err(e)) => Err(PipelineError::Other(e.into())),
        Err(join_err) => Err(PipelineError::Other(anyhow::anyhow!("Task failed: {}", join_err))),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;
    use std::io::Read;
    use tokio::process::Command;
    use tokio::task;
    use tokio::fs::File as TokioFile;
    use tokio::time::{self, Duration};
    use crate::utils::fastx::fastx_generator;
    use crate::config::defs::{RunConfig, StreamDataType};
    use std::sync::Arc;
    use rayon::ThreadPoolBuilder;
    use tokio::sync::{Semaphore, Mutex};
    use std::path::PathBuf;
    use crate::cli::Arguments;
    use tempfile::tempdir;
    use tokio::fs::metadata;
    use std::os::unix::fs::FileTypeExt;
    use crate::utils::fastx::SequenceRecord;
    use tokio::io::AsyncReadExt;


    // Helper function to create a RunConfig for tests
    fn create_test_run_config() -> Arc<RunConfig> {
        let args = Arguments {
            threads: 8, // Laptop
            ..Default::default()
        };
        Arc::new(RunConfig {
            cwd: PathBuf::from("."),
            ram_temp_dir: std::env::temp_dir(),
            out_dir: PathBuf::from("test"),
            args,
            thread_pool: Arc::new(ThreadPoolBuilder::new().num_threads(8).build().unwrap()),
            maximal_semaphore: Arc::new(Semaphore::new(8)),
            base_buffer_size: 5_000_000,
            input_size_mb: 100
        })
    }

    #[tokio::test]
    async fn test_t_junction_zero_streams() -> Result<()> {
        let stream = fastx_generator(10, 143, 35.0, 3.0).map(ParseOutput::Fastq);
        let result = t_junction(
            stream,
            0,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_zero_streams".to_string(),
            None,
        )
            .await;
        assert!(result.is_err());
        let error = result.unwrap_err();
        assert_eq!(
            error.to_string(),
            "No subscribers: cannot process stream in test_t_junction_zero_streams"
        );
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_two_records() -> Result<()> {
        let records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1".to_string(),
                desc: None,
                seq: Arc::new(b"ATCG".to_vec()),
                qual: Arc::new(b"IIII".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
        ];
        let stream = tokio_stream::iter(records);
        let (mut outputs, done_rx) = t_junction(
            stream,
            2,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_two_records".to_string(),
            None,
        )
            .await?;
        let mut output1 = ReceiverStream::new(outputs.pop().unwrap());
        let mut output2 = ReceiverStream::new(outputs.pop().unwrap());
        let mut records1 = Vec::new();
        let mut records2 = Vec::new();
        while let Some(record) = output1.next().await {
            records1.push(record);
        }
        while let Some(record) = output2.next().await {
            records2.push(record);
        }
        assert_eq!(records1.len(), 2);
        assert_eq!(records2.len(), 2);
        if let ParseOutput::Fastq(rec) = &records1[0] {
            assert_eq!(rec.id(), "read1");
        }
        if let ParseOutput::Fastq(rec) = &records2[0] {
            assert_eq!(rec.id(), "read1");
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_long_stream() -> Result<()> {
        let stream = fastx_generator(10_000, 143, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            2,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_long_stream".to_string(),
            None,
        )
            .await?;
        let mut output1 = ReceiverStream::new(outputs.pop().unwrap());
        let mut output2 = ReceiverStream::new(outputs.pop().unwrap());
        let mut records1 = Vec::new();
        let mut records2 = Vec::new();
        while let Some(record) = output1.next().await {
            records1.push(record);
        }
        while let Some(record) = output2.next().await {
            records2.push(record);
        }
        assert_eq!(records1.len(), 10_000);
        assert_eq!(records2.len(), 10_000);
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_ten_thousand_records_ten_streams() -> Result<()> {
        let stream = fastx_generator(10_000, 143, 35.0, 3.0).map(ParseOutput::Fastq);
        let (outputs, done_rx) = t_junction(
            stream,
            10,
            50_000,
            10_000,
            Some(0),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_ten_streams".to_string(),
            None,
        )
            .await?;
        let mut records = vec![Vec::new(); outputs.len()];
        let mut handles = Vec::new();

        for (i, rx) in outputs.into_iter().enumerate() {
            let handle = task::spawn(async move {
                let mut stream = ReceiverStream::new(rx);
                let mut local_records = Vec::new();
                while let Some(record) = stream.next().await {
                    local_records.push(record);
                }
                Ok::<_, anyhow::Error>(local_records)
            });
            handles.push((i, handle));
        }

        for (i, handle) in handles {
            let local_records = handle.await??;
            records[i] = local_records;
        }

        for record_set in &records {
            assert_eq!(record_set.len(), 10_000);
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_empty_stream() -> Result<()> {
        let stream = fastx_generator(0, 50, 35.0, 3.0).map(ParseOutput::Fastq);
        let (outputs, done_rx) = t_junction(
            stream,
            2,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_empty_stream".to_string(),
            None,
        )
            .await?;
        for rx in outputs {
            let mut stream = ReceiverStream::new(rx);
            assert!(stream.next().await.is_none(), "Empty stream should yield no items");
        }
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_junction_single_record() -> Result<()> {
        let stream = fastx_generator(1, 50, 35.0, 3.0).map(ParseOutput::Fastq);
        let (outputs, done_rx) = t_junction(
            stream,
            2,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_single_record".to_string(),
            None,
        )
            .await?;
        let mut handles = Vec::new();
        for rx in outputs {
            handles.push(task::spawn(async move {
                let mut records = Vec::new();
                let mut stream = ReceiverStream::new(rx);
                while let Some(record) = stream.next().await {
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
        let stream = fastx_generator(1_000, 50, 35.0, 3.0).map(ParseOutput::Fastq);
        let (outputs, done_rx) = t_junction(
            stream,
            2,
            500,
            10,
            Some(100),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_slow_consumer".to_string(),
            None,
        )
            .await?;
        let mut handles = Vec::new();
        for (i, rx) in outputs.into_iter().enumerate() {
            let handle = task::spawn(async move {
                let mut records = Vec::new();
                let mut stream = ReceiverStream::new(rx);
                let consumer_id = i;
                while let Some(record) = stream.next().await {
                    records.push(record);
                    if consumer_id == 1 {
                        sleep(Duration::from_millis(5)).await; // Simulate slow consumer
                    }
                }
                eprintln!("Consumer {} collected {} records", consumer_id, records.len());
                Ok::<_, anyhow::Error>(records)
            });
            handles.push(handle);
        }
        let all_records = time::timeout(Duration::from_secs(120), async {
            let mut all_records = Vec::with_capacity(2);
            for (i, handle) in handles.into_iter().enumerate() {
                let records = handle.await??;
                eprintln!("Task {} returned {} records", i, records.len());
                all_records.push(records);
            }
            assert_eq!(all_records.len(), 2, "Expected exactly 2 consumer outputs");
            assert_eq!(
                all_records[0].len(),
                1_000,
                "Consumer 0 should have all 1000 records, got {}",
                all_records[0].len()
            );
            assert_eq!(
                all_records[1].len(),
                1_000,
                "Consumer 1 should have all 1000 records, got {}",
                all_records[1].len()
            );
            Ok::<_, anyhow::Error>(all_records)
        })
            .await
            .map_err(|_| anyhow!("Test timed out after 120 seconds"))??;
        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_t_million_records_ten_streams() -> Result<()> {
        let num_records = 1_000_000;
        let stream = fastx_generator(num_records, 143, 35.0, 3.0).map(ParseOutput::Fastq);
        let (outputs, done_rx) = t_junction(
            stream,
            2,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_million_records_ten_streams".to_string(),
            None,
        )
            .await?;

        let mut handles = Vec::new();
        for rx in outputs {
            let handle = task::spawn(async move {
                let mut records = Vec::new();
                let mut stream = ReceiverStream::new(rx);
                while let Some(record) = stream.next().await {
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

        for i in 0..num_records {
            if let (ParseOutput::Fastq(rec0), ParseOutput::Fastq(rec1)) = (&all_records[0][i], &all_records[1][i]) {
                assert_eq!(
                    rec0.id(),
                    rec1.id(),
                    "Record {} IDs should match",
                    i
                );
                assert_eq!(
                    rec0.seq(),
                    rec1.seq(),
                    "Record {} sequences should match",
                    i
                );
                assert_eq!(
                    rec0.qual(),
                    rec1.qual(),
                    "Record {} quality scores should match",
                    i
                );
            }
        }

        done_rx.await??;
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_valid() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(100, 50, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_t_junction_stream_to_cmd_valid".to_string(),
            None,
        )
            .await?;
        let (child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = {
            let mut guard = child.lock().await;
            guard.stdout.take().unwrap()
        };
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        assert!(!output.is_empty(), "Output should contain data");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_valid_cat() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(2, 10, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_valid_cat".to_string(),
            None,
        )
            .await?;
        let (child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = {
            let mut guard = child.lock().await;
            guard.stdout.take().unwrap()
        };
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        assert!(!output.is_empty(), "Output should contain FASTQ data");
        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1"), "Output should contain first read ID");
        assert!(output_str.contains("@read2"), "Output should contain second read ID");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_valid_parseoutput() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(2, 10, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_valid_parse_output".to_string(),
            None,
        )
            .await?;
        let rx = outputs.pop().unwrap();
        let (tx, rx_parse) = mpsc::channel(100);
        tokio::spawn(async move {
            let mut stream = ReceiverStream::new(rx);
            while let Some(record) = stream.next().await {
                if tx.send(record).await.is_err() {
                    eprintln!("Failed to send ParseOutput");
                    break;
                }
            }
        });
        let (child, task, _err_task) = stream_to_cmd(
            config,
            rx_parse,
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = {
            let mut guard = child.lock().await;
            guard.stdout.take().unwrap()
        };
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        assert!(!output.is_empty(), "Output should contain FASTQ data");
        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1"), "Output should contain first read ID");
        assert!(output_str.contains("@read2"), "Output should contain second read ID");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_invalid_command() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(2, 10, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, _done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_invalid_cmd".to_string(),
            None,
        )
            .await?;
        let result = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "nonexistent_cmd",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await;
        assert!(result.is_err(), "Should fail for invalid command");
        let err = result.unwrap_err();
        assert!(
            err.to_string().contains("Failed to spawn nonexistent_cmd"),
            "Error should mention command name"
        );
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_empty_stream() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(0, 10, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_empty_stream".to_string(),
            None,
        )
            .await?;
        let (child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = {
            let mut guard = child.lock().await;
            guard.stdout.take().unwrap()
        };
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        assert!(output.is_empty(), "Output should be empty for empty stream");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_large_stream() -> Result<()> {
        let config = create_test_run_config();
        let num_records = 10_000;
        let stream = fastx_generator(num_records, 50, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_large_stream".to_string(),
            None,
        )
            .await?;
        let (child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = {
            let mut guard = child.lock().await;
            guard.stdout.take().unwrap()
        };
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;

        let output_str = String::from_utf8_lossy(&output);
        let mut lines = output_str.lines().peekable();
        let mut record_count = 0;

        while let Some(line) = lines.next() {
            if line.starts_with('@') {
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
    async fn test_stream_to_cmd_premature_exit() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(10, 10, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_premature_exit".to_string(),
            None,
        )
            .await?;
        let (child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "head",
            vec!["-n".to_string(), "1".to_string()],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = {
            let mut guard = child.lock().await;
            guard.stdout.take().unwrap()
        };
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??; // Task may error due to BrokenPipe, which is expected for "head -n 1"
        done_rx.await??;
        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1"), "Output should contain first read ID");
        assert!(!output_str.contains("@read2"), "Output should not contain second read ID");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_cmd_resource_cleanup() -> Result<()> {
        let config = create_test_run_config();
        let stream = fastx_generator(5, 10, 35.0, 3.0).map(ParseOutput::Fastq);
        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            50_000,
            10_000,
            Some(1),
            50,
            StreamDataType::IlluminaFastq,
            "test_stream_to_cmd_resource_cleanup".to_string(),
            None,
        )
            .await?;
        let (child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        task.await??;
        done_rx.await??;
        let stdin_closed = {
            let guard = child.lock().await;
            guard.stdin.is_none()
        };
        assert!(stdin_closed, "Stdin should be closed");
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_child_output_fastq() -> Result<()> {
        let mut cmd = Command::new("echo");
        cmd.arg("@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nHHHH\n");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let rx = parse_child_output(&mut child, ChildStream::Stdout, ParseMode::Fastq, 100).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut records = Vec::new();
        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
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
    async fn test_parse_child_output_bytes() -> Result<()> {
        let mut cmd = Command::new("echo");
        cmd.arg("test data");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let rx = parse_child_output(&mut child, ChildStream::Stdout, ParseMode::Bytes, 100).await?;
        let mut stream = ReceiverStream::new(rx);
        while let Some(ParseOutput::Bytes(chunk)) = stream.next().await {
            assert!(String::from_utf8_lossy(&*chunk).contains("test data"));
            break;
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_file_fastq() -> Result<()> {
        let _ = fs::remove_file("stream_to_file_test_illumina.fq");
        let num_records = 2;
        let mut records = fastx_generator(num_records, 150, 30.0, 8.0).map(ParseOutput::Fastq);
        let (tx, rx) = mpsc::channel(1024);

        tokio::spawn(async move {
            while let Some(record) = records.next().await {
                if tx.send(record).await.is_err() {
                    eprintln!("Failed to send record");
                    break;
                }
            }
        });

        let output_path = Path::new("stream_to_file_test_illumina.fq").to_path_buf();
        let write_task = tokio::spawn(stream_to_file(
            rx,
            output_path.clone(),
        ));
        write_task.await??;

        assert!(output_path.exists(), "Output file was not created");

        let file = File::open(&output_path).await?;
        let rx = parse_fastq(file, 1024).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut parsed_records = Vec::new();

        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
            parsed_records.push(record);
        }

        eprintln!("Parsed {} records from file", parsed_records.len());

        assert_eq!(
            parsed_records.len(),
            num_records,
            "Expected {} FASTQ records, found {}",
            num_records,
            parsed_records.len()
        );

        fs::remove_file(output_path)?;
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_file_no_records() -> Result<()> {
        let _ = fs::remove_file("stream_to_file_norecord_test.fq");
        let (tx, rx) = mpsc::channel(1024);
        let mut records = fastx_generator(0, 10, 35.0, 3.0).map(ParseOutput::Fastq);

        tokio::spawn(async move {
            while let Some(record) = records.next().await {
                if tx.send(record).await.is_err() {
                    eprintln!("Failed to send record");
                    break;
                }
            }
        });

        let output_path = Path::new("stream_to_file_norecord_test.fq").to_path_buf();
        let write_task = tokio::spawn(stream_to_file(
            rx,
            output_path.clone(),
        ));
        write_task.await??;

        assert!(
            output_path.exists(),
            "Output file was not created"
        );

        let file = File::open(&output_path).await?;
        let rx = parse_fastq(file, 1024).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut parsed_records = Vec::new();

        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
            parsed_records.push(record);
        }

        eprintln!("Parsed {} records from file", parsed_records.len());

        assert_eq!(
            parsed_records.len(),
            0,
            "Expected 0 FASTQ records, found {}",
            parsed_records.len()
        );

        fs::remove_file(output_path)?;
        Ok(())
    }

    #[tokio::test]
    async fn test_read_child_output_to_vec() -> Result<()> {
        let mut cmd = Command::new("echo");
        cmd.arg("line1\nline2\n\nline3");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "line1");
        assert_eq!(lines[1], "line2");
        assert_eq!(lines[2], "line3");
        Ok(())
    }

    #[tokio::test]
    async fn test_read_child_output_to_vec_stderr() -> Result<()> {
        let mut cmd = Command::new("sh");
        cmd.arg("-c").arg("echo error >&2");
        let mut child = cmd.stderr(std::process::Stdio::piped()).spawn()?;
        let lines = read_child_output_to_vec(&mut child, ChildStream::Stderr).await?;
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "error");
        Ok(())
    }

    #[tokio::test]
    async fn test_spawn_cmd_stderr_verbose() -> Result<()> {
        let config = create_test_run_config();
        let args = vec!["-c".to_string(), "echo error >&2".to_string()];
        let (mut child, stderr_task) = spawn_cmd(config, "sh", args, true).await?;
        let status = child.wait().await?;
        stderr_task.await??;
        assert!(status.success(), "Child process should exit successfully");
        Ok(())
    }

    #[tokio::test]
    async fn test_spawn_cmd_stderr_non_verbose() -> Result<()> {
        let config = create_test_run_config();
        let args = vec!["-c".to_string(), "echo error >&2".to_string()];
        let (mut child, stderr_task) = spawn_cmd(config, "sh", args, false).await?;
        let status = child.wait().await?;
        stderr_task.await??;
        assert!(status.success(), "Child process should exit successfully");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_file_bytes_no_clone() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let data = Arc::new(vec![b'A'; 1_000_000]); // 1MB chunk
        tokio::spawn(async move { tx.send(ParseOutput::Bytes(data)).await.unwrap(); });
        let path = PathBuf::from("test_bytes.fq");
        let task = tokio::spawn(stream_to_file(rx, path.clone())); // Pass raw Receiver
        task.await??;
        let file = std::fs::File::open(&path)?;
        let mut contents = Vec::new();
        std::io::BufReader::new(file).read_to_end(&mut contents)?; // Use std::io::Read
        assert_eq!(contents.len(), 1_000_000);
        std::fs::remove_file(&path)?;
        Ok(())
    }

    #[tokio::test]
    async fn test_create_fifo_success() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test_fifo").to_path_buf();
        create_fifo(&path).await?;
        let meta = metadata(&path).await?;
        assert!(meta.file_type().is_fifo(), "Created file should be a FIFO");
        Ok(())
    }

    #[tokio::test]
    async fn test_create_fifo_already_exists_as_file() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test_fifo").to_path_buf();
        // Create as regular file first
        tokio::fs::File::create(&path).await?;
        let result = create_fifo(&path).await;
        assert!(result.is_err(), "Should error if path exists as non-FIFO");
        if let Err(e) = result {
            assert!(e.to_string().contains("mkfifo failed"), "Error should mention mkfifo failure");
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_create_fifo_invalid_path() -> Result<()> {
        let invalid_path = PathBuf::from("/invalid/nonexistent/dir/test_fifo");
        let result = create_fifo(&invalid_path).await;
        assert!(result.is_err(), "Should error on invalid path");
        if let Err(e) = result {
            assert!(
                e.to_string().contains("Failed to execute mkfifo") || e.to_string().contains("mkfifo failed"),
                "Error should mention mkfifo execution failure"
            );
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_create_fifo_with_ram_temp_dir() -> Result<()> {
        let config = create_test_run_config();
        let fifo_path = config.ram_temp_dir.join("test_ram_fifo").to_path_buf();
        create_fifo(&fifo_path).await?;
        let meta = metadata(&fifo_path).await?;
        assert!(meta.file_type().is_fifo(), "Created file in ram_temp_dir should be a FIFO");
        // Cleanup
        tokio::fs::remove_file(&fifo_path).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_write_to_fifo_success() -> Result<()> {
        let dir = tempdir()?;
        let fifo_path = dir.path().join("test_fifo").to_path_buf();
        create_fifo(&fifo_path).await?;

        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let test_records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1".to_string(),
                desc: None,
                seq: Arc::new(b"ATCG".to_vec()),
                qual: Arc::new(b"IIII".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
        ];

        tokio::spawn(async move {
            for record in test_records {
                tx.send(record).await.unwrap();
            }
        });

        // Spawn reader (cat) to consume FIFO
        let mut reader = Command::new("cat")
            .arg(&fifo_path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat: {}", e))?;

        // Write to FIFO
        let write_handle = tokio::spawn(write_to_fifo(rx, fifo_path.clone()));
        let mut output = Vec::new();
        reader.stdout.take().unwrap().read_to_end(&mut output).await?;
        write_handle.await??;

        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1\nATCG\n+\nIIII\n"), "Should contain read1");
        assert!(output_str.contains("@read2\nGCTA\n+\nHHHH\n"), "Should contain read2");
        tokio::fs::remove_file(&fifo_path).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_write_to_fifo_empty_stream() -> Result<()> {
        let dir = tempdir()?;
        let fifo_path = dir.path().join("test_fifo").to_path_buf();
        create_fifo(&fifo_path).await?;

        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        drop(tx); // Simulate empty stream

        // Spawn reader (cat) to avoid blocking
        let mut reader = Command::new("cat")
            .arg(&fifo_path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat: {}", e))?;

        // Write (should complete immediately)
        let write_handle = tokio::spawn(write_to_fifo(rx, fifo_path.clone()));
        let mut output = Vec::new();
        reader.stdout.take().unwrap().read_to_end(&mut output).await?;
        write_handle.await??;

        // Verify empty output
        assert!(output.is_empty(), "Empty stream should produce no output");
        tokio::fs::remove_file(&fifo_path).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_write_to_fifo_invalid_path() -> Result<()> {
        let invalid_path = PathBuf::from("/blah/blah/flibiddy/blah");
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let write_handle = tokio::spawn(write_to_fifo(rx, invalid_path.clone()));
        let result = write_handle.await?;
        assert!(result.is_err(), "Should error on invalid FIFO path");
        if let Err(e) = result {
            assert!(
                e.to_string().contains("Failed to open FIFO"),
                "Error should mention FIFO open failure"
            );
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_write_to_fifo_with_ram_temp_dir() -> Result<()> {
        let config = create_test_run_config();
        let fifo_path = config.ram_temp_dir.join("test_ram_fifo").to_path_buf();
        create_fifo(&fifo_path).await?;

        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let test_record = ParseOutput::Fastq(SequenceRecord::Fastq {
            id: "read1".to_string(),
            desc: None,
            seq: Arc::new(b"ATCG".to_vec()),
            qual: Arc::new(b"IIII".to_vec()),
        });

        tokio::spawn(async move {
            tx.send(test_record).await.unwrap();
        });

        // Spawn reader
        let mut reader = Command::new("cat")
            .arg(&fifo_path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat: {}", e))?;

        // Write to FIFO
        let write_handle = tokio::spawn(write_to_fifo(rx, fifo_path.clone()));
        let mut output = Vec::new();
        reader.stdout.take().unwrap().read_to_end(&mut output).await?;
        write_handle.await??;

        // Verify
        let output_str = String::from_utf8_lossy(&output);
        assert!(output_str.contains("@read1\nATCG\n+\nIIII\n"), "Should contain read1 in ram_temp_dir");
        tokio::fs::remove_file(&fifo_path).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_paired() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let stream = ReceiverStream::new(rx);

        // Create interleaved FASTQ records (R1, R2, R1, R2)
        let test_records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/1".to_string(),
                desc: None,
                seq: Arc::new(b"ATCG".to_vec()),
                qual: Arc::new(b"IIII".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2/1".to_string(),
                desc: None,
                seq: Arc::new(b"CCCC".to_vec()),
                qual: Arc::new(b"JJJJ".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2/2".to_string(),
                desc: None,
                seq: Arc::new(b"GGGG".to_vec()),
                qual: Arc::new(b"KKKK".to_vec()),
            }),
        ];

        tokio::spawn(async move {
            for record in test_records {
                tx.send(record).await.unwrap();
            }
        });

        // Deinterleave
        let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_paired",
            true,
        )
            .await?;

        // Spawn readers for both FIFOs
        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(&r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle {
            handle.await??;
        }

        let r1_str = String::from_utf8_lossy(&r1_output);
        let r2_str = String::from_utf8_lossy(&r2_output);
        assert!(r1_str.contains("@read1/1\nATCG\n+\nIIII\n"), "R1 should contain read1/1");
        assert!(r1_str.contains("@read2/1\nCCCC\n+\nJJJJ\n"), "R1 should contain read2/1");
        assert!(r2_str.contains("@read1/2\nGCTA\n+\nHHHH\n"), "R2 should contain read1/2");
        assert!(r2_str.contains("@read2/2\nGGGG\n+\nKKKK\n"), "R2 should contain read2/2");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(&r2_fifo).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_single_end() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let stream = ReceiverStream::new(rx);

        // Single-end records
        let test_records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1".to_string(),
                desc: None,
                seq: Arc::new(b"ATCG".to_vec()),
                qual: Arc::new(b"IIII".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
        ];

        tokio::spawn(async move {
            for record in test_records {
                tx.send(record).await.unwrap();
            }
        });

        let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_single",
            false,
        )
            .await?;

        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(&r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle {
            handle.await??;
        }

        let r1_str = String::from_utf8_lossy(&r1_output);
        let r2_str = String::from_utf8_lossy(&r2_output);
        assert!(r1_str.contains("@read1\nATCG\n+\nIIII\n"), "R1 should contain read1");
        assert!(r1_str.contains("@read2\nGCTA\n+\nHHHH\n"), "R1 should contain read2");
        assert!(r2_str.is_empty(), "R2 should be empty for single-end");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(&r2_fifo).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_empty_stream() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        drop(tx); // Close immediately
        let stream = ReceiverStream::new(rx);

        let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_empty",
            true,
        )
            .await?;


        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(&r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;


        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle {
            handle.await??;
        }

        assert!(r1_output.is_empty(), "R1 should be empty for empty stream");
        assert!(r2_output.is_empty(), "R2 should be empty for empty stream");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(&r2_fifo).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_non_fastq() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let stream = ReceiverStream::new(rx);

        tx.send(ParseOutput::Bytes(Arc::new(b"invalid".to_vec()))).await?;

        let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_non_fastq",
            true,
        )
            .await?;

        // Spawn readers to avoid blocking
        let _r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()?;
        let _r2_reader = Command::new("cat")
            .arg(&r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()?;

        // Check for error
        let result = deinterleave_handle.await?;
        assert!(result.is_err(), "Should error on non-FASTQ input");
        if let Err(e) = result {
            assert!(
                e.to_string().contains("Non-FASTQ in deinterleave stream"),
                "Error should mention non-FASTQ input, got: {}",
                e
            );
        }
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle {
            handle.await??;
        }

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(&r2_fifo).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_with_ram_temp_dir() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let stream = ReceiverStream::new(rx);

        let test_records = vec![
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/1".to_string(),
                desc: None,
                seq: Arc::new(b"ATCG".to_vec()),
                qual: Arc::new(b"IIII".to_vec()),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read1/2".to_string(),
                desc: None,
                seq: Arc::new(b"GCTA".to_vec()),
                qual: Arc::new(b"HHHH".to_vec()),
            }),
        ];

        tokio::spawn(async move {
            for record in test_records {
                tx.send(record).await.unwrap();
            }
        });

        let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_ram_paired",
            true,
        )
            .await?;

        // Verify FIFOs exist in ram_temp_dir
        assert!(
            metadata(&r1_fifo).await?.file_type().is_fifo(),
            "R1 FIFO should be in ram_temp_dir"
        );
        assert!(
            metadata(&r2_fifo).await?.file_type().is_fifo(),
            "R2 FIFO should be in ram_temp_dir"
        );

        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(&r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle {
            handle.await??;
        }

        let r1_str = String::from_utf8_lossy(&r1_output);
        let r2_str = String::from_utf8_lossy(&r2_output);
        assert!(r1_str.contains("@read1/1\nATCG\n+\nIIII\n"), "R1 should contain read1/1");
        assert!(r2_str.contains("@read1/2\nGCTA\n+\nHHHH\n"), "R2 should contain read1/2");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(&r2_fifo).await.ok();
        Ok(())
    }

}