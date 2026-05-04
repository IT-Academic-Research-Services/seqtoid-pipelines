use sysinfo::System;
use anyhow::anyhow;
use anyhow::Result;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
use std::io;
use std::os::fd::AsRawFd;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::task::{Context, Poll};
use std::pin::Pin;


use tokio::fs::File;
use log::{self, LevelFilter, debug, info, error, warn};
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt, AsyncBufReadExt, BufReader, BufWriter, ReadBuf};
use tokio::process::{Child, Command};
use tokio::sync::{mpsc, oneshot};
use tokio::time::{Duration, sleep};
use tokio::time::{interval, Interval};
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
use std::os::unix::fs::PermissionsExt;
use tokio::fs::{self, set_permissions};
use std::fs::Permissions;
use futures::stream::BoxStream;
use futures::pin_mut;
use bytes::Bytes;
use rayon::prelude::*;

use crate::utils::fastx::{SequenceRecord, parse_header};
use crate::config::defs::{PipelineError, StreamDataType};
use crate::config::defs::{CoreAllocation, RunConfig, SimdLevel};
use crate::utils::system::{detect_cores_and_load, compute_stream_threads, detect_ram, generate_rng, compute_base_buffer_size, get_ram_temp_dir};


pub trait ToBytes {
    fn to_bytes(&self) -> Result<Bytes>;
}

impl ToBytes for SequenceRecord {
    fn to_bytes(&self) -> Result<Bytes> {
        let mut buffer = Vec::new();

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

        Ok(Bytes::from(buffer))
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ChildStream {
    Stdout,
    Stderr,
}

#[derive(Clone, Copy, Debug)]
pub enum ParseMode {
    Fastq,
    Fasta,
    Bytes,
    Lines,
}

// streams.rs
#[derive(Clone, Debug)]
pub enum ParseOutput {
    Fastq(SequenceRecord),
    Fasta(SequenceRecord),
    Bytes(Bytes),
}

impl ToBytes for ParseOutput {
    fn to_bytes(&self) -> Result<Bytes> {
        match self {
            ParseOutput::Fastq(record) => record.to_bytes(),
            ParseOutput::Fasta(record) => record.to_bytes(),
            ParseOutput::Bytes(bytes) => Ok(bytes.clone()),
        }
    }
}

/// Converts mpsc::Receiver<ParseOutput> output to AsyncRead suitable for parse_fasta / parse_fastq.
pub struct ChannelReader {
    rx: mpsc::Receiver<ParseOutput>,
    buffer: Option<Bytes>,
    offset: usize,
    closed: bool,
}

impl ChannelReader {
    pub fn new(rx: mpsc::Receiver<ParseOutput>) -> Self {
        Self {
            rx,
            buffer: None,
            offset: 0,
            closed: false,
        }
    }
}

impl AsyncRead for ChannelReader {
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        if self.closed {
            return Poll::Ready(Ok(()));
        }

        loop {
            if let Some(bytes) = self.buffer.as_ref().cloned() {
                let remaining = &bytes[self.offset..];
                if remaining.is_empty() {
                    self.buffer = None;
                    self.offset = 0;
                    continue;
                }

                let to_copy = std::cmp::min(remaining.len(), buf.remaining());
                buf.put_slice(&remaining[..to_copy]);
                self.offset += to_copy;

                if self.offset >= bytes.len() {
                    self.buffer = None;
                    self.offset = 0;
                }

                return Poll::Ready(Ok(()));
            }

            match self.rx.poll_recv(cx) {
                Poll::Ready(Some(item)) => match item {
                    ParseOutput::Bytes(bytes) => {
                        self.buffer = Some(bytes);
                        self.offset = 0;
                    }
                    _ => {
                        return Poll::Ready(Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Expected ParseOutput::Bytes, got another variant",
                        )));
                    }
                },
                Poll::Ready(None) => {
                    self.closed = true;
                    return Poll::Ready(Ok(()));
                }
                Poll::Pending => return Poll::Pending,
            }
        }
    }
}

/// Fans oiut data piece by peiecve to donstream consumers
/// NOTE: does NOT give each cocnusmer ALL data. NOT used for "splitting" streams!!
/// if feeding reads a b c to two consumers
/// Task A: recv() → gets `a`
/// Task B: recv() → gets `b`
/// Task A: recv() → gets `c`
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
        available_ram = system.total_memory().max(8_000_000_000) / 4;
    }

    const MAX_PROCESSES: usize = 4;
    const RAM_FRACTION: f64 = 0.5;
    const MIN_BUFFER_PER_STREAM: usize = 10_000;
    const MAX_BUFFER_PER_STREAM: usize = 2_000_000;

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
        warn!("Failed to detect available RAM in {}, using fallback", label);
        (base_buffer_size * n_outputs.max(1) * 2).max(min_buffer_size)
    };
    let buffer_size = (base_buffer_size * n_outputs.max(1))
        .clamp(min_buffer_size, max_buffer_size.min(MAX_BUFFER_PER_STREAM));

    let (done_tx, done_rx) = oneshot::channel::<Result<(), anyhow::Error>>();
    let mut output_txs: Vec<(usize, mpsc::Sender<ParseOutput>)> = Vec::with_capacity(n_outputs);
    let mut output_rxs = Vec::with_capacity(n_outputs);

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
            let mut lowest_capacity_ratio = 1.0f64; // tracks the worst receiver

            for (i, tx) in output_txs.into_iter() {
                let capacity_ratio = tx.capacity() as f64 / buffer_size as f64;
                lowest_capacity_ratio = lowest_capacity_ratio.min(capacity_ratio);

                match tx.try_send(item.clone()) {
                    Ok(()) => active_txs.push((i, tx)),
                    Err(mpsc::error::TrySendError::Full(_)) => {
                        if tx.send(item.clone()).await.is_err() {
                            error!("{}: Receiver {} dropped at item {}", label, i, item_count);
                            dropped_receivers.push(i);
                        } else {
                            active_txs.push((i, tx));
                        }
                    }
                    Err(mpsc::error::TrySendError::Closed(_)) => {
                        error!("{}: Receiver {} dropped at item {}", label, i, item_count);
                        dropped_receivers.push(i);
                    }
                }
            }

            output_txs = active_txs;
            if output_txs.is_empty() {
                error!("{}: All receivers dropped at item {}. Dropped: {:?}", label, item_count, dropped_receivers);
                let _ = done_tx.send(Err(anyhow!("All receivers dropped at item {}", item_count)));
                return;
            }

            // Per-item backpressure: only slow down when the SLOWEST receiver is truly full.
            if lowest_capacity_ratio < 0.05 {
                if item_count % 5000 == 0 {
                    debug!(
                        "{}: Heavy backpressure detected (lowest capacity {:.1}%). Sleep 10ms.",
                        label, lowest_capacity_ratio * 100.0
                    );
                }
                sleep(Duration::from_millis(10)).await;
            } else if lowest_capacity_ratio < 0.20 {
                tokio::task::yield_now().await;
            }

            if item_count % stall_threshold == 0 {
                if let Some(sleep_ms) = stream_sleep_ms {
                    debug!("{}: Periodic stall check sleep {}ms at item {}", label, sleep_ms, item_count);
                    sleep(Duration::from_millis(sleep_ms)).await;
                }
            }

            if last_progress_time.elapsed() > Duration::from_secs(60) {
                debug!("{}: Progress update — processed {} items", label, item_count);
                last_progress_time = Instant::now();
            }
        }

        for (i, _) in output_txs {
            debug!("{}: Closed sender for receiver {} after {} items", label, i, item_count);
        }

        if let Some(n) = notify {
            n.notify_waiters();
        }

        if dropped_receivers.is_empty() {
            let _ = done_tx.send(Ok(()));
        } else {
            let _ = done_tx.send(Err(anyhow!(
                "{}: {} receivers dropped: {:?}",
                label,
                dropped_receivers.len(),
                dropped_receivers
            )));
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
    stderr_log_path: Option<PathBuf>,
) -> Result<(
    Arc<tokio::sync::Mutex<Child>>,
    JoinHandle<Result<(), anyhow::Error>>,
    JoinHandle<Result<(), anyhow::Error>>,
)> {
    let core_allocation = config.get_core_allocation(cmd_tag, None);
    let _permit = if core_allocation == CoreAllocation::Maximal {
        Some(config.maximal_semaphore.clone().acquire_owned().await?)
    } else {
        None
    };

    let writer_capacity = crate::utils::system::compute_buffer_size(
        &config,
        "stream_to_cmd",
        data_type,
        1.8, // 1.0 = normal, 1.5–2.0 for very hot paths later
    );

    let batch_size_bytes = writer_capacity / 4; // keep ~4:1 batching ratio

    let cmd_tag_owned = cmd_tag.to_string();
    let cmd_tag_err_owned = cmd_tag.to_string();
    let stderr_log_path_clone = stderr_log_path.clone();

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
        let mut total_written = 0usize;

        let mut stream = ReceiverStream::new(rx);
        while let Some(item) = stream.next().await {
            match &item {
                ParseOutput::Bytes(bytes) => batch.extend_from_slice(bytes.as_ref()),
                ParseOutput::Fastq(record) => {
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
                        let mut guard = child_clone.lock().await;
                        if let Ok(Some(status)) = guard.try_wait() {
                            if status.success() {
                                debug!(
                                    "Ignoring BrokenPipe in {} after successful exit",
                                    cmd_tag_owned
                                );
                                break;
                            } else {
                                return Err(anyhow!(
                                    "BrokenPipe in {} with failed child exit: {:?}",
                                    cmd_tag_owned,
                                    status
                                ));
                            }
                        } else {
                            return Err(anyhow!(
                                "BrokenPipe in {} before child exit: {}",
                                cmd_tag_owned,
                                e
                            ));
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
                        if !status.success() {
                            return Err(anyhow!("Final BrokenPipe with failed exit: {:?}", status));
                        }
                    } else {
                        return Err(anyhow!("Final BrokenPipe before exit: {}", e));
                    }
                } else {
                    return Err(anyhow!("Final write error: {}", e));
                }
            } else {
                writer.flush().await?;
            }
        }

        writer.flush().await?;
        drop(writer);

        let mut guard = child_clone.lock().await;
        let status = guard.wait().await?;
        if !status.success() {
            return Err(anyhow!("Child {} failed: {:?}", cmd_tag_owned, status));
        }

        debug!(
            "[{} stdin] completed, wrote {} bytes",
            cmd_tag_owned, total_written
        );

        Ok(())
    });

    let stderr_task = tokio::spawn(async move {
        let mut reader = BufReader::with_capacity(1_048_576, stderr);
        let mut buffer = vec![0u8; 8192];

        let mut stderr_file = if let Some(path) = stderr_log_path_clone {
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent).await.map_err(|e| {
                    anyhow!("Failed to create stderr log dir {}: {}", parent.display(), e)
                })?;
            }
            Some(TokioFile::create(&path).await.map_err(|e| {
                anyhow!("Failed to create stderr log file {}: {}", path.display(), e)
            })?)
        } else {
            None
        };

        loop {
            let n = reader
                .read(&mut buffer)
                .await
                .map_err(|e| anyhow!("Failed to read stderr: {}", e))?;

            if n == 0 {
                break;
            }

            if let Some(file) = stderr_file.as_mut() {
                file.write_all(&buffer[..n]).await.map_err(|e| {
                    anyhow!("Failed to write stderr log for {}: {}", cmd_tag_err_owned, e)
                })?;
            }

            if verbose || stderr_log_path.is_some() {
                let stderr_chunk = String::from_utf8_lossy(&buffer[..n]);
                debug!("[{} stderr]: {}", cmd_tag_err_owned, stderr_chunk);
            }
        }

        if let Some(file) = stderr_file.as_mut() {
            file.flush().await.map_err(|e| {
                anyhow!("Failed to flush stderr log for {}: {}", cmd_tag_err_owned, e)
            })?;
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
    stderr_log_path: Option<PathBuf>,
) -> Result<(Child, JoinHandle<Result<(), anyhow::Error>>)> {
    let core_allocation = config.get_core_allocation(cmd_tag, None);
    let _permit = if core_allocation == CoreAllocation::Maximal {
        Some(config.maximal_semaphore.clone().acquire_owned().await?)
    } else {
        None
    };

    let cmd_tag_owned = cmd_tag.to_string();
    eprintln!("[{} cmd]: {}", cmd_tag, args.join(" "));

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
        let stderr_log_path_clone = stderr_log_path.clone();

        tokio::spawn(async move {
            let mut reader = BufReader::with_capacity(1024 * 1024, stderr);
            let mut buffer = vec![0u8; 8192];

            let mut stderr_file = if let Some(path) = stderr_log_path_clone {
                if let Some(parent) = path.parent() {
                    fs::create_dir_all(parent).await.map_err(|e| {
                        anyhow!("Failed to create stderr log dir {}: {}", parent.display(), e)
                    })?;
                }
                Some(TokioFile::create(&path).await.map_err(|e| {
                    anyhow!("Failed to create stderr log file {}: {}", path.display(), e)
                })?)
            } else {
                None
            };

            loop {
                let n = reader.read(&mut buffer).await
                    .map_err(|e| anyhow!("Failed to read stderr: {}", e))?;

                if n == 0 {
                    break;
                }

                if let Some(file) = stderr_file.as_mut() {
                    file.write_all(&buffer[..n]).await.map_err(|e| {
                        anyhow!("Failed to write stderr log for {}: {}", cmd_tag_clone, e)
                    })?;
                }

                if verbose || stderr_log_path.is_some() {
                    let chunk = String::from_utf8_lossy(&buffer[..n]);
                    eprint!("[{} stderr]: {}", cmd_tag_clone, chunk);
                    if !chunk.ends_with('\n') {
                        eprintln!();
                    }
                }
            }

            if let Some(file) = stderr_file.as_mut() {
                file.flush().await.map_err(|e| {
                    anyhow!("Failed to flush stderr log for {}: {}", cmd_tag_clone, e)
                })?;
            }

            Ok(())
        })
    };

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
    config: &RunConfig,
) -> Result<mpsc::Receiver<ParseOutput>> {

    match (stream, mode) {
        // Fastq
        (ChildStream::Stdout, ParseMode::Fastq) => {
            let stdout = child.stdout.take()
                .ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_fastq(stdout, config, StreamDataType::IlluminaFastq).await
        }
        (ChildStream::Stderr, ParseMode::Fastq) => {
            let stderr = child.stderr.take()
                .ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_fastq(stderr, config, StreamDataType::IlluminaFastq).await
        }

        // Fasta
        (ChildStream::Stdout, ParseMode::Fasta) => {
            let stdout = child.stdout.take()
                .ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_fasta(stdout, config, StreamDataType::JustBytes).await
        }
        (ChildStream::Stderr, ParseMode::Fasta) => {
            let stderr = child.stderr.take()
                .ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_fasta(stderr, config, StreamDataType::JustBytes).await
        }

        // Bytes
        (ChildStream::Stdout, ParseMode::Bytes) => {
            let stdout = child.stdout.take()
                .ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_bytes(stdout, config, StreamDataType::JustBytes).await
        }
        (ChildStream::Stderr, ParseMode::Bytes) => {
            let stderr = child.stderr.take()
                .ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_bytes(stderr, config, StreamDataType::JustBytes).await
        }

        // Lines
        (ChildStream::Stdout, ParseMode::Lines) => {
            let stdout = child.stdout.take()
                .ok_or_else(|| anyhow!("Child stdout not available"))?;
            parse_lines(stdout, config, StreamDataType::JustBytes).await
        }
        (ChildStream::Stderr, ParseMode::Lines) => {
            let stderr = child.stderr.take()
                .ok_or_else(|| anyhow!("Child stderr not available"))?;
            parse_lines(stderr, config, StreamDataType::JustBytes).await
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
    config: &RunConfig,           // ← added
    data_type: StreamDataType,    // ← added (IlluminaFastq or OntFastq)
) -> Result<mpsc::Receiver<ParseOutput>> {

    let channel_buffer = crate::utils::system::compute_buffer_size(
        config,
        "parse_fastq",
        data_type,
        0.55,                         // parsing benefits from tighter buffers + backpressure
    );

    let (tx, rx) = mpsc::channel(channel_buffer);

    let mut reader = BufReader::with_capacity(1024 * 1024, reader);
    let mut buffer = String::new();
    let mut count = 0usize;

    tokio::spawn(async move {
        loop {
            buffer.clear();
            let bytes_read = match reader.read_line(&mut buffer).await {
                Ok(n) => n,
                Err(e) => {
                    error!("parse_fastq: Error reading FASTQ line: {}", e);
                    return;
                }
            };
            if bytes_read == 0 {
                debug!("parse_fastq: Reached end of input, processed {} records", count);
                break;
            }

            let id_line = buffer.trim_end();
            if !id_line.starts_with('@') {
                error!("parse_fastq: Invalid FASTQ format: expected '@', got '{}'", id_line);
                return;
            }

            let id = id_line[1..].to_string();
            let desc = None;

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                error!("parse_fastq: Error reading sequence line");
                return;
            }
            let seq = buffer.trim_end().as_bytes().to_vec();
            if seq.is_empty() {
                error!("parse_fastq: Missing sequence");
                return;
            }

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                error!("parse_fastq: Error reading plus line");
                return;
            }
            let plus = buffer.trim_end();
            if plus != "+" {
                error!("parse_fastq: Invalid FASTQ format: expected '+', got '{}'", plus);
                return;
            }

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                error!("parse_fastq: Error reading quality line");
                return;
            }
            let qual = buffer.trim_end().as_bytes().to_vec();
            if qual.len() != seq.len() {
                error!("parse_fastq: Sequence and quality lengths do not match: seq_len={}, qual_len={}", seq.len(), qual.len());
                return;
            }

            let record = SequenceRecord::Fastq {
                id,
                desc,
                seq: Arc::new(seq),
                qual: Arc::new(qual),
            };

            if tx.send(ParseOutput::Fastq(record)).await.is_err() {
                error!("parse_fastq: Receiver dropped after {} records", count);
                break;
            }
            count += 1;
        }
        debug!("parse_fastq: Completed parsing, sent {} FASTQ records", count);
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
    config: &RunConfig,
    data_type: StreamDataType,    // (JustBytes is fine)
) -> Result<mpsc::Receiver<ParseOutput>> {

    let channel_buffer = crate::utils::system::compute_buffer_size(
        config,
        "parse_fasta",
        data_type,
        0.6,                          // parsing path — tighter buffers
    );

    let (tx, rx) = mpsc::channel(channel_buffer);

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
                seq: Arc::new(current_seq),
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
    config: &RunConfig,
    data_type: StreamDataType,
) -> Result<mpsc::Receiver<ParseOutput>> {

    let channel_buffer = crate::utils::system::compute_buffer_size(
        config,
        "parse_bytes",
        data_type,
        0.6,                          // parsing path — tighter buffers for backpressure
    );

    let (tx, rx) = mpsc::channel(channel_buffer);

    tokio::spawn(async move {
        let mut reader = BufReader::with_capacity(1024 * 1024, reader);
        let mut buffer = vec![0u8; 8192];

        loop {
            let bytes_read = match reader.read(&mut buffer).await {
                Ok(n) => n,
                Err(e) => {
                    error!("Error reading bytes: {}", e);
                    return;
                }
            };
            if bytes_read == 0 {
                break;
            }

            if tx.send(ParseOutput::Bytes(Bytes::from(buffer[..bytes_read].to_vec())))
                .await
                .is_err()
            {
                error!("Receiver dropped while sending bytes");
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
    config: &RunConfig,
    data_type: StreamDataType,
) -> Result<mpsc::Receiver<ParseOutput>> {

    let channel_buffer = crate::utils::system::compute_buffer_size(
        config,
        "parse_lines",
        data_type,
        0.6,
    );

    let (tx, rx) = mpsc::channel(channel_buffer);

    let mut reader = BufReader::with_capacity(1024 * 1024, reader);

    tokio::spawn(async move {
        let mut line = String::new();

        loop {
            line.clear();

            let bytes_read = match reader.read_line(&mut line).await {
                Ok(n) => n,
                Err(e) => {
                    error!("Error reading line: {}", e);
                    return;
                }
            };

            if bytes_read == 0 {
                break;
            }

            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }

            let bytes = trimmed.as_bytes().to_vec();

            if tx
                .send(ParseOutput::Bytes(Bytes::from(bytes)))
                .await
                .is_err()
            {
                error!("Receiver dropped while sending line");
                break;
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
pub async fn read_child_output_to_vec(child: &mut Child, stream: ChildStream, config: &RunConfig) -> Result<Vec<String>> {
    let rx = parse_child_output(child, stream, ParseMode::Bytes, &config).await?;
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

pub async fn create_fifo(path: &PathBuf) -> Result<(), PipelineError> {
    // Check for and remove existing FIFO with retry
    let mut exists_as_non_fifo = false;
    for attempt in 1..=3 {
        match fs::metadata(path).await {
            Ok(meta) => {
                if meta.file_type().is_fifo() {
                    debug!("Attempt {} to remove existing FIFO: {}", attempt, path.display());
                    match fs::remove_file(path).await {
                        Ok(_) => {
                            debug!("Successfully removed existing FIFO: {}", path.display());
                            break;
                        }
                        Err(e) => {
                            if attempt < 3 {
                                sleep(Duration::from_millis(100)).await;
                                continue;
                            }
                            return Err(PipelineError::IOError(format!(
                                "Failed to remove existing FIFO {} after {} attempts: {}",
                                path.display(),
                                attempt,
                                e
                            )));
                        }
                    }
                } else {
                    exists_as_non_fifo = true;
                    break;
                }
            }
            Err(_) if attempt == 1 => {
                debug!("No existing FIFO found: {}", path.display());
                break;
            }
            Err(e) => {
                error!("Metadata error on attempt {}: {}", attempt, e);
                if attempt < 3 {
                    sleep(Duration::from_millis(100)).await;
                    continue;
                }
                return Err(PipelineError::IOError(e.to_string()));
            }
        }
    }

    if exists_as_non_fifo {
        return Err(PipelineError::IOError(format!(
            "Path {} exists and is not a FIFO; refusing to remove non-FIFO file",
            path.display()
        )));
    }

    // Create new FIFO
    let status = Command::new("mkfifo")
        .arg(path)
        .status()
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    if !status.success() {
        return Err(PipelineError::IOError(format!("mkfifo failed for {}", path.display())));
    }

    // Set permissions
    set_permissions(path, Permissions::from_mode(0o666))
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

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
    debug!("Opened FIFO for writing: {}", fifo_path.display());

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
            debug!("Wrote {} bytes to FIFO: {}", bytes_written, fifo_path.display());
        }
    }

    writer.flush().await
        .map_err(|e| anyhow!("Flush to FIFO {} failed: {}", fifo_path.display(), e))?;
    debug!("Finished writing {} bytes to FIFO: {}", bytes_written, fifo_path.display());
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
                let bytes_output = ParseOutput::Bytes(Bytes::from(bytes));
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
        debug!("deinterleave_fastq_stream: Completed deinterleaving {} records", record_count);
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
) -> Result<(
    PathBuf,
    Option<PathBuf>,
    JoinHandle<Result<(), anyhow::Error>>,
    JoinHandle<Result<(), anyhow::Error>>,
    Option<JoinHandle<Result<(), anyhow::Error>>>,
), PipelineError> {

    let ram_temp_dir = config.ram_temp_dir.clone();
    let r1_fifo = ram_temp_dir.join(format!("{}_R1.fq", sample_base));
    let r2_fifo = paired.then(|| ram_temp_dir.join(format!("{}_R2.fq", sample_base)));

    // Clean up stale FIFOs
    let fifos_to_clean: Vec<&PathBuf> = std::iter::once(&r1_fifo)
        .chain(r2_fifo.as_ref())
        .collect();

    for fifo in fifos_to_clean {
        for attempt in 1..=3 {
            if metadata(fifo).await.is_ok() {
                match remove_file(fifo).await {
                    Ok(_) => break,
                    Err(e) => {
                        if attempt == 3 {
                            return Err(PipelineError::IOError(format!(
                                "Failed to remove FIFO {} after 3 attempts: {}", fifo.display(), e
                            )));
                        }
                        sleep(Duration::from_millis(100)).await;
                    }
                }
            } else {
                break;
            }
        }
    }

    // Create FIFOs
    create_fifo(&r1_fifo)
        .await
        .map_err(|e| PipelineError::IOError(format!("Failed to create R1 FIFO {}: {}", r1_fifo.display(), e)))?;

    if let Some(ref r2_path) = r2_fifo {
        create_fifo(r2_path)
            .await
            .map_err(|e| PipelineError::IOError(format!("Failed to create R2 FIFO {}: {}", r2_path.display(), e)))?;
    }

    // Dynamic buffer for deinterleave channels
    let channel_buffer = crate::utils::system::compute_buffer_size(
        &config,
        "deinterleave_fastq",
        StreamDataType::IlluminaFastq,
        1.6,                          // hot path for deinterleave
    );

    let (r1_tx, r1_rx) = mpsc::channel::<ParseOutput>(channel_buffer);
    let (r2_tx, r2_rx) = mpsc::channel::<ParseOutput>(channel_buffer);

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
                                .map_err(|_| PipelineError::StreamDataDropped)?;
                        } else {
                            r2_tx.send(ParseOutput::Fastq(record)).await
                                .map_err(|_| PipelineError::StreamDataDropped)?;
                        }
                        is_r1 = !is_r1;
                    } else {
                        r1_tx.send(ParseOutput::Fastq(record)).await
                            .map_err(|_| PipelineError::StreamDataDropped)?;
                    }
                }
                _ => return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in deinterleave stream".to_string()).into()),
            }
        }

        if paired && !is_r1 {
            return Err(PipelineError::InvalidFastqFormat(format!(
                "Incomplete paired-end stream: {} records", record_count
            )).into());
        }

        debug!("deinterleave_fastq_stream_to_fifos: Completed {} records", record_count);
        Ok(())
    });

    // Writer tasks
    let r1_write_handle = tokio::spawn(write_to_fifo(r1_rx, r1_fifo.clone()));
    let r2_write_handle = if paired {
        Some(tokio::spawn(write_to_fifo(r2_rx, r2_fifo.clone().unwrap())))
    } else {
        None
    };

    Ok((
        r1_fifo,
        r2_fifo,
        deinterleave_handle,
        r1_write_handle,
        r2_write_handle,
    ))
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
    debug!("Created R1 file descriptor: {}", r1_fd_str);

    // Hold stdout to prevent premature closure
    let r1_child_stdout = Arc::new(r1_child_stdout); // Keep alive
    let r1_cat_handle: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
        r1_child.wait().await.map_err(|e| anyhow!("R1 cat process failed: {}", e))?;
        debug!("R1 cat process completed");
        Ok(())
    });

    // Spawn R1 writer
    let r1_fifo_clone = r1_fifo.clone(); // Clone for use in closure
    let writers_ready_clone = writers_ready.clone();
    let r1_child_stdout_clone = r1_child_stdout.clone(); // Keep stdout alive
    let r1_handle = tokio::spawn(async move {
        debug!("Opened R1 process substitution for writing");
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
        debug!("Finished writing {} bytes to R1 process substitution", bytes_written);
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
        debug!("Created R2 file descriptor: {}", r2_fd_str);

        let r2_child_stdout = Arc::new(r2_child_stdout); // Keep alive
        let r2_cat_handle: JoinHandle<Result<(), anyhow::Error>> = tokio::spawn(async move {
            r2_child.wait().await.map_err(|e| anyhow!("R2 cat process failed: {}", e))?;
            debug!("R2 cat process completed");
            Ok(())
        });

        let writers_ready_clone = writers_ready.clone();
        let r2_child_stdout_clone = r2_child_stdout.clone(); // Keep stdout alive
        let r2_handle = tokio::spawn(async move {
            debug!("Opened R2 process substitution for writing");
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
            debug!("Finished writing {} bytes to R2 process substitution", bytes_written);
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

pub async fn interleave_fastq_streams(
    r1_rx: mpsc::Receiver<ParseOutput>,
    r2_rx: mpsc::Receiver<ParseOutput>,
    buffer_size: usize,
) -> Result<(mpsc::Receiver<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>), PipelineError> {
    let (inter_tx, inter_rx) = mpsc::channel(buffer_size);
    let inter_task = tokio::spawn(async move {
        let mut r1_stream = ReceiverStream::new(r1_rx);
        let mut r2_stream = ReceiverStream::new(r2_rx);
        let mut record_count = 0;
        while let (Some(r1), Some(r2)) = (r1_stream.next().await, r2_stream.next().await) {
            record_count += 2;
            match (&r1, &r2) {
                (ParseOutput::Fastq(_), ParseOutput::Fastq(_)) => {
                    inter_tx.send(r1).await
                        .map_err(|_| PipelineError::StreamDataDropped)
                        .map_err(Into::<anyhow::Error>::into)?;
                    inter_tx.send(r2).await
                        .map_err(|_| PipelineError::StreamDataDropped)
                        .map_err(Into::<anyhow::Error>::into)?;
                }
                _ => return Err(PipelineError::InvalidFastqFormat("Non-FASTQ in interleave stream".to_string()))
                    .map_err(Into::<anyhow::Error>::into),
            }
        }
        // Check if one stream ended early
        if r1_stream.next().await.is_some() || r2_stream.next().await.is_some() {
            return Err(PipelineError::InvalidFastqFormat(format!(
                "Uneven paired-end streams: {} records", record_count
            )))
                .map_err(Into::<anyhow::Error>::into);
        }
        debug!("interleave_fastq_streams: Interleaved {} records", record_count);
        Ok(())
    });
    Ok((inter_rx, inter_task))
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
                            return Err(anyhow::anyhow!("Failed to send FASTA record"));
                        }
                    }
                    SequenceRecord::Fastq { .. } => {
                        error!("Unexpected FASTQ record in consensus_stats_stream");
                    }
                }
            } else {
                error!("Unexpected non-FASTA item in stream");
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
    let mut leftover: Vec<u8> = Vec::new();

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(rx);

        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(bytes) => {
                    let mut bytes = bytes.to_vec();

                    if !leftover.is_empty() {
                        leftover.append(&mut bytes);
                        bytes = leftover;
                        leftover = Vec::new();
                    }

                    let mut start = 0;
                    for i in 0..bytes.len() {
                        if bytes[i] == b'\n' {
                            let line = Bytes::copy_from_slice(&bytes[start..=i]);
                            if tx.send(ParseOutput::Bytes(line)).await.is_err() {
                                warn!("Warning: Receiver dropped in bytes_to_lines, stopping line processing");
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
            if tx.send(ParseOutput::Bytes(Bytes::from(leftover))).await.is_err() {
                warn!("Warning: Receiver dropped when sending final leftover line");
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



pub fn monitor_stream(
    input: ReceiverStream<ParseOutput>,
    stage_name: &str,
    log_interval: Duration,
) -> ReceiverStream<ParseOutput> {
    let count = Arc::new(AtomicUsize::new(0));
    let stage = stage_name.to_string();
    let count_clone = count.clone();

    let stop_notify = Arc::new(Notify::new());
    let stop_notify_clone = stop_notify.clone();
    let stage_clone = stage.clone();

    tokio::spawn(async move {
        let mut heartbeat: Interval = interval(log_interval);
        let mut last_count: usize = 0;
        let start = Instant::now();

        loop {
            tokio::select! {
                _ = heartbeat.tick() => {
                    let current = count_clone.load(Ordering::Relaxed);
                    let delta = current - last_count;
                    let velocity = (delta as f64) / log_interval.as_secs_f64();

                    if delta == 0 {
                        warn!(
                            "{}: Stalled (0 records in last {:?}, total={})",
                            stage_clone,
                            log_interval,
                            current
                        );
                    } else {
                        info!(
                            "{}: Output velocity: {:.0} records/s ({} total, +{} since last tick, elapsed {:?})",
                            stage_clone,
                            velocity,
                            current,
                            delta,
                            start.elapsed()
                        );
                    }

                    last_count = current;
                }
                _ = stop_notify_clone.notified() => {
                    let final_count = count_clone.load(Ordering::Relaxed);
                    debug!(
                        "{}: Logging task stopped after {} total records in {:?}",
                        stage_clone,
                        final_count,
                        start.elapsed()
                    );
                    break;
                }
            }
        }
    });

    let mut input_receiver = input.into_inner();
    let buf_size = input_receiver.capacity();
    let (forward_tx, forward_rx) = mpsc::channel::<ParseOutput>(buf_size);
    let stop_notify_forward = stop_notify.clone();
    let stage_forward = stage.clone();
    let count_forward = count.clone();

    tokio::spawn(async move {
        struct DropGuard {
            notify: Arc<Notify>,
        }

        impl Drop for DropGuard {
            fn drop(&mut self) {
                self.notify.notify_one();
            }
        }

        let _guard = DropGuard {
            notify: stop_notify,
        };

        let mut forwarded = 0usize;
        let mut last_log = Instant::now();

        while let Some(item) = input_receiver.recv().await {
            count_forward.fetch_add(1, Ordering::Relaxed);
            forwarded += 1;

            if forward_tx.send(item).await.is_err() {
                warn!(
                    "{}: Forward send failed (downstream dropped after {} records)",
                    stage_forward,
                    forwarded
                );
                break;
            }

            if last_log.elapsed() >= Duration::from_secs(10) {
                let current = count_forward.load(Ordering::Relaxed);
                info!(
                    "{}: forwarded {} records total (receiver total={})",
                    stage_forward,
                    forwarded,
                    current
                );
                last_log = Instant::now();
            }
        }

        let final_total = count_forward.load(Ordering::Relaxed);
        info!(
            "{}: Stream complete ({} total records)",
            stage_forward,
            final_total
        );

        drop(forward_tx);
        stop_notify_forward.notify_one();
    });

    ReceiverStream::new(forward_rx)
}


/// Guarantees: zero data drops, respects t_junction backpressure, cleanup_receivers untouched.
///
/// Speed improvements
/// - reuses the input batch buffer allocation instead of `std::mem::take`
/// - reserves space for incoming chunks before extending
/// - surfaces downstream send failure instead of silently ignoring it
pub fn batch_rayon_process<F>(
    config: Arc<RunConfig>,
    mut input: ReceiverStream<ParseOutput>,
    processor: F,
    stage_name: &'static str,
    batch_target_bytes: usize,        // legacy fallback
) -> ReceiverStream<ParseOutput>
where
    F: Fn(Vec<u8>) -> Vec<Vec<u8>> + Send + Sync + 'static,
{
    let (tx, rx) = mpsc::channel(1024);
    let processor = Arc::new(processor);

    async fn send_lines(
        tx: &mpsc::Sender<ParseOutput>,
        lines: Vec<Vec<u8>>,
        stage_name: &'static str,
        label: &'static str,
    ) -> bool {
        for line in lines {
            if tx.send(ParseOutput::Bytes(Bytes::from(line))).await.is_err() {
                warn!("{}: downstream closed while sending {}", stage_name, label);
                return false;
            }
        }
        true
    }

    tokio::spawn(async move {
        let batch_target = crate::utils::system::compute_buffer_size(
            &config,
            "batch_rayon",
            StreamDataType::JustBytes,
            2.0,                          // hot CPU path — larger batches
        ).max(batch_target_bytes);        // respect legacy param if larger

        let mut buf = Vec::with_capacity(batch_target);
        let mut processed = 0u64;

        while let Some(item) = input.next().await {
            let ParseOutput::Bytes(b) = item else { continue; };

            buf.reserve(b.len());
            buf.extend_from_slice(&b);

            if buf.len() >= batch_target {
                let mut batch = Vec::with_capacity(buf.len());
                batch.append(&mut buf);

                let p = Arc::clone(&processor);
                let config_clone = Arc::clone(&config);

                let lines = match tokio::task::spawn_blocking(move || {
                    config_clone.thread_pool.install(|| p(batch))
                }).await {
                    Ok(lines) => lines,
                    Err(join_err) => {
                        error!("{}: batch task panicked: {}", stage_name, join_err);
                        break;
                    }
                };

                let batch_len = lines.len() as u64;
                if !send_lines(&tx, lines, stage_name, "batch").await {
                    break;
                }

                processed += batch_len;
            }
        }

        // Final flush
        if !buf.is_empty() {
            let mut batch = Vec::with_capacity(buf.len());
            batch.append(&mut buf);

            let p = Arc::clone(&processor);
            let config_clone = Arc::clone(&config);

            let lines = match tokio::task::spawn_blocking(move || {
                config_clone.thread_pool.install(|| p(batch))
            }).await {
                Ok(lines) => lines,
                Err(join_err) => {
                    error!("{}: final batch task panicked: {}", stage_name, join_err);
                    return;
                }
            };

            let batch_len = lines.len() as u64;
            let _ = send_lines(&tx, lines, stage_name, "final batch").await;
            processed += batch_len;
        }

        debug!("{}: finished — total {} items (target batch {} bytes)", stage_name, processed, batch_target);
    });

    ReceiverStream::new(rx)
}

/// Router that reads the upstream ONCE and fans it out to N private channels.
/// Guarantees immediate draining + instant EOF propagation. No silent drops.
/// Router that reads the upstream ONCE and fans it out to N private channels.
/// Guarantees immediate draining + instant EOF propagation. No silent drops.
pub async fn fanout_to_channels(
    mut upstream: ReceiverStream<ParseOutput>,
    n_branches: usize,
    label: &str,
    config: &RunConfig,           // ← added
    data_type: StreamDataType,    // ← added
) -> Result<(Vec<mpsc::Receiver<ParseOutput>>, JoinHandle<anyhow::Result<(), anyhow::Error>>), PipelineError> {

    let buffer_per_branch = crate::utils::system::compute_buffer_size(
        config,
        "fanout",                     // phase name
        data_type,
        1.45,                         // good hot-path boost for fanout (not too aggressive)
    );

    let mut txs = Vec::with_capacity(n_branches);
    let mut rxs = Vec::with_capacity(n_branches);

    for _ in 0..n_branches {
        let (tx, rx) = mpsc::channel(buffer_per_branch);
        txs.push(tx);
        rxs.push(rx);
    }

    let label = label.to_string();
    let router: JoinHandle<anyhow::Result<(), anyhow::Error>> = tokio::spawn(async move {
        let mut count = 0usize;
        while let Some(item) = upstream.next().await {
            count += 1;
            // Send to EVERY branch (never drop)
            let sends: Vec<_> = txs.iter().map(|tx| tx.send(item.clone())).collect();
            let _ = futures::future::join_all(sends).await;
        }
        debug!("[fanout {}] finished — forwarded {} items total", label, count);
        Ok(())
    });

    Ok((rxs, router))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;
    use std::io::Read;
    use log::{self, LevelFilter, debug, info, error, warn};
    use tokio::process::Command;
    use tokio::task;
    use tokio::fs::File as TokioFile;
    use tokio::time::{self, Duration};
    use crate::utils::fastx::fastx_generator;
    use crate::config::defs::{GpuDetection, RunConfig, StreamDataType, NRAlignmentBackend, GpuInfo};
    use std::sync::Arc;
    use rayon::ThreadPoolBuilder;
    use tokio::sync::{Semaphore, Mutex};
    use std::path::PathBuf;
    use crate::cli::Arguments;
    use tempfile::tempdir;
    use tokio::fs::metadata;
    use std::os::unix::fs::FileTypeExt;
    use crate::utils::fastx::SequenceRecord;
    use tokio::io::{AsyncReadExt, AsyncWriteExt};
    use tokio::io::duplex;


    // Helper function to create a RunConfig for tests
    // Helper function to create a RunConfig for tests
    fn create_test_run_config() -> Arc<RunConfig> {
        let args = Arguments {
            threads: 8,
            ..Default::default()
        };

        let (total_ram, available_ram) = detect_ram()
            .unwrap_or((16u64 << 30, 8u64 << 30));

        let rng = generate_rng(Some(42));

        let mut run_config = RunConfig {
            cwd: PathBuf::from("."),
            ram_temp_dir: std::env::temp_dir(),
            out_dir: PathBuf::from("test"),
            args,
            thread_pool: Arc::new(ThreadPoolBuilder::new().num_threads(8).build().unwrap()),
            maximal_semaphore: Arc::new(Semaphore::new(8)),
            base_buffer_size: 0,                    // temporary placeholder
            input_size: 100 + 1_048_576,
            physical_cores: 4,
            max_cores: 32,
            available_ram,
            rng,
            log_level: LevelFilter::Debug,
            base_backpressure_pause: 1000,
            simd: SimdLevel::Scalar,
            gpu_info: GpuDetection { count: 0, gpus: vec![] },
            has_gpu: false,
            alignment_backend: NRAlignmentBackend::Diamond,
        };

        // Compute proper buffer size exactly like main.rs
        let sdt = StreamDataType::IlluminaFastq; // default for tests

        let base_buffer_size = crate::utils::system::compute_buffer_size(
            &run_config,
            "test_global_default",
            sdt,
            1.0,
        );

        run_config.base_buffer_size = base_buffer_size;

        Arc::new(run_config)
    }


    #[tokio::test]
    async fn test_to_bytes_formats_fastq_and_fasta() -> Result<()> {
        let fastq = SequenceRecord::Fastq {
            id: "read1".to_string(),
            desc: Some("mate 1".to_string()),
            seq: Arc::new(b"ATCG".to_vec()),
            qual: Arc::new(b"IIII".to_vec()),
        };
        let fasta = SequenceRecord::Fasta {
            id: "contig1".to_string(),
            desc: Some("assembled contig".to_string()),
            seq: Arc::new(b"GATTACA".to_vec()),
        };

        assert_eq!(fastq.to_bytes()?.as_ref(), b"@read1 mate 1\nATCG\n+\nIIII\n");
        assert_eq!(fasta.to_bytes()?.as_ref(), b">contig1 assembled contig\nGATTACA\n");
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_fasta_multiline_and_description() -> Result<()> {
        let config = create_test_run_config();
        let (mut writer, reader) = duplex(64 * 1024);
        let input = b">seq1 first record\nACGT\nTGCA\n>seq2\nNNNN\n";

        let write_task = tokio::spawn(async move {
            writer.write_all(input).await?;
            writer.shutdown().await?;
            Ok::<_, anyhow::Error>(())
        });

        let rx = parse_fasta(reader, &config, StreamDataType::JustBytes).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut records = Vec::new();
        while let Some(item) = stream.next().await {
            records.push(item);
        }

        write_task.await??;
        assert_eq!(records.len(), 2);

        match &records[0] {
            ParseOutput::Fasta(SequenceRecord::Fasta { id, desc, seq }) => {
                assert_eq!(id, "seq1");
                assert_eq!(desc.as_deref(), Some("first record"));
                assert_eq!(seq.as_ref().as_slice(), b"ACGTTGCA");
            }
            other => panic!("Expected first FASTA record, got {:?}", other),
        }

        match &records[1] {
            ParseOutput::Fasta(SequenceRecord::Fasta { id, desc, seq }) => {
                assert_eq!(id, "seq2");
                assert_eq!(desc, &None);
                assert_eq!(seq.as_ref().as_slice(), b"NNNN");
            }
            other => panic!("Expected second FASTA record, got {:?}", other),
        }

        Ok(())
    }

    #[tokio::test]
    async fn test_parse_lines_skips_blanks_and_trims_endings() -> Result<()> {
        let config = create_test_run_config();
        let (mut writer, reader) = duplex(64 * 1024);
        let input = b"line1\n\nline2\r\nline3\n";

        let write_task = tokio::spawn(async move {
            writer.write_all(input).await?;
            writer.shutdown().await?;
            Ok::<_, anyhow::Error>(())
        });

        let rx = parse_lines(reader, &config, StreamDataType::JustBytes).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut lines = Vec::new();
        while let Some(ParseOutput::Bytes(bytes)) = stream.next().await {
            lines.push(String::from_utf8(bytes.to_vec())?);
        }

        write_task.await??;
        assert_eq!(lines, vec!["line1", "line2", "line3"]);
        Ok(())
    }

    #[tokio::test]
    async fn test_bytes_to_lines_reassembles_chunk_boundaries() -> Result<()> {
        let (tx, rx) = mpsc::channel(8);
        let (mut out_rx, task) = bytes_to_lines(rx, 8).await?;

        tx.send(ParseOutput::Bytes(Bytes::from_static(b"abc\nx"))).await?;
        tx.send(ParseOutput::Bytes(Bytes::from_static(b"yz\n\nlast"))).await?;
        drop(tx);

        let mut got = Vec::new();
        while let Some(item) = out_rx.recv().await {
            match item {
                ParseOutput::Bytes(bytes) => got.push(String::from_utf8(bytes.to_vec())?),
                other => panic!("Expected bytes, got {:?}", other),
            }
        }

        task.await??;
        assert_eq!(got, vec!["abc\n", "xyz\n", "\n", "last"]);
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_lines_skips_blank_lines() -> Result<()> {
        let config = create_test_run_config();
        let (mut writer, reader) = tokio::io::duplex(64);

        tokio::spawn(async move {
            use tokio::io::AsyncWriteExt;
            writer.write_all(b"abc\n\nxyz\n\nlast\n").await.unwrap();
        });

        let rx = parse_lines(reader, &config, StreamDataType::JustBytes).await?;
        let mut stream = ReceiverStream::new(rx);

        let mut got = Vec::new();
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(bytes) => got.push(String::from_utf8(bytes.to_vec())?),
                other => panic!("Expected bytes, got {:?}", other),
            }
        }

        assert_eq!(got, vec!["abc", "xyz", "last"]);
        Ok(())
    }

    #[tokio::test]
    async fn test_y_junction_concatenates_streams_in_order() -> Result<()> {
        let (tx1, rx1) = mpsc::channel(8);
        let (tx2, rx2) = mpsc::channel(8);

        tx1.send(ParseOutput::Bytes(Bytes::from_static(b"a1"))).await?;
        tx1.send(ParseOutput::Bytes(Bytes::from_static(b"a2"))).await?;
        drop(tx1);
        tx2.send(ParseOutput::Bytes(Bytes::from_static(b"b1"))).await?;
        drop(tx2);

        let (rx, task) = y_junction(vec![rx1, rx2], 8).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut got = Vec::new();
        while let Some(ParseOutput::Bytes(bytes)) = stream.next().await {
            got.push(String::from_utf8(bytes.to_vec())?);
        }

        task.await??;
        assert_eq!(got, vec!["a1", "a2", "b1"]);
        Ok(())
    }

    #[tokio::test]
    async fn test_convert_fasta_stream_to_sequence_record_filters_non_fasta() -> Result<()> {
        let (tx, rx) = mpsc::channel(8);
        tx.send(ParseOutput::Fasta(SequenceRecord::Fasta {
            id: "seq1".to_string(),
            desc: Some("desc".to_string()),
            seq: Arc::new(b"ACGT".to_vec()),
        })).await?;
        tx.send(ParseOutput::Bytes(Bytes::from_static(b"ignored"))).await?;
        drop(tx);

        let (mut out_rx, task) = convert_fasta_stream_to_sequence_record(rx, 8).await?;
        let mut records = Vec::new();
        while let Some(item) = out_rx.recv().await {
            records.push(item);
        }

        task.await??;
        assert_eq!(records.len(), 1);
        match &records[0] {
            SequenceRecord::Fasta { id, desc, seq } => {
                assert_eq!(id, "seq1");
                assert_eq!(desc.as_deref(), Some("desc"));
                assert_eq!(seq.as_ref().as_slice(), b"ACGT");
            }
            other => panic!("Expected FASTA record, got {:?}", other),
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_fanout_to_channels_duplicates_items_to_all_branches() -> Result<()> {
        let config = create_test_run_config();
        let (up_tx, up_rx) = mpsc::channel(8);
        tokio::spawn(async move {
            up_tx.send(ParseOutput::Bytes(Bytes::from_static(b"first"))).await.unwrap();
            up_tx.send(ParseOutput::Bytes(Bytes::from_static(b"second"))).await.unwrap();
            drop(up_tx);
        });
        let upstream = ReceiverStream::new(up_rx);

        let (receivers, router) = fanout_to_channels(
            upstream,
            3,
            "test_fanout",
            &config,
            StreamDataType::JustBytes,
        ).await?;

        let mut branch_data = Vec::new();
        for rx in receivers {
            let mut stream = ReceiverStream::new(rx);
            let mut items = Vec::new();
            while let Some(ParseOutput::Bytes(bytes)) = stream.next().await {
                items.push(String::from_utf8(bytes.to_vec())?);
            }
            branch_data.push(items);
        }

        router.await??;
        for branch in branch_data {
            assert_eq!(branch, vec!["first", "second"]);
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_parse_child_output_lines() -> Result<()> {
        let config = create_test_run_config();
        let mut cmd = Command::new("sh");
        cmd.arg("-c").arg("printf 'alpha\\n\\n beta\\n gamma\\n'");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let rx = parse_child_output(&mut child, ChildStream::Stdout, ParseMode::Lines, &config).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut got = Vec::new();
        while let Some(ParseOutput::Bytes(bytes)) = stream.next().await {
            got.push(String::from_utf8(bytes.to_vec())?);
        }
        assert_eq!(got, vec!["alpha", " beta", " gamma"]);
        Ok(())
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
                debug!("Consumer {} collected {} records", consumer_id, records.len());
                Ok::<_, anyhow::Error>(records)
            });
            handles.push(handle);
        }
        let all_records = time::timeout(Duration::from_secs(120), async {
            let mut all_records = Vec::with_capacity(2);
            for (i, handle) in handles.into_iter().enumerate() {
                let records = handle.await??;
                debug!("Task {} returned {} records", i, records.len());
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
            None
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
            None
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
                    error!("Failed to send ParseOutput");
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
            None
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
            None
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
            None
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
            None
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
            None
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
            None
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
        let config = create_test_run_config();
        let mut cmd = Command::new("echo");
        cmd.arg("@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nHHHH\n");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let rx = parse_child_output(&mut child, ChildStream::Stdout, ParseMode::Fastq, &config).await?;
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
        let config = create_test_run_config();
        let mut cmd = Command::new("echo");
        cmd.arg("test data");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let rx = parse_child_output(&mut child, ChildStream::Stdout, ParseMode::Bytes, &config).await?;
        let mut stream = ReceiverStream::new(rx);
        while let Some(ParseOutput::Bytes(chunk)) = stream.next().await {
            assert!(String::from_utf8_lossy(&*chunk).contains("test data"));
            break;
        }
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_file_fastq() -> Result<()> {
        let config = create_test_run_config();
        let _ = fs::remove_file("stream_to_file_test_illumina.fq");
        let num_records = 2;
        let mut records = fastx_generator(num_records, 150, 30.0, 8.0).map(ParseOutput::Fastq);
        let (tx, rx) = mpsc::channel(1024);

        tokio::spawn(async move {
            while let Some(record) = records.next().await {
                if tx.send(record).await.is_err() {
                    error!("Failed to send record");
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
        let rx = parse_fastq(file, &config, StreamDataType::IlluminaFastq).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut parsed_records = Vec::new();

        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
            parsed_records.push(record);
        }

        debug!("Parsed {} records from file", parsed_records.len());

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
        let config = create_test_run_config();
        let _ = fs::remove_file("stream_to_file_norecord_test.fq");
        let (tx, rx) = mpsc::channel(1024);
        let mut records = fastx_generator(0, 10, 35.0, 3.0).map(ParseOutput::Fastq);

        tokio::spawn(async move {
            while let Some(record) = records.next().await {
                if tx.send(record).await.is_err() {
                    error!("Failed to send record");
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
        let rx = parse_fastq(file, &config, StreamDataType::IlluminaFastq).await?;
        let mut stream = ReceiverStream::new(rx);
        let mut parsed_records = Vec::new();

        while let Some(ParseOutput::Fastq(record)) = stream.next().await {
            parsed_records.push(record);
        }

        debug!("Parsed {} records from file", parsed_records.len());

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
        let config = create_test_run_config();
        let mut cmd = Command::new("echo");
        cmd.arg("line1\nline2\n\nline3");
        let mut child = cmd.stdout(std::process::Stdio::piped()).spawn()?;
        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout, &config).await?;
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "line1");
        assert_eq!(lines[1], "line2");
        assert_eq!(lines[2], "line3");
        Ok(())
    }

    #[tokio::test]
    async fn test_read_child_output_to_vec_stderr() -> Result<()> {
        let config = create_test_run_config();
        let mut cmd = Command::new("sh");
        cmd.arg("-c").arg("echo error >&2");
        let mut child = cmd.stderr(std::process::Stdio::piped()).spawn()?;
        let lines = read_child_output_to_vec(&mut child, ChildStream::Stderr, &config).await?;
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "error");
        Ok(())
    }

    #[tokio::test]
    async fn test_spawn_cmd_stderr_verbose() -> Result<()> {
        let config = create_test_run_config();
        let args = vec!["-c".to_string(), "echo error >&2".to_string()];
        let (mut child, stderr_task) = spawn_cmd(config, "sh", args, true, None).await?;
        let status = child.wait().await?;
        stderr_task.await??;
        assert!(status.success(), "Child process should exit successfully");
        Ok(())
    }

    #[tokio::test]
    async fn test_spawn_cmd_stderr_non_verbose() -> Result<()> {
        let config = create_test_run_config();
        let args = vec!["-c".to_string(), "echo error >&2".to_string()];
        let (mut child, stderr_task) = spawn_cmd(config, "sh", args, false, None).await?;
        let status = child.wait().await?;
        stderr_task.await??;
        assert!(status.success(), "Child process should exit successfully");
        Ok(())
    }

    #[tokio::test]
    async fn test_stream_to_file_bytes_no_clone() -> Result<()> {
        let (tx, rx) = mpsc::channel(10);
        let data = Bytes::from(vec![b'A'; 1_000_000]);

        tokio::spawn(async move {
            tx.send(ParseOutput::Bytes(data)).await.unwrap();
        });

        let path = PathBuf::from("test_bytes.fq");
        let task = tokio::spawn(stream_to_file(rx, path.clone()));
        task.await??;

        let file = std::fs::File::open(&path)?;
        let mut contents = Vec::new();
        std::io::BufReader::new(file).read_to_end(&mut contents)?;

        assert_eq!(contents.len(), 1_000_000);
        assert!(contents.iter().all(|&b| b == b'A'));

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
            assert!(
                e.to_string().contains("is not a FIFO"),
                "Error should mention non-FIFO file, got: {}",
                e
            );
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
        let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_paired",
            true,
        )
            .await?;

        let r2_fifo = r2_fifo_opt.as_ref().expect("R2 FIFO should exist in paired mode");

        // Spawn readers for both FIFOs
        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle_opt {
            handle.await??;
        }

        let r1_str = String::from_utf8_lossy(&r1_output);
        let r2_str = String::from_utf8_lossy(&r2_output);
        assert!(r1_str.contains("@read1/1\nATCG\n+\nIIII\n"), "R1 should contain read1/1");
        assert!(r1_str.contains("@read2/1\nCCCC\n+\nJJJJ\n"), "R1 should contain read2/1");
        assert!(r2_str.contains("@read1/2\nGCTA\n+\nHHHH\n"), "R2 should contain read1/2");
        assert!(r2_str.contains("@read2/2\nGGGG\n+\nKKKK\n"), "R2 should contain read2/2");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(r2_fifo).await.ok();
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

        let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_single",
            false,
        )
            .await?;

        assert!(r2_fifo_opt.is_none(), "R2 FIFO should not be created in single-end mode");
        assert!(r2_write_handle_opt.is_none(), "R2 writer should not exist in single-end mode");

        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;

        let mut r1_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;

        let r1_str = String::from_utf8_lossy(&r1_output);
        assert!(r1_str.contains("@read1\nATCG\n+\nIIII\n"), "R1 should contain read1");
        assert!(r1_str.contains("@read2\nGCTA\n+\nHHHH\n"), "R1 should contain read2");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_empty_stream() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        drop(tx); // Close immediately
        let stream = ReceiverStream::new(rx);

        let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_empty",
            true,
        )
            .await?;

        let r2_fifo = r2_fifo_opt.as_ref().expect("R2 FIFO should exist");

        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle_opt {
            handle.await??;
        }

        assert!(r1_output.is_empty(), "R1 should be empty for empty stream");
        assert!(r2_output.is_empty(), "R2 should be empty for empty stream");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(r2_fifo).await.ok();
        Ok(())
    }

    #[tokio::test]
    async fn test_deinterleave_fastq_stream_to_fifos_non_fastq() -> Result<()> {
        let config = create_test_run_config();
        let (tx, rx) = mpsc::channel::<ParseOutput>(10);
        let stream = ReceiverStream::new(rx);

        tx.send(ParseOutput::Bytes(Bytes::from_static(b"invalid"))).await?;

        let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_non_fastq",
            true,
        )
            .await?;

        let _ = Command::new("cat").arg(&r1_fifo).spawn()?;
        let _ = Command::new("cat").arg(r2_fifo_opt.as_ref().unwrap()).spawn()?;

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
        if let Some(handle) = r2_write_handle_opt {
            handle.await??;
        }

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(r2_fifo_opt.unwrap()).await.ok();
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

        let (r1_fifo, r2_fifo_opt, deinterleave_handle, r1_write_handle, r2_write_handle_opt) = deinterleave_fastq_stream_to_fifos(
            config.clone(),
            stream,
            "test_ram_paired",
            true,
        )
            .await?;

        let r2_fifo = r2_fifo_opt.as_ref().expect("R2 FIFO should exist");

        // Verify FIFOs exist in ram_temp_dir
        assert!(
            metadata(&r1_fifo).await?.file_type().is_fifo(),
            "R1 FIFO should be in ram_temp_dir"
        );
        assert!(
            metadata(r2_fifo).await?.file_type().is_fifo(),
            "R2 FIFO should be in ram_temp_dir"
        );

        let mut r1_reader = Command::new("cat")
            .arg(&r1_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R1: {}", e))?;
        let mut r2_reader = Command::new("cat")
            .arg(r2_fifo)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn cat for R2: {}", e))?;

        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        r1_reader.stdout.take().unwrap().read_to_end(&mut r1_output).await?;
        r2_reader.stdout.take().unwrap().read_to_end(&mut r2_output).await?;

        deinterleave_handle.await??;
        r1_write_handle.await??;
        if let Some(handle) = r2_write_handle_opt {
            handle.await??;
        }

        let r1_str = String::from_utf8_lossy(&r1_output);
        let r2_str = String::from_utf8_lossy(&r2_output);
        assert!(r1_str.contains("@read1/1\nATCG\n+\nIIII\n"), "R1 should contain read1/1");
        assert!(r2_str.contains("@read1/2\nGCTA\n+\nHHHH\n"), "R2 should contain read1/2");

        tokio::fs::remove_file(&r1_fifo).await.ok();
        tokio::fs::remove_file(r2_fifo).await.ok();
        Ok(())
    }

}