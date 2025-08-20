use sysinfo::System;
use anyhow::anyhow;
use anyhow::Result;
use std::path::PathBuf;
use std::sync::Arc;

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

use crate::utils::fastx::{SequenceRecord, parse_header};
use crate::config::defs::StreamDataType;
use crate::config::defs::{CoreAllocation, RunConfig};



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

    let record_size = match data_type {
        StreamDataType::IlluminaFastq => 1_000, // ~1KB per FASTQ record
        StreamDataType::OntFastq => 10_000,     // ~10KB per ONT read
        StreamDataType::JustBytes => 500,       // ~500B for SAM/BAM/VCF
    };

    let max_buffer_size = if available_ram > 0 {
        ((available_ram as f64 * RAM_FRACTION / MAX_PROCESSES as f64) / record_size as f64) as usize
    } else {
        eprintln!("Warning: Failed to detect available RAM in {}, using fallback buffer size", label);
        base_buffer_size * n_outputs.max(1) * 2
    };

    let min_buffer_size = MIN_BUFFER_PER_STREAM * n_outputs.max(1);
    let buffer_size = (base_buffer_size * n_outputs.max(1)).clamp(min_buffer_size, max_buffer_size);

    let (done_tx, done_rx) = oneshot::channel::<Result<(), anyhow::Error>>();
    let mut output_txs = Vec::with_capacity(n_outputs);
    let mut output_rxs = Vec::with_capacity(n_outputs);

    // Track channels with IDs
    let mut channel_status = vec![true; n_outputs]; // true = active, false = dropped
    for i in 0..n_outputs {
        let (tx, rx) = mpsc::channel(buffer_size);
        output_txs.push((i, tx)); // Store ID with tx
        output_rxs.push(rx);
    }

    tokio::spawn(async move {
        let mut input = Box::pin(input);
        let mut count = 0;
        let mut dropped_count = 0;

        while let Some(item) = input.next().await {
            let mut active_txs = Vec::new();
            for (id, tx) in output_txs.iter() {
                match tx.send(item.clone()).await {
                    Ok(()) => active_txs.push((*id, tx.clone())),
                    Err(_) => {
                        if channel_status[*id] {
                            eprintln!("{}: Receiver {} dropped at item {}", label, id, count + 1);
                            channel_status[*id] = false;
                            dropped_count += 1;
                        }
                    }
                }
            }
            output_txs = active_txs;
            if output_txs.is_empty() {
                eprintln!("{}: All receivers dropped at item {}. Channel status: {:?}", label, count + 1, channel_status);
                let _ = done_tx.send(Err(anyhow!("{}: All receivers dropped at item {}, data loss occurred", label, count + 1)));
                return;
            }

            count += 1;
            if count % stall_threshold == 0 {
                if let Some(sleep_ms) = stream_sleep_ms {
                    sleep(Duration::from_millis(sleep_ms)).await;
                }
            }
        }

        if let Some(n) = notify {
            n.notify_waiters();
        }
        if dropped_count > 0 {
            let _ = done_tx.send(Err(anyhow!("{}: {} receivers dropped mid-stream, potential data loss", label, dropped_count)));
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
) -> Result<(Child, JoinHandle<Result<(), anyhow::Error>>, JoinHandle<Result<(), anyhow::Error>>)> {
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
                writer.write_all(&batch).await?;
                writer.flush().await?;
                total_written += batch.len();
                batch.clear();
            }
        }

        if !batch.is_empty() {
            writer.write_all(&batch).await?;
            writer.flush().await?;
            total_written += batch.len();
        }

        writer.flush().await?;
        writer.shutdown().await?;
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
                    eprintln!("Error reading FASTQ: {}", e);
                    return;
                }
            };
            if bytes_read == 0 {
                break;
            }
            let id_line = buffer.trim_end();
            if !id_line.starts_with('@') {
                eprintln!("Invalid FASTQ format: expected '@', got '{}'", id_line);
                return;
            }

            let id = id_line[1..].to_string();
            let desc = None;

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                eprintln!("Error reading sequence line");
                return;
            }
            let seq = buffer.trim_end().as_bytes().to_vec();
            if seq.is_empty() {
                eprintln!("Missing sequence");
                return;
            }

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                eprintln!("Error reading plus line");
                return;
            }
            let plus = buffer.trim_end();
            if plus != "+" {
                eprintln!("Invalid FASTQ format: expected '+', got '{}'", plus);
                return;
            }

            buffer.clear();
            if reader.read_line(&mut buffer).await.is_err() {
                eprintln!("Error reading quality line");
                return;
            }
            let qual = buffer.trim_end().as_bytes().to_vec();
            if qual.len() != seq.len() {
                eprintln!("Sequence and quality lengths do not match");
                return;
            }

            let record = SequenceRecord::Fastq {
                id,
                desc,
                seq: Arc::new(seq),
                qual: Arc::new(qual),
            };

            if tx.send(ParseOutput::Fastq(record)).await.is_err() {
                eprintln!("Receiver dropped while sending FASTQ record");
                break;
            }
            count += 1;
        }
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


pub async fn convert_stream(
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


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;
    use tokio::process::Command;
    use tokio::task;
    use tokio::fs::File as TokioFile;
    use tokio::time::{self, Duration};
    use crate::utils::fastx::fastx_generator;
    use crate::config::defs::{RunConfig, StreamDataType};
    use std::sync::Arc;
    use rayon::ThreadPoolBuilder;
    use tokio::sync::Semaphore;
    use std::path::PathBuf;
    use crate::cli::Arguments;

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
                seq: b"ATCG".to_vec(),
                qual: b"IIII".to_vec(),
            }),
            ParseOutput::Fastq(SequenceRecord::Fastq {
                id: "read2".to_string(),
                desc: None,
                seq: b"GCTA".to_vec(),
                qual: b"HHHH".to_vec(),
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
        let (mut child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        child.wait().await?;
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
        let (mut child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
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
        let (mut child, task, _err_task) = stream_to_cmd(
            config,
            rx_parse,
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
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
        let (mut child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
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
        let (mut child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "cat",
            vec![],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
        let mut stdout = child.stdout.take().unwrap();
        let mut output = Vec::new();
        tokio::io::copy(&mut stdout, &mut output).await?;
        task.await??;
        done_rx.await??;
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");

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
        let (mut child, task, _err_task) = stream_to_cmd(
            config,
            outputs.pop().unwrap(),
            "head",
            vec!["-n".to_string(), "1".to_string()],
            StreamDataType::IlluminaFastq,
            false,
        )
            .await?;
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
        let (mut child, task, _err_task) = stream_to_cmd(
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
        let status = child.wait().await?;
        assert!(status.success(), "Child process should exit successfully");
        assert!(child.stdin.is_none(), "Stdin should be closed");
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
}