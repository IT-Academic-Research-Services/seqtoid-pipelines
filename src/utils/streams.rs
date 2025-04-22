// src/utils/stream.rs
use std::io::{self, Write};
use std::path::PathBuf;
use anyhow::anyhow;
use tokio::sync::broadcast;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncReadExt, AsyncWriteExt, BufReader, BufWriter};
use tokio::process::{Child, Command};
use tokio::sync::mpsc;
use tokio::time::{sleep, Duration, Instant};
use tokio_stream::{Stream, StreamExt};
use tokio_stream::wrappers::BroadcastStream;
use crate::utils::file::WriteToFile;
use crate::utils::fastx::SequenceRecord;
use tokio::sync::broadcast::error::RecvError;
use tokio::process::ChildStdout;




// ToBytes and implementations (defined locally)
pub trait ToBytes {
    fn to_bytes(&self) -> std::io::Result<Vec<u8>>;
}


/// Trait to convert items to bytes for writing to a process's stdin


// Implementation for generic byte-like types
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

/// Implementation for SequenceRecord
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


/// Generates any number of output streams from a single input stream.
/// The streams are all asynchronous.
///
/// # Arguments
///
/// * `input_rx' - Receiver stream: tokio::mpsc
/// * 'num_streams' - Number of output streams to generate.
/// 
/// # Returns
/// Vector of output streams.
///
pub async fn t_junction<T>(
    input_rx: mpsc::Receiver<T>,
    num_streams: usize,
) -> Vec<BroadcastStream<T>>
where
    T: Clone + Send + 'static,
{
    let (tx, _rx) = broadcast::channel(100000);

    let mut streams = Vec::new();
    for _ in 0..num_streams {
        let rx = tx.subscribe();
        streams.push(BroadcastStream::new(rx));
    }

    tokio::spawn(async move {
        let mut rx = input_rx;
        let mut count = 0;
        let mut last_progress = tokio::time::Instant::now();

        while let Some(item) = rx.recv().await {
            count += 1;
            if last_progress.elapsed() > tokio::time::Duration::from_secs(15) {
                eprintln!("Broadcast stall detected at {} records", count);
                last_progress = tokio::time::Instant::now();
            }

            if tx.send(item).is_err() {
                eprintln!("Broadcast: No active receivers after {} records", count);
            }
            sleep(Duration::from_millis(1)).await;
        }
        eprintln!("Broadcast task completed, sent {} records", count);
    });

    streams
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
pub async fn stream_to_cmd<T>(
    mut input_stream: BroadcastStream<T>,
    cmd: &str,
    args: Vec<String>,
) -> anyhow::Result<Child>
where
    T: ToBytes + Clone + Send + 'static,
    BroadcastStream<T>: Stream<Item = Result<T, RecvError>>,
{
    let cmd = cmd.to_string(); // Clone to owned String
    let mut command = Command::new(&cmd);
    command
        .args(args.as_slice())
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());

    let mut child = command
        .spawn()
        .map_err(|e| anyhow!("Failed to spawn command {}: {}", cmd, e))?;

    let mut stdin = child
        .stdin
        .take()
        .ok_or_else(|| anyhow!("Failed to get stdin for command {}", cmd))?;

    let mut stderr = child
        .stderr
        .take()
        .ok_or_else(|| anyhow!("Failed to get stderr for command {}", cmd))?;

    {
        let cmd = cmd.clone();
        let args = args.clone(); // Clone for stderr task
        tokio::spawn(async move {
            let mut stderr_reader = tokio::io::BufReader::new(stderr);
            let mut line = String::new();
            while let Ok(n) = stderr_reader.read_line(&mut line).await {
                if n == 0 {
                    break;
                }
                eprintln!("{} {} error: {}", cmd, args.join(" "), line.trim());
                line.clear();
            }
        });
    }

    {
        let cmd = cmd.clone();
        let args = args.clone(); // Clone for stdin task
        tokio::spawn(async move {
            let mut stdin_writer = BufWriter::with_capacity(4 * 1024 * 1024, &mut stdin);
            let mut count = 0;
            let mut last_progress = tokio::time::Instant::now();

            while let Some(result) = input_stream.next().await {
                if last_progress.elapsed() > Duration::from_secs(15) {
                    eprintln!("{} {} stall detected at {} records", cmd, args.join(" "), count);
                    last_progress = tokio::time::Instant::now();
                }

                match result {
                    Ok(data) => {
                        let write_result = tokio::task::spawn_blocking(move || data.to_bytes())
                            .await;

                        match write_result {
                            Ok(Ok(bytes)) => {
                                if let Err(e) = stdin_writer.write_all(&bytes).await {
                                    eprintln!("Failed to write to {} {} stdin: {}", cmd, args.join(" "), e);
                                    break;
                                }
                                count += 1;
                            }
                            Ok(Err(e)) => {
                                eprintln!("Failed to convert data to bytes: {} {}: {}", cmd, args.join(" "), e);
                                break;
                            }
                            Err(e) => {
                                eprintln!("Spawn blocking error: {} {}: {}", cmd, args.join(" "), e);
                                break;
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("{} {} stream error after {} records: {:?}", cmd, args.join(" "), count, e);
                    }
                }
            }

            if let Err(e) = stdin_writer.flush().await {
                eprintln!("Failed to flush {} {} stdin: {}", cmd, args.join(" "), e);
            }

            eprintln!("{} {} completed, processed {} records", cmd, args.join(" "), count);
        });
    }

    Ok(child)
}


/// Writes a stream to a file.
/// Capture stdout and return from function.
///
/// # Arguments
///
/// * `input_stream' - tokio ChldStdout
/// * 'output_opatj' - PathBuf for output file

///
/// # Returns
/// Result
pub async fn stream_to_file(
    mut input_stream: ChildStdout,
    output_path: PathBuf,
) -> anyhow::Result<()> {
    let file = File::create(&output_path)
        .await // Add await
        .map_err(|e| anyhow!("Failed to create file {}: {}", output_path.display(), e))?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, &mut input_stream);

    let mut buffer = vec![0; 1024 * 1024]; // 1 MB read buffer
    let mut total_bytes = 0;
    let mut last_progress = Instant::now();

    loop {
        if last_progress.elapsed() > Duration::from_secs(15) {
            eprintln!("File writer stall detected at {} bytes", total_bytes);
            last_progress = Instant::now();
        }

        let bytes_read = reader
            .read(&mut buffer)
            .await
            .map_err(|e| anyhow!("Failed to read from stream: {}", e))?;

        if bytes_read == 0 {
            break;
        }

        writer
            .write_all(&buffer[..bytes_read])
            .await
            .map_err(|e| anyhow!("Failed to write to file {}: {}", output_path.display(), e))?;

        total_bytes += bytes_read;
    }

    writer
        .flush()
        .await
        .map_err(|e| anyhow!("Failed to flush file {}: {}", output_path.display(), e))?;

    eprintln!(
        "File stream for {} completed, wrote {} bytes",
        output_path.display(),
        total_bytes
    );

    Ok(())
}

/// Mostly for testing.
/// Reads child stdout to screen
///
/// # Arguments
///
/// * `child' - Child process from stream_to_cmd.
///
/// # Returns
/// io::Result<()>
pub async fn read_child_stdout(mut child: Child) -> io::Result<()> {
    let mut stdout = child
        .stdout
        .take()
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to open stdout"))?;

    let mut buffer = [0; 1024]; // Adjust buffer size as needed
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
