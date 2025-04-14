// src/utils/stream.rs
use std::fs::File;
use std::io;
use std::io::Write;
use tokio::sync::mpsc;
use tokio::sync::broadcast;
use tokio::process::{Command, Child};
use tokio::io::AsyncWriteExt;
use tokio_stream::StreamExt;
use tokio_stream::wrappers::BroadcastStream;
use tokio_stream::wrappers::errors::BroadcastStreamRecvError;
use crate::utils::fastx::SequenceRecord;
use crate::utils::file::{write_fastq_record, write_fasta_record};


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
    let (tx, _rx) = broadcast::channel(100);

    // Create N broadcast receivers first
    let mut streams = Vec::new();
    for _ in 0..num_streams {
        let rx = tx.subscribe();
        streams.push(BroadcastStream::new(rx));
    }

    // Spawn the task after creating receivers
    tokio::spawn(async move {
        let mut rx = input_rx;
        while let Some(item) = rx.recv().await {
            if tx.send(item).is_err() {
                break; // No active receivers
            }
        }
    });

    streams
}




/// Takes a generic stream from tokio_stream and asynchronously writes to a file.
///
///
/// # Arguments
///
/// * `input_rx' - Receiver stream: tokio::mpsc
/// * 'num_streams' - Number of output streams to generate.
///
/// # Returns
/// Output stream.
pub async fn tee_file <T>(
    input_rx: mpsc::Receiver<T>,
    out_filename: String,
) -> BroadcastStream<T>
where
    T: Clone + Send + 'static + AsRef<[u8]>,
{
    let out_streams = t_junction(input_rx, 2).await;

    let mut streams = out_streams.into_iter();
    let return_stream = streams.next().unwrap();
    let mut file_stream = streams.next().unwrap();

    tokio::spawn(async move {
        let mut file = match File::create(&out_filename) {
            Ok(file) => file,
            Err(e) => {
                eprintln!("Failed to create file {}: {}", out_filename, e);
                return;
            }
        };

        // let mut stream = file_stream;
        while let Some(result) = file_stream.next().await {
            match result {
                Ok(data) => {
                    if let Err(e) = file.write_all(data.as_ref()) {
                        eprintln!("Failed to write data: {}", e);
                        break;
                    }
                    if let Err(e) = file.write_all(b"\n") {
                        eprintln!("Failed to write newline: {}", e);
                        break;
                    }
                }
                Err(err) => match err {
                    BroadcastStreamRecvError::Lagged(skipped) => {
                        eprintln!("File stream lagged, skipped {} items", skipped);
                    }
                },
            }
            
        }
    });

    return_stream
}

/// Takes a stream of FASTX records from seq_io and asynchronously writes to a file.
///
///
/// # Arguments
///
/// * `input_rx' - Receiver stream: tokio::mpsc
/// * 'num_streams' - Number of output streams to generate.
///
/// # Returns
/// Output stream.
pub async fn tee_fastx(
    input_rx: mpsc::Receiver<SequenceRecord>,
    out_filename: String,
) -> BroadcastStream<SequenceRecord>
{
    let out_streams = t_junction(input_rx, 2).await;

    let mut streams = out_streams.into_iter();
    let return_stream = streams.next().unwrap();
    let file_stream = streams.next().unwrap();

    tokio::spawn(async move {
        let mut file = match File::create(&out_filename) {
            Ok(file) => file,
            Err(e) => {
                eprintln!("Failed to create file {}: {}", out_filename, e);
                return;
            }
        };

        let mut stream = file_stream;
        while let Some(result) = stream.next().await {
            match result {
                Ok(record) => {
                    match record {
                        SequenceRecord::Fastq { id, desc, seq, qual } => {
                            if let Err(e) = write_fastq_record(&mut file, &id, desc.as_deref(), &seq, &qual) {
                                eprintln!("Failed to write FASTQ record: {}", e);
                                break;
                            }
                        }
                        SequenceRecord::Fasta { id, desc, seq } => {
                            if let Err(e) = write_fasta_record(&mut file, &id, desc.as_deref(), &seq) {
                                eprintln!("Failed to write FASTA record: {}", e);
                                break;
                            }
                        }
                    }
                }
                Err(err) => match err {
                    BroadcastStreamRecvError::Lagged(skipped) => {
                        eprintln!("File stream lagged, skipped {} items", skipped);
                    }
                },
            }
        }
    });

    return_stream
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
    mut stream: BroadcastStream<T>,
    command: &str,
    args: &[&str],
) -> io::Result<Child>
where
    T: Clone + AsRef<[u8]> + Send + 'static,
{
    let mut child = Command::new(command)
        .args(args)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped()) // Optional
        .stderr(std::process::Stdio::piped()) // Optional
        .spawn()?;
    
    let mut stdin = child.stdin.take().ok_or_else(|| {
        io::Error::new(io::ErrorKind::Other, "Failed to open stdin")
    })?;
    
    tokio::spawn(async move {
        while let Some(result) = stream.next().await {
            match result {
                Ok(data) => {
                    if let Err(e) = stdin.write_all(data.as_ref()).await {
                        eprintln!("Failed to write to stdin: {}", e);
                        break;
                    }
                    if let Err(e) = stdin.write_all(b"\n").await {
                        eprintln!("Failed to write newline: {}", e);
                        break;
                    }
                }
                Err(BroadcastStreamRecvError::Lagged(skipped)) => {
                    eprintln!("Stream lagged, skipped {} items", skipped);
                }
            }
        }
    });

    Ok(child)
}
