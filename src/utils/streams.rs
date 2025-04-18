// src/utils/stream.rs
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::PathBuf;
use tokio::sync::mpsc;
use tokio::sync::broadcast;
use tokio::process::{Command, Child};
use tokio::io::AsyncWriteExt;
use tokio_stream::StreamExt;
use tokio_stream::wrappers::BroadcastStream;
use tokio_stream::wrappers::errors::BroadcastStreamRecvError;
use crate::utils::fastx::SequenceRecord;
use crate::utils::file::{write_fastq_record, write_fasta_record};


/// Trait to abstract writing items to a file
pub trait WriteToFile {
    fn write_to_file(&self, file: &mut File) -> io::Result<()>;
}

/// Implementation for generic byte-like types
impl<T> WriteToFile for T
where
    T: AsRef<[u8]>,
{
    fn write_to_file(&self, file: &mut File) -> io::Result<()> {
        file.write_all(self.as_ref())?;
        file.write_all(b"\n")?;
        Ok(())
    }
}

/// Implementation for SequenceRecord
impl WriteToFile for SequenceRecord {
    fn write_to_file(&self, file: &mut File) -> io::Result<()> {
        match self {
            SequenceRecord::Fastq { id, desc, seq, qual } => {
                write_fastq_record(file, id, desc.as_deref(), seq, qual)
            }
            SequenceRecord::Fasta { id, desc, seq } => {
                write_fasta_record(file, id, desc.as_deref(), seq)
            }
        }
    }
}


/// Trait to convert items to bytes for writing to a process's stdin
pub trait ToBytes {
    fn to_bytes(&self) -> io::Result<Vec<u8>>;
}

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
/// * 'out_filename' - output file.
///
/// # Returns
/// One child stream for continuation of pipeline.
pub async fn tee<T>(
    input_rx: mpsc::Receiver<T>,
    out_path: PathBuf,
) -> BroadcastStream<T>
where
    T: WriteToFile + Clone + Send + 'static,
{
    let out_streams = t_junction(input_rx, 2).await;

    let mut streams = out_streams.into_iter();
    let return_stream = streams.next().unwrap();
    let mut file_stream = streams.next().unwrap();

    tokio::spawn(async move {
        let mut file = match File::create(&out_path) {
            Ok(file) => file,
            Err(e) => {
                eprintln!("Failed to create file {}: {}", out_path.display(), e);
                return;
            }
        };

        while let Some(result) = file_stream.next().await {
            match result {
                Ok(data) => {
                    if let Err(e) = data.write_to_file(&mut file) {
                        eprintln!("Failed to write data to {}: {}", out_path.display(), e);
                        break;
                    }
                }
                Err(BroadcastStreamRecvError::Lagged(skipped)) => {
                    eprintln!("File stream lagged, skipped {} items", skipped);
                }
            }
        }
        eprintln!("File stream for {} completed", out_path.display());
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
    T: ToBytes + Clone + Send + 'static,
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
                    match data.to_bytes() {
                        Ok(bytes) => {
                            if let Err(e) = stdin.write_all(&bytes).await {
                                eprintln!("Failed to write to stdin: {}", e);
                                break;
                            }
                        }
                        Err(e) => {
                            eprintln!("Failed to convert data to bytes: {}", e);
                            break;
                        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::NamedTempFile;
    use tokio::time::{timeout, Duration};
    
    #[tokio::test]
    async fn test_t_junction_splits_stream() {
        let (tx, rx) = mpsc::channel(10);
        let items = vec![b"item1".to_vec(), b"item2".to_vec()];
        
        let items_clone = items.clone();
        tokio::spawn(async move {
            for item in items_clone {
                tx.send(item).await.unwrap();
            }
        });
        
        let mut streams = t_junction(rx, 2).await;
        let mut stream1 = streams.pop().unwrap();
        let mut stream2 = streams.pop().unwrap();
        
        let mut results1 = Vec::new();
        let mut results2 = Vec::new();

        while let Some(result) = timeout(Duration::from_millis(100), stream1.next()).await.unwrap() {
            if let Ok(item) = result {
                results1.push(item);
            }
        }
        while let Some(result) = timeout(Duration::from_millis(100), stream2.next()).await.unwrap() {
            if let Ok(item) = result {
                results2.push(item);
            }
        }
        
        assert_eq!(results1, items);
        assert_eq!(results2, items);
    }
    
    #[tokio::test]
    async fn test_tee_writes_bytes_to_file() {
        let (tx, rx) = mpsc::channel(10);
        let items = vec![b"hello".to_vec(), b"world".to_vec()];
        let temp_file = NamedTempFile::new().unwrap();
        let path = PathBuf::from(temp_file.path());
        
        
        let items_clone = items.clone();
        tokio::spawn(async move {
            for item in items_clone {
                tx.send(item).await.unwrap();
            }
        });
        
        let mut stream = tee(rx, path.clone()).await;
        
        let mut results = Vec::new();
        while let Some(result) = timeout(Duration::from_millis(100), stream.next()).await.unwrap() {
            if let Ok(item) = result {
                results.push(item);
            }
        }
        
        assert_eq!(results, items);
        tokio::time::sleep(Duration::from_millis(100)).await;
        let contents = fs::read_to_string(&path).unwrap();
        assert_eq!(contents, "hello\nworld\n");
    }
    
    #[tokio::test]
    async fn test_tee_writes_fastq_to_file() {
        let (tx, rx) = mpsc::channel(10);
        let record = SequenceRecord::Fastq {
            id: "seq1".to_string(),
            desc: Some("test".to_string()),
            seq: b"ATCG".to_vec(),
            qual: b"IIII".to_vec(),
        };
        let temp_file = NamedTempFile::new().unwrap();
        let path = PathBuf::from(temp_file.path());
        
        let record_clone = record.clone();
        tokio::spawn(async move {
            tx.send(record_clone).await.unwrap();
        });
        
        let mut stream = tee(rx, path.clone()).await;
        
        let mut results = Vec::new();
        while let Some(result) = timeout(Duration::from_millis(100), stream.next()).await.unwrap() {
            if let Ok(item) = result {
                results.push(item);
            }
        }
        
        assert_eq!(results.len(), 1);
        if let SequenceRecord::Fastq { id, desc, seq, qual } = &results[0] {
            assert_eq!(id, "seq1");
            assert_eq!(desc, &Some("test".to_string()));
            assert_eq!(seq, &b"ATCG".to_vec());
            assert_eq!(qual, &b"IIII".to_vec());
        } else {
            panic!("Expected Fastq record");
        }
        
        tokio::time::sleep(Duration::from_millis(100)).await;
        
        let contents = fs::read_to_string(&path).unwrap();
        assert_eq!(contents, "@seq1 test\nATCG\n+\nIIII\n");
    }
    
    #[tokio::test]
    async fn test_stream_to_cmd_feeds_stdin_bytes() {
        let (tx, rx) = mpsc::channel(10);
        let items = vec![b"hello".to_vec(), b"world".to_vec()];
        let stream = t_junction(rx, 1).await.pop().unwrap();
        
        let items_clone = items.clone();
        tokio::spawn(async move {
            for item in items_clone {
                tx.send(item).await.unwrap();
            }
        });
        
        let child = stream_to_cmd(stream, "cat", &[]).await.unwrap();
        let output = child.wait_with_output().await.unwrap();
        
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert_eq!(stdout, "hello\nworld\n");
    }
    
    #[tokio::test]
    async fn test_stream_to_cmd_feeds_stdin_fastq() {
        let (tx, rx) = mpsc::channel(10);
        let record = SequenceRecord::Fastq {
            id: "seq1".to_string(),
            desc: Some("test".to_string()),
            seq: b"ATCG".to_vec(),
            qual: b"IIII".to_vec(),
        };
        let stream = t_junction(rx, 1).await.pop().unwrap();
        
        let record_clone = record.clone();
        tokio::spawn(async move {
            tx.send(record_clone).await.unwrap();
        });
        
        let child = stream_to_cmd(stream, "cat", &[]).await.unwrap();
        let output = child.wait_with_output().await.unwrap();
        
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert_eq!(stdout, "@seq1 test\nATCG\n+\nIIII\n");
    }
}