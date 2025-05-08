use std::fs;
use anyhow::anyhow;
use seqtoid_pipelines::utils::fastx::fastx_generator;
use anyhow::Result;
use std::io::{stderr, Write};
use std::path::Path;
use std::time::{Instant};
use futures::StreamExt;
use sysinfo::{System, Pid, ProcessesToUpdate};
use seqtoid_pipelines::utils::streams::{read_child_output_to_vec, stream_to_cmd, t_junction, parse_child_output, stream_to_file, ChildStream, ParseMode};
use tokio::task::JoinHandle;
use tokio::time::timeout;
use tokio::process::Command;
use tokio::sync::{broadcast};
use tokio_stream::wrappers::BroadcastStream;
use std::sync::{Arc, Mutex};
use tokio::time::{Duration};
use futures::future::join_all;

#[tokio::test]
async fn test_fastx_generator_stress() -> Result<()> {

    let num_reads = vec![10_000, 100_000]; 
    let read_sizes = vec![100, 150, 1000, 5000];
    
    let mut log = std::fs::File::create("fastx_stress.log")?;
    writeln!(
        &mut log,
        "Reads\tSize\tRecords\tTime\tMemory"
    )?;
    log.flush()?;
    
    let mut sys = System::new_all();
    for read_size in &read_sizes {
        for num_read in &num_reads {
            eprintln!(
                "Testing: Reads: {}, Size: {}",
                num_read, read_size
            );
            stderr().flush()?;

            let start = Instant::now();
            
            let memory_before = {
                sys.refresh_memory();
                sys.used_memory()
            };

            // Generate stream and count records
            let stream = fastx_generator(*num_read, *read_size, 35.0, 3.0);
            let record_count = stream
                .fold(0usize, |acc, _| async move { acc + 1 })
                .await;

            let elapsed = start.elapsed();
            let elapsed_secs = elapsed.as_secs_f64();
            let memory_used = {
                    sys.refresh_memory();
                    let memory_after = sys.used_memory();
                    if memory_after >= memory_before {
                        (memory_after - memory_before) / 1024 / 1024 // MB
                    } else {
                        0 // Handle underflow
                    }
                
            };
            
            writeln!(
                &mut log,
                "{}\t{}\t{}\t{}\t{}",
                num_read, read_size, record_count, elapsed_secs, memory_used
            )?;
            log.flush()?;
        }
    }

    Ok(())
}

#[tokio::test]
async fn test_fastx_generator_count() -> Result<()> {
    let num_reads = vec![10_000];
    let read_sizes = vec![50, 100, 1000];
    let mut log = std::fs::File::create("fastx_count.log")?;
    for num_read in num_reads {
        for read_size in &read_sizes {
            writeln!(
                &mut log,"Testing fastx_generator: Reads: {}, Size: {}", num_read, read_size)?;
            stderr().flush()?;
            let stream = fastx_generator(num_read, *read_size, 35.0, 3.0);
            let count = stream.fold(0usize, |acc, _| async move { acc + 1 }).await;
            if count != num_read {
                writeln!(
                    &mut log,"fastx_generator failed: expected {}, got {}", num_read, count)?;
                return Err(anyhow::anyhow!("fastx_generator produced {} records, expected {}", count, num_read));
            }
            writeln!(
                &mut log,"fastx_generator produced {} records", count)?;
            stderr().flush()?;
        }
    }
    Ok(())
}

#[tokio::test]
async fn test_fastx_generator_edge_cases() -> Result<()> {
    // Test zero reads
    let stream = fastx_generator(0, 100, 35.0, 3.0);
    let count = stream.fold(0usize, |acc, _| async move { acc + 1 }).await;
    assert_eq!(count, 0, "Zero reads should produce no records");

    // Test single read
    let stream = fastx_generator(1, 100, 35.0, 3.0);
    let count = stream.fold(0usize, |acc, _| async move { acc + 1 }).await;
    assert_eq!(count, 1, "Single read should produce one record");

    // Test zero read size (if fastx_generator allows it)
    let stream = fastx_generator(1000, 0, 35.0, 3.0);
    let count = stream.fold(0usize, |acc, _| async move { acc + 1 }).await;
    assert_eq!(count, 0, "Zero read size should produce no records");

    Ok(())
}


#[tokio::test]
async fn test_t_junction_stress() -> Result<()> {
    let buffer_sizes = [1000, 10_000, 100_000];
    let stall_thresholds = [100, 1000];
    let sleep_ms_options = [Some(0), Some(10)];
    let backpressure_pause_ms_options = [50, 500];
    let num_reads = 100_000;
    let seq_len = 100;
    let n_outputs = 2;
    let mut log = std::fs::File::create("t_junction_stress.log")?;
    writeln!(
        &mut log,"BufferSize\tStall\tSleep\tBackpressurePause\tStreams\tReads\tSeqLen\tTime\tMemory\tRecords\tSuccess")?;

    for &buffer_size in &buffer_sizes {
        for &stall_threshold in &stall_thresholds {
            for &sleep_ms in &sleep_ms_options {
                for &backpressure_pause_ms in &backpressure_pause_ms_options {
                    eprintln!(
                        "Buffer size: {}  Stall: {}  Sleep: {}  Pause: {}  Streams: {} Reads: {}  Size: {}",
                        buffer_size,
                        stall_threshold,
                        sleep_ms.unwrap_or(0),
                        backpressure_pause_ms,
                        n_outputs,
                        num_reads,
                        seq_len
                    );

                    let stream = fastx_generator(num_reads, seq_len, 30.0, 5.0);
                    let start = Instant::now();
                    let (outputs, done_rx) = t_junction(
                        stream,
                        n_outputs,
                        buffer_size,
                        stall_threshold,
                        sleep_ms,
                        backpressure_pause_ms,
                    )
                        .await?;

                    let mut run_success = true;
                    let record_counts = Arc::new(Mutex::new(vec![0; n_outputs]));
                    for (i, mut output) in outputs.into_iter().enumerate() {
                        let record_counts = Arc::clone(&record_counts);
                        tokio::spawn(async move {
                            while let Some(_item) = output.next().await {
                                {
                                    let mut counts = record_counts.lock().unwrap();
                                    counts[i] += 1;
                                    if counts[i] % 1000 == 0 {
                                        eprintln!("Stream {} processed {} items", i, counts[i]);
                                    }
                                } // Drop MutexGuard here

                            }
                        });
                    }

                    match done_rx.await {
                        Ok(result) => match result {
                            Ok(()) => eprintln!("t_junction completed successfully"),
                            Err(e) => {
                                eprintln!("t_junction failed: {}", e);
                                run_success = false;
                            }
                        },
                        Err(e) => {
                            eprintln!("t_junction task failed to send: {}", e);
                            run_success = false;
                        }
                    }

                    let duration = start.elapsed();
                    let memory = 0; // Update with sysinfo if enabled
                    let record_counts = record_counts.lock().unwrap();
                    eprintln!(
                        "Records: {:?}  Success: {}",
                        *record_counts, run_success
                    );
                    writeln!(
                        &mut log,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}",
                        buffer_size,
                        stall_threshold,
                        sleep_ms.unwrap_or(0),
                        backpressure_pause_ms,
                        n_outputs,
                        num_reads,
                        seq_len,
                        duration.as_secs_f64(),
                        memory,
                        *record_counts,
                        run_success
                    )?;

                    assert_eq!(
                        *record_counts,
                        vec![num_reads; n_outputs],
                        "Incorrect record counts"
                    );

                }
            }
        }
    }
    Ok(())
}


#[tokio::test]
async fn test_t_junction_count() -> Result<()> {
    let num_read = 10_000;
    let read_size = 50;
    let buffer_size = 100_000;
    let stall_threshold = 100;
    let sleep_ms = Some(1);
    let backpressure_pause_ms = 500;
    eprintln!("Testing t_junction: Reads: {}, Size: {}, Buffer: {}, Sleep: {:?}", num_read, read_size, buffer_size, sleep_ms);
    stderr().flush()?;
    let stream = fastx_generator(num_read, read_size, 35.0, 3.0);
    let (mut outputs, done_rx) = t_junction(stream, 2, buffer_size, stall_threshold, sleep_ms, backpressure_pause_ms).await?;
    let mut counts = vec![0usize; 2];
    for (i, mut rx) in outputs.into_iter().enumerate() {
        while let Some(result) = rx.next().await {
            match result {
                Ok(_) => {
                    counts[i] += 1;
                    if counts[i] % 1000 == 0 {
                        eprintln!("t_junction stream {} counted {} records", i, counts[i]);
                        stderr().flush()?;
                    }
                }
                Err(e) => {
                    eprintln!("t_junction stream {} error: {}", i, e);
                    stderr().flush()?;
                }
            }
        }
        eprintln!("t_junction stream {} finished: {} records", i, counts[i]);
        stderr().flush()?;
    }
    match done_rx.await {
        Ok(Ok(())) => eprintln!("t_junction completed successfully"),
        Ok(Err(e)) => eprintln!("t_junction failed: {}", e),
        Err(e) => eprintln!("t_junction send error: {}", e),
    }
    stderr().flush()?;
    eprintln!("t_junction counts: {:?}", counts);
    stderr().flush()?;
    if counts != vec![num_read, num_read] {
        return Err(anyhow!("t_junction produced {:?}, expected [{}, {}]", counts, num_read, num_read));
    }
    Ok(())
}

#[tokio::test]
async fn test_stream_to_cmd_direct() -> Result<()> {
    let num_read = 100;
    let read_size = 50;
    let cmd_tag = "cat";
    let args = vec!["-"];
    eprintln!("Testing stream_to_cmd: Reads: {}, Size: {}, Command: {}", num_read, read_size, cmd_tag);
    stderr().flush()?;

    let stream = fastx_generator(num_read, read_size, 35.0, 3.0);
    let (tx, rx) = broadcast::channel(100_000);
    let rx = BroadcastStream::new(rx);

    tokio::spawn(async move {
        let mut count = 0;
        let mut stream = Box::pin(stream);
        while let Some(item) = stream.next().await {
            count += 1;
            if count % 10 == 0 {
                eprintln!("Sent {} records to broadcast channel", count);
                stderr().flush().ok();
            }
            if tx.send(item).is_err() {
                eprintln!("Broadcast channel dropped");
                break;
            }
        }
        eprintln!("Finished sending {} records", count);
        stderr().flush().ok();
    });

    let (child, inner_task) = stream_to_cmd(rx, cmd_tag, args).await?;
    match timeout(Duration::from_secs(30), inner_task).await {
        Ok(Ok(Ok(()))) => eprintln!("stream_to_cmd completed successfully"),
        Ok(Ok(Err(e))) => {
            eprintln!("stream_to_cmd failed: {}", e);
            return Err(e);
        }
        Ok(Err(e)) => {
            eprintln!("stream_to_cmd join error: {}", e);
            return Err(anyhow!("Join error: {}", e));
        }
        Err(_) => {
            eprintln!("stream_to_cmd timed out after 30s");
            return Err(anyhow!("stream_to_cmd timed out"));
        }
    }
    let output = child.wait_with_output().await?;
    let stderr_output = String::from_utf8_lossy(&output.stderr);
    if !output.status.success() {
        eprintln!("Child process failed: status: {}, stderr: {}", output.status, stderr_output);
        return Err(anyhow!("Child process failed: status: {}, stderr: {}", output.status, stderr_output));
    }
    eprintln!("Child process completed successfully");
    stderr().flush()?;
    Ok(())
}


/// Create multiple streams of fastq data with fatsx_generator
/// The number of records is thus known
/// Spawn N streams and cat the data into files
/// Read each file, and work out if it has the correct number of records
#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn test_stream_to_cmd_stress() -> Result<()> {
    let num_reads = vec![100, 10_000, 100_000];
    let read_sizes = vec![50, 500];
    let stream_nums = vec![1, 5];
    let buffer_sizes = vec![10_000, 100_000];
    let backpressure_pause_ms_options = [50, 500];
    let sleep_ms = vec![0, 1];
    let commands = vec![("cat", vec!["-"])];


    let mut log = std::fs::File::create("stream_to_cmd_stress.log")?;
    writeln!(
        &mut log,
        "Command\tBuffer_Size\tSleep\tBackpressurePause\tStreams\tReads\tSize\tTime\tMemory\tSuccess?"
    )?;
    log.flush()?;

    let mut sys = System::new();
    let pid = Pid::from(std::process::id() as usize);

    for (cmd_tag, args) in &commands {
        for buffer_size in &buffer_sizes {
            for sleep in &sleep_ms {
                for &backpressure_pause_ms in &backpressure_pause_ms_options { 

                for stream_num in &stream_nums {

                    for num_read in &num_reads {
                        for read_size in &read_sizes {
                            let mut run_success = true;
                            let start = Instant::now();
                            sys.refresh_processes(ProcessesToUpdate::Some(&[pid]), true);
                            let memory_before = sys.process(pid).map(|p| p.memory() / 1024 / 1024).unwrap_or(0);

                            eprintln!(
                                "Starting: Command: {}, Buffer: {}, Sleep: {}, Streams: {}, Reads: {}, Size: {}",
                                cmd_tag, buffer_size, sleep, stream_num, num_read, read_size
                            );
                            stderr().flush()?;

                            // Generate stream and split with t_junction
                            let stream = fastx_generator(*num_read, *read_size, 35.0, 3.0);
                            let mut gen_count = 0;
                            let stream = stream.inspect(move |_result| {
                                gen_count += 1;
                            });
                            let (mut outputs, done_rx) = match t_junction(
                                stream,
                                *stream_num,
                                *buffer_size,
                                10000,
                                if *sleep == 0 { None } else { Some(*sleep) }, backpressure_pause_ms
                            )
                                .await {
                                Ok(result) => result,
                                Err(e) => {
                                    eprintln!("t_junction failed to start: {}", e);
                                    run_success = false;
                                    return Err(anyhow!("t_junction failed: {}", e));
                                }
                            };
                            eprintln!("t_junction started with {} streams", *stream_num);
                            stderr().flush()?;
                            
                            let mut tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
                            let mut outfiles = Vec::new();

                            for i in 0..*stream_num {
                                let stream_outfile = format!("stream_to_cmd_stress.{}.log", i);

                                if fs::metadata(&stream_outfile).is_ok() {
                                    fs::remove_file(&stream_outfile)?;
                                }
                                
                                outfiles.push(stream_outfile.clone());
                                let cmd_tag = cmd_tag.to_string();
                                let args = args.iter().map(|&s| s.to_string()).collect::<Vec<_>>();

                                let rx = outputs.pop().ok_or_else(|| anyhow!("Missing stream for cmd"))?;
                                let (mut child, stream_task) = match stream_to_cmd(
                                    rx,
                                    &cmd_tag,
                                    args.iter().map(|s| s.as_str()).collect(),
                                )
                                    .await {
                                    Ok(result) => result,
                                    Err(e) => {
                                        eprintln!("Stream {} failed to spawn {}: {}", i, cmd_tag, e);
                                        run_success = false;
                                        continue;
                                    }
                                };

                                let out_stream = parse_child_output(& mut child, ChildStream::Stdout, ParseMode::Fastq, *buffer_size).await?;
                                let write_task = tokio::spawn(stream_to_file(
                                    out_stream,
                                    Path::new(&stream_outfile).to_path_buf(),
                                ));

                                tasks.push(write_task);
                                tasks.push(stream_task);
                            }

                            let results = join_all(tasks).await;
                            for result in results {
                                result??; // Unwrap JoinHandle and inner Result, propagating errors
                            }

                            done_rx.await??; // Propagates any t_junction errors

                            for outfile in outfiles {
                                
                                let mut cmd = String::new();
                                let mut args_vec: Vec<String> = Vec::new();
                                match *cmd_tag {
                                    "cat" => {
                                        cmd = "wc".to_string();
                                        args_vec.push("-l".to_string());
                                        args_vec.push(outfile.to_string());
                                        
                                    }


                                    _ => {
                                        return Err(anyhow!("Cannot use this command {}", cmd_tag));
                                    }
                                }

                                let mut wc_child = Command::new(&cmd)
                                    .args(&args_vec)
                                    .stdin(std::process::Stdio::piped())
                                    .stdout(std::process::Stdio::piped())
                                    .stderr(std::process::Stdio::piped())
                                    .spawn()
                                    .map_err(|e| anyhow!("Failed to spawn {}: {}.", cmd, e))?;
                                let lines = read_child_output_to_vec(&mut wc_child, ChildStream::Stdout).await?;
                                let first_line = lines
                                    .first()
                                    .ok_or_else(|| anyhow!("No output from fastp -v"))?;

                                let first_cols = first_line.split_whitespace().collect::<Vec<_>>();

                                
                                let line_count : usize = first_cols[0].parse()?;
                                let fastq_count = line_count/4;

                                if fastq_count != *num_read {
                                    run_success = false;
                                }
                                
                                assert_eq!(fastq_count, *num_read);
                            }
                            
                            let elapsed = start.elapsed();
                            let elapsed_secs = elapsed.as_secs_f64();
                            sys.refresh_processes(ProcessesToUpdate::Some(&[pid]), true);
                            let memory_after = sys.process(pid).map(|p| p.memory() / 1024 / 1024).unwrap_or(0);
                            let memory_used = if memory_after >= memory_before {
                                memory_after - memory_before
                            } else {
                                0
                            };
                            
                            writeln!(
                                &mut log,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                cmd_tag,
                                buffer_size,
                                sleep,
                                backpressure_pause_ms,
                                stream_num,
                                num_read,
                                read_size,
                                elapsed_secs,
                                memory_used,
                                run_success
                            )?;
                            log.flush()?;

                            eprintln!("done");
                            stderr().flush()?;

                            if !run_success {
                                return Err(anyhow!("Test failed due to errors in stream_to_cmd or t_junction"));
                            }
                        }
                    }
                }
            }
        }
    }
    }
    Ok(())
}
