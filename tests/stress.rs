use anyhow::anyhow;
use seqtoid_pipelines::utils::fastx::fastx_generator;
use anyhow::Result;
use std::io::{stderr, Write};
use std::time::{Instant};
use futures::StreamExt;
use sysinfo::{System, Pid};
use seqtoid_pipelines::utils::streams::{stream_to_cmd, t_junction};
use tokio::task::JoinHandle;
use tokio::time::timeout;
use tokio::process::Child;
use tokio::sync::broadcast;
use tokio_stream::wrappers::BroadcastStream;
use std::sync::{Arc, Mutex};
use tokio::time::{Duration};


#[tokio::test]
async fn test_fastx_generator_stress() -> Result<()> {

    let num_reads = vec![10_000, 100_000]; 
    let read_sizes = vec![100, 150, 1000, 5000];

    // Log results to file
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


            // Log results
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

    for num_read in num_reads {
        for read_size in &read_sizes {
            eprintln!("Testing fastx_generator: Reads: {}, Size: {}", num_read, read_size);
            stderr().flush()?;
            let stream = fastx_generator(num_read, *read_size, 35.0, 3.0);
            let count = stream.fold(0usize, |acc, _| async move { acc + 1 }).await;
            if count != num_read {
                eprintln!("fastx_generator failed: expected {}, got {}", num_read, count);
                return Err(anyhow::anyhow!("fastx_generator produced {} records, expected {}", count, num_read));
            }
            eprintln!("fastx_generator produced {} records", count);
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
    println!("BufferSize\tStall\tSleep\tBackpressurePause\tStreams\tReads\tSeqLen\tTime\tMemory\tRecords\tSuccess");
    
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
                    println!(
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
                    );

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

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn test_stream_to_cmd_stress() -> Result<()> {
    let num_reads = vec![100];
    let read_sizes = vec![50];
    let stream_nums = vec![1];
    let buffer_sizes = vec![100_000];
    let sleep_ms = vec![0];
    let commands = vec![("cat", vec!["-"])];
    let timeout_secs = 30;
    let backpressure_pause_ms = 500;

    let mut log = std::fs::File::create("stream_to_cmd_stress.log")?;
    writeln!(
        &mut log,
        "Command\tBuffer_Size\tSleep\tStreams\tReads\tSize\tTime\tMemory\tRecords\tSuccess?"
    )?;
    log.flush()?;

    let mut sys = System::new();
    let pid = Pid::from(std::process::id() as usize);

    for (cmd_tag, args) in &commands {
        for buffer_size in &buffer_sizes {
            for sleep in &sleep_ms {
                for stream_num in &stream_nums {
                    for num_read in &num_reads {
                        for read_size in &read_sizes {
                            let mut run_success = true;
                            let start = Instant::now();
                            sys.refresh_process(pid);
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
                                if gen_count % 10 == 0 {
                                    eprintln!("fastx_generator produced {} records", gen_count);
                                    let _ = stderr().flush();
                                }
                            });
                            let (mut outputs, done_rx) = match t_junction(
                                stream,
                                *stream_num + 1,
                                *buffer_size,
                                100,
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
                            eprintln!("t_junction started with {} streams", *stream_num + 1);
                            stderr().flush()?;

                            // Spawn stream_to_cmd for each stream
                            let mut children: Vec<Child> = Vec::new();
                            let mut tasks: Vec<JoinHandle<(Result<(), anyhow::Error>, usize)>> = Vec::new();
                            let mut record_counts = vec![0usize; *stream_num];

                            for i in 0..*stream_num {
                                let cmd_tag = cmd_tag.to_string();
                                let args = args.iter().map(|&s| s.to_string()).collect::<Vec<_>>();
                                let num_read_owned = *num_read;
                                let rx = outputs.pop().ok_or_else(|| anyhow!("Missing stream for cmd"))?;
                                let mut rx_count = outputs.pop().ok_or_else(|| anyhow!("Missing stream for counting"))?;
                                let (child, inner_task) = match stream_to_cmd(
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
                                children.push(child);
                                let task = tokio::spawn(async move {
                                    let mut count = 0;
                                    let result = async {
                                        while let Some(result) = rx_count.next().await {
                                            match result {
                                                Ok(record) => {
                                                    count += 1;
                                                    if count <= 5 {
                                                        eprintln!("Stream {} record {}: {:?}", i, count, record.id());
                                                        stderr().flush()?;
                                                    }
                                                    if count % 10 == 0 {
                                                        eprintln!("Stream {} counted {} records", i, count);
                                                        stderr().flush()?;
                                                    }
                                                }
                                                Err(e) => return Err(anyhow!("Broadcast stream error: {}", e)),
                                            }
                                        }
                                        eprintln!("Stream {} finished counting: {} records", i, count);
                                        stderr().flush()?;
                                        match inner_task.await {
                                            Ok(Ok(())) => {}
                                            Ok(Err(e)) => return Err(anyhow!("stream_to_cmd error: {}", e)),
                                            Err(e) => return Err(anyhow!("stream_to_cmd join error: {}", e)),
                                        }
                                        if count != num_read_owned {
                                            return Err(anyhow!(
                                                "Expected {} records, received {}",
                                                num_read_owned,
                                                count
                                            ));
                                        }
                                        Ok(())
                                    }
                                        .await;
                                    (result, count)
                                });
                                tasks.push(task);
                            }

                            // Wait for tasks with timeout
                            for (i, task) in tasks.into_iter().enumerate() {
                                match timeout(Duration::from_secs(timeout_secs), task).await {
                                    Ok(Ok((Ok(()), count))) => {
                                        record_counts[i] = count;
                                        eprintln!("Stream {} task completed: {} records", i, count);
                                        stderr().flush()?;
                                    }
                                    Ok(Ok((Err(e), count))) => {
                                        eprintln!("Stream {} task failed: {}, received {} records", i, e, count);
                                        record_counts[i] = count;
                                        run_success = false;
                                    }
                                    Ok(Err(e)) => {
                                        eprintln!("Stream {} task join failed: {}", i, e);
                                        run_success = false;
                                    }
                                    Err(_) => {
                                        eprintln!("Stream {} task timed out after {}s", i, timeout_secs);
                                        run_success = false;
                                    }
                                }
                            }

                            // Verify process completion
                            for (i, child) in children.into_iter().enumerate() {
                                match timeout(Duration::from_secs(timeout_secs), child.wait_with_output()).await {
                                    Ok(Ok(output)) => {
                                        let stderr_output = String::from_utf8_lossy(&output.stderr);
                                        if !output.status.success() {
                                            eprintln!(
                                                "Stream {} process {} failed with status: {}, stderr: {}",
                                                i, cmd_tag, output.status, stderr_output
                                            );
                                            run_success = false;
                                        } else {
                                            eprintln!("Stream {} process {} completed successfully", i, cmd_tag);
                                            stderr().flush()?;
                                        }
                                    }
                                    Ok(Err(e)) => {
                                        eprintln!("Stream {} process failed: {}", i, e);
                                        run_success = false;
                                    }
                                    Err(_) => {
                                        eprintln!("Stream {} process timed out after {}s", i, timeout_secs);
                                        run_success = false;
                                    }
                                }
                            }

                            // Check t_junction completion
                            match timeout(Duration::from_secs(timeout_secs), done_rx).await {
                                Ok(Ok(Ok(()))) => {
                                    eprintln!("t_junction completed successfully");
                                    stderr().flush()?;
                                }
                                Ok(Ok(Err(e))) => {
                                    eprintln!("t_junction failed: {}", e);
                                    run_success = false;
                                    return Err(anyhow!("t_junction failed: {}", e));
                                }
                                Ok(Err(e)) => {
                                    eprintln!("t_junction task failed to send: {}", e);
                                    run_success = false;
                                    return Err(anyhow!("t_junction send error: {}", e));
                                }
                                Err(_) => {
                                    eprintln!("t_junction timed out after {}s", timeout_secs);
                                    run_success = false;
                                    return Err(anyhow!("t_junction timed out"));
                                }
                            }

                            let elapsed = start.elapsed();
                            let elapsed_secs = elapsed.as_secs_f64();
                            sys.refresh_process(pid);
                            let memory_after = sys.process(pid).map(|p| p.memory() / 1024 / 1024).unwrap_or(0);
                            let memory_used = if memory_after >= memory_before {
                                memory_after - memory_before
                            } else {
                                0
                            };

                            let record_counts_str = record_counts
                                .iter()
                                .map(|c| c.to_string())
                                .collect::<Vec<_>>()
                                .join(",");
                            writeln!(
                                &mut log,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                cmd_tag,
                                buffer_size,
                                sleep,
                                stream_num,
                                num_read,
                                read_size,
                                elapsed_secs,
                                memory_used,
                                record_counts_str,
                                run_success
                            )?;
                            log.flush()?;

                            eprintln!(
                                "Completed: Command: {}, Buffer: {}, Sleep: {}, Streams: {}, Reads: {}, Size: {}, Success: {}, Records: {}",
                                cmd_tag, buffer_size, sleep, stream_num, num_read, read_size, run_success, record_counts_str
                            );
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
    Ok(())
}