use anyhow::anyhow;
use seqtoid_pipelines::utils::fastx::fastx_generator;
use anyhow::Result;
use std::io::{stderr, Write};
use std::time::{Instant, Duration};
use futures::StreamExt;
use sysinfo::{System, Pid};
use seqtoid_pipelines::utils::streams::{stream_to_cmd, t_junction, ToBytes};
use tokio::task::JoinHandle;
use tokio::time::timeout;
use tokio_stream::wrappers::BroadcastStream;
use tokio::process::{Child, Command};

#[tokio::test]
async fn test_fastx_generator_stress() -> Result<()> {
    // Test parameters
    let num_reads = vec![10_000, 100_000]; // Matches t_junction failing case
    let read_sizes = vec![100, 150, 1000, 5000]; // Illumina and ONT

    // Log results to file
    let mut log = std::fs::File::create("fastx_stress.log")?;
    writeln!(
        &mut log,
        "Reads\tSize\tRecords\tTime\tMemory"
    )?;
    log.flush()?;
    
    // Optional system monitoring
    #[cfg(feature = "sysinfo")]
    let mut sys = System::new_all();

    for read_size in &read_sizes {
        

        for num_read in &num_reads {


            eprintln!(
                "Testing: Reads: {}, Size: {}",
                num_read, read_size
            );
            stderr().flush()?;

            let start = Instant::now();

            #[cfg(feature = "sysinfo")]
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
                #[cfg(feature = "sysinfo")]
                {
                    sys.refresh_memory();
                    let memory_after = sys.used_memory();
                    if memory_after >= memory_before {
                        (memory_after - memory_before) / 1024 / 1024 // MB
                    } else {
                        0 // Handle underflow
                    }
                }
                #[cfg(not(feature = "sysinfo"))]
                0
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
    let buffer_sizes = vec![10_000, 100_000];
    let stall_thresholds = vec![100, 10_000];
    let sleep_ms = vec![0, 1];
    let stream_nums = vec![2, 5];
    let num_reads = vec![10000, 100000];
    let read_sizes  = vec![100, 1000];

    let mut log = std::fs::File::create("stress_test.log")?;
    writeln!(&mut log, "Buffer_Size\tStall\tSleep\tStreams\tReads\tSize\tTime\tMemory\tRecords\tSuccess?")?;
    let mut sys = System::new();
    let pid = Pid::from(std::process::id() as usize);

    for buffer_size in &buffer_sizes {
        for stall in &stall_thresholds {
            for sleep in &sleep_ms {
                for stream_num in &stream_nums {
                    for num_read in &num_reads {
                        for read_size in &read_sizes {
                            let mut run_success = true;
                            let start = Instant::now();
                            sys.refresh_process(pid);
                            let memory_before = sys.process(pid).map(|p| p.memory() / 1024 / 1024).unwrap_or(0);
                            
                            eprintln!("Buffer size: {}  Stall: {}  Sleep: {}  Streams: {} Reads: {}  Size: {}", buffer_size, stall, sleep, stream_num, num_read, read_size);
                            std::io::stderr().flush()?;

                            let stream = fastx_generator(*num_read, *read_size, 35.0, 3.0);
                            let (mut outputs, done_rx) = t_junction(
                                stream,
                                *stream_num,
                                *buffer_size,
                                *stall,
                                if *sleep == 0 { None } else { Some(*sleep) }
                            )
                                .await?;

                            let mut records = vec![Vec::new(); *stream_num];
                            for i in 0..*stream_num {
                                while let Some(result) = outputs[i].next().await {
                                    match result {
                                        Ok(record) => records[i].push(record),
                                        Err(e) => eprintln!("Stream {} error: {}", i, e),
                                    }
                                }
                                eprintln!("Stream {} received {} records", i, records[i].len());
                                stderr().flush()?;
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
                            

                            for i in 0..*stream_num {
                                eprintln!("stream {}", i);
                                stderr().flush()?;
                                if num_read != &records[i].len() {
                                    eprintln!("Stream {} failed: expected {}, got {}", i, num_read, records[i].len());
                                    run_success = false;
                                }
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
                            };
                            stderr().flush()?;

                            let record_counts: String = records
                                .iter()
                                .map(|r| r.len().to_string())
                                .collect::<Vec<_>>()
                                .join(",");
                            
                            writeln!(& mut log, "{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}\t{}\t{}", buffer_size, stall, sleep, stream_num, num_read, read_size, elapsed_secs, memory_used, record_counts, run_success)?;
                            log.flush()?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}


#[tokio::test]
async fn test_stream_to_cmd_stress() -> Result<()> {
    let num_reads = vec![10_000];
    let read_sizes = vec![50, 100, 1000];
    let stream_nums = vec![1];
    let buffer_sizes = vec![10_000, 100_000];
    let sleep_ms = vec![0, 1];
    let commands = vec![
        ("cat", vec!["-"]),
        ("gzip", vec!["-c"]),
    ];
    let timeout_secs = 300;

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
                            let (mut outputs, done_rx) = match t_junction(
                                stream,
                                *stream_num + 1, // Extra stream for counting
                                *buffer_size,
                                100,
                                if *sleep == 0 { None } else { Some(*sleep) },
                            )
                                .await {
                                Ok(result) => result,
                                Err(e) => {
                                    eprintln!("t_junction failed to start: {}", e);
                                    run_success = false;
                                    writeln!(
                                        &mut log,
                                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                        cmd_tag,
                                        buffer_size,
                                        sleep,
                                        stream_num,
                                        num_read,
                                        read_size,
                                        0.0,
                                        0,
                                        "0",
                                        false
                                    )?;
                                    log.flush()?;
                                    continue;
                                }
                            };

                            // Spawn stream_to_cmd for each stream
                            let mut children: Vec<Child> = Vec::new();
                            let mut tasks: Vec<JoinHandle<(Result<(), anyhow::Error>, usize)>> = Vec::new();
                            let mut record_counts = vec![0usize; *stream_num];

                            for i in 0..*stream_num {
                                let cmd_tag = cmd_tag.to_string();
                                let args = args.iter().map(|&s| s.to_string()).collect::<Vec<_>>();
                                let num_read_owned = *num_read; // Copy usize for 'static
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
                                        // Count records from rx_count
                                        while let Some(result) = rx_count.next().await {
                                            match result {
                                                Ok(_) => count += 1,
                                                Err(e) => return Err(anyhow!("Broadcast stream error: {}", e)),
                                            }
                                        }
                                        // Wait for inner_task to complete
                                        inner_task.await??;
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
                            for (i, mut child) in children.into_iter().enumerate() {
                                match timeout(Duration::from_secs(timeout_secs), child.wait_with_output()).await {
                                    Ok(Ok(output)) => {
                                        if !output.status.success() {
                                            eprintln!(
                                                "Stream {} process {} failed with status: {}, stderr: {}",
                                                i, cmd_tag, output.status, String::from_utf8_lossy(&output.stderr)
                                            );
                                            run_success = false;
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
                                Ok(Ok(Ok(()))) => {}
                                Ok(Ok(Err(e))) => {
                                    eprintln!("t_junction failed: {}", e);
                                    run_success = false;
                                }
                                Ok(Err(e)) => {
                                    eprintln!("t_junction task failed to send: {}", e);
                                    run_success = false;
                                }
                                Err(_) => {
                                    eprintln!("t_junction timed out after {}s", timeout_secs);
                                    run_success = false;
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
                        }
                    }
                }
            }
        }
    }
    Ok(())
}