use std::fs;
use anyhow::anyhow;
use seqtoid_pipelines::utils::fastx::{fastx_generator, SequenceRecord};
use anyhow::Result;
use std::io::{stderr, Write};
use std::path::Path;
use std::time::Instant;
use log::{self, LevelFilter, debug, info, error, warn};
use futures::StreamExt;
use sysinfo::{System, Pid, ProcessesToUpdate};
use seqtoid_pipelines::utils::streams::{read_child_output_to_vec, stream_to_cmd, parse_child_output, stream_to_file, ChildStream, ParseMode, ParseOutput, deinterleave_fastq_stream_to_fifos, fanout_to_channels};
use tokio::task::JoinHandle;
use tokio::time::timeout;
use tokio::process::{Child, Command};
use tokio_stream::wrappers::ReceiverStream;
use std::sync::{Arc, Mutex};
use tokio::time::Duration;
use tokio::fs::File as TokioFile;
use futures::future::join_all;
use seqtoid_pipelines::config::defs::{GpuDetection, NRAlignmentBackend, RunConfig, StreamDataType};
use std::path::PathBuf;
use rayon::ThreadPoolBuilder;
use tokio::io::AsyncReadExt;
use tokio::sync::{mpsc, Semaphore};
use seqtoid_pipelines::cli::Arguments;
use seqtoid_pipelines::utils::system::{detect_cores_and_load, detect_ram, generate_rng, compute_buffer_size};
use seqtoid_pipelines::config::defs::SimdLevel;

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
            info!(
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
                &mut log,
                "Testing fastx_generator: Reads: {}, Size: {}",
                num_read, read_size
            )?;
            stderr().flush()?;
            let stream = fastx_generator(num_read, *read_size, 35.0, 3.0);
            let count = stream.fold(0usize, |acc, _| async move { acc + 1 }).await;
            if count != num_read {
                writeln!(
                    &mut log,
                    "fastx_generator failed: expected {}, got {}",
                    num_read, count
                )?;
                return Err(anyhow!(
                    "fastx_generator produced {} records, expected {}",
                    count, num_read
                ));
            }
            writeln!(
                &mut log,
                "fastx_generator produced {} records",
                count
            )?;
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
        run_id: "NULL".to_string(),
        efs_base_dir: PathBuf::from("/dev/null"),
        sample_base: PathBuf::from("stress"),
    };

    // Compute proper buffer size exactly like main.rs
    let sdt = StreamDataType::IlluminaFastq; // default for tests

    let base_buffer_size = compute_buffer_size(
        &run_config,
        "test_global_default",
        sdt,
        1.0,
    );

    run_config.base_buffer_size = base_buffer_size;

    Arc::new(run_config)
}


#[tokio::test]
async fn test_fanout_to_channels_stress() -> Result<()> {
    let num_reads_options = [10_000, 100_000];
    let n_outputs_options = [2, 4];
    let seq_len = 100;
    let config = create_test_run_config();

    let mut log = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open("fanout_stress.log")?;

    writeln!(
        &mut log,
        "Streams\tReads\tSeqLen\tTime\tRecords\tSuccess"
    )?;

    for &num_reads in &num_reads_options {
        for &n_outputs in &n_outputs_options {
            info!(
                "Testing fanout_to_channels: Streams: {} Reads: {}  Size: {}",
                n_outputs,
                num_reads,
                seq_len
            );

            let stream = fastx_generator(num_reads, seq_len, 30.0, 5.0).map(ParseOutput::Fastq);
            let start = Instant::now();
            
            // Re-creating the stream because fanout_to_channels takes ownership of ReceiverStream
            let (tx, rx) = mpsc::channel(1000);
            let stream_input = ReceiverStream::new(rx);
            
            let (outputs, router_handle) = fanout_to_channels(
                stream_input,
                n_outputs,
                "fanout_stress",
                &config,
                StreamDataType::IlluminaFastq
            ).await?;

            let record_counts = Arc::new(Mutex::new(vec![0; n_outputs]));
            let mut handles = Vec::new();

            for (i, rx) in outputs.into_iter().enumerate() {
                let record_counts = Arc::clone(&record_counts);
                let handle = tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(rx);
                    while let Some(_item) = stream.next().await {
                        let mut counts = record_counts.lock().unwrap();
                        counts[i] += 1;
                    }
                });
                handles.push(handle);
            }

            // Feed the input
            let mut input_gen = fastx_generator(num_reads, seq_len, 30.0, 5.0).map(ParseOutput::Fastq);
            while let Some(item) = input_gen.next().await {
                tx.send(item).await.map_err(|e| anyhow!("Failed to send to fanout: {}", e))?;
            }
            drop(tx);

            join_all(handles).await;
            router_handle.await??;

            let duration = start.elapsed();
            let record_counts_final = record_counts.lock().unwrap().clone();
            let run_success = record_counts_final == vec![num_reads; n_outputs];
            
            info!(
                "Records: {:?}  Success: {}",
                record_counts_final, run_success
            );
            writeln!(
                &mut log,
                "{}\t{}\t{}\t{}\t{:?}\t{}",
                n_outputs,
                num_reads,
                seq_len,
                duration.as_secs_f64(),
                record_counts_final,
                run_success
            )?;

            assert!(run_success, "Incorrect record counts: expected {}, got {:?}", num_reads, record_counts_final);
        }
    }
    Ok(())
}

#[tokio::test]
async fn test_fanout_to_channels_count() -> Result<()> {
    let num_read = 10_000;
    let read_size = 50;
    let n_outputs = 2;
    let config = create_test_run_config();

    info!(
        "Testing fanout_to_channels_count: Reads: {}, Size: {}, Streams: {}",
        num_read, read_size, n_outputs
    );

    let (tx, rx) = mpsc::channel(1000);
    let stream_input = ReceiverStream::new(rx);

    let (outputs, router_handle) = fanout_to_channels(
        stream_input,
        n_outputs,
        "fanout_count",
        &config,
        StreamDataType::IlluminaFastq
    ).await?;

    let mut handles = Vec::new();
    for (i, rx) in outputs.into_iter().enumerate() {
        let handle = tokio::spawn(async move {
            let mut stream = ReceiverStream::new(rx);
            let mut count = 0;
            while let Some(_) = stream.next().await {
                count += 1;
            }
            info!("fanout stream {} finished: {} records", i, count);
            Ok::<usize, anyhow::Error>(count)
        });
        handles.push(handle);
    }

    // Feed input
    let mut input_gen = fastx_generator(num_read, read_size, 35.0, 3.0).map(ParseOutput::Fastq);
    while let Some(item) = input_gen.next().await {
        tx.send(item).await?;
    }
    drop(tx);

    let mut counts = Vec::new();
    for handle in handles {
        counts.push(handle.await??);
    }

    router_handle.await??;

    info!("fanout counts: {:?}", counts);
    assert_eq!(counts, vec![num_read; n_outputs]);
    
    Ok(())
}

#[tokio::test]
async fn test_stream_to_cmd_direct() -> Result<()> {
    let config = create_test_run_config();
    let num_read = 100;
    let read_size = 50;
    let cmd_tag = "cat";
    let args = vec!["-".to_string()];
    info!(
        "Testing stream_to_cmd: Reads: {}, Size: {}, Command: {}",
        num_read, read_size, cmd_tag
    );
    stderr().flush()?;

    let (tx, rx) = mpsc::channel(1000);
    let stream_input = ReceiverStream::new(rx);
    let (mut outputs, done_rx) = fanout_to_channels(
        stream_input,
        1,
        "test_stream_to_cmd_direct",
        &config,
        StreamDataType::IlluminaFastq
    ).await?;

    // Feed input in a background task to avoid blocking
    tokio::spawn(async move {
        let mut stream = fastx_generator(num_read, read_size, 35.0, 3.0).map(ParseOutput::Fastq);
        while let Some(item) = stream.next().await {
            let _ = tx.send(item).await;
        }
    });

    let rx = outputs.pop().ok_or_else(|| anyhow!("No output stream"))?;
    let (child, inner_task, _inner_err_task) = stream_to_cmd(
        config,
        rx,
        cmd_tag,
        args,
        StreamDataType::IlluminaFastq,
        false,
        None
    ).await?;

    match timeout(Duration::from_secs(30), inner_task).await {
        Ok(Ok(Ok(()))) => info!("stream_to_cmd completed successfully"),
        Ok(Ok(Err(e))) => {
            error!("stream_to_cmd failed: {}", e);
            return Err(e);
        }
        Ok(Err(e)) => {
            error!("stream_to_cmd join error: {}", e);
            return Err(anyhow!("Join error: {}", e));
        }
        Err(_) => {
            error!("stream_to_cmd timed out after 30s");
            return Err(anyhow!("stream_to_cmd timed out"));
        }
    }

    done_rx.await??;

    let child = Arc::try_unwrap(child)
        .map_err(|_| anyhow!("Multiple references to child remain"))?
        .into_inner();

    let output = child.wait_with_output().await?;
    let stderr_output = String::from_utf8_lossy(&output.stderr);
    if !output.status.success() {
        error!(
            "Child process failed: status: {}, stderr: {}",
            output.status, stderr_output
        );
        return Err(anyhow!(
            "Child process failed: status: {}, stderr: {}",
            output.status, stderr_output
        ));
    }
    info!("Child process completed successfully");
    stderr().flush()?;
    Ok(())
}

#[tokio::test(flavor = "multi_thread", worker_threads = 84)]
async fn test_stream_to_cmd_stress() -> Result<()> {
    let config = create_test_run_config();
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
        "Command\tBuffer_Size\tSleep\tBackpressurePause\tStreams\tReads\tSize\tTime\tMemory_Bytes\tSuccess?"
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
                                let memory_before = sys.process(pid).map(|p| p.memory()).unwrap_or(0);

                                info!(
                                    "Starting: Command: {}, Buffer: {}, Sleep: {}, Streams: {}, Reads: {}, Size: {}",
                                    cmd_tag, buffer_size, sleep, stream_num, num_read, read_size
                                );
                                stderr().flush()?;

                                // Generate stream and fan out
                                let (tx, rx) = mpsc::channel(1000);
                                let stream_input = ReceiverStream::new(rx);
                                let (outputs, done_rx) = match fanout_to_channels(
                                    stream_input,
                                    *stream_num,
                                    "test_stream_to_cmd_stress",
                                    &config,
                                    StreamDataType::IlluminaFastq,
                                )
                                    .await {
                                    Ok(result) => result,
                                    Err(e) => {
                                        error!("fanout_to_channels failed to start: {}", e);
                                        run_success = false;
                                        return Err(anyhow!("fanout_to_channels failed: {}", e));
                                    }
                                };
                                info!("fanout_to_channels started with {} streams", *stream_num);

                                // Feed input in a background task
                                let num_read_val = *num_read;
                                let read_size_val = *read_size;
                                tokio::spawn(async move {
                                    let mut stream = fastx_generator(num_read_val, read_size_val, 35.0, 3.0).map(ParseOutput::Fastq);
                                    while let Some(item) = stream.next().await {
                                        let _ = tx.send(item).await;
                                    }
                                });
                                stderr().flush()?;

                                let mut tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
                                let mut outfiles = Vec::new();

                                for (i, rx) in outputs.into_iter().enumerate() {
                                    let stream_outfile = format!("stream_to_cmd_stress.{}.fq", i);

                                    if fs::metadata(&stream_outfile).is_ok() {
                                        fs::remove_file(&stream_outfile)?;
                                    }

                                    outfiles.push(stream_outfile.clone());
                                    let cmd_tag = cmd_tag.to_string();
                                    let args = args.iter().map(|s| s.to_string()).collect::<Vec<_>>();

                                    let (child, stream_task, _stream_err_task) = match stream_to_cmd(
                                        config.clone(),
                                        rx,
                                        &cmd_tag,
                                        args,
                                        StreamDataType::IlluminaFastq,
                                        false,
                                        None
                                    )
                                        .await {
                                        Ok(result) => result,
                                        Err(e) => {
                                            error!("Stream {} failed to spawn {}: {}", i, cmd_tag, e);
                                            run_success = false;
                                            continue;
                                        }
                                    };

                                    let mut child_guard = child.lock().await;
                                    let out_stream = parse_child_output(
                                        &mut *child_guard,
                                        ChildStream::Stdout,
                                        ParseMode::Fastq,
                                        &config,
                                    ).await?;
                                    drop(child_guard);

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

                                done_rx.await??; // Propagates any fanout errors

                                for outfile in &outfiles {
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
                                    let lines = read_child_output_to_vec(&mut wc_child, ChildStream::Stdout, &config).await?;
                                    let first_line = lines
                                        .first()
                                        .ok_or_else(|| anyhow!("No output from wc"))?;

                                    let first_cols = first_line.split_whitespace().collect::<Vec<_>>();

                                    let line_count: usize = first_cols[0].parse()?;
                                    let fastq_count = line_count / 4;

                                    if fastq_count != *num_read {
                                        run_success = false;
                                    }

                                    assert_eq!(fastq_count, *num_read);
                                }

                                let elapsed = start.elapsed();
                                let elapsed_secs = elapsed.as_secs_f64();
                                sys.refresh_processes(ProcessesToUpdate::Some(&[pid]), true);
                                let memory_after = sys.process(pid).map(|p| p.memory()).unwrap_or(0);
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
                                
                                stderr().flush()?;

                                for outfile in outfiles {
                                    let _ = fs::remove_file(&outfile);
                                }

                                if !run_success {
                                    return Err(anyhow!(
                                        "Test failed due to errors in stream_to_cmd or fanout_to_channels"
                                    ));
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
