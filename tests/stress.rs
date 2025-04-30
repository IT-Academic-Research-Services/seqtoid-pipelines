use seqtoid_pipelines::utils::fastx::{fastx_generator, SequenceRecord};
use anyhow::Result;
use std::io::{stderr, Write};
use std::time::Instant;
use futures::StreamExt;
use sysinfo::System;
use seqtoid_pipelines::utils::streams::t_junction;

#[tokio::test]
async fn test_fastx_generator_stress() -> Result<()> {
    // Test parameters
    let num_reads = vec![10_000, 100_000]; // Matches t_junction failing case
    let read_sizes = vec![100, 150, 1000, 5000]; // Illumina and ONT

    // Log results to file
    let mut log = std::fs::File::create("fastx_stress.log")?;
    writeln!(
        &mut log,
        "Dataset\tReads\tSize\tRecords\tTime\tMemory"
    )?;
    log.flush()?;

    // Optional system monitoring
    #[cfg(feature = "sysinfo")]
    let mut sys = System::new_all();

    for read_size in &read_sizes {
        let dataset_name = if *read_size <= 150 {
            format!("illumina_{}bp", read_size)
        } else {
            format!("ont_{}bp", read_size)
        };
        let max_records = if *read_size >= 1000 { 10_000 } else { 100_000 };

        for num_read in &num_reads {
            if *num_read > max_records {
                continue;
            }

            eprintln!(
                "Testing: Dataset: {}, Reads: {}, Size: {}",
                dataset_name, num_read, read_size
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

            // Verify record count
            assert_eq!(
                *num_read,
                record_count,
                "Dataset: {}, Reads: {}, Size: {}",
                dataset_name, num_read, read_size
            );

            eprintln!(
                "Dataset: {} produced {} records",
                dataset_name, record_count
            );
            stderr().flush()?;

            // Log results
            writeln!(
                &mut log,
                "{}\t{}\t{}\t{}\t{}\t{}",
                dataset_name, num_read, read_size, record_count, elapsed_secs, memory_used
            )?;
            log.flush()?;
        }
    }

    Ok(())
}

#[tokio::test]
async fn test_t_junction_stress() -> Result<()> {
    let buffer_sizes = vec![10_000, 50_000, 100_000, 1_000_000];
    let stall_thresholds = vec![100, 1_000, 10_000, 50_000];
    let sleep_ms = vec![0, 1, 2, 5, 10];
    let stream_nums = vec![2, 5, 10];
    let num_reads = vec![10000, 100000, 1000000];
    let read_sizes  = vec![100, 1000, 10000];

    let mut log = std::fs::File::create("stress_test.log")?;
    writeln!(&mut log, "Buffer size\tStall\tSleep\tStreams\tReads\tSize\tTime\tMemory")?;
    let mut sys = System::new_all();

    for buffer_size in &buffer_sizes {
        for stall in &stall_thresholds {
            for sleep in &sleep_ms {
                for stream_num in &stream_nums {
                    for num_read in &num_reads {
                        for read_size in &read_sizes {
                            let start = Instant::now();
                            sys.refresh_memory();
                            let memory_before = sys.used_memory();
                            eprintln!("Buffer size: {}  Stall: {}  Sleep: {}  Streams: {} Reads: {}  Size: {}", buffer_size, stall, sleep, stream_num, num_read, read_size);
                            std::io::stderr().flush()?;
                            let stream = fastx_generator(*num_read, *read_size, 35.0, 3.0);
                            let (mut outputs, done_rx) = t_junction(
                                stream,
                                *stream_num,
                                *buffer_size,
                                *stall,
                                if *sleep == 0 { None } else { Some(*sleep) },
                            ).await?;

                            let mut records = vec![Vec::new(); *stream_num];
                            for i in 0..*stream_num {
                                while let Some(Ok(record)) = outputs[i].next().await {
                                    records[i].push(record);
                                }
                            }
                            let elapsed = start.elapsed();
                            let elapsed_secs = elapsed.as_secs_f64();
                            sys.refresh_memory();
                            let memory_after = sys.used_memory();
                            let memory_used = if memory_after >= memory_before {
                                (memory_after - memory_before) / 1024 / 1024 // MB
                            } else {
                                0 // Memory usage decreased, likely due to system fluctuations
                            };

                            for i in 0..*stream_num {
                                eprintln!("stream {}", i);
                                std::io::stderr().flush()?;
                                assert_eq!(num_read, &records[i].len());
                            }
                            match done_rx.await {
                                Ok(result) => match result {
                                    Ok(()) => eprintln!("t_junction completed successfully"),
                                    Err(e) => {
                                        eprintln!("t_junction failed: {}", e);
                                        return Err(e);
                                    }
                                },
                                Err(e) => {
                                    eprintln!("t_junction task failed to send: {}", e);
                                    return Err(anyhow::anyhow!("t_junction task failed: {}", e));
                                }
                            };

                            writeln!(& mut log, "{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}", buffer_size, stall, sleep, stream_num, num_read, read_size, elapsed_secs, memory_used)?;
                            log.flush()?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}