use std::path::PathBuf;
use anyhow::{anyhow, Result};
use crate::utils::Arguments;
use crate::utils::fastx::SequenceRecord;
use tokio::io::AsyncBufReadExt;
use tokio_stream::wrappers::{BroadcastStream, ReceiverStream};
use crate::utils::command::generate_cli;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base};
use crate::utils::streams::{t_junction, stream_to_file, stream_bytes_to_file, stream_to_cmd, parse_child_stdout_to_bytes, parse_child_stdout_to_fastq};
use crate::{PIGZ_TAG, FASTP_TAG};

pub async fn run(args: &Arguments) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let cwd = std::env::current_dir()?;

    let file1_path = file_path_manipulator(&PathBuf::from(&args.file1), &cwd, None, None, "");
    eprintln!("{}", file1_path.display());
    let sample_base: String;
    let file1_r1r2 = r1r2_base(&file1_path);
    match file1_r1r2.file_name {
        Some(prefix) => { sample_base = prefix; }
        None => {
            eprintln!("No R1 tag found. Using bare file 1 stem as sample_base.");
            sample_base = file1_path.to_string_lossy().into_owned();
        },
    }

    let file2_path: Option<PathBuf> = match &args.file2 {
        Some(file) => {
            Some(file_path_manipulator(&PathBuf::from(file), &cwd, None, None, ""))
        },
        None => {
            eprintln!("File2 not given");
            None
        },
    };

    let technology = Some(args.technology.clone());

    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(&sample_base), &cwd, None, Some("validated"), "_");
    let rx = read_and_interleave_sequences(file1_path, file2_path, technology, args.max_reads, args.min_read_len, args.max_read_len)?;
    let rx_stream = ReceiverStream::new(rx);

    let (streams, mut done_rx) = t_junction(
        rx_stream,
        2,
        args.stall_threshold * 2,
        Some(args.stream_sleep_ms),
    ).await?;

    let [fastp_stream, pigz_stream]: [BroadcastStream<SequenceRecord>; 2] = streams.try_into().map_err(|_| anyhow!("Expected exactly 2 streams"))?;

    // Pigz stream
    let pigz_args = generate_cli(PIGZ_TAG, &args)?;
    // eprintln!("pigz_args: {:?}", pigz_args);
    let pigz_args: Vec<&str> = pigz_args.iter().map(|s| s.as_str()).collect();
    let mut pigz_child = stream_to_cmd(pigz_stream, PIGZ_TAG, pigz_args).await?;
    let pigz_stdout = pigz_child.stdout.take().ok_or_else(|| anyhow!("Failed to get stdout from pigz"))?;
    let pigz_stream = parse_child_stdout_to_bytes(pigz_stdout).await?;
    let mut pigz_task = tokio::spawn(async move {
        let result = stream_bytes_to_file(pigz_stream, validated_interleaved_file_path.clone()).await;
        result
    });

    // Fastp stream
    let fastp_args = generate_cli(FASTP_TAG, &args)?;
    // eprintln!("fastp_args: {:?}", fastp_args);
    let fastp_args: Vec<&str> = fastp_args.iter().map(|s| s.as_str()).collect();
    let mut fastp_child = stream_to_cmd(fastp_stream, FASTP_TAG, fastp_args).await?;
    let fastp_stdout = fastp_child.stdout.take().ok_or_else(|| anyhow!("Failed to get stdout from fastp"))?;
    let fastp_stderr = fastp_child.stderr.take();
    let (fastp_tx, fastp_rx) = tokio::sync::mpsc::channel(20000); // Use mpsc channel
    let mut fastp_parse_task = tokio::spawn(async move {
        let result = parse_child_stdout_to_fastq(fastp_stdout, fastp_tx).await;

        result
    });
    let mut fastp_write_task = tokio::spawn(async move {
        let result = stream_to_file(fastp_rx, PathBuf::from("test_fastp.fq")).await;
        result
    });
    let mut fastp_stderr_task = tokio::spawn(async move {
        if let Some(fastp_stderr) = fastp_stderr {
            let mut stderr = tokio::io::BufReader::new(fastp_stderr);
            let mut buffer = String::new();
            while let Ok(bytes) = stderr.read_line(&mut buffer).await {
                if bytes == 0 { break; }
                // eprintln!("fastp stderr: {}", buffer.trim());  // turnt his back on once you set up verbosity
                buffer.clear();
            }
        } else {
        }
        Ok::<(), anyhow::Error>(())
    });

    // Wait for tasks and processes with select!
    let mut fastp_done = false;
    let mut pigz_done = false;
    let mut pigz_task_done = false;
    let mut fastp_parse_done = false;
    let mut fastp_write_done = false;
    let mut fastp_stderr_done = false;
    let mut t_junction_done = false;

    loop {
        tokio::select! {
            result = &mut fastp_parse_task, if !fastp_parse_done => {
                let result = result.map_err(|e| anyhow!("Fastp parse task failed: {}", e))?;
                fastp_parse_done = true;
            }
            result = &mut fastp_write_task, if !fastp_write_done => {
                let result = result.map_err(|e| anyhow!("Fastp write task failed: {}", e))?;
                fastp_write_done = true;
            }
            result = &mut fastp_stderr_task, if !fastp_stderr_done => {
                let result = result.map_err(|e| anyhow!("Fastp stderr task failed: {}", e))?;
                fastp_stderr_done = true;
            }
            status = fastp_child.wait(), if !fastp_done && fastp_parse_done && fastp_write_done && fastp_stderr_done => {
                let status = status.map_err(|e| anyhow!("Failed to wait for fastp: {}", e))?;
                if !status.success() {
                    return Err(anyhow!("Fastp exited with non-zero status: {}", status));
                }
                fastp_done = true;
            }
            status = pigz_child.wait(), if !pigz_done => {
                let status = status.map_err(|e| anyhow!("Failed to wait for pigz: {}", e))?;
                if !status.success() {
                    return Err(anyhow!("Pigz exited with non-zero status: {}", status));
                }
                pigz_done = true;
            }
            result = &mut pigz_task, if !pigz_task_done => {
                let result = result.map_err(|e| anyhow!("Pigz task failed: {}", e))?;

                pigz_task_done = true;
            }
            result = &mut done_rx, if !t_junction_done => {

                match result {
                    Ok(()) => t_junction_done = true,
                    Err(e) => return Err(anyhow!("t_junction failed: {}", e)),
                }
            }
            else => break,
        }
        if fastp_done && pigz_done && pigz_task_done && fastp_parse_done && fastp_write_done && fastp_stderr_done && t_junction_done {
            break;
        }
    }

    println!("done!");
    Ok(())
}