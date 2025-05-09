use std::path::PathBuf;
use anyhow::{anyhow, Result};
use crate::utils::Arguments;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::command::generate_cli;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base};
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file};
use crate::utils::defs::{PIGZ_TAG, FASTP_TAG};

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
    let (val_streams, val_done_rx) = t_junction(
        rx_stream,
        2,
        args.buffer_size,
        args.stall_threshold.try_into().unwrap(),
        Some(args.stream_sleep_ms),
        50,
    )
        .await?;

    if val_streams.len() != 2 {
        return Err(anyhow!("Expected exactly 2 streams, got {}", val_streams.len()));
    }
    
    let mut streams_iter = val_streams.into_iter();
    let mut fastp_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing fastp stream"))?;
    let mut pigz_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing file stream"))?;
    
    //Pigz stream to intermediate file output
    let pigz_args = generate_cli(PIGZ_TAG, &args)?;
    let pigz_args: Vec<&str> = pigz_args.iter().map(|s| s.as_str()).collect();
    let (mut pigz_child, _pigz_stream_task) = stream_to_cmd(pigz_stream, PIGZ_TAG, pigz_args, StreamDataType::IlluminaFastq).await?;

    let pigz_out_stream = parse_child_output(
        &mut pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size,
    ).await?;
    let pigz_write_task = tokio::spawn(stream_to_file(
        pigz_out_stream,
        validated_interleaved_file_path,
    ));
    
    
    // Fastp stream
    let fastp_args = generate_cli(FASTP_TAG, &args)?;
    let fastp_args: Vec<&str> = fastp_args.iter().map(|s| s.as_str()).collect();
    let (mut fastp_child, _fastp_stream_task) = stream_to_cmd(fastp_stream, FASTP_TAG, fastp_args, StreamDataType::IlluminaFastq).await?;
    let fastp_out_stream = parse_child_output(
        &mut fastp_child,
        ChildStream::Stdout,
        ParseMode::Fastq,
        args.buffer_size,
    ).await?;
    
    // just for now write out
    let mut fastp_write_task = tokio::spawn(stream_to_file(
        fastp_out_stream,
        PathBuf::from("fastp_out_test.fq"),
    ));


    pigz_write_task.await??;
    fastp_write_task.await??;

    // Check child exit statuses
    let pigz_status = pigz_child.wait().await?;
    if !pigz_status.success() {
        return Err(anyhow!("pigz exited with non-zero status: {}", pigz_status));
    }
    let fastp_status = fastp_child.wait().await?;
    if !fastp_status.success() {
        return Err(anyhow!("fastp exited with non-zero status: {}", fastp_status));
    }

    // Check t_junction completion
    val_done_rx.await??;

    

    Ok(())
}