use std::path::PathBuf;
use anyhow::Result;
use crate::utils::Arguments;
use tokio_stream::wrappers::{BroadcastStream, ReceiverStream};
use crate::utils::command::generate_cli;
use crate::utils::fastx::SequenceRecord;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base};
use crate::utils::streams::{stream_to_cmd, t_junction, stream_to_file, parse_child_stdout_to_bytes};
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
        Some(prefix) => {sample_base = prefix;}
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
    
    
    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(sample_base), &cwd, None, Option::from("validated"), "_");
    let mut rx = read_and_interleave_sequences(file1_path, file2_path, technology, args.max_reads, args.min_read_len, args.max_read_len)?;

    let rx_stream = ReceiverStream::new(rx); // rx is a stream of SequenceRecord, must convert before splitting.

    let (streams, done_rx) = t_junction(
        rx_stream,
        2,
        args.stall_threshold,
        Some(args.stream_sleep_ms),
    )
        .await?;
    
    let [fastp_stream, pigz_stream] = streams
        .try_into()
        .map_err(|_| anyhow::anyhow!("Expected exactly 2 streams"))?;

    // Pigz stream
    let pigz_args = generate_cli(PIGZ_TAG, &args)?;
    let mut pigz_child = stream_to_cmd(pigz_stream, PIGZ_TAG, pigz_args).await?;
    let pigz_stdout = pigz_child
        .stdout
        .take()
        .ok_or_else(|| anyhow::anyhow!("Failed to get stdout from pigz"))?;
    
    let pigz_stream = parse_child_stdout_to_bytes(pigz_stdout).await?;
    stream_to_file(pigz_stream, validated_interleaved_file_path).await?;

    //Fastp stream
    let fastp_args = generate_cli(FASTP_TAG, &args)?;
    let mut fastp_child = stream_to_cmd(fastp_stream, FASTP_TAG, fastp_args).await?;
    let fastp_stdout = fastp_child
        .stdout
        .take()
        .ok_or_else(|| anyhow::anyhow!("Failed to get stdout from fastp"))?;
    let fastp_stream = parse_child_stdout_to_bytes(fastp_stdout).await?;
    stream_to_file(fastp_stream, PathBuf::from("test_fastp.fa")).await?;

    done_rx.await?;
    println!("done!");
    
    // 
    Ok(())

}
