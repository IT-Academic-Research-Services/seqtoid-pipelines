use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::Arguments;
use std::process::Command;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::command::{generate_cli, check_version};
use crate::utils::file::{file_path_manipulator};
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, write_fasta_to_fifo};
use crate::utils::db::write_hdf5_seq_to_fifo;
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG};
use crate::utils::db::{lookup_sequence, load_index, build_new_in_memory_index};

pub async fn run(args: &Arguments) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let cwd = std::env::current_dir()?;
    let verbose = args.verbose.clone();

    
    // Arguments and files check

    let file1_path: PathBuf = match &args.file1 {
        Some(file) => {
            let file1_full_path = file_path_manipulator(&PathBuf::from(file), &cwd, None, None, "");
            if file1_full_path.exists() {
                file1_full_path
            } else {
                return Err(anyhow!("File1 path does not exist: {}", file1_full_path.display()));
            }
        }
        None => {
            return Err(anyhow!("File1 path required"));
        }
    };
    
    
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
    if !file1_path.exists() {
        return Err(anyhow!("Specficied file1 {:?} does not exist.", file1_path));
    }

    let file2_path: Option<PathBuf> = match &args.file2 {
        Some(file) => {
            let file2_full_path = file_path_manipulator(&PathBuf::from(file), &cwd, None, None, "");
            if file2_full_path.exists() {
                Some(file2_full_path)
            } else {
                eprintln!("File2 path does not exist: {}", file2_full_path.display());
                None
            }
        }
        None => {
            None
        }
    };

    let technology = args.technology.clone();
    
    let ref_db = args.ref_db.clone().ok_or_else(|| {
        anyhow!("HDF5 database file must be given (-d).")
    })?;
    let ref_db_path = PathBuf::from(&ref_db);

    

    //*****************
    // Input Validation
    
    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(&sample_base), &cwd, None, Some("validated"), "_");
    let rx = read_and_interleave_sequences(file1_path, file2_path, Some(technology.clone()), args.max_reads, args.min_read_len, args.max_read_len)?;
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
    let pigz_args = generate_cli(PIGZ_TAG, &args, None)?;
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
    let fastp_args = generate_cli(FASTP_TAG, &args, None)?;
    let fastp_args: Vec<&str> = fastp_args.iter().map(|s| s.as_str()).collect();
    let (mut fastp_child, fastp_stream_task) = stream_to_cmd(fastp_stream, FASTP_TAG, fastp_args, StreamDataType::IlluminaFastq).await?;
    let fastp_out_stream = parse_child_output(
        &mut fastp_child,
        ChildStream::Stdout,
        ParseMode::Bytes, // Use Bytes to avoid parsing issues
        args.buffer_size / 4, // Reduce buffer size
    ).await?;
    let mut fastp_out_stream = ReceiverStream::new(fastp_out_stream);
    
    let fastp_err_stream = parse_child_output(
        &mut fastp_child,
        ChildStream::Stderr,
        ParseMode::Bytes,
        args.buffer_size,
    ).await?;
    let fastp_err_task = tokio::spawn(async move {
        let mut err_stream = ReceiverStream::new(fastp_err_stream);
        if verbose {
            while let Some(ParseOutput::Bytes(chunk)) = err_stream.next().await {
                eprintln!("fastp stderr: {}", String::from_utf8_lossy(&chunk));
            }
        } else {
            // Consume the stream without processing
            while err_stream.next().await.is_some() {}
        }
        
        Ok::<(), anyhow::Error>(())
    });


    
    //*****************
    //Fetch Reference from accession
    
    // Retrieve index or create it for host sequence and filter sequence
    let h5_index = if let Some(index_file) = &args.ref_index {
        let index_full_path = file_path_manipulator(&PathBuf::from(index_file), &cwd, None, None, "");
        if index_full_path.exists() {
            load_index(&index_full_path).await?
        } else {
            eprintln!("Index path does not exist: {}", index_full_path.display());

            build_new_in_memory_index(&ref_db_path, &index_full_path).await?
        }
    } else {
        let index_full_path = ref_db_path.with_extension("index.bin");
        eprintln!("No index file provided, creating new index: {}", index_full_path.display());
        build_new_in_memory_index(&ref_db_path, &index_full_path).await?
    };

    
    // Create FIFO pipes
    let ref_temp = NamedTempFile::new()?;
    let ref_pipe_path = ref_temp.path().to_path_buf();
    let query_temp = NamedTempFile::new()?;
    let query_pipe_path = query_temp.path().to_path_buf();


    if ref_pipe_path.exists() {
        std::fs::remove_file(&ref_pipe_path)?;
    }
    Command::new("mkfifo")
        .arg(&ref_pipe_path)
        .status()?;
    if query_pipe_path.exists() {
        std::fs::remove_file(&query_pipe_path)?;
    }
    Command::new("mkfifo")
        .arg(&query_pipe_path)
        .status()?;
    
    
    let host_accession = args.host_accession.clone();
    let host_sequence = args.host_sequence.clone();

    // If the host sequence file is given, load it, if not retrieve it by accession from ref_db
    let seq = match &host_sequence {
        Some(_host_sequence_file) => None,
        None => {
            match &host_accession {
                Some(accession) => {
                    Some(lookup_sequence(&ref_db_path, &h5_index, &accession).await?)
                }
                None => {
                    return Err(anyhow!("Must provide either a host sequence file with --host_sequence or an accession with --host_accession"))
                }
            }
            
        },
        
    };


    // Create FIFO pipe for the fastp output to stream to minimap2
    let query_write_task = tokio::spawn({
        let query_pipe_path = query_pipe_path.clone();
        async move {
            let mut query_file = File::create(&query_pipe_path).await?;
            let mut byte_count = 0;
            while let Some(item) = fastp_out_stream.next().await {
                match item {
                    ParseOutput::Bytes(data) => {
                        query_file.write_all(&data).await?;
                        byte_count += data.len();
                    }
                    _ => return Err(anyhow!("Expected Bytes, got unexpected data")),
                }
            }
            query_file.flush().await?;
            if byte_count == 0 {
                return Err(anyhow!("No data produced by fastp"));
            }
            Ok::<(), anyhow::Error>(())
        }
    });

    // Create FIFO pipe from either the host_sequence of host_accession
    let ref_write_task = tokio::spawn({
        let ref_pipe_path = ref_pipe_path.clone();
        async move {
            
            match seq {
                Some(seq) => {
                    match &host_accession {
                        Some(accession) => {
                            write_hdf5_seq_to_fifo(seq, &accession, &ref_pipe_path).await
                        }
                        None => {
                            return Err(anyhow!("Must provide either a host sequence file with --host_sequence or an accession with --host_accession"))
                        }
                    }
                }
                None => {
                    match &host_sequence {
                        Some(host_sequence_file) => {
                            let host_sequence_path = file_path_manipulator(&PathBuf::from(host_sequence_file), &cwd, None, None, "");
                            tokio::task::spawn_blocking(move || {
                                write_fasta_to_fifo(&host_sequence_path, &ref_pipe_path)
                            }).await?
                        }
                        None => {
                            return Err(anyhow!("Must provide either a host sequence file with --host_sequence or an accession with --host_accession"))
                        }
                    }
                    
    
                }
            }
        }
    });

    
    
    //*****************
    //Host Removal

    let _minimap2_version = match check_version(MINIMAP2_TAG).await {
        Ok(version) => {
            eprintln!("{}", version);
            version
        }
        Err(err) => {
            return Err(anyhow!("Error checking minimap2 version: {}", err));
        }
    };

    let minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(ref_pipe_path.clone(), query_pipe_path.clone())))?;

    let minimap2_args: Vec<&str> = minimap2_args.iter().map(|s| s.as_str()).collect();
    let (mut minimap2_child, minimap2_task) = spawn_cmd(MINIMAP2_TAG, minimap2_args).await?;
    let minimap2_out_stream = parse_child_output(
        &mut minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    let minimap2_write_task = tokio::spawn(stream_to_file(
        minimap2_out_stream,
        PathBuf::from("test_minimap.sam"),
    ));

    let minimap2_err_stream = parse_child_output(
        &mut minimap2_child,
        ChildStream::Stderr,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    let minimap2_err_task = tokio::spawn(async move {
        let mut err_stream = ReceiverStream::new(minimap2_err_stream);
        if verbose {
            while let Some(ParseOutput::Bytes(chunk)) = err_stream.next().await {
                eprintln!("minimap2 stderr: {}", String::from_utf8_lossy(&chunk));
            }
        } else {
            // Consume the stream without processing
            while err_stream.next().await.is_some() {}
        }
        Ok::<(), anyhow::Error>(())
    });


    let _samtools_version = match check_version(SAMTOOLS_TAG).await {
        Ok(version) => {
            eprintln!("{}", version);
            version
        }
        Err(err) => {
            return Err(anyhow!("Error checking samtools version: {}", err));
        }
    };

    //*****************
    // Cleanup, hanging tasks
    


    fastp_stream_task.await??;
    fastp_err_task.await??;
    let fastp_status = fastp_child.wait().await?;
    if fastp_status.success() {
        eprintln!("Fastp exited successfully");
    }
    else {
        return Err(anyhow!("Fastp exited with non-zero status: {}", fastp_status));
    }
    
    //T_junction completion
    val_done_rx.await??;
    eprintln!("Validation t_junction done");

    // Ensure Minimap2 FIFO write tasks complete
    ref_write_task.await??;
    query_write_task.await??;
    eprintln!("Minimap2 fifo tasks done");


    minimap2_err_task.await??;
    minimap2_write_task.await??;
    minimap2_task.await??;
    let minimap2_status = minimap2_child.wait().await?;
    if minimap2_status.success() {
        eprintln!("Minimap2 exited successfully");
    }
    else {
        return Err(anyhow!("Minimap2 exited with non-zero status: {}", minimap2_status));
    }
    
    
    if ref_pipe_path.exists() {
        std::fs::remove_file(&ref_pipe_path)?;
    }
    if query_pipe_path.exists() {
        std::fs::remove_file(&query_pipe_path)?;
    }

    pigz_write_task.await??;
    let pigz_status = pigz_child.wait().await?;
    if pigz_status.success() {
        eprintln!("Pigz exited successfully");
    }
    else {
        return Err(anyhow!("pigz exited with non-zero status: {}", pigz_status));
    }
    
    eprintln!("Finished generating consensus genome");
    
    Ok(())
}