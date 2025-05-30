use std::collections::HashMap;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::{Arguments, Technology};
use std::process::Command;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::command::{generate_cli, check_versions};
use crate::utils::file::{file_path_manipulator, write_parse_output_to_temp};
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, write_fasta_to_fifo};
use crate::utils::db::write_hdf5_seq_to_fifo;
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::db::{lookup_sequence, load_index, build_new_in_memory_index};


const ERCC_FASTA: &str = "ercc_sequences.fasta";

pub async fn run(args: &Arguments) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let cwd = std::env::current_dir()?;
    let verbose = args.verbose.clone();

    //External tools check
    let _tool_versions = check_versions(vec![SAMTOOLS_TAG, MINIMAP2_TAG, FASTP_TAG, SAMTOOLS_TAG]).await?;

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
            let file2_full_path = file_path_manipulator(&PathBuf::from(file), &cwd.clone(), None, None, "");
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

    let ercc_path = file_path_manipulator(&PathBuf::from(ERCC_FASTA), &cwd, None, None, "");
    if !ercc_path.exists() {
        return Err(anyhow!("Specficied ercc {:?} does not exist.", ercc_path));
    }

    //*****************
    // Input Validation
    
    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(&sample_base), &cwd.clone(), None, Some("validated"), "_");
    let rx = read_and_interleave_sequences(file1_path, file2_path, Some(technology.clone()), args.max_reads, args.min_read_len, args.max_read_len)?;
    let val_rx_stream = ReceiverStream::new(rx);
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
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
    let val_fastp_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing fastp stream"))?;
    let val_pigz_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing file stream"))?;

    //Pigz stream to intermediate file output
    let val_pigz_args = generate_cli(PIGZ_TAG, &args, None)?;
    let (mut val_pigz_child, _val_pigz_stream_task) = stream_to_cmd(val_pigz_stream, PIGZ_TAG, val_pigz_args, StreamDataType::IlluminaFastq).await?;

    let val_pigz_out_stream = parse_child_output(
        &mut val_pigz_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size,
    ).await?;
    let val_pigz_write_task = tokio::spawn(stream_to_file(
        val_pigz_out_stream,
        validated_interleaved_file_path,
    ));


    // Fastp stream
    let val_fastp_args = generate_cli(FASTP_TAG, &args, None)?;
    let (mut val_fastp_child, val_fastp_stream_task) = stream_to_cmd(val_fastp_stream, FASTP_TAG, val_fastp_args, StreamDataType::IlluminaFastq).await?;
    let val_fastp_out_stream = parse_child_output(
        &mut val_fastp_child,
        ChildStream::Stdout,
        ParseMode::Bytes, // Use Bytes to avoid parsing issues
        args.buffer_size / 4, // Reduce buffer size
    ).await?;
    let mut val_fastp_out_stream = ReceiverStream::new(val_fastp_out_stream);
    
    let val_fastp_err_stream = parse_child_output(
        &mut val_fastp_child,
        ChildStream::Stderr,
        ParseMode::Bytes,
        args.buffer_size,
    ).await?;
    let val_fastp_err_task = tokio::spawn(async move {
        let mut err_stream = ReceiverStream::new(val_fastp_err_stream);
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
    //Fetch reference
    
    // Retrieve index or create it for host sequence and filter sequence
    let h5_index = if let Some(index_file) = &args.ref_index {
        let index_full_path = file_path_manipulator(&PathBuf::from(index_file), &cwd.clone(), None, None, "");
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



    //*****************
    //Host Removal
    
    // Create FIFO pipes
    let host_ref_temp = NamedTempFile::new()?;
    let host_ref_pipe_path = host_ref_temp.path().to_path_buf();
    
    if host_ref_pipe_path.exists() {
        std::fs::remove_file(&host_ref_pipe_path)?;
    }
    Command::new("mkfifo")
        .arg(&host_ref_pipe_path)
        .status()?;
    
    
    let host_accession = args.host_accession.clone();
    let host_sequence = args.host_sequence.clone();

    // If the host sequence file is given, load it, if not retrieve it by accession from ref_db
    let host_seq = match &host_sequence {
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
    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(val_fastp_out_stream, None).await?;
  

    // Create FIFO pipe from either the host_sequence or host_accession
    let host_ref_write_task = tokio::spawn({
        let cwd = cwd.clone();
        let host_ref_pipe_path = host_ref_pipe_path.clone();
        async move {
            
            match host_seq {
                Some(seq) => {
                    match &host_accession {
                        Some(accession) => {
                            write_hdf5_seq_to_fifo(seq, &accession, &host_ref_pipe_path).await
                        }
                        None => {
                            return Err(anyhow!("Must provide either a host sequence file with --host_sequence or an accession with --host_accession"))
                        }
                    }
                }
                None => {
                    match &host_sequence {
                        Some(host_sequence_file) => {
                            let host_sequence_path = file_path_manipulator(&PathBuf::from(host_sequence_file), &cwd.clone(), None, None, "");
                            tokio::task::spawn_blocking(move || {
                                write_fasta_to_fifo(&host_sequence_path, &host_ref_pipe_path)
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
    
    
    let host_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(host_ref_pipe_path.clone(), host_query_pipe_path.clone())))?;
    
    let (mut host_minimap2_child, host_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, host_minimap2_args, args.verbose).await?;
    let host_minimap2_out_stream = parse_child_output(
        &mut host_minimap2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;

    
    let host_samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([("-f".to_string(), Some("4".to_string())),]), // Require unmapped reads (SAM flag 4). Pass only unmapped reads here
    };
    let host_samtools_args_view = generate_cli(
        SAMTOOLS_TAG,
        &args,
        Some(&host_samtools_config_view),
    )?;

    let (mut host_samtools_child_view, host_samtools_task_view) = stream_to_cmd(host_minimap2_out_stream, SAMTOOLS_TAG, host_samtools_args_view, StreamDataType::JustBytes).await?;
    let host_samtools_out_stream_view = parse_child_output(
        &mut host_samtools_child_view,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    // let host_samtools_out_stream_view = ReceiverStream::new(host_samtools_out_stream_view);


    // //Output to FASTQ through samtools
    let host_samtools_config_fastq = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let host_samtools_args_fastq = generate_cli(
        SAMTOOLS_TAG,
        &args,
        Some(&host_samtools_config_fastq),
    )?;

    let (mut host_samtools_child_fastq, host_samtools_task_fastq) = stream_to_cmd(host_samtools_out_stream_view, SAMTOOLS_TAG, host_samtools_args_fastq, StreamDataType::JustBytes).await?;
    let host_samtools_out_stream_fastq = parse_child_output(
        &mut host_samtools_child_fastq,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    let mut host_samtools_out_stream_fastq = ReceiverStream::new(host_samtools_out_stream_fastq);
    
    let host_samtools_child_fastq_err_stream = parse_child_output(
        &mut host_samtools_child_fastq,
        ChildStream::Stderr,
        ParseMode::Bytes,
        args.buffer_size,
    ).await?;
    let host_samtools_child_fastq_err_task = tokio::spawn(async move {
        let mut err_stream = ReceiverStream::new(host_samtools_child_fastq_err_stream);
        if verbose {
            while let Some(ParseOutput::Bytes(chunk)) = err_stream.next().await {
                eprintln!("samtools fastq stderr: {}", String::from_utf8_lossy(&chunk));
            }
        } else {
            while err_stream.next().await.is_some() {}
        }
        
        Ok::<(), anyhow::Error>(())
    });

    // Split for file write and passing on to next stage
    let no_host_file_path = file_path_manipulator(&PathBuf::from(&sample_base), &cwd.clone(), None, Some("no_host"), "_");
    
    let (host_streams, host_done_rx) = t_junction(
        host_samtools_out_stream_fastq,
        2,
        args.buffer_size,
        args.stall_threshold.try_into().unwrap(),
        Some(args.stream_sleep_ms),
        50,
    )
        .await?;
    
    if host_streams.len() != 2 {
        return Err(anyhow!("Expected exactly 2 streams, got {}", host_streams.len()));
    }
    
    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing output stream"))?;
    let no_host_file_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing file stream"))?;

    let host_samtools_write_task = tokio::spawn(stream_to_file(
        no_host_file_stream,
        PathBuf::from(no_host_file_path),
    ));
    
    let no_host_output_stream = ReceiverStream::new(no_host_output_stream);
    
    
    //*****************
    // Split by Technology

    match technology {
        Technology::Illumina => {
            eprintln!("Illumina");
            //*****************
            // ERCC

            let (ercc_query_write_task, ercc_query_pipe_path) = write_parse_output_to_temp(no_host_output_stream, None).await?;
        
            let ercc_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(ercc_path, ercc_query_pipe_path.clone())))?;
        
            let (mut ercc_minimap2_child, ercc_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, ercc_minimap2_args, args.verbose).await?;
            let ercc_minimap2_out_stream = parse_child_output(
                &mut ercc_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;
        
            let ercc_minimap_write_task = tokio::spawn(stream_to_file(
                ercc_minimap2_out_stream,
                PathBuf::from("test_samtools_ercc.sam"),
            ));
        
        ercc_query_write_task.await??;
        ercc_minimap_write_task.await??;
        ercc_minimap2_err_task.await??;

        } // end tech illumina
        Technology::ONT => {
            return Err(anyhow!("Minion not ready"));
        }

    }

    //*****************
    // Cleanup, hanging tasks

    val_fastp_stream_task.await??;
    val_fastp_err_task.await??;
    let val_fastp_status = val_fastp_child.wait().await?;
    if val_fastp_status.success() {
        eprintln!("Fastp exited successfully");
    }
    else {
        return Err(anyhow!("Fastp exited with non-zero status: {}", val_fastp_status));
    }
    
    //T_junction completion
    val_done_rx.await??;
    eprintln!("Validation t_junction done");

    // Ensure Minimap2 FIFO write tasks complete
    host_ref_write_task.await??;
    host_query_write_task.await??;
    eprintln!("Minimap2 fifo tasks done");





    let minimap2_status = host_minimap2_child.wait().await?;
    if minimap2_status.success() {
        eprintln!("Minimap2 exited successfully");
    }
    else {
        return Err(anyhow!("Minimap2 exited with non-zero status: {}", minimap2_status));
    }
    
    host_samtools_task_view.await??;
    let samtools_status = host_samtools_child_view.wait().await?;
    if samtools_status.success() {
        eprintln!("Samtools exited successfully");
    }
    else {
        return Err(anyhow!("Samtools exited with non-zero status: {}", samtools_status));
    }

    if host_ref_pipe_path.exists() {
        std::fs::remove_file(&host_ref_pipe_path)?;
    }



    

    val_pigz_write_task.await??;
    let pigz_status = val_pigz_child.wait().await?;
    if pigz_status.success() {
        eprintln!("Pigz exited successfully");
    }
    else {
        return Err(anyhow!("pigz exited with non-zero status: {}", pigz_status));
    }

    host_samtools_task_fastq.await??;
    let host_samtools_task_fastq_status = host_samtools_child_fastq.wait().await?;
    if  host_samtools_task_fastq_status.success() {
        eprintln!("Samtools fastq exited successfully");
    }
    else {
        return Err(anyhow!("Samtools fastq exited with non-zero status: {}", host_samtools_task_fastq_status));
    }

    host_samtools_write_task.await??;
    host_samtools_child_fastq_err_task.await??;
    host_minimap2_err_task.await??;
    host_done_rx.await??;
    eprintln!("Host t_junction done");
    
    eprintln!("Finished generating consensus genome");
    
    Ok(())
}