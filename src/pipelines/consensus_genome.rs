use std::collections::HashMap;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::{Arguments, Technology};
use std::process::Command;
use std::time::Instant;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::command::{generate_cli, check_versions};
use crate::utils::file::{extension_remover, file_path_manipulator, write_parse_output_to_temp};
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, write_fasta_to_fifo};
use crate::utils::db::write_hdf5_seq_to_fifo;
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand, KRAKEN2_TAG};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::db::{lookup_sequence, load_index, build_new_in_memory_index, get_index, retrieve_h5_seq};


const ERCC_FASTA: &str = "ercc_sequences.fasta";

pub async fn run(args: &Arguments) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let cwd = std::env::current_dir()?;
    
    //External tools check
    let _tool_versions = check_versions(vec![SAMTOOLS_TAG, MINIMAP2_TAG, FASTP_TAG, SAMTOOLS_TAG, KRAKEN2_TAG]).await?;

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

    let sample_base_buf: PathBuf = PathBuf::from(&sample_base);
    let (no_ext_sample_base_buf, _) = extension_remover(&sample_base_buf);
    let no_ext_sample_base = no_ext_sample_base_buf.to_string_lossy().into_owned();


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
    let (mut val_pigz_child, val_pigz_stream_task, val_pigz_err_task) = stream_to_cmd(val_pigz_stream, PIGZ_TAG, val_pigz_args, StreamDataType::IlluminaFastq, args.verbose).await?;

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
    let (mut val_fastp_child, val_fastp_stream_task, val_fastp_err_task) = stream_to_cmd(val_fastp_stream, FASTP_TAG, val_fastp_args, StreamDataType::IlluminaFastq, args.verbose).await?;
    let val_fastp_out_stream = parse_child_output(
        &mut val_fastp_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    let mut val_fastp_out_stream = ReceiverStream::new(val_fastp_out_stream);
    

    
    //*****************
    //Fetch reference
    let index_start = Instant::now();
    // Retrieve index or create it for host sequence and filter sequence
    let h5_index = get_index(&args).await?;
    println!("Index retrieve time: {} milliseconds.", index_start.elapsed().as_millis());


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
    
    let (host_accession, host_seq) = retrieve_h5_seq(args.host_accession.clone(), args.host_sequence.clone(), Some(&ref_db_path), Some(&h5_index)).await?;
    // Create FIFO pipe from either the host_sequence or host_accession
    let host_ref_write_task = tokio::spawn({
        let host_ref_pipe_path = host_ref_pipe_path.clone();
        async move {
            write_hdf5_seq_to_fifo(&host_seq, &host_accession, &host_ref_pipe_path).await;
        }
    });
    
    // Create FIFO pipe for the fastp output to stream to minimap2
    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(val_fastp_out_stream, None).await?;
    
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

    let (mut host_samtools_child_view, host_samtools_task_view, host_samtools_err_task_view) = stream_to_cmd(host_minimap2_out_stream, SAMTOOLS_TAG, host_samtools_args_view, StreamDataType::JustBytes, args.verbose).await?;
    let host_samtools_out_stream_view = parse_child_output(
        &mut host_samtools_child_view,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    
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

    let (mut host_samtools_child_fastq, host_samtools_task_fastq, host_samtools_err_task_fastq) = stream_to_cmd(host_samtools_out_stream_view, SAMTOOLS_TAG, host_samtools_args_fastq, StreamDataType::JustBytes, args.verbose).await?;
    let host_samtools_out_stream_fastq = parse_child_output(
        &mut host_samtools_child_fastq,
        ChildStream::Stdout,
        ParseMode::Bytes,
        args.buffer_size / 4,
    ).await?;
    let mut host_samtools_out_stream_fastq = ReceiverStream::new(host_samtools_out_stream_fastq);


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

            let (ercc_streams, ercc_done_rx) = t_junction(
                no_host_output_stream,
                2,
                args.buffer_size,
                args.stall_threshold.try_into().unwrap(),
                Some(args.stream_sleep_ms),
                50,
            )
                .await?;

            if ercc_streams.len() != 2 {
                return Err(anyhow!("Expected exactly 2 streams, got {}", ercc_streams.len()));
            }

            let mut streams_iter = ercc_streams.into_iter();
            let ercc_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing ercc stream"))?;
            let ercc_bypass_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing ercc bypass stream"))?;
            let mut ercc_stream = ReceiverStream::new(ercc_stream);



            let (ercc_query_write_task, ercc_query_pipe_path) = write_parse_output_to_temp(ercc_stream, None).await?;
            let ercc_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(ercc_path, ercc_query_pipe_path.clone())))?;

            let (mut ercc_minimap2_child, ercc_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, ercc_minimap2_args, args.verbose).await?;
            let ercc_minimap2_out_stream = parse_child_output(
                &mut ercc_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;
            

            let ercc_samtools_config_view = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::View,
                subcommand_fields: HashMap::from([]),
            };
            let ercc_samtools_args_view = generate_cli(
                SAMTOOLS_TAG,
                &args,
                Some(&ercc_samtools_config_view),
            )?;

            let (mut ercc_samtools_child_view, ercc_samtools_task_view, _ercc_samtools_err_task_view) = stream_to_cmd(ercc_minimap2_out_stream, SAMTOOLS_TAG, ercc_samtools_args_view, StreamDataType::JustBytes, args.verbose).await?;
            let ercc_samtools_out_stream_view = parse_child_output(
                &mut ercc_samtools_child_view,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;
            
            
            let ercc_samtools_config_stats = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::Stats,
                subcommand_fields: HashMap::from([("-".to_string(), None)]),
            };
            let ercc_samtools_args_view = generate_cli(
                SAMTOOLS_TAG,
                &args,
                Some(&ercc_samtools_config_stats),
            )?;
            
            let ercc_stats_file_path = no_ext_sample_base + "_stats.txt";
            let (mut ercc_samtools_child_stats, ercc_samtools_task_stats, _ercc_samtools_err_task_stats) = stream_to_cmd(ercc_samtools_out_stream_view, SAMTOOLS_TAG, ercc_samtools_args_view, StreamDataType::JustBytes, args.verbose).await?;

            let ercc_samtools_out_stream_stats = parse_child_output(
                &mut ercc_samtools_child_stats,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;

            let ercc_stats_write_task = tokio::spawn(stream_to_file(
                ercc_samtools_out_stream_stats,
                PathBuf::from(ercc_stats_file_path),
            ));



            let (filter_align_accession, filter_align_seq) = retrieve_h5_seq(args.ref_accession.clone(), args.ref_sequence.clone(), Some(&ref_db_path), Some(&h5_index)).await?;

            //*****************
            // Filter Reads

            let mut filter_reads_out_stream: ReceiverStream<ParseOutput>;
            
            if args.dont_filter_reads {
                filter_reads_out_stream = ReceiverStream::new(ercc_bypass_stream);
            }
            
            else {
                let mut ercc_bypass_stream = ReceiverStream::new(ercc_bypass_stream);
                let filter_ref_temp = NamedTempFile::new()?;
                let filter_ref_pipe_path = filter_ref_temp.path().to_path_buf();
                
                if filter_ref_pipe_path.exists() {
                    std::fs::remove_file(&filter_ref_pipe_path)?;
                }
                Command::new("mkfifo")
                    .arg(&filter_ref_pipe_path)
                    .status()?;


                let filter_ref_write_task = tokio::spawn({
                    let filter_ref_pipe_path = filter_ref_pipe_path.clone();
                    async move {
                        write_hdf5_seq_to_fifo(&filter_align_seq, &filter_align_accession, &filter_ref_pipe_path).await;
                    }
                });
                let (filter_query_write_task, filter_query_pipe_path) = write_parse_output_to_temp(ercc_bypass_stream, None).await?;

                let filter_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(filter_ref_pipe_path.clone(), filter_query_pipe_path.clone())))?;
                let (mut filter_minimap2_child, filter_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, filter_minimap2_args, args.verbose).await?;
                let filter_minimap2_out_stream = parse_child_output(
                    &mut filter_minimap2_child,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    args.buffer_size / 4,
                ).await?;

                
                //Convert to FASTQ

                let filter_samtools_config_fastq = SamtoolsConfig {
                    subcommand: SamtoolsSubcommand::Fastq,
                    subcommand_fields: HashMap::from([("-".to_string(), None)]),
                };
                let filter_samtools_args_fastq = generate_cli(
                    SAMTOOLS_TAG,
                    &args,
                    Some(&filter_samtools_config_fastq),
                )?;

                let (mut filter_samtools_child_fastq, filter_samtools_task_fastq, filter_samtools_err_task_fastq) = stream_to_cmd(filter_minimap2_out_stream, SAMTOOLS_TAG, filter_samtools_args_fastq, StreamDataType::JustBytes, args.verbose).await?;
                let filter_samtools_out_stream_fastq = parse_child_output(
                    &mut filter_samtools_child_fastq,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    args.buffer_size / 4,
                ).await?;
                // let mut filter_samtools_out_stream_fastq = ReceiverStream::new(filter_samtools_out_stream_fastq);


                // let kraken2_report_path = file_path_manipulator(&PathBuf::from(&no_ext_sample_base_buf), &cwd.clone(), None, Some("kraken2_report.txt"), "_");
                // let kraken2_classified_path = file_path_manipulator(&PathBuf::from(&no_ext_sample_base_buf), &cwd.clone(), None, Some("classified.fq"), "_");
                // let (kraken2_query_write_task, kraken2_query_pipe_path) = write_parse_output_to_temp(ercc_bypass_stream, None).await?;
                // 
                // let filter_reads_kraken2_config = Kraken2Config {
                //     report_path: kraken2_report_path,
                //     classified_path: kraken2_classified_path,
                //     fastq_path: kraken2_query_pipe_path
                // };
                // 
                // let filter_reads_kraken2_args = generate_cli(KRAKEN2_TAG, &args, Some(&filter_reads_kraken2_config))?;
                // eprintln!("Filter reads kraken2 args: {:?}", filter_reads_kraken2_args);
                // 
                // 
                // let (mut filter_kraken2_child, _filter_kraken2_err_task) = spawn_cmd(MINIMAP2_TAG, filter_reads_kraken2_args, args.verbose).await?;
                // let filter_kraken2_out_stream = parse_child_output(
                //     &mut filter_kraken2_child,
                //     ChildStream::Stdout,
                //     ParseMode::Bytes,
                //     args.buffer_size / 4,
                // ).await?;
                

                filter_query_write_task.await??;
                filter_ref_write_task.await?;
                // filter_output_write_task.await?;

                filter_reads_out_stream = ReceiverStream::new(filter_samtools_out_stream_fastq);
            }

            

            let filter_output_write_task = tokio::spawn(stream_to_file(
                filter_reads_out_stream.into_inner(),
                PathBuf::from("test_samtools_filter.fastq"),
            ));

            //*****************
            // Align Reads to Target

            // let align_ref_temp = NamedTempFile::new()?;
            // let align_ref_pipe_path = align_ref_temp.path().to_path_buf();
            // 
            // if align_ref_pipe_path.exists() {
            //     std::fs::remove_file(&align_ref_pipe_path)?;
            // }
            // Command::new("mkfifo")
            //     .arg(&align_ref_pipe_path)
            //     .status()?;
            // 
            // 
            // let align_ref_write_task = tokio::spawn({
            //     let align_ref_pipe_path = align_ref_pipe_path.clone();
            //     async move {
            //         write_hdf5_seq_to_fifo(&filter_align_seq, &filter_align_accession, &filter_ref_pipe_path).await;
            //     }
            // });
            // 
            

            // let (filter_query_write_task, filter_query_pipe_path) = write_parse_output_to_temp(ercc_bypass_stream, None).await?;
            // 
            // let filter_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(filter_ref_pipe_path.clone(), filter_query_pipe_path.clone())))?;
            // let (mut filter_minimap2_child, filter_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, filter_minimap2_args, args.verbose).await?;
            // let filter_minimap2_out_stream = parse_child_output(
            //     &mut filter_minimap2_child,
            //     ChildStream::Stdout,
            //     ParseMode::Bytes,
            //     args.buffer_size / 4,
            // ).await?;
            // 





            ercc_query_write_task.await??;
            // ercc_output_write_task.await??;
            ercc_minimap2_err_task.await??;
            ercc_samtools_task_view.await??;
            ercc_samtools_task_stats.await??;
            ercc_stats_write_task.await??;
            ercc_done_rx.await??;

            
            
            
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
    host_ref_write_task.await?;
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
    host_samtools_err_task_view.await??;
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




    val_pigz_stream_task.await??;
    val_pigz_err_task.await??;
    val_pigz_write_task.await??;
    let pigz_status = val_pigz_child.wait().await?;
    if pigz_status.success() {
        eprintln!("Pigz exited successfully");
    }
    else {
        return Err(anyhow!("pigz exited with non-zero status: {}", pigz_status));
    }

    host_samtools_task_fastq.await??;
    host_samtools_err_task_fastq.await??;

    let host_samtools_task_fastq_status = host_samtools_child_fastq.wait().await?;
    if  host_samtools_task_fastq_status.success() {
        eprintln!("Samtools fastq exited successfully");
    }
    else {
        return Err(anyhow!("Samtools fastq exited with non-zero status: {}", host_samtools_task_fastq_status));
    }

    host_samtools_write_task.await??;

    host_minimap2_err_task.await??;
    host_done_rx.await??;
    eprintln!("Host t_junction done");
    
    eprintln!("Finished generating consensus genome");
    
    Ok(())
}