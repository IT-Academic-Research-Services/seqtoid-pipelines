use std::collections::HashMap;
use std::sync::Arc;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::{Arguments, Technology, TargetType};
use std::process::Command;
use std::time::Instant;
use tokio::fs::File as TokioFile;
use tokio::io::AsyncWriteExt;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::command::{generate_cli, check_versions, check_version};
use crate::utils::file::{extension_remover, file_path_manipulator, write_parse_output_to_temp};
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base, write_fasta_to_fifo, parse_and_filter_fastq_id};
use crate::utils::db::write_hdf5_seq_to_fifo;
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd, parse_bytes, parse_fastq};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand, KRAKEN2_TAG, BCFTOOLS_TAG, BcftoolsSubcommand, IVAR_TAG, IvarSubcommand};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::command::ivar::IvarConfig;
use crate::utils::db::{lookup_sequence, load_index, build_new_in_memory_index, get_index, retrieve_h5_seq};
use tokio::sync::mpsc;
use tokio::time::{timeout, Duration};
use tokio::io::BufWriter;
use futures::future::try_join_all;


const ERCC_FASTA: &str = "ercc_sequences.fasta";

pub async fn run(args: &Arguments) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let cwd = std::env::current_dir()?;

    // Initialize cleanup tasks
    let mut cleanup_tasks: Vec<tokio::task::JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<tokio::sync::oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();
    
    //External tools check
    let mut tool_versions = check_versions(vec![SAMTOOLS_TAG, MINIMAP2_TAG, FASTP_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG]).await?;

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
    cleanup_receivers.push(val_done_rx);

    let mut streams_iter = val_streams.into_iter();
    let val_fastp_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing fastp stream"))?;
    let val_pigz_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing file stream"))?;

    //Pigz stream to intermediate file output
    let val_pigz_args = generate_cli(PIGZ_TAG, &args, None)?;
    let (mut val_pigz_child, val_pigz_stream_task, val_pigz_err_task) = stream_to_cmd(val_pigz_stream, PIGZ_TAG, val_pigz_args, StreamDataType::IlluminaFastq, args.verbose).await?;
    cleanup_tasks.push(val_pigz_stream_task);
    cleanup_tasks.push(val_pigz_err_task);
    
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
    cleanup_tasks.push(val_pigz_write_task);


    // Fastp stream
    let val_fastp_args = generate_cli(FASTP_TAG, &args, None)?;
    let (mut val_fastp_child, val_fastp_stream_task, val_fastp_err_task) = stream_to_cmd(val_fastp_stream, FASTP_TAG, val_fastp_args, StreamDataType::IlluminaFastq, args.verbose).await?;
    cleanup_tasks.push(val_fastp_stream_task);
    cleanup_tasks.push(val_fastp_err_task);
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
    let host_ref_write_task = tokio::spawn({
        let host_ref_pipe_path = host_ref_pipe_path.clone();
        async move {
            let result = write_hdf5_seq_to_fifo(&host_seq, &host_accession, &host_ref_pipe_path).await;
            result
        }
    });
    cleanup_tasks.push(host_ref_write_task);
    
    // Create FIFO pipe for the fastp output to stream to minimap2
    let (host_query_write_task, host_query_pipe_path) = write_parse_output_to_temp(val_fastp_out_stream, None).await?;
    cleanup_tasks.push(host_query_write_task);
    let host_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(host_ref_pipe_path.clone(), host_query_pipe_path.clone())))?;
    let (mut host_minimap2_child, host_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, host_minimap2_args, args.verbose).await?;
    cleanup_tasks.push(host_minimap2_err_task);
    
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
    cleanup_tasks.push(host_samtools_task_view);
    cleanup_tasks.push(host_samtools_err_task_view);
    
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
    cleanup_tasks.push(host_samtools_task_fastq);
    cleanup_tasks.push(host_samtools_err_task_fastq);
    
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
    cleanup_receivers.push(host_done_rx);
    
    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing output stream"))?;
    let no_host_file_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing file stream"))?;
    
    let host_samtools_write_task = tokio::spawn(stream_to_file(
        no_host_file_stream,
        PathBuf::from(no_host_file_path),
    ));
    cleanup_tasks.push(host_samtools_write_task);

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
            cleanup_tasks.push(ercc_query_write_task);
            let ercc_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(ercc_path, ercc_query_pipe_path.clone())))?;
    
            let (mut ercc_minimap2_child, ercc_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, ercc_minimap2_args, args.verbose).await?;
            let ercc_minimap2_out_stream = parse_child_output(
                &mut ercc_minimap2_child,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_minimap2_err_task);
    
            let ercc_samtools_config_view = SamtoolsConfig {
                subcommand: SamtoolsSubcommand::View,
                subcommand_fields: HashMap::from([]),
            };
            let ercc_samtools_args_view = generate_cli(
                SAMTOOLS_TAG,
                &args,
                Some(&ercc_samtools_config_view),
            )?;
    
            let (mut ercc_samtools_child_view, ercc_samtools_task_view, ercc_samtools_err_task_view) = stream_to_cmd(ercc_minimap2_out_stream, SAMTOOLS_TAG, ercc_samtools_args_view, StreamDataType::JustBytes, args.verbose).await?;
            let ercc_samtools_out_stream_view = parse_child_output(
                &mut ercc_samtools_child_view,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_samtools_task_view);
            cleanup_tasks.push(ercc_samtools_err_task_view);
            
            
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
            let (mut ercc_samtools_child_stats, ercc_samtools_task_stats, ercc_samtools_err_task_stats) = stream_to_cmd(ercc_samtools_out_stream_view, SAMTOOLS_TAG, ercc_samtools_args_view, StreamDataType::JustBytes, args.verbose).await?;
    
            let ercc_samtools_out_stream_stats = parse_child_output(
                &mut ercc_samtools_child_stats,
                ChildStream::Stdout,
                ParseMode::Bytes,
                args.buffer_size / 4,
            ).await?;
            cleanup_tasks.push(ercc_samtools_task_stats);
            cleanup_tasks.push(ercc_samtools_err_task_stats);
    
            let ercc_stats_write_task = tokio::spawn(stream_to_file(
                ercc_samtools_out_stream_stats,
                PathBuf::from(ercc_stats_file_path),
            ));
            cleanup_tasks.push(ercc_stats_write_task);
    
    
            //*****************
            // Filter Reads
            
            let (filter_align_accession, filter_align_seq) = retrieve_h5_seq(args.ref_accession.clone(), args.ref_sequence.clone(), Some(&ref_db_path), Some(&h5_index)).await?;
            let filter_align_seq = Arc::new(filter_align_seq);
            let filter_align_accession = Arc::new(filter_align_accession);
            eprintln!("Align accession: {:?}", filter_align_accession);
            
            let mut filter_reads_out_stream: ReceiverStream<ParseOutput>;
            
            if args.dont_filter_reads {
                filter_reads_out_stream = ReceiverStream::new(ercc_bypass_stream);
            }
            // 
            else {
                eprintln!("Filtering");
                let mut ercc_bypass_stream = ReceiverStream::new(ercc_bypass_stream);
                let filter_ref_temp = NamedTempFile::new()?;
                let filter_ref_pipe_path = filter_ref_temp.path().to_path_buf();
            
                if filter_ref_pipe_path.exists() {
                    std::fs::remove_file(&filter_ref_pipe_path)?;
                }
                Command::new("mkfifo")
                    .arg(&filter_ref_pipe_path)
                    .status()?;
            
                let filter_align_seq_clone = Arc::clone(&filter_align_seq);
                let filter_align_accession_clone = Arc::clone(&filter_align_accession);
                eprintln!("Filter accession: {:?}", filter_align_accession_clone);
                let filter_ref_write_task = tokio::spawn({
                    let filter_ref_pipe_path = filter_ref_pipe_path.clone();
                    async move {
                        let result = write_hdf5_seq_to_fifo(&filter_align_seq_clone, &filter_align_accession_clone, &filter_ref_pipe_path).await;
                        result
                    }
                });
                cleanup_tasks.push(filter_ref_write_task);
                
                let (filter_query_write_task, filter_query_pipe_path) = write_parse_output_to_temp(ercc_bypass_stream, None).await?;
                cleanup_tasks.push(filter_query_write_task);
                
                let filter_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(filter_ref_pipe_path.clone(), filter_query_pipe_path.clone())))?;
                let (mut filter_minimap2_child, filter_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, filter_minimap2_args, args.verbose).await?;
                let filter_minimap2_out_stream = parse_child_output(
                    &mut filter_minimap2_child,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    args.buffer_size / 4,
                ).await?;
                cleanup_tasks.push(filter_minimap2_err_task);

                
                
                // Sort output
                let filter_samtools_config_sort = SamtoolsConfig {
                    subcommand: SamtoolsSubcommand::Sort,
                    subcommand_fields: HashMap::from([("-n".to_string(), None), ("-O".to_string(), Some("BAM".to_string())), ("-".to_string(), None)]),
                };
                let filter_samtools_args_sort = generate_cli(
                    SAMTOOLS_TAG,
                    &args,
                    Some(&filter_samtools_config_sort),
                )?;
                eprintln!("Filter samtools: {:?}", filter_samtools_args_sort);
                let (mut filter_samtools_child_sort, filter_samtools_task_sort, filter_samtools_err_task_sort) = stream_to_cmd(filter_minimap2_out_stream, SAMTOOLS_TAG, filter_samtools_args_sort, StreamDataType::JustBytes, args.verbose).await?;
                cleanup_tasks.push(filter_samtools_task_sort);
                cleanup_tasks.push(filter_samtools_err_task_sort);
                let filter_samtools_out_stream_sort = parse_child_output(
                    &mut filter_samtools_child_sort,
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
                
                let (mut filter_samtools_child_fastq, filter_samtools_task_fastq, filter_samtools_err_task_fastq) = stream_to_cmd(filter_samtools_out_stream_sort, SAMTOOLS_TAG, filter_samtools_args_fastq, StreamDataType::JustBytes, args.verbose).await?;
                cleanup_tasks.push(filter_samtools_task_fastq);
                cleanup_tasks.push(filter_samtools_err_task_fastq);
                let filter_samtools_out_stream_fastq = parse_child_output(
                    &mut filter_samtools_child_fastq,
                    ChildStream::Stdout,
                    ParseMode::Bytes,
                    args.buffer_size / 4,
                ).await?;


                // let test_write_task = tokio::spawn(stream_to_file(
                //     filter_samtools_out_stream_fastq,
                //     PathBuf::from("test_filtered.fq"),
                // ));
                // test_write_task.await??;
            
                let mut filter_samtools_out_stream_fastq = ReceiverStream::new(filter_samtools_out_stream_fastq);
            
                let kraken2_report_path = file_path_manipulator(&PathBuf::from(&no_ext_sample_base_buf), &cwd.clone(), None, Some("kraken2_report.txt"), "_");
                let kraken2_classified_temp = NamedTempFile::new()?;
                let kraken2_classified_pipe_path = kraken2_classified_temp.path().to_path_buf();
            
                if kraken2_classified_pipe_path.exists() {
                    std::fs::remove_file(&kraken2_classified_pipe_path)?;
                }
                Command::new("mkfifo")
                    .arg(&kraken2_classified_pipe_path)
                    .status()?
                    .success()
                    .then(|| ())
                    .ok_or_else(|| anyhow!("Failed to create mkfifo for Kraken2 classified output"))?;
            
                // Create a temporary FIFO for kraken2 input
                let (kraken2_query_write_task, kraken2_query_pipe_path) = write_parse_output_to_temp(filter_samtools_out_stream_fastq, None).await?;
                cleanup_tasks.push(kraken2_query_write_task);
            
                let filter_reads_kraken2_config = Kraken2Config {
                    report_path: kraken2_report_path,
                    classified_path: kraken2_classified_pipe_path.clone(),
                    fastq_path: kraken2_query_pipe_path.clone(),
                };
            
                let filter_reads_kraken2_args = generate_cli(KRAKEN2_TAG, &args, Some(&filter_reads_kraken2_config))?;
                eprintln!("kraken2 args: {:?}", filter_reads_kraken2_args);
                let (mut filter_kraken2_child, filter_kraken2_err_task) = spawn_cmd(KRAKEN2_TAG, filter_reads_kraken2_args, args.verbose).await?;
                cleanup_tasks.push(filter_kraken2_err_task);
                
                let kraken2_classified_stream = TokioFile::open(&kraken2_classified_pipe_path).await?;
                eprintln!("Opened kraken2 classified pipe: {}", kraken2_classified_pipe_path.display());
                
                // let (parse_rx, parse_task) = {
                //     let buffer_size = args.buffer_size;
                //     let (tx, rx) = mpsc::channel(buffer_size);
                //     let task = tokio::spawn(async move {
                //         parse_fastq(kraken2_classified_stream, buffer_size)
                //             .await
                //             .map(|_| ())
                //             .map_err(|e| anyhow::anyhow!("parse_fastq failed: {}", e)) // Convert error to anyhow::Error
                //     });
                //     (rx, task)
                // };
                // cleanup_tasks.push(parse_task);
            
                let parse_rx = parse_fastq(kraken2_classified_stream, args.buffer_size).await?; // KEEP THIS until you test the above
                let filter_fn = |id: &str| id.contains("kraken:taxid|2697049"); // TODO no hard code
            
                let (filtered_rx, filter_task) = parse_and_filter_fastq_id(parse_rx, args.buffer_size, filter_fn.clone());
                cleanup_tasks.push(filter_task);
                
                let (parse_output_tx, parse_output_rx) = mpsc::channel(args.buffer_size);
                tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(filtered_rx);
                    while let Some(record) = stream.next().await {
                        if parse_output_tx.send(ParseOutput::Fastq(record)).await.is_err() {
                            eprintln!("Failed to send ParseOutput::Fastq");
                            break;
                        }
                    }
                });
                
            
            //     // filter_reads_out_stream = ReceiverStream::new(parse_output_rx);
                let test_write_task = tokio::spawn(stream_to_file(
                    parse_output_rx,
                    PathBuf::from("test_filter_reads_id.fq"),
                ));
            
                test_write_task.await??;

            }
    
    
    //         //*****************
    //         // Align Reads to Target
    //         
    //         // let align_ref_temp = NamedTempFile::new()?;
    //         // let align_ref_pipe_path = align_ref_temp.path().to_path_buf();
    //         // 
    //         // if align_ref_pipe_path.exists() {
    //         //     std::fs::remove_file(&align_ref_pipe_path)?;
    //         // }
    //         // Command::new("mkfifo")
    //         //     .arg(&align_ref_pipe_path)
    //         //     .status()?;
    //         // 
    //         // let filter_align_seq_clone = Arc::clone(&filter_align_seq);
    //         // let filter_align_accession_clone = Arc::clone(&filter_align_accession);
    //         // let align_ref_write_task = tokio::spawn({
    //         //     let align_ref_pipe_path = align_ref_pipe_path.clone();
    //         //     async move {
    //         //         write_hdf5_seq_to_fifo(&filter_align_seq_clone, &filter_align_accession_clone, &align_ref_pipe_path).await;
    //         //     }
    //         // });
    //         // 
    //         // let (align_query_write_task, align_query_pipe_path) = write_parse_output_to_temp(filter_reads_out_stream, None).await?;
    //         // 
    //         // let align_minimap2_args = generate_cli(MINIMAP2_TAG, &args, Some(&(align_ref_pipe_path.clone(), align_query_pipe_path.clone())))?;
    //         // let (mut align_minimap2_child,  align_minimap2_err_task) = spawn_cmd(MINIMAP2_TAG, align_minimap2_args, args.verbose).await?;
    //         // let align_minimap2_out_stream = parse_child_output(
    //         //     &mut align_minimap2_child,
    //         //     ChildStream::Stdout,
    //         //     ParseMode::Bytes,
    //         //     args.buffer_size / 4,
    //         // ).await?;
    //         // 
    //         // 
    //         // let align_samtools_config_sort = SamtoolsConfig {
    //         //     subcommand: SamtoolsSubcommand::Sort,
    //         //     // subcommand_fields: HashMap::from([("-o".to_string(), Some("test_samtools_align.sam".to_string())),("-".to_string(), None)]),
    //         //     subcommand_fields: HashMap::from([("-".to_string(), None)]),
    //         // };
    //         // let align_samtools_args_sort = generate_cli(
    //         //     SAMTOOLS_TAG,
    //         //     &args,
    //         //     Some(&align_samtools_config_sort),
    //         // )?;
    //         // let (mut align_samtools_child_sort, align_samtools_task_sort, align_samtools_err_task_sort) = stream_to_cmd(align_minimap2_out_stream, SAMTOOLS_TAG, align_samtools_args_sort, StreamDataType::JustBytes, args.verbose).await?;
    //         // let align_samtools_out_stream_sort = parse_child_output(
    //         //     &mut align_samtools_child_sort,
    //         //     ChildStream::Stdout,
    //         //     ParseMode::Bytes,
    //         //     args.buffer_size / 4,
    //         // ).await?;
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    //         // align_query_write_task.await??;
    //         // align_samtools_task_sort.await??;
    //         // align_samtools_err_task_sort.await??;
    //         // test_write_task.await??;
    // 
    //         //*****************
    //         // Make Consensus
    // 
    //         // let align_bam_path = file_path_manipulator(&no_ext_sample_base_buf, &cwd.clone(), None, Some("align.bam"), "_");
    //         // 
    //         // let mut align_samtools_out_stream_sort = ReceiverStream::new(align_samtools_out_stream_sort);
    //         // 
    //         // // Split stream for BAM writing
    //         // let (consensus_bam_streams, consensus_bam_done_rx) = t_junction(
    //         //     align_samtools_out_stream_sort,
    //         //     2,
    //         //     args.buffer_size,
    //         //     args.stall_threshold.try_into().unwrap(),
    //         //     Some(args.stream_sleep_ms),
    //         //     50,
    //         // )
    //         //     .await?;
    //         // 
    //         // if consensus_bam_streams.len() != 2 {
    //         //     return Err(anyhow!("Expected exactly 2 streams, got {}", consensus_bam_streams.len()));
    //         // }
    //         // 
    //         // let mut streams_iter = consensus_bam_streams.into_iter();
    //         // let consensus_bam_output_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing output stream"))?;
    //         // let consensus_bam_file_stream = streams_iter.next().ok_or_else(|| anyhow!("Missing file stream"))?;
    //         // 
    //         // let consensus_bam_write_task = tokio::spawn(stream_to_file(
    //         //     consensus_bam_file_stream,
    //         //     PathBuf::from(align_bam_path),
    //         // ));
    // 
    // 
    //         //Check target type, only allow viral
    //         match args.target_type {
    //             TargetType::Viral => {
    //                
    //                 
    // 
    //             }  // End if viral
    // 
    //             _ => {
    //                 return Err(anyhow!("Only viral consensus target types supported at this time."));
    //             }
    //         }
    // 
    //         //Samtools mpileup
    // 
    //         
    // 
    // 
    // 
    // 
    //         // consensus_bam_write_task.await??;
    //         // consensus_samtools_index_stream_task.await??;
    //         // consensus_index_err_task.await??;
    //         // align_ref_write_task.await?;
    //         // align_query_write_task.await??;
    //         // align_output_write_task.await??;
    //         ercc_query_write_task.await??;
    //         // ercc_output_write_task.await??;
    //         ercc_minimap2_err_task.await??;
    //         ercc_samtools_task_view.await??;
    //         ercc_samtools_task_stats.await??;
    //         ercc_stats_write_task.await??;
    //         ercc_done_rx.await??;
    //         
    // 
    //         
        } // end tech illumina
        Technology::ONT => {
            return Err(anyhow!("Minion not ready"));
        }
    
    }

    //*****************
    // Cleanup, hanging tasks







    // let minimap2_status = host_minimap2_child.wait().await?;
    // if minimap2_status.success() {
    //     eprintln!("Minimap2 exited successfully");
    // }
    // else {
    //     return Err(anyhow!("Minimap2 exited with non-zero status: {}", minimap2_status));
    // }
    //
    // host_samtools_task_view.await??;
    // host_samtools_err_task_view.await??;
    // let samtools_status = host_samtools_child_view.wait().await?;
    // if samtools_status.success() {
    //     eprintln!("Samtools exited successfully");
    // }
    // else {
    //     return Err(anyhow!("Samtools exited with non-zero status: {}", samtools_status));
    // }

    // if host_ref_pipe_path.exists() {
    //     std::fs::remove_file(&host_ref_pipe_path)?;
    // }






    
    // host_samtools_task_fastq.await??;
    // host_samtools_err_task_fastq.await??;

    // let host_samtools_task_fastq_status = host_samtools_child_fastq.wait().await?;
    // if  host_samtools_task_fastq_status.success() {
    //     eprintln!("Samtools fastq exited successfully");
    // }
    // else {
    //     return Err(anyhow!("Samtools fastq exited with non-zero status: {}", host_samtools_task_fastq_status));
    // }

    // host_samtools_write_task.await??;


    let results = try_join_all(cleanup_tasks).await?;
    for result in results {
        result?; // Propagate inner Result errors
    }

    for receiver in cleanup_receivers {
        receiver.await??; 
    }
    
    
    
    eprintln!("Finished generating consensus genome");
    
    Ok(())
}