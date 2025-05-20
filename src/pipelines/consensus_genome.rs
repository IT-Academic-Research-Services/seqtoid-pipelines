use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::Arguments;
use std::process::Command;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;
use tokio_stream::wrappers::ReceiverStream;
use tokio::process::Command as TokioCommand;
use crate::utils::command::{generate_cli, check_version};
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base};
use crate::utils::streams::{t_junction, stream_to_cmd, StreamDataType, parse_child_output, ChildStream, ParseMode, stream_to_file};
use crate::config::defs::{PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG};
use crate::cli::Technology;
use crate::utils::db::{lookup_sequence, load_index, build_new_in_memory_index};

pub async fn run(args: &Arguments) -> Result<()> {
    println!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let cwd = std::env::current_dir()?;

    
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

    let host_accession = args.host_accession.clone().ok_or_else(|| {
        anyhow!("Host accession must be given (-a).")
    })?;
    
    // let ref_accession = match technology {
    //     Technology::Illumina => {
    //         args.ref_accession.clone().ok_or_else(|| {
    //             anyhow::anyhow!("HDF5 database file must be given (-d).")
    //         })?
    //     }
    //     Technology::ONT => {
    //         String::new() // or return an error if ref_accession is required
    //     }
    // };




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

    // NB: the await calls below cause run to pause
    // Await write tasks concurrently
    // let (pigz_write_result, fastp_write_result) = join!(pigz_write_task, fastp_write_task);
    // pigz_write_result??;
    // fastp_write_result??;
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
    
    
    //*****************
    //Fetch Reference from accession
    
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


    let ref_temp = NamedTempFile::new()?;
    let ref_pipe_path = ref_temp.path().to_path_buf();
    let query_temp = NamedTempFile::new()?;
    let query_pipe_path = query_temp.path().to_path_buf();

    #[cfg(unix)]
    {
        Command::new("mkfifo")
            .arg(&ref_pipe_path)
            .status()?;
        Command::new("mkfifo")
            .arg(&query_pipe_path)
            .status()?;
    }
    #[cfg(not(unix))]
    return Err(anyhow!("Named pipes are not supported on non-Unix systems. Only Unix-like systems supported."));
    
    
    let seq = lookup_sequence(&ref_db_path, &h5_index, &host_accession).await?;
    let ref_seq = String::from_utf8(seq)?;
    eprintln!("{}", ref_seq);


    let ref_write_task = tokio::spawn({
        let ref_pipe_path = ref_pipe_path.clone();
        async move {
            let mut ref_file = File::create(&ref_pipe_path).await?;
            ref_file.write_all(format!(">{}\n{}\n", host_accession, ref_seq).as_bytes()).await?;
            ref_file.flush().await?;
            Ok::<(), anyhow::Error>(())
        }
    });

    // let query_write_task = tokio::spawn({
    //     let query_pipe_path = query_pipe_path.clone();
    //     async move {
    //         let mut query_file = File::create(&query_pipe_path).await?;
    //         while let Some(record) = fastp_out_stream.next().await {
    //             query_file.write_all(format!("@{}\n", record.id()).as_bytes()).await?;
    //             query_file.write_all(record.seq()).await?;
    //             query_file.write_all(b"\n+\n").await?;
    //             query_file.write_all(record.qual().unwrap_or(b"").as_bytes()).await?;
    //             query_file.write_all(b"\n").await?;
    //         }
    //         query_file.flush().await?;
    //         Ok::<(), anyhow::Error>(())
    //     }
    // });
    
    //*****************
    //Host Removal

    let minimap2_version = match check_version(MINIMAP2_TAG).await {
        Ok(version) => {
            eprintln!("{}", version);
            version
        }
        Err(err) => {
            return Err(anyhow!("Error checking minimap2 version: {}", err));
        }
    };
    eprintln!("Minimap2 version: {}", minimap2_version);
    
    
    //*****************
    eprintln!("Finished generating consensus genome");
    
    Ok(())
}