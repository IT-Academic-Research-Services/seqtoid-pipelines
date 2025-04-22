use std::io;
use std::path::PathBuf;
use anyhow::Result;
use tokio_stream::StreamExt;
use crate::utils::Arguments;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base};
use crate::utils::streams::{stream_to_cmd, tee};
use crate::utils::command::generate_cli;
use tokio_stream::wrappers::errors::BroadcastStreamRecvError;
use crate::FASTP_TAG;


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
    
    // let start = Instant::now();
    // let fq1count = record_counter(&args.file1);
    // let duration = start.elapsed();
    // println!("Num records: {:?} Time in ms: {}", fq1count,  duration.as_millis());

    // eprintln!("{}", args.quality);



    use tokio_stream::StreamExt;
    use tokio::time::Instant;

    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(sample_base), &cwd, None, Option::from("validated"), "_");
    let mut rx = read_and_interleave_sequences(file1_path, file2_path, technology, args.max_reads, args.min_read_len, args.max_read_len)?;
    let (mut rrx, file_handle) = tee(rx, validated_interleaved_file_path, true).await;

    let start = Instant::now();
    let mut out_counter = 0;
    while let Some(result) = rrx.next().await {
        match result {
            Ok(_) => {
                out_counter += 1;
                if out_counter % 1000 == 0 {
                    eprintln!("Consumed {} records", out_counter);
                }
            }
            Err(BroadcastStreamRecvError::Lagged(skipped)) => {
                eprintln!("Broadcast stream lagged, skipped {} items after {} records", skipped, out_counter);
            }
        }
    }
    eprintln!("out count {}, took {} ms", out_counter, start.elapsed().as_millis());

    // Await the file writer task
    file_handle.await.unwrap();
    eprintln!("File writing completed");

    // let fastp_cmd = generate_cli(FASTP_TAG, &args)?;

    // let rx_child = stream_to_cmd(rrx, fastp_cmd);
    // 
    // while let Some(result) = rrx.next().await {
    //     match result {
    //         Ok(record) => {
    //             eprintln!("id: {}", record.id());
    //             eprintln!("seq: {}", String::from_utf8_lossy(&record.seq()));
    //         }
    //         Err(BroadcastStreamRecvError::Lagged(skipped)) => {
    //             eprintln!("Broadcast stream lagged, skipped {} items", skipped);
    //         }
    //     }
    // }


    // let mut out_counter = 0;
    // 
    // while let Some(record) = rx.recv().await {
    // 
    //     println!("ID: {}", record.id());
    //     // println!("Seq: {:?}", record.seq());
    // 
    //     let seq_len = record.seq().len();
    //     // println!("Sequence length: {}", seq_len);
    // 
    //     let seq_str = String::from_utf8_lossy(&record.seq()).to_string();
    //     // println!("Seq: {}", seq_str);
    // 
    //     // println!("---"); // Separator between records
    //     out_counter += 1;
    // }
    // eprintln!("out coutn {}", out_counter);


    Ok(())

}
