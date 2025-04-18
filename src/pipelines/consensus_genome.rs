use std::path::PathBuf;
use tokio_stream::StreamExt;
use crate::utils::Arguments;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::{read_and_interleave_sequences, r1r2_base};
use crate::utils::streams::tee;
use tokio_stream::wrappers::errors::BroadcastStreamRecvError;


pub async fn run(args: &Arguments) -> anyhow::Result<()> {
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

    let validated_interleaved_file_path = file_path_manipulator(&PathBuf::from(sample_base), &cwd, None, Option::from("validated"), "_");

    


    let rx = read_and_interleave_sequences(file1_path, file2_path, technology, args.max_reads, args.min_read_len, args.max_read_len)?;

    let mut rrx = tee(rx, validated_interleaved_file_path, true).await;
    while let Some(result) = rrx.next().await {
        match result {
            Ok(record) => {
                eprintln!("id: {}", record.id());
                eprintln!("seq: {}", String::from_utf8_lossy(&record.seq()));
            }
            Err(BroadcastStreamRecvError::Lagged(skipped)) => {
                eprintln!("Broadcast stream lagged, skipped {} items", skipped);
            }
        }
    }
    
    

    // // 
    // while let Some(record) = rx.recv().await {
    // 
    //     println!("ID: {}", record.id());
    //     // println!("Seq: {:?}", record.seq());
    // 
    //     let seq_len = record.seq().len();
    //     println!("Sequence length: {}", seq_len);
    // 
    //     let seq_str = String::from_utf8_lossy(&record.seq()).to_string();
    //     println!("Seq: {}", seq_str);
    // 
    //     println!("---"); // Separator between records
    // }

    Ok(())

}
