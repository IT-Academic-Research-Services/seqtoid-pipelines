use std::path::PathBuf;
use crate::utils::Arguments;
use crate::utils::file::file_path_manipulator;
use crate::utils::fastx::read_and_interleave_sequences;


pub async fn run(args: &Arguments) -> anyhow::Result<()> {
    eprintln!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);

    let file1_path = file_path_manipulator(PathBuf::from(&args.file1), None, None).unwrap();
    println!("File 1: {}", file1_path.display());

    let file2_path = match &args.file2 {
        Some(file) => {
            file_path_manipulator(PathBuf::from(file), None, None)
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
    
    let mut rx = read_and_interleave_sequences(file1_path, file2_path, technology, args.max_reads, args.min_read_len, args.max_read_len)?;
    // 
    // 
    while let Some(record) = rx.recv().await {
    
        println!("ID: {}", record.id());
        // println!("Seq: {:?}", record.seq());
    
        let seq_len = record.seq().len();
        println!("Sequence length: {}", seq_len);
    
        let seq_str = String::from_utf8_lossy(&record.seq()).to_string();
        println!("Seq: {}", seq_str);
    
        println!("---"); // Separator between records
    }

    Ok(())

}
