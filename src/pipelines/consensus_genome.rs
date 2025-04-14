
use crate::utils::Arguments;
use crate::utils::fastx::read_and_interleave_sequences;


pub async fn run(args: &Arguments) -> anyhow::Result<()> {
    eprintln!("\n-------------\n Consensus Genome\n-------------\n");
    println!("Running consensus genome with module: {}", args.module);
    println!("File 1: {}", args.file1);
    match &args.file2 {
        Some(file) => println!("File 2: {}", file),
        None => println!("File2 not given"),
    }
    let technology = Some(args.technology.clone());

    // let start = Instant::now();
    // let fq1count = record_counter(&args.file1);
    // let duration = start.elapsed();
    // println!("Num records: {:?} Time in ms: {}", fq1count,  duration.as_millis());
    
    let mut rx = read_and_interleave_sequences(&args.file1, args.file2.as_deref(), technology, args.max_reads)?;
    // 
    // 
    while let Some(record) = rx.recv().await {

        println!("ID: {}", record.id());

        println!("---"); // Separator between records
    }

    Ok(())

}
