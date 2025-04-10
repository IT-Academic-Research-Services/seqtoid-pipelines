mod pipelines;
mod utils;

use std::time::Instant;
use std::env;
use clap::Parser;
use utils::Arguments;

use pipelines::consensus_genome;

#[tokio::main]
async fn main() {
    let run_start = Instant::now();
    eprintln!("\n-------------\n SeqToID\n-------------\n");

    let dir = env::current_dir().unwrap();
    eprintln!("The current directory is {:?}\n", dir);

    let args = Arguments::parse();

    match args.module.as_str() {
        "consensus_genome" => consensus_genome_run(&args).await,
        _ => panic!("Invalid option!"),
    }

    eprintln!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
}

async fn consensus_genome_run(args: &Arguments) {
    consensus_genome::run(args).await.unwrap();
}