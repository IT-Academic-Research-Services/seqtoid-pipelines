mod pipelines;
mod utils;

use std::time::Instant;
use std::env;
use clap::Parser;
use anyhow::Result;
use utils::Arguments;

use pipelines::consensus_genome;



#[tokio::main]
async fn main() -> Result<()> {
    let run_start = Instant::now();
    eprintln!("\n-------------\n SeqToID\n-------------\n");

    let dir = env::current_dir()?;
    eprintln!("The current directory is {:?}\n", dir);

    let args = Arguments::parse();

    if let Err(e) = match args.module.as_str() {
        "consensus_genome" => consensus_genome_run(&args).await,
        _ => Err(anyhow::anyhow!("Invalid module: {}", args.module)),
    } {
        eprintln!("Pipeline failed: {}", e);
        std::process::exit(1);
    }

    eprintln!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
    Ok(())
}

async fn consensus_genome_run(args: &Arguments) -> Result<()> {
    consensus_genome::run(args).await
}