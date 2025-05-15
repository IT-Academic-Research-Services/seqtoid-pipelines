mod pipelines;
mod utils;

use std::time::Instant;
use std::env;
use clap::Parser;
use anyhow::Result;


use crate::cli::{parse, Arguments};
use pipelines::consensus_genome;
use pipelines::db;

mod cli;


#[tokio::main]
async fn main() -> Result<()> {
    let run_start = Instant::now();
    println!("\n-------------\n SeqToID\n-------------\n");

    let dir = env::current_dir()?;
    println!("The current directory is {:?}\n", dir);

    let args = parse();

    if let Err(e) = match args.module.as_str() {
        "consensus_genome" => consensus_genome_run(&args).await,
        "create_db" => create_db_run(&args).await,
        _ => Err(anyhow::anyhow!("Invalid module: {}", args.module)),
    } {
        eprintln!("Pipeline failed: {}", e);
        std::process::exit(1);
    }

    println!("Run complete: {} milliseconds.", run_start.elapsed().as_millis());
    Ok(())
}

async fn consensus_genome_run(args: &Arguments) -> Result<()> {
    consensus_genome::run(args).await
}
async fn create_db_run(args: &Arguments) -> Result<()> {
    db::create_db(args).await
}