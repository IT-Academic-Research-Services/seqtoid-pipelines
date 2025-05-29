use std::fs;
use anyhow::anyhow;
use std::path::PathBuf;
use std::time::Instant;
use crate::cli::Arguments;
use crate::config::defs::{FASTA_EXTS, H5DUMP_TAG};
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::command::check_version;
use crate::utils::file::{file_path_manipulator, scan_files_with_extensions};
use crate::utils::fastx::read_and_interleave_sequences;
use crate::utils::db::{write_sequences_to_hdf5, build_new_in_memory_index, check_db};



/// Creates a new, HDF5-backed database of id's and
/// sequences.
///
/// # Arguments
///
/// * `args` - Parsed CLI arguments.
///
/// # Returns
/// anyhow::Result<()>
///
pub async fn create_db(args: &Arguments) -> anyhow::Result<()> {
    println!("\n-------------\n Create DB\n-------------\n");
    println!("Generating HDF5 DB");
    let start = Instant::now();
    
    let h5v = check_version(H5DUMP_TAG).await?;
    eprintln!("HDF5 version: {}", h5v);

    let technology = Some(args.technology.clone());

    let cwd = std::env::current_dir()?;
    
    let mut all_fastas = Vec::new();


    match &args.file1 {
        Some(file) => {
            let fasta_full_path = file_path_manipulator(&PathBuf::from(file), &cwd, None, None, "");
            if fasta_full_path.exists() {
                all_fastas.push(fasta_full_path.clone());
            } else {
                return Err(anyhow!("File path does not exist: {}", fasta_full_path.display()));
            }
        }
        None => {

        }
    };
    // let ref_accession = match technology {
    if all_fastas.is_empty() {
        let the_fastas = match scan_files_with_extensions(&cwd, FASTA_EXTS)
        {
            Ok(the_fastas) => the_fastas,
            Err(e) => {
                return Err(anyhow!(e))
            }
        };
        

        all_fastas = the_fastas.clone();
    }
    eprintln!("FASTA files for input {:?}", all_fastas);
    
    let h5_path = match &args.ref_db {
        Some(file) => {
            let file_path = file_path_manipulator(&PathBuf::from(file), &cwd, None, None, "");
            
            file_path.with_extension("h5")
        }
        None => {
            file_path_manipulator(&PathBuf::from("db.h5"), &cwd, None, None, "")
        }
    };
    eprintln!("Writing to: {}", h5_path.display());
    let _ = fs::remove_file(&h5_path);
    
    for fasta_path in all_fastas {
        eprintln!("Writing from: {}", fasta_path.display());
        let rx = read_and_interleave_sequences(
            fasta_path,
            None,
            technology.clone(),
            args.max_reads,
            args.min_read_len,
            args.max_read_len,
        )?;
        let mut rx_stream = ReceiverStream::new(rx);
        write_sequences_to_hdf5(&mut rx_stream, &h5_path).await?;
    }
    
    let index_path = h5_path.with_extension("index.bin");
    
    let _index_map = build_new_in_memory_index(&h5_path, &index_path).await?;
    check_db(&h5_path, &index_path, None).await?;

    println!("Created DB File: {}", h5_path.display());
    println!("Created DB Index: {}", index_path.display());
    let elapsed = start.elapsed();
    let elapsed_secs = elapsed.as_secs_f64();
    println!("Time elapsed: {} seconds", elapsed_secs);
    Ok(())
}
