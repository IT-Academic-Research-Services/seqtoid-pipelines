use std::fs;
use std::path::PathBuf;
use std::time::Instant;
use tokio_stream::wrappers::ReceiverStream;
use crate::cli::Arguments;
use crate::utils::fastx::{read_and_interleave_sequences};
use crate::utils::file::{file_path_manipulator, extension_remover};
use crate::utils::command::check_version;
use crate::config::defs::H5DUMP_TAG;
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

    let cwd = std::env::current_dir()?;
    let fasta_path = file_path_manipulator(&PathBuf::from(&args.file1), &cwd, None, None, "");

    let (stem, _extensions) = extension_remover(&fasta_path);
    // let base = stem.to_str().unwrap_or("");

    let technology = Some(args.technology.clone());

    let rx = read_and_interleave_sequences(
        fasta_path,
        None,
        technology,
        args.max_reads,
        args.min_read_len,
        args.max_read_len,
    )?;
    let mut rx_stream = ReceiverStream::new(rx);

    let h5_path = stem.with_extension("h5");
    eprintln!("{}", h5_path.display());
    let _ = fs::remove_file(&h5_path);
    write_sequences_to_hdf5(&mut rx_stream, &h5_path).await?;


    let index_path = stem.with_extension("index.bin");
    let _ = fs::remove_file(&index_path);

    let _index_map = build_new_in_memory_index(&h5_path, &index_path).await?;
    // check_db(hdf5_file_name.as_str(), &index_file_name, None).await?;
    
    println!("Created DB File: {}", h5_path.display());
    println!("Created DB Index: {}", index_path.display());
    let elapsed = start.elapsed();
    let elapsed_secs = elapsed.as_secs_f64();
    println!("Time elapsed: {} seconds", elapsed_secs);
    Ok(())
}
