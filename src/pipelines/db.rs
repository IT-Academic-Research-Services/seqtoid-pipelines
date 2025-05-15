use std::fs;
use std::path::PathBuf;
use std::time::Instant;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::Arguments;
use crate::utils::fastx::{read_and_interleave_sequences};
use crate::utils::file::{file_path_manipulator, extension_remover};
use crate::utils::command::check_version;
use crate::utils::defs::H5DUMP_TAG;
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
    println!("HDF5 version: {}", h5v);

    let cwd = std::env::current_dir()?;
    let fasta_path = file_path_manipulator(&PathBuf::from(&args.file1), &cwd, None, None, "");
    eprintln!("{}", fasta_path.display());

    let (stem, _extensions) = extension_remover(&fasta_path);
    let base = stem.to_str().unwrap_or("");

    let technology = Some(args.technology.clone());
    eprintln!("Starting read_and_interleave_sequences");
    let rx = read_and_interleave_sequences(
        fasta_path,
        None,
        technology,
        args.max_reads,
        args.min_read_len,
        args.max_read_len,
    )?;
    let mut rx_stream = ReceiverStream::new(rx);

    let hdf5_file_name = format!("{}.h5", base);
    let _ = fs::remove_file(&hdf5_file_name);
    eprintln!("Writing to HDF5: {}", hdf5_file_name);
    write_sequences_to_hdf5(&mut rx_stream, &hdf5_file_name).await?;
    eprintln!("Finished writing HDF5");


    let index_file_name = format!("{}.index.bin", base);
    let _ = fs::remove_file(&index_file_name);
    eprintln!("building index map");
    let _index_map = build_new_in_memory_index(&hdf5_file_name, index_file_name.as_str()).await?;
    eprintln!("building index map complete");

    eprintln!("Checking DB");
    check_db(hdf5_file_name.as_str(), &index_file_name, None).await?;
    eprintln!("Checking DB complete");


    let elapsed = start.elapsed();
    let elapsed_secs = elapsed.as_secs_f64();
    println!("Created DB File: {} seconds", elapsed_secs);
    Ok(())
}
