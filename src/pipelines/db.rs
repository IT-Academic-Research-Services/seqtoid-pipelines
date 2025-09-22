use std::fs;
use anyhow::anyhow;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Instant;
use futures::future::try_join_all;
use crate::config::defs::{PipelineError, RunConfig};
use crate::config::defs::{FASTA_EXTS, H5DUMP_TAG};
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;
use crate::cli::Technology;
use crate::utils::command::version_check;
use crate::utils::file::{file_path_manipulator, scan_files_with_extensions, extension_remover};
use crate::utils::fastx::read_fastq;
use crate::utils::db::{write_sequences_to_hdf5, build_new_in_memory_index, check_db};
use crate::utils::streams::{ChildStream, ParseOutput};

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
pub async fn create_db(config: &RunConfig) -> anyhow::Result<()> {
    println!("\n-------------\n Create DB\n-------------\n");
    println!("Generating HDF5 DB");
    let start = Instant::now();
    let mut cleanup_tasks = Vec::new();

    let h5v = version_check(H5DUMP_TAG, vec!["-V"], 0, 2, ChildStream::Stdout).await?;
    eprintln!("HDF5 version: {}", h5v);

    let technology = Some(config.args.technology.clone());

    let cwd = std::env::current_dir()?;

    let mut all_singles: Vec<PathBuf> = Vec::new();
    let mut all_multis: Vec<PathBuf> = Vec::new();

    match &config.args.file1 {
        Some(f1) => {
            let eff_1 = PathBuf::from_str(&f1)?;

            if eff_1.is_dir() {
                let the_multis = match scan_files_with_extensions(&eff_1, FASTA_EXTS) {
                    Ok(the_multis) => the_multis,
                    Err(e) => {
                        return Err(anyhow!(e))
                    }
                };
                all_multis = the_multis.clone();
            } else {
                let fasta_full_path = file_path_manipulator(&PathBuf::from(f1), Some(&cwd), None, None, "");
                if fasta_full_path.exists() {
                    all_multis.push(fasta_full_path.clone());
                } else {
                    return Err(anyhow!("File path does not exist: {}", fasta_full_path.display()));
                }
            }
        }
        None => {}
    };

    match &config.args.file2 {
        Some(f2) => {
            let eff_2 = PathBuf::from_str(&f2)?;
            if eff_2.is_dir() {
                let the_singles = match scan_files_with_extensions(&eff_2, FASTA_EXTS) {
                    Ok(the_singles) => the_singles,
                    Err(e) => {
                        return Err(anyhow!(e))
                    }
                };
                all_singles = the_singles.clone();
            } else {
                let fasta_full_path = file_path_manipulator(&PathBuf::from(f2), Some(&cwd), None, None, "");
                if fasta_full_path.exists() {
                    all_singles.push(fasta_full_path.clone());
                } else {
                    return Err(anyhow!("File path does not exist: {}", fasta_full_path.display()));
                }
            }
        }
        None => {}
    };

    if all_multis.is_empty() && all_singles.is_empty() {
        return Err(anyhow!("No input files found!"))
    }

    eprintln!("Multi-sample FASTA files for input {:?}", all_multis);
    eprintln!("Single-sample FASTA files for input {:?}", all_singles);

    let h5_path = match &config.args.ref_db {
        Some(file) => {
            let file_path = file_path_manipulator(&PathBuf::from(file), Some(&cwd), None, None, "");
            file_path.with_extension("h5")
        }
        None => {
            file_path_manipulator(&PathBuf::from("db.h5"), Some(&cwd), None, None, "")
        }
    };
    eprintln!("Writing to: {}", h5_path.display());
    let _ = fs::remove_file(&h5_path);

    for multi_path in all_multis {
        eprintln!("Writing from: {}", multi_path.display());

        let (rx, count_task) = read_fastq(
            multi_path,
            None,
            technology.clone(),
            config.args.max_reads as u64,
            config.args.min_read_len,
            config.args.max_read_len,
            100_000,
        )?;

        let multi_count_cleanup_task = tokio::spawn(async move {
            if let Err(e) = count_task.await {
                eprintln!("create_db: Warning - read_fastq task failed for multi-sample: {}", e);
            }
            Ok::<(), anyhow::Error>(())
        });
        cleanup_tasks.push(multi_count_cleanup_task);
        let mut rx_stream = ReceiverStream::new(rx);
        write_sequences_to_hdf5(&mut rx_stream, &h5_path, None).await?;
    }

    for single_path in all_singles {
        eprintln!("Writing from: {}", single_path.display());

        let sample_base_buf: PathBuf = PathBuf::from(&single_path);
        let (no_ext_sample_base_buf, _) = extension_remover(&sample_base_buf);
        let no_ext_sample_base = no_ext_sample_base_buf.to_string_lossy().into_owned();
        eprintln!("ID in h5 will be: {}", no_ext_sample_base);

        let (rx, count_task) = read_fastq(
            single_path.clone(),
            None,
            technology.clone(),
            config.args.max_reads as u64,
            config.args.min_read_len,
            config.args.max_read_len,
            100_000
        )?;

        let single_count_cleanup_task = tokio::spawn(async move {
            if let Err(e) = count_task.await {
                eprintln!("create_db: Warning - read_fastq task failed for single-sample {}: {}", single_path.display(), e);
            }
            Ok::<(), anyhow::Error>(())
        });
        cleanup_tasks.push(single_count_cleanup_task);

        let mut rx_stream = ReceiverStream::new(rx);
        write_sequences_to_hdf5(&mut rx_stream, &h5_path, Some(no_ext_sample_base)).await?;
    }

    let index_path = h5_path.with_extension("index.bin");

    let _index_map = build_new_in_memory_index(&h5_path, &index_path).await?;
    check_db(&h5_path, &index_path, Option::from(None)).await?;

    println!("Created DB File: {}", h5_path.display());
    println!("Created DB Index: {}", index_path.display());

    // Cleanup
    for task in cleanup_tasks {
        if let Err(e) = task.await {
            eprintln!("create_db: Cleanup task failed: {}", e);
        }
    }

    let elapsed = start.elapsed();
    let elapsed_secs = elapsed.as_secs_f64();
    println!("Time elapsed: {} seconds", elapsed_secs);
    Ok(())
}
