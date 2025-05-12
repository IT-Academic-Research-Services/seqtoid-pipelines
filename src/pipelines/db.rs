use std::path::PathBuf;
use futures::StreamExt;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::Arguments;
use crate::utils::fastx::{read_and_interleave_sequences, SequenceRecord};
use crate::utils::file::file_path_manipulator;
use hdf5_metno::filters::{blosc_set_nthreads, Blosc, BloscShuffle};
use hdf5_metno::{File, Result, Extent};
use hdf5_metno::types::{VarLenArray, VarLenUnicode};
use tokio::task;

pub async fn create_db(args: &Arguments) -> anyhow::Result<()> {
    println!("\n-------------\n Create DB\n-------------\n");
    println!("Generating HDF5 DB");

    let cwd = std::env::current_dir()?;
    let fasta_path = file_path_manipulator(&PathBuf::from(&args.file1), &cwd, None, None, "");
    eprintln!("{}", fasta_path.display());

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

    let hdf5_file_name = args
        .out_file
        .as_ref()
        .cloned()
        .unwrap_or_else(|| "hdf5_out.h5".to_string());

    write_sequences_to_hdf5(&mut rx_stream, &hdf5_file_name, args.threads).await?;
    Ok(())
}

async fn write_sequences_to_hdf5(
    rx_stream: &mut ReceiverStream<SequenceRecord>,
    hdf5_file_name: &str,
    threads: usize,
) -> anyhow::Result<()> {
    let hdf5_file = File::create(hdf5_file_name)?;
    let hdf5_group = hdf5_file.create_group("db")?;

    blosc_set_nthreads(threads.min(255) as u8);
    let chunk_size = 1000;

    // Create extensible dataset for sequences (VarLenArray<u8>)
    let seq_dataset = hdf5_group
        .new_dataset::<VarLenArray<u8>>()
        .shape([Extent::resizable(0)])
        .chunk([chunk_size])
        .blosc(Blosc::BloscLZ, 9, BloscShuffle::Byte)
        .create("sequences")?;

    // Create extensible dataset for IDs (VarLenUnicode)
    let id_dataset = hdf5_group
        .new_dataset::<VarLenUnicode>()
        .shape([Extent::resizable(0)])
        .chunk([chunk_size])
        .blosc(Blosc::BloscLZ, 9, BloscShuffle::Byte)
        .create("id")?;

    let mut seq_buffer: Vec<VarLenArray<u8>> = Vec::new();
    let mut id_buffer: Vec<VarLenUnicode> = Vec::new();

    while let Some(record) = rx_stream.next().await {

        seq_buffer.push(record.seq().into()); // Convert &[u8] to VarLenArray<u8>
        id_buffer.push(unsafe { VarLenUnicode::from_str_unchecked(record.id()) });

        if seq_buffer.len() >= chunk_size {
            let count = seq_buffer.len();
            write_chunk_async(
                seq_dataset.clone(),
                id_dataset.clone(),
                seq_buffer,
                id_buffer,
                count,
            )
                .await?;
            seq_buffer = Vec::new();
            id_buffer = Vec::new();
        }
    }

    if !seq_buffer.is_empty() {
        let count = seq_buffer.len();
        write_chunk_async(seq_dataset, id_dataset, seq_buffer, id_buffer, count).await?;
    }

    Ok(())
}

async fn write_chunk_async(
    seq_dataset: hdf5_metno::Dataset,
    id_dataset: hdf5_metno::Dataset,
    seq_buffer: Vec<VarLenArray<u8>>,
    id_buffer: Vec<VarLenUnicode>,
    count: usize,
) -> anyhow::Result<()> {
    task::spawn_blocking(move || write_chunk(&seq_dataset, &id_dataset, &seq_buffer, &id_buffer, count))
        .await??;
    Ok(())
}

fn write_chunk(
    seq_dataset: &hdf5_metno::Dataset,
    id_dataset: &hdf5_metno::Dataset,
    seq_buffer: &[VarLenArray<u8>],
    id_buffer: &[VarLenUnicode],
    count: usize,
) -> Result<()> {
    let current_size = seq_dataset.shape()[0];
    seq_dataset.resize([current_size + count])?;
    id_dataset.resize([current_size + count])?;

    seq_dataset.write_slice(seq_buffer, current_size..current_size + count)?;
    id_dataset.write_slice(id_buffer, current_size..current_size + count)?;
    Ok(())
}