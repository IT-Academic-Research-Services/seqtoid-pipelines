use std::path::PathBuf;
use hdf5_metno::{Extent, File, H5Type};
use hdf5_metno::types::{FixedAscii, VarLenArray};
use tokio::task;
use tokio_stream::wrappers::ReceiverStream;
use crate::utils::fastx::SequenceRecord;
use fxhash::FxHashMap as HashMap;
use futures::StreamExt;
use crate::cli::Technology;

const CHUNK_SIZE: usize = 1000;
const TEST_FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/test_7_nt.fa");
const TEST_FASTA_TOOLONG_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/test_7_nt_long_id.fa");
const TEST_FASTA_ID: &str = "test7";
const TEST_FASTA_SEQ: &str = "ACGT";

#[derive(H5Type, Clone, PartialEq)]
#[repr(C)]
struct IndexEntry {
    id: FixedAscii<24>,
    index: u64,
}

/// Asynchronous fxn for writing to an hdf5 file.
///
/// # Arguments
///
/// * `rx_stream` - ReceiverStream<SequenceRecord> from read_and_interleave_sequences or similar.
/// * `hdf_file_name' - HDF5 file to write to.
/// * `threads' - Number of Blosc compression threads.
///
/// # Returns
/// anyhow::Result<()>
///
pub async fn write_sequences_to_hdf5(
    rx_stream: &mut ReceiverStream<SequenceRecord>,
    hdf5_file_name: &PathBuf,
) -> anyhow::Result<()> {
    let hdf5_file = File::create(hdf5_file_name)?;
    let hdf5_group = hdf5_file.create_group("db")?;

    let seq_dataset = hdf5_group
        .new_dataset::<VarLenArray<u8>>()
        .shape([Extent::resizable(0)])
        .chunk([CHUNK_SIZE])
        .shuffle().deflate(6)
        .create("sequences")?;

    let id_dataset = hdf5_group
        .new_dataset::<FixedAscii<24>>()
        .shape([Extent::resizable(0)])
        .chunk([CHUNK_SIZE])
        .shuffle().deflate(6)
        .create("id")?;

    let index_dataset = hdf5_group
        .new_dataset::<IndexEntry>()
        .shape([Extent::resizable(0)])
        .chunk([CHUNK_SIZE])
        .shuffle().deflate(6)
        .create("index")?;

    let mut seq_buffer: Vec<VarLenArray<u8>> = Vec::new();
    let mut id_buffer: Vec<FixedAscii<24>> = Vec::new();
    let mut index_buffer: Vec<IndexEntry> = Vec::new();
    let mut global_index: u64 = 0;

    while let Some(record) = rx_stream.next().await {
        let accession = record
            .id()
            .split_whitespace()
            .next()
            .ok_or_else(|| anyhow::anyhow!("Invalid header: {}", record.id()))?;
        if accession.len() > 23 {
            eprintln!(
                "ERROR: Skipping ID '{}' ({} bytes), exceeds 23-byte limit",
                accession,
                accession.len()
            );
            continue;
        }
        let mut accession_bytes = accession.as_bytes().to_vec();
        accession_bytes.push(0); // Add null terminator
        let id = FixedAscii::from_ascii(&accession_bytes)
            .map_err(|e| anyhow::anyhow!("Invalid ASCII in '{}': {}", accession, e))?;
        seq_buffer.push(record.seq().into());
        id_buffer.push(id.clone());
        index_buffer.push(IndexEntry { id, index: global_index });
        global_index += 1;

        if seq_buffer.len() >= CHUNK_SIZE {
            let count = seq_buffer.len();
            index_buffer.sort_by(|a, b| a.id.as_str().cmp(b.id.as_str()));
            write_chunk_async(
                seq_dataset.clone(),
                id_dataset.clone(),
                index_dataset.clone(),
                seq_buffer,
                id_buffer,
                index_buffer,
                count,
            )
                .await?;
            seq_buffer = Vec::new();
            id_buffer = Vec::new();
            index_buffer = Vec::new();
        }
    }

    if !seq_buffer.is_empty() {
        let count = seq_buffer.len();
        index_buffer.sort_by(|a, b| a.id.as_str().cmp(b.id.as_str()));
        write_chunk_async(
            seq_dataset,
            id_dataset,
            index_dataset,
            seq_buffer,
            id_buffer,
            index_buffer,
            count,
        )
            .await?;
    }

    Ok(())
}

/// Asynchronous caller of write_chunk.
///
/// # Arguments
///
/// * `seq_dataset` - Dataset of sequence data.
/// * `id_dataset' - Dataset of id data.
/// * `index_dataset' - "dataset of id:sequence associations.
/// * `seq_buffer` - Buffer of data to be written to the sequence dataset.
/// * `id_buffer` - Buffer of data to be written to the id dataset.
/// * `index_buffer` - Buffer of data to be written to the index dataset.
/// * 'count' - Passed to use for resizing.
///
/// # Returns
/// anyhow::Result<()>
///
async fn write_chunk_async(
    seq_dataset: hdf5_metno::Dataset,
    id_dataset: hdf5_metno::Dataset,
    index_dataset: hdf5_metno::Dataset,
    seq_buffer: Vec<VarLenArray<u8>>,
    id_buffer: Vec<FixedAscii<24>>,
    index_buffer: Vec<IndexEntry>,
    count: usize,
) -> anyhow::Result<()> {
    task::spawn_blocking(move || {
        write_chunk(
            &seq_dataset,
            &id_dataset,
            &index_dataset,
            &seq_buffer,
            &id_buffer,
            &index_buffer,
            count,
        )
    })
        .await??;
    Ok(())
}

/// Writes new data to the HDF5 datasets.
///
/// # Arguments
///
/// * `seq_dataset` - Dataset of sequence data.
/// * `id_dataset' - Dataset of id data.
/// * `index_dataset' - "dataset of id:sequence associations.
/// * `seq_buffer` - Buffer of data to be written to the sequence dataset.
/// * `id_buffer` - Buffer of data to be written to the id dataset.
/// * `index_buffer` - Buffer of data to be written to the index dataset.
/// * 'count' - Passed to use for resizing.
///
/// # Returns
/// anyhow::Result<()>
///
fn write_chunk(
    seq_dataset: &hdf5_metno::Dataset,
    id_dataset: &hdf5_metno::Dataset,
    index_dataset: &hdf5_metno::Dataset,
    seq_buffer: &[VarLenArray<u8>],
    id_buffer: &[FixedAscii<24>],
    index_buffer: &[IndexEntry],
    count: usize,
) -> hdf5_metno::Result<()> {
    let current_size = seq_dataset.shape()[0];
    seq_dataset.resize([current_size + count])?;
    id_dataset.resize([current_size + count])?;
    index_dataset.resize([current_size + count])?;

    seq_dataset.write_slice(seq_buffer, current_size..current_size + count)?;
    id_dataset.write_slice(id_buffer, current_size..current_size + count)?;
    index_dataset.write_slice(index_buffer, current_size..current_size + count)?;
    Ok(())
}

/// Determines if a file path is a FASTA, FASTQ, or neither.
///
/// # Arguments
///
/// * `path` - Header line of a FASTX record.
///
/// # Returns
/// Result<String>. Ok fastq or fasta, or err.
///
pub async fn build_new_in_memory_index(h5_path: &PathBuf, cache_file_name: &PathBuf) -> anyhow::Result<HashMap<[u8; 24], u64>> {
    eprintln!("Building in-memory index for: {}", h5_path.display());
    let file = File::open(h5_path)?;
    let group = file.group("db")?;
    let index_dataset = group.dataset("index")?;

    let index_len = index_dataset.shape()[0];

    let mut index_map: HashMap<[u8; 24], u64> = HashMap::default();
    index_map.reserve(index_len);

    for start in (0..index_len).step_by(CHUNK_SIZE) {
        let end = (start + CHUNK_SIZE).min(index_len);
        let entries: Vec<IndexEntry> = index_dataset.read_slice(start..end)?.to_vec();
        for entry in entries {
            let id_bytes = entry.id.as_str().as_bytes();
            let mut id_bytes_array = [0u8; 24];
            id_bytes_array[..id_bytes.len()].copy_from_slice(&id_bytes[..id_bytes.len()]);
            index_map.insert(id_bytes_array, entry.index);
        }
    }

    let config = bincode::config::standard();
    let (serialized_data, index_map) = tokio::task::spawn_blocking({
        move || {
            let serialized = bincode::serde::encode_to_vec(&index_map, config)?;
            Ok::<_, bincode::error::EncodeError>((serialized, index_map))
        }
    }).await??;
    tokio::fs::write(&cache_file_name, serialized_data).await?;

    eprintln!("In-memory index built: {} entries", index_map.len());
    Ok(index_map)
}


/// Loads an index file associated with an HDF5 file.
///
/// # Arguments
///
/// * `index_file_name` - Name of index file.
///
/// # Returns
/// anyhow::Result<HashMap<[u8; 24], u64>> the index
///
pub async fn load_index(index_path: &PathBuf) -> anyhow::Result<HashMap<[u8; 24], u64>> {
    let config = bincode::config::standard();
    let data = tokio::fs::read(index_path).await?;
    let (index_map, _): (HashMap<[u8; 24], u64>, usize) = bincode::serde::decode_from_slice(&data, config)?;
    eprintln!("Loaded index from {}: {} entries", index_path.display(), index_map.len());
    Ok(index_map)
}

/// Database integrity check function
///
/// # Arguments
///
/// * `h5_file_name` - Name of index file.
/// * `index_file_name` - Name of index file.
/// * 'target_id` - Accession ID to test.
///
/// # Returns
/// anyhow::Result<HashMap<[u8; 24], u64>> the index
///
pub async fn check_db(h5_path: &PathBuf, index_path: &PathBuf, target_id: Option<&str>) -> anyhow::Result<()> {
    println!("Checking HDF5 file: {}", h5_path.display());
    let file = File::open(h5_path)?;
    let group = file.group("db")?;
    let id_dataset = group.dataset("id")?;
    let seq_dataset = group.dataset("sequences")?;
    let index_dataset = group.dataset("index")?;

    let id_len = id_dataset.shape()[0];
    let seq_len = seq_dataset.shape()[0];
    let index_len = index_dataset.shape()[0];
    if id_len != seq_len || id_len != index_len {
        return Err(anyhow::anyhow!(
            "Mismatched dataset lengths: id={} seq={} index={}",
            id_len,
            seq_len,
            index_len
        ));
    }
    println!("Dataset sizes: {} records", id_len);

    if let Some(id) = target_id {
        if id.len() > 23 {
            return Err(anyhow::anyhow!(
                "Target ID '{}' ({} bytes) exceeds 23-byte limit",
                id,
                id.len()
            ));
        }

        let index_map = load_index(index_path).await?;
        let seq = lookup_sequence(h5_path, &index_map, &id.to_string()).await?;
        println!("Found sequence for ID '{}': {:?}", id, seq);
    }

    Ok(())
}


/// Retrieve a sequence from an H5 file.
///
/// # Arguments
///
/// * `h5_file_name` - Name of index file.
/// * 'target_id` - Accession ID to test.
/// * `index_map' - Index associated with H5 file.
///
/// # Returns
/// anyhow::Result<HashMap<[u8; 24], u64>> the index
///
pub async fn lookup_sequence(h5_path: &PathBuf,  index_map: &HashMap<[u8; 24], u64>, target_id: &String) -> anyhow::Result<Vec<u8>> {
    if target_id.len() > 23 {
        return Err(anyhow::anyhow!(
            "Target ID '{}' ({} bytes) exceeds 23-byte limit",
            target_id,
            target_id.len()
        ));
    }
    let file = File::open(h5_path)?;
    let group = file.group("db")?;
    let seq_dataset = group.dataset("sequences")?;

    let mut id_bytes = [0u8; 24];
    if target_id.len() > 23 {
        return Err(anyhow::anyhow!(
        "Target ID '{}' ({} bytes) exceeds 23-byte limit",
        target_id,
        target_id.len()
    ));
    }
    id_bytes[..target_id.len()].copy_from_slice(target_id.as_bytes());
    let index = index_map
        .get(&id_bytes)
        .ok_or_else(|| anyhow::anyhow!("ID '{}' not found in dataset", target_id))?;
    
    if *index as usize >= seq_dataset.shape()[0] {
        return Err(anyhow::anyhow!(
            "Index {} out of bounds for dataset size {}",
            *index,
            seq_dataset.shape()[0]
        ));
    }

    let seq: VarLenArray<u8> = seq_dataset.read_slice::<VarLenArray<u8>, std::ops::Range<usize>, ndarray::Ix1>(*index as usize..*index as usize + 1)?[0].clone();
    Ok(seq.to_vec())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use super::*;
    use anyhow::anyhow;
    use tempfile::NamedTempFile;
    use crate::utils::fastx::{fastx_generator, read_and_interleave_sequences};
    use crate::utils::file::extension_remover;
    use crate::utils::streams::t_junction;

    #[tokio::test]
    async fn test_create_db_illumina_small() -> anyhow::Result<()> {
        let temp_file = NamedTempFile::new().unwrap();
        let hdf5_path = PathBuf::from(temp_file.path());
        let stream = fastx_generator(100, 150, 35.0, 3.0);

        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            100_000,
            100,
            Some(0),
            50
        ).await?;

        let rx = outputs.pop().ok_or_else(|| anyhow!("No output stream"))?;
        let mut rx_stream = ReceiverStream::new(rx);
        let result = write_sequences_to_hdf5(&mut rx_stream, &hdf5_path).await;
        assert!(result.is_ok());

        let file = File::open(hdf5_path).unwrap();
        let group = file.group("db").unwrap();
        let seq_dataset = group.dataset("sequences").unwrap();
        let id_dataset = group.dataset("id").unwrap();
        let index_dataset = group.dataset("index").unwrap();

        assert_eq!(seq_dataset.shape()[0], 100);
        assert_eq!(id_dataset.shape()[0], 100);
        assert_eq!(index_dataset.shape()[0], 100);

        let ids: Vec<FixedAscii<24>> = id_dataset.read_slice(0..100).unwrap().to_vec();
        assert_eq!(ids[0].as_str(), "read1");

        Ok(())
    }

    #[tokio::test]
    async fn test_create_db_ont_small() -> anyhow::Result<()> {
        let temp_file = NamedTempFile::new().unwrap();
        let hdf5_path = PathBuf::from(temp_file.path());
        let stream = fastx_generator(100, 5000, 35.0, 3.0);

        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            100_000,
            100,
            Some(0),
            50
        ).await?;

        let rx = outputs.pop().ok_or_else(|| anyhow!("No output stream"))?;
        let mut rx_stream = ReceiverStream::new(rx);
        let result = write_sequences_to_hdf5(&mut rx_stream, &hdf5_path).await;
        assert!(result.is_ok());

        let file = File::open(hdf5_path).unwrap();
        let group = file.group("db").unwrap();
        let seq_dataset = group.dataset("sequences").unwrap();
        let id_dataset = group.dataset("id").unwrap();
        let index_dataset = group.dataset("index").unwrap();

        assert_eq!(seq_dataset.shape()[0], 100);
        assert_eq!(id_dataset.shape()[0], 100);
        assert_eq!(index_dataset.shape()[0], 100);

        let ids: Vec<FixedAscii<24>> = id_dataset.read_slice(0..100).unwrap().to_vec();
        assert_eq!(ids[0].as_str(), "read1");

        Ok(())
    }

    #[tokio::test]
    async fn test_read_known_file_to_db() -> anyhow::Result<()> {
        let path = PathBuf::from(TEST_FASTA_PATH);
        assert!(path.exists(), "Test file {} does not exist", TEST_FASTA_PATH);

        let temp_file = NamedTempFile::new().unwrap();
        let hdf5_path = PathBuf::from(temp_file.path());

        let rx = read_and_interleave_sequences(
            PathBuf::from(TEST_FASTA_PATH),
            None,
            Some(Technology::Illumina),
            500000000,
            None,
            None,
        )?;

        let mut rx_stream = ReceiverStream::new(rx);
        let result = write_sequences_to_hdf5(&mut rx_stream, &hdf5_path).await;
        assert!(result.is_ok());

        let file = File::open(hdf5_path).unwrap();
        let group = file.group("db").unwrap();
        let seq_dataset = group.dataset("sequences").unwrap();
        let id_dataset = group.dataset("id").unwrap();
        let index_dataset = group.dataset("index").unwrap();

        assert_eq!(seq_dataset.shape()[0], 7);
        assert_eq!(id_dataset.shape()[0], 7);
        assert_eq!(index_dataset.shape()[0], 7);

        let ids: Vec<FixedAscii<24>> = id_dataset.read_slice(0..7).unwrap().to_vec();
        let seqs: Vec<VarLenArray<u8>> = seq_dataset.read_slice(0..7).unwrap().to_vec();
        assert_eq!(ids[0].as_str(), "NC_049488.1");
        assert_eq!(ids[6].as_str(), "test7");
        assert_eq!(seqs[6].to_vec(), b"ACGT");

        Ok(())
    }
    
    #[tokio::test]
    async fn test_seq_lookup() -> anyhow::Result<()> {

        let path = PathBuf::from(TEST_FASTA_PATH);
        assert!(path.exists(), "Test file {} does not exist", TEST_FASTA_PATH);

        let temp_file = NamedTempFile::new().unwrap();
        let hdf5_path = PathBuf::from(temp_file.path());

        let rx = read_and_interleave_sequences(
            PathBuf::from(TEST_FASTA_PATH),
            None,
            Some(Technology::Illumina),
            500000000,
            None,
            None,
        )?;

        let mut rx_stream = ReceiverStream::new(rx);
        let result = write_sequences_to_hdf5(&mut rx_stream, &hdf5_path).await;
        assert!(result.is_ok());

        let file = File::open(&hdf5_path).unwrap();
        let group = file.group("db").unwrap();
        let seq_dataset = group.dataset("sequences").unwrap();
        let id_dataset = group.dataset("id").unwrap();
        let index_dataset = group.dataset("index").unwrap();
        let ids: Vec<FixedAscii<24>> = id_dataset.read_slice(0..7).unwrap().to_vec();
        let seqs: Vec<VarLenArray<u8>> = seq_dataset.read_slice(0..7).unwrap().to_vec();
        assert_eq!(ids[6].as_str(), "test7");
        
        let index_path = hdf5_path.with_extension("index.bin");
        let _ = fs::remove_file(&index_path);
        
        let index_map = build_new_in_memory_index(&hdf5_path, &index_path).await?;
        
        let seq = lookup_sequence(&hdf5_path,  &index_map, &TEST_FASTA_ID.to_string()).await?;
        let seq_str = std::str::from_utf8(&seq).unwrap();
        assert_eq!(seq_str, TEST_FASTA_SEQ);
        
        Ok(())
        
    }

    #[tokio::test]
    async fn test_load_index() -> anyhow::Result<()> {
        let temp_file = NamedTempFile::new().unwrap();
        let hdf5_path = PathBuf::from(temp_file.path());
        let stream = fastx_generator(100, 150, 35.0, 3.0);

        let (mut outputs, done_rx) = t_junction(
            stream,
            1,
            100_000,
            100,
            Some(0),
            50
        ).await?;

        let rx = outputs.pop().ok_or_else(|| anyhow!("No output stream"))?;
        let mut rx_stream = ReceiverStream::new(rx);
        let result = write_sequences_to_hdf5(&mut rx_stream, &hdf5_path).await;
        assert!(result.is_ok());

        let (stem, _extensions) = extension_remover(&PathBuf::from(&hdf5_path));
        let base = stem.to_str().unwrap_or("");

        let index_path = hdf5_path.with_extension("index.bin");
        let _ = fs::remove_file(&index_path);

        let index_map_from_build = build_new_in_memory_index(&hdf5_path, &index_path).await?;

        let index_map_from_file = load_index(&index_path).await?;
        
        assert_eq!(index_map_from_build.len(), index_map_from_file.len());

        for (key, value) in index_map_from_build.iter() {
            match index_map_from_file.get(key) {
                Some(&map2_value) => {
                    assert_eq!(map2_value, *value);

                }
                None => assert!(false), // Key missing
            }
        }
        
        Ok(())
    }

    #[tokio::test]
    async fn test_read_too_long_to_db() -> anyhow::Result<()> {
        
        let path = PathBuf::from(TEST_FASTA_TOOLONG_PATH);
        
        assert!(path.exists(), "Test file {} does not exist", TEST_FASTA_TOOLONG_PATH);

        let temp_file = NamedTempFile::new().unwrap();
        let hdf5_path = PathBuf::from(temp_file.path());

        let rx = read_and_interleave_sequences(
            PathBuf::from(TEST_FASTA_TOOLONG_PATH),
            None,
            Some(Technology::Illumina),
            500000000,
            None,
            None,
        )?;

        let mut rx_stream = ReceiverStream::new(rx);
        let result = write_sequences_to_hdf5(&mut rx_stream, &hdf5_path).await;
        assert!(result.is_ok());

        let file = File::open(&hdf5_path).unwrap();
        let group = file.group("db").unwrap();
        let seq_dataset = group.dataset("sequences").unwrap();
        let id_dataset = group.dataset("id").unwrap();
        let index_dataset = group.dataset("index").unwrap();

        assert_eq!(seq_dataset.shape()[0], 6);
        assert_eq!(id_dataset.shape()[0], 6);
        assert_eq!(index_dataset.shape()[0], 6);

        let ids: Vec<FixedAscii<24>> = id_dataset.read_slice(0..6).unwrap().to_vec();
        let seqs: Vec<VarLenArray<u8>> = seq_dataset.read_slice(0..6).unwrap().to_vec();

        assert_eq!(ids[0].as_str(), "NC_039199.1"); // NB: This is actually the second entry, 
        // so this assertion implies the first entry was skipped for excessive length
        assert_eq!(ids[5].as_str(), "test7");
        assert_eq!(seqs[5].to_vec(), b"ACGT");

        Ok(())
    }

}