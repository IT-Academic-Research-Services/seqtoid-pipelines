use std::path::PathBuf;
use std::sync::Arc;
use futures::future::try_join_all;
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use crate::config::defs::{PipelineError, RunConfig, StreamDataType, ReadStats};
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file};
use crate::utils::fastx::{raw_read_count, read_fastq};
use crate::utils::streams::{t_junction, ParseOutput, join_with_error_handling};




/// Called read_fastq in single or paired FASTQ's and streams interleaved output
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `file1_path` - Path to R1 or single ended FASTQ.
/// * `file2_path` - Optional path to R2.
///  * `sample_base_buf` - PathBuf of sample basename.
/// * `out_dir` - Base dir for output files.
///
/// # Returns
/// ParseOutput validated interleaved FASTQ stream, vecs of cleaup takss and recoevers, handle for raw and validated read counts.
async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> anyhow::Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>>, JoinHandle<anyhow::Result<u64, anyhow::Error>>, JoinHandle<anyhow::Result<ReadStats, anyhow::Error>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();


    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated.fq.gz"),
        "_",
    );

    let raw_count_task = raw_read_count(file1_path.clone(), file2_path.clone());


    let (rx, val_count_task) = read_fastq(
        file1_path,
        file2_path,
        Some(config.args.technology.clone()),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        config.base_buffer_size,
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);

    // Split the byte stream for fastp and write
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,
        "validate_input".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(val_done_rx);

    let mut val_streams_iter = val_streams.into_iter();
    let val_ercc_bowtie2_filter_out_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let val_file_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;


    let val_file_write_task = write_byte_stream_to_file(
        &validated_interleaved_file_path,
        ReceiverStream::new(val_file_stream),
        Some(config.base_buffer_size),
    )
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    cleanup_tasks.push(val_file_write_task);

    Ok((ReceiverStream::new(val_ercc_bowtie2_filter_out_stream),  cleanup_tasks, cleanup_receivers, raw_count_task, val_count_task))
}



/// Run function for Short Read mNGS pipelines
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
///
/// # Returns
/// Result<(), PipelineError>
pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();
    let mut cleanup_tasks: Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>> = Vec::new();

    let (file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;


    // Input Validation
    let (val_ercc_bowtie2_filter_out_stream, validate_cleanup_tasks, validate_cleanup_receivers, raw_count_task, val_count_task) = validate_input(
        config.clone(),
        file1_path,
        file2_path,
        no_ext_sample_base_buf.clone(),
        &out_dir,
    )
        .await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);


    let raw_count = join_with_error_handling(raw_count_task).await?;
    println!("Processed {} raw reads (additive from R1 and R2 if paired)", raw_count);

    let stats = join_with_error_handling(val_count_task).await?;
    println!("Processed {} validated, {} undersized, {} oversized reads",
             stats.validated, stats.undersized, stats.oversized);

    // Cleanup
    let results = try_join_all(cleanup_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }
    for receiver in cleanup_receivers {
        receiver
            .await
            .map_err(|e| PipelineError::Other(e.into()))?
            .map_err(|e| PipelineError::Other(e))?;
    }

    println!("Finished short read mNGS.");
    Ok(())
}
