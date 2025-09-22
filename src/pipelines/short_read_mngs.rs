use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use futures::future::try_join_all;
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use crate::config::defs::{PipelineError, RunConfig, StreamDataType, ReadStats, MINIMAP2_TAG, BOWTIE2_TAG, SAMTOOLS_TAG, FASTP_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, SEQKIT_TAG, QUAST_TAG, HISAT2_TAG, SamtoolsSubcommand};
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file};
use crate::utils::fastx::{raw_read_count, read_fastq};
use crate::utils::streams::{t_junction, ParseOutput, join_with_error_handling, stream_to_cmd, parse_child_output, ChildStream, ParseMode, stream_to_file, read_child_output_to_vec};
use crate::utils::command::bowtie2::{Bowtie2Config, bowtie2_index_prep};
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::command::samtools::SamtoolsConfig;


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


/// Pre-filters ERCC reads from the input stream
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns
///
async fn ercc_bowtie2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>, // FASTQ raw byte stream from
    bt2_index_path: PathBuf,
    paired: bool,

) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {

    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let ercc_bt2_config_view = Bowtie2Config{
        bt2_index_path: bt2_index_path.clone(),
        option_fields: HashMap::from([("--very-sensitive-local" .to_string(), None)]),
    };

    let ercc_bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&ercc_bt2_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut ercc_bt2_child, ercc_bt2_stream_task, ercc_bt2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        BOWTIE2_TAG,
        ercc_bt2_args,
        StreamDataType::JustBytes, // FASTQ bytes
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(ercc_bt2_stream_task);
    cleanup_tasks.push(ercc_bt2_err_task);

    let ercc_bt2_out_stream = {
        let mut guard = ercc_bt2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BOWTIE2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let ercc_samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-u".to_string(), None), // Uncompressed
            ("-O".to_string(), Some("bam".to_string())), // BAM output
            ("-".to_string(), None)
        ]),
    };
    let ercc_samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&ercc_samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut ercc_samtools_sort_child, ercc_samtools_sort_task, ercc_samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        ercc_bt2_out_stream,
        SAMTOOLS_TAG,
        ercc_samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(ercc_samtools_sort_task);
    cleanup_tasks.push(ercc_samtools_sort_err_task);

    let ercc_samtools_sort_out_stream = {
        let mut guard = ercc_samtools_sort_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };


    let (ercc_bam_streams, ercc_bam_done_rx) = t_junction(
        ReceiverStream::new(ercc_samtools_sort_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,  // Semantically FASTQ bytes
        "ercc_bam_split".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(ercc_bam_done_rx);


    let mut ercc_bam_streams_iter = ercc_bam_streams.into_iter();
    let ercc_bam_count_stream = ercc_bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let ercc_unmapped_stream = ercc_bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    //Count mapped reads to ERCC index
    let mapped_count_flag = if paired { "-F13".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None), // Count reads, don't pass anything else
            (mapped_count_flag, None),
            ("-".to_string(), None),
        ]),
    };
    let samtools_count_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_count_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut count_child_arc, count_stream_task, count_err_task) = stream_to_cmd(
        config.clone(),
        ercc_bam_count_stream,
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (ercc_count_tx, ercc_count_rx) = oneshot::channel::<u64>();

    let count_future = tokio::spawn(async move {
        let mut count_child = count_child_arc;
        let mut guard = count_child.lock().await;
        let count_lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        let _ = ercc_count_tx.send(mapped_count); // Send just u64
        Ok(())
    });
    cleanup_tasks.push(count_future);


    //Unmapped goes to samtools fastq for output stream,
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),  // Filter unmapped directly
            ("-".to_string(), None), // stdin/stdout
            // Add "-n" for name sorting if needed; "-1/-2" for paired FASTQ output
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        ercc_unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,  // Input: sorted BAM bytes
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let ercc_unmapped_fastq_stream= {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,  // Output: FASTQ records
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    Ok((ReceiverStream::new(ercc_unmapped_fastq_stream), ercc_count_rx , cleanup_tasks, cleanup_receivers))
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

    // External tools check
    check_versions(vec![
    BOWTIE2_TAG,
    HISAT2_TAG
    ])
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let (file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;

    let paired = file2_path.is_some();

    // Input Validation
    let (val_out_stream, validate_cleanup_tasks, validate_cleanup_receivers, raw_count_task, val_count_task) = validate_input(
        config.clone(),
        file1_path,
        file2_path,
        no_ext_sample_base_buf.clone(),
        &out_dir,
    )
        .await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);


    let ercc_bt2_index_path = bowtie2_index_prep(&config.args.ercc_bowtie2_index, &cwd)?;


    let (ercc_bt2_out_stream, ercc_count_rx, ercc_bt2_cleanup_tasks, ercc_bt2_cleanup_receivers) = ercc_bowtie2_filter(config.clone(), val_out_stream, ercc_bt2_index_path, paired).await?;
    cleanup_tasks.extend(ercc_bt2_cleanup_tasks);
    cleanup_receivers.extend(ercc_bt2_cleanup_receivers);


    let test_write_task = tokio::spawn(stream_to_file(
        ercc_bt2_out_stream.into_inner(),
        PathBuf::from("ercc-test.fq"),
    ));
    test_write_task.await;


    let ercc_mapped_count = ercc_count_rx.await.map_err(|e| PipelineError::Other(anyhow::anyhow!("ERCC count receiver failed: {}", e)))?;


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
