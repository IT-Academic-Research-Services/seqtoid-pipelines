use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::fs::File;
use std::io::Write;
use regex::Regex;
use anyhow::{anyhow, Result};
use futures::future::try_join_all;
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;
use crate::config::defs::{PipelineError, RunConfig, StreamDataType, ReadStats, MINIMAP2_TAG, BOWTIE2_TAG, SAMTOOLS_TAG, FASTP_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, SEQKIT_TAG, QUAST_TAG, HISAT2_TAG, SamtoolsSubcommand, KALLISTO_TAG, KallistoSubcommand};
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file};
use crate::utils::fastx::{raw_read_count, read_fastq, stream_record_counter};
use crate::utils::streams::{t_junction, ParseOutput, join_with_error_handling, stream_to_cmd, parse_child_output, ChildStream, ParseMode, stream_to_file, read_child_output_to_vec};
use crate::utils::command::bowtie2::{Bowtie2Config, bowtie2_index_prep};
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::kallisto::KallistoConfig;

#[derive(Debug)]
pub struct KallistoResults {
    pub ercc_counts: Vec<(String, f64)>, // (target_id, est_counts)
    pub transcript_to_gene: Vec<(String, String)>, // (transcript_id, gene_id)
}

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

    Ok((ReceiverStream::new(ercc_unmapped_fastq_stream), ercc_count_rx, cleanup_tasks, cleanup_receivers))
}



/// QC's input stream using FASTP
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns
///
async fn fastp_qc(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
) -> Result<(ReceiverStream<ParseOutput>,  Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, oneshot::Receiver<Result<u64, anyhow::Error>>), PipelineError>{
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let qc_fastp_config_view = FastpConfig {

        //These default QC thresholds are loosely based on the pre-2022 CZI pipeline using PriceSeq & LZW
        command_fields: HashMap::from([
            ("--dont_eval_duplication".to_string(), None),
            ("--length_required".to_string(), Some("35".to_string())),
            ("--qualified_quality_phred".to_string(), Some("17".to_string())),
            ("--unqualified_percent_limit".to_string(), Some("15".to_string())),
            ("--n_base_limit".to_string(), Some("15".to_string())),
            ("--complexity_threshold".to_string(), Some("60".to_string())),
        ]),
    };

    let qc_fastp_args = generate_cli(FASTP_TAG, &config, Some(&qc_fastp_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut qc_fastp_child, qc_fastp_stream_task, qc_fastp_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        FASTP_TAG,
        qc_fastp_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: FASTP_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(qc_fastp_stream_task);
    cleanup_tasks.push(qc_fastp_err_task);

    let qc_fastp_out_stream = {
        let mut guard = qc_fastp_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: FASTP_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Tee the byte stream for counting
    let tee_count_input = ReceiverStream::new(qc_fastp_out_stream);
    let (mut tee_count_streams, tee_count_done_rx) = t_junction(
        tee_count_input,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "fastp_count".to_string(),
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "t_junction".to_string(),
            error: e.to_string(),
        })?;
    cleanup_receivers.push(tee_count_done_rx);

    let mut tee_count_streams_iter = tee_count_streams.into_iter();
    let qc_out_stream = tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let qc_count_stream = tee_count_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let (count_tx, count_result_rx) = oneshot::channel::<Result<u64, anyhow::Error>>();

    let count_task = tokio::spawn(async move {
        match stream_record_counter(qc_count_stream, false).await {
            Ok(count) => {
                eprintln!("Fastp output reads: {}", count);
                let _ = count_tx.send(Ok(count)); // Send the count result
                Ok(())
            }
            Err(e) => {
                let _ = count_tx.send(Err(e));
                Err(anyhow!("Failed to count fastp output reads"))
            },
        }
    });
    cleanup_tasks.push(count_task);

    Ok((ReceiverStream::new(qc_out_stream), cleanup_tasks, cleanup_receivers, count_result_rx))
}

/// Runs Kallisto
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
///
/// # Returns
///
async fn kallisto_quant(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    output_dir: PathBuf,
    paired: bool, // Ignored since we force single-end mode
    sample_base_buf: PathBuf,
) -> Result<
    (
        oneshot::Receiver<Vec<(String, f64)>>, // ERCC counts
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let sample_base = sample_base_buf.file_stem().unwrap_or_default().to_string_lossy().to_string();

    // Prepare KallistoConfig (always single-end mode)
    let kallisto_config = KallistoConfig {
        subcommand: KallistoSubcommand::Quant,
        subcommand_fields: HashMap::from([]),
        output_dir: output_dir.clone(),
        reproducible: false,
    };

    let kallisto_args = generate_cli(KALLISTO_TAG, &*config, Some(&kallisto_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    eprintln!("Kallisto QUANT: {:?}", kallisto_args);

    // Use stream_to_cmd for stdin streaming (treats interleaved as single-end)
    let (kallisto_child, stream_task, err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        KALLISTO_TAG,
        kallisto_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(stream_task);
    cleanup_tasks.push(err_task);
    let kallisto_child_arc = kallisto_child;

    // Capture kallisto stdout (plaintext abundance.tsv)
    let kallisto_out_stream = {
        let mut guard = kallisto_child_arc.lock().await;
        parse_child_output(
            &mut *guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: KALLISTO_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Split abundance stream for saving and ERCC filtering
    let (abundance_streams, abundance_done_rx) = t_junction(
        ReceiverStream::new(kallisto_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "kallisto_abundance_split".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;
    cleanup_receivers.push(abundance_done_rx);

    let mut abundance_streams_iter = abundance_streams.into_iter();
    let save_stream = abundance_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let ercc_stream = abundance_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let abundance_tsv_path = file_path_manipulator(
        &sample_base_buf,
        Some(&output_dir),
        None,
        Some("reads_per_transcript.kallisto.tsv"),
        "_",
    );
    let save_task = tokio::spawn(async move {
        stream_to_file(save_stream, abundance_tsv_path.clone()).await?;
        // Verify output file exists
        if !abundance_tsv_path.exists() {
            return Err(anyhow!("Kallisto failed to produce abundance.tsv"));
        }
        Ok(())
    });
    cleanup_tasks.push(save_task);

    // Process ERCC counts
    let (ercc_tx, ercc_rx) = oneshot::channel::<Vec<(String, f64)>>();
    let ercc_regex = Regex::new(r"^ERCC-").map_err(|e| PipelineError::Other(anyhow!("Regex error: {}", e)))?;
    let ercc_counts_path = file_path_manipulator(
        &sample_base_buf,
        Some(&output_dir),
        None,
        Some("ERCC_counts.tsv"),
        "_",
    );

    let ercc_task = tokio::spawn(async move {
        let mut ercc_counts = vec![];
        let mut ercc_file = File::create(&ercc_counts_path)
            .map_err(|e| anyhow!("Failed to create ERCC_counts.tsv: {}", e))?;
        writeln!(ercc_file, "target_id\test_counts")
            .map_err(|e| anyhow!("Failed to write ERCC_counts.tsv header: {}", e))?;
        let header = r"target_id\tlength\teff_length\est_counts\ttpm";
        let mut ercc_stream = ReceiverStream::new(ercc_stream);
        while let Some(ParseOutput::Bytes(line)) = ercc_stream.next().await {
            let line = String::from_utf8_lossy(&line).to_string();
            if line.trim() == header {
                continue; // Skip header
            }
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() >= 5 && ercc_regex.is_match(fields[0]) {
                let est_counts: f64 = fields[3].parse()
                    .map_err(|e| anyhow!("Failed to parse est_counts: {}", e))?;
                writeln!(ercc_file, "{}\t{}", fields[0], est_counts)
                    .map_err(|e| anyhow!("Failed to write ERCC counts: {}", e))?;
                ercc_counts.push((fields[0].to_string(), est_counts));
            }
        }
        ercc_tx.send(ercc_counts)
            .map_err(|_| anyhow!("Failed to send ERCC counts"))?;
        Ok(())
    });
    cleanup_tasks.push(ercc_task);

    Ok((
        ercc_rx,
        cleanup_tasks,
        cleanup_receivers,
    ))
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
    eprintln!("start");
    // External tools check
    check_versions(vec![
    BOWTIE2_TAG,
    HISAT2_TAG,
    KALLISTO_TAG
    ])
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    eprintln!("end");

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


    let (qc_fastp_out_stream, qc_cleanup_tasks, qc_cleanup_receivers, qc_count_result_rx) = fastp_qc(config.clone(), ercc_bt2_out_stream).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);
    cleanup_receivers.extend(qc_cleanup_receivers);


    // Split here as kallisto does not pass its input stream back out:
    let (mut kallisto_streams,  kallisto_split_done_rx) = t_junction(
        qc_fastp_out_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "kallisto_split".to_string(),
        None,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "t_junction".to_string(),
            error: e.to_string(),
        })?;
    cleanup_receivers.push(kallisto_split_done_rx);

    let mut kallisto_streams_iter = kallisto_streams.into_iter();
    let kallisto_bypass_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let kallisto_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let kallisto_stream = ReceiverStream::new(kallisto_stream);

    let (kallisto_ercc_rx, kallisto_cleanup_tasks, kallisto_cleanup_receivers) =
        kallisto_quant(
            config.clone(),
            kallisto_stream,
            out_dir.clone(),
            paired,
            no_ext_sample_base_buf.clone(),
        )
            .await?;
    cleanup_tasks.extend(kallisto_cleanup_tasks);
    cleanup_receivers.extend(kallisto_cleanup_receivers);

    // Test write out for the main stream until pipeline construction complete
    let test_write_task = tokio::spawn(stream_to_file(
        kallisto_bypass_stream,
        PathBuf::from("qc-test.fq"),
    ));
    test_write_task.await;


    // Retrieve ERCC counts
    let kallisto_ercc_counts = kallisto_ercc_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC counts receiver failed: {}", e)))?;
    eprintln!("ERCC counts: {:?}", kallisto_ercc_counts);


    let ercc_mapped_count = ercc_count_rx.await.map_err(|e| PipelineError::Other(anyhow::anyhow!("ERCC count receiver failed: {}", e)))?;


    let raw_count = join_with_error_handling(raw_count_task).await?;
    println!("Processed {} raw reads (additive from R1 and R2 if paired)", raw_count);

    let stats = join_with_error_handling(val_count_task).await?;
    println!("Processed {} validated, {} undersized, {} oversized reads",
             stats.validated, stats.undersized, stats.oversized);


    let _qc_fastp_read_count = match qc_count_result_rx.await {
        Ok(Ok(count)) => {
            eprintln!("Received fastp read count: {}", count);
            count
        }
        Ok(Err(e)) => {
            eprintln!("Error getting fastp read count: {}", e);
            0
        }
        Err(e) => {
            eprintln!("Count receiver dropped: {}", e);
            0
        }
    };


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
