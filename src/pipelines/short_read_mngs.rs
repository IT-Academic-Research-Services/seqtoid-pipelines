use crate::utils::fastx::SequenceRecord;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::fs::File;
use std::io::Write;
use std::io::BufReader;
use std::io::BufRead;
use regex::Regex;
use anyhow::{anyhow, Result};
use futures::future::try_join_all;
use tokio::fs;
use tokio::time::{sleep, Duration, timeout};
use tokio::sync::oneshot;
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;
use tokio::sync::Notify;
use tokio::fs::OpenOptions as TokioOpenOptions;
use crate::config::defs::{PipelineError, RunConfig, StreamDataType, ReadStats, MINIMAP2_TAG, BOWTIE2_TAG, SAMTOOLS_TAG, FASTP_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, MAFFT_TAG, SEQKIT_TAG, QUAST_TAG, HISAT2_TAG, SamtoolsSubcommand, KALLISTO_TAG, KallistoSubcommand, STAR_TAG};
use crate::utils::file::{file_path_manipulator, validate_file_inputs, write_byte_stream_to_file};
use crate::utils::fastx::{raw_read_count, read_fastq, stream_record_counter};
use crate::utils::streams::{t_junction, ParseOutput, join_with_error_handling, stream_to_cmd, parse_child_output, ChildStream, ParseMode, stream_to_file, read_child_output_to_vec, spawn_cmd, parse_fastq, ChannelReader, write_to_fifo};
use crate::utils::command::bowtie2::{Bowtie2Config, bowtie2_index_prep};
use crate::utils::command::{check_versions, generate_cli};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::fastp::FastpConfig;
use crate::utils::command::kallisto::KallistoConfig;
use crate::utils::command::hisat2::{Hisat2Config, hisat2_index_prep};
use crate::utils::streams::{deinterleave_fastq_stream_to_fifos};
use crate::utils::command::star::{StarConfig, star_index_prep};


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

/// bowtie2 filter function where the passing stream contains the unmapped reads
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream
/// * `bt2_index_path` - Path to Bowtie2 index.
/// * `paired` - Whether the input is paired-end.
/// * `bowtie2_options` - Additional Bowtie2 options as a HashMap (e.g., HashMap::from([("--very-sensitive-local".to_string(), None)])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - unmapped FASTQ stream.
/// - Optional receiver for the total mapped count (u64) if `count_mapped` is true.
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.
async fn bowtie2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    bt2_index_path: PathBuf,
    paired: bool,
    bowtie2_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // BT2
    let bt2_config_view = Bowtie2Config {
        bt2_index_path: bt2_index_path.clone(),
        option_fields: bowtie2_options,
    };

    let bt2_args = generate_cli(BOWTIE2_TAG, &config, Some(&bt2_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bt2_child, bt2_stream_task, bt2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),
        BOWTIE2_TAG,
        bt2_args,
        StreamDataType::JustBytes, // FASTQ bytes
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BOWTIE2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bt2_stream_task);
    cleanup_tasks.push(bt2_err_task);

    let bt2_out_stream = {
        let mut guard = bt2_child.lock().await;
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

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted (required for paired-end fastq extraction)
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        bt2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
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

    // Determine number of streams for t_junction (unmapped + optional BAM output + mapped count)
    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 }; // Unmapped + count + optional BAM

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = if num_tees > 1 {
        t_junction(
            bam_rx_stream,
            num_tees,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes, // BAM bytes
            "bowtie2_bam_split".to_string(),
            None,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?
    } else {
        // Single stream case (shouldn't happen since we always need unmapped + count)
        (vec![bam_rx_stream.into_inner()], oneshot::channel::<Result<(), anyhow::Error>>().1)
    };
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads (always performed)
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let mapped_flag = if paired { "-F13".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None), // Count
            (mapped_flag, None),
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
        bam_count_stream,
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

    let (count_tx, count_rx) = oneshot::channel::<u64>();

    let count_future = tokio::spawn(async move {
        let mut count_child = count_child_arc;
        let mut guard = count_child.lock().await;
        let count_lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        let _ = count_tx.send(mapped_count);
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Unmapped stream
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Extract unmapped FASTQ
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    Ok((ReceiverStream::new(unmapped_fastq_stream), count_rx, cleanup_tasks, cleanup_receivers))
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
async fn kallisto_quant(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    out_dir: PathBuf,
    paired: bool,
    sample_base_buf: PathBuf,
) -> Result<(oneshot::Receiver<KallistoResults>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let sample_base = sample_base_buf.file_stem().unwrap_or_default().to_string_lossy().to_string();
    let kallisto_out_dir = out_dir.join("kallisto");
    fs::create_dir_all(&kallisto_out_dir).await
        .map_err(|e| PipelineError::IOError(format!("Failed to create Kallisto output directory: {}", e)))?;

    let kallisto_index = config.args.kallisto_index.clone()
        .ok_or(PipelineError::MissingArgument("kallisto_index required".to_string()))?;
    let (ercc_tx, ercc_rx) = oneshot::channel::<KallistoResults>();

    // Convert raw byte stream to FASTQ
    let byte_rx = input_stream.into_inner();
    let byte_reader = ChannelReader::new(byte_rx);
    let fastq_rx = parse_fastq(byte_reader, config.base_buffer_size).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: "parse_fastq".to_string(),
            error: e.to_string(),
        })?;
    let fastq_stream = ReceiverStream::new(fastq_rx);

    eprintln!("Deinterleaving FASTQ stream for sample: {}", sample_base);
    let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        fastq_stream,
        &sample_base,
        paired,
    ).await?;

    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle {
        cleanup_tasks.push(r2_handle);
    }


    let kallisto_config = if paired {
        KallistoConfig {
            subcommand: KallistoSubcommand::Quant,
            subcommand_fields: HashMap::from([
                ("R1".to_string(), Some(r1_fifo.to_string_lossy().to_string())),
                ("R2".to_string(), Some(r2_fifo.to_string_lossy().to_string())),
            ]),
            output_dir: kallisto_out_dir.clone(),
            reproducible: false
        }
    }
    else {
        KallistoConfig {
            subcommand: KallistoSubcommand::Quant,
            subcommand_fields: HashMap::from([
                ("--single".to_string(), None),
                ("-l".to_string(), Some("200".to_string())),
                ("-s".to_string(), Some("20".to_string())),
                ("R1".to_string(), Some(r1_fifo.to_string_lossy().to_string())),
            ]),
            output_dir: kallisto_out_dir.clone(),
            reproducible: false
        }
    };

    let kallisto_args = generate_cli(KALLISTO_TAG, &config, Some(&kallisto_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;

    eprintln!("Spawning Kallisto with args: {:?}", kallisto_args);

    let (mut kallisto_child, kallisto_err_task) = spawn_cmd(
        config.clone(),
        KALLISTO_TAG,
        kallisto_args,
        config.args.verbose,
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KALLISTO_TAG.to_string(),
            error: e.to_string(),
        })?;

    cleanup_tasks.push(kallisto_err_task);

    // Await Kallisto exit
    let r1_fifo_clone = r1_fifo.clone();
    let r2_fifo_clone = r2_fifo.clone();
    let kallisto_exit_task = tokio::spawn(async move {
        let status = kallisto_child.wait().await
            .map_err(|e| anyhow!("Kallisto wait failed: {}", e))?;
        if !status.success() {
            return Err(anyhow!("Kallisto exited with code: {:?}", status.code()));
        }
        // Clean up FIFOs after Kallisto exits
        tokio::fs::remove_file(&r1_fifo_clone).await
            .map_err(|e| anyhow!("Failed to remove R1 FIFO {}: {}", r1_fifo_clone.display(), e))?;
        if paired {
            tokio::fs::remove_file(&r2_fifo_clone).await
                .map_err(|e| anyhow!("Failed to remove R2 FIFO {}: {}", r2_fifo_clone.display(), e))?;
        }
        Ok(())
    });
    cleanup_tasks.push(kallisto_exit_task);

    // Process Kallisto output (abundance.tsv)
    let kallisto_results_task = tokio::spawn(async move {
        let abundance_path = kallisto_out_dir.join("abundance.tsv");
        eprintln!("Starting to read abundance.tsv from: {}", abundance_path.display());
        let mut ercc_counts = Vec::new();
        let mut transcript_to_gene = Vec::new();

        match std::fs::File::open(&abundance_path) {
            Ok(file) => {
                let reader = BufReader::new(file);
                for line in reader.lines().skip(1) { // Skip header
                    let line = line.map_err(|e| anyhow!("Failed to read abundance.tsv: {}", e))?;
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() >= 5 {
                        let target_id = fields[0].to_string();
                        let est_counts = fields[4].parse::<f64>()
                            .map_err(|e| anyhow!("Failed to parse est_counts: {}", e))?;
                        ercc_counts.push((target_id.clone(), est_counts));
                        transcript_to_gene.push((target_id, format!("gene_{}", fields[0])));
                    }
                }
                eprintln!("Finished reading abundance.tsv: {} transcripts", ercc_counts.len());
            }
            Err(e) => {
                eprintln!("No abundance.tsv found at {}: {}, returning empty results", abundance_path.display(), e);
            }
        }

        Ok::<KallistoResults, anyhow::Error>(KallistoResults {
            ercc_counts,
            transcript_to_gene,
        })
    });


    let ercc_send_task = tokio::spawn(async move {
        let results = kallisto_results_task.await
            .map_err(|e| anyhow!("Kallisto results task failed: {}", e))??;
        ercc_tx.send(results)
            .map_err(|_| anyhow!("Failed to send Kallisto results"))?;
        eprintln!("Sent Kallisto results");
        Ok(())
    });
    cleanup_tasks.push(ercc_send_task);

    Ok((ercc_rx, cleanup_tasks, cleanup_receivers))
}
/// HISAT2 filter
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream.
/// * `hisat2_index_path` - Path to HISAT2 index.
/// * `paired` - Whether the input is paired-end.
/// * `hisat2_options` - Additional HISAT2 options as a HashMap (e.g., HashMap::from([("--no-spliced-alignment".to_string(), None)])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - Unmapped FASTQ stream.
/// - Optional receiver for the total mapped count (u64) if count is needed.
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.
async fn hisat2_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    hisat2_index_path: PathBuf,
    paired: bool,
    hisat2_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Create FIFOs for input (deinterleaved for paired, single FIFO for single-end)
    let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        input_stream,
        "hisat2_input",
        paired,
    ).await?;

    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle {
        cleanup_tasks.push(r2_handle);
    }

    // HISAT2 configuration with FIFOs
    let hisat2_config = Hisat2Config {
        hisat2_index_path: hisat2_index_path.clone(),
        option_fields: hisat2_options,
        r1_fifo: r1_fifo.clone(),
        r2_fifo: if paired { Some(r2_fifo.clone()) } else { None },
    };

    let hisat2_args = generate_cli(HISAT2_TAG, &config, Some(&hisat2_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: HISAT2_TAG.to_string(),
            error: e.to_string(),
        })?;

    eprintln!("hisat2 args {:?}", hisat2_args);

    let (mut hisat2_child, hisat2_err_task) = spawn_cmd(
        config.clone(),
        HISAT2_TAG,
        hisat2_args,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: HISAT2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(hisat2_err_task);

    let hisat2_out_stream = parse_child_output(
        &mut hisat2_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: HISAT2_TAG.to_string(),
            error: e.to_string(),
        })?;

    // Sort output to BAM (name-sorted for fastq extraction)
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted
            ("-u".to_string(), None), // Uncompressed BAM
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None), // Output to stdout
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        hisat2_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
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

    // Determine number of streams for t_junction (unmapped + count + optional BAM)
    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = if num_tees > 1 {
        t_junction(
            bam_rx_stream,
            num_tees,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes,
            "hisat2_bam_split".to_string(),
            None,
        )
            .await
            .map_err(|_| PipelineError::StreamDataDropped)?
    } else {
        (vec![bam_rx_stream.into_inner()], oneshot::channel::<Result<(), anyhow::Error>>().1)
    };
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        )
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let mapped_flag = if paired { "-F13".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None), // Count
            (mapped_flag, None),
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
        bam_count_stream,
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

    let (count_tx, count_rx) = oneshot::channel::<u64>();

    let count_future = tokio::spawn(async move {
        let mut count_child = count_child_arc;
        let mut guard = count_child.lock().await;
        let count_lines = read_child_output_to_vec(&mut guard, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        let _ = count_tx.send(mapped_count);
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Unmapped stream
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Extract unmapped FASTQ
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Cleanup FIFOs
    let r1_fifo_cleanup = r1_fifo.clone();
    let r2_fifo_cleanup = r2_fifo.clone();
    let fifo_cleanup_task = tokio::spawn(async move {
        tokio::fs::remove_file(&r1_fifo_cleanup).await.ok();
        if paired {
            tokio::fs::remove_file(&r2_fifo_cleanup).await.ok();
        }
        Ok(())
    });
    cleanup_tasks.push(fifo_cleanup_task);

    Ok((ReceiverStream::new(unmapped_fastq_stream), count_rx, cleanup_tasks, cleanup_receivers))
}


/// STAR filter
///
/// # Arguments
///
/// * `config` - RunConfig struct from main.
/// * `input_stream` - Raw byte FASTQ stream (interleaved).
/// * `star_index_dir` - Path to STAR index directory.
/// * `paired` - Whether the input is paired-end.
/// * `star_options` - Additional STAR options as a HashMap (e.g., HashMap::from([("--outFilterMultimapNmax".to_string(), Some("5".to_string()))])).
/// * `output_bam_path` - Optional path to save the aligned BAM file (name-sorted).
///
/// # Returns
/// Tuple:
/// - unmapped FASTQ stream.
/// - Receiver for the total mapped count (u64).
/// - Vector of cleanup tasks.
/// - Vector of cleanup receivers.
async fn star_filter(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    star_index_path: PathBuf,
    paired: bool,
    star_options: HashMap<String, Option<String>>,
    output_bam_path: Option<PathBuf>,
) -> Result<(ReceiverStream<ParseOutput>, oneshot::Receiver<u64>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Deinterleave to FIFOs
    let (r1_fifo, r2_fifo, deinterleave_handle, r1_write_handle, r2_write_handle) = deinterleave_fastq_stream_to_fifos(
        config.clone(),
        input_stream,
        "star_filter",
        paired,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: "deinterleave".to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(deinterleave_handle);
    cleanup_tasks.push(r1_write_handle);
    if let Some(r2_handle) = r2_write_handle {
        cleanup_tasks.push(r2_handle);
    }

    // STAR config
    let star_config = StarConfig {
        star_index_dir: star_index_prep(&star_index_path, &std::env::current_dir().unwrap())?,
        option_fields: star_options,
        r1_fifo: r1_fifo.clone(),
        r2_fifo: if paired { Some(r2_fifo.clone()) } else { None },
    };

    let star_args = generate_cli(STAR_TAG, &config, Some(&star_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: STAR_TAG.to_string(),
            error: e.to_string(),
        })?;

    // Spawn STAR (no stdin stream, uses FIFOs)
    let (mut star_child, star_err_task) = spawn_cmd(
        config.clone(),
        STAR_TAG,
        star_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: STAR_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(star_err_task);

    let star_out_stream = parse_child_output(
        &mut star_child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        config.base_buffer_size,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: STAR_TAG.to_string(),
        error: e.to_string(),
    })?;

    // Sort, output uncompressed BAM
    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-n".to_string(), None), // Name-sorted (required for paired-end fastq extraction)
            ("-u".to_string(), None),
            ("-O".to_string(), Some("bam".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        star_out_stream,
        SAMTOOLS_TAG,
        samtools_sort_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(samtools_sort_task);
    cleanup_tasks.push(samtools_sort_err_task);

    let samtools_sort_out_stream = {
        let mut guard = samtools_sort_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?
    };

    // number of streams for t_junction (unmapped + count + optional BAM)
    let num_tees = 2 + if output_bam_path.is_some() { 1 } else { 0 };

    let bam_rx_stream = ReceiverStream::new(samtools_sort_out_stream);

    let (bam_streams, bam_done_rx) = if num_tees > 1 {
        t_junction(
            bam_rx_stream,
            num_tees,
            config.base_buffer_size,
            config.args.stall_threshold,
            None,
            100,
            StreamDataType::JustBytes,
            "star_bam_split".to_string(),
            None,
        ).await.map_err(|_| PipelineError::StreamDataDropped)?
    } else {
        (vec![bam_rx_stream.into_inner()], oneshot::channel::<Result<(), anyhow::Error>>().1)
    };
    cleanup_receivers.push(bam_done_rx);

    let mut bam_streams_iter = bam_streams.into_iter();

    // Optional: Write BAM to file
    if let Some(bam_path) = output_bam_path {
        let stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
        let bam_write_task = write_byte_stream_to_file(
            &bam_path,
            ReceiverStream::new(stream),
            Some(config.base_buffer_size),
        ).await.map_err(|e| PipelineError::IOError(e.to_string()))?;
        cleanup_tasks.push(bam_write_task);
    }

    // Count total mapped reads
    let bam_count_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let mapped_flag = if paired { "-F13".to_string() } else { "-F4".to_string() };
    let samtools_count_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-c".to_string(), None), // Count
            (mapped_flag, None),
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
        bam_count_stream,
        SAMTOOLS_TAG,
        samtools_count_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(count_stream_task);
    cleanup_tasks.push(count_err_task);

    let (count_tx, count_rx) = oneshot::channel::<u64>();

    let count_future = tokio::spawn(async move {
        let mut count_child = count_child_arc.lock().await;
        let count_lines = read_child_output_to_vec(&mut count_child, ChildStream::Stdout).await?;
        let mapped_count: u64 = count_lines.get(0).unwrap_or(&"0".to_string()).trim().parse()?;
        let _ = count_tx.send(mapped_count);
        Ok(())
    });
    cleanup_tasks.push(count_future);

    // Unmapped stream
    let unmapped_stream = bam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Extract unmapped FASTQ
    let unmapped_flag = if paired { "-f13".to_string() } else { "-f4".to_string() };
    let samtools_fastq_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([
            (unmapped_flag, None),
            ("-".to_string(), None), // Output to stdout (interleaved for paired)
        ]),
    };
    let samtools_fastq_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_fastq_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut fastq_child, fastq_stream_task, fastq_err_task) = stream_to_cmd(
        config.clone(),
        unmapped_stream,
        SAMTOOLS_TAG,
        samtools_fastq_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: SAMTOOLS_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(fastq_stream_task);
    cleanup_tasks.push(fastq_err_task);

    let unmapped_fastq_stream = {
        let mut guard = fastq_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Fastq,
            config.base_buffer_size,
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?
    };

    // Cleanup FIFOs
    let r1_fifo_cleanup = r1_fifo.clone();
    let r2_fifo_cleanup = r2_fifo.clone();
    let fifo_cleanup_task = tokio::spawn(async move {
        tokio::fs::remove_file(&r1_fifo_cleanup).await.ok();
        if paired {
            tokio::fs::remove_file(&r2_fifo_cleanup).await.ok();
        }
        Ok(())
    });
    cleanup_tasks.push(fifo_cleanup_task);

    Ok((ReceiverStream::new(unmapped_fastq_stream), count_rx, cleanup_tasks, cleanup_receivers))
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
    check_versions(vec![BOWTIE2_TAG, STAR_TAG, KALLISTO_TAG])
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    // Check required files
    let host_bowtie2_index: String = config.args.host_bowtie2_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_bowtie2_index is required".to_string()))?;
    let host_star_index: String = config.args.host_star_index.clone()
        .ok_or_else(|| PipelineError::MissingArgument("host_star_index is required".to_string()))?;

    let (file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;
    let paired = file2_path.is_some();

    // Input Validation
    let (val_out_stream, validate_cleanup_tasks, validate_cleanup_receivers, raw_count_task, val_count_task) = validate_input(
        config.clone(),
        file1_path,
        file2_path,
        no_ext_sample_base_buf.clone(),
        &out_dir,
    ).await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);

    // ERCC bt2 filtering and count
    let ercc_bt2_index_path = bowtie2_index_prep(&config.args.ercc_bowtie2_index, &cwd)?;
    let ercc_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
    let (ercc_bt2_out_stream, ercc_count_rx, ercc_bt2_cleanup_tasks, ercc_bt2_cleanup_receivers) = bowtie2_filter(
        config.clone(),
        val_out_stream,
        ercc_bt2_index_path,
        paired,
        ercc_bt2_options,
        None,
    ).await?;
    cleanup_tasks.extend(ercc_bt2_cleanup_tasks);
    cleanup_receivers.extend(ercc_bt2_cleanup_receivers);

    let (qc_fastp_out_stream, qc_cleanup_tasks, qc_cleanup_receivers, qc_count_result_rx) = fastp_qc(
        config.clone(),
        ercc_bt2_out_stream,
    ).await?;
    cleanup_tasks.extend(qc_cleanup_tasks);
    cleanup_receivers.extend(qc_cleanup_receivers);

    // Split for Kallisto and Bowtie2
    let (mut kallisto_streams, kallisto_split_done_rx) = t_junction(
        qc_fastp_out_stream,
        2,
        config.base_buffer_size * 10,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "kallisto_split".to_string(),
        None,
    ).await?;
    cleanup_receivers.push(kallisto_split_done_rx);

    let mut kallisto_streams_iter = kallisto_streams.into_iter();
    let kallisto_bypass_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let kallisto_stream = kallisto_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let kallisto_stream = ReceiverStream::new(kallisto_stream);

    let (kallisto_ercc_rx, kallisto_cleanup_tasks, kallisto_cleanup_receivers) = kallisto_quant(
        config.clone(),
        kallisto_stream,
        out_dir.clone(),
        paired,
        no_ext_sample_base_buf.clone(),
    ).await?;
    cleanup_tasks.extend(kallisto_cleanup_tasks);
    cleanup_receivers.extend(kallisto_cleanup_receivers);

    // Host filtering: bt2
    let host_bt2_index_path = bowtie2_index_prep(host_bowtie2_index, &cwd)?;
    let host_bt2_options = HashMap::from([("--very-sensitive-local".to_string(), None)]);
    let (host_bt2_out_stream, host_bt2_count_rx, host_bt2_cleanup_tasks, host_bt2_cleanup_receivers) = bowtie2_filter(
        config.clone(),
        ReceiverStream::new(kallisto_bypass_stream),
        host_bt2_index_path,
        paired,
        host_bt2_options,
        None,
    ).await?;
    cleanup_tasks.extend(host_bt2_cleanup_tasks);
    cleanup_receivers.extend(host_bt2_cleanup_receivers);

    // Host filtering: STAR
    let host_star_index_path = star_index_prep(host_star_index, &cwd)?;
    let host_star_options = HashMap::from([]);
    let (host_star_out_stream, host_star_count_rx, host_star_cleanup_tasks, host_star_cleanup_receivers) = star_filter(
        config.clone(),
        host_bt2_out_stream,
        host_star_index_path,
        paired,
        host_star_options,
        None,
    ).await?;
    cleanup_tasks.extend(host_star_cleanup_tasks);
    cleanup_receivers.extend(host_star_cleanup_receivers);

    // Test write out for the main stream until pipeline construction complete
    let test_write_task = tokio::spawn(stream_to_file(
        host_star_out_stream.into_inner(),
        PathBuf::from("test.fq"),
    ));
    test_write_task.await;

    // Results retrieval
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

    let kallisto_ercc_counts = kallisto_ercc_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC counts receiver failed: {}", e)))?;
    // eprintln!("ERCC counts: {:?}", kallisto_ercc_counts);

    let ercc_mapped_count = ercc_count_rx.await
        .map_err(|e| PipelineError::Other(anyhow!("ERCC count receiver failed: {}", e)))?;

    let host_bt2_counts = host_bt2_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host bt2 counts receiver failed: {}", e)))?;
    eprintln!("Host bt2 counts: {:?}", host_bt2_counts);

    let host_star_counts = host_star_count_rx
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Host star counts receiver failed: {}", e)))?;
    eprintln!("Host star counts: {:?}", host_star_counts);

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
