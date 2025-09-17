use tokio::io::{BufWriter, AsyncWriteExt};
use std::fs::File;
use std::sync::Arc;
use std::io::Write;
use std::collections::HashMap;
use tokio_stream::StreamExt;
use crate::utils::streams::ParseOutput;
use std::path::PathBuf;
use anyhow::{anyhow, Result};
use tempfile::NamedTempFile;
use crate::cli::Technology;
use tokio::task::JoinHandle;
use std::time::Instant;
use tokio::fs;
use tokio::fs::File as TokioFile;
use tokio_stream::wrappers::ReceiverStream;
use tokio::sync::mpsc::Receiver;
use serde::Serialize;
use tokio::sync::{mpsc, oneshot};
use futures::future::try_join_all;
use fxhash::FxHashMap as FxHashMap;
use crate::utils::command::{generate_cli, check_versions};
use crate::utils::file::{extension_remover, file_path_manipulator, write_parse_output_to_temp_fifo, write_vecu8_to_file, validate_file_inputs};
use crate::utils::fastx::{read_fastq, r1r2_base, parse_and_filter_fastq_id, concatenate_paired_reads, parse_byte_stream_to_fastq};
use crate::utils::streams::{t_junction, stream_to_cmd, parse_child_output, ChildStream, ParseMode, stream_to_file, spawn_cmd, parse_fastq, parse_bytes, y_junction};
use crate::config::defs::{PipelineError, StreamDataType, PIGZ_TAG, FASTP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, SamtoolsSubcommand, KRAKEN2_TAG, BCFTOOLS_TAG, BcftoolsSubcommand, MAFFT_TAG, QUAST_TAG, SEQKIT_TAG, SeqkitSubcommand};
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::kraken2::Kraken2Config;
use crate::utils::command::bcftools::BcftoolsConfig;
use crate::utils::command::seqkit::SeqkitConfig;
use crate::utils::db::{get_index, retrieve_h5_seq};
use crate::config::defs::RunConfig;
use crate::utils::command::quast::QuastConfig;
use crate::utils::stats::{parse_samtools_stats, parse_samtools_depth, compute_depth_stats, parse_seqkit_stats, parse_ercc_stats, compute_allele_counts, compute_coverage_bins};
use crate::utils::vcf::count_variants_from_bcftools_stats;
use crate::utils::plotting::plot_depths;


#[derive(Serialize)]
struct Stats {
    sample_name: String,
    depth_avg: f64,
    depth_q25: f64,
    depth_q50: f64,
    depth_q75: f64,
    depth_frac_above_10x: f64,
    depth_frac_above_25x: f64,
    depth_frac_above_50x: f64,
    depth_frac_above_100x: f64,
    allele_counts: HashMap<char, u64>,
    total_reads: u64,
    mapped_reads: u64,
    mapped_paired: Option<u64>,
    paired_inward: Option<u64>,
    paired_outward: Option<u64>,
    paired_other_orientation: Option<u64>,
    ercc_mapped_reads: Option<u64>,
    ercc_mapped_paired: Option<u64>,
    ref_snps: u64,
    ref_mnps: u64,
    ref_indels: u64,
    n_actg: u64,
    n_missing: u64,
    n_gap: u64,
    n_ambiguous: u64,
    coverage_breadth: f64,
    max_aligned_length: usize,
    total_length: usize,
    coverage_bin_size: f64,
    coverage: Vec<(usize, f64, f64, u8, u8)>,
}


async fn prepare_reference_and_index<'a>(
    config: &RunConfig,
    ref_db_path: Option<PathBuf>,
    ram_temp_dir: &PathBuf,
    h5_index: Option<&'a FxHashMap<[u8; 24], u64>>,
    accession: Option<String>,
    sequence: Option<String>,
    index_path: Option<String>,
    ref_type: &str,
) -> Result<
    (
        Option<PathBuf>,             // FASTA path (None if only index provided)
        PathBuf,                     // Index path (.mmi)
        Option<NamedTempFile>,       // FASTA temp file (None if only index)
        Option<NamedTempFile>,       // Index temp file (if created/copied)
        Vec<JoinHandle<Result<(), anyhow::Error>>>, // Tasks
    ),
    PipelineError,
> {
    let mut tasks = Vec::new();
    let mut ref_fasta_path: Option<PathBuf> = None;
    let mut ref_temp: Option<NamedTempFile> = None;

    // Check if only index is provided
    let (index_path, index_temp) = match index_path {
        Some(index) => {
            let index_path = PathBuf::from(&index);
            if !index_path.exists() {
                return Err(PipelineError::FileNotFound(index_path));
            }
            if index_path.extension().map_or(true, |ext| ext != "mmi") {
                return Err(PipelineError::InvalidConfig(format!(
                    "{} index must have .mmi extension: {}",
                    ref_type,
                    index_path.display()
                )));
            }
            // Copy index to ram_temp_dir
            let index_temp = NamedTempFile::with_suffix_in(".mmi", ram_temp_dir)
                .map_err(|e| PipelineError::Other(e.into()))?;
            let index_temp_path = index_temp.path().to_path_buf();
            let index_temp_path_clone = index_temp_path.clone();
            let ref_type_owned = ref_type.to_string();
            let copy_task = tokio::spawn(async move {
                fs::copy(&index_path, &index_temp_path_clone)
                    .await
                    .map_err(|e| anyhow!("Failed to copy {} index: {}", ref_type_owned, e))?;
                Ok(())
            });
            tasks.push(copy_task);
            (index_temp_path, Some(index_temp))
        }
        None => {
            // Require sequence or accession for FASTA
            let (_accession, seq) = retrieve_h5_seq(
                accession.clone(),
                sequence.clone(),
                ref_db_path.as_ref(),
                h5_index,
            )
                .await
                .map_err(|e| PipelineError::ReferenceRetrievalFailed(e.to_string()))?;

            let ref_temp_file = NamedTempFile::with_suffix_in(".fasta", ram_temp_dir)
                .map_err(|e| PipelineError::Other(e.into()))?;
            ref_fasta_path = Some(ref_temp_file.path().to_path_buf());
            ref_temp = Some(ref_temp_file);
            let ref_write_task = write_vecu8_to_file(Arc::new(seq.clone()), ref_fasta_path.as_ref().unwrap(), config.base_buffer_size)
                .await
                .map_err(|e| PipelineError::Other(e.into()))?;
            tasks.push(ref_write_task);

            // Check for .mmi in sequence directory
            let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
            let sequence_path_buf = file_path_manipulator(&PathBuf::from(ref_fasta_path.as_ref().unwrap()), Some(&cwd), None, None, "");
            let (no_ext_sequence_path, _) = extension_remover(&sequence_path_buf);
            let mut index_path = no_ext_sequence_path.with_extension("mmi");
            if !index_path.exists() {
                // No .mmi found; create one
                eprintln!("No {} index found at {}; creating new index", ref_type, index_path.display());
                let index_temp = NamedTempFile::with_suffix_in(".mmi", ram_temp_dir)
                    .map_err(|e| PipelineError::Other(e.into()))?;
                index_path = index_temp.path().to_path_buf();
                let minimap2_args = vec![
                    "-d".to_string(),
                    index_path.to_string_lossy().to_string(),
                    ref_fasta_path.as_ref().unwrap().to_string_lossy().to_string(),
                ];
                let (mut child, err_task) = spawn_cmd(
                    Arc::new(config.clone()),
                    MINIMAP2_TAG,
                    minimap2_args,
                    config.args.verbose,
                )
                    .await
                    .map_err(|e| PipelineError::ToolExecution {
                        tool: MINIMAP2_TAG.to_string(),
                        error: e.to_string(),
                    })?;
                let index_task = tokio::spawn(async move {
                    let status = child.wait().await?;
                    if !status.success() {
                        return Err(anyhow!("minimap2 index creation failed with exit code: {:?}", status.code()));
                    }
                    Ok(())
                });
                tasks.push(index_task);
                tasks.push(err_task);
                (index_path, Some(index_temp))
            } else {
                // Found .mmi; copy to ram_temp_dir
                eprintln!("Found {} index at {}; copying to temp dir", ref_type, index_path.display());
                let index_temp = NamedTempFile::with_suffix_in(".mmi", ram_temp_dir)
                    .map_err(|e| PipelineError::Other(e.into()))?;
                let index_temp_path = index_temp.path().to_path_buf();
                let index_temp_path_clone = index_temp_path.clone();
                let ref_type_owned = ref_type.to_string();
                let copy_task = tokio::spawn(async move {
                    fs::copy(&index_path, &index_temp_path_clone)
                        .await
                        .map_err(|e| anyhow!("Failed to copy {} index: {}", ref_type_owned, e))?;
                    Ok(())
                });
                tasks.push(copy_task);
                (index_temp_path, Some(index_temp))
            }
        }
    };

    Ok((ref_fasta_path, index_path, ref_temp, index_temp, tasks))
}


async fn validate_input(
    config: Arc<RunConfig>,
    file1_path: PathBuf,
    file2_path: Option<PathBuf>,
    sample_base_buf: PathBuf,
    out_dir: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let validated_interleaved_file_path = file_path_manipulator(
        &PathBuf::from(&sample_base_buf),
        Some(out_dir),
        None,
        Some("validated.fq.gz"),
        "_",
    );
    
    let rx = read_fastq(
        file1_path,
        file2_path,
        Some(config.args.technology.clone()),
        config.args.max_reads as u64,
        config.args.min_read_len,
        config.args.max_read_len,
        config.base_buffer_size,  // Use as chunk_size
    )
        .map_err(|e| PipelineError::InvalidFastqFormat(e.to_string()))?;

    let val_rx_stream = ReceiverStream::new(rx);

    // Split the byte stream for fastp and quast write
    let (val_streams, val_done_rx) = t_junction(
        val_rx_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,  // Semantically FASTQ bytes
        "validate_input".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    if val_streams.len() < 2 {
        return Err(PipelineError::EmptyStream);
    }

    // Consume Vec to get Receivers
    let mut val_streams_iter = val_streams.into_iter();
    let val_fastp_out_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let val_quast_out_stream = val_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Write validated interleaved FASTQ to file using byte chunks
    let val_quast_write_task = tokio::spawn(async move {
        let mut file = BufWriter::new(
            TokioFile::create(&validated_interleaved_file_path)
                .await
                .map_err(|e| anyhow!("Failed to create validated FASTQ file: {}", e))?,
        );
        let mut stream = ReceiverStream::new(val_quast_out_stream);
        let mut total_bytes = 0;

        while let Some(ParseOutput::Bytes(bytes)) = stream.next().await {
            total_bytes += bytes.len() as u64;
            file.write_all(&bytes)
                .await
                .map_err(|e| anyhow!("Failed to write to validated FASTQ: {}", e))?;
        }
        file.flush()
            .await
            .map_err(|e| anyhow!("Failed to flush validated FASTQ: {}", e))?;

        // Log for correctness (no silent drops)
        // eprintln!("validate_input: Wrote {} bytes to {}", total_bytes, validated_interleaved_file_path.display());
        Ok(())
    });

    Ok((
        ReceiverStream::new(val_fastp_out_stream),
        vec![val_quast_write_task],
        vec![val_done_rx]
    ))
}



async fn align_to_host(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>, // FASTQ raw byte stream from fastp
    host_index_path: PathBuf,  // minimap2 .mmi index file
    no_host_file_path: PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, ReceiverStream<ParseOutput>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    // Generate minimap2 args for host alignment
    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(host_index_path)))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut minimap2_child, minimap2_stream_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),  // Convert to Receiver<ParseOutput> for streaming
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::JustBytes,  // FASTQ bytes
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_stream_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    // Samtools view to filter unmapped reads (-f 4), output as uncompressed BAM (-b -u)
    let samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-f".to_string(), Some("4".to_string())), // Unmapped reads
            ("--no-PG".to_string(), None),
            ("-h".to_string(), None),   // Include header
            ("-b".to_string(), None),   // BAM output
            ("-u".to_string(), None),   // Uncompressed BAM
        ]),
    };
    let samtools_args_view = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_child_view, samtools_task_view, samtools_err_task_view) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
        SAMTOOLS_TAG,
        samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_view);
    cleanup_tasks.push(samtools_err_task_view);

    let samtools_out_stream_view = {
        let mut guard = samtools_child_view.lock().await;
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

    // Samtools fastq to convert uncompressed BAM to FASTQ bytes
    let samtools_config_fastq = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Fastq,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let samtools_args_fastq = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_fastq))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_child_fastq, samtools_task_fastq, samtools_err_task_fastq) = stream_to_cmd(
        config.clone(),
        samtools_out_stream_view,
        SAMTOOLS_TAG,
        samtools_args_fastq,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_fastq);
    cleanup_tasks.push(samtools_err_task_fastq);

    let samtools_out_stream_fastq = {
        let mut guard = samtools_child_fastq.lock().await;
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
    let samtools_out_stream_fastq = ReceiverStream::new(samtools_out_stream_fastq);

    // Split FASTQ byte stream for output, file write, count, and check
    let (host_streams, host_done_rx) = t_junction(
        samtools_out_stream_fastq,
        4,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "align_to_host".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(host_done_rx);

    // Consume Vec to get Receivers
    let mut streams_iter = host_streams.into_iter();
    let no_host_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_count_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let no_host_check_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Check task for host-removed empty FASTQ (adapted for bytes)
    let (check_tx, mut check_rx) = mpsc::channel(1);
    let check_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(no_host_check_stream);
        let mut count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(bytes) => {
                    if !bytes.is_empty() {
                        count = 1;
                        // Spawn background drain for remaining items to avoid blocking
                        tokio::spawn(async move {
                            while stream.next().await.is_some() {}
                            Ok::<(), anyhow::Error>(())
                        });
                        break; // Exit loop early to unblock sender
                    }
                }
                _ => return Err(anyhow!("Unexpected item type in no_host_check_stream")),
            }
        }
        check_tx
            .send(count)
            .await
            .map_err(|e| anyhow!("Failed to send check count: {}", e))?;
        Ok(())
    });
    cleanup_tasks.push(check_task);

    let check_count = check_rx.recv().await.ok_or(PipelineError::EmptyStream)?;
    if check_count == 0 {
        let drain_rxs = vec![no_host_output_stream, no_host_file_stream, no_host_count_stream];
        let drain_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = drain_rxs
            .into_iter()
            .map(|rx| {
                tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(rx);
                    while stream.next().await.is_some() {}
                    Ok(())
                })
            })
            .collect();
        cleanup_tasks.extend(drain_tasks);
        return Err(PipelineError::EmptyStream);
    }

    // Pigz for compression
    let pigz_args = generate_cli(PIGZ_TAG, &config, None)
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (mut pigz_child, pigz_stream_task, pigz_err_task) = stream_to_cmd(
        config.clone(),
        no_host_file_stream,
        PIGZ_TAG,
        pigz_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(pigz_stream_task);
    cleanup_tasks.push(pigz_err_task);

    let pigz_out_stream = {
        let mut guard = pigz_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: PIGZ_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    let pigz_write_task = tokio::spawn(stream_to_file(pigz_out_stream, no_host_file_path));
    cleanup_tasks.push(pigz_write_task);

    Ok((ReceiverStream::new(no_host_output_stream), ReceiverStream::new(no_host_count_stream), cleanup_tasks, cleanup_receivers))
}



async fn process_ercc(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>, // FASTQ byte stream
    ercc_index_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base: &str,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = Vec::new();
    let mut cleanup_receivers = Vec::new();

    let ercc_stats_file_path = file_path_manipulator(
        &PathBuf::from(no_ext_sample_base),
        Some(out_dir),
        None,
        Some("ercc_stats.txt"),
        "_",
    );

    let (ercc_streams, ercc_done_rx) = t_junction(
        input_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "process_ercc_bypass".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    let mut streams_iter = ercc_streams.into_iter();
    let bypass_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let alignment_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    cleanup_receivers.push(ercc_done_rx);

    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&ercc_index_path))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut minimap2_child, minimap2_stream_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        alignment_stream, // Use alignment_stream, not input_stream
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::JustBytes, // FASTQ bytes
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_stream_task);
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };


    // Split SAM stream for check and further processing
    let (sam_streams, sam_done_rx) = t_junction(
        ReceiverStream::new(minimap2_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "process_ercc_empty_check".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(sam_done_rx);
    let mut sam_streams_iter = sam_streams.into_iter();
    let sam_check_stream = sam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_view_stream = sam_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    // Check task for empty ERCC alignments
    let (count_tx, mut count_rx) = mpsc::channel(1);
    let count_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(sam_check_stream);
        let mut count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(bytes) => {
                    let line_str = String::from_utf8_lossy(&bytes);
                    if !line_str.starts_with('@') && !line_str.trim().is_empty() {
                        count = 1;
                        // Spawn background drain to avoid blocking
                        tokio::spawn(async move {
                            while stream.next().await.is_some() {}

                            Ok::<(), anyhow::Error>(())
                        });
                        break; // Exit early to unblock sender
                    }
                }
                _ => return Err(anyhow!("Unexpected item type in sam_check_stream")),
            }
        }
        count_tx
            .send(count)
            .await
            .map_err(|e| anyhow!("Failed to send SAM alignment count: {}", e))?;
        Ok::<(), anyhow::Error>(())
    });
    cleanup_tasks.push(count_task);

    let alignment_count = count_rx.recv().await.ok_or(PipelineError::EmptyStream)?;

    if alignment_count == 0 {
        // Case of no ERCC spike-ins
        let (dummy_tx, dummy_rx) = mpsc::channel(config.base_buffer_size);
        drop(dummy_tx); // Ensure stream is empty but valid
        let drain_task = tokio::spawn(async move {
            let mut stream = ReceiverStream::new(bypass_output_stream);
            while stream.next().await.is_some() {}
            Ok(())
        });
        cleanup_tasks.push(drain_task);

        let zero_stats = HashMap::from([
            ("ercc_mapped_reads".to_string(), 0),
            ("ercc_mapped_paired".to_string(), 0),
        ]);
        let stats_write_task = tokio::spawn(async move {
            let mut file = File::create(&ercc_stats_file_path)?;
            writeln!(file, "ercc_mapped_reads: 0\nercc_mapped_paired: 0")?;
            Ok(())
        });
        cleanup_tasks.push(stats_write_task);
        return Ok((
            ReceiverStream::new(dummy_rx),
            Some(tokio::spawn(async { Ok(zero_stats) })),
            cleanup_tasks,
            cleanup_receivers,
        ));
    }

    // Non-zero case: process stats
    let samtools_config_view = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("--no-PG".to_string(), None),
            ("-b".to_string(), None), // BAM output
            ("-u".to_string(), None), // Uncompressed BAM
        ]),
    };
    let samtools_args_view = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_config_view))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_child_view, samtools_task_view, samtools_err_task_view) = stream_to_cmd(
        config.clone(),
        sam_view_stream,
        SAMTOOLS_TAG,
        samtools_args_view,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_task_view);
    cleanup_tasks.push(samtools_err_task_view);

    let samtools_out_stream_view = {
        let mut guard = samtools_child_view.lock().await;
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


    let stats_samtools_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Stats,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let stats_samtools_args = generate_cli(SAMTOOLS_TAG, &config, Some(&stats_samtools_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut stats_samtools_child, stats_samtools_task, stats_samtools_err_task) = stream_to_cmd(
        config.clone(),
        samtools_out_stream_view,
        SAMTOOLS_TAG,
        stats_samtools_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(stats_samtools_task);
    cleanup_tasks.push(stats_samtools_err_task);

    let stats_samtools_out_stream = {
        let mut guard = stats_samtools_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Lines,
            config.base_buffer_size / 2,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };


    let (stats_streams, stats_done_rx) = t_junction(
        ReceiverStream::new(stats_samtools_out_stream),
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "process_ercc_stats".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(stats_done_rx);
    let mut stats_streams_iter = stats_streams.into_iter();
    let stats_file_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let stats_parse_stream = stats_streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let stats_write_task = tokio::spawn(stream_to_file(stats_file_stream, ercc_stats_file_path));
    cleanup_tasks.push(stats_write_task);

    let ercc_stats_task = Some(tokio::spawn(parse_ercc_stats(stats_parse_stream)));

    Ok((
        ReceiverStream::new(bypass_output_stream),
        ercc_stats_task,
        cleanup_tasks,
        cleanup_receivers,
    ))
}


/// Filters input FASTQ stream through Kraken2
///
/// # Arguments
///
/// * `config` - reading stream
/// * `input_stream` - ReceiverStream<ParseOutput>,
/// * 'targer_ref_path' -
///
/// # Returns
/// samtools_sort_out_stream: Result<(ReceiverStream<ParseOutput>, <--- FASTQ uncompressed
/// cleanup_tasks
/// quast_write_tasks
///
async fn filter_with_kraken(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,
    target_ref_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
    target_taxid: &str,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks = vec![];

    // Parse the byte stream into Fastq records
    let (parse_rx, parse_task) = parse_byte_stream_to_fastq(
        input_stream.into_inner(), // Convert ReceiverStream to Receiver
        config.base_buffer_size,
        config.args.stall_threshold,
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    cleanup_tasks.push(parse_task);

    // Concatenate paired-end reads
    let (concat_stream, concat_task) = concatenate_paired_reads(
        ReceiverStream::new(parse_rx),
        config.base_buffer_size,
        config.args.stall_threshold,
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    cleanup_tasks.push(concat_task);


    let final_compressed_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("classified_filtered.fq.gz"),
        "_",
    );

    // Write input stream to temp FIFO for Kraken2
    let (kraken_query_write_task, kraken_query_pipe_path) = write_parse_output_to_temp_fifo(
        concat_stream,
        Some(config.base_buffer_size),
        Some(".fq"),
        Some(&config.ram_temp_dir),
    )
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    cleanup_tasks.push(kraken_query_write_task);

    // Create temp FIFO for Kraken2 classified output
    let mut builder = tempfile::Builder::new();
    builder.suffix(".fq");
    let temp_name = builder
        .tempfile_in(&config.ram_temp_dir)
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create temp file in {}: {}", config.ram_temp_dir.display(), e)))?;
    let kraken2_pipe_path = temp_name.path().to_path_buf();
    temp_name
        .close()
        .map_err(|e| PipelineError::Other(anyhow!("Failed to close temp file: {}", e)))?;

    let _ = tokio::fs::remove_file(&kraken2_pipe_path).await;
    tokio::process::Command::new("mkfifo")
        .arg(&kraken2_pipe_path)
        .status()
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to create named pipe: {}", e)))?;

    // Run Kraken2 in single-end mode
    let kraken2_config = Kraken2Config {
        report_path: file_path_manipulator(
            no_ext_sample_base_buf,
            Some(out_dir),
            None,
            Some("kraken2_report.txt"),
            "_",
        ),
        classified_path: kraken2_pipe_path.clone(),
        fastq_path: kraken_query_pipe_path.clone(),
    };
    let kraken2_args = generate_cli(KRAKEN2_TAG, &config, Some(&kraken2_config))?;

    let (mut kraken2_child, kraken2_err_task) = spawn_cmd(config.clone(), KRAKEN2_TAG, kraken2_args, config.args.verbose)
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: KRAKEN2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(kraken2_err_task);

    let kraken2_wait_task = tokio::spawn(async move {
        let status = kraken2_child
            .wait()
            .await
            .map_err(|e| anyhow!("Failed to wait on kraken2: {}", e))?;
        if !status.success() {
            return Err(anyhow!("kraken2 exited with non-zero status"));
        }
        Ok(())
    });
    cleanup_tasks.push(kraken2_wait_task);

    // Stream classified FASTQ from named pipe
    let kraken2_file = tokio::fs::File::open(&kraken2_pipe_path)
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to open named pipe: {}", e)))?;
    let parse_rx = parse_fastq(kraken2_file, config.base_buffer_size)
        .await
        .map_err(|e| PipelineError::Other(anyhow!("Failed to parse FASTQ from named pipe: {}", e)))?;

    // Filter using parse_and_filter_fastq_id
    let pattern = format!("kraken:taxid|{}", target_taxid);
    let (filtered_rx, filter_task) = parse_and_filter_fastq_id(
        parse_rx,
        config.base_buffer_size,
        pattern,
    );
    cleanup_tasks.push(filter_task);


    let (parse_output_tx, parse_output_rx) = mpsc::channel(config.base_buffer_size);
    let conversion_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(filtered_rx);
        while let Some(record) = stream.next().await {
            parse_output_tx
                .send(ParseOutput::Fastq(record))
                .await
                .map_err(|e| anyhow!("Failed to send ParseOutput: {}", e))?;
        }
        Ok(())
    });

    cleanup_tasks.push(conversion_task);



    let filtered_stream = ReceiverStream::new(parse_output_rx);

    // Before t_junction, split for check
    let (check_tx, mut check_rx) = mpsc::channel(1);
    let (forward_tx, forward_rx) = mpsc::channel(config.base_buffer_size);
    let check_split_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(filtered_stream.into_inner());
        let mut found = false;
        while let Some(item) = stream.next().await {
            if matches!(item, ParseOutput::Fastq(_)) {
                found = true;
            }
            forward_tx
                .send(item)
                .await
                .map_err(|e| anyhow!("Failed to forward item: {}", e))?;
            if found {
                break; // Exit after first Fastq record
            }
        }
        // eprintln!("check_split_task found {} Fastq record", if found { "a" } else { "no" });
        check_tx
            .send(if found { 1 } else { 0 })
            .await
            .map_err(|e| anyhow!("Failed to send check count: {}", e))?;
        // Continue forwarding remaining records
        while let Some(item) = stream.next().await {
            forward_tx
                .send(item)
                .await
                .map_err(|e| anyhow!("Failed to forward item: {}", e))?;
        }
        Ok(())
    });
    cleanup_tasks.push(check_split_task);

    let filtered_stream = ReceiverStream::new(forward_rx);

    // Split stream for output and compression
    let (kraken_streams, kraken_done_rx) = t_junction(
        filtered_stream,
        2, // Only two streams now
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::IlluminaFastq,
        "filter_reads_output".to_string(),
        None,
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    let mut cleanup_receivers = vec![kraken_done_rx];
    let mut streams_iter = kraken_streams.into_iter();
    let kraken_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let kraken_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let check_count = check_rx.recv().await.ok_or_else(|| {
        eprintln!("check_rx.recv failed");
        PipelineError::EmptyStream
    })?;


    let pigz_args = generate_cli(PIGZ_TAG, &config, None)?;
    let (mut pigz_child, pigz_stream_task, pigz_err_task) = stream_to_cmd(
        config.clone(),
        kraken_file_stream,
        PIGZ_TAG,
        pigz_args,
        StreamDataType::IlluminaFastq,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: PIGZ_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(pigz_stream_task);
    cleanup_tasks.push(pigz_err_task);

    let pigz_out_stream = {
        let mut guard = pigz_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout, ParseMode::Bytes, config.base_buffer_size)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: PIGZ_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    let pigz_write_task = tokio::spawn(stream_to_file(pigz_out_stream, final_compressed_path));
    cleanup_tasks.push(pigz_write_task);

    let cleanup_pipe_task = tokio::spawn(async move {
        tokio::fs::remove_file(&kraken2_pipe_path)
            .await
            .map_err(|e| anyhow!("Failed to remove named pipe: {}", e))?;
        Ok(())
    });
    cleanup_tasks.push(cleanup_pipe_task);

    Ok((ReceiverStream::new(kraken_output_stream), cleanup_tasks, cleanup_receivers))
}


/// Runs minimap2 to align to the target reference, then sorts.
///
/// # Arguments
///
/// * `config` - reading stream
/// * `input_stream` - stream buffer size
/// * 'targer_ref_path' -
///
/// # Returns
/// samtools_sort_out_stream: Result<(ReceiverStream<ParseOutput>, <- uncomrpessed BAM
/// cleanup_tasks
/// quast_write_tasks
///
async fn align_to_target(
    config: Arc<RunConfig>,
    input_stream: ReceiverStream<ParseOutput>,  // FASTQ SequenceRecord stream
    target_index_path: PathBuf,  // Minimap2 .mmi index
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, Vec<JoinHandle<Result<(), anyhow::Error>>>, PathBuf), PipelineError> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];
    let mut quast_write_tasks = vec![];

    
    let minimap2_args = generate_cli(MINIMAP2_TAG, &config, Some(&(target_index_path)))
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut minimap2_child, minimap2_stream_task, minimap2_err_task) = stream_to_cmd(
        config.clone(),
        input_stream.into_inner(),  // Convert to Receiver<ParseOutput> for streaming
        MINIMAP2_TAG,
        minimap2_args,
        StreamDataType::JustBytes,  // FASTQ bytes, not structured records
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MINIMAP2_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(minimap2_err_task);

    let minimap2_out_stream = {
        let mut guard = minimap2_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MINIMAP2_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_sort_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Sort,
        subcommand_fields: HashMap::from([
            ("-u".to_string(), None), // Uncompressed
            ("-O".to_string(), Some("bam".to_string())), // BAM output
            ("-".to_string(), None)
        ]),
    };
    let samtools_sort_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_sort_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut samtools_sort_child, samtools_sort_task, samtools_sort_err_task) = stream_to_cmd(
        config.clone(),
        minimap2_out_stream,
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

    let samtools_sort_out_stream = ReceiverStream::new(samtools_sort_out_stream);

    let align_bam_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("target_aligned.bam"),
        "_",
    );

    let (sam_streams, sam_done_rx) = t_junction(
        samtools_sort_out_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "align_to_target".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(sam_done_rx);
    let mut streams_iter = sam_streams.into_iter();
    let sam_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let sam_check_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    //  check stream convert to SAM for non-empty chk
    let sam_view_args = generate_cli(SAMTOOLS_TAG, &config, Some(&SamtoolsConfig {
        subcommand: SamtoolsSubcommand::View,
        subcommand_fields: HashMap::from([
            ("-h".to_string(), None),
            ("-".to_string(), None),
        ]),
    }))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut sam_view_child, sam_view_task, sam_view_err_task) = stream_to_cmd(
        config.clone(),
        sam_check_stream,
        SAMTOOLS_TAG,
        sam_view_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(sam_view_task);
    cleanup_tasks.push(sam_view_err_task);

    let sam_check_stream = {
        let mut guard = sam_view_child.lock().await;
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

    //check for alignment records (not just ehader)
    let (check_tx, mut check_rx) = mpsc::channel(1);
    let check_task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(sam_check_stream);
        let mut alignment_count = 0;
        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Bytes(line) => {
                    let line_str = String::from_utf8_lossy(&line);
                    if !line_str.starts_with('@') && !line_str.trim().is_empty() {
                        alignment_count = 1;
                        // Spawn background drain to avoid blocking
                        tokio::spawn(async move {
                            while stream.next().await.is_some() {}
                            Ok::<(), anyhow::Error>(())
                        });
                        break; // Exit early to unblock sender
                    }
                }
                _ => return Err(anyhow!("Unexpected item type in sam_check_stream")),
            }
        }
        check_tx
            .send(alignment_count)
            .await
            .map_err(|e| anyhow!("Failed to send alignment count: {}", e))?;
        Ok::<(), anyhow::Error>(())
    });
    cleanup_tasks.push(check_task);

    let alignment_count = check_rx
        .recv()
        .await
        .ok_or(PipelineError::EmptyStream)?;
    if alignment_count == 0 {
        let drain_rxs = vec![sam_output_stream, sam_file_stream];
        let drain_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = drain_rxs
            .into_iter()
            .map(|rx| {
                tokio::spawn(async move {
                    let mut stream = ReceiverStream::new(rx);
                    while stream.next().await.is_some() {}
                    Ok(())
                })
            })
            .collect();
        cleanup_tasks.extend(drain_tasks);
        return Err(PipelineError::Other(anyhow!("No alignment records in BAM stream")));
    }

    let bam_write_task = tokio::spawn(stream_to_file(sam_file_stream, align_bam_path.clone()));
    quast_write_tasks.push(bam_write_task);

    Ok((
        ReceiverStream::new(sam_output_stream),
        cleanup_tasks,
        cleanup_receivers,
        quast_write_tasks,
        align_bam_path
    ))
}

async fn generate_consensus(
    config: Arc<RunConfig>,
    bam_stream: ReceiverStream<ParseOutput>, // Uncompressed BAM byte stream
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, ReceiverStream<ParseOutput>, PathBuf, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>, Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];
    let mut quast_write_tasks = vec![];

    // Add primer trimming with samtools ampliconclip
    let primer_bed_path = config.args.primer_bed_path.clone().ok_or_else(|| PipelineError::InvalidConfig("Primer BED file not provided".to_string()))?;
    let trim_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Ampliconclip,
        subcommand_fields: HashMap::from([
            ("-".to_string(), None),
            ("-b".to_string(), Some(primer_bed_path)),
            ("--both-ends".to_string(), None),
            ("-o".to_string(), Some("-".to_string())),
        ]),
    };
    let trim_args = generate_cli(SAMTOOLS_TAG, &config, Some(&trim_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;


    let (mut trim_child, trim_task, trim_err_task) = stream_to_cmd(
        config.clone(),
        bam_stream.into_inner(),
        SAMTOOLS_TAG,
        trim_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(trim_task);
    cleanup_tasks.push(trim_err_task);

    let trimmed_bam_stream = {
        let mut guard = trim_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Bytes, config.base_buffer_size)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_consensus_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Consensus,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let samtools_consensus_args = generate_cli(SAMTOOLS_TAG, &config, Some(&samtools_consensus_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (samtools_consensus_child, samtools_consensus_task, samtools_consensus_err_task) = stream_to_cmd(
        config.clone(),
        ReceiverStream::new(trimmed_bam_stream).into_inner(),
        SAMTOOLS_TAG,
        samtools_consensus_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: SAMTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(samtools_consensus_task);
    cleanup_tasks.push(samtools_consensus_err_task);

    let samtools_consensus_out_stream = {
        let mut guard = samtools_consensus_child.lock().await;
        parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Fasta, config.base_buffer_size)
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: SAMTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let samtools_consensus_out_stream = ReceiverStream::new(samtools_consensus_out_stream);

    let (consensus_streams, consensus_done_rx) = t_junction(
        samtools_consensus_out_stream,
        3,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "generate_consensus".to_string(),
        None,
    ).await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(consensus_done_rx);
    let mut streams_iter = consensus_streams.into_iter();
    let consensus_realign_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let consensus_stats_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let consensus_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let consensus_file_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("consensus.fa"),
        "_",
    );

    let consensus_write_task = tokio::spawn(stream_to_file(
        consensus_file_stream,
        consensus_file_path.clone(),
    ));
    quast_write_tasks.push(consensus_write_task);

    Ok((
        ReceiverStream::new(consensus_realign_stream),
        ReceiverStream::new(consensus_stats_stream),
        consensus_file_path,
        cleanup_tasks,
        cleanup_receivers,
        quast_write_tasks,
    ))
}

async fn call_variants(
    config: Arc<RunConfig>,
    bam_stream: ReceiverStream<ParseOutput>,  //Uncompressed BAM byte stream
    target_ref_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(ReceiverStream<ParseOutput>, PathBuf, Vec<JoinHandle<Result<(), anyhow::Error>>>, Vec<oneshot::Receiver<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];
    let mut cleanup_receivers = vec![];

    let bcftools_mpileup_config = BcftoolsConfig {
        subcommand: BcftoolsSubcommand::Mpileup,
        subcommand_fields: HashMap::from([
            ("-f".to_string(), Some(target_ref_path.to_string_lossy().into_owned())),
            ("-Ou".to_string(), None),  // Uncompressed BCF
            ("-".to_string(), None),]),
    };
    let bcftools_mpileup_args = generate_cli(BCFTOOLS_TAG, &config, Some(&bcftools_mpileup_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bcftools_mpileup_child, bcftools_mpileup_task, bcftools_mpileup_err_task) = stream_to_cmd(
        config.clone(),
        bam_stream.into_inner(),
        BCFTOOLS_TAG,
        bcftools_mpileup_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bcftools_mpileup_task);
    cleanup_tasks.push(bcftools_mpileup_err_task);

    let bcftools_mpileup_out_stream = {
        let mut guard = bcftools_mpileup_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BCFTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };

    let bcftools_call_config = BcftoolsConfig {
        subcommand: BcftoolsSubcommand::Call,
        subcommand_fields: HashMap::from([
            ("--ploidy".to_string(), Some("1".to_string())),
            ("-m".to_string(), None),
            ("-v".to_string(), None),
            ("-P".to_string(), Some(config.args.bcftools_call_theta.to_string())),
            ("-Ou".to_string(), None),  // Uncompressed BCF
            ("-".to_string(), None),
        ])
    };
    let bcftools_call_args = generate_cli(BCFTOOLS_TAG, &config, Some(&bcftools_call_config))
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;

    let (mut bcftools_call_child, bcftools_call_task, bcftools_call_err_task) = stream_to_cmd(
        config.clone(),
        bcftools_mpileup_out_stream,
        BCFTOOLS_TAG,
        bcftools_call_args,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: BCFTOOLS_TAG.to_string(),
            error: e.to_string(),
        })?;
    cleanup_tasks.push(bcftools_call_task);
    cleanup_tasks.push(bcftools_call_err_task);

    // Parse the bcftools call output
    let bcftools_call_out_stream = {
        let mut guard = bcftools_call_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        )
            .await
            .map_err(|e| PipelineError::ToolExecution {
                tool: BCFTOOLS_TAG.to_string(),
                error: e.to_string(),
            })?
    };
    let bcftools_call_out_stream = ReceiverStream::new(bcftools_call_out_stream);

    let (vcf_streams, vcf_done_rx) = t_junction(
        bcftools_call_out_stream,
        2,
        config.base_buffer_size,
        config.args.stall_threshold,
        None,
        100,
        StreamDataType::JustBytes,
        "call_variants".to_string(),
        None
    )
        .await
        .map_err(|_| PipelineError::StreamDataDropped)?;

    cleanup_receivers.push(vcf_done_rx);

    let mut streams_iter = vcf_streams.into_iter();
    let vcf_output_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;
    let vcf_file_stream = streams_iter.next().ok_or(PipelineError::EmptyStream)?;

    let vcf_file_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("variants.bcf"),
        "_",
    );

    let vcf_write_task = tokio::spawn(stream_to_file(
        vcf_file_stream,
        vcf_file_path.clone(),
    ));
    cleanup_tasks.push(vcf_write_task);

    Ok((ReceiverStream::new(vcf_output_stream), vcf_file_path, cleanup_tasks, cleanup_receivers))
}


async fn realign_consensus_to_ref(
    config: Arc<RunConfig>,
    consensus_realign_stream: Receiver<ParseOutput>,  // FASTA
    target_ref_fasta_path: PathBuf,
    out_dir: &PathBuf,
    no_ext_sample_base_buf: &PathBuf,
) -> Result<(Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let reference_file = TokioFile::open(&target_ref_fasta_path).await
        .map_err(|e| PipelineError::Other(e.into()))?;
    let reference_rx = parse_bytes(reference_file, config.base_buffer_size).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let realign_streams = vec![consensus_realign_stream, reference_rx];
    let (combined_rx, combined_task) = y_junction(realign_streams, config.base_buffer_size).await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let realign_consensus_path = file_path_manipulator(
        no_ext_sample_base_buf,
        Some(out_dir),
        None,
        Some("consensus_realigned.fa"),
        "_"
    );
    let realign_mafft_args = generate_cli(MAFFT_TAG, &config, None)
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;
    let (
        mut realign_consensus_mafft_child,
        realign_consensus_mafft_task,
        realign_consensus_mafft_err_task
    ) = stream_to_cmd(config.clone(),
                      combined_rx,
                      MAFFT_TAG,
                      realign_mafft_args,
                      StreamDataType::JustBytes,
                      config.args.verbose
    ).await
        .map_err(|e| PipelineError::ToolExecution {
            tool: MAFFT_TAG.to_string(),
            error: e.to_string(),
        })?;

    let mut cleanup_tasks = vec![combined_task, realign_consensus_mafft_task, realign_consensus_mafft_err_task];

    let consensus_samtools_out_stream = {
        let mut guard = realign_consensus_mafft_child.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Bytes,
            config.base_buffer_size,
        ).await
            .map_err(|e| PipelineError::ToolExecution {
                tool: MAFFT_TAG.to_string(),
                error: e.to_string(),
            })?
    };


    let realign_consensus_write_task = tokio::spawn(stream_to_file(
        consensus_samtools_out_stream,
        realign_consensus_path,
    ));
    cleanup_tasks.push(realign_consensus_write_task);



    Ok(cleanup_tasks)
}


async fn calculate_statistics(
    config: Arc<RunConfig>,
    no_ext_sample_base: &str,
    consensus_bam_stats_stream: Option<ReceiverStream<ParseOutput>>,  //Uncompressed BAM byte stream
    consensus_bam_depth_stream: Option<ReceiverStream<ParseOutput>>,  //Uncompressed BAM byte stream
    no_host_seqkit_out_stream_stats: Receiver<ParseOutput>,
    ercc_stats_task: Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>>,
    consensus_stats_stream: Option<ReceiverStream<ParseOutput>>,
    call_bcftools_stats_stream: Option<ReceiverStream<ParseOutput>>, // Uncompressed BCF stream (-Ou)
    out_dir: &PathBuf,
    technology: Technology,
) -> Result<(), anyhow::Error> {
    let mut local_cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();


    let stats_samtools_config_stats = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Stats,
        subcommand_fields: HashMap::from([("-".to_string(), None)]),
    };
    let stats_samtools_args_stats = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&stats_samtools_config_stats),
    )?;

    let stats_samtools_out_stream_stats = match consensus_bam_stats_stream {
        Some(stream) => {
            let (mut stats_samtools_child_stats, stats_samtools_task_stats, stats_samtools_err_task_stats) = stream_to_cmd(config.clone(),
                                                                                                                           stream.into_inner(),
                                                                                                                           SAMTOOLS_TAG,
                                                                                                                           stats_samtools_args_stats.clone(),
                                                                                                                           StreamDataType::JustBytes,
                                                                                                                           config.args.verbose,
            ).await?;
            local_cleanup_tasks.push(stats_samtools_task_stats);
            local_cleanup_tasks.push(stats_samtools_err_task_stats);
            {
                let mut guard = stats_samtools_child_stats.lock().await;
                parse_child_output(&mut guard, ChildStream::Stdout, ParseMode::Lines, config.base_buffer_size / 2).await?
            }
        }
        None => return Err(anyhow!("consensus_bam_stats_stream is not available")),
    };
    let samtools_stats_out = parse_samtools_stats(stats_samtools_out_stream_stats).await?;

    let depth_samtools_config = SamtoolsConfig {
        subcommand: SamtoolsSubcommand::Depth,
        subcommand_fields: HashMap::from([
            ("-aa".to_string(), None),
            ("-d".to_string(), Some("0".to_string())),
            ("-".to_string(), None),
        ]),
    };
    let depth_samtools_args = generate_cli(
        SAMTOOLS_TAG,
        &config,
        Some(&depth_samtools_config),
    )?;

    let depth_samtools_out_stream = match consensus_bam_depth_stream {
        Some(stream) => {
            let (mut depth_samtools_child, depth_samtools_task, depth_samtools_err_task) = stream_to_cmd(config.clone(),
                                                                                                         stream.into_inner(),
                                                                                                         SAMTOOLS_TAG,
                                                                                                         depth_samtools_args,
                                                                                                         StreamDataType::JustBytes,
                                                                                                         config.args.verbose,
            ).await?;
            local_cleanup_tasks.push(depth_samtools_task);
            local_cleanup_tasks.push(depth_samtools_err_task);
            {
                let mut guard = depth_samtools_child.lock().await;
                parse_child_output(
                    &mut guard,
                    ChildStream::Stdout,
                    ParseMode::Lines,
                    config.base_buffer_size / 2,
                ).await?
            }
        }
        None => {
            return Err(anyhow!("consensus_bam_depth_stream is not available"));
        }
    };

    let depth_map = parse_samtools_depth(depth_samtools_out_stream).await?;
    if depth_map.is_empty() {
        return Err(anyhow!("No depth data found"));
    }
    let first_chr = depth_map.keys().next().ok_or_else(|| anyhow!("No chromosomes found"))?.clone();
    let first_chr_depth_map = depth_map.get(&first_chr).unwrap();
    let depths: Vec<u32> = first_chr_depth_map.values().copied().collect();
    let samtools_depth_stats = compute_depth_stats(&depths)?;

    let depth_plot_path = file_path_manipulator(
        &PathBuf::from(no_ext_sample_base),
        Some(out_dir),
        None,
        Some("depth.png"),
        "_"
    );

    plot_depths(&first_chr_depth_map, no_ext_sample_base, &depth_plot_path)?;

    let seqkit_stats = parse_seqkit_stats(no_host_seqkit_out_stream_stats).await?;

    let ercc_stats = if let Technology::Illumina = technology {
        match ercc_stats_task {
            Some(task) => match task.await? {
                Ok(stats) => stats,
                Err(e) => return Err(anyhow!("Failed to parse ERCC stats: {}", e)),
            },
            None => return Err(anyhow!("ERCC stats task not initialized for Illumina technology")),
        }
    } else {
        HashMap::new()
    };

    let allele_counts = if let Some(stream) = consensus_stats_stream {
        compute_allele_counts(stream.into_inner()).await?
    } else {
        HashMap::new()
    };


    let (ref_snps, ref_mnps, ref_indels) = if let Some(stream) = call_bcftools_stats_stream {
        let bcftools_stats_args = vec!["stats".to_string(), "-".to_string()];
        let (mut bcftools_stats_child, bcftools_stats_stream_task, bcftools_stats_err_task) = stream_to_cmd(
            config.clone(),
            stream.into_inner(),
            BCFTOOLS_TAG,
            bcftools_stats_args,
            StreamDataType::JustBytes,
            config.args.verbose,
        ).await?;

        local_cleanup_tasks.push(bcftools_stats_stream_task);
        local_cleanup_tasks.push(bcftools_stats_err_task);

        // Parse the stdout of bcftools stats as lines
        let bcftools_stats_out_rx = {
            let mut guard = bcftools_stats_child.lock().await;
            parse_child_output(
                &mut guard,
                ChildStream::Stdout,
                ParseMode::Lines,
                config.base_buffer_size,
            ).await?
        };

        count_variants_from_bcftools_stats(bcftools_stats_out_rx).await?
    } else {
        (0, 0, 0)
    };



    let n_actg = allele_counts.iter().filter(|&(k, _)| "ACTGU".contains(*k)).map(|(_, &v)| v).sum::<u64>();
    let n_missing = allele_counts.get(&'N').copied().unwrap_or(0);
    let n_gap = allele_counts.get(&'-').copied().unwrap_or(0);
    let n_ambiguous = allele_counts.iter().filter(|&(k, _)| !"ACTGUN-".contains(*k)).map(|(_, &v)| v).sum::<u64>();

    let coverage_breadth = if !depths.is_empty() { depths.iter().filter(|&&d| d > 0).count() as f64 / depths.len() as f64 } else { 0.0 };
    let max_aligned_length = depths.len();
    let total_length = depths.len();
    let (coverage_bin_size, coverage) = if !depths.is_empty() { compute_coverage_bins(&depths, 500) } else { (0.0, Vec::new()) };

    let stats = Stats {
        sample_name: no_ext_sample_base.to_string(),
        depth_avg: samtools_depth_stats.get("depth_avg").copied().unwrap_or(0.0),
        depth_q25: samtools_depth_stats.get("depth_q.25").copied().unwrap_or(0.0),
        depth_q50: samtools_depth_stats.get("depth_q.5").copied().unwrap_or(0.0),
        depth_q75: samtools_depth_stats.get("depth_q.75").copied().unwrap_or(0.0),
        depth_frac_above_10x: samtools_depth_stats.get("depth_frac_above_10x").copied().unwrap_or(0.0),
        depth_frac_above_25x: samtools_depth_stats.get("depth_frac_above_25x").copied().unwrap_or(0.0),
        depth_frac_above_50x: samtools_depth_stats.get("depth_frac_above_50x").copied().unwrap_or(0.0),
        depth_frac_above_100x: samtools_depth_stats.get("depth_frac_above_100x").copied().unwrap_or(0.0),
        allele_counts,
        total_reads: seqkit_stats.get("num_seqs").and_then(|s| s.parse::<u64>().ok()).unwrap_or(0),
        mapped_reads: samtools_stats_out.get("reads mapped").and_then(|s| s.parse::<u64>().ok()).unwrap_or(0),
        mapped_paired: samtools_stats_out.get("reads mapped and paired").and_then(|s| s.parse::<u64>().ok()),
        paired_inward: samtools_stats_out.get("inward oriented pairs").and_then(|s| s.parse::<u64>().ok()).map(|v| v * 2),
        paired_outward: samtools_stats_out.get("outward oriented pairs").and_then(|s| s.parse::<u64>().ok()).map(|v| v * 2),
        paired_other_orientation: samtools_stats_out.get("pairs with other orientation").and_then(|s| s.parse::<u64>().ok()).map(|v| v * 2),
        ercc_mapped_reads: ercc_stats.get("ercc_mapped_reads").copied(),
        ercc_mapped_paired: ercc_stats.get("ercc_mapped_paired").copied(),
        ref_snps,
        ref_mnps,
        ref_indels,
        n_actg,
        n_missing,
        n_gap,
        n_ambiguous,
        coverage_breadth,
        max_aligned_length,
        total_length,
        coverage_bin_size,
        coverage,
    };

    let stats_file_path = out_dir.join(format!("{}_stats.json", no_ext_sample_base));
    let mut stats_file = File::create(&stats_file_path)?;
    serde_json::to_writer_pretty(&mut stats_file, &stats)?;

    try_join_all(local_cleanup_tasks).await?.into_iter().collect::<Result<Vec<_>, _>>()?;

    Ok(())
}

async fn evaluate_assembly(
    config: Arc<RunConfig>,
    target_ref_fasta_path: PathBuf,
    align_bam_path: PathBuf,
    consensus_file_path: PathBuf,
) -> Result<(Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    let mut cleanup_tasks = vec![];
    let quast_config = QuastConfig {
        ref_fasta: target_ref_fasta_path.to_string_lossy().into_owned(),
        ref_bam: align_bam_path.to_string_lossy().into_owned(),
        assembly_fasta: consensus_file_path.to_string_lossy().into_owned(),
    };

    let assembly_eval_quast_args = generate_cli(
        QUAST_TAG,
        &*config,
        Some(&quast_config),
    ).map_err(|e| PipelineError::ToolExecution {
        tool: QUAST_TAG.to_string(),
        error: e.to_string(),
    })?;

    let (mut assembly_eval_quast_child, assembly_eval_quast_err_task) = spawn_cmd(
        config.clone(),
        QUAST_TAG,
        assembly_eval_quast_args,
        config.args.verbose,
    ).await.map_err(|e| PipelineError::ToolExecution {
        tool: QUAST_TAG.to_string(),
        error: e.to_string(),
    })?;
    cleanup_tasks.push(assembly_eval_quast_err_task);

    // Await child exit in a task to keep it alive
    let quast_wait_task = tokio::spawn(async move {
        let status = assembly_eval_quast_child.wait().await?;
        if !status.success() {
            return Err(anyhow!("QUAST exited with non-zero status: {}", status));
        }
        Ok(())
    });
    cleanup_tasks.push(quast_wait_task);


    Ok(cleanup_tasks)
}





pub async fn run(config: Arc<RunConfig>) -> Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let ram_temp_dir = config.ram_temp_dir.clone();
    let out_dir = config.out_dir.clone();
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let mut cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut cleanup_receivers: Vec<oneshot::Receiver<Result<(), anyhow::Error>>> = Vec::new();
    let mut quast_write_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut stats_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();
    let mut ercc_stats_task: Option<JoinHandle<Result<HashMap<String, u64>, anyhow::Error>>> = None;

    let mut target_ref_fasta_path: Option<PathBuf> = None;
    let mut consensus_file_path: Option<PathBuf> = None;
    let mut align_bam_path: Option<PathBuf> = None;
    let mut align_sam_stats_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut align_sam_depth_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut consensus_stats_stream: Option<ReceiverStream<ParseOutput>> = None;
    let mut call_bcftools_stats_stream: Option<ReceiverStream<ParseOutput>> = None;

    // External tools check
    check_versions(vec![
        SAMTOOLS_TAG,
        MINIMAP2_TAG,
        FASTP_TAG,
        SAMTOOLS_TAG,
        KRAKEN2_TAG,
        BCFTOOLS_TAG,
        MAFFT_TAG,
        SEQKIT_TAG,
        QUAST_TAG,
    ])
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    let (file1_path, file2_path, no_ext_sample_base_buf, no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;

    let technology = config.args.technology.clone();

    let ref_db_path: Option<PathBuf> = config.args.ref_db.as_ref().map(PathBuf::from);


    // Retrieve Index
    let index_start = Instant::now();
    let h5_index = get_index(&config.args)
        .await
        .map_err(|e| PipelineError::ReferenceRetrievalFailed(e.to_string()))?;
    if h5_index.is_some() {
        println!("Index retrieve time: {} milliseconds.", index_start.elapsed().as_millis());
    }

    let (host_ref_fasta_path, host_ref_index_path, host_ref_temp, host_index_temp, host_ref_tasks) = prepare_reference_and_index(
        &config,
        ref_db_path.clone(),
        &ram_temp_dir,
        h5_index.as_ref(),
        config.args.host_accession.clone(),
        config.args.host_sequence.clone(),
        config.args.host_index.clone(),
        "host",
    )
        .await?;
    cleanup_tasks.extend(host_ref_tasks);
    if let Some(temp) = host_ref_temp {
        temp_files.push(temp);
    }
    if let Some(index_temp) = host_index_temp {
        temp_files.push(index_temp);
    }

    let (target_ref_fasta_path_inner, target_ref_index_path, target_ref_temp, target_index_temp, target_ref_tasks) = prepare_reference_and_index(
        &config,
        ref_db_path.clone(),
        &ram_temp_dir,
        h5_index.as_ref(),
        config.args.target_accession.clone(),
        config.args.target_sequence.clone(),
        config.args.target_index.clone(),
        "target",
    )
        .await?;
    cleanup_tasks.extend(target_ref_tasks);
    if let Some(temp) = target_ref_temp {
        temp_files.push(temp);
    }
    if let Some(index_temp) = target_index_temp {
        temp_files.push(index_temp);
    }
    let target_fasta = target_ref_fasta_path_inner.ok_or_else(|| {
        PipelineError::InvalidConfig("Target FASTA required for filter_with_kraken, call_variants, and evaluate_assembly".to_string())
    })?;
    target_ref_fasta_path = Some(target_fasta.clone());


    // Input Validation
    let (val_fastp_out_stream, validate_cleanup_tasks, validate_cleanup_receivers) = validate_input(
        config.clone(),
        file1_path,
        file2_path,
        no_ext_sample_base_buf.clone(),
        &out_dir,
    )
        .await?;
    cleanup_tasks.extend(validate_cleanup_tasks);
    cleanup_receivers.extend(validate_cleanup_receivers);


    // Host Removal
    let no_host_file_path = file_path_manipulator(
        &no_ext_sample_base_buf,
        Some(&out_dir),
        None,
        Some("no_host.fq.gz"),
        "_",
    );

    let (no_host_output_stream, no_host_seqkit_out_stream_stats, no_host_cleanup_tasks, no_host_cleanup_receivers) = align_to_host(
        config.clone(),
        val_fastp_out_stream,
        host_ref_index_path,
        no_host_file_path,
    )
        .await?;
    cleanup_tasks.extend(no_host_cleanup_tasks);
    cleanup_receivers.extend(no_host_cleanup_receivers);

    // Counting stats for the host-removed reads
    let stats_seqkit_config_stats = SeqkitConfig {
        subcommand: SeqkitSubcommand::Stats,
        subcommand_fields: HashMap::from([]),
    };
    let stats_seqkit_args_stats = generate_cli(SEQKIT_TAG, &config, Some(&stats_seqkit_config_stats))?;
    let no_host_seqkit_out_stream_stats = no_host_seqkit_out_stream_stats.into_inner();
    let (mut no_host_seqkit_child_stats, no_host_seqkit_task_stats, no_host_seqkit_err_task_stats) = stream_to_cmd(
        config.clone(),
        no_host_seqkit_out_stream_stats,
        SEQKIT_TAG,
        stats_seqkit_args_stats,
        StreamDataType::JustBytes,
        config.args.verbose,
    )
        .await?;
    let no_host_seqkit_out_stream_stats = {
        let mut guard = no_host_seqkit_child_stats.lock().await;
        parse_child_output(
            &mut guard,
            ChildStream::Stdout,
            ParseMode::Lines,
            config.base_buffer_size / 2,
        )
            .await?
    };
    stats_tasks.push(no_host_seqkit_task_stats);
    cleanup_tasks.push(no_host_seqkit_err_task_stats);

    // Split by Technology
    match technology {
        Technology::Illumina => {
            eprintln!("Technology: Illumina");
            let (ercc_fasta_path, ercc_index_path, ercc_ref_temp, ercc_index_temp, ercc_ref_tasks) = prepare_reference_and_index(
                &config,
                None, // No HDF5 for ERCC
                &ram_temp_dir,
                None, // No HDF5 index
                None, // No accession
                config.args.ercc_sequence.clone(), // Use ercc_sequences as sequence path
                config.args.ercc_index.clone(),
                "ercc",
            )
                .await?;
            cleanup_tasks.extend(ercc_ref_tasks);
            if let Some(temp) = ercc_ref_temp {
                temp_files.push(temp);
            }
            if let Some(index_temp) = ercc_index_temp {
                temp_files.push(index_temp);
            }

            // ERCC
            let (no_host_ercc_stream, ercc_stats_out_task, mut ercc_cleanup_tasks, mut ercc_cleanup_receivers) = process_ercc(
                config.clone(),
                no_host_output_stream,
                ercc_index_path,
                &out_dir,
                &no_ext_sample_base,
            )
                .await?;
            ercc_stats_task = ercc_stats_out_task;
            cleanup_tasks.append(&mut ercc_cleanup_tasks);
            cleanup_receivers.append(&mut ercc_cleanup_receivers);

            // Filter Reads
            let (filter_reads_out_stream, filter_reads_cleanup_tasks, filter_reads_cleanup_receivers) = filter_with_kraken(
                config.clone(),
                no_host_ercc_stream,
                target_fasta.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
                config.args.target_taxid.as_ref().expect("target_taxid must be set"),
            )
                .await?;
            cleanup_tasks.extend(filter_reads_cleanup_tasks);
            cleanup_receivers.extend(filter_reads_cleanup_receivers);

            // Align Reads to Target
            let (sam_output_stream, align_cleanup_tasks, align_cleanup_receivers, align_quast_tasks, align_bam_path_inner) = align_to_target(
                config.clone(),
                filter_reads_out_stream,
                target_ref_index_path,
                &out_dir,
                &no_ext_sample_base_buf,
            )
                .await?;
            cleanup_tasks.extend(align_cleanup_tasks);
            cleanup_receivers.extend(align_cleanup_receivers);
            quast_write_tasks.extend(align_quast_tasks);
            align_bam_path = Some(align_bam_path_inner);

            // Split SAM streams for bypass, stats, etc
            let (align_sam_streams, align_sam_done_rx) = t_junction(
                sam_output_stream,
                4,
                config.base_buffer_size,
                config.args.stall_threshold,
                None,
                100,
                StreamDataType::JustBytes,
                "pipeline_aligned_sam_split".to_string(),
                None,
            )
                .await
                .map_err(|_| PipelineError::StreamDataDropped)?;
            cleanup_receivers.push(align_sam_done_rx);
            let mut align_streams_iter = align_sam_streams.into_iter();
            let align_sam_output_stream = align_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
            let align_sam_call_stream = align_streams_iter.next().ok_or(PipelineError::EmptyStream)?;
            align_sam_stats_stream = Some(ReceiverStream::new(align_streams_iter.next().ok_or(PipelineError::EmptyStream)?));
            align_sam_depth_stream = Some(ReceiverStream::new(align_streams_iter.next().ok_or(PipelineError::EmptyStream)?));

            let align_sam_output_stream = ReceiverStream::new(align_sam_output_stream);
            let align_sam_call_stream = ReceiverStream::new(align_sam_call_stream);

            // Make Consensus
            let (consensus_realign_stream, consensus_stats_stream_rx, consensus_file_path_x, consensus_cleanup_tasks, consensus_cleanup_receivers, consensus_quast_tasks) =
                generate_consensus(config.clone(), align_sam_output_stream, &out_dir, &no_ext_sample_base_buf).await?;
            consensus_file_path = Some(consensus_file_path_x);
            cleanup_tasks.extend(consensus_cleanup_tasks);
            cleanup_receivers.extend(consensus_cleanup_receivers);
            quast_write_tasks.extend(consensus_quast_tasks);
            consensus_stats_stream = Some(consensus_stats_stream_rx);

            // Call Variants
            let (call_bcftools_stats_stream_out, _, call_cleanup_tasks, call_cleanup_receivers) =
                call_variants(config.clone(), align_sam_call_stream, target_fasta.clone(), &out_dir, &no_ext_sample_base_buf).await?;
            cleanup_tasks.extend(call_cleanup_tasks);
            cleanup_receivers.extend(call_cleanup_receivers);
            call_bcftools_stats_stream = Some(call_bcftools_stats_stream_out);

            // Realign Consensus to Ref
            let realign_cleanup_tasks = realign_consensus_to_ref(
                config.clone(),
                consensus_realign_stream.into_inner(),
                target_fasta.clone(),
                &out_dir,
                &no_ext_sample_base_buf,
            )
                .await?;
            cleanup_tasks.extend(realign_cleanup_tasks);
        }
        Technology::ONT => {
            eprintln!("Technology: ONT not ready");
        }
    }

    // Join stats tasks
    let results = try_join_all(stats_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }

    // Calculate Statistics
    calculate_statistics(
        config.clone(),
        &no_ext_sample_base,
        align_sam_stats_stream,
        align_sam_depth_stream,
        no_host_seqkit_out_stream_stats,
        ercc_stats_task,
        consensus_stats_stream,
        call_bcftools_stats_stream,
        &out_dir,
        technology,
    )
        .await?;

    // Assembly Evaluation
    let results = try_join_all(quast_write_tasks)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;
    for result in results {
        result.map_err(|e| PipelineError::Other(e))?;
    }

    if let (Some(target_ref_fasta_path), Some(align_bam_path), Some(consensus_file_path)) = (target_ref_fasta_path, align_bam_path, consensus_file_path) {
        evaluate_assembly(config, target_ref_fasta_path, align_bam_path, consensus_file_path).await?;
    }


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
    drop(temp_files);

    println!("Finished generating consensus genome");
    Ok(())
}