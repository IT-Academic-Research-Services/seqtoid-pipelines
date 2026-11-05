//! AMR Pipeline

use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::anyhow;
use log::{error, info, warn};
use tempfile::TempDir;
use bytes::Bytes;
use dashmap::DashMap;
use tokio::sync::{mpsc, oneshot};
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;

use crate::config::defs::{ClusterInfo, PairingMode, PipelineError, RunConfig, StreamDataType};
use crate::pipelines::short_read_mngs::{
    kallisto_results,
    run_host_filter,
    HostFilterCompletion,
};
use crate::utils::fastx::{parse_byte_stream_to_fastq, read_fastq, write_fasta_stream_to_file, SequenceRecord};
use crate::utils::file::{resolve_existing_input_path, validate_file_inputs, write_count_file};
use crate::utils::streams::{join_with_error_handling, ParseOutput};
use crate::utils::subsample::load_duplicate_clusters_csv;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AmrInputMode {
    Raw,
    Filtered,
}

impl AmrInputMode {
    fn from_config(config: &RunConfig) -> Result<Self, PipelineError> {
        match config.args.amr_input_mode.as_str() {
            "raw" => Ok(Self::Raw),
            "filtered" => Ok(Self::Filtered),
            other => Err(PipelineError::InvalidConfig(format!(
                "Invalid AMR input mode '{}'; expected 'raw' or 'filtered'",
                other
            ))),
        }
    }
}

struct RunRedupHandle {
    /// One live FASTA stream per input end, in WDL `redups_1.fa`,
    /// `redups_2.fa` order.
    redups_streams: Vec<ReceiverStream<ParseOutput>>,

    /// Completes after every per-end producer has reached EOF.
    /// Callers must install consumers for `redups_streams` before awaiting it.
    completion_task: JoinHandle<Result<(), anyhow::Error>>,
}

/// Native streaming equivalent of the WDL RunRedup task.
///
/// WDL semantics:
///   1. `subsampled_reads` form the prefix of each `redups_N.fa`.
///   2. Only a duplicate cluster whose representative survived subsampling is
///      eligible for re-expansion.
///   3. The representative itself is already present in the subsampled FASTA;
///      append only its non-representative cluster members.
///   4. Appended records are selected from the corresponding host-filtered,
///      pre-dedup FASTQ and converted to FASTA.
///
/// Rust keeps the cluster relationship in memory, so the WDL-only
/// `cluster_sizes` / `clusters.csv` parsing and temporary ID files are not
/// reproduced here.
async fn run_redup(
    config: Arc<RunConfig>,
    host_filtered_streams: Vec<ReceiverStream<ParseOutput>>,
    subsampled_streams: Vec<ReceiverStream<ParseOutput>>,
    duplicate_clusters_task: JoinHandle<Result<Arc<DashMap<String, ClusterInfo>>, anyhow::Error>>,
) -> Result<RunRedupHandle, PipelineError> {
    if host_filtered_streams.is_empty() || host_filtered_streams.len() > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRedup expected one or two host-filtered streams, got {}",
            host_filtered_streams.len()
        )));
    }

    if host_filtered_streams.len() != subsampled_streams.len() {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRedup input cardinality mismatch: host_filtered_streams={}, subsampled_streams={}",
            host_filtered_streams.len(),
            subsampled_streams.len()
        )));
    }

    let part_count = subsampled_streams.len();
    let mut output_streams = Vec::with_capacity(part_count);
    let mut prefix_tasks = Vec::with_capacity(part_count);
    let mut host_buffer_tasks = Vec::with_capacity(part_count);

    for (part_index, (mut host_stream, mut subsampled_stream)) in host_filtered_streams
        .into_iter()
        .zip(subsampled_streams.into_iter())
        .enumerate()
    {
        let (tx, rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);
        output_streams.push(ReceiverStream::new(rx));

        // Drain the pre-dedup stream immediately. It is upstream of dedup and
        // therefore cannot be left unread while we wait for subsampling or the
        // final cluster map. Preserve stream order so the later duplicate suffix
        // matches `seqkit grep ... | seqkit fq2fa`.
        let host_buffer_task = tokio::spawn(async move {
            let mut records = Vec::<SequenceRecord>::new();

            while let Some(item) = host_stream.next().await {
                let ParseOutput::Fastq(record) = item else {
                    return Err(anyhow!(
                        "RunRedup expected FASTQ on host-filtered stream for part {}",
                        part_index + 1
                    ));
                };

                records.push(record);
            }

            Ok::<Vec<SequenceRecord>, anyhow::Error>(records)
        });
        host_buffer_tasks.push(host_buffer_task);

        // The WDL output begins with the complete subsampled input. Forward it
        // immediately as FASTA while recording which representatives survived
        // subsampling. The sender is returned so the coordinator can append the
        // duplicate suffix after the cluster map becomes available.
        let prefix_task = tokio::spawn(async move {
            let mut passed_representatives = HashSet::<String>::new();

            while let Some(item) = subsampled_stream.next().await {
                let ParseOutput::Fastq(record) = item else {
                    return Err(anyhow!(
                        "RunRedup expected FASTQ on subsampled stream for part {}",
                        part_index + 1
                    ));
                };

                passed_representatives.insert(record.id().to_string());

                let fasta = SequenceRecord::Fasta {
                    id: record.id().to_string(),
                    desc: None,
                    seq: record.seq(),
                };

                tx.send(ParseOutput::Fasta(fasta))
                    .await
                    .map_err(|_| anyhow!(
                        "RunRedup downstream receiver dropped for part {}",
                        part_index + 1
                    ))?;
            }

            Ok::<_, anyhow::Error>((passed_representatives, tx))
        });
        prefix_tasks.push(prefix_task);
    }

    let completion_task = tokio::spawn(async move {
        // Every live input already has a consumer installed above, so it is now
        // safe to await dedup finalization.
        let duplicate_clusters = duplicate_clusters_task
            .await
            .map_err(|e| anyhow!("RunRedup duplicate-cluster task panicked: {}", e))??;

        for (part_index, (prefix_task, host_buffer_task)) in prefix_tasks
            .into_iter()
            .zip(host_buffer_tasks.into_iter())
            .enumerate()
        {
            let (passed_representatives, tx) = prefix_task
                .await
                .map_err(|e| anyhow!(
                    "RunRedup prefix task for part {} panicked: {}",
                    part_index + 1,
                    e
                ))??;

            let host_records = host_buffer_task
                .await
                .map_err(|e| anyhow!(
                    "RunRedup host-buffer task for part {} panicked: {}",
                    part_index + 1,
                    e
                ))??;

            let mut duplicate_member_ids = HashSet::<String>::new();

            for entry in duplicate_clusters.iter() {
                let rep_id = entry.key();
                let info = entry.value();

                if info.size <= 1 || !passed_representatives.contains(rep_id.as_str()) {
                    continue;
                }

                for member_id in &info.members {
                    if member_id != rep_id {
                        duplicate_member_ids.insert(member_id.clone());
                    }
                }
            }

            // Append selected pre-dedup records in their original stream order.
            for record in host_records {
                if !duplicate_member_ids.contains(record.id()) {
                    continue;
                }

                let fasta = SequenceRecord::Fasta {
                    id: record.id().to_string(),
                    desc: None,
                    seq: record.seq(),
                };

                tx.send(ParseOutput::Fasta(fasta))
                    .await
                    .map_err(|_| anyhow!(
                        "RunRedup downstream receiver dropped for part {}",
                        part_index + 1
                    ))?;
            }
        }

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunRedupHandle {
        redups_streams: output_streams,
        completion_task,
    })
}




async fn fastq_file_to_record_stream(
    path: PathBuf,
    config: Arc<RunConfig>,
    tag: &str,
) -> Result<(ReceiverStream<ParseOutput>, JoinHandle<Result<(), anyhow::Error>>), PipelineError> {
    let (byte_rx, read_task) = read_fastq(
        path,
        None,
        PairingMode::Strict,
        None,
        u64::MAX,
        None,
        None,
        tag,
        &config,
    )
        .map_err(PipelineError::Other)?;

    let (record_rx, parse_task) = parse_byte_stream_to_fastq(
        byte_rx,
        config.base_buffer_size,
        config.args.stall_threshold,
    )
        .await
        .map_err(PipelineError::Other)?;

    let completion_task = tokio::spawn(async move {
        read_task
            .await
            .map_err(|e| anyhow!("FASTQ file reader task panicked: {}", e))??;
        parse_task
            .await
            .map_err(|e| anyhow!("FASTQ byte parser task panicked: {}", e))??;
        Ok(())
    });

    Ok((ReceiverStream::new(record_rx), completion_task))
}

pub async fn run(
    config: Arc<RunConfig>,
) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir()
        .map_err(|e| PipelineError::Other(e.into()))?;

    let out_dir = config.out_dir.clone();

    info!("Starting AMR pipeline.");

    match AmrInputMode::from_config(&config)? {
        AmrInputMode::Raw => {
            info!(
                "[AMR] raw input selected; running shared host filtering"
            );

            // ------------------------------------------------------------
            // Raw-input setup
            // ------------------------------------------------------------

            let (
                file1_path,
                file2_path,
                sample_base_buf,
                _sample_base,
            ) = validate_file_inputs(
                &config,
                &cwd,
                true,
            )
                .await?;

            let paired = file2_path.is_some();

            let mut final_temp_dirs: Vec<TempDir> = Vec::new();

            let mut cleanup_tasks:
                Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> =
                Vec::new();

            let mut cleanup_receivers:
                Vec<
                    oneshot::Receiver<
                        anyhow::Result<(), anyhow::Error>
                    >
                > =
                Vec::new();

            // ------------------------------------------------------------
            // Shared host filtering
            //
            // AMR retains the exact host-filtered stream immediately
            // before deduplication for later AMR processing.
            // ------------------------------------------------------------

            let host_filter = run_host_filter(
                config.clone(),
                &cwd,
                &out_dir,
                &file1_path,
                file2_path.clone(),
                &sample_base_buf,
                paired,
                true,
                &mut final_temp_dirs,
                &mut cleanup_tasks,
                &mut cleanup_receivers,
            )
                .await?;

            let HostFilterCompletion {
                mut outputs,

                val_count_task,
                qc_count_result_rx,
                ercc_count_rx,

                ercc_bt2_bam_write_handle,
                host_bt2_bam_write_handle,
                optional_human_bam_write_handle,

                host_bt2_count_rx,

                ercc_bt2_bam_path: _,
                host_hisat2_count_rx,
                host_bt2_bam_path: _,

                kallisto_exit_task,
                kallisto_ercc_tx,
                kallisto_ercc_rx,

                dedup_count_rx,
                subsample_count_rx,
            } = host_filter;

            // ------------------------------------------------------------
            // Acknowledge AMR preprocessing outputs.
            // ------------------------------------------------------------
            if outputs.pre_dedup_r1.is_none() {
                return Err(PipelineError::Other(anyhow!(
        "AMR requested pre-dedup stream but R1 was not returned"
    )));
            }

            if paired && outputs.pre_dedup_r2.is_none() {
                return Err(PipelineError::Other(anyhow!(
        "AMR paired-end run requested pre-dedup stream but R2 was not returned"
    )));
            }

            info!(
    "[AMR] host-filter streaming outputs ready: \
     pre_dedup_r1=true, pre_dedup_r2={}, \
     non_host_r1=true, non_host_r2={}, \
     duplicate_clusters=pending",
    outputs.pre_dedup_r2.is_some(),
    outputs.non_host_r2.is_some(),
);

            // Install RunRedup consumers immediately so every bounded host-filter
            // output has a downstream reader before we await host-filter completion.
            let mut pre_dedup_streams = vec![
                outputs.pre_dedup_r1.take().ok_or_else(|| {
                    PipelineError::Other(anyhow!("AMR raw RunRedup missing pre-dedup R1 stream"))
                })?,
            ];
            if let Some(r2) = outputs.pre_dedup_r2.take() {
                pre_dedup_streams.push(r2);
            }

            let mut subsampled_streams = vec![outputs.non_host_r1];
            if let Some(r2) = outputs.non_host_r2 {
                subsampled_streams.push(r2);
            }

            let redup = run_redup(
                config.clone(),
                pre_dedup_streams,
                subsampled_streams,
                outputs.duplicate_clusters,
            )
                .await?;

            // Current AMR staging terminates RunRedup at its WDL-equivalent files.
            // These writers are installed immediately; downstream stages can later
            // consume redup.redups_streams directly instead without changing RunRedup.
            let mut redup_write_tasks = Vec::new();
            for (idx, stream) in redup.redups_streams.into_iter().enumerate() {
                redup_write_tasks.push(write_fasta_stream_to_file(
                    stream,
                    out_dir.join(format!("redups_{}.fa", idx + 1)),
                    config.clone(),
                    StreamDataType::JustBytes,
                    &format!("AMR RunRedup part {}", idx + 1),
                ));
            }
            let redup_completion_task = redup.completion_task;
            // ------------------------------------------------------------
            // Drain validation / QC results.
            // ------------------------------------------------------------

            let validation_and_qc_fut = async move {
                let stats =
                    join_with_error_handling(val_count_task).await?;

                info!(
                    "[AMR] validation: raw_reads={}, raw_bases={}, \
                     validated_reads={}, validated_bases={}, \
                     undersized={}, oversized={}",
                    stats.raw_reads,
                    stats.raw_bases,
                    stats.validated,
                    stats.validated_bases,
                    stats.undersized,
                    stats.oversized,
                );

                let qc_fastp_read_count =
                    match qc_count_result_rx.await {
                        Ok(Ok(count)) => {
                            info!(
                                "[AMR] fastp output reads: {}",
                                count
                            );
                            count
                        }

                        Ok(Err(e)) => {
                            error!(
                                "[AMR] error getting fastp \
                                 read count: {}",
                                e
                            );
                            0
                        }

                        Err(e) => {
                            error!(
                                "[AMR] fastp count receiver \
                                 dropped: {}",
                                e
                            );
                            0
                        }
                    };

                let ercc_filtered_count =
                    ercc_count_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "ERCC count receiver failed: {}",
                                e
                            ))
                        })?;

                Ok::<_, PipelineError>((
                    stats,
                    ercc_filtered_count,
                    qc_fastp_read_count,
                ))
            };

            // ------------------------------------------------------------
            // Drain BAM writers and host counts.
            // ------------------------------------------------------------

            let host_counts_fut = async move {
                ercc_bt2_bam_write_handle
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "ERCC BAM writer task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                host_bt2_bam_write_handle
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "Host Bowtie2 BAM writer task \
                             panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                if let Some(handle) =
                    optional_human_bam_write_handle
                {
                    handle
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Human BAM writer task \
                                 panicked: {}",
                                e
                            ))
                        })?
                        .map_err(PipelineError::Other)?;
                }

                let host_bt2_count =
                    host_bt2_count_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Host Bowtie2 count receiver \
                                 failed: {}",
                                e
                            ))
                        })?;

                let host_hisat2_count =
                    host_hisat2_count_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Host HISAT2 count receiver \
                                 failed: {}",
                                e
                            ))
                        })?;

                info!(
                    "[AMR] host filtering counts: \
                     bowtie2={}, hisat2={}",
                    host_bt2_count,
                    host_hisat2_count
                );

                Ok::<_, PipelineError>((
                    host_bt2_count,
                    host_hisat2_count,
                ))
            };

            // ------------------------------------------------------------
            // Drain Kallisto.
            // ------------------------------------------------------------

            let out_dir_kallisto = out_dir.clone();

            let kallisto_fut = async move {
                let kallisto_exit =
                    join_with_error_handling(
                        kallisto_exit_task
                    )
                        .await;

                match kallisto_exit {
                    Ok(_) => {
                        info!(
                            "[AMR] Kallisto completed successfully"
                        );
                    }

                    Err(e) => {
                        warn!(
                            "[AMR] Kallisto failed \
                             (non-fatal, possibly no ERCC \
                             spike-in): {}",
                            e
                        );
                    }
                }

                let results_task =
                    kallisto_results(
                        out_dir_kallisto.join("kallisto"),
                        kallisto_ercc_tx,
                    )
                        .await?;

                join_with_error_handling(results_task).await?;

                let results =
                    kallisto_ercc_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Kallisto results receiver \
                                 failed: {}",
                                e
                            ))
                        })?;

                info!(
                    "[AMR] Kallisto results drained: \
                     ercc_entries={}, transcript_mappings={}",
                    results.ercc_counts.len(),
                    results.transcript_to_gene.len(),
                );

                Ok::<(), PipelineError>(())
            };

            // ------------------------------------------------------------
            // Drain dedup / subsampling counts.
            // ------------------------------------------------------------

            let dedup_and_subsample_fut = async move {
                let unique_units =
                    dedup_count_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Dedup count receiver failed: {}",
                                e
                            ))
                        })?;

                let unique_reads =
                    unique_units * if paired { 2 } else { 1 };

                let subsample_units =
                    subsample_count_rx
                        .await
                        .map_err(|e| {
                            PipelineError::Other(anyhow!(
                                "Subsample count receiver \
                                 failed: {}",
                                e
                            ))
                        })?;

                let subsample_reads =
                    subsample_units * if paired { 2 } else { 1 };

                info!(
                    "[AMR] dedup/subsample: \
                     unique_reads={}, subsampled_reads={}",
                    unique_reads,
                    subsample_reads
                );

                Ok::<_, PipelineError>((
                    unique_reads,
                    subsample_reads,
                ))
            };

            // ------------------------------------------------------------
            // Complete all host-filter result branches.
            // ------------------------------------------------------------

            let (
                validation_and_qc,
                host_counts,
                (),
                dedup_and_subsample,
            ) = tokio::try_join!(
                validation_and_qc_fut,
                host_counts_fut,
                kallisto_fut,
                dedup_and_subsample_fut,
            )?;

            let (
                validation_stats,
                ercc_filtered_count,
                qc_fastp_read_count,
            ) = validation_and_qc;

            let (
                host_bt2_filtered_count,
                host_hisat2_filtered_count,
            ) = host_counts;

            let (
                unique_reads,
                subsample_reads,
            ) = dedup_and_subsample;

            // ------------------------------------------------------------
            // Write host-filter result artifacts.
            //
            // These exist only for raw AMR input because only raw input
            // executes the host-filter stage.
            // ------------------------------------------------------------

            write_count_file(
                &out_dir,
                "fastqs.count",
                validation_stats.raw_reads,
            )
                .await?;

            write_count_file(
                &out_dir,
                "validate_input_out.count",
                validation_stats.validated,
            )
                .await?;

            let validation_summary = serde_json::json!({
                "<50": validation_stats.validation_too_short,
                "50-500": validation_stats.validation_normal,
                "500-10000": validation_stats.validation_long,
                "10000+": validation_stats.validation_too_long,
            });

            let validation_summary_path =
                out_dir.join("validate_input_summary.json");

            tokio::fs::write(
                &validation_summary_path,
                serde_json::to_vec(
                    &validation_summary
                )
                    .map_err(|e| {
                        PipelineError::Other(e.into())
                    })?,
            )
                .await
                .map_err(|e| {
                    PipelineError::IOError(format!(
                        "Failed to write {}: {}",
                        validation_summary_path.display(),
                        e
                    ))
                })?;

            write_count_file(
                &out_dir,
                "fastp_out.count",
                qc_fastp_read_count,
            )
                .await?;

            write_count_file(
                &out_dir,
                "bowtie2_ercc_filtered_out.count",
                ercc_filtered_count,
            )
                .await?;

            write_count_file(
                &out_dir,
                "bowtie2_host_filtered_out.count",
                host_bt2_filtered_count,
            )
                .await?;

            write_count_file(
                &out_dir,
                "hisat2_host_filtered_out.count",
                host_hisat2_filtered_count,
            )
                .await?;

            write_count_file(
                &out_dir,
                "czid_dedup_out.count",
                unique_reads,
            )
                .await?;

            write_count_file(
                &out_dir,
                "subsampled_out.count",
                subsample_reads,
            )
                .await?;

            redup_completion_task
                .await
                .map_err(|e| PipelineError::Other(anyhow!("RunRedup task panicked: {}", e)))?
                .map_err(PipelineError::Other)?;

            for task in redup_write_tasks {
                task.await
                    .map_err(|e| PipelineError::Other(anyhow!("RunRedup writer task panicked: {}", e)))?
                    .map_err(PipelineError::Other)?;
            }

            // ------------------------------------------------------------
            // Drain cleanup owned by the host-filter stage.
            // ------------------------------------------------------------

            info!(
                "[AMR] draining {} cleanup tasks and \
                 {} cleanup receivers",
                cleanup_tasks.len(),
                cleanup_receivers.len()
            );

            for task in cleanup_tasks {
                task
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "AMR host-filter cleanup task \
                             panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;
            }

            for rx in cleanup_receivers {
                rx.await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "AMR host-filter cleanup receiver \
                             dropped: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;
            }

            // Keep host-filter temp directories alive until all
            // host-filter outputs and cleanup work have completed.
            info!(
                "[AMR] raw-input host filtering complete: \
                 {} retained temp directories can now be released",
                final_temp_dirs.len()
            );
        }

        AmrInputMode::Filtered => {
            // ------------------------------------------------------------
            // Already-filtered input.
            //
            // This corresponds to the WDL filtered_sample branch.
            // There is no host-filter invocation in that branch, so:
            //
            //   - run_host_filter() is not called
            //   - HostFilterCompletion does not exist
            //   - validation/QC/ERCC/host/dedup/subsample receivers
            //     do not exist
            //   - no host-filter count/report artifacts are written
            //
            // Later AMR stages will begin directly from the supplied
            // filtered-sample inputs.
            // ------------------------------------------------------------

            info!(
                "[AMR] filtered input selected; \
                 skipping host filtering"
            );


            let (file1_path, file2_path, _, _) =
                validate_file_inputs(&config, &cwd, true).await?;

            let pre_dedup_r1 = config
                .args
                .amr_pre_dedup_file1
                .as_deref()
                .ok_or_else(|| PipelineError::InvalidConfig(
                    "--amr-pre-dedup-file1 is required for --amr-input-mode filtered".to_string(),
                ))?;
            let pre_dedup_r1 = resolve_existing_input_path(pre_dedup_r1, &cwd, true)?;

            let pre_dedup_r2 = match config.args.amr_pre_dedup_file2.as_deref() {
                Some(path) => Some(resolve_existing_input_path(path, &cwd, true)?),
                None => None,
            };

            if file2_path.is_some() != pre_dedup_r2.is_some() {
                return Err(PipelineError::InvalidConfig(
                    "Filtered AMR RunRedup requires file2 and amr-pre-dedup-file2 to either both be present or both be absent".to_string(),
                ));
            }

            let clusters_path = config
                .args
                .amr_duplicate_clusters
                .as_deref()
                .ok_or_else(|| PipelineError::InvalidConfig(
                    "--amr-duplicate-clusters is required for --amr-input-mode filtered".to_string(),
                ))?;
            let clusters_path = resolve_existing_input_path(clusters_path, &cwd, true)?;

            let mut input_tasks = Vec::new();

            let (subsampled_r1, task) =
                fastq_file_to_record_stream(file1_path, config.clone(), "amr_filtered_subsampled_r1").await?;
            input_tasks.push(task);
            let mut subsampled_streams = vec![subsampled_r1];

            if let Some(path) = file2_path {
                let (stream, task) =
                    fastq_file_to_record_stream(path, config.clone(), "amr_filtered_subsampled_r2").await?;
                subsampled_streams.push(stream);
                input_tasks.push(task);
            }

            let (pre_r1, task) =
                fastq_file_to_record_stream(pre_dedup_r1, config.clone(), "amr_filtered_pre_dedup_r1").await?;
            input_tasks.push(task);
            let mut pre_dedup_streams = vec![pre_r1];

            if let Some(path) = pre_dedup_r2 {
                let (stream, task) =
                    fastq_file_to_record_stream(path, config.clone(), "amr_filtered_pre_dedup_r2").await?;
                pre_dedup_streams.push(stream);
                input_tasks.push(task);
            }

            let redup = run_redup(
                config.clone(),
                pre_dedup_streams,
                subsampled_streams,
                load_duplicate_clusters_csv(clusters_path),
            )
                .await?;

            let mut redup_write_tasks = Vec::new();
            for (idx, stream) in redup.redups_streams.into_iter().enumerate() {
                redup_write_tasks.push(write_fasta_stream_to_file(
                    stream,
                    out_dir.join(format!("redups_{}.fa", idx + 1)),
                    config.clone(),
                    StreamDataType::JustBytes,
                    &format!("AMR filtered RunRedup part {}", idx + 1),
                ));
            }

            redup.completion_task
                .await
                .map_err(|e| PipelineError::Other(anyhow!("Filtered RunRedup task panicked: {}", e)))?
                .map_err(PipelineError::Other)?;

            for task in input_tasks {
                task.await
                    .map_err(|e| PipelineError::Other(anyhow!("Filtered FASTQ adapter task panicked: {}", e)))?
                    .map_err(PipelineError::Other)?;
            }

            for task in redup_write_tasks {
                task.await
                    .map_err(|e| PipelineError::Other(anyhow!("Filtered RunRedup writer task panicked: {}", e)))?
                    .map_err(PipelineError::Other)?;
            }
        }
    }



    info!(
        "Finished AMR host-filter staging."
    );

    Ok(())
}
