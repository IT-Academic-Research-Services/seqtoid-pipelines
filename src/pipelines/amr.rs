//! AMR Pipeline

use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::anyhow;
use log::{error, info, warn};
use tempfile::TempDir;
use bytes::Bytes;
use dashmap::DashMap;
use tokio::fs::File as TokioFile;
use tokio::io::BufReader;
use tokio::sync::{mpsc, oneshot};
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;

use crate::config::defs::{ClusterInfo, PipelineError, RunConfig, StreamDataType};
use crate::pipelines::short_read_mngs::{
    kallisto_results,
    run_host_filter,
    HostFilterCompletion,
};
use crate::utils::fastx::SequenceRecord;
use crate::utils::file::{validate_file_inputs, write_count_file};
use crate::utils::streams::{join_with_error_handling, parse_fastq, ParseOutput};


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
    host_filtered_reads: Vec<PathBuf>,
    subsampled_reads: Vec<PathBuf>,
    duplicate_clusters: Arc<DashMap<String, ClusterInfo>>,
) -> Result<RunRedupHandle, PipelineError> {
    if host_filtered_reads.is_empty() || host_filtered_reads.len() > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRedup expected one or two host-filtered read files, got {}",
            host_filtered_reads.len()
        )));
    }

    if host_filtered_reads.len() != subsampled_reads.len() {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRedup input cardinality mismatch: host_filtered_reads={}, subsampled_reads={}",
            host_filtered_reads.len(),
            subsampled_reads.len()
        )));
    }

    let mut output_streams = Vec::with_capacity(subsampled_reads.len());
    let mut producer_tasks = Vec::with_capacity(subsampled_reads.len());

    for (part_index, (host_filtered_path, subsampled_path)) in host_filtered_reads
        .into_iter()
        .zip(subsampled_reads.into_iter())
        .enumerate()
    {
        let (tx, rx) = mpsc::channel::<ParseOutput>(config.base_buffer_size);
        output_streams.push(ReceiverStream::new(rx));

        let config_for_task = config.clone();
        let clusters = duplicate_clusters.clone();
        let task = tokio::spawn(async move {
            // -------------------------------------------------------------
            // Prefix: subsampled FASTQ -> FASTA, preserving file order.
            // -------------------------------------------------------------
            let subsampled_file = TokioFile::open(&subsampled_path)
                .await
                .map_err(|e| anyhow!(
                    "RunRedup failed to open subsampled input {}: {}",
                    subsampled_path.display(), e
                ))?;

            let subsampled_rx = parse_fastq(
                BufReader::new(subsampled_file),
                &config_for_task,
                StreamDataType::IlluminaFastq,
            )
                .await
                .map_err(|e| anyhow!(
                "RunRedup failed to parse subsampled input {}: {}",
                subsampled_path.display(), e
            ))?;

            let mut subsampled_stream = ReceiverStream::new(subsampled_rx);
            let mut passed_representatives = HashSet::<String>::new();

            while let Some(item) = subsampled_stream.next().await {
                let ParseOutput::Fastq(record) = item else {
                    return Err(anyhow!(
                        "RunRedup expected parsed FASTQ from {}",
                        subsampled_path.display()
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

            // WDL builds one passed_filters set across all ends. For strict
            // paired-end data both ends carry the same read IDs, so the per-end
            // set is equivalent and avoids coupling the bounded output streams
            // behind a cross-stream barrier.
            let mut duplicate_member_ids = HashSet::<String>::new();

            for entry in clusters.iter() {
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

            if duplicate_member_ids.is_empty() {
                return Ok::<(), anyhow::Error>(());
            }

            // -------------------------------------------------------------
            // Suffix: scan this pre-dedup host-filtered end in its original
            // order, select duplicate members, and convert them to FASTA.
            // This matches `seqkit grep ... | seqkit fq2fa` ordering.
            // -------------------------------------------------------------
            let host_file = TokioFile::open(&host_filtered_path)
                .await
                .map_err(|e| anyhow!(
                    "RunRedup failed to open host-filtered input {}: {}",
                    host_filtered_path.display(), e
                ))?;

            let host_rx = parse_fastq(
                BufReader::new(host_file),
                &config_for_task,
                StreamDataType::IlluminaFastq,
            )
                .await
                .map_err(|e| anyhow!(
                "RunRedup failed to parse host-filtered input {}: {}",
                host_filtered_path.display(), e
            ))?;

            let mut host_stream = ReceiverStream::new(host_rx);

            while let Some(item) = host_stream.next().await {
                let ParseOutput::Fastq(record) = item else {
                    return Err(anyhow!(
                        "RunRedup expected parsed FASTQ from {}",
                        host_filtered_path.display()
                    ));
                };

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

            Ok::<(), anyhow::Error>(())
        });

        producer_tasks.push(task);
    }

    let completion_task = tokio::spawn(async move {
        for task in producer_tasks {
            task.await
                .map_err(|e| anyhow!("RunRedup producer task panicked: {}", e))??;
        }
        Ok::<(), anyhow::Error>(())
    });

    Ok(RunRedupHandle {
        redups_streams: output_streams,
        completion_task,
    })
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
                outputs,

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

            let pre_dedup_r1 =
                outputs
                    .pre_dedup_r1
                    .as_ref()
                    .ok_or_else(|| {
                        PipelineError::Other(anyhow!(
                            "AMR requested pre-dedup retention \
                             but R1 was not returned"
                        ))
                    })?;

            if paired && outputs.pre_dedup_r2.is_none() {
                return Err(PipelineError::Other(anyhow!(
                    "AMR paired-end run requested pre-dedup \
                     retention but R2 was not returned"
                )));
            }

            info!(
                "[AMR] host-filter outputs ready: \
                 pre_dedup_r1={}, pre_dedup_r2={:?}, \
                 non_host_r1={}, non_host_r2={:?}, \
                 duplicate_clusters={}",
                pre_dedup_r1.display(),
                outputs
                    .pre_dedup_r2
                    .as_ref()
                    .map(|p| p.display().to_string()),
                outputs.non_host_r1.display(),
                outputs
                    .non_host_r2
                    .as_ref()
                    .map(|p| p.display().to_string()),
                outputs.duplicate_clusters.len(),
            );

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
        }
    }



    info!(
        "Finished AMR host-filter staging."
    );

    Ok(())
}