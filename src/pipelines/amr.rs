//! AMR Pipeline

use std::collections::{HashMap, HashSet};
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

use crate::config::defs::{ClusterInfo, PairingMode, PipelineError, RunConfig, StreamDataType, SPADES_TAG};
use crate::pipelines::short_read_mngs::{
    kallisto_results,
    run_host_filter,
    HostFilterCompletion,
};
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::generate_cli;
use crate::utils::fastx::{parse_byte_stream_to_fastq, read_fasta, read_fastq, write_fasta_stream_to_file, SequenceRecord};
use crate::utils::file::{choose_temp_dir, resolve_existing_input_path, validate_file_inputs, write_count_file};
use crate::utils::streams::{interleave_fastq_streams, join_with_error_handling, parse_child_output, spawn_cmd, stream_to_file, ChildStream, ParseMode, ParseOutput};
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
    mut host_filtered_stream: ReceiverStream<ParseOutput>,
    mut subsampled_stream: ReceiverStream<ParseOutput>,
    paired: bool,
    duplicate_clusters_task: JoinHandle<
        Result<Arc<DashMap<String, ClusterInfo>>, anyhow::Error>
    >,
) -> Result<RunRedupHandle, PipelineError> {
    let part_count = if paired { 2 } else { 1 };

    let mut output_streams = Vec::with_capacity(part_count);
    let mut output_txs = Vec::with_capacity(part_count);

    for _ in 0..part_count {
        let (tx, rx) =
            mpsc::channel::<ParseOutput>(config.base_buffer_size);

        output_txs.push(tx);
        output_streams.push(ReceiverStream::new(rx));
    }

    // Drain the pre-dedup stream immediately. It is upstream of dedup and
    // therefore cannot be left unread while dedup/subsampling finalize.
    //
    // Keep records separated by end so duplicate members can later be
    // appended to the corresponding redups_N.fa output.
    let host_buffer_task = tokio::spawn(async move {
        let mut records_by_part =
            (0..part_count)
                .map(|_| Vec::<SequenceRecord>::new())
                .collect::<Vec<_>>();

        let mut record_index = 0usize;

        while let Some(item) = host_filtered_stream.next().await {
            let ParseOutput::Fastq(record) = item else {
                return Err(anyhow!(
                    "RunRedup expected FASTQ on interleaved host-filtered stream"
                ));
            };

            let part_index = if paired {
                record_index % 2
            } else {
                0
            };

            records_by_part[part_index].push(record);
            record_index += 1;
        }

        if paired && record_index % 2 != 0 {
            return Err(anyhow!(
                "RunRedup received incomplete paired pre-dedup stream: {} records",
                record_index
            ));
        }

        Ok::<Vec<Vec<SequenceRecord>>, anyhow::Error>(records_by_part)
    });

    // Forward the complete subsampled input immediately as FASTA while
    // recording which deduplicated representatives survived subsampling.
    let prefix_task = tokio::spawn(async move {
        let mut passed_representatives = HashSet::<String>::new();
        let mut record_index = 0usize;

        while let Some(item) = subsampled_stream.next().await {
            let ParseOutput::Fastq(record) = item else {
                return Err(anyhow!(
                    "RunRedup expected FASTQ on interleaved subsampled stream"
                ));
            };

            let part_index = if paired {
                record_index % 2
            } else {
                0
            };

            passed_representatives.insert(record.id().to_string());

            let fasta = SequenceRecord::Fasta {
                id: record.id().to_string(),
                desc: None,
                seq: record.seq(),
            };

            output_txs[part_index]
                .send(ParseOutput::Fasta(fasta))
                .await
                .map_err(|_| {
                    anyhow!(
                        "RunRedup downstream receiver dropped for part {}",
                        part_index + 1
                    )
                })?;

            record_index += 1;
        }

        if paired && record_index % 2 != 0 {
            return Err(anyhow!(
                "RunRedup received incomplete paired subsampled stream: {} records",
                record_index
            ));
        }

        Ok::<_, anyhow::Error>((
            passed_representatives,
            output_txs,
        ))
    });

    let completion_task = tokio::spawn(async move {
        // Both live inputs already have consumers installed above.
        let duplicate_clusters = duplicate_clusters_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunRedup duplicate-cluster task panicked: {}",
                    e
                )
            })??;

        let (passed_representatives, output_txs) = prefix_task
            .await
            .map_err(|e| {
                anyhow!("RunRedup prefix task panicked: {}", e)
            })??;

        let host_records_by_part = host_buffer_task
            .await
            .map_err(|e| {
                anyhow!("RunRedup host-buffer task panicked: {}", e)
            })??;

        let mut duplicate_member_ids =
            HashSet::<String>::new();

        for entry in duplicate_clusters.iter() {
            let rep_id = entry.key();
            let info = entry.value();

            if info.size <= 1
                || !passed_representatives.contains(rep_id.as_str())
            {
                continue;
            }

            for member_id in &info.members {
                if member_id != rep_id {
                    duplicate_member_ids.insert(member_id.clone());
                }
            }
        }

        for (part_index, host_records) in
            host_records_by_part.into_iter().enumerate()
        {
            let tx = &output_txs[part_index];

            // Preserve original pre-dedup order within each end.
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
                    .map_err(|_| {
                        anyhow!(
                            "RunRedup downstream receiver dropped for part {}",
                            part_index + 1
                        )
                    })?;
            }
        }

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunRedupHandle {
        redups_streams: output_streams,
        completion_task,
    })
}

#[derive(Debug)]
struct RunSpadesHandle {
    contigs_path: PathBuf,
    scaffolds_path: PathBuf,
    completion_task: JoinHandle<Result<(), anyhow::Error>>,
}

/// Native AMR equivalent of the WDL RunSpades task.
///
/// RunRedup produces live FASTA streams, but SPAdes inspects and reopens its
/// read inputs during a run. Materializing the inputs is therefore a genuine
/// external-tool boundary rather than a pipeline checkpoint. The per-end
/// writers are installed before this function returns so RunRedup is never
/// left blocked on an unconsumed bounded output.
async fn run_spades_amr(
    config: Arc<RunConfig>,
    reduplicated_streams: Vec<ReceiverStream<ParseOutput>>,
) -> Result<RunSpadesHandle, PipelineError> {
    if reduplicated_streams.is_empty() || reduplicated_streams.len() > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunSpades expected one or two reduplicated FASTA streams, got {}",
            reduplicated_streams.len()
        )));
    }

    let input_temp_dir = choose_temp_dir(
        config.input_size,
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        2,
        false,
    )
        .await?;

    let spades_out_dir = config.out_dir.join("spades");
    let contigs_path = spades_out_dir.join("contigs.fasta");
    let scaffolds_path = spades_out_dir.join("scaffolds.fasta");

    let mut input_paths = Vec::with_capacity(reduplicated_streams.len());
    let mut input_write_tasks = Vec::with_capacity(reduplicated_streams.len());

    for (idx, stream) in reduplicated_streams.into_iter().enumerate() {
        let path = input_temp_dir
            .path()
            .join(format!("redups_{}.fa", idx + 1));

        input_write_tasks.push(write_fasta_stream_to_file(
            stream,
            path.clone(),
            config.clone(),
            StreamDataType::JustBytes,
            &format!("AMR RunSpades input part {}", idx + 1),
        ));
        input_paths.push(path);
    }

    let completion_config = config.clone();
    let completion_contigs = contigs_path.clone();
    let completion_scaffolds = scaffolds_path.clone();

    let completion_task = tokio::spawn(async move {
        // Keep the scratch directory alive for the entire SPAdes invocation.
        let _input_temp_dir = input_temp_dir;

        for (idx, task) in input_write_tasks.into_iter().enumerate() {
            task.await
                .map_err(|e| anyhow!(
                    "RunSpades input writer {} panicked: {}",
                    idx + 1,
                    e
                ))??;
        }

        // A RunRedup/input-adapter failure is an upstream pipeline failure.
        // From this point onward, match WDL RunSpades: assembly failures are
        // non-fatal and produce compatibility marker files.
        async fn write_failure_outputs(
            contigs_path: &PathBuf,
            scaffolds_path: &PathBuf,
        ) -> Result<(), anyhow::Error> {
            if let Some(parent) = contigs_path.parent() {
                tokio::fs::create_dir_all(parent).await?;
            }
            tokio::fs::write(contigs_path, b";ASSEMBLY FAILED\n").await?;
            tokio::fs::write(scaffolds_path, b";ASSEMBLY FAILED\n").await?;
            Ok(())
        }

        if spades_out_dir.exists() {
            tokio::fs::remove_dir_all(&spades_out_dir).await?;
        }

        let mut options = HashMap::new();
        options.insert("--only-assembler".to_string(), None);

        let spades_config = SpadesConfig {
            r1_path: input_paths[0].clone(),
            r2_path_opt: input_paths.get(1).cloned(),
            outdir_path: spades_out_dir.clone(),
            tempdir_path: Some(_input_temp_dir.path().to_path_buf()),
            option_fields: options,
        };

        let spades_args = match generate_cli(
            SPADES_TAG,
            &completion_config,
            Some(&spades_config),
        ) {
            Ok(args) => args,
            Err(e) => {
                warn!("[AMR] RunSpades command generation failed (non-fatal): {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        };

        info!("[AMR] RunSpades args: {:?}", spades_args);

        let (mut child, stderr_task) = match spawn_cmd(
            completion_config.clone(),
            SPADES_TAG,
            spades_args,
            completion_config.args.verbose,
            None,
        )
            .await
        {
            Ok(v) => v,
            Err(e) => {
                warn!("[AMR] RunSpades spawn failed (non-fatal): {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        };

        // Drain stdout concurrently so a verbose SPAdes process cannot block on
        // its pipe. Preserve it as a useful sidecar log.
        let stdout_task = match parse_child_output(
            &mut child,
            ChildStream::Stdout,
            ParseMode::Lines,
            &completion_config,
        )
            .await
        {
            Ok(stdout_rx) => Some(tokio::spawn(stream_to_file(
                stdout_rx,
                completion_config.out_dir.join("amr_spades_stdout.log"),
            ))),
            Err(e) => {
                warn!("[AMR] RunSpades stdout capture failed (non-fatal): {}", e);
                None
            }
        };

        let status = child.wait().await;
        let stderr_result = stderr_task.await;
        if let Some(task) = stdout_task {
            let _ = task.await;
        }

        let spades_ok = matches!(status, Ok(ref s) if s.success())
            && matches!(stderr_result, Ok(Ok(())));

        if !spades_ok || !completion_contigs.is_file() {
            warn!(
                "[AMR] RunSpades assembly failed or contigs.fasta is missing; writing WDL-compatible marker"
            );
            write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
            return Ok(());
        }

        let metadata = tokio::fs::metadata(&completion_contigs).await?;
        if metadata.len() == 0 {
            warn!("[AMR] RunSpades produced empty contigs.fasta; writing failure marker");
            write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
            return Ok(());
        }

        // Match `seqtk seq -L min_contig_length`: read the completed SPAdes
        // FASTA and rewrite contigs.fasta with records at or above the configured
        // minimum length. Use a separate path to avoid truncating the input.
        let unfiltered_contigs = spades_out_dir.join("contigs_unfiltered.fasta");
        if let Err(e) = tokio::fs::rename(&completion_contigs, &unfiltered_contigs).await {
            warn!("[AMR] RunSpades could not stage contigs for filtering: {}", e);
            write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
            return Ok(());
        }

        let filtered_rx = match read_fasta(
            unfiltered_contigs.clone(),
            u64::MAX,
            Some(completion_config.args.min_contig_length),
            None,
            &completion_config,
        ) {
            Ok(rx) => rx,
            Err(e) => {
                warn!("[AMR] RunSpades contig filtering read failed: {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        };

        let filter_write_task = write_fasta_stream_to_file(
            ReceiverStream::new(filtered_rx),
            completion_contigs.clone(),
            completion_config.clone(),
            StreamDataType::JustBytes,
            "AMR RunSpades contig filtering",
        );

        match filter_write_task.await {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                warn!("[AMR] RunSpades contig filtering write failed: {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
            Err(e) => {
                warn!("[AMR] RunSpades contig filtering writer panicked: {}", e);
                write_failure_outputs(&completion_contigs, &completion_scaffolds).await?;
                return Ok(());
            }
        }

        let _ = tokio::fs::remove_file(unfiltered_contigs).await;

        info!(
            "[AMR] RunSpades complete: contigs={}, scaffolds={}",
            completion_contigs.display(),
            completion_scaffolds.display()
        );

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunSpadesHandle {
        contigs_path,
        scaffolds_path,
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
            if outputs.pre_dedup_interleaved.is_none() {
                return Err(PipelineError::Other(anyhow!(
        "AMR requested retained interleaved pre-dedup stream but it was not returned"
    )));
            }

            if outputs.subsampled_interleaved.is_none() {
                return Err(PipelineError::Other(anyhow!(
        "AMR requested retained interleaved subsampled stream but it was not returned"
    )));
            }

            info!(
    "[AMR] host-filter streaming outputs ready: \
     pre_dedup_interleaved=true, \
     subsampled_interleaved=true, \
     non_host_r1=true, non_host_r2={}, \
     duplicate_clusters=pending",
    outputs.non_host_r2.is_some(),
);



            // Install RunRedup consumers immediately so every bounded host-filter
            // output has a downstream reader before we await host-filter completion.
            let pre_dedup_stream = outputs.pre_dedup_interleaved.take().ok_or_else(|| {
                PipelineError::Other(anyhow!(
        "AMR raw RunRedup missing retained interleaved pre-dedup stream"
    ))
            })?;

            let subsampled_stream = outputs.subsampled_interleaved.take().ok_or_else(|| {
                PipelineError::Other(anyhow!(
        "AMR raw RunRedup missing retained interleaved subsampled stream"
    ))
            })?;

            let redup = run_redup(
                config.clone(),
                pre_dedup_stream,
                subsampled_stream,
                paired,
                outputs.duplicate_clusters,
            )
                .await?;




            // Raw WDL topology is RunRedupRaw -> RunSpades. Install the SPAdes
            // adapter immediately so every bounded RunRedup output has a consumer
            // before any host-filter completion work is awaited.
            let spades = run_spades_amr(
                config.clone(),
                redup.redups_streams,
            )
                .await?;
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

            spades.completion_task
                .await
                .map_err(|e| PipelineError::Other(anyhow!("RunSpades task panicked: {}", e)))?
                .map_err(PipelineError::Other)?;

            info!(
                "[AMR] raw contigs ready: {} (scaffolds: {})",
                spades.contigs_path.display(),
                spades.scaffolds_path.display()
            );

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

            let paired = file2_path.is_some();

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

            let subsampled_stream = if let Some(path) = file2_path {
                let (subsampled_r2, task) =
                    fastq_file_to_record_stream(path, config.clone(), "amr_filtered_subsampled_r2").await?;
                input_tasks.push(task);

                let (interleaved_rx, interleave_task) = interleave_fastq_streams(
                    subsampled_r1.into_inner(),
                    subsampled_r2.into_inner(),
                    config.base_buffer_size * 4,
                )
                    .await?;
                input_tasks.push(interleave_task);
                ReceiverStream::new(interleaved_rx)
            } else {
                subsampled_r1
            };

            let (pre_r1, task) =
                fastq_file_to_record_stream(
                    pre_dedup_r1,
                    config.clone(),
                    "amr_filtered_pre_dedup_r1",
                )
                    .await?;

            input_tasks.push(task);

            let pre_dedup_stream = if let Some(path) = pre_dedup_r2 {
                let (pre_r2, task) =
                    fastq_file_to_record_stream(
                        path,
                        config.clone(),
                        "amr_filtered_pre_dedup_r2",
                    )
                        .await?;

                input_tasks.push(task);

                let (interleaved_rx, interleave_task) =
                    interleave_fastq_streams(
                        pre_r1.into_inner(),
                        pre_r2.into_inner(),
                        config.base_buffer_size * 4,
                    )
                        .await?;

                input_tasks.push(interleave_task);

                ReceiverStream::new(interleaved_rx)
            } else {
                pre_r1
            };

            let redup = run_redup(
                config.clone(),
                pre_dedup_stream,
                subsampled_stream,
                paired,
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
