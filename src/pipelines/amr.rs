//! AMR Pipeline

use std::sync::Arc;

use log::{error, info, warn};
use anyhow::anyhow;
use tempfile::TempDir;
use tokio::sync::oneshot;
use tokio::task::JoinHandle;

use crate::config::defs::{PipelineError, RunConfig};
use crate::pipelines::short_read_mngs::{
    run_host_filter,
    HostFilterCompletion,
    kallisto_results,
};
use crate::utils::file::{validate_file_inputs, write_count_file};
use crate::utils::streams::join_with_error_handling;



pub async fn run(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;
    let out_dir = config.out_dir.clone();

    let mut final_temp_dirs: Vec<TempDir> = Vec::new();

    let mut cleanup_tasks:
        Vec<JoinHandle<anyhow::Result<(), anyhow::Error>>> =
        Vec::new();

    let mut cleanup_receivers:
        Vec<oneshot::Receiver<anyhow::Result<(), anyhow::Error>>> =
        Vec::new();

    info!("Starting AMR pipeline.");
    info!("AMR Milestone B: running shared host filtering only.");

    // ------------------------------------------------------------
    // Input setup
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

    // ------------------------------------------------------------
    // Shared host filtering
    //
    // AMR differs from ordinary short-read mNGS here:
    // retain the exact host-filtered stream immediately before
    // deduplication.
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
    // Acknowledge the actual AMR preprocessing outputs.
    //
    // run_host_filter() has already joined dedup finalization and
    // pre-dedup retention before returning, so these are complete
    // files / structures here.
    // ------------------------------------------------------------

    let pre_dedup_r1 =
        outputs.pre_dedup_r1.as_ref().ok_or_else(|| {
            PipelineError::Other(anyhow!(
                "AMR requested pre-dedup retention but R1 was not returned"
            ))
        })?;

    if paired && outputs.pre_dedup_r2.is_none() {
        return Err(PipelineError::Other(anyhow!(
            "AMR paired-end run requested pre-dedup retention but R2 was not returned"
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
    // Drain validation / QC counts.
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

        // Match short-read behavior: QC counting failure is logged
        // and represented as zero rather than killing the run.
        let qc_fastp_read_count =
            match qc_count_result_rx.await {
                Ok(Ok(count)) => {
                    info!("[AMR] fastp output reads: {}", count);
                    count
                }

                Ok(Err(e)) => {
                    error!(
                        "[AMR] error getting fastp read count: {}",
                        e
                    );
                    0
                }

                Err(e) => {
                    error!(
                        "[AMR] fastp count receiver dropped: {}",
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
    //
    // For Milestone B we only need to guarantee these side branches
    // complete. We do not yet need downstream AMR processing from
    // these BAMs.
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
                    "Host Bowtie2 BAM writer task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;

        if let Some(handle) = optional_human_bam_write_handle {
            handle
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Human BAM writer task panicked: {}",
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
                        "Host Bowtie2 count receiver failed: {}",
                        e
                    ))
                })?;

        let host_hisat2_count =
            host_hisat2_count_rx
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "Host HISAT2 count receiver failed: {}",
                        e
                    ))
                })?;

        info!(
            "[AMR] host filtering counts: bowtie2={}, hisat2={}",
            host_bt2_count,
            host_hisat2_count
        );

        Ok::<_, PipelineError>((
            host_bt2_count,
            host_hisat2_count,
        ))
    };

    // ------------------------------------------------------------
    // Finish Kallisto.
    //
    // This branch is upstream host-filter bookkeeping even though
    // AMR itself does not consume Kallisto results.
    //
    // IMPORTANT: keep the receiver alive until the producer path is
    // settled.
    // ------------------------------------------------------------

    let out_dir_kallisto = out_dir.clone();

    let kallisto_fut = async move {
        let kallisto_exit =
            join_with_error_handling(kallisto_exit_task).await;

        match kallisto_exit {
            Ok(_) => {
                info!("[AMR] Kallisto completed successfully");
            }

            Err(e) => {
                warn!(
                "[AMR] Kallisto failed \
                 (non-fatal, possibly no ERCC spike-in): {}",
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
                "Kallisto results receiver failed: {}",
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
    // Dedup / subsampling counts.
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
                        "Subsample count receiver failed: {}",
                        e
                    ))
                })?;

        let subsample_reads =
            subsample_units * if paired { 2 } else { 1 };

        info!(
            "[AMR] dedup/subsample: unique_reads={}, subsampled_reads={}",
            unique_reads,
            subsample_reads
        );

        Ok::<_, PipelineError>((
            unique_reads,
            subsample_reads,
        ))
    };

    // Drain all independent HostFilterCompletion branches concurrently.
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
    // Persist the same basic host-filter count artifacts as mNGS.
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
        serde_json::to_vec(&validation_summary)
            .map_err(|e| PipelineError::Other(e.into()))?,
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
    // Drain general cleanup owned by run_host_filter().
    // ------------------------------------------------------------

    info!(
        "[AMR] draining {} cleanup tasks and {} cleanup receivers",
        cleanup_tasks.len(),
        cleanup_receivers.len()
    );

    for task in cleanup_tasks {
        task
            .await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "AMR host-filter cleanup task panicked: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;
    }

    for rx in cleanup_receivers {
        rx.await
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "AMR host-filter cleanup receiver dropped: {}",
                    e
                ))
            })?
            .map_err(PipelineError::Other)?;
    }

    // Keep final_temp_dirs alive through all of the above. They are
    // intentionally dropped only now, after AMR has acknowledged the
    // retained pre-dedup/non-host outputs.
    info!(
        "[AMR] host-filter-only milestone complete: \
         {} retained temp directories can now be released",
        final_temp_dirs.len()
    );

    info!("Finished AMR (cleanup completed, some non-fatal warnings may have occurred).");
    Ok(())
}