//! Aligning functions
use std::path::PathBuf;
use std::sync::Arc;

use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation;
use tokio::sync::mpsc::Receiver;

use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use futures::stream::select_all;
use log::{debug, info, warn};
use tokio::fs;
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::sync::{mpsc, oneshot};
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;

use crate::config::defs::{RunConfig, PipelineError, MINIMAP2_TAG};
use crate::utils::fastx::SequenceRecord;
use crate::utils::command::{generate_cli};
use crate::utils::command::minimap2::Minimap2Config;
use crate::utils::streams::{ParseOutput, spawn_cmd};



/// Minimap2 mapping preset used with a prebuilt sharded reference.
///
/// IMPORTANT:
/// The `.mmi` files supplied to `sharded_minimap2_align()` must have been
/// built with index parameters compatible with the selected preset.
///
/// In particular, the current NT reference sets are expected to be:
///
/// - ShortRead: k=14, w=8
/// - MapOnt:    k=15, w=10
///
/// Those indexing parameters are intentionally NOT supplied here. They are
/// properties of the prebuilt reference assets and cannot be meaningfully
/// changed when mapping against an existing `.mmi`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Minimap2Preset {
    ShortRead,
    MapOnt,
}

impl Minimap2Preset {
    fn as_str(self) -> &'static str {
        match self {
            Self::ShortRead => "sr",
            Self::MapOnt => "map-ont",
        }
    }
}


/// Align contigs using rust-bio
///
/// # Arguments
///
/// * `contigs_rx` - stream of SequenceRecord
/// * `reference` - vec of ascii reference
///
/// # Returns
/// Result wrapping a Vec of tuple (position, Alignment Operations)
pub async fn align_contigs_to_reference(
    mut contigs_rx: Receiver<SequenceRecord>,
    reference: &[u8],
) -> anyhow::Result<Vec<(u64, Vec<AlignmentOperation>)>> {
    let mut alignments = Vec::new();
    let mut aligner = Aligner::new(-5, -1, |a, b| if a == b { 2 } else { -2 });

    while let Some(record) = contigs_rx.recv().await {
        if let SequenceRecord::Fasta { seq, .. } = record {
            let alignment = aligner.global(&seq, reference);
            let aligned_length = (alignment.yend - alignment.ystart) as u64;
            alignments.push((aligned_length, alignment.operations));
        }
    }
    if alignments.is_empty() {
        return Err(anyhow!("No valid FASTA records for alignment"));
    }
    Ok(alignments)
}




/// Configuration for sharded minimap2 alignment.
///
/// Every `.mmi` file in `index_dir` is treated as an independent reference
/// shard. The complete query input is aligned against every shard.
///
/// `query2` is intended for paired short-read input. Long-read/FASTA callers
/// normally provide only `query1`.
#[derive(Debug, Clone)]
pub struct ShardedMinimap2Config {
    pub index_dir: PathBuf,
    pub preset: Minimap2Preset,
    pub query1: PathBuf,
    pub query2: Option<PathBuf>,
}

/// Align a complete query set independently against a directory of prebuilt
/// minimap2 `.mmi` reference shards.
///
/// The function:
///
/// 1. discovers all `.mmi` files in `index_dir`;
/// 2. sorts them largest-first for better worker load balancing;
/// 3. estimates safe minimap2 concurrency from available RAM;
/// 4. creates a fixed pool of minimap2 workers;
/// 5. assigns reference shards round-robin to those workers;
/// 6. runs standard minimap2 against each shard;
/// 7. streams PAF records from all workers into one merged output stream.
///
/// The returned stream contains:
///
///     ParseOutput::Bytes
///
/// where each item is one newline-terminated PAF record.
///
/// No global primary/secondary reconciliation is performed here. This
/// intentionally preserves the architecture currently used by the Rust
/// short-read NT pipeline: all shard-produced candidate alignments are
/// retained for downstream M8 conversion and global hit processing.
///
/// # Failure semantics
///
/// Failure of ANY minimap2 shard is fatal. A nonzero minimap2 exit status is
/// returned from that worker as an error rather than being logged and ignored.
///
/// Worker and merger tasks are returned to the caller so that the pipeline's
/// normal cleanup/error-propagation machinery can await them.
pub async fn sharded_minimap2_align(
    config: Arc<RunConfig>,
    alignment: ShardedMinimap2Config,
) -> Result<
    (
        ReceiverStream<ParseOutput>,
        Vec<JoinHandle<Result<(), anyhow::Error>>>,
        Vec<oneshot::Receiver<Result<(), anyhow::Error>>>,
    ),
    PipelineError,
> {
    let mut cleanup_tasks: Vec<JoinHandle<Result<(), anyhow::Error>>> = Vec::new();

    // Retained for API compatibility with the pipeline cleanup model.
    // This utility currently creates no oneshot-based cleanup operations.
    let cleanup_receivers: Vec<
        oneshot::Receiver<Result<(), anyhow::Error>>,
    > = Vec::new();

    let base_buffer_size = config.base_buffer_size;
    let channel_buffer = (base_buffer_size * 8).max(1024);

    // =====================================================================
    // 1. Discover reference shards.
    // =====================================================================

    let mut chunk_paths: Vec<PathBuf> = Vec::new();
    let mut max_size_bytes: u64 = 0;

    let mut entries = fs::read_dir(&alignment.index_dir)
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Cannot read minimap2 split-index directory {}: {}",
                alignment.index_dir.display(),
                e
            ))
        })?;

    while let Some(entry) = entries
        .next_entry()
        .await
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "Failed reading minimap2 split-index directory {}: {}",
                alignment.index_dir.display(),
                e
            ))
        })?
    {
        let path = entry.path();

        if !path.is_file() {
            continue;
        }

        if !path
            .extension()
            .map_or(false, |ext| ext == "mmi")
        {
            continue;
        }

        if let Ok(meta) = entry.metadata().await {
            max_size_bytes = max_size_bytes.max(meta.len());
        }

        chunk_paths.push(path);
    }

    if chunk_paths.is_empty() {
        return Err(PipelineError::Other(anyhow!(
            "No .mmi reference shards found in {}",
            alignment.index_dir.display()
        )));
    }

    // Largest-first makes the round-robin assignment less likely to leave
    // one worker with the largest shard at the very end of the stage.
    //
    // Use the metadata already available from the filesystem here. Failure
    // to stat an individual shard is treated as size zero for ordering only;
    // minimap2 itself will subsequently fail loudly if the file is unusable.
    chunk_paths.sort_by_key(|path| {
        std::cmp::Reverse(
            path.metadata()
                .map(|meta| meta.len())
                .unwrap_or(0),
        )
    });

    let num_chunks = chunk_paths.len();

    info!(
        "[sharded_minimap2_align] found {} minimap2 shards in {} \
         (largest-first, preset={})",
        num_chunks,
        alignment.index_dir.display(),
        alignment.preset.as_str(),
    );

    // =====================================================================
    // 2. Determine worker concurrency.
    //
    // These values intentionally preserve the resource model proven by the
    // existing short-read NT implementation.
    // =====================================================================

    const TARGET_RAM_FRACTION: f64 = 0.55;
    const HARD_MAX_CONCURRENCY: usize = 14;
    const MIN_CONCURRENCY: usize = 4;
    const MIN_THREADS_PER_JOB: usize = 6;
    const MAX_THREADS_PER_JOB: usize = 12;

    let largest_gb = if max_size_bytes > 0 {
        ((max_size_bytes + (1 << 30) - 1) / (1 << 30)) as f64
    } else {
        // Conservative fallback inherited from the existing implementation.
        70.0
    };

    // Existing empirically tuned estimate from the short-read NT stage.
    let est_gb_per_job = (largest_gb * 0.85)
        .max(35.0)
        .min(65.0);

    info!(
        "[sharded_minimap2_align] largest shard: {:.1} GiB; \
         estimating {:.0} GiB RAM/minimap2 job",
        largest_gb,
        est_gb_per_job,
    );

    let available_ram_gb =
        config.available_ram as f64 / 1_000_000_000.0;

    let effective_min_concurrency = if available_ram_gb < 64.0 {
        1
    } else {
        MIN_CONCURRENCY
    };

    let ram_based = (
        available_ram_gb
            * TARGET_RAM_FRACTION
            / est_gb_per_job
    )
        .floor() as usize;

    let max_useful_concurrency =
        HARD_MAX_CONCURRENCY.min(num_chunks).max(1);

    let min_concurrency =
        effective_min_concurrency.min(max_useful_concurrency);

    let mut concurrency = ram_based
        .clamp(min_concurrency, max_useful_concurrency);

    // Additional safety valve inherited from the existing implementation.
    let estimated_peak = concurrency as f64 * est_gb_per_job;

    if estimated_peak > available_ram_gb * 0.82 {
        let new_concurrency = concurrency
            .saturating_sub(1)
            .max(min_concurrency);

        if new_concurrency < concurrency {
            warn!(
                "[sharded_minimap2_align] estimated peak {:.0} GiB exceeds \
                 82% of available RAM ({:.0} GiB); reducing concurrency \
                 {} -> {}",
                estimated_peak,
                available_ram_gb,
                concurrency,
                new_concurrency,
            );

            concurrency = new_concurrency;
        }
    }

    let threads_per_job = if available_ram_gb < 64.0 {
        config.max_cores.min(6).max(1)
    } else {
        (config.max_cores / concurrency.max(1))
            .clamp(
                MIN_THREADS_PER_JOB,
                MAX_THREADS_PER_JOB,
            )
    };

    // Recalculate after the safety adjustment.
    let estimated_peak =
        concurrency as f64 * est_gb_per_job;

    info!(
        concat!(
            "[sharded_minimap2_align] minimap2 stage\n",
            "  preset: {}\n",
            "  shards: {}\n",
            "  workers: {}\n",
            "  threads/job: {}\n",
            "  estimated peak RAM: ~{:.0} GiB ({:.0}% of available)\n",
            "  largest shard: {:.1} GiB\n",
            "  estimated RAM/job: {:.0} GiB"
        ),
        alignment.preset.as_str(),
        num_chunks,
        concurrency,
        threads_per_job,
        estimated_peak,
        if available_ram_gb > 0.0 {
            (estimated_peak / available_ram_gb * 100.0).round()
        } else {
            0.0
        },
        largest_gb,
        est_gb_per_job,
    );

    // =====================================================================
    // 3. Global merged PAF stream.
    // =====================================================================

    let (merged_tx, merged_rx) =
        mpsc::channel::<ParseOutput>(channel_buffer * 4);

    // Each worker receives:
    //
    //     (reference shard path, display name)
    //
    // and owns one output channel.
    let mut worker_job_txs:
        Vec<mpsc::Sender<(PathBuf, String)>> =
        Vec::with_capacity(concurrency);

    let mut worker_output_streams:
        Vec<ReceiverStream<ParseOutput>> =
        Vec::with_capacity(concurrency);

    let mut worker_handles:
        Vec<JoinHandle<Result<(), anyhow::Error>>> =
        Vec::with_capacity(concurrency);

    // =====================================================================
    // 4. Spawn fixed minimap2 worker pool.
    // =====================================================================

    for worker_idx in 0..concurrency {
        let (job_tx, mut job_rx) =
            mpsc::channel::<(PathBuf, String)>(
                num_chunks.max(1),
            );

        let (out_tx, out_rx) =
            mpsc::channel::<ParseOutput>(channel_buffer);

        worker_job_txs.push(job_tx);
        worker_output_streams.push(
            ReceiverStream::new(out_rx),
        );

        let config = config.clone();

        let query1 = alignment.query1.clone();
        let query2 = alignment.query2.clone();
        let preset = alignment.preset;

        let handle = tokio::spawn(async move {
            while let Some((chunk_mmi, chunk_name)) =
                job_rx.recv().await
            {
                let mut mm_config =
                    Minimap2Config::new()
                        .option(
                            "-t",
                            threads_per_job.to_string(),
                        )
                        .flag("-c")
                        .option("-x", preset.as_str())
                        .option("--secondary", "yes")
                        .positional(
                            chunk_mmi
                                .to_string_lossy()
                                .to_string(),
                        )
                        .positional(
                            query1
                                .to_string_lossy()
                                .to_string(),
                        );

                if let Some(query2_path) = &query2 {
                    mm_config = mm_config.positional(
                        query2_path
                            .to_string_lossy()
                            .to_string(),
                    );
                }

                let args = generate_cli(
                    MINIMAP2_TAG,
                    &config,
                    Some(&mm_config),
                )
                    .with_context(|| {
                        format!(
                            "Failed to generate minimap2 args \
                         for shard {}",
                            chunk_name
                        )
                    })?;

                debug!(
                    "[sharded_minimap2_align] worker {} \
                     launching minimap2 shard {} \
                     (preset={})",
                    worker_idx,
                    chunk_name,
                    preset.as_str(),
                );

                let (mut child, stderr_task) =
                    spawn_cmd(
                        config.clone(),
                        MINIMAP2_TAG,
                        args,
                        config.args.verbose,
                        None,
                    )
                        .await
                        .with_context(|| {
                            format!(
                                "Failed to spawn minimap2 \
                             for shard {}",
                                chunk_name
                            )
                        })?;

                // ---------------------------------------------------------
                // Stream PAF stdout.
                //
                // We intentionally retain the existing line-oriented PAF
                // behavior so downstream paf_to_m8 receives complete PAF
                // records as ParseOutput::Bytes.
                // ---------------------------------------------------------

                let mut downstream_closed = false;

                let stdout = child
                    .stdout
                    .take()
                    .ok_or_else(|| {
                        anyhow!(
                            "minimap2 stdout was not piped \
                             for shard {}",
                            chunk_name
                        )
                    })?;

                let mut reader = BufReader::new(stdout);
                let mut line = String::new();

                loop {
                    line.clear();

                    let n = reader
                        .read_line(&mut line)
                        .await
                        .with_context(|| {
                            format!(
                                "Failed reading minimap2 \
                                 stdout for shard {}",
                                chunk_name
                            )
                        })?;

                    if n == 0 {
                        break;
                    }

                    let trimmed = line.trim_end();

                    if trimmed.is_empty()
                        || trimmed.starts_with('#')
                    {
                        continue;
                    }

                    if downstream_closed {
                        // Continue draining stdout so minimap2 cannot block
                        // on a full pipe. We still need its exit status.
                        continue;
                    }

                    let mut bytes =
                        trimmed.as_bytes().to_vec();
                    bytes.push(b'\n');

                    if out_tx
                        .send(ParseOutput::Bytes(
                            Bytes::from(bytes),
                        ))
                        .await
                        .is_err()
                    {
                        downstream_closed = true;

                        debug!(
                            "[sharded_minimap2_align] \
                             downstream closed while worker {} \
                             was processing {}; continuing to \
                             drain minimap2 stdout",
                            worker_idx,
                            chunk_name,
                        );
                    }
                }

                // ---------------------------------------------------------
                // Always wait for minimap2 AND its stderr task.
                // ---------------------------------------------------------

                let status = child
                    .wait()
                    .await
                    .with_context(|| {
                        format!(
                            "Failed waiting for minimap2 \
                             shard {}",
                            chunk_name
                        )
                    })?;

                stderr_task
                    .await
                    .with_context(|| {
                        format!(
                            "Failed joining minimap2 stderr \
                             task for shard {}",
                            chunk_name
                        )
                    })??;

                // A missing shard must never silently turn into missing
                // candidate hits.
                if !status.success() {
                    return Err(anyhow!(
                        "minimap2 failed for shard {} \
                         with status {}",
                        chunk_name,
                        status
                    ));
                }

                info!(
                    "[sharded_minimap2_align] worker {} \
                     completed shard {}",
                    worker_idx,
                    chunk_name,
                );
            }

            Ok::<(), anyhow::Error>(())
        });

        worker_handles.push(handle);
    }

    // =====================================================================
    // 5. Merge worker PAF streams.
    //
    // select_all() intentionally permits records from different shards to
    // interleave. No downstream stage should depend on shard ordering.
    // =====================================================================

    let merger_handle = tokio::spawn({
        let mut merged_stream =
            select_all(worker_output_streams);

        async move {
            while let Some(item) =
                merged_stream.next().await
            {
                if merged_tx.send(item).await.is_err() {
                    // The computational consumer has disappeared.
                    // Dropping merged_tx will propagate closure toward the
                    // worker output channels.
                    break;
                }
            }

            Ok::<(), anyhow::Error>(())
        }
    });

    cleanup_tasks.push(merger_handle);

    // =====================================================================
    // 6. Assign shards to workers.
    //
    // chunk_paths is already largest-first. Round-robin assignment spreads
    // the largest shards across workers rather than concentrating them.
    // =====================================================================

    for (idx, chunk_mmi) in
        chunk_paths.into_iter().enumerate()
    {
        let chunk_name = chunk_mmi
            .file_name()
            .map(|name| {
                name.to_string_lossy().to_string()
            })
            .unwrap_or_else(|| {
                format!("chunk_{}", idx)
            });

        let worker_idx =
            idx % worker_job_txs.len();

        worker_job_txs[worker_idx]
            .send((chunk_mmi, chunk_name))
            .await
            .map_err(|_| {
                PipelineError::Other(anyhow!(
                    "minimap2 worker {} closed before \
                     shard {} could be scheduled",
                    worker_idx,
                    idx
                ))
            })?;
    }

    // Closing all job senders tells workers to exit after their assigned
    // shards have completed.
    drop(worker_job_txs);

    cleanup_tasks.extend(worker_handles);

    info!(
        "[sharded_minimap2_align] launched {} workers \
         across {} reference shards; merged PAF stream ready",
        concurrency,
        num_chunks,
    );

    Ok((
        ReceiverStream::new(merged_rx),
        cleanup_tasks,
        cleanup_receivers,
    ))
}

