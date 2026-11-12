//! AMR Pipeline

use crate::pipelines::short_read_mngs::HostFilterOutputs;
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

use crate::config::defs::{ClusterInfo, PairingMode, PipelineError, RunConfig, StreamDataType,
                          SPADES_TAG, KMA_TAG, SAMTOOLS_TAG, BAMTOOLS_TAG, SamtoolsSubcommand, BamtoolsSubcommand};
use crate::pipelines::short_read_mngs::{
    kallisto_results,
    run_host_filter,
    HostFilterCompletion,
};
use crate::utils::command::kma::KmaConfig;
use crate::utils::command::samtools::SamtoolsConfig;
use crate::utils::command::bamtools::BamtoolsConfig;
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::generate_cli;
use crate::utils::fastx::{parse_byte_stream_to_fastq, read_fasta, read_fastq, write_fasta_stream_to_file, SequenceRecord};
use crate::utils::file::{choose_temp_dir, resolve_existing_input_path, validate_file_inputs, write_count_file, write_byte_stream_to_file};
use crate::utils::sambam::parse_sam_record;
use crate::utils::streams::{
    bytes_to_lines,
    fanout_to_channels,
    interleave_fastq_streams,
    join_with_error_handling,
    parse_bytes,
    parse_child_output,
    spawn_cmd,
    stream_to_file,
    stream_to_cmd,
    ChildStream,
    ParseMode,
    ParseOutput,
};
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

#[derive(Debug)]
struct RunKmaHandle {
    /// Complete line-oriented KMA SAM for the native AMR parser.
    sam_stream: ReceiverStream<ParseOutput>,

    /// Raw KMA SAM bytes, unchanged, for direct streaming into samtools fixmate.
    samtools_stream: mpsc::Receiver<ParseOutput>,

    /// Completes when every raw KMA stdout byte has reached both branches.
    sam_fanout_done: oneshot::Receiver<Result<(), anyhow::Error>>,

    /// Completes when the parser branch has converted raw byte chunks to lines.
    sam_lines_task: JoinHandle<Result<(), anyhow::Error>>,

    completion_task: JoinHandle<Result<(), anyhow::Error>>,

    // Keep KMA sidecars alive for downstream AMR post-processing.
    output_temp_dir: TempDir,
    output_prefix: PathBuf,
}

#[derive(Debug, Clone)]
struct KmaReference {
    name: String,
    length: u64,
}

#[derive(Debug, Clone)]
struct KmaAlignment {
    qname: String,
    flag: u16,
    rname: String,
    pos: u64,
    mapq: u8,
    cigar: String,
    rnext: String,
    seq: String,
}

#[derive(Debug)]
struct KmaSamResults {
    references: HashMap<String, KmaReference>,
    alignments: Vec<KmaAlignment>,

    /// Number of non-header SAM records emitted by KMA.
    total_records: u64,

    /// KMA records that are not usable alignments despite potentially
    /// misleading FLAG values, e.g. RNAME="*", POS=0, CIGAR="*".
    invalid_alignment_records: u64,
}

#[derive(Debug, Clone)]
struct RgiReadCounts {
    reference_length: u64,
    mapped: u64,
    unmapped: u64,
    all: u64,
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

async fn run_kma_amr(
    config: Arc<RunConfig>,
    reduplicated_streams: Vec<ReceiverStream<ParseOutput>>,
    kma_db: PathBuf,
) -> Result<RunKmaHandle, PipelineError> {
    let stream_count = reduplicated_streams.len();

    if stream_count == 0 || stream_count > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunKMA expected one or two reduplicated FASTA streams, got {}",
            stream_count
        )));
    }

    let paired = stream_count == 2;

    // ------------------------------------------------------------
    // Adapt RunRedup output to KMA stdin.
    //
    // Single-end:
    //
    //     R1 -> KMA -i --
    //
    // Paired-end:
    //
    //     R1 --\
    //            -> R1,R2,R1,R2,... -> KMA -int --
    //     R2 --/
    //
    // RunRedup emits ParseOutput::Fasta records. Keep them structured
    // until stream_to_cmd serializes them to KMA stdin.
    // ------------------------------------------------------------

    let (kma_input_rx, input_adapter_task) = if paired {
        let mut streams = reduplicated_streams.into_iter();

        let mut r1 = streams
            .next()
            .expect("paired RunKMA input missing R1 stream");

        let mut r2 = streams
            .next()
            .expect("paired RunKMA input missing R2 stream");

        let (tx, rx) =
            mpsc::channel::<ParseOutput>(config.base_buffer_size);

        let task = tokio::spawn(async move {
            let mut pairs = 0u64;

            loop {
                let left = r1.next().await;
                let right = r2.next().await;

                match (left, right) {
                    (
                        Some(ParseOutput::Fasta(r1_record)),
                        Some(ParseOutput::Fasta(r2_record)),
                    ) => {
                        tx.send(ParseOutput::Fasta(r1_record))
                            .await
                            .map_err(|_| {
                                anyhow!(
                                    "RunKMA paired input receiver dropped on R1"
                                )
                            })?;

                        tx.send(ParseOutput::Fasta(r2_record))
                            .await
                            .map_err(|_| {
                                anyhow!(
                                    "RunKMA paired input receiver dropped on R2"
                                )
                            })?;

                        pairs += 1;
                    }

                    (None, None) => {
                        break;
                    }

                    (None, Some(_)) | (Some(_), None) => {
                        return Err(anyhow!(
                            "RunKMA R1/R2 streams have unequal record counts \
                             after {} complete pairs",
                            pairs
                        ));
                    }

                    (Some(left), Some(right)) => {
                        return Err(anyhow!(
                            "RunKMA expected paired FASTA records, got {:?} and {:?}",
                            left,
                            right
                        ));
                    }
                }
            }

            info!(
                "[AMR] RunKMA prepared {} interleaved read pairs",
                pairs
            );

            Ok::<(), anyhow::Error>(())
        });

        (rx, task)
    } else {
        let mut streams = reduplicated_streams.into_iter();

        let mut r1 = streams
            .next()
            .expect("single-end RunKMA input missing R1 stream");

        let (tx, rx) =
            mpsc::channel::<ParseOutput>(config.base_buffer_size);

        let task = tokio::spawn(async move {
            let mut reads = 0u64;

            while let Some(item) = r1.next().await {
                let ParseOutput::Fasta(record) = item else {
                    return Err(anyhow!(
                        "RunKMA expected single-end FASTA record, got {:?}",
                        item
                    ));
                };

                tx.send(ParseOutput::Fasta(record))
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "RunKMA single-end input receiver dropped"
                        )
                    })?;

                reads += 1;
            }

            info!(
                "[AMR] RunKMA prepared {} single-end reads",
                reads
            );

            Ok::<(), anyhow::Error>(())
        });

        (rx, task)
    };

    // ------------------------------------------------------------
    // KMA requires an output prefix even though SAM is emitted on
    // stdout.
    //
    // Keep the TempDir in RunKmaHandle so KMA sidecars remain alive
    // for downstream RunRgiBwtKma processing.
    // ------------------------------------------------------------

    let kma_temp_dir =
        tempfile::TempDir::new_in(&config.ram_temp_dir)
            .map_err(|e| {
                PipelineError::Other(anyhow!(
                    "Failed to create RunKMA temp directory: {}",
                    e
                ))
            })?;

    let kma_output_prefix =
        kma_temp_dir.path().join("read_amr");

    let kma_config = KmaConfig::amr_stdin(
        kma_db,
        kma_output_prefix.clone(),
        paired,
    );

    let kma_args = generate_cli(
        KMA_TAG,
        &config,
        Some(&kma_config),
    )
        .map_err(PipelineError::Other)?;

    info!(
        "[AMR] RunKMA mode={}, args={:?}",
        if paired {
            "paired-interleaved-stdin"
        } else {
            "single-end-stdin"
        },
        kma_args
    );

    let stderr_path =
        config.out_dir.join("amr_kma_stderr.log");

    // ------------------------------------------------------------
    // Feed structured FASTA directly to KMA.
    // ------------------------------------------------------------

    let (child, stdin_task, stderr_task) =
        stream_to_cmd(
            config.clone(),
            kma_input_rx,
            KMA_TAG,
            kma_args,
            StreamDataType::IlluminaFastq,
            config.args.verbose,
            Some(stderr_path),
        )
            .await
            .map_err(PipelineError::Other)?;

    // ------------------------------------------------------------
    // Take KMA stdout without retaining the child mutex.
    //
    // stream_to_cmd's stdin task ultimately waits on the child.
    // Holding the child lock while consuming stdout would deadlock.
    // ------------------------------------------------------------

    let stdout = {
        let mut child_guard = child.lock().await;

        child_guard
            .stdout
            .take()
            .ok_or_else(|| {
                PipelineError::Other(anyhow!(
                    "RunKMA child stdout not available"
                ))
            })?
    };

    // ------------------------------------------------------------
    // KMA `-sam` emits serialized SAM on stdout.
    //
    // IMPORTANT:
    //
    // Keep this as RAW BYTES before fanout. This preserves the exact
    // SAM byte stream, including line endings, for samtools.
    //
    //                         ┌-> bytes_to_lines()
    //                         │      -> native SAM parser
    // KMA -> parse_bytes() -> fanout
    //                         │
    //                         └-> samtools fixmate stdin
    //
    // There is deliberately no raw-SAM file boundary here.
    // ------------------------------------------------------------

    let sam_bytes_rx = parse_bytes(
        stdout,
        &config,
        StreamDataType::JustBytes,
    )
        .await
        .map_err(PipelineError::Other)?;

    let (mut sam_branches, sam_fanout_done) =
        fanout_to_channels(
            ReceiverStream::new(sam_bytes_rx),
            2,
            "AMR KMA SAM bytes",
            &config,
            StreamDataType::JustBytes,
        )
            .await?;

    // ------------------------------------------------------------
    // Branch 1: native AMR parser.
    //
    // KMA stdout was chunked arbitrarily by parse_bytes(), so convert
    // the byte stream back into complete SAM lines before handing it
    // to consume_kma_sam().
    //
    // bytes_to_lines() preserves the newline in each emitted line.
    // consume_kma_sam() therefore needs to trim the line ending before
    // parsing (see small change below).
    // ------------------------------------------------------------

    let native_sam_bytes_rx = sam_branches.remove(0);

    let (sam_lines_rx, sam_lines_task) =
        bytes_to_lines(
            native_sam_bytes_rx,
            config.base_buffer_size,
        )
            .await
            .map_err(PipelineError::Other)?;

    let sam_stream =
        ReceiverStream::new(sam_lines_rx);

    // ------------------------------------------------------------
    // Branch 2: exact raw SAM bytes for samtools fixmate.
    //
    // Do NOT parse, rewrite, materialize, or otherwise modify this
    // branch. stream_to_cmd() can feed these ParseOutput::Bytes
    // directly into samtools stdin.
    // ------------------------------------------------------------

    let samtools_stream =
        sam_branches.remove(0);

    // ------------------------------------------------------------
    // Own KMA-side producer/process completion in one task.
    //
    // kma_temp_dir is deliberately NOT moved into this task. The
    // returned RunKmaHandle owns it so KMA sidecars survive process
    // completion.
    //
    // Neither SAM branch is consumed here. The caller must install
    // both downstream consumers immediately:
    //
    //   sam_stream      -> consume_kma_sam()
    //   samtools_stream -> samtools fixmate
    //
    // This ensures KMA stdout cannot backpressure indefinitely.
    // ------------------------------------------------------------

    let completion_task = tokio::spawn(async move {
        input_adapter_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunKMA input adapter task panicked: {}",
                    e
                )
            })??;

        stdin_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunKMA stdin task panicked: {}",
                    e
                )
            })??;

        stderr_task
            .await
            .map_err(|e| {
                anyhow!(
                    "RunKMA stderr task panicked: {}",
                    e
                )
            })??;

        info!(
            "[AMR] RunKMA process completed ({})",
            if paired {
                "paired-end"
            } else {
                "single-end"
            }
        );

        Ok::<(), anyhow::Error>(())
    });

    Ok(RunKmaHandle {
        sam_stream,
        samtools_stream,
        sam_fanout_done,
        sam_lines_task,
        completion_task,
        output_temp_dir: kma_temp_dir,
        output_prefix: kma_output_prefix,
    })
}

fn parse_kma_sq_header(line: &str) -> anyhow::Result<Option<KmaReference>> {
    if !line.starts_with("@SQ\t") {
        return Ok(None);
    }

    let mut name: Option<&str> = None;
    let mut length: Option<u64> = None;

    for field in line.split('\t').skip(1) {
        if let Some(value) = field.strip_prefix("SN:") {
            name = Some(value);
        } else if let Some(value) = field.strip_prefix("LN:") {
            length = Some(
                value
                    .parse::<u64>()
                    .map_err(|e| {
                        anyhow!(
                            "Invalid @SQ LN value '{}' in KMA SAM header '{}': {}",
                            value,
                            line,
                            e
                        )
                    })?,
            );
        }
    }

    let name = name.ok_or_else(|| {
        anyhow!("KMA @SQ header missing SN field: {}", line)
    })?;

    let length = length.ok_or_else(|| {
        anyhow!("KMA @SQ header missing LN field: {}", line)
    })?;

    Ok(Some(KmaReference {
        name: name.to_string(),
        length,
    }))
}

fn is_valid_kma_alignment(record: &crate::utils::sambam::SamRecordRef<'_>) -> bool {
    !record.flag_is_unmapped()
        && record.rname != "*"
        && record.pos != 0
        && record.cigar != "*"
}

async fn consume_kma_sam(
    mut sam_stream: ReceiverStream<ParseOutput>,
) -> anyhow::Result<KmaSamResults> {
    let mut references = HashMap::<String, KmaReference>::new();
    let mut alignments = Vec::<KmaAlignment>::new();

    let mut total_records = 0u64;
    let mut invalid_alignment_records = 0u64;

    while let Some(item) = sam_stream.next().await {
        let ParseOutput::Bytes(line_bytes) = item else {
            return Err(anyhow!(
                "RunKMA expected line-oriented Bytes output"
            ));
        };

        let line = std::str::from_utf8(&line_bytes)
            .map_err(|e| {
                anyhow!(
            "RunKMA emitted non-UTF8 SAM line: {}",
            e
        )
            })?
            .trim_end_matches(['\r', '\n']);

        if line.starts_with('@') {
            if let Some(reference) = parse_kma_sq_header(line)? {
                references.insert(
                    reference.name.clone(),
                    reference,
                );
            }

            continue;
        }

        total_records += 1;

        let record = parse_sam_record(line)?;

        if !is_valid_kma_alignment(&record) {
            invalid_alignment_records += 1;
            continue;
        }

        alignments.push(KmaAlignment {
            qname: record.qname.to_string(),
            flag: record.flag,
            rname: record.rname.to_string(),
            pos: record.pos,
            mapq: record.mapq,
            cigar: record.cigar.to_string(),
            rnext: record.rnext.to_string(),
            seq: record.seq.to_string(),
        });
    }

    info!(
        "[AMR] RunKMA SAM parsed: {} records, {} valid mapped alignments, \
         {} non-alignment records, {} reference sequences",
        total_records,
        alignments.len(),
        invalid_alignment_records,
        references.len(),
    );

    Ok(KmaSamResults {
        references,
        alignments,
        total_records,
        invalid_alignment_records,
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
            // Split host-filter outputs by ownership.
            //
            // Raw AMR consumes:
            //
            //   * pre_dedup_interleaved
            //   * subsampled_interleaved
            //   * duplicate_clusters
            //
            // run_host_filter() also creates the normal short-read
            // deinterleaved non_host_r1/non_host_r2 outputs. AMR does not
            // use those outputs, but they MUST still be drained. Otherwise
            // their bounded channels can fill and backpressure the shared
            // subsampled fanout.
            // ------------------------------------------------------------

            let HostFilterOutputs {
                non_host_r1,
                non_host_r2,
                subsampled_interleaved,
                duplicate_clusters,
                pre_dedup_interleaved,
            } = outputs;

            // ------------------------------------------------------------
            // Drain the unused normal short-read non-host outputs.
            //
            // Install these consumers immediately, before awaiting any
            // host-filter completion work.
            // ------------------------------------------------------------

            let mut unused_non_host_r1 = non_host_r1;

            cleanup_tasks.push(tokio::spawn(async move {
                while unused_non_host_r1.next().await.is_some() {}

                Ok::<(), anyhow::Error>(())
            }));

            let has_non_host_r2 = non_host_r2.is_some();

            if let Some(mut unused_non_host_r2) = non_host_r2 {
                cleanup_tasks.push(tokio::spawn(async move {
                    while unused_non_host_r2.next().await.is_some() {}

                    Ok::<(), anyhow::Error>(())
                }));
            }

            // ------------------------------------------------------------
            // Validate and take the AMR-specific retained streams.
            // ------------------------------------------------------------

            let pre_dedup_stream =
                pre_dedup_interleaved.ok_or_else(|| {
                    PipelineError::Other(anyhow!(
                        "AMR raw RunRedup missing retained interleaved pre-dedup stream"
                    ))
                })?;

            let subsampled_stream =
                subsampled_interleaved.ok_or_else(|| {
                    PipelineError::Other(anyhow!(
                        "AMR raw RunRedup missing retained interleaved subsampled stream"
                    ))
                })?;

            info!(
                "[AMR] host-filter streaming outputs ready: \
                 pre_dedup_interleaved=true, \
                 subsampled_interleaved=true, \
                 non_host_r1=draining, non_host_r2={}, \
                 duplicate_clusters=pending",
                has_non_host_r2,
            );

            // ------------------------------------------------------------
            // Install RunRedup consumers immediately.
            //
            // Both retained streams are structured interleaved FASTQ:
            //
            //   pre_dedup_stream
            //       exact host-filtered records immediately before dedup
            //
            //   subsampled_stream
            //       deduplicated + subsampled representatives
            //
            // duplicate_clusters completes only after dedup reaches EOF.
            // run_redup() installs consumers for both live streams before
            // awaiting that task.
            // ------------------------------------------------------------

            let redup = run_redup(
                config.clone(),
                pre_dedup_stream,
                subsampled_stream,
                paired,
                duplicate_clusters,
            )
                .await?;

            // ------------------------------------------------------------
            // Raw WDL topology:
            //
            //     RunRedupRaw -> RunSpades
            //
            // Install the SPAdes consumer immediately so every bounded
            // RunRedup output has a downstream reader before any
            // host-filter completion work is awaited.
            // ------------------------------------------------------------


            let mut spades_streams = Vec::new();
            let mut kma_streams = Vec::new();
            let mut redup_fanout_tasks = Vec::new();

            for (idx, redup_stream) in
                redup.redups_streams.into_iter().enumerate()
            {
                let (mut branches, done_rx) = fanout_to_channels(
                    redup_stream,
                    2,
                    &format!("AMR redups part {}", idx + 1),
                    &config,
                    StreamDataType::IlluminaFastq,
                )
                    .await?;

                spades_streams.push(
                    ReceiverStream::new(branches.remove(0))
                );

                kma_streams.push(
                    ReceiverStream::new(branches.remove(0))
                );

                redup_fanout_tasks.push(done_rx);
            }



            let spades = run_spades_amr(
                config.clone(),
                spades_streams,
            )
                .await?;

            let kma = run_kma_amr(
                config.clone(),
                kma_streams,
                PathBuf::from(
                    "/data/refs/card/rgi_db/localDB/bwt/card_reference/kma"
                ),
            )
                .await?;
            let RunKmaHandle {
                sam_stream: kma_sam_stream,
                samtools_stream: kma_samtools_stream,
                sam_fanout_done,
                sam_lines_task,
                completion_task: kma_completion_task,
                output_temp_dir: kma_output_temp_dir,
                output_prefix: kma_output_prefix,
            } = kma;
            // Keep KMA sidecars alive through downstream AMR post-processing.
            final_temp_dirs.push(kma_output_temp_dir);

            let redup_completion_task = redup.completion_task;

            // ------------------------------------------------------------
            // Parse KMA SAM stdout directly into the native AMR
            // representation.
            //
            // This consumer is installed immediately after RunKMA returns so
            // KMA stdout can never block on the bounded line channel.
            // ------------------------------------------------------------

            let kma_sam_task = tokio::spawn(
                consume_kma_sam(kma_sam_stream)
            );

            // ------------------------------------------------------------
            // Stream raw KMA SAM directly through the same samtools fixmate
            // transformation used by RGI.
            //
            // RGI:
            //
            //   samtools fixmate
            //       --verbosity 0
            //       -z all
            //       --threads N
            //       -m
            //       input.sam
            //       -O bam
            //       output.bam
            //
            // Here both file operands are "-" so the transformation remains
            // fully streaming:
            //
            //   KMA stdout -> samtools fixmate stdin -> BAM stdout
            // ------------------------------------------------------------

            let fixmate_threads =
                config.thread_allocation(SAMTOOLS_TAG, Some("fixmate"));

            let fixmate_config =
                SamtoolsConfig::new(SamtoolsSubcommand::Fixmate)
                    .option("--verbosity", "0")
                    .option("-z", "all")
                    .option("--threads", fixmate_threads.to_string())
                    .flag("-m")
                    .option("-O", "bam")
                    .positional("-")
                    .positional("-");

            let fixmate_args = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&fixmate_config),
            )
                .map_err(PipelineError::Other)?;

            info!(
    "[AMR] RunRgiBwtKma samtools fixmate args: {:?}",
    fixmate_args
);

            let (
                fixmate_child,
                fixmate_stdin_task,
                fixmate_stderr_task,
            ) = stream_to_cmd(
                config.clone(),
                kma_samtools_stream,
                SAMTOOLS_TAG,
                fixmate_args,
                StreamDataType::JustBytes,
                config.args.verbose,
                None,
            )
                .await
                .map_err(PipelineError::Other)?;

            // ------------------------------------------------------------
            // Stream fixmate BAM stdout directly into samtools sort.
            //
            // RGI:
            //
            //   samtools sort
            //       --threads N
            //       -T <prefix>.sorted
            //       -o <sorted.bam>
            //       <unsorted.bam>
            //
            // We preserve the sort operation but keep the pipeline streaming:
            //
            //   fixmate stdout -> samtools sort stdin -> sorted BAM stdout
            //
            // No intermediate unsorted BAM is materialized.
            // ------------------------------------------------------------

            let fixmate_stdout = {
                let mut child_guard = fixmate_child.lock().await;

                child_guard
                    .stdout
                    .take()
                    .ok_or_else(|| {
                        PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools fixmate stdout not available"
            ))
                    })?
            };

            let fixmate_bam_rx = parse_bytes(
                fixmate_stdout,
                &config,
                StreamDataType::JustBytes,
            )
                .await
                .map_err(PipelineError::Other)?;

            let sort_temp_dir = choose_temp_dir(
                config.input_size,
                &config.ram_temp_dir,
                &config.args.nvme_scratch,
                4,
                false,
            )
                .await?;

            let sort_temp_prefix =
                sort_temp_dir.path().join("sr_amr_report.sorted");

            let sort_threads =
                config.thread_allocation(SAMTOOLS_TAG, Some("sort"));

            let sort_config =
                SamtoolsConfig::new(SamtoolsSubcommand::Sort)
                    .option("--threads", sort_threads.to_string())
                    .option(
                        "-T",
                        sort_temp_prefix.to_string_lossy().to_string(),
                    )
                    .option("-O", "bam")
                    .positional("-");

            let sort_args = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&sort_config),
            )
                .map_err(PipelineError::Other)?;

            info!(
    "[AMR] RunRgiBwtKma samtools sort args: {:?}",
    sort_args
);

            let (
                sort_child,
                sort_stdin_task,
                sort_stderr_task,
            ) = stream_to_cmd(
                config.clone(),
                fixmate_bam_rx,
                SAMTOOLS_TAG,
                sort_args,
                StreamDataType::JustBytes,
                config.args.verbose,
                None,
            )
                .await
                .map_err(PipelineError::Other)?;



            // ------------------------------------------------------------
            // Stream the coordinate-sorted BAM through the same predicate-free
            // bamtools filter transformation used by RGI, then re-emit it with
            // samtools before materializing the final length_100 BAM:
            //
            //   sort_stdout
            //       -> parse_bytes
            //       -> stream_to_cmd(BAMTOOLS_TAG, ["filter"])
            //       -> parse_bytes
            //       -> stream_to_cmd(SAMTOOLS_TAG, ["view", "-b", "-"])
            //       -> parse_bytes
            //       -> write_byte_stream_to_file(length_100.bam)
            //
            // BamTools 2.5.2 accepts BAM on stdin and emits the same decoded
            // alignments as its file-based filter invocation, but its stdout BAM
            // omits the BGZF EOF marker. `samtools view -b -` re-emits a properly
            // finalized BAM while preserving the alignment stream. This lets us
            // avoid both RGI's sorted.temp.bam materialization and its unused
            // intermediate BAI; the only BAM file boundary is the durable
            // sr_amr_report.sorted.length_100.bam required downstream.
            // ------------------------------------------------------------

            let sort_stdout = {
                let mut child_guard = sort_child.lock().await;

                child_guard
                    .stdout
                    .take()
                    .ok_or_else(|| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools sort stdout not available"
                        ))
                    })?
            };

            let sorted_bam_rx = parse_bytes(
                sort_stdout,
                &config,
                StreamDataType::JustBytes,
            )
                .await
                .map_err(PipelineError::Other)?;

            let bamtools_config =
                BamtoolsConfig::new(BamtoolsSubcommand::Filter);

            let bamtools_args = generate_cli(
                BAMTOOLS_TAG,
                &config,
                Some(&bamtools_config),
            )
                .map_err(PipelineError::Other)?;

            info!(
                "[AMR] RunRgiBwtKma bamtools filter args: {:?}",
                bamtools_args
            );

            let (
                bamtools_child,
                bamtools_stdin_task,
                bamtools_stderr_task,
            ) = stream_to_cmd(
                config.clone(),
                sorted_bam_rx,
                BAMTOOLS_TAG,
                bamtools_args,
                StreamDataType::JustBytes,
                config.args.verbose,
                None,
            )
                .await
                .map_err(PipelineError::Other)?;

            let bamtools_stdout = {
                let mut child_guard = bamtools_child.lock().await;

                child_guard
                    .stdout
                    .take()
                    .ok_or_else(|| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma bamtools filter stdout not available"
                        ))
                    })?
            };

            let bamtools_bam_rx = parse_bytes(
                bamtools_stdout,
                &config,
                StreamDataType::JustBytes,
            )
                .await
                .map_err(PipelineError::Other)?;

            let finalize_config =
                SamtoolsConfig::new(SamtoolsSubcommand::View)
                    .flag("-b")
                    .positional("-");

            let finalize_args = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&finalize_config),
            )
                .map_err(PipelineError::Other)?;

            info!(
                "[AMR] RunRgiBwtKma samtools BAM finalizer args: {:?}",
                finalize_args
            );

            let (
                finalize_child,
                finalize_stdin_task,
                finalize_stderr_task,
            ) = stream_to_cmd(
                config.clone(),
                bamtools_bam_rx,
                SAMTOOLS_TAG,
                finalize_args,
                StreamDataType::JustBytes,
                config.args.verbose,
                None,
            )
                .await
                .map_err(PipelineError::Other)?;

            let finalize_stdout = {
                let mut child_guard = finalize_child.lock().await;

                child_guard
                    .stdout
                    .take()
                    .ok_or_else(|| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools BAM finalizer stdout not available"
                        ))
                    })?
            };

            let final_bam_rx = parse_bytes(
                finalize_stdout,
                &config,
                StreamDataType::JustBytes,
            )
                .await
                .map_err(PipelineError::Other)?;

            let sorted_bam_path =
                out_dir.join("sr_amr_report.sorted.length_100.bam");

            let sorted_bam_write_task = write_byte_stream_to_file(
                &sorted_bam_path,
                ReceiverStream::new(final_bam_rx),
                config.clone(),
                StreamDataType::JustBytes,
                "RunRgiBwtKma length_100 BAM",
                false,
            )
                .await
                .map_err(PipelineError::Other)?;

            // Keep samtools sort scratch alive until sort has completed.
            final_temp_dirs.push(sort_temp_dir);


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

            // ------------------------------------------------------------
            // Complete the streaming
            // KMA -> fixmate -> sort -> bamtools -> samtools view -> BAM chain.
            //
            // All stdout consumers have already been installed above, so it
            // is now safe to await producer/process completion. The final BAM
            // is not considered ready until every stdin/process task has
            // completed successfully and samtools view stdout has been fully
            // materialized.
            // ------------------------------------------------------------

            let rgi_bwt_bam_fut = async move {
                let (
                    fixmate_stdin_result,
                    sort_stdin_result,
                    bamtools_stdin_result,
                    finalize_stdin_result,
                    sorted_bam_write_result,
                ) = tokio::join!(
                    fixmate_stdin_task,
                    sort_stdin_task,
                    bamtools_stdin_task,
                    finalize_stdin_task,
                    sorted_bam_write_task,
                );

                fixmate_stdin_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools fixmate stdin task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                sort_stdin_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools sort stdin task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                bamtools_stdin_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma bamtools filter stdin task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                finalize_stdin_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools BAM finalizer stdin task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                sorted_bam_write_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma length_100 BAM writer task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                let (
                    fixmate_stderr_result,
                    sort_stderr_result,
                    bamtools_stderr_result,
                    finalize_stderr_result,
                ) = tokio::join!(
                    fixmate_stderr_task,
                    sort_stderr_task,
                    bamtools_stderr_task,
                    finalize_stderr_task,
                );

                fixmate_stderr_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools fixmate stderr task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                sort_stderr_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools sort stderr task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                bamtools_stderr_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma bamtools filter stderr task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                finalize_stderr_result
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools BAM finalizer stderr task panicked: {}",
                            e
                        ))
                    })?
                    .map_err(PipelineError::Other)?;

                info!(
                    "[AMR] RunRgiBwtKma length_100 BAM ready: {}",
                    sorted_bam_path.display()
                );

                Ok::<PathBuf, PipelineError>(sorted_bam_path)
            };

            let (
                validation_and_qc,
                host_counts,
                (),
                dedup_and_subsample,
                sorted_bam_path,
            ) = tokio::try_join!(
                validation_and_qc_fut,
                host_counts_fut,
                kallisto_fut,
                dedup_and_subsample_fut,
                rgi_bwt_bam_fut,
            )?;

            // ------------------------------------------------------------
            // Index the completed coordinate-sorted AMR BAM.
            //
            // rgi_bwt_bam_fut does not return until the final
            // sr_amr_report.sorted.length_100.bam writer has reached EOF.
            // samtools index therefore starts only after the durable BAM
            // file boundary is complete.
            // ------------------------------------------------------------

            let sorted_bai_path =
                PathBuf::from(format!("{}.bai", sorted_bam_path.display()));

            let index_threads =
                config.thread_allocation(SAMTOOLS_TAG, Some("index"));

            let index_config =
                SamtoolsConfig::new(SamtoolsSubcommand::Index)
                    .option("-@", index_threads.to_string())
                    .flag("-b")
                    .positional(
                        sorted_bam_path
                            .to_string_lossy()
                            .to_string(),
                    )
                    .positional(
                        sorted_bai_path
                            .to_string_lossy()
                            .to_string(),
                    );

            let index_args = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&index_config),
            )
                .map_err(PipelineError::Other)?;

            info!(
    "[AMR] RunRgiBwtKma samtools index args: {:?}",
    index_args
);

            let (mut index_child, index_stderr_task) =
                spawn_cmd(
                    config.clone(),
                    SAMTOOLS_TAG,
                    index_args,
                    config.args.verbose,
                    None,
                )
                    .await
                    .map_err(PipelineError::Other)?;

            let index_status =
                index_child
                    .wait()
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                "RunRgiBwtKma samtools index wait failed: {}",
                e
            ))
                    })?;

            index_stderr_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunRgiBwtKma samtools index stderr task panicked: {}",
            e
        ))
                })?
                .map_err(PipelineError::Other)?;

            if !index_status.success() {
                return Err(PipelineError::Other(anyhow!(
        "RunRgiBwtKma samtools index failed with exit status {}",
        index_status
    )));
            }

            if !sorted_bai_path.is_file() {
                return Err(PipelineError::Other(anyhow!(
        "RunRgiBwtKma samtools index completed but BAI was not created: {}",
        sorted_bai_path.display()
    )));
            }

            info!(
    "[AMR] RunRgiBwtKma BAM index ready: {}",
    sorted_bai_path.display()
);

            // ------------------------------------------------------------
            // Pull per-reference alignment counts with samtools idxstats.
            //
            // RGI materializes this stdout as sr_amr_report.temp.txt and
            // get_reads_count() later keeps only rows with mapped > 0.
            // There is no external-tool requirement for that text-file
            // boundary, so consume idxstats stdout directly and preserve the
            // same mapped-reference semantics in memory.
            // ------------------------------------------------------------

            let idxstats_config =
                SamtoolsConfig::new(SamtoolsSubcommand::Idxstats)
                    .positional(
                        sorted_bam_path
                            .to_string_lossy()
                            .to_string(),
                    );

            let idxstats_args = generate_cli(
                SAMTOOLS_TAG,
                &config,
                Some(&idxstats_config),
            )
                .map_err(PipelineError::Other)?;

            info!(
                "[AMR] RunRgiBwtKma samtools idxstats args: {:?}",
                idxstats_args
            );

            let (mut idxstats_child, idxstats_stderr_task) =
                spawn_cmd(
                    config.clone(),
                    SAMTOOLS_TAG,
                    idxstats_args,
                    config.args.verbose,
                    None,
                )
                    .await
                    .map_err(PipelineError::Other)?;

            let mut idxstats_stdout = parse_child_output(
                &mut idxstats_child,
                ChildStream::Stdout,
                ParseMode::Lines,
                &config,
            )
                .await
                .map_err(PipelineError::Other)?;

            let mut rgi_read_counts = HashMap::<String, RgiReadCounts>::new();

            while let Some(item) = idxstats_stdout.recv().await {
                let ParseOutput::Bytes(line_bytes) = item else {
                    return Err(PipelineError::Other(anyhow!(
                        "RunRgiBwtKma samtools idxstats expected line-oriented Bytes output"
                    )));
                };

                let line = std::str::from_utf8(&line_bytes)
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma samtools idxstats emitted non-UTF8 output: {}",
                            e
                        ))
                    })?
                    .trim_end_matches(['\r', '\n']);

                if line.is_empty() {
                    continue;
                }

                let fields: Vec<&str> = line.split('\t').collect();

                if fields.len() != 4 {
                    return Err(PipelineError::Other(anyhow!(
                        "RunRgiBwtKma samtools idxstats expected 4 tab-separated fields, got {}: {}",
                        fields.len(),
                        line
                    )));
                }

                let reference_length = fields[1]
                    .parse::<u64>()
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma invalid idxstats reference length '{}': {}",
                            fields[1],
                            e
                        ))
                    })?;

                let mapped = fields[2]
                    .parse::<u64>()
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma invalid idxstats mapped count '{}': {}",
                            fields[2],
                            e
                        ))
                    })?;

                let unmapped = fields[3]
                    .parse::<u64>()
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                            "RunRgiBwtKma invalid idxstats unmapped count '{}': {}",
                            fields[3],
                            e
                        ))
                    })?;

                // Match RGI get_reads_count(): retain only references with at
                // least one mapped read. The '*' idxstats row therefore drops
                // out naturally because its mapped count is zero.
                if mapped > 0 {
                    rgi_read_counts.insert(
                        fields[0].to_string(),
                        RgiReadCounts {
                            reference_length,
                            mapped,
                            unmapped,
                            all: mapped + unmapped,
                        },
                    );
                }
            }

            let idxstats_status = idxstats_child
                .wait()
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "RunRgiBwtKma samtools idxstats wait failed: {}",
                        e
                    ))
                })?;

            idxstats_stderr_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
                        "RunRgiBwtKma samtools idxstats stderr task panicked: {}",
                        e
                    ))
                })?
                .map_err(PipelineError::Other)?;

            if !idxstats_status.success() {
                return Err(PipelineError::Other(anyhow!(
                    "RunRgiBwtKma samtools idxstats failed with exit status {}",
                    idxstats_status
                )));
            }

            info!(
                "[AMR] RunRgiBwtKma idxstats parsed: {} references with mapped reads",
                rgi_read_counts.len()
            );

            for (reference, counts) in &rgi_read_counts {
                info!(
                    "[AMR] RGI idxstats: reference={} length={} mapped={} unmapped={} all={}",
                    reference,
                    counts.reference_length,
                    counts.mapped,
                    counts.unmapped,
                    counts.all,
                );
            }


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
            // Complete RunRedup and both downstream consumers.
            // ------------------------------------------------------------

            redup_completion_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunRedup task panicked: {}",
            e
        ))
                })?
                .map_err(PipelineError::Other)?;

            // Every RunRedup item must have reached both SPAdes and KMA.
            for (idx, done_rx) in redup_fanout_tasks.into_iter().enumerate() {
                done_rx
                    .await
                    .map_err(|e| {
                        PipelineError::Other(anyhow!(
                "RunRedup fanout {} completion receiver dropped: {}",
                idx + 1,
                e
            ))
                    })?
                    .map_err(PipelineError::Other)?;
            }

            // KMA stdout must be drained concurrently with KMA execution.
            // At this point wait for both sides to finish.
            let kma_sam_results = kma_sam_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunKMA SAM parser task panicked: {}",
            e
        ))
                })?
                .map_err(PipelineError::Other)?;

            sam_lines_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunKMA SAM line adapter task panicked: {}",
            e
        ))
                })?
                .map_err(PipelineError::Other)?;

            sam_fanout_done
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunKMA SAM fanout completion receiver dropped: {}",
            e
        ))
                })?
                .map_err(PipelineError::Other)?;



            info!(
        "[AMR] RunKMA native results: records={}, mapped={}, references={}",
        kma_sam_results.total_records,
        kma_sam_results.alignments.len(),
        kma_sam_results.references.len(),
    );

            for alignment in &kma_sam_results.alignments {
                info!(
        "[AMR] KMA mapped: qname={} flag={} rname={} pos={} mapq={} cigar={}",
        alignment.qname,
        alignment.flag,
        alignment.rname,
        alignment.pos,
        alignment.mapq,
        alignment.cigar,
    );
            }



            kma_completion_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunKMA completion task panicked: {}",
            e
        ))
                })?
                .map_err(PipelineError::Other)?;

            spades.completion_task
                .await
                .map_err(|e| {
                    PipelineError::Other(anyhow!(
            "RunSpades task panicked: {}",
            e
        ))
                })?
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
