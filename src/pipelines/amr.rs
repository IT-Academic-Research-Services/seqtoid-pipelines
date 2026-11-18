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
                          SPADES_TAG, RgiSubcommand, RGI_TAG};
use crate::pipelines::short_read_mngs::{
    kallisto_results,
    run_host_filter,
    HostFilterCompletion,
};
use crate::utils::command::spades::SpadesConfig;
use crate::utils::command::rgi::RgiConfig;
use crate::utils::command::generate_cli;
use crate::utils::fastx::{parse_byte_stream_to_fastq, read_fasta, read_fastq, write_fasta_stream_to_file, SequenceRecord};
use crate::utils::file::{resolve_existing_input_path, validate_file_inputs, write_count_file};
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

#[derive(Debug)]
struct AmrPreparedInputs {
    /// Sample basename returned by validate_file_inputs for the selected input mode.
    sample_base: String,
    /// WDL-equivalent `non_host_reads_fa`: durable RunRedup FASTA files.
    non_host_reads_fa: Vec<PathBuf>,
    /// WDL-equivalent `contigs_fa`: raw SPAdes contigs or filtered-sample contigs.
    contigs_fa: PathBuf,
}

struct RunRedupHandle {
    /// One live FASTA stream per input end, in WDL `redups_1.fa`,
    /// `redups_2.fa` order.
    redups_streams: Vec<ReceiverStream<ParseOutput>>,

    /// Completes after every per-end producer has reached EOF.
    /// Callers must install consumers for `redups_streams` before awaiting it.
    completion_task: JoinHandle<Result<(), anyhow::Error>>,
}

/// Resolve the prepared localDB from the existing native-KMA prefix argument.
/// The prefix itself need not be a file (KMA appends index suffixes).
async fn resolve_amr_rgi_root(config: &RunConfig) -> Result<PathBuf, PipelineError> {
    let prefix = config.args.amr_kma_db.as_deref().ok_or_else(|| {
        PipelineError::InvalidConfig(
            "--amr-kma-db is required; expected <root>/localDB/bwt/card_reference/kma".into(),
        )
    })?;
    let prefix = PathBuf::from(prefix);
    let reference_dir = prefix.parent().ok_or_else(|| {
        PipelineError::InvalidConfig("Invalid --amr-kma-db prefix".into())
    })?;
    let reference_dir = tokio::fs::canonicalize(reference_dir).await
        .map_err(|e| PipelineError::IOError(format!(
            "Cannot resolve --amr-kma-db parent {}: {}", reference_dir.display(), e,
        )))?;
    let bwt_dir = reference_dir.parent();
    let local_db = bwt_dir.and_then(|p| p.parent());
    if prefix.file_name() != Some(std::ffi::OsStr::new("kma"))
        || reference_dir.file_name() != Some(std::ffi::OsStr::new("card_reference"))
        || bwt_dir.and_then(|p| p.file_name()) != Some(std::ffi::OsStr::new("bwt"))
        || local_db.and_then(|p| p.file_name()) != Some(std::ffi::OsStr::new("localDB"))
    {
        return Err(PipelineError::InvalidConfig(
            "--amr-kma-db must identify <root>/localDB/bwt/card_reference/kma; RGI selects that localDB, not an arbitrary KMA prefix".into(),
        ));
    }
    let local_db = local_db.unwrap(); // checked above
    for name in ["card.json", "card_reference.fasta", "loaded_databases.json",
        "61mer_database.json", "amr_61mer.txt", "index-for-model-sequences.txt"] {
        let path = local_db.join(name);
        if !tokio::fs::metadata(&path).await
            .map_err(|e| PipelineError::IOError(format!("RGI reference {}: {}", path.display(), e)))?
            .is_file()
        {
            return Err(PipelineError::InvalidConfig(format!(
                "RGI reference is not a regular file: {}", path.display(),
            )));
        }
    }
    // RGI reads installed localDB/card.json, not the original source path.
    // If supplied, reject a different source rather than silently ignoring it.
    if let Some(source) = config.args.amr_card_json.as_deref() {
        let source_bytes = tokio::fs::read(source).await
            .map_err(|e| PipelineError::IOError(format!("CARD source {}: {}", source, e)))?;
        let installed_bytes = tokio::fs::read(local_db.join("card.json")).await
            .map_err(|e| PipelineError::IOError(e.to_string()))?;
        let source_json: serde_json::Value = serde_json::from_slice(&source_bytes)
            .map_err(|e| PipelineError::Other(anyhow!("Invalid source CARD JSON: {}", e)))?;
        let installed_json: serde_json::Value = serde_json::from_slice(&installed_bytes)
            .map_err(|e| PipelineError::Other(anyhow!("Invalid installed CARD JSON: {}", e)))?;
        if source_json != installed_json {
            return Err(PipelineError::InvalidConfig(format!(
                "--amr-card-json differs from {}", local_db.join("card.json").display(),
            )));
        }
    }
    Ok(local_db.parent().unwrap().to_path_buf())
}

/// AMR-only launcher:
/// stderr goes directly to a durable file, so it cannot fill a pipe.
/// The returned completed task preserves the callers' existing join contract.
async fn spawn_amr_rgi(
    root: &std::path::Path,
    args: Vec<String>,
    stderr_log: PathBuf,
) -> anyhow::Result<(tokio::process::Child, JoinHandle<anyhow::Result<()>>)> {
    let stderr_file = tokio::fs::File::create(&stderr_log).await?.into_std().await;
    let child = tokio::process::Command::new(RGI_TAG)
        .args(args)
        .current_dir(root)
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::from(stderr_file))
        .kill_on_drop(true)
        .spawn()?;
    Ok((child, tokio::spawn(async { Ok(()) })))
}

#[derive(Debug)]
struct RgiBwtKmaOutputs {
    /// Owns the returned scratch paths. Retain until all downstream
    /// readers and artifact-packaging tasks have completed.
    _work_dir: TempDir,

    output_kma: Vec<PathBuf>,

    kma_amr_results: PathBuf,
    artifacts_mapping_stats: PathBuf,
    gene_mapping_data: PathBuf,
    overall_mapping_stats: PathBuf,
    reference_mapping_stats: PathBuf,
    output_sorted_length_100: PathBuf,
    output_sorted_length_100_bai: PathBuf,
    kma_amr_results_json: PathBuf,
}

/// WDL RunRgiBwtKma: one RGI invocation with file-backed outputs.
///
/// The output prefix is placed inside choose_temp_dir. This does not
/// redirect every internal
/// temporary file RGI or its subprocesses might create.
///
/// Diagnostic logs are durable files in config.out_dir. Computational
/// outputs remain temporary and are owned by the returned TempDir.
async fn run_rgi_bwt_kma(
    config: Arc<RunConfig>,
    non_host_reads_fa: &[PathBuf],
    rgi_db_root: &std::path::Path,
) -> Result<RgiBwtKmaOutputs, PipelineError> {
    if !(1..=2).contains(&non_host_reads_fa.len()) {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRgiBwtKma requires one or two FASTA files; got {}",
            non_host_reads_fa.len(),
        )));
    }

    let mut reads = Vec::with_capacity(non_host_reads_fa.len());
    let mut input_bytes = 0u64;

    for path in non_host_reads_fa {
        let absolute = tokio::fs::canonicalize(path)
            .await
            .map_err(|e| {
                PipelineError::IOError(format!(
                    "RunRgiBwtKma cannot resolve input {}: {}",
                    path.display(),
                    e,
                ))
            })?;

        let metadata = tokio::fs::metadata(&absolute)
            .await
            .map_err(|e| {
                PipelineError::IOError(format!(
                    "RunRgiBwtKma cannot stat input {}: {}",
                    absolute.display(),
                    e,
                ))
            })?;

        if !metadata.is_file() {
            return Err(PipelineError::InvalidConfig(format!(
                "RunRgiBwtKma input is not a regular file: {}",
                absolute.display(),
            )));
        }

        input_bytes = input_bytes.saturating_add(metadata.len());
        reads.push(absolute);
    }

    // Heuristic, not a quota: SAM/BAM and KMA outputs may exceed the
    // input FASTA volume. choose_temp_dir applies headroom and fallback.
    let work_dir = crate::utils::file::choose_temp_dir(
        input_bytes.saturating_mul(4).max(1),
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    )
        .await?;

    let work_path = tokio::fs::canonicalize(work_dir.path())
        .await
        .map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiBwtKma cannot resolve scratch directory {}: {}",
                work_dir.path().display(),
                e,
            ))
        })?;

    tokio::fs::create_dir_all(&config.out_dir)
        .await
        .map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiBwtKma cannot create log directory {}: {}",
                config.out_dir.display(),
                e,
            ))
        })?;

    let stdout_log_path = config.out_dir.join("amr_rgi_bwt_stdout.log");
    let stderr_log_path = config.out_dir.join("amr_rgi_bwt_stderr.log");
    let prefix = work_path.join("sr_amr_report");

    let mut rgi_config = RgiConfig::new(RgiSubcommand::Bwt)
        .option("-1", reads[0].to_string_lossy().into_owned());

    if let Some(r2) = reads.get(1) {
        rgi_config =
            rgi_config.option("-2", r2.to_string_lossy().into_owned());
    }

    // Preserve the WDL command and defaults.
    rgi_config = rgi_config
        .flag("--local")
        .option("-a", "kma")
        .option("-o", prefix.to_string_lossy().into_owned())
        .flag("--clean");

    let args = generate_cli(RGI_TAG, &config, Some(&rgi_config))?;

    info!(
        "[RunRgiBwtKma] rgi {:?}; output directory={}; \
         stdout={}; stderr={}",
        args,
        work_path.display(),
        stdout_log_path.display(),
        stderr_log_path.display(),
    );

    let (mut child, stderr_task) = spawn_amr_rgi(
        rgi_db_root,
        args,
        stderr_log_path.clone(),
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!("RunRgiBwtKma spawn failed: {}", e),
        })?;

    // stderr is redirected to its log. Install a stdout consumer before
    // waiting so neither pipe can backpressure the running process.
    let stdout_rx = match parse_child_output(
        &mut child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        &config,
    )
        .await
    {
        Ok(rx) => rx,
        Err(e) => {
            // Do not leave the process running while work_dir is dropped.
            let kill_result = child.kill().await;
            let wait_result = child.wait().await;
            let stderr_result = stderr_task.await;

            if let Err(kill_error) = kill_result {
                warn!(
                    "[RunRgiBwtKma] cleanup kill failed: {}",
                    kill_error,
                );
            }
            if let Err(wait_error) = wait_result {
                warn!(
                    "[RunRgiBwtKma] cleanup wait failed: {}",
                    wait_error,
                );
            }
            if !matches!(&stderr_result, Ok(Ok(()))) {
                warn!(
                    "[RunRgiBwtKma] stderr capture during cleanup: {:?}",
                    stderr_result,
                );
            }

            return Err(PipelineError::ToolExecution {
                tool: RGI_TAG.to_string(),
                error: format!(
                    "RunRgiBwtKma stdout capture setup failed: {}; \
                     stderr log: {}",
                    e,
                    stderr_log_path.display(),
                ),
            });
        }
    };

    let stdout_task = tokio::spawn(stream_to_file(
        stdout_rx,
        stdout_log_path.clone(),
    ));

    // Join all three, rather than returning early and detaching a writer.
    let (status_result, stdout_result, stderr_result) =
        tokio::join!(child.wait(), stdout_task, stderr_task);

    // Report secondary capture failures even when the process itself
    // failed. The process error remains the primary failure.
    if !matches!(&stdout_result, Ok(Ok(()))) {
        warn!(
            "[RunRgiBwtKma] stdout capture failed: {:?}; log={}",
            stdout_result,
            stdout_log_path.display(),
        );
    }

    if !matches!(&stderr_result, Ok(Ok(()))) {
        warn!(
            "[RunRgiBwtKma] stderr capture failed: {:?}; log={}",
            stderr_result,
            stderr_log_path.display(),
        );
    }

    let status =
        status_result.map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiBwtKma wait failed: {}; stdout={}; stderr={}",
                e,
                stdout_log_path.display(),
                stderr_log_path.display(),
            ),
        })?;

    if !status.success() {
        return Err(PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiBwtKma exited with {}; stdout={}; stderr={}",
                status,
                stdout_log_path.display(),
                stderr_log_path.display(),
            ),
        });
    }

    stdout_result
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma stdout writer task failed: {}; log={}",
                e,
                stdout_log_path.display(),
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma stdout write failed: {}; log={}",
                e,
                stdout_log_path.display(),
            ))
        })?;

    stderr_result
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma stderr task failed: {}; log={}",
                e,
                stderr_log_path.display(),
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiBwtKma stderr capture failed: {}; log={}",
                e,
                stderr_log_path.display(),
            ))
        })?;

    let output = |suffix: &str| {
        work_path.join(format!("sr_amr_report.{}", suffix))
    };

    let required = [
        output("allele_mapping_data.txt"),
        output("artifacts_mapping_stats.txt"),
        output("gene_mapping_data.txt"),
        output("overall_mapping_stats.txt"),
        output("reference_mapping_stats.txt"),
        output("sorted.length_100.bam"),
        output("sorted.length_100.bam.bai"),
        output("allele_mapping_data.json"),
    ];

    for path in &required {
        let metadata = tokio::fs::metadata(path)
            .await
            .map_err(|e| {
                PipelineError::IOError(format!(
                    "RunRgiBwtKma missing/unreadable required output \
                     {}: {}; stderr={}",
                    path.display(),
                    e,
                    stderr_log_path.display(),
                ))
            })?;

        if !metadata.is_file() {
            return Err(PipelineError::IOError(format!(
                "RunRgiBwtKma output is not a regular file: {}",
                path.display(),
            )));
        }

        // Existence/type validation only. Do not impose a blanket
        // nonempty-file requirement on no-hit reports.
    }

    // WDL glob("sr_amr_report*"), scoped to this fresh output directory.
    let mut output_kma = Vec::new();
    let mut entries = tokio::fs::read_dir(&work_path)
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;

    while let Some(entry) = entries
        .next_entry()
        .await
        .map_err(|e| PipelineError::IOError(e.to_string()))?
    {
        if entry
            .file_name()
            .to_string_lossy()
            .starts_with("sr_amr_report")
            && entry
            .file_type()
            .await
            .map_err(|e| PipelineError::IOError(e.to_string()))?
            .is_file()
        {
            output_kma.push(entry.path());
        }
    }

    output_kma.sort();

    let [
    kma_amr_results,
    artifacts_mapping_stats,
    gene_mapping_data,
    overall_mapping_stats,
    reference_mapping_stats,
    output_sorted_length_100,
    output_sorted_length_100_bai,
    kma_amr_results_json,
    ] = required;

    Ok(RgiBwtKmaOutputs {
        _work_dir: work_dir,
        output_kma,
        kma_amr_results,
        artifacts_mapping_stats,
        gene_mapping_data,
        overall_mapping_stats,
        reference_mapping_stats,
        output_sorted_length_100,
        output_sorted_length_100_bai,
        kma_amr_results_json,
    })
}
#[derive(Debug)]
struct RgiMainOutputs {
    /// Owns returned scratch paths until downstream readers/packaging finish.
    _work_dir: TempDir,
    output_main: Vec<PathBuf>,
    output_json: PathBuf,
    main_amr_results: PathBuf,
}

/// WDL RunRgiMain: run RGI main on assembled contigs, or synthesize the
/// WDL-compatible empty outputs when SPAdes emitted its failure marker.
async fn run_rgi_main(
    config: Arc<RunConfig>,
    contigs_fa: &PathBuf,
    rgi_db_root: &std::path::Path,
) -> Result<RgiMainOutputs, PipelineError> {
    use tokio::io::{AsyncBufReadExt, BufReader};

    let contigs = tokio::fs::canonicalize(contigs_fa).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiMain cannot resolve input {}: {}", contigs_fa.display(), e,
        )))?;
    let metadata = tokio::fs::metadata(&contigs).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiMain cannot stat input {}: {}", contigs.display(), e,
        )))?;
    if !metadata.is_file() {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRgiMain input is not a regular file: {}", contigs.display(),
        )));
    }

    let work_dir = crate::utils::file::choose_temp_dir(
        metadata.len().saturating_mul(4).max(1),
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    ).await?;
    let work_path = tokio::fs::canonicalize(work_dir.path()).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiMain cannot resolve scratch directory {}: {}",
            work_dir.path().display(), e,
        )))?;

    tokio::fs::create_dir_all(&config.out_dir).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiMain cannot create log directory {}: {}",
            config.out_dir.display(), e,
        )))?;

    let prefix = work_path.join("contig_amr_report");
    let output_json = work_path.join("contig_amr_report.json");
    let main_amr_results = work_path.join("contig_amr_report.txt");

    let first_line = {
        let file = tokio::fs::File::open(&contigs).await
            .map_err(|e| PipelineError::IOError(format!(
                "RunRgiMain cannot open input {}: {}", contigs.display(), e,
            )))?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        reader.read_line(&mut line).await
            .map_err(|e| PipelineError::IOError(format!(
                "RunRgiMain cannot read input {}: {}", contigs.display(), e,
            )))?;
        line.trim_end_matches(['\r', '\n']).to_owned()
    };

    if first_line == ";ASSEMBLY FAILED" {
        // Exact bytes of the supplied WDL /tmp/empty-main-header.txt.
        const EMPTY_MAIN_HEADER: &str = "ORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_type\tSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug Class\tResistance Mechanism\tAMR Gene Family\tPredicted_DNA\tPredicted_Protein\tCARD_Protein_Sequence\tPercentage Length of Reference Sequence\tID\tModel_ID\tHit_Start\tHit_End\tNudged\tNote\r\n";
        tokio::try_join!(
            tokio::fs::write(&output_json, b"{}\n"),
            tokio::fs::write(&main_amr_results, EMPTY_MAIN_HEADER.as_bytes()),
        ).map_err(|e| PipelineError::IOError(format!(
            "RunRgiMain failed writing assembly-failure outputs: {}", e,
        )))?;
    } else {
        let stdout_log_path = config.out_dir.join("amr_rgi_main_stdout.log");
        let stderr_log_path = config.out_dir.join("amr_rgi_main_stderr.log");

        let rgi_config = RgiConfig::new(RgiSubcommand::Main)
            .flag("--local")
            .option("-i", contigs.to_string_lossy().into_owned())
            .option("-o", prefix.to_string_lossy().into_owned())
            .option("-t", "contig")
            .option("-a", "BLAST")
            .flag("--include_nudge")
            .flag("--clean");
        let args = generate_cli(RGI_TAG, &config, Some(&rgi_config))?;

        info!(
            "[RunRgiMain] rgi {:?}; output directory={}; stdout={}; stderr={}",
            args, work_path.display(), stdout_log_path.display(), stderr_log_path.display(),
        );

        let (mut child, stderr_task) = spawn_amr_rgi(
            rgi_db_root, args, stderr_log_path.clone(),
        ).await.map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!("RunRgiMain spawn failed: {}", e),
        })?;

        let stdout_rx = match parse_child_output(
            &mut child, ChildStream::Stdout, ParseMode::Bytes, &config,
        ).await {
            Ok(rx) => rx,
            Err(e) => {
                let _ = child.kill().await;
                let _ = child.wait().await;
                let _ = stderr_task.await;
                return Err(PipelineError::ToolExecution {
                    tool: RGI_TAG.to_string(),
                    error: format!(
                        "RunRgiMain stdout capture setup failed: {}; stderr log: {}",
                        e, stderr_log_path.display(),
                    ),
                });
            }
        };
        let stdout_task = tokio::spawn(stream_to_file(stdout_rx, stdout_log_path.clone()));
        let (status_result, stdout_result, stderr_result) =
            tokio::join!(child.wait(), stdout_task, stderr_task);

        let status = status_result.map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiMain wait failed: {}; stdout={}; stderr={}",
                e, stdout_log_path.display(), stderr_log_path.display(),
            ),
        })?;
        if !status.success() {
            return Err(PipelineError::ToolExecution {
                tool: RGI_TAG.to_string(),
                error: format!(
                    "RunRgiMain exited with {}; stdout={}; stderr={}",
                    status, stdout_log_path.display(), stderr_log_path.display(),
                ),
            });
        }
        stdout_result
            .map_err(|e| PipelineError::Other(anyhow!(
                "RunRgiMain stdout writer task failed: {}; log={}", e, stdout_log_path.display(),
            )))?
            .map_err(|e| PipelineError::Other(anyhow!(
                "RunRgiMain stdout write failed: {}; log={}", e, stdout_log_path.display(),
            )))?;
        stderr_result
            .map_err(|e| PipelineError::Other(anyhow!(
                "RunRgiMain stderr task failed: {}; log={}", e, stderr_log_path.display(),
            )))?
            .map_err(|e| PipelineError::Other(anyhow!(
                "RunRgiMain stderr capture failed: {}; log={}", e, stderr_log_path.display(),
            )))?;
    }

    for path in [&output_json, &main_amr_results] {
        let metadata = tokio::fs::metadata(path).await.map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiMain missing/unreadable required output {}: {}", path.display(), e,
            ))
        })?;
        if !metadata.is_file() {
            return Err(PipelineError::IOError(format!(
                "RunRgiMain output is not a regular file: {}", path.display(),
            )));
        }
    }

    // WDL glob("contig_amr_report*"), scoped to this fresh output directory.
    let mut output_main = Vec::new();
    let mut entries = tokio::fs::read_dir(&work_path).await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    while let Some(entry) = entries.next_entry().await
        .map_err(|e| PipelineError::IOError(e.to_string()))? {
        if entry.file_name().to_string_lossy().starts_with("contig_amr_report")
            && entry.file_type().await
            .map_err(|e| PipelineError::IOError(e.to_string()))?.is_file()
        {
            output_main.push(entry.path());
        }
    }
    output_main.sort();

    Ok(RgiMainOutputs {
        _work_dir: work_dir,
        output_main,
        output_json,
        main_amr_results,
    })
}

#[derive(Debug)]
struct GeneCoverageOutputs {
    _work_dir: TempDir,
    output_gene_coverage: PathBuf,
}

// Read quoted TSV fields without adding a new crate dependency. RGI ORF IDs
// and hit IDs are kept verbatim for the JSON lookup.
fn gene_coverage_tsv_records(input: &str) -> anyhow::Result<Vec<Vec<String>>> {
    let mut records = Vec::new();
    let mut record = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = input.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' if quoted => {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    field.push('"');
                } else {
                    quoted = false;
                }
            }
            '"' if field.is_empty() => quoted = true,
            '\t' if !quoted => record.push(std::mem::take(&mut field)),
            '\n' | '\r' if !quoted => {
                if c == '\r' && chars.peek() == Some(&'\n') { chars.next(); }
                record.push(std::mem::take(&mut field));
                if record.len() != 1 || !record[0].is_empty() {
                    records.push(std::mem::take(&mut record));
                } else { record.clear(); }
            }
            _ => field.push(c),
        }
    }
    if quoted { return Err(anyhow!("MakeGeneCoverage: unterminated quoted TSV field")); }
    if !field.is_empty() || !record.is_empty() {
        record.push(field);
        records.push(record);
    }
    Ok(records)
}

/// WDL MakeGeneCoverage: union reference coordinates per (ARO, Model_ID).
/// Reverse Python ranges include start and exclude end, hence [end+1,start+1).
fn gene_coverage_table(tsv: &str, json: &serde_json::Value) -> anyhow::Result<String> {
    use std::collections::BTreeMap;
    let records = gene_coverage_tsv_records(tsv)?;
    let header = records.first().ok_or_else(|| anyhow!("MakeGeneCoverage: missing TSV header"))?;
    let required = ["ORF_ID", "ID", "ARO", "Model_ID", "Hit_Start", "Hit_End",
        "Percentage Length of Reference Sequence"];
    let columns: Vec<usize> = required.iter().map(|name| {
        header.iter().position(|s| s == name)
            .ok_or_else(|| anyhow!("MakeGeneCoverage: missing column {}", name))
    }).collect::<anyhow::Result<_>>()?;
    let mut groups: BTreeMap<(i64, i64), Vec<(i64, i64)>> = BTreeMap::new();
    let mut lengths = HashMap::new();
    for (index, row) in records.iter().enumerate().skip(1) {
        let get = |col: usize| -> anyhow::Result<&str> {
            row.get(columns[col]).map(String::as_str)
                .ok_or_else(|| anyhow!("MakeGeneCoverage: short TSV row {}", index + 1))
        };
        let number = |col: usize| -> anyhow::Result<i64> {
            get(col)?.parse::<i64>().map_err(|e| anyhow!(
                "MakeGeneCoverage: row {}, invalid {}: {}", index + 1, required[col], e))
        };
        let aro = number(2)?;
        let model = number(3)?;
        let orf = get(0)?;
        let hit = get(1)?;
        let sequence = json.get(orf).and_then(|v| v.get(hit))
            .and_then(|v| v.get("dna_sequence_from_broadstreet"))
            .and_then(serde_json::Value::as_str)
            .ok_or_else(|| anyhow!(
                "MakeGeneCoverage: missing DNA sequence for ORF {:?}, hit {:?}", orf, hit))?;
        // Preserve WDL's last-row-wins length map, keyed only by Model_ID.
        lengths.insert(model, sequence.chars().count());
        let start = number(4)?;
        let end = number(5)?;
        if start == end {
            return Err(anyhow!("MakeGeneCoverage: zero-length hit at row {} (WDL range step undefined)", index + 1));
        }
        let interval = if start < end { (start, end) } else {
            (end.checked_add(1).ok_or_else(|| anyhow!("coordinate overflow"))?,
             start.checked_add(1).ok_or_else(|| anyhow!("coordinate overflow"))?)
        };
        groups.entry((aro, model)).or_default().push(interval);
    }
    let mut output = String::from("ID\tgene_coverage_bps\tdb_seq_length\tgene_coverage_perc\n");
    for ((aro, model), mut intervals) in groups {
        intervals.sort_unstable();
        let mut merged: Option<(i64, i64)> = None;
        let mut covered = 0i128;
        for (start, end) in intervals {
            match merged {
                Some((lo, hi)) if start <= hi => merged = Some((lo, hi.max(end))),
                Some((lo, hi)) => {
                    covered += i128::from(hi) - i128::from(lo);
                    merged = Some((start, end));
                }
                None => merged = Some((start, end)),
            }
        }
        if let Some((lo, hi)) = merged { covered += i128::from(hi) - i128::from(lo); }
        let length = lengths[&model];
        if length == 0 { return Err(anyhow!("MakeGeneCoverage: zero reference length for model {}", model)); }
        let scaled = (covered as f64 / length as f64) * 100.0 * 10000.0;
        // NumPy round uses ties-to-even; do not cap coverage at 100 percent.
        let percentage = scaled.round_ties_even() / 10000.0;
        let mut formatted = format!("{:.4}", percentage);
        while formatted.ends_with('0') { formatted.pop(); }
        if formatted.ends_with('.') { formatted.push('0'); }
        output.push_str(&format!("{}\t{}\t{}\t{}\n", aro, covered, length, formatted));
    }
    Ok(output)
}

async fn make_gene_coverage(
    config: Arc<RunConfig>,
    main: &RgiMainOutputs,
) -> Result<GeneCoverageOutputs, PipelineError> {
    let tsv_path = main.main_amr_results.clone();
    let json_path = main.output_json.clone();
    let input_bytes = tokio::fs::metadata(&json_path).await
        .map_err(|e| PipelineError::IOError(e.to_string()))?.len();
    let work_dir = crate::utils::file::choose_temp_dir(
        input_bytes.max(1), &config.ram_temp_dir, &config.args.nvme_scratch, 4, false,
    ).await?;
    let output_gene_coverage = work_dir.path().join("gene_coverage.tsv");
    let output = tokio::task::spawn_blocking(move || -> anyhow::Result<String> {
        let tsv = std::fs::read_to_string(&tsv_path)
            .map_err(|e| anyhow!("MakeGeneCoverage: {}: {}", tsv_path.display(), e))?;
        let file = std::fs::File::open(&json_path)
            .map_err(|e| anyhow!("MakeGeneCoverage: {}: {}", json_path.display(), e))?;
        let json = serde_json::from_reader(std::io::BufReader::new(file))?;
        gene_coverage_table(&tsv, &json)
    }).await.map_err(|e| PipelineError::Other(anyhow!("MakeGeneCoverage worker: {}", e)))?
        .map_err(PipelineError::Other)?;
    tokio::fs::write(&output_gene_coverage, output).await
        .map_err(|e| PipelineError::IOError(format!("MakeGeneCoverage output: {}", e)))?;
    Ok(GeneCoverageOutputs { _work_dir: work_dir, output_gene_coverage })
}

#[derive(Debug)]
struct RgiKmerBwtOutputs {
    /// Owns returned scratch paths until downstream readers/packaging finish.
    _work_dir: TempDir,
    output_kmer_bwt: Vec<PathBuf>,
    sr_species_allele: PathBuf,
    kma_species_calling: PathBuf,
}

async fn run_rgi_kmer_bwt(
    config: Arc<RunConfig>,
    bwt: &RgiBwtKmaOutputs,
    rgi_db_root: &std::path::Path,
) -> Result<RgiKmerBwtOutputs, PipelineError> {
    let bam = tokio::fs::canonicalize(&bwt.output_sorted_length_100)
        .await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiKmerBwt cannot resolve input {}: {}",
            bwt.output_sorted_length_100.display(), e,
        )))?;
    let metadata = tokio::fs::metadata(&bam).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiKmerBwt cannot stat input {}: {}", bam.display(), e,
        )))?;
    if !metadata.is_file() {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRgiKmerBwt input is not a regular file: {}", bam.display(),
        )));
    }
    // Heuristic, not a quota; choose_temp_dir applies headroom/fallback.
    let work_dir = crate::utils::file::choose_temp_dir(
        metadata.len().saturating_mul(4).max(1),
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    ).await?;
    let work_path = tokio::fs::canonicalize(work_dir.path()).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiKmerBwt cannot resolve scratch directory {}: {}",
            work_dir.path().display(), e,
        )))?;

    tokio::fs::create_dir_all(&config.out_dir)
        .await
        .map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiKmerBwt cannot create log directory {}: {}",
                config.out_dir.display(),
                e,
            ))
        })?;

    let stdout_log_path = config.out_dir.join("amr_rgi_kmer_bwt_stdout.log");
    let stderr_log_path = config.out_dir.join("amr_rgi_kmer_bwt_stderr.log");
    let prefix = work_path.join("sr_species_report");

    let rgi_config = RgiConfig::new(RgiSubcommand::KmerQuery)
        .flag("--bwt")
        .flag("--local")
        .option("-k", "61")
        .option("-i", bam.to_string_lossy().into_owned())
        .option("--output", prefix.to_string_lossy().into_owned());

    let args = generate_cli(RGI_TAG, &config, Some(&rgi_config))?;

    info!(
        "[RunRgiKmerBwt] rgi {:?}; output directory={}; \
         stdout={}; stderr={}",
        args,
        work_path.display(),
        stdout_log_path.display(),
        stderr_log_path.display(),
    );

    let (mut child, stderr_task) = spawn_amr_rgi(
        rgi_db_root,
        args,
        stderr_log_path.clone(),
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!("RunRgiKmerBwt spawn failed: {}", e),
        })?;

    // stderr is redirected to its log. Install a stdout consumer before
    // waiting so neither pipe can backpressure the running process.
    let stdout_rx = match parse_child_output(
        &mut child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        &config,
    )
        .await
    {
        Ok(rx) => rx,
        Err(e) => {
            // Do not leave the process running while work_dir is dropped.
            let kill_result = child.kill().await;
            let wait_result = child.wait().await;
            let stderr_result = stderr_task.await;

            if let Err(kill_error) = kill_result {
                warn!(
                    "[RunRgiKmerBwt] cleanup kill failed: {}",
                    kill_error,
                );
            }
            if let Err(wait_error) = wait_result {
                warn!(
                    "[RunRgiKmerBwt] cleanup wait failed: {}",
                    wait_error,
                );
            }
            if !matches!(&stderr_result, Ok(Ok(()))) {
                warn!(
                    "[RunRgiKmerBwt] stderr capture during cleanup: {:?}",
                    stderr_result,
                );
            }

            return Err(PipelineError::ToolExecution {
                tool: RGI_TAG.to_string(),
                error: format!(
                    "RunRgiKmerBwt stdout capture setup failed: {}; \
                     stderr log: {}",
                    e,
                    stderr_log_path.display(),
                ),
            });
        }
    };

    let stdout_task = tokio::spawn(stream_to_file(
        stdout_rx,
        stdout_log_path.clone(),
    ));

    // Join all three, rather than returning early and detaching a writer.
    let (status_result, stdout_result, stderr_result) =
        tokio::join!(child.wait(), stdout_task, stderr_task);

    // Report secondary capture failures even when the process itself
    // failed. The process error remains the primary failure.
    if !matches!(&stdout_result, Ok(Ok(()))) {
        warn!(
            "[RunRgiKmerBwt] stdout capture failed: {:?}; log={}",
            stdout_result,
            stdout_log_path.display(),
        );
    }

    if !matches!(&stderr_result, Ok(Ok(()))) {
        warn!(
            "[RunRgiKmerBwt] stderr capture failed: {:?}; log={}",
            stderr_result,
            stderr_log_path.display(),
        );
    }

    let status =
        status_result.map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiKmerBwt wait failed: {}; stdout={}; stderr={}",
                e,
                stdout_log_path.display(),
                stderr_log_path.display(),
            ),
        })?;

    if !status.success() {
        return Err(PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiKmerBwt exited with {}; stdout={}; stderr={}",
                status,
                stdout_log_path.display(),
                stderr_log_path.display(),
            ),
        });
    }

    stdout_result
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerBwt stdout writer task failed: {}; log={}",
                e,
                stdout_log_path.display(),
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerBwt stdout write failed: {}; log={}",
                e,
                stdout_log_path.display(),
            ))
        })?;

    stderr_result
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerBwt stderr task failed: {}; log={}",
                e,
                stderr_log_path.display(),
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerBwt stderr capture failed: {}; log={}",
                e,
                stderr_log_path.display(),
            ))
        })?;

    let sr_species_allele =
        work_path.join("sr_species_report_61mer_analysis.allele.txt");
    let kma_species_calling =
        work_path.join("sr_species_report_61mer_analysis.gene.txt");
    for path in [&sr_species_allele, &kma_species_calling] {
        let metadata = tokio::fs::metadata(path).await.map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiKmerBwt missing/unreadable required output {}: {}; stderr={}",
                path.display(), e, stderr_log_path.display(),
            ))
        })?;
        if !metadata.is_file() {
            return Err(PipelineError::IOError(format!(
                "RunRgiKmerBwt output is not a regular file: {}", path.display(),
            )));
        }
        // Do not impose a blanket nonempty-file requirement on reports.
    }

    // WDL glob("sr_species_report*"), scoped to this invocation.
    let mut output_kmer_bwt = Vec::new();
    let mut entries = tokio::fs::read_dir(&work_path).await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    while let Some(entry) = entries.next_entry().await
        .map_err(|e| PipelineError::IOError(e.to_string()))? {
        if entry.file_name().to_string_lossy().starts_with("sr_species_report")
            && entry.file_type().await
            .map_err(|e| PipelineError::IOError(e.to_string()))?.is_file()
        {
            output_kmer_bwt.push(entry.path());
        }
    }
    output_kmer_bwt.sort();
    Ok(RgiKmerBwtOutputs {
        _work_dir: work_dir,
        output_kmer_bwt,
        sr_species_allele,
        kma_species_calling,
    })
}

#[derive(Debug)]
struct RgiKmerMainOutputs {
    /// Owns returned scratch paths until downstream readers/packaging finish.
    _work_dir: TempDir,
    output_kmer_main: Vec<PathBuf>,
    species_calling: PathBuf,
}

async fn run_rgi_kmer_main(
    config: Arc<RunConfig>,
    main: &RgiMainOutputs,
    rgi_db_root: &std::path::Path,
) -> Result<RgiKmerMainOutputs, PipelineError> {
    let input_json = tokio::fs::canonicalize(&main.output_json)
        .await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiKmerMain cannot resolve input {}: {}",
            main.output_json.display(), e,
        )))?;
    let metadata = tokio::fs::metadata(&input_json).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiKmerMain cannot stat input {}: {}", input_json.display(), e,
        )))?;
    if !metadata.is_file() {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRgiKmerMain input is not a regular file: {}", input_json.display(),
        )));
    }
    // Heuristic, not a quota; choose_temp_dir applies headroom/fallback.
    let work_dir = crate::utils::file::choose_temp_dir(
        metadata.len().saturating_mul(4).max(1),
        &config.ram_temp_dir,
        &config.args.nvme_scratch,
        4,
        false,
    ).await?;
    let work_path = tokio::fs::canonicalize(work_dir.path()).await
        .map_err(|e| PipelineError::IOError(format!(
            "RunRgiKmerMain cannot resolve scratch directory {}: {}",
            work_dir.path().display(), e,
        )))?;

    tokio::fs::create_dir_all(&config.out_dir)
        .await
        .map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiKmerMain cannot create log directory {}: {}",
                config.out_dir.display(),
                e,
            ))
        })?;

    let stdout_log_path = config.out_dir.join("amr_rgi_kmer_main_stdout.log");
    let stderr_log_path = config.out_dir.join("amr_rgi_kmer_main_stderr.log");
    let prefix = work_path.join("contig_species_report");

    let rgi_config = RgiConfig::new(RgiSubcommand::KmerQuery)
        .flag("--rgi")
        .flag("--local")
        .option("-k", "61")
        .option("-i", input_json.to_string_lossy().into_owned())
        .option("--output", prefix.to_string_lossy().into_owned());

    let args = generate_cli(RGI_TAG, &config, Some(&rgi_config))?;

    info!(
        "[RunRgiKmerMain] rgi {:?}; output directory={}; \
         stdout={}; stderr={}",
        args,
        work_path.display(),
        stdout_log_path.display(),
        stderr_log_path.display(),
    );

    let (mut child, stderr_task) = spawn_amr_rgi(
        rgi_db_root,
        args,
        stderr_log_path.clone(),
    )
        .await
        .map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!("RunRgiKmerMain spawn failed: {}", e),
        })?;

    // stderr is redirected to its log. Install a stdout consumer before
    // waiting so neither pipe can backpressure the running process.
    let stdout_rx = match parse_child_output(
        &mut child,
        ChildStream::Stdout,
        ParseMode::Bytes,
        &config,
    )
        .await
    {
        Ok(rx) => rx,
        Err(e) => {
            // Do not leave the process running while work_dir is dropped.
            let kill_result = child.kill().await;
            let wait_result = child.wait().await;
            let stderr_result = stderr_task.await;

            if let Err(kill_error) = kill_result {
                warn!(
                    "[RunRgiKmerMain] cleanup kill failed: {}",
                    kill_error,
                );
            }
            if let Err(wait_error) = wait_result {
                warn!(
                    "[RunRgiKmerMain] cleanup wait failed: {}",
                    wait_error,
                );
            }
            if !matches!(&stderr_result, Ok(Ok(()))) {
                warn!(
                    "[RunRgiKmerMain] stderr capture during cleanup: {:?}",
                    stderr_result,
                );
            }

            return Err(PipelineError::ToolExecution {
                tool: RGI_TAG.to_string(),
                error: format!(
                    "RunRgiKmerMain stdout capture setup failed: {}; \
                     stderr log: {}",
                    e,
                    stderr_log_path.display(),
                ),
            });
        }
    };

    let stdout_task = tokio::spawn(stream_to_file(
        stdout_rx,
        stdout_log_path.clone(),
    ));

    // Join all three, rather than returning early and detaching a writer.
    let (status_result, stdout_result, stderr_result) =
        tokio::join!(child.wait(), stdout_task, stderr_task);

    // Report secondary capture failures even when the process itself
    // failed. The process error remains the primary failure.
    if !matches!(&stdout_result, Ok(Ok(()))) {
        warn!(
            "[RunRgiKmerMain] stdout capture failed: {:?}; log={}",
            stdout_result,
            stdout_log_path.display(),
        );
    }

    if !matches!(&stderr_result, Ok(Ok(()))) {
        warn!(
            "[RunRgiKmerMain] stderr capture failed: {:?}; log={}",
            stderr_result,
            stderr_log_path.display(),
        );
    }

    let status =
        status_result.map_err(|e| PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiKmerMain wait failed: {}; stdout={}; stderr={}",
                e,
                stdout_log_path.display(),
                stderr_log_path.display(),
            ),
        })?;

    if !status.success() {
        return Err(PipelineError::ToolExecution {
            tool: RGI_TAG.to_string(),
            error: format!(
                "RunRgiKmerMain exited with {}; stdout={}; stderr={}",
                status,
                stdout_log_path.display(),
                stderr_log_path.display(),
            ),
        });
    }

    stdout_result
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerMain stdout writer task failed: {}; log={}",
                e,
                stdout_log_path.display(),
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerMain stdout write failed: {}; log={}",
                e,
                stdout_log_path.display(),
            ))
        })?;

    stderr_result
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerMain stderr task failed: {}; log={}",
                e,
                stderr_log_path.display(),
            ))
        })?
        .map_err(|e| {
            PipelineError::Other(anyhow!(
                "RunRgiKmerMain stderr capture failed: {}; log={}",
                e,
                stderr_log_path.display(),
            ))
        })?;

    let species_calling =
        work_path.join("contig_species_report_61mer_analysis_rgi_summary.txt");
    for path in [&species_calling] {
        let metadata = tokio::fs::metadata(path).await.map_err(|e| {
            PipelineError::IOError(format!(
                "RunRgiKmerMain missing/unreadable required output {}: {}; stderr={}",
                path.display(), e, stderr_log_path.display(),
            ))
        })?;
        if !metadata.is_file() {
            return Err(PipelineError::IOError(format!(
                "RunRgiKmerMain output is not a regular file: {}", path.display(),
            )));
        }
        // Do not impose a blanket nonempty-file requirement on reports.
    }

    // WDL glob("contig_species_report*"), scoped to this invocation.
    let mut output_kmer_main = Vec::new();
    let mut entries = tokio::fs::read_dir(&work_path).await
        .map_err(|e| PipelineError::IOError(e.to_string()))?;
    while let Some(entry) = entries.next_entry().await
        .map_err(|e| PipelineError::IOError(e.to_string()))? {
        if entry.file_name().to_string_lossy().starts_with("contig_species_report")
            && entry.file_type().await
            .map_err(|e| PipelineError::IOError(e.to_string()))?.is_file()
        {
            output_kmer_main.push(entry.path());
        }
    }
    output_kmer_main.sort();
    Ok(RgiKmerMainOutputs {
        _work_dir: work_dir,
        output_kmer_main,
        species_calling,
    })
}

// Native RunResultsPerSample. No dataframe runtime: indexed outer joins retain
// duplicate-key Cartesian products; sets are sorted for reproducible output.
type AmrResultRow = HashMap<String, Option<String>>;

#[derive(Clone, Debug)]
struct AmrResultTable {
    columns: Vec<String>,
    rows: Vec<AmrResultRow>,
}

fn amr_result_na(s: &str) -> bool {
    matches!(s, "" | "#N/A" | "#N/A N/A" | "#NA" | "-1.#IND" | "-1.#QNAN" |
        "-NaN" | "-nan" | "1.#IND" | "1.#QNAN" | "<NA>" | "N/A" | "NA" |
        "NULL" | "NaN" | "None" | "n/a" | "nan" | "null")
}

impl AmrResultTable {
    fn read(path: &std::path::Path) -> anyhow::Result<Self> {
        use anyhow::Context;
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)
            .with_context(|| format!("RunResultsPerSample: reading {}", path.display()))?;
        let columns: Vec<String> = reader.headers()?.iter().map(str::to_owned).collect();
        if columns.is_empty() { return Err(anyhow!("Missing TSV header: {}", path.display())); }
        let mut seen = HashSet::new();
        for c in &columns {
            if !seen.insert(c) { return Err(anyhow!("Duplicate column {:?} in {}", c, path.display())); }
        }
        let mut rows = Vec::new();
        for record in reader.records() {
            let record = record.with_context(|| format!("Invalid TSV in {}", path.display()))?;
            rows.push(columns.iter().cloned().zip(record.iter().map(|v| {
                if amr_result_na(v) { None } else { Some(v.to_owned()) }
            })).collect());
        }
        Ok(Self { columns, rows })
    }

    fn require(&self, names: &[&str]) -> anyhow::Result<()> {
        for name in names {
            if !self.columns.iter().any(|c| c == name) {
                return Err(anyhow!("RunResultsPerSample: missing column {:?}", name));
            }
        }
        Ok(())
    }

    fn clean_sort(&mut self, name: &str) -> anyhow::Result<()> {
        self.require(&[name])?;
        for row in &mut self.rows {
            if let Some(Some(v)) = row.get_mut(name) { *v = v.trim().to_owned(); }
        }
        self.sort(name);
        Ok(())
    }

    fn sort(&mut self, name: &str) {
        self.rows.sort_by(|a, b| match (amr_cell(a, name), amr_cell(b, name)) {
            (Some(a), Some(b)) => a.cmp(b),
            (None, None) => std::cmp::Ordering::Equal,
            (None, _) => std::cmp::Ordering::Greater,
            (_, None) => std::cmp::Ordering::Less,
        });
    }

    fn suffix(&mut self, suffix: &str) {
        self.columns = self.columns.iter().map(|c| format!("{}{}", c, suffix)).collect();
        for row in &mut self.rows {
            *row = std::mem::take(row).into_iter().map(|(k, v)| (format!("{}{}", k, suffix), v)).collect();
        }
    }

    fn coalesce(&mut self, dest: &str, left: &str, right: &str) -> anyhow::Result<()> {
        self.require(&[left, right])?;
        self.columns.push(dest.to_owned());
        for row in &mut self.rows {
            let value = row.get(left).cloned().flatten().or_else(|| row.get(right).cloned().flatten());
            row.insert(dest.to_owned(), value);
        }
        Ok(())
    }

    // Reproduce pandas astype(str) after an outer join: integer columns with
    // missing values are promoted to float, and missing cells become "nan".
    fn accession_strings(&mut self, column: &str) -> anyhow::Result<()> {
        self.require(&[column])?;
        let missing = self.rows.iter().any(|r| amr_cell(r, column).is_none());
        let numeric = self.rows.iter().filter_map(|r| amr_cell(r, column))
            .all(|s| s.parse::<f64>().is_ok());
        let floating = missing || self.rows.iter().filter_map(|r| amr_cell(r, column))
            .any(|s| s.contains('.') || s.contains('e') || s.contains('E'));
        for row in &mut self.rows {
            let value = match amr_cell(row, column) {
                None => "nan".to_owned(),
                Some(v) if numeric && floating => amr_float(v.parse::<f64>()?),
                Some(v) => v.to_owned(),
            };
            row.insert(column.to_owned(), Some(value));
        }
        Ok(())
    }

    fn outer(&self, other: &Self, left: &[&str], right: &[&str]) -> anyhow::Result<Self> {
        self.require(left)?;
        other.require(right)?;
        if left.len() != right.len() { return Err(anyhow!("Join key length mismatch")); }
        if other.columns.iter().any(|c| self.columns.contains(c)) {
            return Err(anyhow!("RunResultsPerSample: overlapping join columns"));
        }
        let key = |row: &AmrResultRow, cols: &[&str]| -> Vec<Option<String>> {
            cols.iter().map(|c| row.get(*c).cloned().flatten()).collect()
        };
        let mut index: HashMap<Vec<Option<String>>, Vec<usize>> = HashMap::new();
        for (i, row) in other.rows.iter().enumerate() {
            index.entry(key(row, right)).or_default().push(i);
        }
        let mut used = vec![false; other.rows.len()];
        let mut rows = Vec::new();
        for row in &self.rows {
            if let Some(matches) = index.get(&key(row, left)) {
                for &i in matches {
                    let mut merged = row.clone();
                    merged.extend(other.rows[i].clone());
                    rows.push(merged);
                    used[i] = true;
                }
            } else { rows.push(row.clone()); }
        }
        for (i, row) in other.rows.iter().enumerate() {
            if !used[i] { rows.push(row.clone()); }
        }
        let columns = self.columns.iter().chain(&other.columns).cloned().collect();
        Ok(Self { columns, rows })
    }

    fn write(&self, path: &std::path::Path, selected: Option<&[&str]>) -> anyhow::Result<()> {
        let columns: Vec<&str> = selected.map(|v| v.to_vec())
            .unwrap_or_else(|| self.columns.iter().map(String::as_str).collect());
        self.require(&columns)?;
        let mut writer = csv::WriterBuilder::new().delimiter(b'\t')
            .terminator(csv::Terminator::Any(b'\n')).from_path(path)?;
        writer.write_record(&columns)?;
        for row in &self.rows {
            writer.write_record(columns.iter().map(|c| amr_cell(row, c).unwrap_or("")))?;
        }
        writer.flush()?;
        Ok(())
    }
}

fn amr_cell<'a>(row: &'a AmrResultRow, column: &str) -> Option<&'a str> {
    row.get(column).and_then(|v| v.as_deref())
}

fn amr_float(value: f64) -> String {
    let mut s = value.to_string();
    if !s.contains('.') && !s.contains('e') && value.is_finite() { s.push_str(".0"); }
    s
}

fn amr_unique(rows: &[&AmrResultRow], columns: &[&str]) -> std::collections::BTreeSet<String> {
    rows.iter().flat_map(|row| columns.iter().filter_map(move |c| {
        amr_cell(row, c).filter(|s| !amr_result_na(s)).map(str::to_owned)
    })).collect()
}

fn amr_join(values: &std::collections::BTreeSet<String>, separator: &str) -> Option<String> {
    if values.is_empty() { None } else { Some(values.iter().cloned().collect::<Vec<_>>().join(separator)) }
}

fn amr_numbers(rows: &[&AmrResultRow], column: &str) -> anyhow::Result<Vec<f64>> {
    let mut numbers = Vec::new();
    for value in amr_unique(rows, &[column]) {
        let n = value.parse::<f64>().map_err(|e| anyhow!("RunResultsPerSample: invalid {} {:?}: {}", column, value, e))?;
        if n.is_finite() && !numbers.contains(&n) { numbers.push(n); }
    }
    Ok(numbers)
}

fn amr_max(rows: &[&AmrResultRow], column: &str) -> anyhow::Result<Option<String>> {
    Ok(amr_numbers(rows, column)?.into_iter().reduce(f64::max).map(amr_float))
}

const AMR_BIG_COLUMNS: &[&str] = &[
    "ARO_overall", "Gene_Family_overall", "Drug_Class_overall", "Resistance_Mechanism_overall",
    "model_overall", "All Mapped Reads_kma_amr", "Percent Coverage_kma_amr", "Depth_kma_amr",
    "CARD*kmer Prediction_kma_sp", "Cut_Off_contig_amr",
    "Percentage Length of Reference Sequence_contig_amr", "Best_Identities_contig_amr",
    "CARD*kmer Prediction_contig_sp",
];
const AMR_PRIMARY_COLUMNS: &[&str] = &[
    "gene_name", "sample_name", "gene_family", "aro_accession", "drug_class", "high_level_drug_class",
    "resistance_mechanism", "model_type", "num_contigs", "cutoff", "contig_coverage_breadth",
    "contig_percent_id", "contig_species", "num_reads", "read_gene_id", "read_coverage_breadth",
    "read_coverage_depth", "read_species",
];

/// Pure file-backed transformation, also usable for fixture/differential tests.
fn amr_results_per_sample_files(
    inputs: &[PathBuf; 5], ontology_path: &std::path::Path, sample: &str,
    directory: &std::path::Path,
) -> anyhow::Result<()> {
    let mut main = AmrResultTable::read(&inputs[0])?;
    let mut main_sp = AmrResultTable::read(&inputs[1])?;
    let mut kma = AmrResultTable::read(&inputs[2])?;
    let mut kma_sp = AmrResultTable::read(&inputs[3])?;
    main.clean_sort("Best_Hit_ARO")?;
    main_sp.clean_sort("Best_Hit_ARO")?;
    kma.clean_sort("ARO Term")?;
    kma_sp.clean_sort("ARO term")?;
    main.suffix("_contig_amr");
    main_sp.suffix("_contig_sp");
    main.require(&["Contig_contig_amr"])?;
    for row in &mut main.rows {
        if let Some(Some(v)) = row.get_mut("Contig_contig_amr") { *v = v.trim().to_owned(); }
    }
    let mut a = main.outer(&main_sp,
                           &["Best_Hit_ARO_contig_amr", "Contig_contig_amr"],
                           &["Best_Hit_ARO_contig_sp", "Contig_contig_sp"])?;
    a.coalesce("ARO_contig", "Best_Hit_ARO_contig_amr", "Best_Hit_ARO_contig_sp")?;
    a.write(&directory.join("merge_a.tsv"), None)?;
    a.accession_strings("ARO_contig_amr")?;
    kma.suffix("_kma_amr");
    kma_sp.suffix("_kma_sp");
    let mut b = kma.outer(&kma_sp, &["ARO Term_kma_amr"], &["ARO term_kma_sp"])?;
    b.accession_strings("ARO Accession_kma_amr")?;
    b.require(&["Reference Model Type_kma_amr"])?;
    b.rows.retain(|r| !matches!(amr_cell(r, "Reference Model Type_kma_amr"),
        Some("protein variant model" | "rRNA gene variant model")));
    b.coalesce("ARO_kma", "ARO Term_kma_amr", "ARO term_kma_sp")?;
    b.write(&directory.join("merge_b.tsv"), None)?;
    let mut x = a.outer(&b, &["ARO_contig_amr"], &["ARO Accession_kma_amr"])?;
    for (dest, left, right) in [
        ("ARO_accession", "ARO_contig_amr", "ARO Accession_kma_amr"),
        ("ARO_overall", "ARO_contig", "ARO_kma"),
        ("Gene_Family_overall", "AMR Gene Family_contig_amr", "AMR Gene Family_kma_amr"),
        ("Drug_Class_overall", "Drug Class_contig_amr", "Drug Class_kma_amr"),
        ("model_overall", "Model_type_contig_amr", "Reference Model Type_kma_amr"),
        ("Resistance_Mechanism_overall", "Resistance Mechanism_contig_amr", "Resistance Mechanism_kma_amr"),
    ] { x.coalesce(dest, left, right)?; }
    x.write(&directory.join("comprehensive_AMR_metrics.tsv"), None)?;
    // WDL rereads this file with pandas, turning literal "nan" back into NA.
    for row in &mut x.rows {
        for value in row.values_mut() {
            if value.as_deref().map(amr_result_na).unwrap_or(false) { *value = None; }
        }
        if let Some(value) = amr_cell(row, "ARO_contig_amr") {
            let n: f64 = value.parse()?;
            row.insert("ARO_contig_amr".into(), Some(format!("{:.0}", n)));
        }
    }
    x.require(AMR_BIG_COLUMNS)?;
    let mut big = x.clone();
    big.sort("Gene_Family_overall");
    big.write(&directory.join("bigtable_report.tsv"), Some(AMR_BIG_COLUMNS))?;
    if x.rows.is_empty() {
        std::fs::write(directory.join("primary_AMR_report.tsv"), b"")?;
        return Ok(());
    }
    let coverage = AmrResultTable::read(&inputs[4])?;
    coverage.require(&["ID", "gene_coverage_perc"])?;
    let ontology: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(ontology_path)?))?;
    if !ontology.is_object() { return Err(anyhow!("CARD ontology must be a JSON object")); }
    x.require(&["Nudged_contig_amr", "Reference Sequence_kma_amr"])?;
    let mut groups: std::collections::BTreeMap<String, Vec<&AmrResultRow>> = std::collections::BTreeMap::new();
    for row in &x.rows {
        if let Some(name) = amr_cell(row, "ARO_overall") { groups.entry(name.to_owned()).or_default().push(row); }
    }
    let mut primary = AmrResultTable { columns: AMR_PRIMARY_COLUMNS.iter().map(|s| (*s).to_owned()).collect(), rows: Vec::new() };
    for (name, rows) in groups {
        let mut out = AmrResultRow::new();
        out.insert("gene_name".into(), Some(name.clone()));
        out.insert("sample_name".into(), Some(sample.to_owned()));
        for (dest, cols) in [
            ("gene_family", ["AMR Gene Family_contig_amr", "AMR Gene Family_kma_amr"]),
            ("drug_class", ["Drug Class_contig_amr", "Drug Class_kma_amr"]),
            ("resistance_mechanism", ["Resistance Mechanism_contig_amr", "Resistance Mechanism_kma_amr"]),
        ] { out.insert(dest.into(), amr_join(&amr_unique(&rows, &cols), ";")); }
        let accessions = amr_unique(&rows, &["ARO_accession"]).into_iter().map(|s| {
            let n: f64 = s.parse()?;
            if !n.is_finite() { return Err(anyhow!("Invalid ARO accession {:?}", s)); }
            Ok(format!("ARO:{:.0}", n.trunc()))
        }).collect::<anyhow::Result<std::collections::BTreeSet<_>>>()?;
        out.insert("aro_accession".into(), amr_join(&accessions, ";"));
        let mut high = std::collections::BTreeSet::new();
        for dc in amr_unique(&rows, &["Drug Class_contig_amr", "Drug Class_kma_amr"]) {
            for class in dc.split("; ") {
                if let Some(entry) = ontology.get(class) {
                    if let Some(classes) = entry.get("highLevelDrugClasses") {
                        let classes = classes.as_array().ok_or_else(|| anyhow!("Invalid highLevelDrugClasses for {}", class))?;
                        for c in classes {
                            high.insert(c.as_str().ok_or_else(|| anyhow!("Non-string drug class"))?.to_owned());
                        }
                    }
                }
            }
        }
        out.insert("high_level_drug_class".into(), amr_join(&high, "; "));
        let mut cutoffs = std::collections::BTreeSet::new();
        for row in &rows {
            if let Some(c) = amr_cell(row, "Cut_Off_contig_amr") {
                let nudged = amr_cell(row, "Nudged_contig_amr").map(|s| s.eq_ignore_ascii_case("true")).unwrap_or(false);
                cutoffs.insert(if c == "Strict" && nudged { "Nudged".to_owned() } else { c.to_owned() });
            }
        }
        out.insert("cutoff".into(), amr_join(&cutoffs, ";"));
        let models: Vec<String> = amr_unique(&rows, &["Model_type_contig_amr", "Reference Model Type_kma_amr"])
            .into_iter().map(|s| { let mut words: Vec<&str> = s.split(' ').collect(); words.pop(); words.join(" ") }).collect();
        out.insert("model_type".into(), if models.is_empty() { None } else { Some(models.join(";")) });
        for (dest, col) in [("read_coverage_breadth", "Percent Coverage_kma_amr"),
            ("read_coverage_depth", "Depth_kma_amr"), ("num_reads", "All Mapped Reads_kma_amr")] {
            out.insert(dest.into(), amr_max(&rows, col)?);
        }
        let gene_id = amr_join(&amr_unique(&rows, &["ARO_contig_amr"]), ";");
        let breadth = match gene_id {
            None => None,
            Some(id) => {
                let row = coverage.rows.iter().find(|r| amr_cell(r, "ID") == Some(id.as_str()))
                    .ok_or_else(|| anyhow!("RunResultsPerSample: no gene coverage for {:?} (gene {:?})", id, name))?;
                amr_cell(row, "gene_coverage_perc").map(str::to_owned)
            }
        };
        out.insert("contig_coverage_breadth".into(), breadth);
        out.insert("num_contigs".into(), Some(amr_unique(&rows, &["Contig_contig_amr"]).len().to_string()));
        out.insert("read_gene_id".into(), amr_join(&amr_unique(&rows, &["Reference Sequence_kma_amr"]), ";"));
        let ids = amr_numbers(&rows, "Best_Identities_contig_amr")?;
        out.insert("contig_percent_id".into(), if ids.is_empty() { None } else { Some(amr_float(ids.iter().sum::<f64>() / ids.len() as f64)) });
        for (dest, col) in [("contig_species", "CARD*kmer Prediction_contig_sp"), ("read_species", "CARD*kmer Prediction_kma_sp")] {
            out.insert(dest.into(), Some(amr_join(&amr_unique(&rows, &[col]), " ").unwrap_or_default()));
        }
        // Preserve validation performed by the WDL's otherwise-unused species
        // dictionary; its contents do not contribute to the output columns.
        for prediction in amr_unique(&rows, &["CARD*kmer Prediction_contig_sp", "CARD*kmer Prediction_kma_sp"]) {
            if prediction.contains("Unknown") { continue; }
            for item in prediction.split(';') {
                if item.trim().contains(':') {
                    let count = item.split(':').nth(1).unwrap_or("").trim();
                    count.parse::<i64>().map_err(|e| anyhow!("Invalid species count {:?}: {}", count, e))?;
                }
            }
        }
        if out.get("drug_class").and_then(|v| v.as_ref()).is_some() { primary.rows.push(out); }
    }
    primary.write(&directory.join("primary_AMR_report.tsv"), None)?;
    Ok(())
}

#[derive(Debug)]
struct ResultsPerSampleOutputs {
    _work_dir: TempDir,
    artifacts: Vec<PathBuf>,
    final_summary: PathBuf,
    synthesized_report: PathBuf,
}

async fn run_results_per_sample(
    config: Arc<RunConfig>, main: &RgiMainOutputs, main_species: &RgiKmerMainOutputs,
    kma: &RgiBwtKmaOutputs, kma_species: &RgiKmerBwtOutputs, coverage: &GeneCoverageOutputs,
    ontology: &std::path::Path, sample: &str,
) -> Result<ResultsPerSampleOutputs, PipelineError> {
    let inputs = [main.main_amr_results.clone(), main_species.species_calling.clone(),
        kma.kma_amr_results.clone(), kma_species.kma_species_calling.clone(), coverage.output_gene_coverage.clone()];
    let mut bytes = 0u64;
    for path in &inputs {
        bytes = bytes.saturating_add(tokio::fs::metadata(path).await
            .map_err(|e| PipelineError::IOError(format!("RunResultsPerSample {}: {}", path.display(), e)))?.len());
    }
    let work_dir = crate::utils::file::choose_temp_dir(bytes.saturating_mul(4).max(1),
                                                       &config.ram_temp_dir, &config.args.nvme_scratch, 4, false).await?;
    let directory = work_dir.path().to_path_buf();
    let worker_dir = directory.clone();
    let ontology = ontology.to_path_buf();
    let sample = sample.to_owned();
    tokio::task::spawn_blocking(move || amr_results_per_sample_files(&inputs, &ontology, &sample, &worker_dir))
        .await.map_err(|e| PipelineError::Other(anyhow!("RunResultsPerSample worker: {}", e)))?
        .map_err(PipelineError::Other)?;
    let artifacts = ["merge_a.tsv", "merge_b.tsv", "comprehensive_AMR_metrics.tsv", "bigtable_report.tsv", "primary_AMR_report.tsv"]
        .iter().map(|name| directory.join(name)).collect();
    Ok(ResultsPerSampleOutputs { _work_dir: work_dir, artifacts,
        final_summary: directory.join("comprehensive_AMR_metrics.tsv"),
        synthesized_report: directory.join("primary_AMR_report.tsv") })
}

#[derive(Debug)]
struct TsvToSamOutputs {
    _work_dir: TempDir,
    output_sorted: PathBuf,
    output_sorted_bai: PathBuf,
}

/// WDL re.sub(r'_\d*$', '', name), including removal of a trailing underscore.
fn amr_bam_contig_name(name: &str) -> &str {
    match name.rsplit_once('_') {
        Some((prefix, suffix)) if suffix.bytes().all(|b| b.is_ascii_digit()) => prefix,
        _ => name,
    }
}

/// Write the WDL's synthetic alignments as SAM. samtools performs only encoding
/// and indexing, not alignment or sorting. Duplicate summary rows are retained.
fn amr_write_contig_sam(
    contigs: &std::path::Path, summary: &std::path::Path, sam: &std::path::Path,
) -> anyhow::Result<()> {
    use std::io::{BufRead, Write};
    let mut fasta = std::io::BufReader::new(std::fs::File::open(contigs)?);
    let mut first = String::new();
    fasta.read_line(&mut first)?;
    let failed = first == ";ASSEMBLY FAILED\n";
    let mut rows = Vec::<(String, String)>::new();
    let mut references = Vec::<String>::new();
    let mut seen = HashSet::new();
    if !failed {
        let table = AmrResultTable::read(summary)?;
        table.require(&["ARO_contig_amr", "Contig_contig_amr"])?;
        for row in &table.rows {
            let Some(name) = amr_cell(row, "Contig_contig_amr") else { continue; };
            let accession = amr_cell(row, "ARO_contig_amr")
                .ok_or_else(|| anyhow!("tsvToSam: contig {:?} has no ARO", name))?;
            let n: f64 = accession.parse()?;
            if !n.is_finite() { return Err(anyhow!("tsvToSam: invalid ARO {:?}", accession)); }
            let aro = format!("ARO:{:.0}", n.trunc());
            if seen.insert(aro.clone()) { references.push(aro.clone()); }
            rows.push((aro, amr_bam_contig_name(name).to_owned()));
        }
    }
    let wanted: HashSet<&str> = rows.iter().map(|(_, n)| n.as_str()).collect();
    let mut sequences = HashMap::<String, String>::new();
    if !failed && !wanted.is_empty() {
        // Fetch equivalent: FASTA IDs end at the first whitespace. Keep only
        // requested contigs, and retain the first occurrence of duplicate IDs.
        let mut name = String::new();
        let mut sequence = String::new();
        let mut keep = false;
        for line in std::iter::once(Ok(first)).chain(fasta.lines()) {
            let line = line?;
            if let Some(header) = line.strip_prefix('>') {
                if keep { sequences.entry(name.clone()).or_insert(std::mem::take(&mut sequence)); }
                sequence.clear();
                name = header.split_whitespace().next().unwrap_or("").to_owned();
                keep = wanted.contains(name.as_str());
            } else if keep {
                sequence.extend(line.chars().filter(|c| !c.is_whitespace()));
            }
        }
        if keep { sequences.entry(name).or_insert(sequence); }
    }
    let mut writer = std::io::BufWriter::new(std::fs::File::create(sam)?);
    if references.is_empty() { references.push("NoGenes".to_owned()); }
    for reference in &references { writeln!(writer, "@SQ\tSN:{}\tLN:100", reference)?; }
    for (aro, name) in rows {
        let sequence = sequences.get(&name)
            .ok_or_else(|| anyhow!("tsvToSam: contig {:?} not found in {}", name, contigs.display()))?;
        // pysam defaults: flag=0, MAPQ=0, no CIGAR/mate/qualities. Its
        // reference_start=1 is zero-based, hence SAM POS=2, not POS=1.
        writeln!(writer, "{}\t0\t{}\t2\t0\t*\t*\t0\t0\t{}\t*",
                 name, aro, if sequence.is_empty() { "*" } else { sequence })?;
    }
    writer.flush()?;
    Ok(())
}

async fn tsv_to_sam(
    config: Arc<RunConfig>, contigs: &std::path::Path, results: &ResultsPerSampleOutputs,
) -> Result<TsvToSamOutputs, PipelineError> {
    let bytes = tokio::fs::metadata(contigs).await
        .map_err(|e| PipelineError::IOError(format!("tsvToSam input: {}", e)))?.len();
    let work_dir = crate::utils::file::choose_temp_dir(bytes.saturating_mul(4).max(1),
                                                       &config.ram_temp_dir, &config.args.nvme_scratch, 4, false).await?;
    let sam = work_dir.path().join("contig_amr_report.sam");
    let output_sorted = work_dir.path().join("contig_amr_report.sorted.bam");
    let output_sorted_bai = work_dir.path().join("contig_amr_report.sorted.bam.bai");
    let input = contigs.to_path_buf();
    let summary = results.final_summary.clone();
    let worker_sam = sam.clone();
    tokio::task::spawn_blocking(move || amr_write_contig_sam(&input, &summary, &worker_sam))
        .await.map_err(|e| PipelineError::Other(anyhow!("tsvToSam worker: {}", e)))?
        .map_err(PipelineError::Other)?;
    // Match WDL: despite the filename, it does not sort records. Do not silently
    // reorder them here; index failure is an error, just as with pysam.index.
    for (stage, args) in [
        ("encode", vec![std::ffi::OsString::from("view"), "--no-PG".into(), "-b".into(), "-o".into(),
                        output_sorted.as_os_str().to_owned(), sam.as_os_str().to_owned()]),
        ("index", vec![std::ffi::OsString::from("index"), output_sorted.as_os_str().to_owned(),
                       output_sorted_bai.as_os_str().to_owned()]),
    ] {
        let stdout = config.out_dir.join(format!("amr_tsv_to_sam_{}_stdout.log", stage));
        let stderr = config.out_dir.join(format!("amr_tsv_to_sam_{}_stderr.log", stage));
        let out_file = tokio::fs::File::create(&stdout).await.map_err(|e| PipelineError::IOError(e.to_string()))?.into_std().await;
        let err_file = tokio::fs::File::create(&stderr).await.map_err(|e| PipelineError::IOError(e.to_string()))?.into_std().await;
        let status = tokio::process::Command::new("samtools").args(args)
            .stdin(std::process::Stdio::null()).stdout(std::process::Stdio::from(out_file))
            .stderr(std::process::Stdio::from(err_file)).kill_on_drop(true).status().await
            .map_err(|e| PipelineError::Other(anyhow!("tsvToSam samtools {}: {}", stage, e)))?;
        if !status.success() { return Err(PipelineError::Other(anyhow!(
            "tsvToSam samtools {} exited {}; stderr={}", stage, status, stderr.display()))); }
    }
    for path in [&output_sorted, &output_sorted_bai] {
        if !tokio::fs::metadata(path).await.map_err(|e| PipelineError::IOError(e.to_string()))?.is_file() {
            return Err(PipelineError::Other(anyhow!("tsvToSam missing {}", path.display())));
        }
    }
    Ok(TsvToSamOutputs { _work_dir: work_dir, output_sorted, output_sorted_bai })
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
/// The RunRedup streams have already been materialized at the AMR external-tool
/// boundary. SPAdes consumes those durable files, which are also retained as the
/// downstream `non_host_reads_fa` inputs.
async fn run_spades_amr(
    config: Arc<RunConfig>,
    input_paths: Vec<PathBuf>,
) -> Result<RunSpadesHandle, PipelineError> {
    if input_paths.is_empty() || input_paths.len() > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunSpades expected one or two reduplicated FASTA files, got {}",
            input_paths.len()
        )));
    }

    for path in &input_paths {
        if !path.is_file() {
            return Err(PipelineError::FileNotFound(path.clone()));
        }
    }

    let spades_out_dir = config.out_dir.join("spades");
    let contigs_path = spades_out_dir.join("contigs.fasta");
    let scaffolds_path = spades_out_dir.join("scaffolds.fasta");

    let completion_config = config.clone();
    let completion_contigs = contigs_path.clone();
    let completion_scaffolds = scaffolds_path.clone();

    let completion_task = tokio::spawn(async move {
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
            tempdir_path: None,
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
            ParseMode::Bytes,
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




async fn materialize_redups(
    config: Arc<RunConfig>,
    reduplicated_streams: Vec<ReceiverStream<ParseOutput>>,
    out_dir: &PathBuf,
    tag: &str,
) -> Result<(Vec<PathBuf>, Vec<JoinHandle<Result<(), anyhow::Error>>>), PipelineError> {
    if reduplicated_streams.is_empty() || reduplicated_streams.len() > 2 {
        return Err(PipelineError::InvalidConfig(format!(
            "RunRedup expected one or two output streams, got {}",
            reduplicated_streams.len()
        )));
    }

    let mut paths = Vec::with_capacity(reduplicated_streams.len());
    let mut tasks = Vec::with_capacity(reduplicated_streams.len());

    for (idx, stream) in reduplicated_streams.into_iter().enumerate() {
        let path = out_dir.join(format!("redups_{}.fa", idx + 1));
        tasks.push(write_fasta_stream_to_file(
            stream,
            path.clone(),
            config.clone(),
            StreamDataType::JustBytes,
            &format!("{} part {}", tag, idx + 1),
        ));
        paths.push(path);
    }

    Ok((paths, tasks))
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

/// Schedule durable artifact writers without blocking downstream analysis.
/// The caller retains the source TempDir owners until cleanup_tasks are joined.
fn schedule_rgi_artifact_copies(
    config: Arc<RunConfig>,
    artifacts: &[PathBuf],
    cleanup_tasks: &mut Vec<JoinHandle<anyhow::Result<()>>>,
) {
    for source in artifacts {
        let source = source.clone();
        let config = config.clone();
        cleanup_tasks.push(tokio::spawn(async move {
            let name = source.file_name().ok_or_else(|| anyhow!(
                "RGI artifact has no filename: {}", source.display(),
            ))?;
            let destination = config.out_dir.join(name);
            let writer = crate::utils::file::copy_file_streaming(
                &source, &destination, Some(config), None,
            ).await?;
            writer.await.map_err(|e| anyhow!(
                "RGI artifact writer panicked for {}: {}", source.display(), e,
            ))??;
            info!("[AMR] preserved RGI artifact: {}", destination.display());
            Ok(())
        }));
    }
}

pub async fn run(
    config: Arc<RunConfig>,
) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir()
        .map_err(|e| PipelineError::Other(e.into()))?;
    let mut cleanup_tasks: Vec<JoinHandle<anyhow::Result<()>>> = Vec::new();
    let out_dir = config.out_dir.clone();

    info!("Starting AMR pipeline.");
    // Fail fast on report-only configuration before expensive analysis starts.
    let report_ontology = config.args.amr_card_ontology.as_deref().ok_or_else(||
        PipelineError::InvalidConfig("--amr-card-ontology is required for AMR reporting (ontology.json, not card.json)".into()))?;
    let report_ontology = resolve_existing_input_path(report_ontology, &cwd, true)?;
    let rgi_db_root = resolve_amr_rgi_root(&config).await?;
    info!("[AMR] prepared RGI database: {}", rgi_db_root.display());

    let prepared = match AmrInputMode::from_config(&config)? {
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
                sample_base,
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

            let (redup_paths, redup_write_tasks) = materialize_redups(
                config.clone(),
                redup.redups_streams,
                &out_dir,
                "AMR raw RunRedup",
            )
                .await?;

            let redup_completion_task = redup.completion_task;

            // ------------------------------------------------------------
            // Drain validation / QC results.
            // ------------------------------------------------------------





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
                    .map_err(|e| PipelineError::Other(anyhow!("Raw RunRedup writer task panicked: {}", e)))?
                    .map_err(PipelineError::Other)?;
            }

            let spades = run_spades_amr(config.clone(), redup_paths.clone()).await?;
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

            AmrPreparedInputs {
                sample_base,
                non_host_reads_fa: redup_paths,
                contigs_fa: spades.contigs_path,
            }
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


            let (file1_path, file2_path, _, sample_base) =
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

            let filtered_contigs = config
                .args
                .amr_filtered_contigs
                .as_deref()
                .ok_or_else(|| PipelineError::InvalidConfig(
                    "--amr-filtered-contigs is required for --amr-input-mode filtered".to_string(),
                ))?;
            let filtered_contigs = resolve_existing_input_path(filtered_contigs, &cwd, true)?;

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



            let (redup_paths, redup_write_tasks) = materialize_redups(
                config.clone(),
                redup.redups_streams,
                &out_dir,
                "AMR filtered RunRedup",
            )
                .await?;

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

            AmrPreparedInputs {
                sample_base,
                non_host_reads_fa: redup_paths,
                contigs_fa: filtered_contigs,
            }
        }
    };
    info!(
    "[AMR] prepared common inputs: reads={:?}, contigs={}",
    prepared.non_host_reads_fa,
    prepared.contigs_fa.display(),
);

    // Keep this output owner in scope through all downstream RGI consumers.
    let rgi_bwt_kma = run_rgi_bwt_kma(
        config.clone(),
        &prepared.non_host_reads_fa,
        &rgi_db_root,
    ).await?;
    info!(
        "[AMR] RunRgiBwtKma complete: {} artifacts, allele report={}",
        rgi_bwt_kma.output_kma.len(),
        rgi_bwt_kma.kma_amr_results.display(),
    );
    schedule_rgi_artifact_copies(
        config.clone(), &rgi_bwt_kma.output_kma, &mut cleanup_tasks,
    );

    // Do not return early on query failure: first finish the BWT writers
    // while their source directory owner is still alive.
    let kmer_result = run_rgi_kmer_bwt(
        config.clone(), &rgi_bwt_kma, &rgi_db_root,
    ).await;
    if let Ok(rgi_kmer_bwt) = &kmer_result {
        info!(
            "[AMR] RunRgiKmerBwt complete: {} artifacts, allele={}, gene={}",
            rgi_kmer_bwt.output_kmer_bwt.len(),
            rgi_kmer_bwt.sr_species_allele.display(),
            rgi_kmer_bwt.kma_species_calling.display(),
        );
        schedule_rgi_artifact_copies(
            config.clone(), &rgi_kmer_bwt.output_kmer_bwt, &mut cleanup_tasks,
        );
    }

    // RunRgiMain is independent of the BWT k-mer result. Preserve that work
    // even if the k-mer query failed, and defer error propagation until after
    // all independently scheduled artifact writers have been drained.
    let main_result = run_rgi_main(
        config.clone(), &prepared.contigs_fa, &rgi_db_root,
    ).await;
    if let Ok(rgi_main) = &main_result {
        info!(
            "[AMR] RunRgiMain complete: {} artifacts, main={}, json={}",
            rgi_main.output_main.len(),
            rgi_main.main_amr_results.display(),
            rgi_main.output_json.display(),
        );
        schedule_rgi_artifact_copies(
            config.clone(), &rgi_main.output_main, &mut cleanup_tasks,
        );
    }

    // Consume scratch results directly; publication copies are independent.
    // Defer errors until all already-scheduled writers have drained.
    let coverage_result = match &main_result {
        Ok(main) => Some(make_gene_coverage(config.clone(), main).await),
        Err(_) => None,
    };
    if let Some(Ok(coverage)) = &coverage_result {
        info!("[AMR] MakeGeneCoverage complete: {}", coverage.output_gene_coverage.display());
        schedule_rgi_artifact_copies(
            config.clone(), std::slice::from_ref(&coverage.output_gene_coverage),
            &mut cleanup_tasks,
        );
    }

    // Run against the main JSON directly, reusing the prepared local database.
    // Keep the owner alive until publication and downstream consumers finish.
    let main_kmer_result = match &main_result {
        Ok(main) => Some(run_rgi_kmer_main(config.clone(), main, &rgi_db_root).await),
        Err(_) => None,
    };
    if let Some(Ok(main_kmer)) = &main_kmer_result {
        info!(
            "[AMR] RunRgiKmerMain complete: {} artifacts, species={}",
            main_kmer.output_kmer_main.len(),
            main_kmer.species_calling.display(),
        );
        schedule_rgi_artifact_copies(
            config.clone(), &main_kmer.output_kmer_main, &mut cleanup_tasks,
        );
    }

    let results_result = match (&main_result, &main_kmer_result, &kmer_result, &coverage_result) {
        (Ok(main), Some(Ok(main_species)), Ok(kma_species), Some(Ok(coverage))) => {
            Some(run_results_per_sample(config.clone(), main, main_species,
                                        &rgi_bwt_kma, kma_species, coverage, &report_ontology, &prepared.sample_base).await)
        }
        _ => None,
    };
    if let Some(Ok(results)) = &results_result {
        info!("[AMR] RunResultsPerSample complete: summary={}, primary={}",
            results.final_summary.display(), results.synthesized_report.display());
        schedule_rgi_artifact_copies(config.clone(), &results.artifacts, &mut cleanup_tasks);
    }

    let sam_result = match &results_result {
        Some(Ok(results)) => Some(tsv_to_sam(config.clone(), &prepared.contigs_fa, results).await),
        _ => None,
    };
    if let Some(Ok(sam)) = &sam_result {
        info!("[AMR] tsvToSam complete: bam={}, bai={}", sam.output_sorted.display(), sam.output_sorted_bai.display());
        schedule_rgi_artifact_copies(config.clone(), &[sam.output_sorted.clone(), sam.output_sorted_bai.clone()], &mut cleanup_tasks);
    }

    // Future downstream analysis belongs ABOVE final cleanup. It reads the
    // scratch paths, never waits for publication copies, and must retain owners.

    let mut publication_error = None;
    for task in cleanup_tasks {
        let result = match task.await {
            Ok(result) => result,
            Err(e) => Err(anyhow!("RGI artifact cleanup task panicked: {}", e)),
        };
        if let Err(e) = result {
            error!("[AMR] RGI artifact publication failed: {}", e);
            if publication_error.is_none() {
                publication_error = Some(e);
            }
        }
    }

    // Propagate analysis failures only after their independent writers finish.
    let rgi_kmer_bwt = kmer_result?;
    let rgi_main = main_result?;
    let gene_coverage = coverage_result
        .ok_or_else(|| PipelineError::Other(anyhow!("MakeGeneCoverage was not scheduled")))??;
    let rgi_kmer_main = main_kmer_result
        .ok_or_else(|| PipelineError::Other(anyhow!("RunRgiKmerMain was not scheduled")))??;
    let results = results_result
        .ok_or_else(|| PipelineError::Other(anyhow!("RunResultsPerSample was not scheduled")))??;
    let sam = sam_result.ok_or_else(|| PipelineError::Other(anyhow!("tsvToSam was not scheduled")))??;
    if let Some(e) = publication_error {
        return Err(PipelineError::Other(e));
    }
    drop(sam);
    drop(results);
    drop(rgi_kmer_main);
    drop(gene_coverage);
    drop(rgi_main);
    drop(rgi_kmer_bwt);
    drop(rgi_bwt_kma);

    Ok(())
}



#[cfg(test)]
mod amr_results_tests {
    use super::*;

    fn table(columns: &[&str], rows: &[&[&str]]) -> AmrResultTable {
        AmrResultTable {
            columns: columns.iter().map(|s| (*s).to_owned()).collect(),
            rows: rows.iter().map(|row| columns.iter().zip(row.iter()).map(|(k,v)|
                ((*k).to_owned(), if v.is_empty() { None } else { Some((*v).to_owned()) })).collect()).collect(),
        }
    }

    #[test]
    fn outer_join_preserves_cartesian_duplicates_and_unmatched_rows() {
        let a = table(&["a", "left"], &[&["x", "1"], &["x", "2"], &["z", "3"]]);
        let b = table(&["b", "right"], &[&["x", "4"], &["x", "5"], &["y", "6"]]);
        let joined = a.outer(&b, &["a"], &["b"]).unwrap();
        assert_eq!(joined.rows.len(), 6);
        assert_eq!(joined.rows.iter().filter(|r| amr_cell(r, "a") == Some("x")).count(), 4);
    }

    #[test]
    fn accession_stringification_retains_pandas_nan_and_float_promotion() {
        let mut t = table(&["id"], &[&["3001"], &[""]]);
        t.accession_strings("id").unwrap();
        assert_eq!(amr_cell(&t.rows[0], "id"), Some("3001.0"));
        assert_eq!(amr_cell(&t.rows[1], "id"), Some("nan"));
    }

    fn fixture(root: &std::path::Path, populated: bool, variant: bool) -> [PathBuf; 5] {
        let inputs = ["main.tsv", "main_sp.tsv", "kma.tsv", "kma_sp.tsv", "coverage.tsv"].map(|s| root.join(s));
        let main_cols = ["Best_Hit_ARO", "Contig", "ARO", "AMR Gene Family", "Drug Class", "Model_type",
            "Resistance Mechanism", "Cut_Off", "Nudged", "Percentage Length of Reference Sequence", "Best_Identities"];
        let main_sp_cols = ["Best_Hit_ARO", "Contig", "CARD*kmer Prediction"];
        let kma_cols = ["ARO Term", "ARO Accession", "Reference Model Type", "AMR Gene Family", "Drug Class",
            "Resistance Mechanism", "All Mapped Reads", "Percent Coverage", "Depth", "Reference Sequence"];
        let kma_sp_cols = ["ARO term", "CARD*kmer Prediction"];
        let main_row: &[&str] = &[" gene ", " ctg ", "3001", "family", "class", "protein homolog model", "mechanism", "Strict", "True", "90", "95"];
        let main_sp_row: &[&str] = &["gene", "ctg", "Species A: 2"];
        let kma_row: &[&str] = &["gene", "3001", if variant { "protein variant model" } else { "protein homolog model" },
            "family", "class", "mechanism", "20", "80", "3", "ref"];
        let kma_sp_row: &[&str] = &["gene", "Species A: 3"];
        let main_rows = [main_row];
        let main_sp_rows = [main_sp_row];
        let kma_rows = [kma_row];
        let kma_sp_rows = [kma_sp_row];
        let coverage_row: &[&str] = &["3001", "90.0"];
        let coverage_rows = [coverage_row];

        let data = [
            table(
                &main_cols,
                if populated { &main_rows } else { &[] },
            ),
            table(
                &main_sp_cols,
                if populated { &main_sp_rows } else { &[] },
            ),
            table(
                &kma_cols,
                if populated { &kma_rows } else { &[] },
            ),
            table(
                &kma_sp_cols,
                if populated { &kma_sp_rows } else { &[] },
            ),
            table(
                &["ID", "gene_coverage_perc"],
                if populated { &coverage_rows } else { &[] },
            ),
        ];
        for (t, path) in data.iter().zip(&inputs) { t.write(path, None).unwrap(); }
        std::fs::write(root.join("ontology.json"), r#"{"class":{"highLevelDrugClasses":["high class"]}}"#).unwrap();
        inputs
    }

    #[test]
    fn populated_reports_include_ontology_nudge_and_metrics() {
        let dir = tempfile::tempdir().unwrap();
        let inputs = fixture(dir.path(), true, false);
        amr_results_per_sample_files(&inputs, &dir.path().join("ontology.json"), "sample", dir.path()).unwrap();
        let report = AmrResultTable::read(&dir.path().join("primary_AMR_report.tsv")).unwrap();
        assert_eq!(report.rows.len(), 1);
        let r = &report.rows[0];
        assert_eq!(amr_cell(r, "gene_name"), Some("gene"));
        assert_eq!(amr_cell(r, "cutoff"), Some("Nudged"));
        assert_eq!(amr_cell(r, "high_level_drug_class"), Some("high class"));
        assert_eq!(amr_cell(r, "num_reads"), Some("20.0"));
        assert_eq!(amr_cell(r, "contig_coverage_breadth"), Some("90.0"));
    }

    #[test]
    fn empty_primary_is_zero_bytes_as_in_wdl() {
        let dir = tempfile::tempdir().unwrap();
        let inputs = fixture(dir.path(), false, false);
        amr_results_per_sample_files(&inputs, &dir.path().join("missing_ontology.json"), "sample", dir.path()).unwrap();
        assert_eq!(std::fs::metadata(dir.path().join("primary_AMR_report.tsv")).unwrap().len(), 0);
        assert!(dir.path().join("comprehensive_AMR_metrics.tsv").is_file());
    }

    #[test]
    fn kma_variant_rows_are_removed_before_final_merge() {
        let dir = tempfile::tempdir().unwrap();
        let inputs = fixture(dir.path(), true, true);
        amr_results_per_sample_files(&inputs, &dir.path().join("ontology.json"), "sample", dir.path()).unwrap();
        let b = AmrResultTable::read(&dir.path().join("merge_b.tsv")).unwrap();
        assert!(b.rows.is_empty());
        let primary = AmrResultTable::read(&dir.path().join("primary_AMR_report.tsv")).unwrap();
        assert_eq!(primary.rows.len(), 1);
        assert_eq!(amr_cell(&primary.rows[0], "num_reads"), None);
    }
}

#[cfg(test)]
mod gene_coverage_tests {
    use super::*;
    const HEADER: &str = "ORF_ID\tID\tARO\tModel_ID\tHit_Start\tHit_End\tPercentage Length of Reference Sequence\n";

    #[test]
    fn empty_report() {
        assert_eq!(gene_coverage_table(HEADER, &serde_json::json!({})).unwrap(),
                   "ID\tgene_coverage_bps\tdb_seq_length\tgene_coverage_perc\n");
    }

    #[test]
    fn overlapping_forward_and_reverse_ranges() {
        let tsv = format!("{}a\th\t10\t20\t0\t5\t0\nb\th\t10\t20\t7\t3\t0\n", HEADER);
        let json = serde_json::json!({
            "a": {"h": {"dna_sequence_from_broadstreet": "AAAAAAAAAA"}},
            "b": {"h": {"dna_sequence_from_broadstreet": "AAAAAAAAAA"}}
        });
        // range(0,5) U range(7,3,-1) = {0,1,2,3,4,5,6,7}.
        assert!(gene_coverage_table(&tsv, &json).unwrap().ends_with("10\t8\t10\t80.0\n"));
    }

    #[test]
    fn missing_json_hit_is_an_error() {
        let tsv = format!("{}a\th\t10\t20\t0\t5\t0\n", HEADER);
        assert!(gene_coverage_table(&tsv, &serde_json::json!({})).is_err());
    }
}

#[cfg(test)]
mod tsv_to_sam_tests {
    use super::*;
    #[test]
    fn strips_only_terminal_underscore_digits() {
        assert_eq!(amr_bam_contig_name("NODE_cov_1.23_4"), "NODE_cov_1.23");
        assert_eq!(amr_bam_contig_name("NODE_cov_1.23"), "NODE_cov_1.23");
        assert_eq!(amr_bam_contig_name("NODE_"), "NODE");
    }
    #[test]
    fn assembly_failure_has_dummy_reference() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("contigs.fa");
        let sam = dir.path().join("out.sam");
        std::fs::write(&fa, ";ASSEMBLY FAILED\n").unwrap();
        amr_write_contig_sam(&fa, &dir.path().join("absent.tsv"), &sam).unwrap();
        assert_eq!(std::fs::read_to_string(sam).unwrap(), "@SQ\tSN:NoGenes\tLN:100\n");
    }
    #[test]
    fn preserves_duplicates_and_pysam_position() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("contigs.fa");
        let summary = dir.path().join("summary.tsv");
        let sam = dir.path().join("out.sam");
        std::fs::write(&fa, ">NODE_cov_1.2 description\nACGT\n").unwrap();
        std::fs::write(&summary, "ARO_contig_amr\tContig_contig_amr\n3001.0\tNODE_cov_1.2_1\n3001.0\tNODE_cov_1.2_1\n\t\n").unwrap();
        amr_write_contig_sam(&fa, &summary, &sam).unwrap();
        let text = std::fs::read_to_string(sam).unwrap();
        assert_eq!(text.lines().count(), 3);
        assert!(text.contains("NODE_cov_1.2\t0\tARO:3001\t2\t0\t*\t*\t0\t0\tACGT\t*"));
    }
}
