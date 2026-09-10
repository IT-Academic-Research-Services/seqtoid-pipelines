//! Worker-side execution of one distributed NR work unit.
//!
//! This module executes a single claimed work unit using one of the supported
//! NR backends:
//!
//!   MMseqs CPU
//!   MMseqs GPU
//!   Diamond
//!
//! Queue discovery and durable claim coordination are intentionally outside
//! this module.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use log::{debug, info};
use tokio::fs;
use tokio::io::{AsyncBufReadExt, BufReader};

use crate::config::defs::{DiamondSubcommand, MMSEQS_TAG, DIAMOND_TAG};
use crate::utils::command::diamond::{
    generate_diamond_args,
    DiamondConfig,
    DiamondExecutionConfig,
};
use crate::utils::command::mmseqs::{
    generate_mmseqs_args,
    MmseqsBackend,
    MmseqsConfig,
    MmseqsExecutionConfig,
    MmseqsSubcommand,
};
use crate::utils::fastx::write_combined_fastq;
use crate::utils::streams::spawn_external_cmd;
use crate::utils::work_units::{
    WorkUnit,
    WorkUnitResult,
};

/// Backend selected for one worker execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkerBackend {
    MmseqsCpu,
    MmseqsGpu,
    Diamond,
}

/// Configuration required by the worker executor.
///
/// This is deliberately independent of the launch-node `RunConfig`.
#[derive(Debug, Clone)]
pub struct WorkerExecutorConfig {
    /// EC2 instance ID identifying this worker.
    pub worker_id: String,

    /// Local ephemeral scratch directory.
    pub scratch_dir: PathBuf,

    /// Durable EFS result directory.
    pub results_dir: PathBuf,

    /// Number of threads assigned to the execution.
    pub threads: usize,

    /// Maximum worker cores, used by the shared NUMA policy.
    pub max_cores: usize,

    /// Worker-local MMseqs database.
    pub mmseqs_db: Option<PathBuf>,

    /// Worker-local Diamond database prefix.
    pub diamond_db: Option<PathBuf>,
}

/// Executes one distributed work unit.
#[derive(Debug, Clone)]
pub struct WorkerExecutor {
    config: WorkerExecutorConfig,
    backend: WorkerBackend,
}

impl WorkerExecutor {
    pub fn new(
        config: WorkerExecutorConfig,
        backend: WorkerBackend,
    ) -> Self {
        Self { config, backend }
    }

    pub fn worker_id(&self) -> &str {
        &self.config.worker_id
    }

    /// Claims and executes one work unit, persisting every lifecycle
    /// transition to the WorkUnit JSON on EFS.
    pub async fn claim_and_execute(
        &self,
        work_unit_path: &Path,
        work_unit: &mut WorkUnit,
    ) -> Result<WorkUnitResult> {
        let attempt = work_unit
            .claim(self.worker_id())
            .map_err(|e| anyhow!(
                "failed to claim work unit {}: {}",
                work_unit.id(),
                e
            ))?;

        self.persist_work_unit(work_unit_path, work_unit).await?;

        if let Err(err) = self.validate_backend_reference(work_unit).await {
            let reason = err.to_string();
            work_unit
                .fail(self.worker_id(), attempt, reason.clone(), true)
                .map_err(|state_err| anyhow!(
                    "work unit {} reference validation failed with '{}', and FAILED transition also failed: {}",
                    work_unit.id(), reason, state_err
                ))?;
            self.persist_work_unit(work_unit_path, work_unit).await?;
            return Err(anyhow!(
                "work unit {} attempt {} failed: {}",
                work_unit.id(), attempt, reason
            ));
        }

        work_unit
            .start(self.worker_id(), attempt)
            .map_err(|e| anyhow!(
                "failed to start work unit {}: {}",
                work_unit.id(),
                e
            ))?;

        self.persist_work_unit(work_unit_path, work_unit).await?;

        info!(
            "[worker:{}] START work unit={} attempt={} backend={:?}",
            self.worker_id(),
            work_unit.id(),
            attempt,
            self.backend
        );

        match self.execute_attempt(work_unit, attempt).await {
            Ok(result) => {
                work_unit
                    .complete(self.worker_id(), attempt, result.clone())
                    .map_err(|e| anyhow!(
                        "execution succeeded but completion transition failed for {}: {}",
                        work_unit.id(), e
                    ))?;

                self.persist_work_unit(work_unit_path, work_unit).await?;
                self.publish_completion_metadata(work_unit).await?;

                info!(
                    "[worker:{}] DONE work unit={} attempt={} result={} bytes={} rows={}",
                    self.worker_id(), work_unit.id(), attempt,
                    result.result_path.display(), result.result_bytes, result.result_rows
                );
                Ok(result)
            }

            Err(err) => {
                let reason = err.to_string();
                work_unit
                    .fail(self.worker_id(), attempt, reason.clone(), true)
                    .map_err(|state_err| anyhow!(
                        "work unit {} failed with '{}', and FAILED transition also failed: {}",
                        work_unit.id(), reason, state_err
                    ))?;

                self.persist_work_unit(work_unit_path, work_unit).await?;
                let _ = self.publish_completion_metadata(work_unit).await;

                Err(anyhow!(
                    "work unit {} attempt {} failed: {}",
                    work_unit.id(), attempt, reason
                ))
            }
        }
    }

    /// Persist the current WorkUnit state atomically: write a complete JSON
    /// document beside the original file, then rename it over the original.
    async fn persist_work_unit(
        &self,
        work_unit_path: &Path,
        work_unit: &WorkUnit,
    ) -> Result<()> {
        let payload = serde_json::to_vec_pretty(work_unit)
            .context("failed to serialize WorkUnit")?;

        let parent = work_unit_path.parent().ok_or_else(|| {
            anyhow!("work unit path has no parent: {}", work_unit_path.display())
        })?;

        fs::create_dir_all(parent).await?;

        let temp_path = work_unit_path.with_extension("json.tmp");
        fs::write(&temp_path, payload).await.with_context(|| {
            format!("failed to write temporary WorkUnit {}", temp_path.display())
        })?;

        fs::rename(&temp_path, work_unit_path).await.with_context(|| {
            format!(
                "failed to publish WorkUnit state {} -> {}",
                temp_path.display(), work_unit_path.display()
            )
        })?;

        Ok(())
    }

    async fn execute_attempt(
        &self,
        work_unit: &WorkUnit,
        attempt: u32,
    ) -> Result<WorkUnitResult> {
        self.validate_input(work_unit).await?;

        let attempt_dir = self
            .config
            .scratch_dir
            .join("seqtoid-worker")
            .join(&work_unit.run_id)
            .join(&work_unit.sample_id)
            .join(format!("{:08}", work_unit.chunk_id))
            .join(format!("attempt_{:03}", attempt));

        fs::create_dir_all(&attempt_dir)
            .await
            .with_context(|| {
                format!(
                    "failed to create worker scratch directory {}",
                    attempt_dir.display()
                )
            })?;

        match self.backend {
            WorkerBackend::MmseqsCpu => {
                self.execute_mmseqs(
                    work_unit,
                    attempt,
                    &attempt_dir,
                    MmseqsBackend::Cpu,
                )
                    .await
            }

            WorkerBackend::MmseqsGpu => {
                self.execute_mmseqs(
                    work_unit,
                    attempt,
                    &attempt_dir,
                    MmseqsBackend::Gpu,
                )
                    .await
            }

            WorkerBackend::Diamond => {
                self.execute_diamond(
                    work_unit,
                    attempt,
                    &attempt_dir,
                )
                    .await
            }
        }
    }

    async fn execute_mmseqs(
        &self,
        work_unit: &WorkUnit,
        attempt: u32,
        attempt_dir: &Path,
        backend: MmseqsBackend,
    ) -> Result<WorkUnitResult> {
        let target_db = self
            .config
            .mmseqs_db
            .clone()
            .ok_or_else(|| {
                anyhow!(
                    "MMseqs backend selected but no worker-local MMseqs DB configured"
                )
            })?;

        let query_fastq = attempt_dir.join("query.fastq");
        let query_db = attempt_dir.join("queryDB");
        let result_db = attempt_dir.join("resultDB");
        let tmp_dir = attempt_dir.join("tmp");
        let output_m8 = attempt_dir.join("result.m8");

        fs::create_dir_all(&tmp_dir).await?;

        // Preserve the existing pipeline boundary:
        // paired R1/R2 -> combined FASTQ -> MMseqs createdb.
        write_combined_fastq(
            work_unit.r1_path.clone(),
            work_unit.r2_path.clone(),
            &query_fastq,
        )
            .await
            .context("failed to create combined FASTQ for MMseqs")?;

        let execution = MmseqsExecutionConfig {
            threads: self.config.threads,
            default_target_db: Some(target_db.clone()),
        };

        // ------------------------------------------------------------
        // createdb
        // ------------------------------------------------------------
        let createdb_config = MmseqsConfig {
            subcommand: MmseqsSubcommand::Createdb,
            backend: MmseqsBackend::Cpu,
            input: Some(query_fastq),
            target_db: None,
            result_db: None,
            output: Some(query_db.clone()),
            tmp_dir: None,
            threads: Some(self.config.threads),
            sensitivity: None,
            search_type: None,
            max_seqs: None,
            prefilter_mode: None,
            db_load_mode: None,
            alignment_mode: None,
            index_subset: None,
            format_output: None,
            cuda_visible_devices: None,
            option_fields: std::collections::HashMap::new(),
            gpu_server: false,
        };

        let args = generate_mmseqs_args(
            &execution,
            &createdb_config,
        )?;

        self.run_external(
            MMSEQS_TAG,
            args,
            attempt_dir,
            "createdb",
        )
            .await?;

        self.require_mmseqs_db(
            &query_db,
            "MMseqs query database",
        )
            .await?;

        // ------------------------------------------------------------
        // search
        // ------------------------------------------------------------
        let search_config = MmseqsConfig {
            subcommand: MmseqsSubcommand::Search,
            backend,
            input: Some(query_db.clone()),
            target_db: Some(target_db.clone()),
            result_db: Some(result_db.clone()),
            output: None,
            tmp_dir: Some(tmp_dir.clone()),
            threads: Some(self.config.threads),

            // Same production MMseqs parameters used by the existing
            // single-machine implementation.
            sensitivity: Some("5.7".to_string()),
            search_type: Some("3".to_string()),
            max_seqs: Some("1000".to_string()),
            prefilter_mode: Some("0".to_string()),
            db_load_mode: Some("2".to_string()),
            alignment_mode: Some("3".to_string()),
            index_subset: None,
            format_output: None,
            cuda_visible_devices: None,
            option_fields: std::collections::HashMap::from([
                (
                    "-e".to_string(),
                    Some("0.001".to_string()),
                ),
                (
                    "--min-seq-id".to_string(),
                    Some("0.25".to_string()),
                ),
            ]),
            gpu_server: false,
        };

        let args = generate_mmseqs_args(
            &execution,
            &search_config,
        )?;

        self.run_external(
            MMSEQS_TAG,
            args,
            attempt_dir,
            "search",
        )
            .await?;

        self.require_mmseqs_db(
            &result_db,
            "MMseqs search result database",
        )
            .await?;

        // ------------------------------------------------------------
        // convertalis
        // ------------------------------------------------------------
        let convert_config = MmseqsConfig {
            subcommand: MmseqsSubcommand::ConvertAlis,
            backend: MmseqsBackend::Cpu,
            input: Some(query_db),
            target_db: Some(target_db),
            result_db: Some(result_db),
            output: Some(output_m8.clone()),
            tmp_dir: None,
            threads: None,
            sensitivity: None,
            search_type: None,
            max_seqs: None,
            prefilter_mode: None,
            db_load_mode: None,
            alignment_mode: None,
            index_subset: None,
            format_output: Some(
                "query,target,pident,alnlen,mismatch,gapopen,\
                 qstart,qend,tstart,tend,evalue,bits"
                    .replace(' ', ""),
            ),
            cuda_visible_devices: None,
            option_fields: std::collections::HashMap::new(),
            gpu_server: false,
        };

        let args = generate_mmseqs_args(
            &execution,
            &convert_config,
        )?;

        self.run_external(
            MMSEQS_TAG,
            args,
            attempt_dir,
            "convertalis",
        )
            .await?;

        self.publish_result(
            work_unit,
            attempt,
            &output_m8,
        )
            .await
    }

    async fn execute_diamond(
        &self,
        work_unit: &WorkUnit,
        attempt: u32,
        attempt_dir: &Path,
    ) -> Result<WorkUnitResult> {
        let diamond_db = self
            .config
            .diamond_db
            .clone()
            .ok_or_else(|| {
                anyhow!(
                    "Diamond backend selected but no worker-local Diamond DB configured"
                )
            })?;

        let execution = DiamondExecutionConfig {
            threads: self.config.threads,
        };

        let mut output_parts = Vec::new();

        // Diamond operates independently on R1 and R2 when paired-end input
        // is provided. The resulting m8 files are concatenated into the
        // durable per-chunk result.
        let r1_out = attempt_dir.join("diamond_R1.m8");

        let r1_config = DiamondConfig {
            subcommand: DiamondSubcommand::Blastx,
            db: diamond_db.clone(),
            r1_path: Some(work_unit.r1_path.clone()),
            r2_path: None,
            subcommand_fields: std::collections::HashMap::from([
                (
                    "--query".to_string(),
                    Some(
                        work_unit
                            .r1_path
                            .to_string_lossy()
                            .to_string(),
                    ),
                ),
                (
                    "--out".to_string(),
                    Some(
                        r1_out
                            .to_string_lossy()
                            .to_string(),
                    ),
                ),
                (
                    "-f".to_string(),
                    Some("6".to_string()),
                ),
            ]),
        };

        let args = generate_diamond_args(
            &execution,
            &r1_config,
        )?;

        self.run_external(
            DIAMOND_TAG,
            args,
            attempt_dir,
            "blastx_R1",
        )
            .await?;

        self.require_path(
            &r1_out,
            "Diamond R1 m8",
        )
            .await?;

        output_parts.push(r1_out);

        if let Some(r2) = &work_unit.r2_path {
            let r2_out = attempt_dir.join("diamond_R2.m8");

            let r2_config = DiamondConfig {
                subcommand: DiamondSubcommand::Blastx,
                db: diamond_db,
                r1_path: Some(r2.clone()),
                r2_path: None,
                subcommand_fields: std::collections::HashMap::from([
                    (
                        "--query".to_string(),
                        Some(
                            r2.to_string_lossy().to_string(),
                        ),
                    ),
                    (
                        "--out".to_string(),
                        Some(
                            r2_out
                                .to_string_lossy()
                                .to_string(),
                        ),
                    ),
                    (
                        "-f".to_string(),
                        Some("6".to_string()),
                    ),
                ]),
            };

            let args = generate_diamond_args(
                &execution,
                &r2_config,
            )?;

            self.run_external(
                DIAMOND_TAG,
                args,
                attempt_dir,
                "blastx_R2",
            )
                .await?;

            self.require_path(
                &r2_out,
                "Diamond R2 m8",
            )
                .await?;

            output_parts.push(r2_out);
        }

        let merged_m8 = attempt_dir.join("result.m8");

        let mut merged =
            fs::File::create(&merged_m8)
                .await
                .context("failed to create merged Diamond m8")?;

        for part in output_parts {
            let mut input =
                fs::File::open(&part).await?;

            tokio::io::copy(
                &mut input,
                &mut merged,
            )
                .await?;
        }

        drop(merged);

        self.publish_result(
            work_unit,
            attempt,
            &merged_m8,
        )
            .await
    }

    /// Shared external command execution.
    ///
    /// Uses the repository-wide spawn_external_cmd() implementation so the
    /// worker and launch-node pipeline share the same process-launching and
    /// NUMA behavior.
    async fn run_external(
        &self,
        cmd_tag: &str,
        args: Vec<String>,
        work_dir: &Path,
        label: &str,
    ) -> Result<()> {
        let stderr_log = work_dir.join(
            format!("{}.stderr.log", label),
        );

        debug!(
            "[worker:{}] spawning {}: {:?}",
            self.worker_id(),
            cmd_tag,
            args
        );

        let (mut child, stderr_task) =
            spawn_external_cmd(
                cmd_tag,
                args,
                self.config.max_cores,
                false,
                Some(stderr_log.clone()),
            )
                .await
                .with_context(|| {
                    format!(
                        "failed to spawn {} for {}",
                        cmd_tag, label
                    )
                })?;

        let status = child.wait().await
            .with_context(|| {
                format!(
                    "failed waiting for {} {}",
                    cmd_tag, label
                )
            })?;

        stderr_task
            .await
            .context("stderr handler task failed")??;

        if !status.success() {
            return Err(anyhow!(
                "{} {} exited with status {:?}; stderr={}",
                cmd_tag,
                label,
                status.code(),
                stderr_log.display()
            ));
        }

        info!(
            "[worker:{}] {} {} completed successfully",
            self.worker_id(),
            cmd_tag,
            label
        );

        Ok(())
    }

    async fn validate_input(
        &self,
        work_unit: &WorkUnit,
    ) -> Result<()> {
        if !work_unit.r1_path.is_file() {
            return Err(anyhow!(
                "work unit R1 does not exist: {}",
                work_unit.r1_path.display()
            ));
        }

        match (
            work_unit.paired_end,
            work_unit.r2_path.as_ref(),
        ) {
            (true, Some(r2)) if r2.is_file() => {}

            (true, None) => {
                return Err(anyhow!(
                    "work unit {} is paired-end but has no R2 path",
                    work_unit.id()
                ));
            }

            (true, Some(r2)) => {
                return Err(anyhow!(
                    "work unit R2 does not exist: {}",
                    r2.display()
                ));
            }

            (false, Some(_)) => {
                return Err(anyhow!(
                    "work unit {} is single-end but has an R2 path",
                    work_unit.id()
                ));
            }

            (false, None) => {}
        }

        Ok(())
    }

    async fn validate_backend_reference(&self, work_unit: &WorkUnit) -> Result<()> {
        let (db, version_path) = match self.backend {
            WorkerBackend::MmseqsCpu => (
                self.config.mmseqs_db.as_ref().ok_or_else(|| {
                    anyhow!("MMseqs CPU backend selected without worker-local MMseqs DB")
                })?,
                self.config.scratch_dir.join("refs/mmseqs/.reference_version"),
            ),
            WorkerBackend::MmseqsGpu => (
                self.config.mmseqs_db.as_ref().ok_or_else(|| {
                    anyhow!("MMseqs GPU backend selected without worker-local MMseqs DB")
                })?,
                self.config.scratch_dir.join("refs/mmseqs-gpu/.reference_version"),
            ),
            WorkerBackend::Diamond => (
                self.config.diamond_db.as_ref().ok_or_else(|| {
                    anyhow!("Diamond backend selected without worker-local Diamond DB")
                })?,
                self.config.scratch_dir.join("refs/diamond/.reference_version"),
            ),
        };

        if !db.exists() {
            return Err(anyhow!(
                "worker-local reference DB does not exist: {}",
                db.display()
            ));
        }

        let local_version = fs::read_to_string(&version_path)
            .await
            .with_context(|| {
                format!(
                    "failed to read worker-local reference version {}",
                    version_path.display()
                )
            })?
            .trim()
            .to_string();

        if local_version.is_empty() {
            return Err(anyhow!(
                "worker-local reference version is empty: {}",
                version_path.display()
            ));
        }

        if local_version != work_unit.reference_version {
            return Err(anyhow!(
                "reference version mismatch for {}: work unit requires {}, worker has {}",
                work_unit.id(), work_unit.reference_version, local_version
            ));
        }

        info!(
            "[worker:{}] reference version validated: {}",
            self.worker_id(), local_version
        );

        Ok(())
    }

    async fn require_path(
        &self,
        path: &Path,
        description: &str,
    ) -> Result<()> {
        if !path.exists() {
            return Err(anyhow!(
                "{} was not produced: {}",
                description,
                path.display()
            ));
        }

        Ok(())
    }

    /// Validate an MMseqs database prefix.
    ///
    /// MMseqs databases may be represented by a prefix plus sidecar files and
    /// split parts, so the prefix itself does not necessarily exist as a file.
    async fn require_mmseqs_db(
        &self,
        path: &Path,
        description: &str,
    ) -> Result<()> {
        if path.is_file() {
            return Ok(());
        }

        let dbtype = path.with_extension("dbtype");
        if !dbtype.is_file() {
            return Err(anyhow!(
                "{} was not produced: missing MMseqs dbtype file {}",
                description,
                dbtype.display()
            ));
        }

        if path.with_extension("index").is_file() {
            return Ok(());
        }

        let parent = path.parent().unwrap_or_else(|| Path::new("."));
        let prefix = path.file_name().and_then(|v| v.to_str()).unwrap_or_default();
        let numbered_prefix = format!("{}.", prefix);

        let mut entries = fs::read_dir(parent).await.with_context(|| {
            format!("failed to inspect MMseqs database directory {}", parent.display())
        })?;

        while let Some(entry) = entries.next_entry().await? {
            let name = entry.file_name();
            let name = name.to_string_lossy();
            if let Some(suffix) = name.strip_prefix(&numbered_prefix) {
                if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
                    return Ok(());
                }
            }
        }

        Err(anyhow!(
            "{} was not produced: MMseqs database prefix {} has no .index or numbered split parts",
            description,
            path.display()
        ))
    }

    /// Validates and atomically publishes a completed m8 result.
    async fn publish_result(
        &self,
        work_unit: &WorkUnit,
        attempt: u32,
        temporary_m8: &Path,
    ) -> Result<WorkUnitResult> {
        let validation = validate_m8(temporary_m8).await?;

        fs::create_dir_all(&self.config.results_dir).await?;

        let final_name = format!(
            "chunk_{:08}_attempt_{:03}.m8",
            work_unit.chunk_id,
            attempt
        );

        let final_path = self.config.results_dir.join(final_name);

        let temp_name = format!(
            "chunk_{:08}_attempt_{:03}.m8.tmp",
            work_unit.chunk_id,
            attempt
        );

        let temp_path = self.config.results_dir.join(temp_name);

        fs::copy(temporary_m8, &temp_path)
            .await
            .with_context(|| {
                format!(
                    "failed to copy m8 {} -> {}",
                    temporary_m8.display(),
                    temp_path.display()
                )
            })?;

        // Publication is now atomic because both paths are on EFS.
        fs::rename(&temp_path, &final_path)
            .await
            .with_context(|| {
                format!(
                    "failed to publish m8 {} -> {}",
                    temp_path.display(),
                    final_path.display()
                )
            })?;

        let metadata = fs::metadata(&final_path).await?;

        Ok(WorkUnitResult {
            worker_id: self.worker_id().to_string(),
            attempt,
            result_path: final_path,
            result_bytes: metadata.len(),
            result_rows: validation.rows,
            checksum: None,
        })
    }

    async fn publish_completion_metadata(
        &self,
        work_unit: &WorkUnit,
    ) -> Result<()> {
        fs::create_dir_all(&self.config.results_dir)
            .await?;

        let metadata_name = format!(
            "chunk_{:08}_attempt_{:03}.json",
            work_unit.chunk_id,
            work_unit.attempt
        );

        let metadata_path =
            self.config.results_dir.join(metadata_name);

        let temporary_path =
            metadata_path.with_extension("json.tmp");

        let payload =
            serde_json::to_vec_pretty(work_unit)
                .context(
                    "failed to serialize work-unit metadata",
                )?;

        fs::write(
            &temporary_path,
            payload,
        )
            .await?;

        fs::rename(
            &temporary_path,
            &metadata_path,
        )
            .await
            .with_context(|| {
                format!(
                    "failed to publish completion metadata {} -> {}",
                    temporary_path.display(),
                    metadata_path.display()
                )
            })?;

        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
struct M8Validation {
    rows: u64,
}

/// Validate the generated 12-column tabular m8.
///
/// A zero-row result is valid.
async fn validate_m8(path: &Path) -> Result<M8Validation> {
    let file = fs::File::open(path)
        .await
        .with_context(|| {
            format!(
                "failed to open generated m8 {}",
                path.display()
            )
        })?;

    let mut reader =
        BufReader::new(file).lines();

    let mut rows = 0u64;

    while let Some(line) =
        reader.next_line().await?
    {
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> =
            line.split('\t').collect();

        if fields.len() != 12 {
            return Err(anyhow!(
                "invalid m8 row in {}: expected 12 columns, found {}",
                path.display(),
                fields.len()
            ));
        }

        rows += 1;
    }

    Ok(M8Validation { rows })
}