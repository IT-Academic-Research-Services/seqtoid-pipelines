/// Functions and structs for working with creating command-line arguments

use anyhow::{anyhow, Result};
use log::{self, LevelFilter, debug, info, error, warn};
use num_cpus;
use tokio::process::Command;
use futures::future::try_join_all;
use crate::config::defs::{RunConfig, PipelineError, TOOL_VERSIONS, FASTP_TAG, PIGZ_TAG, H5DUMP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, IVAR_TAG, MUSCLE_TAG, MAFFT_TAG, QUAST_TAG, NUCMER_TAG, SHOW_COORDS_TAG, SEQKIT_TAG, BOWTIE2_TAG, HISAT2_TAG, KALLISTO_TAG, STAR_TAG, FASTA_EXTS, CZID_DEDUP_TAG, DIAMOND_TAG, SPADES_TAG, BLASTN_TAG, BLASTX_TAG, MAKEBLASTDB_TAG, SORT_TAG, MMSEQS_TAG};
use crate::cli::Arguments;
use crate::utils::streams::{read_child_output_to_vec, ChildStream};
use std::path::PathBuf;
use std::sync::Arc;
use tokio::fs;
use tokio::task::JoinHandle;
use tempfile::{NamedTempFile, TempDir};
use fxhash::FxHashMap;
use crate::utils::file::{write_vecu8_to_file, extension_remover};
use crate::utils::streams::{spawn_cmd};


pub trait ArgGenerator {
    fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>>;
}

/// Checks the version of a CLI external tool.
/// NB: ONLY keeps tracks of major plus minor, i.e. 2.1, 4.3, not a point release like 3.4.5
/// Strips all non digit parts of the version.
/// # Arguments
///
/// * `command_tag` - Executable name of tool.
/// * `version_args` - CLI args necessary to get tool to output version.
/// * `version_line` - Which output line contains the version.
/// * `version_column` - Which output column contains the version.
/// * `child_stream` - Is the verion on stdout or stderr.
///
/// # Returns
/// Result<f32>major/minor version number
pub async fn version_check(
    command_tag: &str,
    version_args: Vec<&str>,
    version_line: usize,
    version_column: usize,
    child_stream: ChildStream,
    version_file: Option<PathBuf>,
    config: &RunConfig
) -> Result<f32> {
    let cmd_tag_owned = command_tag.to_string();
    debug!("Running command: {}", &cmd_tag_owned);
    let args: Vec<&str> = version_args;

    let mut child = Command::new(&cmd_tag_owned)
        .args(&args)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow!("Failed to spawn {}. Is it installed? Error: {}.", cmd_tag_owned.clone(), e))?;

    let lines = read_child_output_to_vec(&mut child, child_stream, &config).await?;

    if let Some(file_path) = version_file {
        if let Some(parent) = file_path.parent() {
            fs::create_dir_all(parent).await
                .map_err(|e| anyhow!("Failed to create directory for version file: {}", e))?;
        }
        let full_output = lines.join("\n");
        fs::write(&file_path, full_output.as_bytes()).await
            .map_err(|e| anyhow!("Failed to write version file {:?}: {}", file_path, e))?;
        info!("Wrote version output to {:?}", file_path);
    }

    let line_w_version = lines
        .get(version_line)
        .ok_or_else(|| anyhow!("No line {} in {} version output", version_line, cmd_tag_owned.clone()))?;

    let version_string = line_w_version
        .split_whitespace()
        .nth(version_column)
        .ok_or_else(|| anyhow!("Invalid {} version output: {}", cmd_tag_owned.clone(), line_w_version))?;

    if version_string.is_empty() {
        return Err(anyhow!("Empty version number string in {} version output: {}", cmd_tag_owned.clone(), line_w_version));
    }

    let version_parts: Vec<_> = version_string.split(".").collect();
    let major_version = version_parts[0];
    let major_version_digits: String = major_version.chars().filter(|c| c.is_digit(10)).collect();
    let mut minor_version = "0";
    if version_parts.len() > 1 {
        minor_version = version_parts[1];
    }
    let minor_version_digits: String = minor_version.chars().filter(|c| c.is_digit(10)).collect();
    let version_string_formatted = format!("{}.{}", major_version_digits, minor_version_digits);
    let version = version_string_formatted.parse::<f32>()?;
    Ok(version)
}

/// Prepends `numactl --interleave=all` for large multi-socket Linux machines (EPYC).
/// Does nothing on MacBooks, small machines, or non-Linux.
pub fn prepend_numactl_if_beneficial(config: &RunConfig, mut args: Vec<String>) -> Vec<String> {
    if cfg!(target_os = "linux") && config.max_cores >= 64 {
        let mut numa_args = vec![
            "--interleave=all".to_string(),
            MMSEQS_TAG.to_string(),   // or whatever binary you're calling
        ];
        numa_args.extend(args);
        debug!("Prepended numactl --interleave=all for large EPYC machine");
        numa_args
    } else {
        args
    }
}

pub mod fastp {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use log::{self, LevelFilter, debug, info, error, warn};
    use crate::config::defs::{FASTP_TAG,  RunConfig};
    use crate::utils::streams::{ChildStream};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::file::file_path_manipulator;

    #[derive(Debug)]
    pub struct FastpConfig {
        pub command_fields: HashMap<String, Option<String>>,
        pub paired: bool,
    }
    pub struct FastpArgGenerator;

    pub async fn fastp_presence_check(config: &RunConfig) -> Result<f32> {
        let version = version_check(FASTP_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for FastpArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            debug!("Allocating {} threads for fastp", RunConfig::thread_allocation(run_config, FASTP_TAG, None));

            let config = extra
                .and_then(|e| e.downcast_ref::<FastpConfig>())
                .ok_or_else(|| anyhow!("FASTP requires a FaspConfig as extra argument"))?;

            let args = &run_config.args;
            let mut args_vec: Vec<String> = Vec::new();
            args_vec.push("--stdin".to_string());
            args_vec.push("--stdout".to_string());
            args_vec.push("--interleaved_in".to_string());
            // args_vec.push("-q".to_string());
            // args_vec.push(args.quality.to_string());

            let json_out = run_config.out_dir.join("fastp.json");
            args_vec.push("-j".to_string());
            args_vec.push(json_out.as_os_str().to_str().unwrap().to_string());
            let html_out = run_config.out_dir.join("fastp.html");
            args_vec.push("-h".to_string());
            args_vec.push(html_out.as_os_str().to_str().unwrap().to_string());
            args_vec.push("-w".to_string());
            args_vec.push(RunConfig::thread_allocation(run_config, FASTP_TAG, None).to_string());

            if let Some(adapter_fasta) = &args.adapter_fasta {
                let cwd = std::env::current_dir()?;
                let adapter_path = file_path_manipulator(&PathBuf::from(adapter_fasta), Some(&cwd.clone()), None, None, "");
                if !adapter_path.exists() {
                    return Err(anyhow!("Adapter FASTA file does not exist: {}", adapter_path.display()));
                }
                args_vec.push("--adapter_fasta".to_string());
                args_vec.push(adapter_path.to_string_lossy().into_owned());
            }

            if args.adapter_fasta.is_some() && config.paired {
                args_vec.push("--detect_adapter_for_pe".to_string());
            }

            for (key, value) in config.command_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            Ok(args_vec)
        }
    }
}

mod pigz {
    use crate::config::defs::{FASTP_TAG, PIGZ_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    pub struct PigzArgGenerator;

    pub async fn pigz_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(PIGZ_TAG,vec!["--version"], 0, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for PigzArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, _extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let mut args_vec: Vec<String> = Vec::new();
            args_vec.push("-c".to_string());
            args_vec.push("-p".to_string());
            args_vec.push(RunConfig::thread_allocation(run_config, PIGZ_TAG, None).to_string());
            Ok(args_vec)
        }
    }
}

mod h5dump {
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{H5DUMP_TAG, RunConfig};
    use crate::utils::command::version_check;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub async fn h5dump_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(H5DUMP_TAG,vec!["-V"], 0, 2 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }
}

pub mod minimap2 {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use std::sync::Arc;
    use tokio::fs;
    use tokio::task::JoinHandle;
    use tempfile::{NamedTempFile, TempDir};
    use log::{self, LevelFilter, debug, info, error, warn};
    use anyhow::{anyhow, Result};
    use crate::config::defs::{RunConfig, PipelineError, MINIMAP2_TAG, FASTA_EXTS};
    use crate::utils::file::available_space_for_path;
    use crate::utils::command::{ArgGenerator, version_check, spawn_cmd};
    use crate::utils::streams::ChildStream;

    pub struct Minimap2Config {
        pub minimap2_index_path: PathBuf,
        pub r1_path: Option<PathBuf>,             // ← NEW: optional query file (None = stdin)
        pub r2_path: Option<PathBuf>,             // ← NEW: optional query file (None = stdin)
        pub option_fields: HashMap<String, Option<String>>,
        pub num_threads: Option<usize>,
    }

    pub struct Minimap2ArgGenerator;

    pub async fn minimap2_presence_check(config: &RunConfig, version_file: Option<PathBuf>) -> Result<f32> {
        let version = version_check(MINIMAP2_TAG, vec!["--version"], 0, 0, ChildStream::Stdout, version_file, &config).await?;
        Ok(version)
    }

    pub async fn minimap2_index_prep(
        config: &RunConfig,
        ram_temp_dir: &PathBuf,
        sequence: Option<String>,
        index_path: Option<String>,
        ref_type: &str,
    ) -> Result<
        (
            Option<PathBuf>,             // FASTA path
            PathBuf,                     // index path (.mmi)
            Option<NamedTempFile>,       // FASTA temp
            Option<NamedTempFile>,       // Index temp
            Option<TempDir>,             // temp dir for splits
            Vec<JoinHandle<Result<(), anyhow::Error>>>,
        ),
        PipelineError,
    > {
        let mut tasks = Vec::new();
        let mut ref_fasta_path: Option<PathBuf> = None;
        let mut ref_temp: Option<NamedTempFile> = None;
        let mut final_index_path: PathBuf = PathBuf::new();
        let mut index_temp_file: Option<NamedTempFile> = None;
        let mut index_temp_dir: Option<TempDir> = None;

        match index_path {
            // pre-existing .mmi index, maybe w split parts
            Some(index) => {
                let orig_index_path = PathBuf::from(&index);
                if !orig_index_path.exists() {
                    return Err(PipelineError::FileNotFound(orig_index_path));
                }
                if orig_index_path.extension().map_or(true, |ext| ext != "mmi") {
                    return Err(PipelineError::InvalidConfig(format!(
                        "{} index must have .mmi extension: {}",
                        ref_type,
                        orig_index_path.display()
                    )));
                }

                let orig_dir = orig_index_path.parent().ok_or_else(|| {
                    PipelineError::InvalidConfig("Invalid index path".to_string())
                })?;
                let basename = orig_index_path.file_stem().ok_or_else(|| {
                    PipelineError::InvalidConfig("Invalid index filename".to_string())
                })?.to_str().ok_or_else(|| {
                    PipelineError::InvalidConfig("Invalid UTF-8 in filename".to_string())
                })?;

                // find split parts: {basename}.part_001.mmi, etc.
                let mut split_files: Vec<PathBuf> = Vec::new();
                let split_prefix = format!("{}.part_", basename);
                let mut dir_entries = fs::read_dir(orig_dir).await
                    .map_err(|e| PipelineError::Other(e.into()))?;
                while let Some(entry) = dir_entries.next_entry().await
                    .map_err(|e| PipelineError::Other(e.into()))? {
                    let path = entry.path();
                    if path.is_file() {
                        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                            if name.starts_with(&split_prefix) && name.ends_with(".mmi") {
                                split_files.push(path);
                            }
                        }
                    }
                }

                let is_split = !split_files.is_empty();

                let mut files_to_copy: Vec<(PathBuf, String)> = Vec::new();
                let mmi_name = orig_index_path.file_name().unwrap().to_string_lossy().to_string();
                files_to_copy.push((orig_index_path.clone(), mmi_name.clone()));
                for split in &split_files {
                    let name = split.file_name().unwrap().to_string_lossy().to_string();
                    files_to_copy.push((split.clone(), name));
                }

                // calc size of all index files
                let mut total_size: u64 = 0;
                for (p, _) in &files_to_copy {
                    let meta = fs::metadata(p).await
                        .map_err(|e| PipelineError::Other(e.into()))?;
                    total_size += meta.len();
                }

                let buffer = 10 * 1024 * 1024; // 10 megs buffer
                let avail = available_space_for_path(ram_temp_dir)
                    .await
                    .map_err(|e| PipelineError::Other(e))?;
                let enough_space = avail >= total_size + buffer;

                if enough_space {
                    if is_split {
                        let temp_dir = TempDir::new_in(ram_temp_dir)
                            .map_err(|e| PipelineError::Other(e.into()))?;
                        let temp_index_path = temp_dir.path().join(&mmi_name);

                        for (src, name) in files_to_copy {
                            let dest = temp_dir.path().join(&name);
                            let src_cl = src.clone();
                            let dest_cl = dest.clone();
                            let ref_type_cl = ref_type.to_string();
                            let task = tokio::spawn(async move {
                                fs::copy(&src_cl, &dest_cl).await
                                    .map_err(|e| anyhow!("Failed to copy {} index file {}: {}", ref_type_cl, name, e))?;
                                Ok(())
                            });
                            tasks.push(task);
                        }

                        final_index_path = temp_index_path;
                        index_temp_dir = Some(temp_dir);
                    } else {
                        let temp_file = NamedTempFile::with_suffix_in(".mmi", ram_temp_dir)
                            .map_err(|e| PipelineError::Other(e.into()))?;
                        let temp_index_path = temp_file.path().to_path_buf();
                        let temp_index_path_cl = temp_index_path.clone();
                        let ref_type_cl = ref_type.to_string();
                        let task = tokio::spawn(async move {
                            fs::copy(&orig_index_path, &temp_index_path_cl).await
                                .map_err(|e| anyhow!("Failed to copy {} index: {}", ref_type_cl, e))?;
                            Ok(())
                        });
                        tasks.push(task);
                        final_index_path = temp_index_path;
                        index_temp_file = Some(temp_file);
                    }
                } else {
                    info!(
                        "Not enough space in RAM temp dir for {} index (need: {} bytes, have: {} bytes). Using original path.",
                        ref_type, total_size, avail
                    );
                    final_index_path = orig_index_path;
                }
            }

            // or build index from FASTA
            None => {
                let sequence_path = sequence.ok_or_else(|| {
                    PipelineError::InvalidConfig("No FASTA or index provided for minimap2".to_string())
                })?;
                let sequence_path = PathBuf::from(&sequence_path);
                if !sequence_path.exists() {
                    return Err(PipelineError::FileNotFound(sequence_path));
                }
                let is_fasta = sequence_path
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .map(|ext| FASTA_EXTS.iter().any(|&f| ext.eq_ignore_ascii_case(f)))
                    .unwrap_or(false);
                if !is_fasta {
                    return Err(PipelineError::InvalidConfig(format!(
                        "Sequence file must have a FASTA extension ({}): {}",
                        FASTA_EXTS.join(", "),
                        sequence_path.display()
                    )));
                }

                let temp_fasta = NamedTempFile::with_suffix_in(".fasta", ram_temp_dir)
                    .map_err(|e| PipelineError::Other(e.into()))?;
                let temp_fasta_path = temp_fasta.path().to_path_buf();
                let temp_fasta_path_cl = temp_fasta_path.clone();
                let ref_type_cl = ref_type.to_string();
                let copy_task = tokio::spawn(async move {
                    fs::copy(&sequence_path, &temp_fasta_path_cl).await
                        .map_err(|e| anyhow!("Failed to copy {} FASTA: {}", ref_type_cl, e))?;
                    Ok(())
                });
                tasks.push(copy_task);
                ref_fasta_path = Some(temp_fasta_path.clone());
                ref_temp = Some(temp_fasta);

                let index_temp = NamedTempFile::with_suffix_in(".mmi", ram_temp_dir)
                    .map_err(|e| PipelineError::Other(e.into()))?;
                let index_temp_path = index_temp.path().to_path_buf();

                let minimap2_args = vec![
                    "-d".to_string(),
                    index_temp_path.to_string_lossy().to_string(),
                    temp_fasta_path.to_string_lossy().to_string(),
                ];

                let (mut child, err_task) = spawn_cmd(
                    Arc::new(config.clone()),
                    MINIMAP2_TAG,
                    minimap2_args,
                    config.args.verbose,
                    None
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
                index_task.await.map_err(|e| PipelineError::Other(e.into()))??;
                debug!("{} mmi created in RAM", ref_type);
                tasks.push(err_task);

                final_index_path = index_temp_path;
                index_temp_file = Some(index_temp);
            }
        }

        Ok((
            ref_fasta_path,
            final_index_path,
            ref_temp,
            index_temp_file,
            index_temp_dir,
            tasks,
        ))
    }

    impl ArgGenerator for Minimap2ArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<Minimap2Config>())
                .ok_or_else(|| anyhow!("Minimap2 requires a Minimap2Config as extra argument"))?;

            let index_path = config.minimap2_index_path.clone();
            if !index_path.exists() {
                return Err(anyhow!("Minimap2 index file does not exist: {}", index_path.display()));
            }
            if index_path.extension().map_or(true, |ext| ext != "mmi") {
                return Err(anyhow!("Minimap2 index must have .mmi extension, got: {}", index_path.display()));
            }

            let mut args_vec: Vec<String> = Vec::new();

            let mut threads = if let Some(override_threads) = config.num_threads {
                override_threads
            } else {
                run_config.thread_allocation(MINIMAP2_TAG, None)
            };
            
            // let num_cores: usize = RunConfig::thread_allocation(run_config, MINIMAP2_TAG, None);
            args_vec.push("-t".to_string());
            args_vec.push(threads.to_string());

            for (key, value) in config.option_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            // Index path
            args_vec.push(index_path.to_string_lossy().to_string());

            // Input: file if provided, otherwise stdin
            if let Some(r1_path) = &config.r1_path {
                args_vec.push(r1_path.to_string_lossy().to_string());
                if let Some(r2_path) = &config.r2_path { args_vec.push(r2_path.to_string_lossy().to_string());}

            } else {
                args_vec.push("-".to_string()); // stdin
            }

            Ok(args_vec)
        }
    }
}

pub mod samtools {
    use std::collections::HashMap;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{SAMTOOLS_TAG, SamtoolsSubcommand, RunConfig};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::{version_check, ArgGenerator};

    #[derive(Debug)]
    pub struct SamtoolsConfig {
        pub subcommand: SamtoolsSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct SamtoolsArgGenerator;

    pub async fn samtools_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(SAMTOOLS_TAG, vec!["--version"], 0, 1, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for SamtoolsArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let config = extra
                .and_then(|e| e.downcast_ref::<SamtoolsConfig>())
                .ok_or_else(|| anyhow!("Samtools requires a SamtoolsConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand {
                SamtoolsSubcommand::View => {
                    args_vec.push("view".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("view")).to_string());
                    args_vec.push("--no-PG".to_string());
                }
                SamtoolsSubcommand::Fastq => {
                    args_vec.push("fastq".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("fastq")).to_string());
                    args_vec.push("-c".to_string());
                    args_vec.push("6".to_string());
                    args_vec.push("-n".to_string());
                }
                SamtoolsSubcommand::Stats => {
                    args_vec.push("stats".to_string());
                }
                SamtoolsSubcommand::Sort => {
                    args_vec.push("sort".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("sort")).to_string());
                    args_vec.push("-m".to_string());
                    args_vec.push("2G".to_string());
                }
                SamtoolsSubcommand::Index => {
                    args_vec.push("index".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("index")).to_string());
                }
                SamtoolsSubcommand::Mpileup => {
                    args_vec.push("mpileup".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("mpileup")).to_string());
                }
                SamtoolsSubcommand::Consensus => {
                    args_vec.push("consensus".to_string());
                    args_vec.push("-f".to_string());
                    args_vec.push("fasta".to_string());
                    args_vec.push("-m".to_string());
                    args_vec.push("simple".to_string());
                    args_vec.push("--min-BQ".to_string());
                    args_vec.push(args.quality.to_string());
                    args_vec.push("-d".to_string());
                    args_vec.push(args.min_depth.to_string());
                }
                SamtoolsSubcommand::Depth => {
                    args_vec.push("depth".to_string());
                }
                SamtoolsSubcommand::Ampliconclip => {
                    args_vec.push("ampliconclip".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("ampliconclip")).to_string());
                }
                SamtoolsSubcommand::Fixmate => {
                    args_vec.push("fixmate".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("fixmate")).to_string());
                }
            }
            for (key, value) in config.subcommand_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            Ok(args_vec)
        }
    }
}

pub mod bcftools {
    use std::collections::HashMap;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{BcftoolsSubcommand, BCFTOOLS_TAG, RunConfig, SAMTOOLS_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[derive(Debug)]
    pub struct BcftoolsConfig {
        pub subcommand: BcftoolsSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct BcftoolsArgGenerator;

    pub async fn bcftools_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(BCFTOOLS_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for BcftoolsArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let config = extra
                .and_then(|e| e.downcast_ref::<BcftoolsConfig>())
                .ok_or_else(|| anyhow!("Bcftools requires a BcftoolsConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand {
                BcftoolsSubcommand::Consensus => {
                    args_vec.push("consensus".to_string());
                }
                BcftoolsSubcommand::Call => {
                    args_vec.push("call".to_string());
                }
                BcftoolsSubcommand::Mpileup => {
                    args_vec.push("mpileup".to_string());
                    args_vec.push("-a".to_string());
                    args_vec.push("AD".to_string()); // include allele depth (AD) for all positions, including those with zero coverage, mimicking -aa.
                    args_vec.push("-d".to_string());
                    args_vec.push("10000".to_string()); // max depth essentially without limit
                    args_vec.push("-L".to_string());
                    args_vec.push("100000000".to_string()); // max per-file depth essentially without limit
                    args_vec.push("-Q".to_string());
                    args_vec.push(args.quality.to_string()); // skip bases with baseQ/BAQ smaller than INT [13]
                    args_vec.push("--threads".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, BCFTOOLS_TAG, Some("mpileup")).to_string());

                }
                BcftoolsSubcommand::View => {
                    args_vec.push("view".to_string());
                }
            }

            for (key, value) in config.subcommand_fields.iter() {
                args_vec.push(format!("{}", key));
                match value {
                    Some(v) => args_vec.push(format!("{}", v)),
                    None => { },
                }
            }

            Ok(args_vec)
        }
    }
}

pub mod kraken2 {
    use anyhow::anyhow;
    use std::path::PathBuf;
    use tokio::process::Command;
    use crate::config::defs::{KRAKEN2_TAG, RunConfig};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::file::file_path_manipulator;

    #[derive(Debug)]
    pub struct Kraken2Config {
        pub report_path: PathBuf,
        pub classified_path: PathBuf,
        pub fastq_path: PathBuf,
    }

    pub struct Kraken2ArgGenerator;

    pub async fn kraken2_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(KRAKEN2_TAG,vec!["--version"], 0, 2 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for Kraken2ArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let config = extra
                .and_then(|e| e.downcast_ref::<Kraken2Config>())
                .ok_or_else(|| anyhow!("Kraken2 requires a Kraken2Config as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();
            let cwd = std::env::current_dir()?;
            match &args.kraken_db{
                Some(db) => {
                    let kraken2_db_path = file_path_manipulator(&PathBuf::from(db), Some(&cwd), None, None, "");
                    if !kraken2_db_path.exists() || !kraken2_db_path.is_dir() {
                        return Err(anyhow!("Kraken2 database path does not exist or is not a directory: {:?}", kraken2_db_path));
                    }

                    args_vec.push("--db".to_string());
                    args_vec.push(kraken2_db_path.to_string_lossy().to_string());
                }
                None => {
                    return Err(anyhow!("No kraken_db specified"));
                }
            }

            let num_cores: usize = RunConfig::thread_allocation(run_config, KRAKEN2_TAG, None);
            args_vec.push("--threads".to_string());
            args_vec.push(num_cores.to_string());
            args_vec.push("--report".to_string());
            args_vec.push(config.report_path.to_string_lossy().to_string());
            args_vec.push("--classified-out".to_string());
            args_vec.push(config.classified_path.to_string_lossy().to_string());
            args_vec.push("--output".to_string());
            args_vec.push("-".to_string()); // "-" will suppress normal output

            args_vec.push(config.fastq_path.to_string_lossy().to_string());

            Ok(args_vec)
        }
    }
}

pub mod ivar {
    use std::collections::HashMap;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::RunConfig;
    use crate::config::defs::{IVAR_TAG, IvarSubcommand, IVAR_QUAL_THRESHOLD, IVAR_FREQ_THRESHOLD};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[derive(Debug)]
    pub struct IvarConfig {
        pub subcommand: IvarSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct IvarArgGenerator;

    pub async fn ivar_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(IVAR_TAG,vec!["version"], 0, 2 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for IvarArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let config = extra
                .and_then(|e| e.downcast_ref::<IvarConfig>())
                .ok_or_else(|| anyhow!("IVar requires a IvarConfig as extra argument"))?;
            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand {
                IvarSubcommand::Consensus => {
                    args_vec.push("consensus".to_string());
                    args_vec.push("-q".to_string());
                    args_vec.push(IVAR_QUAL_THRESHOLD.to_string());
                    args_vec.push("-t".to_string());
                    args_vec.push(IVAR_FREQ_THRESHOLD.to_string());
                    args_vec.push("-m".to_string());
                    args_vec.push(args.min_depth.to_string());
                    args_vec.push("-n N".to_string());
                }
            }

            for (key, value) in config.subcommand_fields.iter() {
                args_vec.push(format!("{}", key));
                match value {
                    Some(v) => args_vec.push(format!("{}", v)),
                    None => { },
                }
            }
            Ok(args_vec)
        }
    }
}

pub mod muscle {
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{RunConfig,  MUSCLE_TAG};
    use crate::utils::command::version_check;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub struct MuscleArgGenerator;
    pub async fn muscle_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(MUSCLE_TAG,vec!["-version"], 0, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }
}

pub mod mafft {
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{MAFFT_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub struct MafftArgGenerator;
    pub async fn mafft_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(MAFFT_TAG,vec!["--version"], 0, 0 , ChildStream::Stderr, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for MafftArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, _extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let mut args_vec: Vec<String> = Vec::new();

            let num_cores: usize = RunConfig::thread_allocation(run_config, MAFFT_TAG, None);
            args_vec.push("--auto".to_string());
            args_vec.push("--thread".to_string());
            args_vec.push(num_cores.to_string());
            args_vec.push("-".to_string());

            Ok(args_vec)
        }
    }
}

pub mod quast {
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{QUAST_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub struct QuastArgGenerator;

    pub struct QuastConfig {
        pub ref_fasta: String,
        pub ref_bam: String,
        pub assembly_fasta: String,
    }

    pub async fn quast_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(QUAST_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for QuastArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let config = extra
                .and_then(|e| e.downcast_ref::<QuastConfig>())
                .ok_or_else(|| anyhow!("Quast requires a QuastConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            let num_cores: usize = RunConfig::thread_allocation(run_config, QUAST_TAG, None);
            args_vec.push("--min-contig".to_string());
            args_vec.push("0".to_string());
            args_vec.push("-t".to_string());
            args_vec.push(num_cores.to_string());
            args_vec.push("-o".to_string());

            let quast_out = run_config.out_dir.join("quast");
            args_vec.push(quast_out.as_os_str().to_str().unwrap().to_string());
            args_vec.push("-r".to_string());
            args_vec.push(config.ref_fasta.to_string());
            args_vec.push("--ref-bam".to_string());
            args_vec.push(config.ref_bam.to_string());
            args_vec.push(config.assembly_fasta.to_string());

            Ok(args_vec)
        }
    }
}

pub mod nucmer {
    use anyhow::anyhow;
    use std::path::PathBuf;
    use tokio::process::Command;
    use crate::config::defs::{NUCMER_DELTA, NUCMER_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub struct NucmerArgGenerator;

    pub struct NucmerConfig {
        pub ref_fasta: PathBuf,
        pub assembly_fasta: PathBuf,
    }

    pub async fn nucmer_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(NUCMER_TAG,vec!["--version"], 0, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for NucmerArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<NucmerConfig>())
                .ok_or_else(|| anyhow!("Nucmer requires a NucmerConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            let num_cores: usize = RunConfig::thread_allocation(run_config, NUCMER_TAG, None);
            args_vec.push("-t".to_string());
            args_vec.push(num_cores.to_string());
            args_vec.push("--prefix=alignment".to_string());
            args_vec.push(config.ref_fasta.to_string_lossy().to_string());
            args_vec.push(config.assembly_fasta.to_string_lossy().to_string());

            Ok(args_vec)
        }
    }
}

pub mod show_coords {
    use anyhow::anyhow;
    use crate::config::defs::RunConfig;
    use crate::config::defs::{SHOW_COORDS_TAG, NUCMER_DELTA};
    use crate::utils::command::ArgGenerator;

    pub struct ShowCoordsArgGenerator;

    impl ArgGenerator for ShowCoordsArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let mut args_vec: Vec<String> = Vec::new();

            args_vec.push("-r".to_string());
            args_vec.push("-c".to_string());
            args_vec.push(NUCMER_DELTA.to_string());

            Ok(args_vec)
        }
    }
}

pub mod seqkit {
    use std::collections::HashMap;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::RunConfig;
    use crate::config::defs::{SeqkitSubcommand, SEQKIT_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[derive(Debug)]
    pub struct SeqkitConfig {
        pub subcommand: SeqkitSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct SeqkitArgGenerator;

    pub async fn seqkit_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(SEQKIT_TAG,vec!["--help"], 1, 1 , ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for SeqkitArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<SeqkitConfig>())
                .ok_or_else(|| anyhow!("Seqkit requires a SeqkitConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand { 
                SeqkitSubcommand::Stats => {
                    args_vec.push("stats".to_string());

                }
                
                SeqkitSubcommand::Rmdup => {
                    args_vec.push("rmdup".to_string());

                }
                
            }

            for (key, value) in config.subcommand_fields.iter() {
                args_vec.push(format!("{}", key));
                match value {
                    Some(v) => args_vec.push(format!("{}", v)),
                    None => { },
                }
            }

            args_vec.push("-".to_string());

            Ok(args_vec)
        }
    }
}

pub mod bowtie2 {
    use std::path::{Path, PathBuf};
    use std::collections::HashMap;
    use anyhow::{anyhow, Result};
    use std::fs::{self, DirEntry};
    use std::process::Command;
    use log::{self, LevelFilter, debug, info, error, warn};
    use crate::config::defs::{RunConfig, BOWTIE2_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    pub struct Bowtie2Config {
        pub bt2_index_path: PathBuf,
        pub paired: bool,
        pub r1_path: Option<PathBuf>,
        pub r2_path: Option<PathBuf>,
        pub option_fields: HashMap<String, Option<String>>,
    }

    pub struct Bowtie2ArgGenerator;

    pub async fn bowtie2_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(BOWTIE2_TAG, vec!["-h"], 0, 3, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    /// Prepares the Bowtie2 index by handling directory, tar, or tar.gz inputs.
    /// Returns the PathBuf to the basename (e.g., /path/to/ercc for files ercc.1.bt2, etc.).
    pub fn bowtie2_index_prep(input_path: impl AsRef<Path>, cwd: &PathBuf) -> Result<PathBuf> {
        let input = input_path.as_ref();

        // Generate unique unpack directory
        let unique_id = if input.is_file() {
            input.file_stem().unwrap_or_default().to_string_lossy().to_string()
        } else {
            input.file_name().unwrap_or_default().to_string_lossy().to_string()
        };
        let unpack_dir = cwd.join(format!("unpacked_index_bt2_{}", unique_id));

        if unpack_dir.exists() {
            fs::remove_dir_all(&unpack_dir)?;
        }
        fs::create_dir_all(&unpack_dir)?;
        debug!("Bowtie2: Created unpack directory: {}", unpack_dir.display());

        let mut index_dir = unpack_dir.clone();

        if input.is_dir() {
            index_dir = input.to_path_buf();
        } else if input.exists() {
            let ext = input.extension().map(|s| s.to_str().unwrap_or("")).unwrap_or("");
            let stem_ext = input.file_stem().and_then(|s| Path::new(s).extension()).map(|s| s.to_str().unwrap_or("")).unwrap_or("");

            let unpack_cmd = if ext == "gz" && stem_ext == "tar" {
                debug!("Bowtie2: Unpacking GZIP TAR: {}", input.display());
                Command::new("tar")
                    .arg("xzf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else if ext == "tar" {
                debug!("Bowtie2: Unpacking TAR: {}", input.display());
                Command::new("tar")
                    .arg("xf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else {
                return Err(anyhow!("Unsupported Bowtie2 index file format: {}", input.display()));
            };

            if !unpack_cmd.status.success() {
                return Err(anyhow!("Failed to unpack Bowtie2 index: {:?}", String::from_utf8_lossy(&unpack_cmd.stderr)));
            }
        } else {
            return Err(anyhow!("Bowtie2 index input not found: {}", input.display()));
        }

        // Strip single top-level subdir (e.g., "ercc/")
        let entries: Vec<DirEntry> = fs::read_dir(&index_dir)?.collect::<Result<_, _>>()?;
        let dirs: Vec<DirEntry> = entries.into_iter().filter(|e| e.file_type().map(|ft| ft.is_dir()).unwrap_or(false)).collect();

        if dirs.len() == 1 {
            let subdir = dirs[0].path();
            if fs::read_dir(&subdir)?.next().is_some() {
                index_dir = subdir;
            }
        }
        debug!("Bowtie2 index directory: {}", index_dir.display());

        // Validate
        let extensions = vec!["bt2".to_string(), "bt2l".to_string()];
        let required_suffixes = vec!["1", "2", "3", "4", "rev.1", "rev.2"];

        for suffix in &required_suffixes {
            let mut found = false;
            for ext in &extensions {
                let pattern = format!("{}.{}", suffix, ext);
                let mut candidates = Vec::new();
                for entry in fs::read_dir(&index_dir)? {
                    let path = entry?.path();
                    if path.is_file() && path.to_str().unwrap_or("").ends_with(&pattern) {
                        candidates.push(path);
                    }
                }
                if !candidates.is_empty() {
                    found = true;
                    break;
                }
            }
            if !found {
                return Err(anyhow!("Missing Bowtie2 index file for suffix '{}' in: {}", suffix, index_dir.display()));
            }
        }

        // Derive basename using suffix "1", avoiding "rev"
        let pattern_suffix = "1";
        let mut basename = None;
        for ext in &extensions {
            let pattern = format!("{}.{}", pattern_suffix, ext);
            let mut candidates = Vec::new();
            for entry in fs::read_dir(&index_dir)? {
                let path = entry?.path();
                if path.is_file() && path.to_str().unwrap_or("").ends_with(&pattern) && !path.to_str().unwrap_or("").contains("rev") {
                    candidates.push(path);
                }
            }
            if let Some(file) = candidates.first() {
                let file_name = file.file_name().ok_or(anyhow!("Invalid file path"))?.to_str().ok_or(anyhow!("Invalid UTF-8"))?;
                let stripped = file_name.strip_suffix(&pattern).ok_or(anyhow!("Suffix mismatch in Bowtie2 file: {}", file_name))?;
                basename = Some(stripped.trim_end_matches('.').to_string());
                break;
            }
        }
        let basename = basename.ok_or(anyhow!("No Bowtie2 index file with suffix '1' avoiding 'rev' found in: {}", index_dir.display()))?;

        let final_path = index_dir.join(&basename);
        debug!("Final Bowtie2 index path: {}", final_path.display());
        Ok(final_path)
    }

    impl ArgGenerator for Bowtie2ArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<Bowtie2Config>())
                .ok_or_else(|| anyhow!("Bowtie2 requires a Bowtie2Config as extra argument"))?;

            let mut args: Vec<String> = Vec::new();

            for (key, value) in &config.option_fields {
                args.push(key.clone());
                if let Some(v) = value {
                    args.push(v.clone());
                }
            }

            // Index
            args.push("-x".to_string());
            args.push(config.bt2_index_path.to_string_lossy().to_string());

            // Threads
            args.push("-p".to_string());
            let num_cores = RunConfig::thread_allocation(run_config, BOWTIE2_TAG, None);
            args.push(num_cores.to_string());

            // Input handling
            match (config.r1_path.as_ref(), config.r2_path.as_ref()) {
                (Some(r1), Some(r2)) => {
                    // Paired files
                    args.push("-1".to_string());
                    args.push(r1.to_string_lossy().to_string());
                    args.push("-2".to_string());
                    args.push(r2.to_string_lossy().to_string());
                }
                (Some(r1), None) => {
                    // Single-end file
                    args.push("-U".to_string());
                    args.push(r1.to_string_lossy().to_string());
                }
                (None, _) => {
                    // Stdin mode
                    if config.paired {
                        args.push("--interleaved".to_string());
                    } else {
                        args.push("-U".to_string());
                    }
                    args.push("-".to_string());  // stdin
                }
            }

            // Safety check: paired flag should match input
            if config.paired && config.r1_path.is_none() && config.r2_path.is_none() {
                debug!("Bowtie2: using --interleaved for paired stdin input");
            } else if config.paired && (config.r1_path.is_none() || config.r2_path.is_none()) {
                return Err(anyhow!("Bowtie2 paired mode requires both R1 and R2 paths (or stdin)"));
            }

            debug!("Bowtie2 final args: {}", args.join(" "));

            Ok(args)
        }
    }
}

pub mod hisat2 {
    use std::path::{Path, PathBuf};
    use std::collections::HashMap;
    use anyhow::{anyhow, Result};
    use log::{self, LevelFilter, debug, info, error, warn};
    use std::fs::{self, DirEntry};
    use std::process::Command;
    use crate::config::defs::{RunConfig, HISAT2_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    #[derive(Debug)]
    pub struct Hisat2Config {
        pub hisat2_index_path: PathBuf,
        pub option_fields: HashMap<String, Option<String>>,
        pub r1_path: String,
        pub r2_path: Option<String>,
    }

    pub struct Hisat2ArgGenerator;

    pub async fn hisat2_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(HISAT2_TAG, vec!["--version"], 0, 2, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    /// Prepares the HISAT2 index by handling directory, tar, or tar.gz inputs.
    /// Returns the PathBuf to the basename (e.g., /path/to/human for files human.1.ht2, etc.).
    pub fn hisat2_index_prep(input_path: impl AsRef<Path>, cwd: &PathBuf) -> Result<PathBuf> {
        let input = input_path.as_ref();

        // Generate unique unpack directory
        let unique_id = if input.is_file() {
            input.file_stem().unwrap_or_default().to_string_lossy().to_string()
        } else {
            input.file_name().unwrap_or_default().to_string_lossy().to_string()
        };
        let unpack_dir = cwd.join(format!("unpacked_index_ht2_{}", unique_id));

        // Clean unpack_dir
        if unpack_dir.exists() {
            fs::remove_dir_all(&unpack_dir)?;
        }
        fs::create_dir_all(&unpack_dir)?;
        debug!("HISAT2: Created unpack directory: {}", unpack_dir.display());

        let mut index_dir = unpack_dir.clone();

        // Handle input type
        if input.is_dir() {
            index_dir = input.to_path_buf();
        } else if input.exists() {
            let ext = input.extension().map(|s| s.to_str().unwrap_or("")).unwrap_or("");
            let stem_ext = input.file_stem().and_then(|s| Path::new(s).extension()).map(|s| s.to_str().unwrap_or("")).unwrap_or("");

            let unpack_cmd = if ext == "gz" && stem_ext == "tar" {
                debug!("HISAT2: Unpacking GZIP TAR: {}", input.display());
                Command::new("tar")
                    .arg("xzf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else if ext == "tar" {
                debug!("HISAT2: Unpacking TAR: {}", input.display());
                Command::new("tar")
                    .arg("xf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else {
                return Err(anyhow!("Unsupported HISAT2 index file format: {}", input.display()));
            };

            if !unpack_cmd.status.success() {
                return Err(anyhow!("Failed to unpack HISAT2 index: {:?}", String::from_utf8_lossy(&unpack_cmd.stderr)));
            }
        } else {
            return Err(anyhow!("HISAT2 index input not found: {}", input.display()));
        }

        // Strip single top-level subdir (e.g., "human/")
        let entries: Vec<DirEntry> = fs::read_dir(&index_dir)?.collect::<Result<_, _>>()?;
        let dirs: Vec<DirEntry> = entries.into_iter().filter(|e| e.file_type().map(|ft| ft.is_dir()).unwrap_or(false)).collect();

        if dirs.len() == 1 {
            let subdir = dirs[0].path();
            if fs::read_dir(&subdir)?.next().is_some() {
                index_dir = subdir;
            }
        }

        // Validate required files
        let extensions = vec!["ht2".to_string(), "ht2l".to_string()];
        let required_suffixes = vec!["1", "2", "3", "4", "5", "6", "7", "8"];

        for suffix in &required_suffixes {
            let mut found = false;
            for ext in &extensions {
                let pattern = format!("{}.{}", suffix, ext);
                let mut candidates = Vec::new();
                for entry in fs::read_dir(&index_dir)? {
                    let path = entry?.path();
                    if path.is_file() && path.to_str().unwrap_or("").ends_with(&pattern) {
                        candidates.push(path);
                    }
                }
                if !candidates.is_empty() {
                    found = true;
                    break;
                }
            }
            if !found {
                return Err(anyhow!("Missing HISAT2 index file for suffix '{}' in: {}", suffix, index_dir.display()));
            }
        }

        // Derive basename using suffix "1", avoiding any unexpected patterns
        let pattern_suffix = "1";
        let mut basename = None;
        for ext in &extensions {
            let pattern = format!("{}.{}", pattern_suffix, ext);
            let mut candidates = Vec::new();
            for entry in fs::read_dir(&index_dir)? {
                let path = entry?.path();
                if path.is_file() && path.to_str().unwrap_or("").ends_with(&pattern) {
                    candidates.push(path);
                }
            }
            if let Some(file) = candidates.first() {
                let file_name = file.file_name().ok_or(anyhow!("Invalid file path"))?.to_str().ok_or(anyhow!("Invalid UTF-8"))?;
                let stripped = file_name.strip_suffix(&pattern).ok_or(anyhow!("Suffix mismatch in HISAT2 file: {}", file_name))?;
                basename = Some(stripped.trim_end_matches('.').to_string());
                break;
            }
        }
        let basename = basename.ok_or(anyhow!("No HISAT2 index file with suffix '1' found in: {}", index_dir.display()))?;

        let final_path = index_dir.join(&basename);
        debug!("Final HISAT2 index path: {}", final_path.display());
        Ok(final_path)
    }

    impl ArgGenerator for Hisat2ArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<Hisat2Config>())
                .ok_or_else(|| anyhow!("HISAT2 requires a Hisat2Config as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            for (key, value) in config.option_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            args_vec.push("-x".to_string());
            args_vec.push(config.hisat2_index_path.to_string_lossy().to_string());
            args_vec.push("-p".to_string());
            let num_cores: usize = RunConfig::thread_allocation(run_config, HISAT2_TAG, None);
            args_vec.push(num_cores.to_string());

            if config.r2_path.is_some() {
                args_vec.push("-1".to_string());
                args_vec.push(config.r1_path.clone());
                args_vec.push("-2".to_string());
                args_vec.push(config.r2_path.clone().unwrap());
            } else {
                args_vec.push("-U".to_string());
                args_vec.push(config.r1_path.clone());
            }


            Ok(args_vec)
        }
    }
}




pub mod kallisto {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::anyhow;
    use crate::config::defs::{RunConfig, KALLISTO_TAG, KallistoSubcommand};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    pub struct KallistoConfig {
        pub subcommand: KallistoSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
        pub output_dir: PathBuf,
        pub reproducible: bool,
    }

    pub struct KallistoArgGenerator;

    pub async fn kallisto_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(KALLISTO_TAG, vec!["version"], 0, 2, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for KallistoArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let config = extra
                .and_then(|e| e.downcast_ref::<KallistoConfig>())
                .ok_or_else(|| anyhow!("Kallisto requires a KallistoConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand {
                KallistoSubcommand::Index => {
                    args_vec.push("index".to_string());
                }
                KallistoSubcommand::Quant => {
                    args_vec.push("quant".to_string());
                    args_vec.push("-i".to_string());
                    args_vec.push(args.kallisto_index.clone().expect("Must provide kallisto index with --kallisto-index arg.").to_string());
                    args_vec.push("-o".to_string());
                    args_vec.push(config.output_dir.to_string_lossy().to_string());
                    args_vec.push("--plaintext".to_string());
                    if config.reproducible {
                        args_vec.push("--threads".to_string());
                        args_vec.push("1".to_string());
                        args_vec.push("--seed".to_string());
                        args_vec.push("42".to_string());
                    } else {
                        args_vec.push("--threads".to_string());
                        let num_cores: usize = RunConfig::thread_allocation(run_config, KALLISTO_TAG, None);
                        args_vec.push(num_cores.to_string());
                    }

                    // Add non-file options from subcommand_fields (e.g., --single, -l, -s)
                    for (key, value) in config.subcommand_fields.iter() {
                        if key != "R1" && key != "R2" {
                            args_vec.push(key.clone());
                            if let Some(v) = value {
                                args_vec.push(v.clone());
                            }
                        }
                    }

                    // Add R1 and R2 paths in order as positional arguments
                    if let Some(r1_path) = config.subcommand_fields.get("R1").and_then(|v| v.as_ref()) {
                        args_vec.push(r1_path.clone());
                    }
                    if let Some(r2_path) = config.subcommand_fields.get("R2").and_then(|v| v.as_ref()) {
                        args_vec.push(r2_path.clone());
                    }
                }
            }

            Ok(args_vec)
        }
    }
}

pub mod star {
    use std::collections::HashMap;
    use std::path::{Path, PathBuf};
    use log::{self, LevelFilter, debug, info, error, warn};
    use anyhow::{anyhow, Result};
    use std::fs::{self, DirEntry};
    use std::process::Command;
    use crate::config::defs::{RunConfig, STAR_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    pub struct StarConfig {
        pub star_index_dir: PathBuf,
        pub option_fields: HashMap<String, Option<String>>,
        pub r1_fifo: PathBuf,
        pub r2_fifo: Option<PathBuf>,
    }

    pub struct StarArgGenerator;

    pub async fn star_presence_check(config: &RunConfig)-> anyhow::Result<f32> {
        let version = version_check(STAR_TAG, vec!["--version"], 0, 0, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    /// Prepares the STAR index by handling directory, tar, or tar.gz inputs.
    /// Returns the PathBuf to the index directory (e.g., /path/to/star_human).
    pub fn star_index_prep(input_path: impl AsRef<Path>, cwd: &PathBuf) -> Result<PathBuf> {
        let input = input_path.as_ref();

        // Generate unique unpack directory
        let unique_id = if input.is_file() {
            input.file_stem().unwrap_or_default().to_string_lossy().to_string()
        } else {
            input.file_name().unwrap_or_default().to_string_lossy().to_string()
        };
        let unpack_dir = cwd.join(format!("unpacked_index_star_{}", unique_id));

        // Clean unpack_dir
        if unpack_dir.exists() {
            fs::remove_dir_all(&unpack_dir)?;
        }
        fs::create_dir_all(&unpack_dir)?;
        debug!("STAR: Created unpack directory: {}", unpack_dir.display());

        let mut index_dir = unpack_dir.clone();

        // Handle input type
        if input.is_dir() {
            index_dir = input.to_path_buf();
        } else if input.exists() {
            let ext = input.extension().map(|s| s.to_str().unwrap_or("")).unwrap_or("");
            let stem_ext = input.file_stem().and_then(|s| Path::new(s).extension()).map(|s| s.to_str().unwrap_or("")).unwrap_or("");

            let unpack_cmd = if ext == "gz" && stem_ext == "tar" {
                debug!("STAR: Unpacking GZIP TAR: {}", input.display());
                Command::new("tar")
                    .arg("xzf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else if ext == "tar" {
                debug!("STAR: Unpacking TAR: {}", input.display());
                Command::new("tar")
                    .arg("xf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else {
                return Err(anyhow!("Unsupported STAR index file format: {}", input.display()));
            };

            if !unpack_cmd.status.success() {
                return Err(anyhow!("Failed to unpack STAR index: {:?}", String::from_utf8_lossy(&unpack_cmd.stderr)));
            }
        } else {
            return Err(anyhow!("STAR index input not found: {}", input.display()));
        }

        // Strip single top-level subdir (e.g., "star_human/")
        let entries: Vec<DirEntry> = fs::read_dir(&index_dir)?.collect::<Result<_, _>>()?;
        let dirs: Vec<DirEntry> = entries.into_iter().filter(|e| e.file_type().map(|ft| ft.is_dir()).unwrap_or(false)).collect();

        if dirs.len() == 1 {
            let subdir = dirs[0].path();
            if fs::read_dir(&subdir)?.next().is_some() {
                index_dir = subdir;
            }
        }

        // Validate required files
        let required_files = vec![
            "SA".to_string(),
            "SAindex".to_string(),
            "Genome".to_string(),
            "genomeParameters.txt".to_string(),
        ];

        for file in &required_files {
            let path = index_dir.join(file);
            if !path.exists() {
                return Err(anyhow!("Missing STAR index file '{}' in: {}", file, index_dir.display()));
            }
        }

        Ok(index_dir)
    }

    impl ArgGenerator for StarArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<StarConfig>())
                .ok_or_else(|| anyhow!("STAR requires a StarConfig as extra argument"))?;

            // Validate FIFO paths
            if !config.r1_fifo.exists() {
                return Err(anyhow!("R1 FIFO does not exist: {}", config.r1_fifo.display()));
            }
            if let Some(r2_fifo) = &config.r2_fifo {
                if !r2_fifo.exists() {
                    return Err(anyhow!("R2 FIFO does not exist: {}", r2_fifo.display()));
                }
            }

            let mut args_vec: Vec<String> = Vec::new();

            args_vec.push("--genomeDir".to_string());
            args_vec.push(config.star_index_dir.to_string_lossy().to_string());

            args_vec.push("--runThreadN".to_string());
            let num_cores: usize = std::cmp::min(128, RunConfig::thread_allocation(run_config, STAR_TAG, None));
            args_vec.push(num_cores.to_string());

            args_vec.push("--readFilesIn".to_string());
            args_vec.push(config.r1_fifo.to_string_lossy().to_string());
            if let Some(r2_fifo) = &config.r2_fifo {
                args_vec.push(r2_fifo.to_string_lossy().to_string());
            }

            args_vec.push("--readFilesCommand".to_string());
            args_vec.push("cat".to_string());

            args_vec.push("--outSAMtype".to_string());
            args_vec.push("SAM".to_string());

            args_vec.push("--outStd".to_string());
            args_vec.push("SAM".to_string());

            args_vec.push("--outFileNamePrefix".to_string());
            args_vec.push("/tmp/star_".to_string());

            for (key, value) in config.option_fields.iter() {
                args_vec.push(format!("{}", key));
                if let Some(v) = value {
                    args_vec.push(format!("{}", v));
                }
            }

            Ok(args_vec)
        }
    }
}

pub mod czid_dedup {

    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{RunConfig, CZID_DEDUP_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    #[derive(Debug)]
    pub struct CzidDedupConfig {
        pub input_paths: Vec<PathBuf>,  // FIFOs or files for inputs (1 for single, 2 for paired)
        pub output_paths: Vec<PathBuf>, // FIFOs or files for outputs (matching inputs)
        pub prefix_length: Option<u32>,
        pub cluster_output: Option<PathBuf>,
        pub cluster_size_output: Option<PathBuf>,
    }

    pub struct CzidDedupArgGenerator;

    pub async fn czid_dedup_presence_check(config: &RunConfig)-> Result<f32> {
        let version = version_check(CZID_DEDUP_TAG, vec!["--help"], 0, 1, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for CzidDedupArgGenerator {
        fn generate_args(&self, _run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<CzidDedupConfig>())
                .ok_or_else(|| anyhow!("czid-dedup requires CzidDedupConfig as extra argument"))?;

            if config.input_paths.len() != config.output_paths.len() {
                return Err(anyhow!("Number of inputs must match outputs for czid-dedup"));
            }

            let mut args_vec: Vec<String> = Vec::new();

            for input in &config.input_paths {
                args_vec.push("--inputs".to_string());
                args_vec.push(input.to_string_lossy().to_string());
            }

            for output in &config.output_paths {
                args_vec.push("--deduped-outputs".to_string());
                args_vec.push(output.to_string_lossy().to_string());
            }

            if let Some(len) = config.prefix_length {
                args_vec.push("--prefix-length".to_string());
                args_vec.push(len.to_string());
            }

            if let Some(cluster) = &config.cluster_output {
                args_vec.push("--cluster-output".to_string());
                args_vec.push(cluster.to_string_lossy().to_string());
            }

            if let Some(size) = &config.cluster_size_output {
                args_vec.push("--cluster-size-output".to_string());
                args_vec.push(size.to_string_lossy().to_string());
            }

            Ok(args_vec)
        }
    }
}


pub mod diamond {

    const DIAMOND_TEMP:u64 = 20 * 1024 * 1024 * 1024; // 20 GB

    use std::collections::HashMap;
    use std::path::PathBuf;
    use std::sync::Arc;
    use std::path::Path;
    use regex::Regex;
    use log::{debug, info, error, warn};
    use tokio::fs::{self, DirEntry};
    use tokio::task::JoinHandle;
    use tokio::process::Command;
    use tempfile::{NamedTempFile, TempDir};
    use anyhow::{anyhow, Result as AnyhowResult};
    use crate::config::defs::{DIAMOND_TAG, DiamondSubcommand, RunConfig, PipelineError};
    use crate::utils::file::{available_space_for_path, choose_temp_dir};
    use crate::utils::system::detect_ram;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream, spawn_cmd};
    use crate::utils::command::{version_check, ArgGenerator};

    #[derive(Debug)]
    pub struct DiamondConfig {
        pub subcommand: DiamondSubcommand,
        pub db: PathBuf,
        pub r1_path: Option<PathBuf>,
        pub r2_path: Option<PathBuf>,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct DiamondArgGenerator;

    pub async fn diamond_presence_check(config: &RunConfig, version_file: Option<PathBuf>) -> anyhow::Result<f32> {
        let version = version_check(DIAMOND_TAG, vec!["help"], 0, 1, ChildStream::Stdout, version_file, &config).await?;
        Ok(version)
    }

    pub async fn diamond_index_prep(
        index_path: Option<String>,
        ref_type: &str,
    ) -> Result<
        (
            PathBuf,                     // DB prefix path (without .dmnd)
            Vec<JoinHandle<Result<(), anyhow::Error>>>,
        ),
        PipelineError,
    > {
        let dmnd_path = index_path.ok_or_else(|| {
            PipelineError::MissingArgument(format!(
                "{}: --diamond-db <PATH> is required (must be a .dmnd file)",
                ref_type
            ))
        })?;

        let path = PathBuf::from(&dmnd_path);

        if !path.exists() {
            return Err(PipelineError::FileNotFound(path));
        }
        if path.extension() != Some(std::ffi::OsStr::new("dmnd")) {
            return Err(PipelineError::InvalidConfig(format!(
                "{} index must be a .dmnd file: {}",
                ref_type,
                path.display()
            )));
        }

        // Return the prefix without .dmnd
        let prefix = path.with_extension("");

        Ok((prefix, vec![]))
    }

    impl ArgGenerator for DiamondArgGenerator {
        fn generate_args(
            &self,
            run_config: &RunConfig,
            extra: Option<&dyn std::any::Any>,
        ) -> AnyhowResult<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<DiamondConfig>())
                .ok_or_else(|| anyhow!("Diamond requires DiamondConfig in extra"))?;

            let mut args_vec = vec![];

            match config.subcommand {
                DiamondSubcommand::Blastx => args_vec.push("blastx".to_string()),
            }

            args_vec.push("-d".to_string());
            args_vec.push(config.db.to_string_lossy().to_string());


            let threads = run_config.thread_allocation(DIAMOND_TAG, None);
            args_vec.push("--threads".to_string());
            args_vec.push(threads.to_string());

            for (key, value) in &config.subcommand_fields {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            Ok(args_vec)
        }
    }

    pub async fn compute_optimal_block_size(run_config: &RunConfig) -> AnyhowResult<f64> {
        let (_, available_ram) = detect_ram()?;
        let available_ram_gb = available_ram as f64 / 1_073_741_824.0;

        let db_path = run_config.args.diamond_db.as_deref()
            .ok_or_else(|| anyhow!("--diamond-db required"))?;
        let db_path = if db_path.ends_with(".dmnd") {
            db_path.to_string()
        } else {
            format!("{}.dmnd", db_path)
        };

        let db_stats = get_diamond_db_stats(&db_path).await
            .unwrap_or((0, 300_000_000_000)); // fallback NR size

        let total_letters_billions = db_stats.1 as f64 / 1e9;

        let ram_factor = if available_ram_gb >= 1400.0 {
            9.0
        } else if available_ram_gb >= 900.0 {
            10.0
        } else if available_ram_gb >= 128.0 {
            12.0
        } else {
            15.0
        };

        let mut block_size = available_ram_gb / ram_factor;

        block_size = block_size.min(total_letters_billions);

        block_size = block_size.max(6.0);
        block_size = block_size.min(200.0);

        let scratch_path = PathBuf::from(run_config.args.nvme_scratch.as_deref().unwrap_or("."));
        let scratch_avail = available_space_for_path(&scratch_path).await.unwrap_or(0);
        let db_size_gb = std::fs::metadata(&db_path)
            .map(|m| m.len() as f64 / 1_073_741_824.0)
            .unwrap_or(50.0);

        let scratch_avail_gib = scratch_avail as f64 / 1_073_741_824.0;

        let estimated_scratch_gib = block_size * 3.0 + 20.0;

        if scratch_avail_gib < estimated_scratch_gib * 1.2 {
            warn!(
            "Low scratch space — need {:.1} GiB, have {:.1} GiB; reducing block size by 50%",
            estimated_scratch_gib, scratch_avail_gib
        );
            block_size *= 0.5;
        } else if scratch_avail_gib < estimated_scratch_gib * 1.5 {
            warn!(
            "Low scratch space — need {:.1} GiB, have {:.1} GiB; reducing block size by 30%",
            estimated_scratch_gib, scratch_avail_gib
        );
            block_size *= 0.7;
        }

        info!(
        "Diamond block size: {:.1} (RAM {:.0} GB, DB ~{:.0}B letters, scratch {:.1} GB free)",
        block_size, available_ram_gb, total_letters_billions, scratch_avail as f64 / 1e9 / 1_073_741_824.0
    );

        Ok(block_size)
    }

    async fn get_diamond_db_stats(db_path: &str) -> AnyhowResult<(u64, u64)> {
        let output = Command::new("diamond")
            .args(["dbinfo", "--db", db_path])
            .output()
            .await
            .map_err(|e| anyhow!("Failed to spawn diamond dbinfo: {}", e))?;

        if !output.status.success() {
            return Err(anyhow!(
            "diamond dbinfo failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        debug!("Diamond dbinfo output: {}", stdout);

        let seq_re = Regex::new(r"Sequences\s+(\d+)")?;
        let letters_re = Regex::new(r"Letters\s+(\d+)")?;

        let sequences = seq_re
            .captures(&stdout)
            .and_then(|c| c.get(1))
            .and_then(|m| m.as_str().parse::<u64>().ok())
            .ok_or_else(|| anyhow!("Failed to parse sequences from: {}", stdout))?;

        let letters = letters_re
            .captures(&stdout)
            .and_then(|c| c.get(1))
            .and_then(|m| m.as_str().parse::<u64>().ok())
            .ok_or_else(|| anyhow!("Failed to parse letters from: {}", stdout))?;

        Ok((sequences, letters))
    }


}


pub mod spades {
    // command.rs spades mod - no changes needed (presence check was not failing; kept as is with -v)
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{SPADES_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;
    use log::{debug, info, error, warn};

    pub struct SpadesConfig {
        pub r1_path: PathBuf,                // Always required (R1 or single-end)
        pub r2_path_opt: Option<PathBuf>,    // Some(R2) for paired, None for single-end
        pub outdir_path: PathBuf,
        pub option_fields: HashMap<String, Option<String>>,  // e.g., {"--only-assembler": None}
    }

    pub struct SpadesArgGenerator;

    pub async fn spades_presence_check(config: &RunConfig, version_file: Option<PathBuf>) -> Result<f32> {
        let version = version_check(SPADES_TAG, vec!["-v"], 0, 3, ChildStream::Stdout, version_file, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for SpadesArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<SpadesConfig>())
                .ok_or_else(|| anyhow!("SPAdes requires a SpadesConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            // Input flags: Separate R1/R2 for paired, single for unpaired
            if let Some(r2_path) = &config.r2_path_opt {
                args_vec.push("--pe1-1".to_string());
                args_vec.push(config.r1_path.to_string_lossy().to_string());
                args_vec.push("--pe1-2".to_string());
                args_vec.push(r2_path.to_string_lossy().to_string());
            } else {
                args_vec.push("-s".to_string());
                args_vec.push(config.r1_path.to_string_lossy().to_string());
            }

            // SPAdes scales decently up to ~80 threads on large machines
            // Beyond that, memory contention + scheduling overhead usually hurts
            //no minimum, to allow smaller machines to run
            let num_cores: usize = RunConfig::thread_allocation(run_config, SPADES_TAG, None).min(80);
            args_vec.push("-t".to_string());
            args_vec.push(num_cores.to_string());

            // Memory
            // Scale down aggressively on smaller machines
            let available_gb = (run_config.available_ram as f64 / 1_000_000_000.0) as u64;
            let spades_gb = match available_gb {
                0..=128 => available_gb / 2,
                129..=512 => available_gb / 2,
                _ => 1000,
            }.max(8);  // minimum 8 GB to avoid SPAdes complaining

            info!(
            "SPAdes memory allocation: {} GB (available: {} GB, capped for safety)",
            spades_gb, available_gb
        );
            args_vec.push("-m".to_string());
            args_vec.push(spades_gb.to_string());

            // Custom flags from config (e.g., --only-assembler)
            for (key, value) in config.option_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            // Output dir
            args_vec.push("-o".to_string());
            args_vec.push(config.outdir_path.to_string_lossy().to_string());

            Ok(args_vec)
        }
    }
}

pub mod makeblastdb {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{MAKEBLASTDB_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    #[derive(Debug)]
    pub struct MakeblastdbConfig {
        pub input: PathBuf,  // reference_fasta
        pub dbtype: String,  // "nucl" or "prot"
        pub output: PathBuf, // blast_index_path (without extension)
        pub option_fields: HashMap<String, Option<String>>,
    }

    pub struct MakeblastdbArgGenerator;

    pub async fn makeblastdb_presence_check(config: &RunConfig)-> Result<f32> {
        let version = version_check(MAKEBLASTDB_TAG, vec!["-version"], 0, 1, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for MakeblastdbArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<MakeblastdbConfig>())
                .ok_or_else(|| anyhow!("makeblastdb requires MakeblastdbConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            args_vec.push("-in".to_string());
            args_vec.push(config.input.to_string_lossy().to_string());

            args_vec.push("-dbtype".to_string());
            args_vec.push(config.dbtype.clone());

            args_vec.push("-out".to_string());
            args_vec.push(config.output.to_string_lossy().to_string());

            for (key, value) in config.option_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            Ok(args_vec)
        }
    }
}

pub mod blastn {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{BLASTN_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    #[derive(Debug)]
    pub struct BlastnConfig {
        pub query: PathBuf,          // assembled_contig
        pub db: PathBuf,             // blast_index_path (without extension)
        pub outfmt: String,          // e.g., "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        pub evalue: f64,             // e.g. 1e-10
        pub max_target_seqs: u32,    // e.g. 5000
        pub option_fields: HashMap<String, Option<String>>,
    }

    pub struct BlastnArgGenerator;

    pub async fn blastn_presence_check(config: &RunConfig)-> Result<f32> {
        let version = version_check(BLASTN_TAG, vec!["-version"], 0, 1, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for BlastnArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<BlastnConfig>())
                .ok_or_else(|| anyhow!("blastn requires BlastnConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            args_vec.push("-query".to_string());
            args_vec.push(config.query.to_string_lossy().to_string());

            args_vec.push("-db".to_string());
            args_vec.push(config.db.to_string_lossy().to_string());


            args_vec.push("-outfmt".to_string());
            args_vec.push(config.outfmt.clone());

            args_vec.push("-evalue".to_string());
            args_vec.push(config.evalue.to_string());

            args_vec.push("-max_target_seqs".to_string());
            args_vec.push(config.max_target_seqs.to_string());

            let num_cores: usize = RunConfig::thread_allocation(run_config, BLASTN_TAG, None);
            args_vec.push("-num_threads".to_string());
            args_vec.push(num_cores.to_string());

            for (key, value) in config.option_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            Ok(args_vec)
        }
    }


}

pub mod blastx {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{BLASTX_TAG, RunConfig};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    #[derive(Debug)]
    pub struct BlastxConfig {
        pub query: PathBuf,          // assembled_contig
        pub db: PathBuf,             // blast_index_path (without extension)
        pub outfmt: u32,             // 6 (standard tabular)
        pub evalue: f64,             // e.g. 1e-10
        pub num_alignments: u32,     // e.g. 5
        pub option_fields: HashMap<String, Option<String>>,
    }

    pub struct BlastxArgGenerator;

    pub async fn blastx_presence_check(config: &RunConfig)-> Result<f32> {
        let version = version_check(BLASTX_TAG, vec!["-version"], 0, 1, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for BlastxArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<BlastxConfig>())
                .ok_or_else(|| anyhow!("blastx requires BlastxConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            args_vec.push("-query".to_string());
            args_vec.push(config.query.to_string_lossy().to_string());

            args_vec.push("-db".to_string());
            args_vec.push(config.db.to_string_lossy().to_string());


            args_vec.push("-outfmt".to_string());
            args_vec.push(config.outfmt.to_string());

            args_vec.push("-evalue".to_string());
            args_vec.push(config.evalue.to_string());

            args_vec.push("-num_alignments".to_string());
            args_vec.push(config.num_alignments.to_string());

            let num_cores: usize = RunConfig::thread_allocation(run_config, BLASTX_TAG, None);
            args_vec.push("-num_threads".to_string());
            args_vec.push(num_cores.to_string());

            for (key, value) in config.option_fields.iter() {
                args_vec.push(key.clone());
                if let Some(v) = value {
                    args_vec.push(v.clone());
                }
            }

            Ok(args_vec)
        }
    }
}

// ────────────────────────────────────────────────────────────────
// GNU sort
// ────────────────────────────────────────────────────────────────
pub mod sort {
    use std::collections::HashMap;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{RunConfig, PipelineError};
    use crate::utils::command::{ArgGenerator, version_check};
    use crate::utils::streams::ChildStream;

    #[derive(Debug)]
    pub struct SortConfig {
        pub key: String,                    // e.g. "-k1,1"
        pub parallel: Option<usize>,
        pub buffer_size: Option<String>,    // e.g. "50%"
        pub temp_dir: Option<String>,       // -T /path
        pub output: Option<String>,         // -o sorted.m8
        pub extra_fields: HashMap<String, Option<String>>,
        pub input: Option<String>,
    }

    pub struct SortArgGenerator;

    pub async fn sort_presence_check(config: &RunConfig)-> Result<f32> {
        // GNU sort doesn't have a clean --version that always works the same way,
        // but "sort --version" is reliable on all modern coreutils.
        let version = version_check("sort", vec!["--version"], 0, 1, ChildStream::Stdout, None, &config).await?;
        Ok(version)
    }

    impl ArgGenerator for SortArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<SortConfig>())
                .ok_or_else(|| anyhow!("sort requires a SortConfig as extra argument"))?;

            let mut args: Vec<String> = vec![];

            // Core sorting key (required for m8: sort by read ID = column 1)
            args.push(config.key.clone());

            // Parallelism — default to full allocation unless overridden
            let threads = config.parallel
                .unwrap_or_else(|| run_config.thread_allocation("sort", None));
            if threads > 1 {
                args.push("--parallel".to_string());
                args.push(threads.to_string());
            }

            // Memory buffer
            if let Some(buf) = &config.buffer_size {
                args.push("-S".to_string());
                args.push(buf.clone());
            } else {
                args.push("-S".to_string());
                args.push("50%".to_string());
            }

            // Temp directory (critical for large sorts on NVMe)
            if let Some(td) = &config.temp_dir {
                args.push("-T".to_string());
                args.push(td.clone());
            }

            // Output file (-o)
            if let Some(out) = &config.output {
                args.push("-o".to_string());
                args.push(out.clone());
            }

            // Any extra flags the caller wants
            for (key, value) in &config.extra_fields {
                args.push(key.clone());
                if let Some(v) = value {
                    args.push(v.clone());
                }
            }

            // Input
            if let Some(input_file) = &config.input{
                args.push(input_file.clone());
            }

            Ok(args)
        }
    }
}



pub mod mmseqs {
    use std::any::Any;
    use std::collections::HashMap;
    use std::path::{Path, PathBuf};

    use anyhow::{anyhow, Result};
    use log::{debug, info, warn};
    use tokio::process::{Child, Command};

    use crate::config::defs::{MMSEQS_TAG, RunConfig};
    use crate::utils::command::ArgGenerator;

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum MmseqsBackend {
        Cpu,
        Gpu,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum MmseqsSubcommand {
        Version,
        Createdb,
        MakePaddedSeqDb,
        CreateIndex,
        EasySearch,
        Search,
        ConvertAlis,
        GpuServer,
    }

    #[derive(Debug, Clone)]
    pub struct MmseqsConfig {
        pub subcommand: MmseqsSubcommand,
        pub backend: MmseqsBackend,

        pub input: Option<PathBuf>,
        pub target_db: Option<PathBuf>,
        pub result_db: Option<PathBuf>,
        pub output: Option<PathBuf>,
        pub tmp_dir: Option<PathBuf>,

        pub threads: Option<usize>,
        pub sensitivity: Option<String>,

        pub search_type: Option<String>,
        pub max_seqs: Option<String>,
        pub prefilter_mode: Option<String>,
        pub db_load_mode: Option<String>,
        pub alignment_mode: Option<String>,
        pub index_subset: Option<String>,

        pub format_output: Option<String>,
        pub cuda_visible_devices: Option<String>,

        /// Extra raw args appended verbatim. Use this for rare flags that do not
        /// deserve a dedicated field yet.
        pub option_fields: HashMap<String, Option<String>>,

        /// If true, generate client-side flags for MMseqs GPU server mode.
        /// This is only meaningful for Search / ConvertAlis workflows, not GpuServer itself.
        pub gpu_server: bool,
    }

    impl Default for MmseqsConfig {
        fn default() -> Self {
            Self {
                subcommand: MmseqsSubcommand::EasySearch,
                backend: MmseqsBackend::Cpu,
                input: None,
                target_db: None,
                result_db: None,
                output: None,
                tmp_dir: None,
                threads: None,
                sensitivity: None,
                search_type: None,
                max_seqs: None,
                prefilter_mode: None,
                db_load_mode: None,
                alignment_mode: None,
                index_subset: None,
                format_output: None,
                cuda_visible_devices: None,
                option_fields: HashMap::new(),
                gpu_server: false,
            }
        }
    }

    pub struct MmseqsArgGenerator;

    fn required_path<'a>(value: &'a Option<PathBuf>, what: &str) -> Result<&'a PathBuf> {
        value
            .as_ref()
            .ok_or_else(|| anyhow!("MMseqs requires {}", what))
    }

    fn resolve_mmseqs_db(run_config: &RunConfig, explicit: Option<&PathBuf>) -> Result<PathBuf> {
        if let Some(p) = explicit {
            return Ok(p.clone());
        }

        let db = run_config
            .args
            .mmseqs_db
            .as_ref()
            .ok_or_else(|| anyhow!("MMseqs requires --mmseqs-db"))?;

        let path = PathBuf::from(db);
        if !path.exists() {
            return Err(anyhow!("MMseqs DB path does not exist: {}", path.display()));
        }

        Ok(path)
    }

    fn validate_cpu_index(db_prefix: &Path) -> Result<()> {
        let candidates = [
            db_prefix.with_extension("index"),
            db_prefix.with_extension("idx"),
            db_prefix.with_extension("lookup"),
            db_prefix.with_extension("source"),
        ];

        if candidates.iter().any(|p| p.exists()) {
            return Ok(());
        }

        warn!(
            "MMseqs CPU backend: no obvious index sidecar found near {}. \
             Proceeding anyway, but repeated searches will be slower.",
            db_prefix.display()
        );

        Ok(())
    }

    fn push_threads(args: &mut Vec<String>, run_config: &RunConfig, threads: Option<usize>) {
        let n = threads.unwrap_or_else(|| run_config.thread_allocation(MMSEQS_TAG, None));
        args.push("--threads".to_string());
        args.push(n.to_string());
    }

    fn push_optional_arg(args: &mut Vec<String>, flag: &str, value: &Option<String>) {
        if let Some(v) = value {
            args.push(flag.to_string());
            args.push(v.clone());
        }
    }

    fn push_option_fields(args: &mut Vec<String>, fields: &HashMap<String, Option<String>>) {
        for (key, value) in fields {
            args.push(key.clone());
            if let Some(v) = value {
                args.push(v.clone());
            }
        }
    }

    fn push_gpu_client_flags(args: &mut Vec<String>, config: &MmseqsConfig) {
        if config.backend != MmseqsBackend::Gpu {
            return;
        }

        args.push("--gpu".to_string());
        args.push("1".to_string());
        debug!("MMseqs2: GPU client mode enabled (--gpu 1)");

        if config.gpu_server {
            args.push("--gpu-server".to_string());
            args.push("1".to_string());

            args.push("--db-load-mode".to_string());
            args.push(config.db_load_mode.clone().unwrap_or_else(|| "2".to_string()));

            if config.prefilter_mode.is_none() {
                args.push("--prefilter-mode".to_string());
                args.push("1".to_string());
            }

            debug!("MMseqs2: GPU server client mode enabled (--gpu-server 1 --db-load-mode 2)");
        }
    }

    fn push_search_common_flags(
        args: &mut Vec<String>,
        run_config: &RunConfig,
        config: &MmseqsConfig,
        is_easy_search: bool,
    ) {
        push_threads(args, run_config, config.threads);

        if let Some(st) = &config.search_type {
            args.push("--search-type".to_string());
            args.push(st.clone());
        }

        if let Some(ms) = &config.max_seqs {
            args.push("--max-seqs".to_string());
            args.push(ms.clone());
        }

        if let Some(pm) = &config.prefilter_mode {
            args.push("--prefilter-mode".to_string());
            args.push(pm.clone());
        }

        if config.backend == MmseqsBackend::Cpu {
            if let Some(s) = &config.sensitivity {
                args.push("-s".to_string());
                args.push(s.clone());
            }
        } else if !is_easy_search {
            if let Some(s) = &config.sensitivity {
                args.push("-s".to_string());
                args.push(s.clone());
            }
        }

        push_option_fields(args, &config.option_fields);
    }

    pub async fn mmseqs_presence_check(config: &RunConfig)-> Result<f32> {
        let output = Command::new(MMSEQS_TAG)
            .arg("version")
            .output()
            .await
            .map_err(|e| anyhow!("Failed to spawn mmseqs. Is it installed? Error: {}", e))?;

        if !output.status.success() {
            return Err(anyhow!(
                "mmseqs version check failed with exit code {:?}",
                output.status.code()
            ));
        }

        Ok(0.0)
    }

    pub async fn spawn_gpuserver(
        target_db: &Path,
        cuda_visible_devices: Option<&str>,
    ) -> Result<Child> {
        if !target_db.exists() {
            return Err(anyhow!(
                "MMseqs GPU server target DB does not exist: {}",
                target_db.display()
            ));
        }

        let mut cmd = Command::new(MMSEQS_TAG);
        cmd.arg("gpuserver");
        cmd.arg(target_db);

        if let Some(devices) = cuda_visible_devices {
            cmd.env("CUDA_VISIBLE_DEVICES", devices);
        }

        // Important: current MMseqs2 release notes say gpuserver no longer accepts --gpu.
        // The client gets the GPU flags; the server is started without --gpu.
        cmd.stdout(std::process::Stdio::null());
        cmd.stderr(std::process::Stdio::inherit());

        let child = cmd
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn mmseqs gpuserver: {}", e))?;

        Ok(child)
    }

    pub async fn stop_gpuserver(child: &mut Child) -> Result<()> {
        match child.try_wait() {
            Ok(Some(_)) => return Ok(()),
            Ok(None) => {}
            Err(e) => return Err(anyhow!("Failed to query gpuserver state: {}", e)),
        }

        child
            .kill()
            .await
            .map_err(|e| anyhow!("Failed to stop mmseqs gpuserver: {}", e))?;

        let _ = child.wait().await;
        Ok(())
    }

    pub fn spawn_with_numa(
        config: &RunConfig,
        program: &str,
        mut args: Vec<String>,
        label: &str,
    ) -> Result<Command> {
        let mut cmd = if cfg!(target_os = "linux") && config.max_cores >= 64 {
            // Dual-socket EPYC (r6id, r8id, 256-core nodes) → use interleave
            let mut c = Command::new("numactl");
            c.arg("--interleave=all");           // Best for large memory workloads
            c.arg(program);
            c.args(&args);
            debug!("{}: Launched with numactl --interleave=all", label);
            c
        } else {
            // MacBook, small machines, or non-Linux → plain command
            let mut c = Command::new(program);
            c.args(&args);
            c
        };

        Ok(cmd)
    }

    impl ArgGenerator for MmseqsArgGenerator {
        fn generate_args(
            &self,
            run_config: &RunConfig,
            extra: Option<&dyn Any>,
        ) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<MmseqsConfig>())
                .ok_or_else(|| anyhow!("MMseqs requires a MmseqsConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand {
                MmseqsSubcommand::Version => {
                    args_vec.push("version".to_string());
                }

                MmseqsSubcommand::Createdb => {
                    let input = required_path(&config.input, "input FASTA/FASTQ for createdb")?;
                    let output = required_path(&config.output, "output DB for createdb")?;

                    args_vec.push("createdb".to_string());
                    args_vec.push(input.to_string_lossy().to_string());
                    args_vec.push(output.to_string_lossy().to_string());

                    // If the caller wants to create a GPU DB directly, pass --gpu via option_fields.
                    push_option_fields(&mut args_vec, &config.option_fields);
                }

                MmseqsSubcommand::MakePaddedSeqDb => {
                    let input = required_path(&config.input, "input DB for makepaddedseqdb")?;
                    let output = required_path(&config.output, "output GPU DB for makepaddedseqdb")?;

                    args_vec.push("makepaddedseqdb".to_string());
                    args_vec.push(input.to_string_lossy().to_string());
                    args_vec.push(output.to_string_lossy().to_string());

                    push_option_fields(&mut args_vec, &config.option_fields);
                }

                MmseqsSubcommand::CreateIndex => {
                    let target_db = resolve_mmseqs_db(run_config, config.target_db.as_ref())?;
                    let tmp_dir = required_path(&config.tmp_dir, "tmp dir for createindex")?;

                    args_vec.push("createindex".to_string());
                    args_vec.push(target_db.to_string_lossy().to_string());
                    args_vec.push(tmp_dir.to_string_lossy().to_string());

                    if let Some(subset) = &config.index_subset {
                        args_vec.push("--index-subset".to_string());
                        args_vec.push(subset.clone());
                    }

                    push_option_fields(&mut args_vec, &config.option_fields);
                }

                MmseqsSubcommand::GpuServer => {
                    let target_db = resolve_mmseqs_db(run_config, config.target_db.as_ref())?;

                    args_vec.push("gpuserver".to_string());
                    args_vec.push(target_db.to_string_lossy().to_string());

                    // Do NOT add --gpu here. Current MMseqs2 release notes say gpuserver
                    // no longer accepts --gpu; the client side gets the GPU flags.
                    if let Some(devs) = &config.cuda_visible_devices {
                        args_vec.push("--cuda-visible-devices".to_string());
                        args_vec.push(devs.clone());
                    }

                    push_option_fields(&mut args_vec, &config.option_fields);
                }

                MmseqsSubcommand::EasySearch => {
                    let query = required_path(&config.input, "query FASTA/FASTQ for easy-search")?;
                    let target_db = resolve_mmseqs_db(run_config, config.target_db.as_ref())?;
                    let output = required_path(&config.output, "output m8 for easy-search")?;
                    let tmp_dir = required_path(&config.tmp_dir, "tmp dir for easy-search")?;

                    args_vec.push("easy-search".to_string());
                    push_gpu_client_flags(&mut args_vec, config);
                    push_threads(&mut args_vec, run_config, config.threads);

                    if config.backend == MmseqsBackend::Cpu {
                        push_optional_arg(&mut args_vec, "-s", &config.sensitivity);
                    }

                    if let Some(st) = &config.search_type {
                        args_vec.push("--search-type".to_string());
                        args_vec.push(st.clone());
                    }

                    if let Some(ms) = &config.max_seqs {
                        args_vec.push("--max-seqs".to_string());
                        args_vec.push(ms.clone());
                    }

                    if let Some(pm) = &config.prefilter_mode {
                        args_vec.push("--prefilter-mode".to_string());
                        args_vec.push(pm.clone());
                    }

                    args_vec.push(query.to_string_lossy().to_string());
                    args_vec.push(target_db.to_string_lossy().to_string());
                    args_vec.push(output.to_string_lossy().to_string());
                    args_vec.push(tmp_dir.to_string_lossy().to_string());

                    args_vec.push("--format-output".to_string());
                    args_vec.push(config.format_output.clone().unwrap_or_else(|| {
                        "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
                            .to_string()
                    }));

                    push_option_fields(&mut args_vec, &config.option_fields);
                }

                MmseqsSubcommand::Search => {
                    let query = required_path(&config.input, "query DB/FASTA for search")?;
                    let target_db = resolve_mmseqs_db(run_config, config.target_db.as_ref())?;
                    let result_db = required_path(&config.result_db, "result DB for search")?;
                    let tmp_dir = required_path(&config.tmp_dir, "tmp dir for search")?;

                    args_vec.push("search".to_string());
                    push_gpu_client_flags(&mut args_vec, config);
                    push_threads(&mut args_vec, run_config, config.threads);

                    if let Some(s) = &config.sensitivity {
                        args_vec.push("-s".to_string());
                        args_vec.push(s.clone());
                    }

                    // Production-safe default: search does not match easy-search identity semantics
                    // unless alignment-mode is set. We default to 3 so pident is real rather than
                    // estimated, unless the caller overrides it.
                    args_vec.push("--alignment-mode".to_string());
                    args_vec.push(
                        config
                            .alignment_mode
                            .clone()
                            .unwrap_or_else(|| "3".to_string()),
                    );

                    if let Some(st) = &config.search_type {
                        args_vec.push("--search-type".to_string());
                        args_vec.push(st.clone());
                    }

                    if let Some(ms) = &config.max_seqs {
                        args_vec.push("--max-seqs".to_string());
                        args_vec.push(ms.clone());
                    }

                    if let Some(pm) = &config.prefilter_mode {
                        args_vec.push("--prefilter-mode".to_string());
                        args_vec.push(pm.clone());
                    }

                    if config.gpu_server && config.backend == MmseqsBackend::Gpu {
                        if config.prefilter_mode.is_none() {
                            args_vec.push("--prefilter-mode".to_string());
                            args_vec.push("1".to_string());
                        }
                    }

                    args_vec.push(query.to_string_lossy().to_string());
                    args_vec.push(target_db.to_string_lossy().to_string());
                    args_vec.push(result_db.to_string_lossy().to_string());
                    args_vec.push(tmp_dir.to_string_lossy().to_string());

                    push_option_fields(&mut args_vec, &config.option_fields);
                }

                MmseqsSubcommand::ConvertAlis => {
                    let query_db = required_path(&config.input, "query DB for convertalis")?;
                    let target_db = required_path(&config.target_db, "target DB for convertalis")?;
                    let result_db = required_path(&config.result_db, "result DB for convertalis")?;
                    let output = required_path(&config.output, "output m8 for convertalis")?;

                    args_vec.push("convertalis".to_string());
                    args_vec.push(query_db.to_string_lossy().to_string());
                    args_vec.push(target_db.to_string_lossy().to_string());
                    args_vec.push(result_db.to_string_lossy().to_string());
                    args_vec.push(output.to_string_lossy().to_string());

                    // Keep the same text columns you are already using unless overridden.
                    args_vec.push("--format-output".to_string());
                    args_vec.push(config.format_output.clone().unwrap_or_else(|| {
                        "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
                            .to_string()
                    }));

                    push_option_fields(&mut args_vec, &config.option_fields);
                }
            }

            if matches!(config.subcommand, MmseqsSubcommand::EasySearch | MmseqsSubcommand::Search)
                && config.backend == MmseqsBackend::Cpu
            {
                let target = resolve_mmseqs_db(run_config, config.target_db.as_ref())?;
                validate_cpu_index(&target)?;
            }

            if config.backend == MmseqsBackend::Gpu
                && matches!(config.subcommand, MmseqsSubcommand::EasySearch | MmseqsSubcommand::Search)
            {
                debug!(
                    "MMseqs GPU client selected; GPU-compatible DB should be padded/indexed already"
                );
            }

            debug!("MMseqs argv: {:?}", args_vec);
            Ok(args_vec)
        }
    }

    /// Convenience helper for the new GPU server workflow.
    ///
    /// Typical usage:
    /// 1) createdb / makepaddedseqdb / createindex --index-subset 2
    /// 2) spawn_gpuserver(...)
    /// 3) run Search with backend=Gpu and gpu_server=true
    /// 4) run ConvertAlis if you want m8 output
    pub async fn build_gpu_server_client_args(
        run_config: &RunConfig,
        config: &MmseqsConfig,
    ) -> Result<Vec<String>> {
        let generator = MmseqsArgGenerator;
        generator.generate_args(run_config, Some(config as &dyn Any))
    }

    pub fn mmseqs_gpu_server_note() -> &'static str {
        "Use gpuserver without --gpu, then use search with --gpu 1 --gpu-server 1 --db-load-mode 2."
    }

    pub fn mmseqs_cpu_production_note() -> &'static str {
        "Use search for production control; add --alignment-mode 3 to preserve real pident, or keep easy-search for compatibility."
    }

    pub fn mmseqs_gpu_db_prep_note() -> &'static str {
        "GPU DB prep: createdb -> makepaddedseqdb -> createindex --index-subset 2."
    }

    pub fn default_easy_search_format_output() -> &'static str {
        "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    }

    pub fn default_search_alignment_mode() -> &'static str {
        "3"
    }

    pub fn default_gpu_server_db_load_mode() -> &'static str {
        "2"
    }

    pub fn default_gpu_server_prefilter_mode() -> &'static str {
        "1"
    }

    pub fn default_gpu_createindex_subset() -> &'static str {
        "2"
    }

    pub fn recommended_cuda_visible_devices_for_server(all_visible: bool) -> Option<String> {
        if all_visible {
            None
        } else {
            Some("0".to_string())
        }
    }

    pub fn example_gpu_server_workflow(
        query_db: PathBuf,
        target_db: PathBuf,
        result_db: PathBuf,
        out_m8: PathBuf,
        tmp_dir: PathBuf,
    ) -> (MmseqsConfig, MmseqsConfig, MmseqsConfig) {
        let search_cfg = MmseqsConfig {
            subcommand: MmseqsSubcommand::Search,
            backend: MmseqsBackend::Gpu,
            input: Some(query_db),
            target_db: Some(target_db.clone()),
            result_db: Some(result_db),
            output: None,
            tmp_dir: Some(tmp_dir),
            threads: None,
            sensitivity: None,
            search_type: None,
            max_seqs: None,
            prefilter_mode: None,
            db_load_mode: Some(default_gpu_server_db_load_mode().to_string()),
            alignment_mode: Some(default_search_alignment_mode().to_string()),
            index_subset: None,
            format_output: None,
            cuda_visible_devices: None,
            option_fields: HashMap::new(),
            gpu_server: true,
        };

        let convert_cfg = MmseqsConfig {
            subcommand: MmseqsSubcommand::ConvertAlis,
            backend: MmseqsBackend::Cpu,
            input: None,
            target_db: Some(target_db),
            result_db: search_cfg.result_db.clone(),
            output: Some(out_m8),
            tmp_dir: None,
            threads: None,
            sensitivity: None,
            search_type: None,
            max_seqs: None,
            prefilter_mode: None,
            db_load_mode: None,
            alignment_mode: None,
            index_subset: None,
            format_output: Some(default_easy_search_format_output().to_string()),
            cuda_visible_devices: None,
            option_fields: HashMap::new(),
            gpu_server: false,
        };

        let server_cfg = MmseqsConfig {
            subcommand: MmseqsSubcommand::GpuServer,
            backend: MmseqsBackend::Gpu,
            input: None,
            target_db: search_cfg.target_db.clone(),
            result_db: None,
            output: None,
            tmp_dir: None,
            threads: None,
            sensitivity: None,
            search_type: None,
            max_seqs: None,
            prefilter_mode: None,
            db_load_mode: None,
            alignment_mode: None,
            index_subset: None,
            format_output: None,
            cuda_visible_devices: None,
            option_fields: HashMap::new(),
            gpu_server: false,
        };

        (server_cfg, search_cfg, convert_cfg)
    }

    pub async fn mmseqs_createdb(
        input: &Path,
        output: &Path,
    ) -> Result<Command> {
        let mut cmd = Command::new(MMSEQS_TAG);
        cmd.arg("createdb");
        cmd.arg(input);
        cmd.arg(output);
        Ok(cmd)
    }

    pub async fn mmseqs_makepaddedseqdb(
        input_db: &Path,
        output_db: &Path,
    ) -> Result<Command> {
        let mut cmd = Command::new(MMSEQS_TAG);
        cmd.arg("makepaddedseqdb");
        cmd.arg(input_db);
        cmd.arg(output_db);
        Ok(cmd)
    }

    pub async fn mmseqs_createindex(
        target_db: &Path,
        tmp_dir: &Path,
        index_subset: Option<&str>,
    ) -> Result<Command> {
        let mut cmd = Command::new(MMSEQS_TAG);
        cmd.arg("createindex");
        cmd.arg(target_db);
        cmd.arg(tmp_dir);

        if let Some(subset) = index_subset {
            cmd.arg("--index-subset");
            cmd.arg(subset);
        }

        Ok(cmd)
    }

    pub async fn mmseqs_convertalis(
        query_db: &Path,
        target_db: &Path,
        result_db: &Path,
        output_m8: &Path,
        format_output: Option<&str>,
    ) -> Result<Command> {
        let mut cmd = Command::new(MMSEQS_TAG);
        cmd.arg("convertalis");
        cmd.arg(query_db);
        cmd.arg(target_db);
        cmd.arg(result_db);
        cmd.arg(output_m8);

        if let Some(fmt) = format_output {
            cmd.arg("--format-output");
            cmd.arg(fmt);
        }

        Ok(cmd)
    }

    pub fn is_gpu_server_client(config: &MmseqsConfig) -> bool {
        config.backend == MmseqsBackend::Gpu && config.gpu_server
    }

    pub fn uses_search_command(config: &MmseqsConfig) -> bool {
        matches!(
            config.subcommand,
            MmseqsSubcommand::EasySearch | MmseqsSubcommand::Search
        )
    }

    pub fn uses_convertalis_command(config: &MmseqsConfig) -> bool {
        matches!(config.subcommand, MmseqsSubcommand::ConvertAlis)
    }

    pub fn uses_server_command(config: &MmseqsConfig) -> bool {
        matches!(config.subcommand, MmseqsSubcommand::GpuServer)
    }

    pub fn uses_db_prep_command(config: &MmseqsConfig) -> bool {
        matches!(
            config.subcommand,
            MmseqsSubcommand::Createdb
                | MmseqsSubcommand::MakePaddedSeqDb
                | MmseqsSubcommand::CreateIndex
        )
    }

    pub fn maybe_warn_on_easy_search_cpu(config: &MmseqsConfig) {
        if matches!(config.subcommand, MmseqsSubcommand::EasySearch)
            && config.backend == MmseqsBackend::Cpu
        {
            info!(
                "MMseqs CPU easy-search selected. This preserves current behavior."
            );
        }
    }
}





pub fn generate_cli(tool: &str, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> Result<Vec<String>> {
    let generator: Box<dyn ArgGenerator> = match tool {
        FASTP_TAG => Box::new(fastp::FastpArgGenerator),
        PIGZ_TAG => Box::new(pigz::PigzArgGenerator),
        MINIMAP2_TAG => Box::new(minimap2::Minimap2ArgGenerator),
        SAMTOOLS_TAG => Box::new(samtools::SamtoolsArgGenerator),
        KRAKEN2_TAG => Box::new(kraken2::Kraken2ArgGenerator),
        BCFTOOLS_TAG => Box::new(bcftools::BcftoolsArgGenerator),
        MAFFT_TAG => Box::new(mafft::MafftArgGenerator),
        NUCMER_TAG => Box::new(nucmer::NucmerArgGenerator),
        SHOW_COORDS_TAG => Box::new(show_coords::ShowCoordsArgGenerator),
        QUAST_TAG => Box::new(quast::QuastArgGenerator),
        SEQKIT_TAG => Box::new(seqkit::SeqkitArgGenerator),
        H5DUMP_TAG => return Err(anyhow!("h5dump argument generation not implemented")),
        BOWTIE2_TAG => Box::new(bowtie2::Bowtie2ArgGenerator),
        HISAT2_TAG => Box::new(hisat2::Hisat2ArgGenerator),
        KALLISTO_TAG => Box::new(kallisto::KallistoArgGenerator),
        STAR_TAG => Box::new(star::StarArgGenerator),
        CZID_DEDUP_TAG => Box::new(czid_dedup::CzidDedupArgGenerator),
        DIAMOND_TAG => Box::new(diamond::DiamondArgGenerator),
        SPADES_TAG => Box::new(spades::SpadesArgGenerator),
        BLASTN_TAG => Box::new(blastn::BlastnArgGenerator),
        BLASTX_TAG => Box::new(blastx::BlastxArgGenerator),
        MAKEBLASTDB_TAG => {Box::new(makeblastdb::MakeblastdbArgGenerator)},
        SORT_TAG => Box::new(sort::SortArgGenerator),
        MMSEQS_TAG => Box::new(mmseqs::MmseqsArgGenerator),
        _ => return Err(anyhow!("Unknown tool: {}", tool)),
    };

    generator.generate_args(run_config, extra)
}

pub async fn check_versions(tools: Vec<&str>, out_dir: &PathBuf, config: &RunConfig) -> Result<()> {
    let assembly_dir = out_dir.join("assembly");
    fs::create_dir_all(&assembly_dir).await
        .map_err(|e| anyhow!("Failed to create assembly directory: {}", e))?;

    let out_dir_cloned = out_dir.clone();
    let assembly_dir_cloned = assembly_dir.clone();

    let checks = tools.into_iter().map(move |tool| {
        let out_dir = out_dir_cloned.clone();
        let assembly_dir = assembly_dir_cloned.clone();

        async move {
            let version_file = match tool {
                SPADES_TAG => Some(assembly_dir.join("assembly_version.txt")),
                MINIMAP2_TAG => Some(out_dir.join("minimap2_version.txt")),
                DIAMOND_TAG => Some(out_dir.join("diamond_version.txt")),
                _ => None,
            };

            let version = match tool {
                FASTP_TAG => fastp::fastp_presence_check(&config).await,
                H5DUMP_TAG => h5dump::h5dump_presence_check(&config).await,
                MINIMAP2_TAG => minimap2::minimap2_presence_check(&config, version_file).await,
                SAMTOOLS_TAG => samtools::samtools_presence_check(&config).await,
                KRAKEN2_TAG => kraken2::kraken2_presence_check(&config).await,
                BCFTOOLS_TAG => bcftools::bcftools_presence_check(&config).await,
                IVAR_TAG => ivar::ivar_presence_check(&config).await,
                MUSCLE_TAG => muscle::muscle_presence_check(&config).await,
                MAFFT_TAG => mafft::mafft_presence_check(&config).await,
                QUAST_TAG => quast::quast_presence_check(&config).await,
                NUCMER_TAG => nucmer::nucmer_presence_check(&config).await,
                SEQKIT_TAG => seqkit::seqkit_presence_check(&config).await,
                BOWTIE2_TAG => bowtie2::bowtie2_presence_check(&config).await,
                HISAT2_TAG => hisat2::hisat2_presence_check(&config).await,
                KALLISTO_TAG => kallisto::kallisto_presence_check(&config).await,
                STAR_TAG => star::star_presence_check(&config).await,
                CZID_DEDUP_TAG => czid_dedup::czid_dedup_presence_check(&config).await,
                DIAMOND_TAG => diamond::diamond_presence_check(&config, version_file).await,
                SPADES_TAG => spades::spades_presence_check(&config, version_file).await,
                BLASTN_TAG => blastn::blastn_presence_check(&config).await,
                BLASTX_TAG => blastx::blastx_presence_check(&config).await,
                MAKEBLASTDB_TAG => makeblastdb::makeblastdb_presence_check(&config).await,
                SORT_TAG => sort::sort_presence_check(&config).await,
                MMSEQS_TAG => mmseqs::mmseqs_presence_check(&config).await,
                _ => return Err(anyhow!("Unknown tool: {}", tool)),
            }?;
            Ok((tool.to_string(), version))
        }
    });

    let results = try_join_all(checks).await?;
    let mut failed_tools = Vec::new();

    for (tool, version) in &results {
        if let Some(&min_version) = TOOL_VERSIONS.get(tool.as_str()) {
            if *version < min_version {
                failed_tools.push(format!(
                    "{} (version found: {}, required: >= {})",
                    tool, version, min_version
                ));
            }
        } else {
            warn!("Warning: No minimum version specified for tool: {}", tool);
        }
    }

    if !failed_tools.is_empty() {
        return Err(anyhow!(
            "The following tools have versions below the minimum required: {}",
            failed_tools.join(", ")
        ));
    }

    Ok(())
}