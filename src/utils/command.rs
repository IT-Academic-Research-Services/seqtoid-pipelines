/// Functions and structs for working with creating command-line arguments

use anyhow::{anyhow, Result};
use num_cpus;
use tokio::process::Command;
use futures::future::try_join_all;
use crate::config::defs::{RunConfig, PipelineError, TOOL_VERSIONS, FASTP_TAG, PIGZ_TAG, H5DUMP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, IVAR_TAG, MUSCLE_TAG, MAFFT_TAG, QUAST_TAG, NUCMER_TAG, SHOW_COORDS_TAG, SEQKIT_TAG, BOWTIE2_TAG, HISAT2_TAG, KALLISTO_TAG, STAR_TAG, FASTA_EXTS};
use crate::cli::Arguments;
use crate::utils::streams::{read_child_output_to_vec, ChildStream};
use std::path::PathBuf;
use std::sync::Arc;
use tokio::fs;
use tokio::task::JoinHandle;
use tempfile::NamedTempFile;
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
pub async fn version_check(command_tag: &str, version_args: Vec<&str>, version_line: usize, version_column: usize, child_stream: ChildStream) -> Result<f32> {
    let cmd_tag_owned = command_tag.to_string();
    // eprintln!("Running command: {}", &cmd_tag_owned);
    let args: Vec<&str> = version_args;

    let mut child = Command::new(&cmd_tag_owned)
        .args(&args)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow!("Failed to spawn {}. Is it installed? Error: {}.", cmd_tag_owned.clone(), e))?;

    let lines = read_child_output_to_vec(&mut child, child_stream).await?;
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

    let version_parts:  Vec<_>  = version_string.split(".").collect();
    let major_version = version_parts[0];
    let major_version_digits: String = major_version.chars().filter(|c| c.is_digit(10)).collect();
    let mut minor_version = "0";
    if version_parts.len() > 1 {
        minor_version = version_parts[1];  // To avoid problems with things like version 2.15.1, only going to track major/minor a la 2.15
    }
    let minor_version_digits: String = minor_version.chars().filter(|c| c.is_digit(10)).collect();
    let version_string_formatted = format!("{}.{}", major_version_digits, minor_version_digits);
    let version = version_string_formatted.parse::<f32>()?;
    Ok(version)
}

pub mod fastp {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use crate::config::defs::{FASTP_TAG,  RunConfig};
    use crate::utils::streams::{ChildStream};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::file::file_path_manipulator;

    #[derive(Debug)]
    pub struct FastpConfig {
        pub command_fields: HashMap<String, Option<String>>,
    }
    pub struct FastpArgGenerator;

    pub async fn fastp_presence_check() -> Result<f32> {
        let version = version_check(FASTP_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout).await?;
        Ok(version)
    }

    impl ArgGenerator for FastpArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            // eprintln!("Allocating {} threads for fastp", RunConfig::thread_allocation(run_config, FASTP_TAG, None));

            let config = extra
                .and_then(|e| e.downcast_ref::<FastpConfig>())
                .ok_or_else(|| anyhow!("FASTP requires a FaspConfig as extra argument"))?;

            let args = &run_config.args;
            let mut args_vec: Vec<String> = Vec::new();
            args_vec.push("--stdin".to_string());
            args_vec.push("--stdout".to_string());
            args_vec.push("--interleaved_in".to_string());
            args_vec.push("-q".to_string());
            args_vec.push(args.quality.to_string());

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

    pub async fn pigz_presence_check() -> anyhow::Result<f32> {
        let version = version_check(PIGZ_TAG,vec!["--version"], 0, 1 , ChildStream::Stdout).await?;
        Ok(version)
    }

    impl ArgGenerator for PigzArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, _extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            // eprintln!("Allocating {} threads for pigz", RunConfig::thread_allocation(run_config, PIGZ_TAG, None));
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
    use crate::config::defs::{H5DUMP_TAG, PIGZ_TAG};
    use crate::utils::command::version_check;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub async fn h5dump_presence_check() -> anyhow::Result<f32> {
        let version = version_check(H5DUMP_TAG,vec!["-V"], 0, 2 , ChildStream::Stdout).await?;
        Ok(version)
    }
}


pub mod minimap2 {
    use std::collections::HashMap;
    use super::*;
    use tokio::fs;

    pub struct Minimap2Config {
        pub minimap2_index_path: PathBuf,
        pub option_fields: HashMap<String, Option<String>>,
    }
    pub struct Minimap2ArgGenerator;

    pub async fn minimap2_presence_check() -> anyhow::Result<f32> {
        let version = version_check(MINIMAP2_TAG, vec!["--version"], 0, 0, ChildStream::Stdout).await?;
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
                    .map(|ext| FASTA_EXTS.iter().any(|&fasta_ext| ext.eq_ignore_ascii_case(fasta_ext)))
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
                let temp_fasta_path_clone = temp_fasta_path.clone();
                let ref_type_owned = ref_type.to_string();
                let copy_task = tokio::spawn(async move {
                    fs::copy(&sequence_path, &temp_fasta_path_clone)
                        .await
                        .map_err(|e| anyhow!("Failed to copy {} FASTA: {}", ref_type_owned, e))?;
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
                eprintln!("{} mmi created", ref_type);
                tasks.push(err_task);
                (index_temp_path, Some(index_temp))
            }
        };

        Ok((ref_fasta_path, index_path, ref_temp, index_temp, tasks))
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
            let num_cores: usize = RunConfig::thread_allocation(run_config, MINIMAP2_TAG, None);
            args_vec.push("-t".to_string());
            args_vec.push(num_cores.to_string());

            for (key, value) in config.option_fields.iter() {
                args_vec.push(format!("{}", key));
                if let Some(v) = value {
                    args_vec.push(format!("{}", v));
                }
            }
            args_vec.push(index_path.to_string_lossy().to_string());
            args_vec.push("-".to_string()); // Query from stdin

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

    pub async fn samtools_presence_check() -> anyhow::Result<f32> {
        let version = version_check(SAMTOOLS_TAG, vec!["--version"], 0, 1, ChildStream::Stdout).await?;
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

    pub async fn bcftools_presence_check() -> anyhow::Result<f32> {
        let version = version_check(BCFTOOLS_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout).await?;
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

    pub async fn kraken2_presence_check() -> anyhow::Result<f32> {
        let version = version_check(KRAKEN2_TAG,vec!["--version"], 0, 2 , ChildStream::Stdout).await?;
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

    pub async fn ivar_presence_check() -> anyhow::Result<f32> {
        let version = version_check(IVAR_TAG,vec!["version"], 0, 2 , ChildStream::Stdout).await?;
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
    use crate::config::defs::{IVAR_TAG, MUSCLE_TAG};
    use crate::utils::command::version_check;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub struct MuscleArgGenerator;
    pub async fn muscle_presence_check() -> anyhow::Result<f32> {
        let version = version_check(MUSCLE_TAG,vec!["-version"], 0, 1 , ChildStream::Stdout).await?;
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
    pub async fn mafft_presence_check() -> anyhow::Result<f32> {
        let version = version_check(MAFFT_TAG,vec!["--version"], 0, 0 , ChildStream::Stderr).await?;
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

    pub async fn quast_presence_check() -> anyhow::Result<f32> {
        let version = version_check(QUAST_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout).await?;
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

    pub async fn nucmer_presence_check() -> anyhow::Result<f32> {
        let version = version_check(NUCMER_TAG,vec!["--version"], 0, 1 , ChildStream::Stdout).await?;
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

    pub async fn seqkit_presence_check() -> anyhow::Result<f32> {
        let version = version_check(SEQKIT_TAG,vec!["--help"], 1, 1 , ChildStream::Stdout).await?;
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
                    args_vec.push("-".to_string());
                }
                SeqkitSubcommand::Grep => {
                    args_vec.push("grep".to_string());
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

pub mod bowtie2 {
    use std::path::{Path, PathBuf};
    use std::collections::HashMap;
    use anyhow::{anyhow, Result};
    use std::fs::{self, DirEntry};
    use std::process::Command;
    use crate::config::defs::{RunConfig, BOWTIE2_TAG};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::streams::ChildStream;

    pub struct Bowtie2Config {
        pub bt2_index_path: PathBuf,
        pub option_fields: HashMap<String, Option<String>>,
    }

    pub struct Bowtie2ArgGenerator;

    pub async fn bowtie2_presence_check() -> anyhow::Result<f32> {
        let version = version_check(BOWTIE2_TAG, vec!["-h"], 0, 3, ChildStream::Stdout).await?;
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
        // eprintln!("Bowtie2: Created unpack directory: {}", unpack_dir.display());

        let mut index_dir = unpack_dir.clone();

        if input.is_dir() {
            index_dir = input.to_path_buf();
        } else if input.exists() {
            let ext = input.extension().map(|s| s.to_str().unwrap_or("")).unwrap_or("");
            let stem_ext = input.file_stem().and_then(|s| Path::new(s).extension()).map(|s| s.to_str().unwrap_or("")).unwrap_or("");

            let unpack_cmd = if ext == "gz" && stem_ext == "tar" {
                eprintln!("Bowtie2: Unpacking GZIP TAR: {}", input.display());
                Command::new("tar")
                    .arg("xzf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else if ext == "tar" {
                eprintln!("Bowtie2: Unpacking TAR: {}", input.display());
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
                eprintln!("Bowtie2: Stripped to subdir: {}", subdir.display());
                index_dir = subdir;
            }
        }
        // eprintln!("Bowtie2 index directory: {}", index_dir.display());

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
        eprintln!("Bowtie2 derived basename: {}", basename);

        let final_path = index_dir.join(&basename);
        eprintln!("Final Bowtie2 index path: {}", final_path.display());
        Ok(final_path)
    }

    impl ArgGenerator for Bowtie2ArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<Bowtie2Config>())
                .ok_or_else(|| anyhow!("Bowtie2 requires a Bowtie2Config as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            for (key, value) in config.option_fields.iter() {
                args_vec.push(format!("{}", key));
                if let Some(v) = value {
                    args_vec.push(format!("{}", v));
                }
            }

            args_vec.push("-x".to_string());
            args_vec.push(config.bt2_index_path.to_string_lossy().to_string());
            args_vec.push("-p".to_string());
            let num_cores: usize = RunConfig::thread_allocation(run_config, BOWTIE2_TAG, None);
            args_vec.push(num_cores.to_string());
            args_vec.push("--interleaved".to_string());
            args_vec.push("-".to_string());

            Ok(args_vec)
        }
    }
}

pub mod hisat2 {
    use std::path::{Path, PathBuf};
    use std::collections::HashMap;
    use anyhow::{anyhow, Result};
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

    pub async fn hisat2_presence_check() -> anyhow::Result<f32> {
        let version = version_check(HISAT2_TAG, vec!["--version"], 0, 2, ChildStream::Stdout).await?;
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
        eprintln!("HISAT2: Created unpack directory: {}", unpack_dir.display());

        let mut index_dir = unpack_dir.clone();

        // Handle input type
        if input.is_dir() {
            index_dir = input.to_path_buf();
        } else if input.exists() {
            let ext = input.extension().map(|s| s.to_str().unwrap_or("")).unwrap_or("");
            let stem_ext = input.file_stem().and_then(|s| Path::new(s).extension()).map(|s| s.to_str().unwrap_or("")).unwrap_or("");

            let unpack_cmd = if ext == "gz" && stem_ext == "tar" {
                eprintln!("HISAT2: Unpacking GZIP TAR: {}", input.display());
                Command::new("tar")
                    .arg("xzf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else if ext == "tar" {
                eprintln!("HISAT2: Unpacking TAR: {}", input.display());
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
                eprintln!("HISAT2: Stripped to subdir: {}", subdir.display());
                index_dir = subdir;
            }
        }
        eprintln!("HISAT2 index directory: {}", index_dir.display());

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
        eprintln!("HISAT2 derived basename: {}", basename);

        let final_path = index_dir.join(&basename);
        eprintln!("Final HISAT2 index path: {}", final_path.display());
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

            args_vec.push("-S".to_string());
            args_vec.push("-".to_string());

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

    pub async fn kallisto_presence_check() -> anyhow::Result<f32> {
        let version = version_check(KALLISTO_TAG, vec!["version"], 0, 2, ChildStream::Stdout).await?;
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

    pub async fn star_presence_check() -> anyhow::Result<f32> {
        let version = version_check(STAR_TAG, vec!["--version"], 0, 0, ChildStream::Stdout).await?;
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
        eprintln!("STAR: Created unpack directory: {}", unpack_dir.display());

        let mut index_dir = unpack_dir.clone();

        // Handle input type
        if input.is_dir() {
            index_dir = input.to_path_buf();
        } else if input.exists() {
            let ext = input.extension().map(|s| s.to_str().unwrap_or("")).unwrap_or("");
            let stem_ext = input.file_stem().and_then(|s| Path::new(s).extension()).map(|s| s.to_str().unwrap_or("")).unwrap_or("");

            let unpack_cmd = if ext == "gz" && stem_ext == "tar" {
                eprintln!("STAR: Unpacking GZIP TAR: {}", input.display());
                Command::new("tar")
                    .arg("xzf")
                    .arg(input)
                    .arg("-C")
                    .arg(&unpack_dir)
                    .output()?
            } else if ext == "tar" {
                eprintln!("STAR: Unpacking TAR: {}", input.display());
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
                eprintln!("STAR: Stripped to subdir: {}", subdir.display());
                index_dir = subdir;
            }
        }
        eprintln!("STAR index directory: {}", index_dir.display());

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
        _ => return Err(anyhow!("Unknown tool: {}", tool)),
    };

    generator.generate_args(run_config, extra)
}

pub async fn check_versions(tools: Vec<&str>) -> Result<()> {
    let checks = tools.into_iter().map(|tool| async move {
        let version = match tool {
            FASTP_TAG => fastp::fastp_presence_check().await,
            H5DUMP_TAG => h5dump::h5dump_presence_check().await,
            MINIMAP2_TAG => minimap2::minimap2_presence_check().await,
            SAMTOOLS_TAG => samtools::samtools_presence_check().await,
            KRAKEN2_TAG => kraken2::kraken2_presence_check().await,
            BCFTOOLS_TAG => bcftools::bcftools_presence_check().await,
            IVAR_TAG => ivar::ivar_presence_check().await,
            MUSCLE_TAG => muscle::muscle_presence_check().await,
            MAFFT_TAG => mafft::mafft_presence_check().await,
            QUAST_TAG => quast::quast_presence_check().await,
            NUCMER_TAG => nucmer::nucmer_presence_check().await,
            SEQKIT_TAG => seqkit::seqkit_presence_check().await,
            BOWTIE2_TAG => bowtie2::bowtie2_presence_check().await,
            HISAT2_TAG => hisat2::hisat2_presence_check().await,
            KALLISTO_TAG => kallisto::kallisto_presence_check().await,
            STAR_TAG => star::star_presence_check().await,

            _ => return Err(anyhow!("Unknown tool: {}", tool)),
        }?;
        Ok((tool.to_string(), version))
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
            println!("Warning: No minimum version specified for tool: {}", tool);
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