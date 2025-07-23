/// Functions and structs for working with creating command-line arguments

use std::collections::HashMap;
use anyhow::{anyhow, Result};
use num_cpus;
use tokio::process::Command;
use futures::future::try_join_all;
use crate::config::defs::{RunConfig, TOOL_VERSIONS, FASTP_TAG, PIGZ_TAG, H5DUMP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, IVAR_TAG, MUSCLE_TAG, MAFFT_TAG, QUAST_TAG, NUCMER_TAG, SHOW_COORDS_TAG, SEQKIT_TAG};
use crate::cli::Arguments;
use crate::utils::streams::{read_child_output_to_vec, ChildStream};

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

mod fastp {
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use tokio::process::Command;
    use crate::config::defs::{FASTP_TAG, KRAKEN2_TAG, RunConfig};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::{version_check, ArgGenerator};
    use crate::utils::file::file_path_manipulator;

    pub struct FastpArgGenerator;

    pub async fn fastp_presence_check() -> Result<f32> {
        let version = version_check(FASTP_TAG,vec!["-v"], 0, 1 , ChildStream::Stdout).await?;
        Ok(version)
    }

    impl ArgGenerator for FastpArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, _extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let mut args_vec: Vec<String> = Vec::new();
            args_vec.push("--stdin".to_string());
            args_vec.push("--stdout".to_string());
            args_vec.push("--interleaved_in".to_string());
            args_vec.push("-q".to_string());
            args_vec.push(args.quality.to_string());
            args_vec.push("-w".to_string());
            args_vec.push(RunConfig::thread_allocation(run_config, FASTP_TAG, None).to_string());

            if let Some(adapter_fasta) = &args.adapter_fasta {
                let cwd = std::env::current_dir()?;
                let adapter_path = file_path_manipulator(&PathBuf::from(adapter_fasta), &cwd.clone(), None, None, "");
                if !adapter_path.exists() {
                    return Err(anyhow!("Adapter FASTA file does not exist: {}", adapter_path.display()));
                }
                args_vec.push("--adapter_fasta".to_string());
                args_vec.push(adapter_path.to_string_lossy().into_owned());
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

mod minimap2 {
    use std::path::PathBuf;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::cli::{Technology};
    use crate::config::defs::{MINIMAP2_TAG, RunConfig};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::{version_check, ArgGenerator};

    pub struct Minimap2ArgGenerator;
    pub async fn minimap2_presence_check() -> anyhow::Result<f32> {
        let version = version_check(MINIMAP2_TAG,vec!["--version"], 0, 0 , ChildStream::Stdout).await?;
        Ok(version)
    }

    impl ArgGenerator for Minimap2ArgGenerator {
        fn generate_args(&self, run_config: &RunConfig, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let args = &run_config.args;
            let paths = extra
                .and_then(|e| e.downcast_ref::<(PathBuf, PathBuf)>())
                .ok_or_else(|| anyhow!("Minimap2 requires (ref_pipe_path, query_pipe_path) as extra arguments"))?;

            let (ref_pipe_path, query_pipe_path) = paths;

            let mut args_vec: Vec<String> = Vec::new();

            let num_cores: usize = RunConfig::thread_allocation(run_config, MINIMAP2_TAG, None);
            args_vec.push("-t".to_string());
            args_vec.push(num_cores.to_string());

            args_vec.push("-ax".to_string());
            let technology = args.technology.clone();
            match technology {
                Technology::Illumina => {
                    args_vec.push("sr".to_string());
                }
                Technology::ONT => {
                    args_vec.push("map-ont".to_string());
                }
            }

            args_vec.push(ref_pipe_path.to_string_lossy().to_string());
            args_vec.push(query_pipe_path.to_string_lossy().to_string());

            Ok(args_vec)
        }
    }

    pub fn arg_generator(run_config: &RunConfig, ref_pipe_path: &PathBuf, query_pipe_path: &PathBuf) -> Vec<String> {
        let args = &run_config.args;
        let mut args_vec: Vec<String> = Vec::new();

        let num_cores: usize = RunConfig::thread_allocation(run_config, MINIMAP2_TAG, None);
        args_vec.push("-t".to_string());
        args_vec.push(num_cores.to_string());

        let technology = args.technology.clone();
        match technology {
            Technology::Illumina => {
                args_vec.push("-ax sr".to_string());
            }
            Technology::ONT => {
                args_vec.push("-ax map-ont".to_string());
            }
        }

        args_vec.push(ref_pipe_path.to_string_lossy().to_string());
        args_vec.push(query_pipe_path.to_string_lossy().to_string());

        args_vec
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
        let version = version_check(SAMTOOLS_TAG,vec!["--version"], 0, 1 , ChildStream::Stdout).await?;
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
                }
                SamtoolsSubcommand::Index => {
                    args_vec.push("index".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(RunConfig::thread_allocation(run_config, SAMTOOLS_TAG, Some("index")).to_string());
                }
                SamtoolsSubcommand::Mpileup => {
                    args_vec.push("mpileup".to_string());
                    args_vec.push("-A".to_string()); // do not discard anomalous read pairs
                    args_vec.push("-d".to_string()); // max depth zero
                    args_vec.push("0".to_string());
                    args_vec.push("-Q".to_string()); // skip bases with baseQ/BAQ smaller than INT [13]
                    args_vec.push("0".to_string());
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

pub mod bcftools {
    use std::collections::HashMap;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{BcftoolsSubcommand, BCFTOOLS_TAG, RunConfig};
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
                    args_vec.push("-O".to_string());
                    args_vec.push("b".to_string());  // Compressed BCF output to save space in stream
                }
                BcftoolsSubcommand::Mpileup => {
                    args_vec.push("mpileup".to_string());
                    args_vec.push("-a".to_string());
                    args_vec.push("AD".to_string()); // include allele depth (AD) for all positions, including those with zero coverage, mimicking -aa.
                    args_vec.push("-d".to_string());
                    args_vec.push("100000000".to_string()); // max depth essentially without limit
                    args_vec.push("-L".to_string());
                    args_vec.push("100000000".to_string()); // max per-file depth essentially without limit
                    args_vec.push("-Q".to_string());
                    args_vec.push(args.quality.to_string()); // skip bases with baseQ/BAQ smaller than INT [13]
                    args_vec.push("-O".to_string());
                    args_vec.push("b".to_string());  // Compressed BCF output to save space in stream
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
                    let kraken2_db_path = file_path_manipulator(&PathBuf::from(db), &cwd, None, None, "");
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
            args_vec.push("quast".to_string());
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