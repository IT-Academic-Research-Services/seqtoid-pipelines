use std::collections::HashMap;
/// Functions and structs for working with creating command-line arguments

use anyhow::{anyhow, Result};
use num_cpus;
use crate::config::defs::{FASTP_TAG, PIGZ_TAG, H5DUMP_TAG, MINIMAP2_TAG, SAMTOOLS_TAG, KRAKEN2_TAG, BCFTOOLS_TAG, IVAR_TAG, MUSCLE_TAG};
use crate::cli::Arguments;


pub trait ArgGenerator {
    fn generate_args(&self, args: &Arguments, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>>; 
}

mod fastp {
    use std::path::PathBuf;
    use anyhow::{anyhow, Result};
    use tokio::process::Command;
    use crate::cli::Arguments;
    use crate::config::defs::FASTP_TAG;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::ArgGenerator;
    use crate::utils::file::file_path_manipulator;

    pub struct FastpArgGenerator;
    
    pub async fn fastp_presence_check() -> Result<String> {
        let args: Vec<&str> = vec!["-v"];

        let cmd_tag_owned = FASTP_TAG.to_string();
        let mut child = Command::new(FASTP_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn {}: {}. Is fastp installed?", cmd_tag_owned, e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from fastp -v"))?;
        let version = first_line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| anyhow!("Invalid fastp -v output: {}", first_line))?
            .to_string();
        if version.is_empty() {
            return Err(anyhow!("Empty version number in fastp -v output: {}", first_line));
        }
        Ok(version)
    }
    
    impl ArgGenerator for FastpArgGenerator {
        fn generate_args(&self, args: &Arguments, _extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let mut args_vec: Vec<String> = Vec::new();
            args_vec.push("--stdin".to_string());
            args_vec.push("--stdout".to_string());
            args_vec.push("--interleaved_in".to_string());
            args_vec.push("-q".to_string());
            args_vec.push(args.quality.to_string());
            args_vec.push("-w".to_string());
            args_vec.push(args.threads.to_string());

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
    use crate::cli::Arguments;
    use crate::utils::command::ArgGenerator;
    pub struct PigzArgGenerator;

    impl ArgGenerator for PigzArgGenerator {
        fn generate_args(&self, args: &Arguments, _extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let mut args_vec: Vec<String> = Vec::new();
            args_vec.push("-c".to_string());
            args_vec.push("-p".to_string());
            args_vec.push(args.threads.to_string());
            Ok(args_vec)
        }
    }
}

mod h5dump {
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::config::defs::{H5DUMP_TAG};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    
    pub async fn h5dump_presence_check() -> anyhow::Result<String> {
        let args: Vec<&str> = vec!["-V"];


        let mut child = Command::new(H5DUMP_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is hddump installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from h5dump -V"))?;
        let version = first_line
            .split_whitespace()
            .nth(2)
            .ok_or_else(|| anyhow!("Invalid h5dump -V output: {}", first_line))?
            .to_string();
        if version.is_empty() {
            return Err(anyhow!("Empty version number in h5dump -V output: {}", first_line));
        }
        Ok(version)
    }
}

mod minimap2 {
    use std::path::PathBuf;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::cli::{Arguments, Technology};
    use crate::config::defs::MINIMAP2_TAG;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::ArgGenerator;

    pub struct Minimap2ArgGenerator;
    pub async fn minimap2_presence_check() -> anyhow::Result<String> {
        let args: Vec<&str> = vec!["-V"];

        let mut child = Command::new(MINIMAP2_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is minimap2 installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from minimap2 -V"))?;
        let version = first_line
            .split_whitespace()
            .nth(0)
            .ok_or_else(|| anyhow!("Invalid minimap2 -V output: {}", first_line))?
            .to_string();
        if version.is_empty() {
            return Err(anyhow!("Empty version number in minimap2 -V output: {}", first_line));
        }
        Ok(version)
        
    }

    impl ArgGenerator for Minimap2ArgGenerator {
        fn generate_args(&self, args: &Arguments, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let paths = extra
                .and_then(|e| e.downcast_ref::<(PathBuf, PathBuf)>())
                .ok_or_else(|| anyhow!("Minimap2 requires (ref_pipe_path, query_pipe_path) as extra arguments"))?;

            let (ref_pipe_path, query_pipe_path) = paths;
            
            let mut args_vec: Vec<String> = Vec::new();

            let num_cores: usize = match args.limit_align_threads {
                true => args.threads,
                false => num_cpus::get()-1,
            };
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

    pub fn arg_generator(args: &Arguments, ref_pipe_path: &PathBuf, query_pipe_path: &PathBuf) -> Vec<String> {
        let mut args_vec: Vec<String> = Vec::new();

        let num_cores : usize = match &args.limit_align_threads {
            true => {
                args.threads
            }
            false => {
                num_cpus::get()
            }
        };

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
    use crate::cli::{Arguments};
    use crate::config::defs::{SAMTOOLS_TAG, SamtoolsSubcommand};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::ArgGenerator;
    
    #[derive(Debug)]
    pub struct SamtoolsConfig {
        pub subcommand: SamtoolsSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }
    
    pub struct SamtoolsArgGenerator;




    pub async fn samtools_presence_check() -> anyhow::Result<String> {
        let args: Vec<&str> = vec!["--version"];

        let mut child = Command::new(SAMTOOLS_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is samtools installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from samtools --version"))?;
        let version = first_line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| anyhow!("Invalid samtools --version output: {}", first_line))?
            .to_string();
        if version.is_empty() {
            return Err(anyhow!("Empty version number in samtools --version output: {}", first_line));
        }
        Ok(version)
    }

    impl ArgGenerator for SamtoolsArgGenerator {
        fn generate_args(&self, args: &Arguments, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<SamtoolsConfig>())
                .ok_or_else(|| anyhow!("Samtools requires a SamtoolsConfig as extra argument"))?;

            let mut args_vec: Vec<String> = Vec::new();

            match config.subcommand {
                SamtoolsSubcommand::View => {
                    args_vec.push("view".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(args.threads.to_string());
                    args_vec.push("--no-PG".to_string());
                
                    }
                SamtoolsSubcommand::Fastq => {
                    args_vec.push("fastq".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(args.threads.to_string());
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
                    args_vec.push(args.threads.to_string());
                    }
                SamtoolsSubcommand::Index => {
                    args_vec.push("index".to_string());
                    args_vec.push("-@".to_string());
                    args_vec.push(args.threads.to_string());
                }
                SamtoolsSubcommand::Mpileup => {
                    args_vec.push("mpileup".to_string());
                    args_vec.push("-A".to_string()); // do not discard anomalous read pairs
                    args_vec.push("-d".to_string()); // max depth zero
                    args_vec.push("0".to_string());
                    args_vec.push("-Q".to_string()); // skip bases with baseQ/BAQ smaller than INT [13]
                    args_vec.push("0".to_string());
                    // args_vec.push(args.quality.to_string());

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
                    // args_vec.push("-l".to_string());
                    // args_vec.push("50".to_string());
                    

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
    use crate::cli::Arguments;
    use crate::config::defs::{BcftoolsSubcommand,  BCFTOOLS_TAG};
    use crate::utils::command::ArgGenerator;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[derive(Debug)]
    pub struct BcftoolsConfig {
        pub subcommand: BcftoolsSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct BcftoolsArgGenerator;

    pub async fn bcftools_presence_check() -> anyhow::Result<String> {
        let args: Vec<&str> = vec!["-v"];

        let mut child = Command::new(BCFTOOLS_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is samtools installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from samtools --version"))?;
        let version = first_line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| anyhow!("Invalid samtools --version output: {}", first_line))?
            .to_string();
        if version.is_empty() {
            return Err(anyhow!("Empty version number in samtools --version output: {}", first_line));
        }
        Ok(version)
    }

    impl ArgGenerator for BcftoolsArgGenerator {
        fn generate_args(&self, args: &Arguments, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
            let config = extra
                .and_then(|e| e.downcast_ref::<BcftoolsConfig>())
                .ok_or_else(|| anyhow!("Samtools requires a BcftoolsConfig as extra argument"))?;

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
    use crate::cli::Arguments;
    use crate::config::defs::{KRAKEN2_TAG};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};
    use crate::utils::command::ArgGenerator;
        use crate::utils::file::file_path_manipulator;

        #[derive(Debug)]
    pub struct Kraken2Config {
        pub report_path: PathBuf, 
        pub classified_path: PathBuf,
        pub fastq_path: PathBuf,  
    }

    pub struct Kraken2ArgGenerator;

    pub async fn kraken2_presence_check() -> anyhow::Result<String> {

        let args: Vec<&str> = vec!["--version"];

        let mut child = Command::new(KRAKEN2_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is samtools installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from samtools --version"))?;

        let version = first_line
            .split_whitespace()
            .nth(2)
            .ok_or_else(|| anyhow!("Invalid samtools --version output: {}", first_line))?
            .to_string();

        if version.is_empty() {
            return Err(anyhow!("Empty version number in samtools --version output: {}", first_line));
        }
        Ok(version)
    }
    
    impl ArgGenerator for Kraken2ArgGenerator {
        fn generate_args(&self, args: &Arguments, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
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
            
            let num_cores: usize = match args.limit_align_threads {
                true => args.threads,
                false => num_cpus::get()-1,
            };
            args_vec.push("--threads".to_string());
            args_vec.push(num_cores.to_string());
            args_vec.push("--report".to_string());
            args_vec.push( config.report_path.to_string_lossy().to_string());
            args_vec.push("--classified-out".to_string());
            args_vec.push(config.classified_path.to_string_lossy().to_string());
            args_vec.push("--output".to_string());
            args_vec.push("-".to_string()); // "-" will suppress normal output
            // args_vec.push("--memory-mapping".to_string());
            // args_vec.push("--gzip-compressed".to_string());

            args_vec.push(config.fastq_path.to_string_lossy().to_string()); // input, should be a mkfifo

            Ok(args_vec)
        }
    }
}

pub mod ivar {
    use std::collections::HashMap;
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::cli::Arguments;
    use crate::config::defs::{IVAR_TAG, IvarSubcommand, IVAR_QUAL_THRESHOLD, IVAR_FREQ_THRESHOLD};
    use crate::utils::command::ArgGenerator;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[derive(Debug)]
    pub struct IvarConfig {
        pub subcommand: IvarSubcommand,
        pub subcommand_fields: HashMap<String, Option<String>>,
    }

    pub struct IvarArgGenerator;

    pub async fn ivar_presence_check() -> anyhow::Result<String> {
        let args: Vec<&str> = vec!["version"];

        let mut child = Command::new(IVAR_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is samtools installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from samtools --version"))?;

        let version = first_line
            .split_whitespace()
            .nth(2)
            .ok_or_else(|| anyhow!("Invalid samtools --version output: {}", first_line))?
            .to_string();

        if version.is_empty() {
            return Err(anyhow!("Empty version number in samtools --version output: {}", first_line));
        }
        Ok(version)
    }

    impl ArgGenerator for IvarArgGenerator {
        fn generate_args(&self, args: &Arguments, extra: Option<&dyn std::any::Any>) -> anyhow::Result<Vec<String>> {
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
    use crate::config::defs::MUSCLE_TAG;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    pub struct MuscleArgGenerator;
    pub async fn muscle_presence_check() -> anyhow::Result<String> {

        let args: Vec<&str> = vec!["-version"];
        
        let mut child = Command::new(MUSCLE_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is muscle installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from muscle -version"))?;

        let version = first_line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| anyhow!("Invalid muscle -versionoutput: {}", first_line))?
            .to_string();

        if version.is_empty() {
            return Err(anyhow!("Empty version number in muscle -version output: {}", first_line));
        }
        Ok(version)
    }
}

pub fn generate_cli(tool: &str, args: &Arguments, extra: Option<&dyn std::any::Any>) -> Result<Vec<String>> {
    let generator: Box<dyn ArgGenerator> = match tool {
        FASTP_TAG => Box::new(fastp::FastpArgGenerator),
        PIGZ_TAG => Box::new(pigz::PigzArgGenerator),
        MINIMAP2_TAG => Box::new(minimap2::Minimap2ArgGenerator),
        SAMTOOLS_TAG => Box::new(samtools::SamtoolsArgGenerator),
        KRAKEN2_TAG => Box::new(kraken2::Kraken2ArgGenerator),
        BCFTOOLS_TAG => Box::new(bcftools::BcftoolsArgGenerator),
        H5DUMP_TAG => return Err(anyhow!("h5dump argument generation not implemented")),
        _ => return Err(anyhow!("Unknown tool: {}", tool)),
    };

    generator.generate_args(args, extra)
}


pub async fn check_version(tool: &str) -> Result<String> {
    let version = match tool {
        FASTP_TAG => fastp::fastp_presence_check().await,
        H5DUMP_TAG => h5dump::h5dump_presence_check().await,
        MINIMAP2_TAG => {minimap2::minimap2_presence_check().await},
        SAMTOOLS_TAG => {samtools::samtools_presence_check().await},
        KRAKEN2_TAG => kraken2::kraken2_presence_check().await,
        BCFTOOLS_TAG => bcftools::bcftools_presence_check().await,
        IVAR_TAG => ivar::ivar_presence_check().await,
        MUSCLE_TAG => muscle::muscle_presence_check().await,
        _ => return Err(anyhow!("Unknown tool: {}", tool)),
    };
    Ok(version?)
}


pub async fn check_versions(tools: Vec<&str>) -> Result<HashMap<String, String>>{
    let mut versions: HashMap<String, String> = HashMap::new();
    for tool in tools {
        let _tool_version = match check_version(tool).await {
            Ok(version) => {
                versions.insert(tool.to_string(), version);
            }
            Err(err) => {
                return Err(anyhow!("Cannot find external tool in path: {}. Error: {}", tool, err));
            }
        };
    }
    Ok(versions)
}