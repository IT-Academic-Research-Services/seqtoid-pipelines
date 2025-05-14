/// Functions and structs for working with creating command-line arguments

use anyhow::{anyhow, Result};
use crate::utils::defs::{FASTP_TAG, PIGZ_TAG, H5DUMP_TAG};
use crate::utils::Arguments;



mod fastp {
    use anyhow::{anyhow, Result};
    use tokio::process::Command;
    use crate::utils::Arguments;
    use crate::utils::defs::FASTP_TAG;
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[allow(dead_code)]
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

    pub fn arg_generator(args: &Arguments) -> Vec<String> {
        let mut args_vec: Vec<String> = Vec::new();
        args_vec.push("--stdin".to_string());
        args_vec.push("--stdout".to_string());
        args_vec.push("--interleaved_in".to_string());
        args_vec.push("-q".to_string());
        args_vec.push(args.quality.to_string());
        args_vec.push("-w".to_string());
        args_vec.push(args.threads.to_string());
        args_vec
    }
}

mod pigz {
    use crate::utils::Arguments;
    
    pub fn arg_generator(args: &Arguments) -> Vec<String> {
        let mut args_vec: Vec<String> = Vec::new();
        args_vec.push("-c".to_string());
        args_vec.push("-p".to_string());
        args_vec.push(args.threads.to_string());
        
        args_vec
    }
}

mod h5dump {
    use anyhow::anyhow;
    use tokio::process::Command;
    use crate::utils::defs::{H5DUMP_TAG};
    use crate::utils::streams::{read_child_output_to_vec, ChildStream};

    #[allow(dead_code)]
    pub async fn h5dump_presence_check() -> anyhow::Result<String> {
        let args: Vec<&str> = vec!["-V"];


        let mut child = Command::new(H5DUMP_TAG)
            .args(&args)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| anyhow!("Failed to spawn: {}. Is fastp installed?",  e))?;

        let lines = read_child_output_to_vec(&mut child, ChildStream::Stdout).await?;
        let first_line = lines
            .first()
            .ok_or_else(|| anyhow!("No output from h5dump -V"))?;
        let version = first_line
            .split_whitespace()
            .nth(2)
            .ok_or_else(|| anyhow!("Invalid fastp -v output: {}", first_line))?
            .to_string();
        if version.is_empty() {
            return Err(anyhow!("Empty version number in fastp -v output: {}", first_line));
        }
        Ok(version)
    }
}

pub fn generate_cli(tool: &str, args: &Arguments) -> Result<Vec<String>> {
    

    let cmd = match tool {
        FASTP_TAG => fastp::arg_generator(&args),
        PIGZ_TAG => pigz::arg_generator(&args),
        _ => return Err(anyhow::anyhow!("Unknown tool: {}", tool)),
    };
    
    Ok(cmd)
}


pub async fn check_version(tool: &str) -> Result<String> {
    let version = match tool {
        FASTP_TAG => fastp::fastp_presence_check().await,
        H5DUMP_TAG => h5dump::h5dump_presence_check().await,
        _ => return Err(anyhow!("Unknown tool: {}", tool)),
    };
    Ok(version?)
}
