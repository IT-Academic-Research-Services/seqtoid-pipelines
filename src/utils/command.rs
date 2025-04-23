/// Functions and structs for working with creating command-line arguments

use anyhow::Result;
use crate::{FASTP_TAG, PIGZ_TAG};
use crate::utils::Arguments;


// fn fastp_presence_check() -> bool {
// 
//     let mut cmd = String::new();
//     cmd.push_str(FASTP_TAG);
// 
//     cmd.push_str(" -v");
//     
// 
//     
// }

mod fastp {
    use crate::utils::Arguments;

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

        args_vec.push("-p".to_string());
        args_vec.push(args.threads.to_string());
        
        args_vec
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
