/// Functions and structs for working with creating command-line arguments


use crate::FASTP_TAG;
use std::collections::{HashMap, HashSet};
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

fn fastp_arg_generator(args: &Arguments) -> Vec<String> {
    
    
    let mut args = Vec::new();
    args.push("--stdin");
    args.push("--stdout");
    args.push("--interleaved_in");
    args.push("-q");
    args.push(args.quality.to_string().as_str());
    args.push("-w");
    args.push(args.threads.to_string().as_str());

    
    args 
}



pub fn generate_cli(tool: &str, args: &Arguments) -> Result<String, String> {
    

    let cmd = match tool {
        FASTP_TAG => fastp_arg_generator(&args),
        _ => return Err(format!("Unknown tool: {}", tool)),
    };
    
    eprintln!("cmd = {:?}", cmd);
    
    Ok(cmd)
}
