/// Functions and structs for working with creating command-line arguments


use crate::FASTP_TAG;
use std::collections::{HashMap, HashSet};
use lazy_static::lazy_static;
use crate::utils::Arguments;


lazy_static! {
    static ref FASTP_ARGS: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        m.insert("verbose", "verbose");
        m.insert("threads", "thread");
        m
    };
}

lazy_static! {
    static ref FASTP_DEFAULTS: HashSet<&'static str> = {
        let mut s = HashSet::new();
        s.insert("stdin");
        s.insert("interleaved_in");
        s
    };
}

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

fn fastp_arg_generator(args: &Arguments) -> String {
    
    
    let mut cmd = String::new();
    cmd.push_str(FASTP_TAG);
    cmd.push_str(" --stdin --stdout --interleaved_in");
    cmd.push_str(" -q");
    cmd.push_str(args.quality.to_string().as_str());
    cmd.push_str(" -w");
    cmd.push_str(args.threads.to_string().as_str());
    
    cmd
}

fn arg_selector(tool: &str) -> Result<(&'static HashMap<&'static str, &'static str>, &'static HashSet<&'static str>), String> {

    match tool {
        FASTP_TAG => Ok((&FASTP_ARGS, &FASTP_DEFAULTS)),
        _ => Err(format!("Unknown tool: {}", tool)),
    }
}

pub fn generate_cli(tool: &str, args: &Arguments) -> Result<String, String> {
    

    let cmd = match tool {
        FASTP_TAG => fastp_arg_generator(args),
        _ => return Err(format!("Unknown tool: {}", tool)),
    };
    
    eprintln!("cmd = {:?}", cmd);
    
    Ok(cmd)
}
