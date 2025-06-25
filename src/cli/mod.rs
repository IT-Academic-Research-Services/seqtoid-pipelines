use clap::Parser;
pub mod args;
pub use args::{Arguments, Technology, TargetType};

pub fn parse() -> Arguments {
    Arguments::parse()
}