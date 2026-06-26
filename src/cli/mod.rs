use clap::Parser;
pub mod args;
pub use args::{Arguments, TargetType, Technology};

pub fn parse() -> Arguments {
    Arguments::parse()
}