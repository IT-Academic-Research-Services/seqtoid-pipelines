use clap::Parser;
pub mod args;
pub use args::{Arguments, Technology};

pub fn parse() -> Arguments {
    Arguments::parse()
}