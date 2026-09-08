// src/lib.rs
pub mod config;
pub mod utils;
pub mod pipelines;
pub mod cli;
pub mod distributed;

pub use cli::{Arguments, Technology};
