//! Functions and structs for handling AWS workers

use anyhow::{anyhow, Context, Result};
use aws_sdk_ec2::{
    types::{Filter, Instance},
    Client as Ec2Client,
};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use tokio::fs;


pub const WORKER_ROLE_TAG: &str = "Role";
pub const WORKER_ROLE_VALUE: &str = "seqtoid-nr-cpu-worker";
pub const WORKER_BACKEND_TAG: &str = "Backend";
pub const WORKER_BACKEND_VALUE: &str = "cpu";
pub const WORKER_REFERENCE_SET_TAG: &str = "ReferenceSet";
pub const WORKER_REFERENCE_SET_VALUE: &str = "phase2";
pub const WORKER_ENVIRONMENT_TAG: &str = "Environment";
pub const WORKER_ENVIRONMENT_VALUE: &str = "dev";