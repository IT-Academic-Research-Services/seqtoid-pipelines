//! Functions and structs for handling AWS workers

use anyhow::{anyhow, Context, Result};
use aws_sdk_ec2::{
    types::{Filter, Instance},
    Client as Ec2Client,
};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use tokio::fs;
