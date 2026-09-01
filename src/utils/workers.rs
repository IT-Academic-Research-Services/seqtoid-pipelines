//! Functions and structs for handling AWS workers

use anyhow::{anyhow, Context, Result};
use aws_sdk_ec2::{
    types::{Filter, Instance},
    Client as Ec2Client,
};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use tokio::fs;

use crate::config::defs::{PipelineError};


pub const WORKER_ROLE_TAG: &str = "Role";
pub const WORKER_ROLE_VALUE: &str = "seqtoid-nr-cpu-worker";
pub const WORKER_BACKEND_TAG: &str = "Backend";
pub const WORKER_BACKEND_VALUE: &str = "cpu";
pub const WORKER_REFERENCE_SET_TAG: &str = "ReferenceSet";
pub const WORKER_REFERENCE_SET_VALUE: &str = "phase2";
pub const WORKER_ENVIRONMENT_TAG: &str = "Environment";
pub const WORKER_ENVIRONMENT_VALUE: &str = "dev";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Worker {
    pub instance_id: String,
    pub private_ip: String,
    pub instance_type: String,
    pub availability_zone: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct WorkerStatus {
    pub state: Option<String>,
    pub worker_ready: Option<bool>,
}

#[derive(Clone)]
pub struct WorkerManager {
    ec2: Ec2Client,
    workers_root: PathBuf,
}

impl WorkerManager {
    pub async fn new(workers_root: PathBuf) -> Result<Self> {
        let sdk_config =
            aws_config::load_defaults(aws_config::BehaviorVersion::latest()).await;

        Ok(Self {
            ec2: Ec2Client::new(&sdk_config),
            workers_root,
        })
    }

    /// Discovers currently running EC2 instances configured as distributed workers.
    ///
    /// Queries EC2 for running instances matching the configured worker role,
    /// backend, reference set, and environment tags. Returned workers are sorted
    /// by instance ID for deterministic ordering.
    ///
    /// # Arguments
    ///
    /// * `self` - Worker manager containing the configured EC2 client
    ///
    /// # Returns
    /// Vector of running, correctly tagged workers
    pub async fn discover_running_workers(&self) -> Result<Vec<Worker>> {
        let filters = vec![
            Filter::builder()
                .name("instance-state-name")
                .values("running")
                .build(),
            Filter::builder()
                .name(format!("tag:{}", WORKER_ROLE_TAG))
                .values(WORKER_ROLE_VALUE)
                .build(),
            Filter::builder()
                .name(format!("tag:{}", WORKER_BACKEND_TAG))
                .values(WORKER_BACKEND_VALUE)
                .build(),
            Filter::builder()
                .name(format!("tag:{}", WORKER_REFERENCE_SET_TAG))
                .values(WORKER_REFERENCE_SET_VALUE)
                .build(),
            Filter::builder()
                .name(format!("tag:{}", WORKER_ENVIRONMENT_TAG))
                .values(WORKER_ENVIRONMENT_VALUE)
                .build(),
        ];

        let mut next_token: Option<String> = None;
        let mut workers = Vec::new();

        loop {
            let mut request = self.ec2
                .describe_instances()
                .set_filters(Some(filters.clone()));

            if let Some(token) = next_token.clone() {
                request = request.next_token(token);
            }

            let response = request
                .send()
                .await
                .context("EC2 DescribeInstances failed")?;

            for reservation in response.reservations() {
                for instance in reservation.instances() {
                    if let Some(worker) = Self::instance_to_worker(instance)? {
                        workers.push(worker);
                    }
                }
            }

            next_token = response.next_token().map(ToOwned::to_owned);
            if next_token.is_none() {
                break;
            }
        }

        workers.sort_by(|a, b| a.instance_id.cmp(&b.instance_id));

        Ok(workers)
    }

    /// Discovers running workers that currently report themselves as READY.
    ///
    /// Checks the status file for each running worker and retains workers whose
    /// state is READY and whose worker_ready flag is true.
    ///
    /// # Arguments
    ///
    /// * `self` - Worker manager containing the EC2 client and worker status directory
    ///
    /// # Returns
    /// Vector of running workers currently reporting READY
    pub async fn discover_ready_workers(&self) -> Result<Vec<Worker>> {
        let running = self.discover_running_workers().await?;
        let mut ready = Vec::with_capacity(running.len());

        for worker in running {
            let status_path = self
                .workers_root
                .join(&worker.instance_id)
                .join("status.json");

            match Self::read_status(&status_path).await {
                Ok(status) => {
                    let is_ready =
                        status.worker_ready.unwrap_or(false)
                            && status.state.as_deref() == Some("READY");

                    if is_ready {
                        ready.push(worker);
                    }
                }

                Err(err) => {
                    log::debug!(
                        "Worker {} is running but has no usable status.json: {}",
                        worker.instance_id,
                        err
                    );
                }
            }
        }

        ready.sort_by(|a, b| a.instance_id.cmp(&b.instance_id));

        Ok(ready)
    }

    /// Discovers and returns the requested number of READY workers.
    ///
    /// Fails if fewer READY workers are available than requested.
    ///
    /// # Arguments
    ///
    /// * `self` - Worker manager containing the EC2 client and worker status directory
    /// * `requested` - Number of READY workers required
    ///
    /// # Returns
    /// Vector containing the requested number of READY workers
    pub async fn require_ready_workers(&self, requested: usize) -> Result<Vec<Worker>> {
        if requested == 0 {
            return Err(anyhow!(
                "--distributed-workers must be greater than zero"
            ));
        }

        let workers = self.discover_ready_workers().await?;

        if workers.len() < requested {
            let running = self.discover_running_workers().await?;

            return Err(anyhow!(
                "Insufficient READY MMseqs workers: requested {}, found {} READY, {} running tagged workers",
                requested,
                workers.len(),
                running.len()
            ));
        }

        Ok(workers.into_iter().take(requested).collect())
    }

    /// Converts an EC2 instance returned by worker discovery into a Worker.
    ///
    /// Instances without an instance ID are ignored. A running worker without a
    /// private IP address is treated as an error.
    ///
    /// # Arguments
    ///
    /// * `instance` - EC2 instance returned by DescribeInstances
    ///
    /// # Returns
    /// Worker representation of the EC2 instance, or None if no instance ID exists
    fn instance_to_worker(instance: &Instance) -> Result<Option<Worker>> {
        let Some(instance_id) = instance.instance_id() else {
            return Ok(None);
        };

        let Some(private_ip) = instance.private_ip_address() else {
            return Err(anyhow!(
                "Running worker {} has no private IP address",
                instance_id
            ));
        };

        let instance_type = instance
            .instance_type()
            .map(|value| value.as_str().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let availability_zone = instance
            .placement()
            .and_then(|placement| placement.availability_zone())
            .map(ToOwned::to_owned);

        Ok(Some(Worker {
            instance_id: instance_id.to_string(),
            private_ip: private_ip.to_string(),
            instance_type,
            availability_zone,
        }))
    }

    /// Reads and parses a worker status file.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the worker status.json file
    ///
    /// # Returns
    /// Parsed WorkerStatus
    async fn read_status(path: &Path) -> Result<WorkerStatus> {
        let bytes = fs::read(path)
            .await
            .with_context(|| format!("failed to read {}", path.display()))?;

        serde_json::from_slice(&bytes)
            .with_context(|| format!("invalid worker status JSON: {}", path.display()))
    }
}


/// Selects the requested number of workers from a discovered worker inventory.
///
/// Sorts workers by instance ID to ensure deterministic selection and fails
/// if fewer workers are available than requested.
///
/// # Arguments
///
/// * `workers` - Running workers returned by worker discovery
/// * `requested` - Number of distributed workers requested for the run
///
/// # Returns
/// Vector containing the selected workers
pub fn select_workers(
    mut workers: Vec<Worker>,
    requested: usize,
) -> Result<Vec<Worker>, PipelineError> {
    if requested == 0 {
        return Err(PipelineError::InvalidConfig(
            "--distributed-workers must be greater than zero".to_string(),
        ));
    }

    // Do not depend on EC2 response ordering.
    workers.sort_by(|a, b| a.instance_id.cmp(&b.instance_id));

    if workers.len() < requested {
        return Err(PipelineError::InvalidConfig(format!(
            "Insufficient distributed workers: requested {}, found {} running tagged workers",
            requested,
            workers.len()
        )));
    }

    Ok(workers.into_iter().take(requested).collect())
}


#[cfg(test)]
mod tests {
    use super::*;

    fn worker(instance_id: &str) -> Worker {
        Worker {
            instance_id: instance_id.to_string(),
            private_ip: "10.140.1.1".to_string(),
            instance_type: "r6id.32xlarge".to_string(),
            availability_zone: Some("us-west-2b".to_string()),
        }
    }

    #[test]
    fn select_workers_returns_requested_workers_in_deterministic_order() {
        let workers = vec![
            worker("i-004"),
            worker("i-002"),
            worker("i-001"),
            worker("i-003"),
        ];

        let selected = select_workers(workers, 3).unwrap();

        let ids: Vec<&str> = selected
            .iter()
            .map(|worker| worker.instance_id.as_str())
            .collect();

        assert_eq!(ids, vec!["i-001", "i-002", "i-003"]);
    }

    #[test]
    fn select_workers_fails_when_insufficient_workers_exist() {
        let workers = vec![
            worker("i-001"),
            worker("i-002"),
        ];

        let result = select_workers(workers, 4);

        assert!(matches!(
            result,
            Err(PipelineError::InvalidConfig(_))
        ));
    }

    #[test]
    fn select_workers_rejects_zero_workers() {
        let workers = vec![worker("i-001")];

        let result = select_workers(workers, 0);

        assert!(matches!(
            result,
            Err(PipelineError::InvalidConfig(_))
        ));
    }
}