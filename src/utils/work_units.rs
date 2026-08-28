//! Defines the distributed NR work-unit contract and lifecycle state machine.

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Lifecycle state of a distributed NR work unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum WorkUnitState {
    Available,
    Claimed,
    Running,
    Done,
    Failed,
}

/// Metadata describing a completed distributed NR result.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct WorkUnitResult {
    /// Worker instance that produced the result.
    pub worker_id: String,

    /// Attempt number that produced the result.
    pub attempt: u32,

    /// Path to the durable m8 result.
    pub result_path: PathBuf,

    /// Size of the m8 result in bytes.
    pub result_bytes: u64,

    /// Number of m8 alignment rows produced.
    ///
    /// This may legitimately be zero if the query completed successfully
    /// but produced no alignments.
    pub result_rows: u64,

    /// Optional checksum of the completed result.
    pub checksum: Option<String>,
}

/// Metadata describing a failed distributed NR attempt.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct WorkUnitFailure {
    /// Worker instance on which the failed attempt ran.
    pub worker_id: String,

    /// Attempt number that failed.
    pub attempt: u32,

    /// Human-readable failure description.
    pub reason: String,

    /// Whether the work unit may be returned to AVAILABLE for retry.
    pub retryable: bool,
}

/// Error returned when an invalid work-unit state transition is attempted.
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum WorkUnitTransitionError {
    #[error("invalid work-unit state transition: {state:?} -> {requested}")]
    InvalidTransition {
        state: WorkUnitState,
        requested: &'static str,
    },

    #[error("work unit is owned by worker {expected}, not {actual}")]
    WrongWorker { expected: String, actual: String },

    #[error("work unit attempt mismatch: expected {expected}, got {actual}")]
    AttemptMismatch { expected: u32, actual: u32 },

    #[error("work unit has no failed retry available")]
    RetryNotAvailable,

    #[error("work unit attempt counter exhausted")]
    AttemptExhausted,
}

/// Identifies one deterministic unit of distributed NR work.
///
/// A work unit represents one complete FASTQ chunk. For paired-end input,
/// `expected_input_count` counts logical read pairs; for single-end input,
/// it counts logical reads.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct WorkUnit {
    /// Pipeline run identifier containing this work unit.
    pub run_id: String,

    /// Sample identifier.
    pub sample_id: String,

    /// Deterministic chunk identifier within the sample.
    pub chunk_id: u64,

    /// Whether the work unit contains paired-end input.
    pub paired_end: bool,

    /// Number of logical query units expected in this chunk.
    ///
    /// For paired-end input this is the number of read pairs.
    /// For single-end input this is the number of reads.
    pub expected_input_count: u64,

    /// Path to the chunk R1 FASTQ.
    pub r1_path: PathBuf,

    /// Path to the chunk R2 FASTQ, if paired-end.
    pub r2_path: Option<PathBuf>,

    /// Reference set/version required to execute the work unit.
    pub reference_version: String,

    /// Number of the current execution attempt.
    ///
    /// The first claim is attempt 1. A retry receives the next attempt number.
    pub attempt: u32,

    /// Current lifecycle state.
    pub state: WorkUnitState,

    /// Worker currently owning the claimed/running attempt.
    pub claimed_by: Option<String>,

    /// Successful result metadata, present only when state is DONE.
    pub result: Option<WorkUnitResult>,

    /// Failure metadata, present only when state is FAILED.
    pub failure: Option<WorkUnitFailure>,
}

impl WorkUnit {
    /// Creates a new work unit in the AVAILABLE state.
    ///
    /// # Arguments
    ///
    /// * `run_id` - Pipeline run identifier
    /// * `sample_id` - Sample identifier
    /// * `chunk_id` - Deterministic chunk identifier
    /// * `paired_end` - Whether the chunk contains paired-end reads
    /// * `expected_input_count` - Number of logical reads or read pairs
    /// * `reference_version` - Required MMseqs reference version
    /// * `r1_path` - Chunk R1 FASTQ path
    /// * `r2_path` - Optional chunk R2 FASTQ path
    ///
    /// # Returns
    /// Newly created AVAILABLE work unit
    pub fn new(
        run_id: impl Into<String>,
        sample_id: impl Into<String>,
        chunk_id: u64,
        paired_end: bool,
        expected_input_count: u64,
        reference_version: impl Into<String>,
        r1_path: PathBuf,
        r2_path: Option<PathBuf>,
    ) -> Self {
        Self {
            run_id: run_id.into(),
            sample_id: sample_id.into(),
            chunk_id,
            paired_end,
            expected_input_count,
            r1_path,
            r2_path,
            reference_version: reference_version.into(),
            attempt: 0,
            state: WorkUnitState::Available,
            claimed_by: None,
            result: None,
            failure: None,
        }
    }

    /// Returns the deterministic identity of this work unit.
    ///
    /// # Returns
    ///
    /// String containing run, sample, and chunk identifiers
    pub fn id(&self) -> String {
        format!(
            "{}/{}:{:08}",
            self.run_id, self.sample_id, self.chunk_id
        )
    }

    /// Claims an AVAILABLE work unit for a worker.
    ///
    /// Claiming increments the attempt number. The attempt number therefore
    /// identifies the specific execution instance of this work unit.
    ///
    /// # Arguments
    ///
    /// * `worker_id` - Worker instance claiming the work unit
    ///
    /// # Returns
    /// Current attempt number
    pub fn claim(&mut self, worker_id: impl Into<String>) -> Result<u32, WorkUnitTransitionError> {
        if self.state != WorkUnitState::Available {
            return Err(WorkUnitTransitionError::InvalidTransition {
                state: self.state,
                requested: "CLAIMED",
            });
        }

        let next_attempt = self
            .attempt
            .checked_add(1)
            .ok_or(WorkUnitTransitionError::AttemptExhausted)?;

        self.attempt = next_attempt;
        self.state = WorkUnitState::Claimed;
        self.claimed_by = Some(worker_id.into());
        self.result = None;
        self.failure = None;

        Ok(self.attempt)
    }

    /// Marks a claimed work unit as RUNNING.
    ///
    /// # Arguments
    ///
    /// * `worker_id` - Worker executing the claimed attempt
    /// * `attempt` - Attempt number returned by `claim`
    ///
    /// # Returns
    /// Result indicating whether the transition succeeded
    pub fn start(
        &mut self,
        worker_id: &str,
        attempt: u32,
    ) -> Result<(), WorkUnitTransitionError> {
        self.verify_owner(worker_id, attempt)?;

        if self.state != WorkUnitState::Claimed {
            return Err(WorkUnitTransitionError::InvalidTransition {
                state: self.state,
                requested: "RUNNING",
            });
        }

        self.state = WorkUnitState::Running;

        Ok(())
    }

    /// Marks a running work unit as DONE and records its durable result metadata.
    ///
    /// Completion is accepted only from the worker and attempt that currently
    /// own the work unit. This prevents a late result from an older retry from
    /// overwriting a newer attempt.
    ///
    /// # Arguments
    ///
    /// * `worker_id` - Worker completing the attempt
    /// * `attempt` - Attempt number being completed
    /// * `result` - Durable result metadata
    ///
    /// # Returns
    /// Result indicating whether the transition succeeded
    pub fn complete(
        &mut self,
        worker_id: &str,
        attempt: u32,
        result: WorkUnitResult,
    ) -> Result<(), WorkUnitTransitionError> {
        self.verify_owner(worker_id, attempt)?;

        if self.state != WorkUnitState::Running {
            return Err(WorkUnitTransitionError::InvalidTransition {
                state: self.state,
                requested: "DONE",
            });
        }

        if result.worker_id != worker_id {
            return Err(WorkUnitTransitionError::WrongWorker {
                expected: worker_id.to_string(),
                actual: result.worker_id,
            });
        }

        if result.attempt != attempt {
            return Err(WorkUnitTransitionError::AttemptMismatch {
                expected: attempt,
                actual: result.attempt,
            });
        }

        self.state = WorkUnitState::Done;
        self.result = Some(result);
        self.failure = None;

        Ok(())
    }

    /// Marks a claimed or running work unit as FAILED.
    ///
    /// # Arguments
    ///
    /// * `worker_id` - Worker that owned the attempt
    /// * `attempt` - Attempt number that failed
    /// * `reason` - Human-readable failure description
    /// * `retryable` - Whether this failure may later be retried
    ///
    /// # Returns
    /// Result indicating whether the transition succeeded
    pub fn fail(
        &mut self,
        worker_id: &str,
        attempt: u32,
        reason: impl Into<String>,
        retryable: bool,
    ) -> Result<(), WorkUnitTransitionError> {
        self.verify_owner(worker_id, attempt)?;

        match self.state {
            WorkUnitState::Claimed | WorkUnitState::Running => {}
            _ => {
                return Err(WorkUnitTransitionError::InvalidTransition {
                    state: self.state,
                    requested: "FAILED",
                });
            }
        }

        self.state = WorkUnitState::Failed;
        self.result = None;
        self.failure = Some(WorkUnitFailure {
            worker_id: worker_id.to_string(),
            attempt,
            reason: reason.into(),
            retryable,
        });

        Ok(())
    }

    /// Returns a failed retryable work unit to AVAILABLE.
    ///
    /// The current attempt number is retained as historical state. The next
    /// claim increments it before execution begins.
    ///
    /// # Returns
    /// Result indicating whether the retry transition succeeded
    pub fn retry(&mut self) -> Result<(), WorkUnitTransitionError> {
        if self.state != WorkUnitState::Failed {
            return Err(WorkUnitTransitionError::InvalidTransition {
                state: self.state,
                requested: "AVAILABLE",
            });
        }

        let failure = self
            .failure
            .as_ref()
            .ok_or(WorkUnitTransitionError::RetryNotAvailable)?;

        if !failure.retryable {
            return Err(WorkUnitTransitionError::RetryNotAvailable);
        }

        self.state = WorkUnitState::Available;
        self.claimed_by = None;
        self.result = None;

        Ok(())
    }

    /// Verifies that a worker owns the current attempt.
    ///
    /// # Arguments
    ///
    /// * `worker_id` - Worker claiming ownership
    /// * `attempt` - Attempt being validated
    ///
    /// # Returns
    /// Result indicating whether ownership matches
    fn verify_owner(
        &self,
        worker_id: &str,
        attempt: u32,
    ) -> Result<(), WorkUnitTransitionError> {
        let expected_worker = self
            .claimed_by
            .as_deref()
            .unwrap_or("unclaimed");

        if expected_worker != worker_id {
            return Err(WorkUnitTransitionError::WrongWorker {
                expected: expected_worker.to_string(),
                actual: worker_id.to_string(),
            });
        }

        if self.attempt != attempt {
            return Err(WorkUnitTransitionError::AttemptMismatch {
                expected: self.attempt,
                actual: attempt,
            });
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn work_unit() -> WorkUnit {
        WorkUnit::new(
            "run-001",
            "sample-001",
            17,
            true,
            123_456,
            "phase2-mmseqs-cpu-2026-08-28",
            PathBuf::from("/efs/runs/run-001/sample-001/chunks/chunk_000017_R1.fastq"),
            Some(PathBuf::from(
                "/efs/runs/run-001/sample-001/chunks/chunk_000017_R2.fastq",
            )),
        )
    }

    fn result(worker_id: &str, attempt: u32) -> WorkUnitResult {
        WorkUnitResult {
            worker_id: worker_id.to_string(),
            attempt,
            result_path: PathBuf::from(
                "/efs/runs/run-001/sample-001/results/chunk_000017_attempt_001.m8",
            ),
            result_bytes: 1024,
            result_rows: 0,
            checksum: Some("abc123".to_string()),
        }
    }

    #[test]
    fn new_work_unit_starts_available() {
        let unit = work_unit();

        assert_eq!(unit.state, WorkUnitState::Available);
        assert_eq!(unit.attempt, 0);
        assert!(unit.claimed_by.is_none());
        assert!(unit.result.is_none());
        assert!(unit.failure.is_none());
    }

    #[test]
    fn lifecycle_available_claimed_running_done() {
        let mut unit = work_unit();

        let attempt = unit.claim("i-worker-001").unwrap();

        assert_eq!(attempt, 1);
        assert_eq!(unit.state, WorkUnitState::Claimed);
        assert_eq!(unit.claimed_by.as_deref(), Some("i-worker-001"));

        unit.start("i-worker-001", attempt).unwrap();

        assert_eq!(unit.state, WorkUnitState::Running);

        unit.complete("i-worker-001", attempt, result("i-worker-001", attempt))
            .unwrap();

        assert_eq!(unit.state, WorkUnitState::Done);
        assert!(unit.result.is_some());
        assert!(unit.failure.is_none());
    }

    #[test]
    fn cannot_claim_non_available_work_unit() {
        let mut unit = work_unit();

        unit.claim("i-worker-001").unwrap();

        let result = unit.claim("i-worker-002");

        assert!(matches!(
            result,
            Err(WorkUnitTransitionError::InvalidTransition {
                state: WorkUnitState::Claimed,
                requested: "CLAIMED"
            })
        ));
    }

    #[test]
    fn wrong_worker_cannot_start_or_complete_attempt() {
        let mut unit = work_unit();

        let attempt = unit.claim("i-worker-001").unwrap();

        let start = unit.start("i-worker-002", attempt);
        assert!(matches!(
            start,
            Err(WorkUnitTransitionError::WrongWorker { .. })
        ));

        unit.start("i-worker-001", attempt).unwrap();

        let complete = unit.complete("i-worker-002", attempt, result("i-worker-002", attempt));
        assert!(matches!(
            complete,
            Err(WorkUnitTransitionError::WrongWorker { .. })
        ));

        assert_eq!(unit.state, WorkUnitState::Running);
    }

    #[test]
    fn stale_attempt_cannot_complete_after_retry() {
        let mut unit = work_unit();

        let first_attempt = unit.claim("i-worker-001").unwrap();
        unit.start("i-worker-001", first_attempt).unwrap();
        unit.fail(
            "i-worker-001",
            first_attempt,
            "worker terminated",
            true,
        )
            .unwrap();

        unit.retry().unwrap();

        let second_attempt = unit.claim("i-worker-002").unwrap();
        assert_eq!(second_attempt, 2);

        unit.start("i-worker-002", second_attempt).unwrap();

        let stale_completion =
            unit.complete("i-worker-001", first_attempt, result("i-worker-001", first_attempt));

        assert!(stale_completion.is_err());
        assert_eq!(unit.state, WorkUnitState::Running);
        assert_eq!(unit.attempt, 2);
        assert_eq!(unit.claimed_by.as_deref(), Some("i-worker-002"));
    }

    #[test]
    fn retry_increments_attempt_on_next_claim() {
        let mut unit = work_unit();

        let first_attempt = unit.claim("i-worker-001").unwrap();
        unit.start("i-worker-001", first_attempt).unwrap();
        unit.fail(
            "i-worker-001",
            first_attempt,
            "MMseqs exited non-zero",
            true,
        )
            .unwrap();

        assert_eq!(unit.state, WorkUnitState::Failed);
        assert_eq!(unit.attempt, 1);

        unit.retry().unwrap();

        assert_eq!(unit.state, WorkUnitState::Available);
        assert_eq!(unit.attempt, 1);
        assert!(unit.claimed_by.is_none());

        let second_attempt = unit.claim("i-worker-002").unwrap();

        assert_eq!(second_attempt, 2);
        assert_eq!(unit.attempt, 2);
        assert_eq!(unit.claimed_by.as_deref(), Some("i-worker-002"));
    }

    #[test]
    fn non_retryable_failure_cannot_be_requeued() {
        let mut unit = work_unit();

        let attempt = unit.claim("i-worker-001").unwrap();
        unit.start("i-worker-001", attempt).unwrap();
        unit.fail(
            "i-worker-001",
            attempt,
            "invalid FASTQ chunk",
            false,
        )
            .unwrap();

        assert_eq!(unit.state, WorkUnitState::Failed);

        let retry = unit.retry();

        assert!(matches!(
            retry,
            Err(WorkUnitTransitionError::RetryNotAvailable)
        ));
    }

    #[test]
    fn work_unit_id_is_deterministic() {
        let unit = work_unit();

        assert_eq!(unit.id(), "run-001/sample-001:00000017");
    }

    #[test]
    fn work_unit_serializes_state_and_required_metadata() {
        let mut unit = work_unit();

        let attempt = unit.claim("i-worker-001").unwrap();
        unit.start("i-worker-001", attempt).unwrap();

        let json = serde_json::to_value(&unit).unwrap();

        assert_eq!(json["run_id"], "run-001");
        assert_eq!(json["sample_id"], "sample-001");
        assert_eq!(json["chunk_id"], 17);
        assert_eq!(json["paired_end"], true);
        assert_eq!(json["expected_input_count"], 123456);
        assert_eq!(json["reference_version"], "phase2-mmseqs-cpu-2026-08-28");
        assert_eq!(json["attempt"], 1);
        assert_eq!(json["state"], "RUNNING");
        assert_eq!(json["claimed_by"], "i-worker-001");
    }
}