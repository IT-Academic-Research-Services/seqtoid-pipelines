use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::anyhow;
use log::{info, warn};
use rand::prelude::*;
use rand::prelude::SliceRandom;
use tokio::sync::mpsc;
use tokio::task::JoinHandle;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt;

use crate::config::defs::{ClusterInfo, PipelineError, RunConfig};
use crate::utils::fastx::SequenceRecord;
use crate::utils::streams::ParseOutput;


/// Represents a sequence record with an associated weight for sampling.
#[derive(Clone)]
struct WeightedSampleItem {
    /// Random key for sampling.
    key: f64,
    /// The sequence record.
    record: SequenceRecord,
    /// The weight (e.g., duplicate count) of this record.
    weight: u64,
}


#[allow(dead_code)]
/// Subsamples the input FASTQ stream using weights from deduplication clusters.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `uniques_stream`: receiver for unique FASTQ records
/// * `duplicate_clusters`: map of cluster information
/// * `seed`: random seed for reproducibility
/// * `max_subsample`: maximum number of reads to keep
///
/// # Returns
///
/// Result containing:
/// - stream of subsampled FASTQ records
/// - receiver for the count of subsampled reads
/// - handle for the subsampling task
pub async fn subsample_weighted(
    config: Arc<RunConfig>,
    uniques_stream: mpsc::Receiver<ParseOutput>,
    duplicate_clusters: Arc<HashMap<String, ClusterInfo>>,
    seed: u64,
    max_subsample: u64,
) -> anyhow::Result<
    (
        ReceiverStream<ParseOutput>,
        tokio::sync::oneshot::Receiver<u64>,
        JoinHandle<anyhow::Result<(), anyhow::Error>>,
    ),
    PipelineError,
> {
    if max_subsample == 0 {
        let (tx, rx) = mpsc::channel(1);
        drop(tx);
        let (count_tx, count_rx) = tokio::sync::oneshot::channel();
        let _ = count_tx.send(0);
        return Ok((
            ReceiverStream::new(rx),
            count_rx,
            tokio::spawn(async { Ok(()) }),
        ));
    }

    let mut rng = config.rng.clone();

    let mut reservoir: Vec<WeightedSampleItem> = Vec::with_capacity(max_subsample as usize);
    let mut stream = ReceiverStream::new(uniques_stream);

    while let Some(item) = stream.next().await {
        let record = match item {
            ParseOutput::Fastq(rec) => rec,
            _ => continue,
        };

        let rep_id = record.id().to_string();
        let weight = duplicate_clusters
            .get(&rep_id)
            .map_or(1u64, |c| c.size as u64);

        let key = rng.random::<f64>().powf(1.0 / weight as f64);

        if reservoir.len() < max_subsample as usize {
            reservoir.push(WeightedSampleItem {
                key,
                record,
                weight,
            });
        } else if let Some(min_item) = reservoir
            .iter()
            .min_by(|a, b| a.key.partial_cmp(&b.key).unwrap())
        {
            if key > min_item.key {
                // Replace the smallest key
                if let Some(idx) = reservoir.iter().position(|r| r.key == min_item.key) {
                    reservoir[idx] = WeightedSampleItem {
                        key,
                        record,
                        weight,
                    };
                }
            }
        }
    }

    // Sort by key for consistent output order (same as original)
    reservoir.sort_by(|a, b| a.key.partial_cmp(&b.key).unwrap_or(Ordering::Equal));

    let count = reservoir.len() as u64;
    let (count_tx, count_rx) = tokio::sync::oneshot::channel();
    count_tx
        .send(count)
        .map_err(|_| PipelineError::Other(anyhow!("Count send failed")))?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let send_task = tokio::spawn(async move {
        for item in reservoir {
            tx.send(ParseOutput::Fastq(item.record)).await?;
        }
        Ok(())
    });

    Ok((ReceiverStream::new(rx), count_rx, send_task))
}

/// Subsamples the input FASTQ stream uniformly.
///
/// # Arguments
///
/// * `config`: the run configuration
/// * `uniques_rx`: receiver for unique FASTQ records
/// * `max_subsample`: maximum number of reads to keep
/// * `paired`: whether the input is paired-end
///
/// # Returns
///
/// Result containing:
/// - stream of subsampled FASTQ records
/// - receiver for the count of subsampled reads
/// - handle for the subsampling task
pub async fn subsample_uniform(
    config: Arc<RunConfig>,
    uniques_rx: mpsc::Receiver<ParseOutput>,
    max_subsample: u64,
    paired: bool,
) -> anyhow::Result<(
    ReceiverStream<ParseOutput>,
    tokio::sync::oneshot::Receiver<u64>,
    JoinHandle<anyhow::Result<()>>,
)> {
    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let (count_tx, count_rx) = tokio::sync::oneshot::channel();

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(uniques_rx);
        let mut records = Vec::new();
        let mut total = 0u64;

        while let Some(item) = stream.next().await {
            match &item {
                ParseOutput::Fastq(_) => {
                    records.push(item);
                    total += 1;
                }
                other => {
                    return Err(anyhow!(
                        "subsample_uniform: unexpected non-FASTQ item: {:?}",
                        other
                    ));
                }
            }
        }

        // Safety check: paired streams must be even-length
        if paired && total % 2 != 0 {
            return Err(anyhow!("Odd number of FASTQ records in paired mode"));
        }

        let num_units = if paired { total / 2 } else { total };
        let subsample_units = (max_subsample).min(num_units);
        let subsample_size = subsample_units * if paired { 2 } else { 1 };

        if count_tx.send(total).is_err() {
            warn!("Subsample count send failed (receiver might have dropped)");
        }

        let mut rng = config.rng.clone();

        if paired {
            let pairs: Vec<(ParseOutput, ParseOutput)> = records
                .chunks_exact(2)
                .map(|chunk| (chunk[0].clone(), chunk[1].clone()))
                .collect();

            let mut indices: Vec<usize> = (0..pairs.len()).collect();
            indices.shuffle(&mut rng);

            let sampled_indices: HashSet<usize> =
                indices.into_iter().take(subsample_units as usize).collect();

            for (i, (r1, r2)) in pairs.into_iter().enumerate() {
                if sampled_indices.contains(&i) {
                    tx.send(r1).await.map_err(|_| anyhow!("Send R1 failed"))?;
                    tx.send(r2).await.map_err(|_| anyhow!("Send R2 failed"))?;
                }
            }
        } else {
            let mut indices: Vec<usize> = (0..records.len()).collect();
            indices.shuffle(&mut rng);

            let sampled_indices: HashSet<usize> =
                indices.into_iter().take(subsample_units as usize).collect();

            for (i, item) in records.into_iter().enumerate() {
                if sampled_indices.contains(&i) {
                    tx.send(item).await.map_err(|_| anyhow!("Send failed"))?;
                }
            }
        }

        info!(
            "Subsample distributor processed {} uniques, sampled {}",
            total, subsample_size
        );

        Ok(())
    });

    Ok((ReceiverStream::new(rx), count_rx, task))
}

/// Uniformly subsample single-end FASTQ records using reservoir sampling.
///
/// Given X input FASTQ records and a requested maximum sample size Y:
///
///   - if X <= Y, every input record is emitted
///   - if X > Y, exactly Y records are selected
///   - every input record has probability Y / X of being selected
///   - sampling is without replacement
///
/// Memory usage is O(Y), rather than O(X).
///
/// The final reservoir is sorted by original stream position before emission,
/// so selected reads are emitted in the same relative order they appeared in
/// the input stream.
async fn subsample_uniform_reservoir(
    config: Arc<RunConfig>,
    input_rx: mpsc::Receiver<ParseOutput>,
    max_subsample: u64,
) -> anyhow::Result<(
    ReceiverStream<ParseOutput>,
    tokio::sync::oneshot::Receiver<u64>,
    JoinHandle<anyhow::Result<()>>,
)> {
    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let (count_tx, count_rx) = tokio::sync::oneshot::channel();

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(input_rx);
        let mut rng = config.rng.clone();

        // Each reservoir entry retains its original input position so that
        // the final selected records can be restored to input order.
        let mut reservoir: Vec<(u64, ParseOutput)> =
            Vec::with_capacity(max_subsample as usize);

        let mut total = 0u64;

        while let Some(item) = stream.next().await {
            match item {
                ParseOutput::Fastq(_) => {}
                other => {
                    return Err(anyhow!(
                        "subsample_uniform_reservoir: unexpected non-FASTQ item: {:?}",
                        other
                    ));
                }
            }

            let input_index = total;
            total += 1;

            if reservoir.len() < max_subsample as usize {
                // The first Y records initially fill the reservoir.
                reservoir.push((input_index, item));
            } else if max_subsample > 0 {
                // We have now seen `total` records.
                //
                // Choose uniformly from [0, total).
                let j = rng.random_range(0..total);

                // If j falls inside the reservoir, replace that entry.
                //
                // Probability of accepting this record:
                //
                //     Y / total
                //
                // which is exactly what is required for uniform reservoir
                // sampling.
                if j < max_subsample {
                    reservoir[j as usize] = (input_index, item);
                }
            }
        }

        // Reservoir replacement destroys input ordering. Restore it.
        reservoir.sort_unstable_by_key(|(input_index, _)| *input_index);

        let sampled = reservoir.len() as u64;

        if count_tx.send(sampled).is_err() {
            warn!(
                "subsample_uniform_reservoir: sample count receiver dropped"
            );
        }

        for (_, item) in reservoir {
            tx.send(item)
                .await
                .map_err(|_| anyhow!(
                    "subsample_uniform_reservoir: output receiver dropped"
                ))?;
        }

        info!(
            "Uniform reservoir subsampling processed {} reads, sampled {}",
            total,
            sampled
        );

        Ok(())
    });

    Ok((
        ReceiverStream::new(rx),
        count_rx,
        task,
    ))
}