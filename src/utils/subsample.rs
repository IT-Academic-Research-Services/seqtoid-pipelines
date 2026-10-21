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

        // Safety check: paired streams must be even-length.
        if paired && total % 2 != 0 {
            return Err(anyhow!(
                "Odd number of FASTQ records in paired mode"
            ));
        }

        // Sampling units are:
        //   SE: individual reads
        //   PE: read pairs
        let num_units = if paired { total / 2 } else { total };

        let subsample_units = max_subsample.min(num_units);

        // Number of physical FASTQ records that will be emitted.
        let subsample_size =
            subsample_units * if paired { 2 } else { 1 };

        // Report the number of sampled LOGICAL units.
        //
        // SE: sampled reads
        // PE: sampled pairs
        //
        // The caller can convert this to physical FASTQ records with:
        //
        //   count * if paired { 2 } else { 1 }
        if count_tx.send(subsample_units).is_err() {
            warn!(
                "Subsample count send failed \
                 (receiver might have dropped)"
            );
        }

        let mut rng = config.rng.clone();

        if paired {
            let pairs: Vec<(ParseOutput, ParseOutput)> = records
                .chunks_exact(2)
                .map(|chunk| (chunk[0].clone(), chunk[1].clone()))
                .collect();

            let mut indices: Vec<usize> =
                (0..pairs.len()).collect();

            indices.shuffle(&mut rng);

            let sampled_indices: HashSet<usize> = indices
                .into_iter()
                .take(subsample_units as usize)
                .collect();

            for (i, (r1, r2)) in pairs.into_iter().enumerate() {
                if sampled_indices.contains(&i) {
                    tx.send(r1)
                        .await
                        .map_err(|_| anyhow!("Send R1 failed"))?;

                    tx.send(r2)
                        .await
                        .map_err(|_| anyhow!("Send R2 failed"))?;
                }
            }
        } else {
            let mut indices: Vec<usize> =
                (0..records.len()).collect();

            indices.shuffle(&mut rng);

            let sampled_indices: HashSet<usize> = indices
                .into_iter()
                .take(subsample_units as usize)
                .collect();

            for (i, item) in records.into_iter().enumerate() {
                if sampled_indices.contains(&i) {
                    tx.send(item)
                        .await
                        .map_err(|_| anyhow!("Send failed"))?;
                }
            }
        }

        info!(
            "Subsample distributor processed {} FASTQ records \
             ({} logical units), sampled {} FASTQ records \
             ({} logical units)",
            total,
            num_units,
            subsample_size,
            subsample_units
        );

        Ok(())
    });

    Ok((
        ReceiverStream::new(rx),
        count_rx,
        task,
    ))
}

/// Uniformly subsample FASTQ records using reservoir sampling.
///
/// Sampling units:
///   - single-end: one FASTQ record
///   - paired-end: two consecutive FASTQ records (R1, R2)
///
/// Given X logical input units and requested maximum Y:
///
///   - if X <= Y, every unit is emitted
///   - if X > Y, exactly Y units are selected
///   - every input unit has probability Y / X of being selected
///   - sampling is without replacement
///   - paired reads are always sampled and emitted together
///
/// Memory usage is O(Y), rather than O(X).
///
/// Selected units are sorted by their original stream position before
/// emission, preserving their relative input order.
pub async fn subsample_uniform_reservoir(
    config: Arc<RunConfig>,
    input_rx: mpsc::Receiver<ParseOutput>,
    max_subsample: u64,
    paired: bool,
) -> anyhow::Result<(
    ReceiverStream<ParseOutput>,
    tokio::sync::oneshot::Receiver<u64>,
    JoinHandle<anyhow::Result<()>>,
)> {
    struct SampleUnit {
        input_index: u64,
        records: Vec<ParseOutput>,
    }

    let reservoir_capacity = usize::try_from(max_subsample)
        .map_err(|_| {
            anyhow!(
                "subsample_uniform_reservoir: max_subsample {} exceeds usize capacity",
                max_subsample
            )
        })?;

    let (tx, rx) = mpsc::channel(config.base_buffer_size);
    let (count_tx, count_rx) = tokio::sync::oneshot::channel();

    let task = tokio::spawn(async move {
        let mut stream = ReceiverStream::new(input_rx);
        let mut rng = config.rng.clone();

        let mut reservoir: Vec<SampleUnit> =
            Vec::with_capacity(reservoir_capacity);

        let mut total_records = 0u64;
        let mut total_units = 0u64;

        loop {
            // ---------------------------------------------------------
            // Construct one logical sampling unit.
            //
            // SE: one record
            // PE: two consecutive records
            // ---------------------------------------------------------

            let first = match stream.next().await {
                Some(item) => item,
                None => break,
            };

            let first_record = match first {
                ParseOutput::Fastq(record) => ParseOutput::Fastq(record),

                other => {
                    return Err(anyhow!(
                        "subsample_uniform_reservoir: \
                         unexpected non-FASTQ item: {:?}",
                        other
                    ));
                }
            };

            total_records += 1;

            let records = if paired {
                let second = stream.next().await.ok_or_else(|| {
                    anyhow!(
                        "subsample_uniform_reservoir: \
                         odd number of FASTQ records in paired mode"
                    )
                })?;

                let second_record = match second {
                    ParseOutput::Fastq(record) => {
                        ParseOutput::Fastq(record)
                    }

                    other => {
                        return Err(anyhow!(
                            "subsample_uniform_reservoir: \
                             unexpected non-FASTQ R2 item: {:?}",
                            other
                        ));
                    }
                };

                total_records += 1;

                vec![first_record, second_record]
            } else {
                vec![first_record]
            };

            // `input_index` is the logical-unit position, not the physical
            // FASTQ-record position.
            let input_index = total_units;

            total_units += 1;

            let unit = SampleUnit {
                input_index,
                records,
            };

            // ---------------------------------------------------------
            // Algorithm R reservoir sampling.
            //
            // First Y units fill the reservoir.
            //
            // For every subsequent unit i, choose uniformly from all
            // units seen so far. If the chosen position lies within the
            // Y reservoir slots, replace that slot.
            // ---------------------------------------------------------

            if reservoir.len() < reservoir_capacity {
                reservoir.push(unit);
            } else if reservoir_capacity > 0 {
                let j = rng.random_range(0..total_units);

                if j < max_subsample {
                    reservoir[j as usize] = unit;
                }
            }
        }

        // -------------------------------------------------------------
        // Reservoir replacement does not preserve stream ordering.
        // Restore the selected units to their original relative order.
        // -------------------------------------------------------------

        reservoir.sort_unstable_by_key(|unit| unit.input_index);

        let sampled_units = reservoir.len() as u64;

        // Preserve the same logical-unit counting convention as the
        // existing subsample_uniform:
        //
        //   SE -> number of sampled reads
        //   PE -> number of sampled pairs
        if count_tx.send(sampled_units).is_err() {
            warn!(
                "subsample_uniform_reservoir: \
                 sample count receiver dropped"
            );
        }

        // -------------------------------------------------------------
        // Emit selected units.
        //
        // In PE mode both records are emitted consecutively, so a pair
        // can never be split by the sampler.
        // -------------------------------------------------------------

        for unit in reservoir {
            for record in unit.records {
                tx.send(record)
                    .await
                    .map_err(|_| {
                        anyhow!(
                            "subsample_uniform_reservoir: \
                             output receiver dropped"
                        )
                    })?;
            }
        }

        let sampled_records =
            sampled_units * if paired { 2 } else { 1 };

        info!(
            "Uniform reservoir subsampling processed {} records \
             ({} logical units), sampled {} records \
             ({} logical units)",
            total_records,
            total_units,
            sampled_records,
            sampled_units
        );

        Ok(())
    });

    Ok((
        ReceiverStream::new(rx),
        count_rx,
        task,
    ))
}