use anyhow::anyhow;
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::Aligner;
use tokio::sync::mpsc::Receiver;
use crate::utils::fastx::SequenceRecord;


/// Align contigs using rust-bio
///
/// # Arguments
///
/// * `contigs_rx` - stream of SequenceRecord
/// * `reference` - vec of ascii reference
///
/// # Returns
/// Result wrapping a Vec of tuple (position, Alignment Operations)
pub async fn align_contigs_to_reference(
    mut contigs_rx: Receiver<SequenceRecord>,
    reference: &[u8],
) -> anyhow::Result<Vec<(u64, Vec<AlignmentOperation>)>> {
    let mut alignments = Vec::new();
    let mut aligner = Aligner::new(-5, -1, |a, b| if a == b { 2 } else { -2 });

    while let Some(record) = contigs_rx.recv().await {
        if let SequenceRecord::Fasta { seq, .. } = record {
            let alignment = aligner.global(&seq, reference);
            let aligned_length = (alignment.yend - alignment.ystart) as u64;
            alignments.push((aligned_length, alignment.operations));
        }
    }
    if alignments.is_empty() {
        return Err(anyhow!("No valid FASTA records for alignment"));
    }
    Ok(alignments)
}