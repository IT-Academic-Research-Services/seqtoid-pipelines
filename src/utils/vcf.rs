use anyhow::anyhow;
use tokio::sync::mpsc::Receiver;
use tokio_stream::wrappers::ReceiverStream;
use tokio_stream::StreamExt; // Import StreamExt for next()
use crate::utils::streams::ParseOutput;
use std::collections::HashSet;

/// Compute allel counts from a stream
///
/// # Arguments
///
/// * `rx` - A Receiver stream of  ParseMode::Fasta,
///
/// # Returns
///
/// result: tuple of (snps, mnps, indels) as (u64, u64, u64)
pub async fn parse_vcf_stream(mut stream: Receiver<ParseOutput>) -> anyhow::Result<(u64, u64, u64)> {
    let mut snps = 0;
    let mut mnps = 0;
    let mut indels = 0;
    let mut stream = ReceiverStream::new(stream); // Create stream once

    while let Some(item) = stream.next().await {
        match item {
            ParseOutput::Bytes(bytes) => {
                let lines = String::from_utf8_lossy(&bytes);
                for line in lines.lines() {
                    if !line.starts_with('#') {
                        let fields: Vec<&str> = line.split('\t').collect();
                        if fields.len() >= 5 {
                            let ref_allele = fields[3];
                            let alt_alleles = fields[4].split(',');
                            let allele_lens: HashSet<usize> = [ref_allele.len()]
                                .iter()
                                .copied()
                                .chain(alt_alleles.map(|a| a.len()))
                                .collect();
                            if allele_lens.len() > 1 {
                                indels += 1;
                            } else {
                                let l = allele_lens.into_iter().next().unwrap();
                                if l == 1 {
                                    snps += 1;
                                } else {
                                    mnps += 1;
                                }
                            }
                        } else {
                            return Err(anyhow!(
                                "Invalid VCF line format: expected at least 5 fields, found {}",
                                fields.len()
                            ));
                        }
                    }
                }
            }
            _ => return Err(anyhow!("Unexpected ParseOutput variant in VCF stream")),
        }
    }
    Ok((snps, mnps, indels))
}

