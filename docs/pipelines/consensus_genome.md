### Consensus Genome (`--module consensus_genome`)
This workflow generates a consensus genome from paired or single-end FASTQ reads, typically for viral or microbial targets in metagenomic samples.

**Steps**:
1. **Input Validation**: Reads and interleaves FASTQ files (R1/R2), filters by length, and counts records. Outputs validated interleaved FASTQ.
2. **Host Removal**: Aligns reads to a host reference (e.g., human genome) using `minimap2` and filters unmapped reads with `samtools`.
3. **ERCC Processing** (Illumina only): Optional alignment to ERCC spike-ins for quality control, with mapped read stats.
4. **Kraken Filtering**: Classifies reads with `kraken2` against a target taxonomy ID, retaining relevant reads.
5. **Target Alignment**: Aligns filtered reads to a target reference using `minimap2` and sorts/indexes BAM with `samtools`.
6. **Consensus Generation**: Calls variants with `bcftools`, generates consensus FASTA, and realigns to reference with `mafft`.
7. **Variant Calling**: Produces VCF with `bcftools mpileup` and stats.
8. **Statistics**: Computes read counts, mapping rates, depth distributions, allele frequencies, and coverage bins using `samtools stats/depth`, `seqkit`, and custom parsers.
9. **Assembly Evaluation**: Runs `quast` on consensus vs. reference, including alignment BAM.

**Inputs**:
- FASTQ files (`--file1`, optional `--file2` for paired-end).
- Reference sequences or accessions (e.g., `--target_accession`, `--host_sequence`).
- Optional: HDF5 database for reference retrieval, Kraken DB path.

**Outputs** (in `--out_dir` or auto-generated dir):
- Validated/no-host/filtered FASTQ files.
- BAM alignments, VCF variants, consensus FASTA.
- JSON stats (e.g., `stats.json` with depth, coverage, variants).
- QUAST reports for assembly quality.
- Depth plots (via custom plotting utils).

**Supported Technologies**: Illumina (fully implemented); ONT (in progress, with FASTQ parsing and basic streaming).

## Data Flow

```mermaid
flowchart TD
    A[Start: FASTQ Inputs R1/R2] --> B{Technology?}
    B -->|Illumina| C[Input Validation: Read, Interleave, Filter, Count, Write Validated FQ.gz]
    B -->|ONT| Z[ONT Pipeline - In Progress]
    C --> D[Host Removal: Minimap2 Align to Host, Samtools Filter Unmapped, Convert to FASTQ]
    D --> E[ERCC Processing: Minimap2 Align to ERCC, Samtools Stats for QC]
    E --> F[Kraken Filtering: Kraken2 Classify, Filter by Target TaxID, Output Filtered FASTQ]
    F --> G[Target Alignment: Minimap2 Align to Target, Samtools View Mapped, Sort, Index BAM]
    G --> H[Consensus Generation: Samtools Mpileup, Bcftools Call, Bcftools Consensus → FASTA]
    G --> I[Variant Calling: Bcftools Mpileup, Call → VCF, Stats]
    H --> J[Realign Consensus: Mafft Align Consensus to Ref]
    G --> K[Statistics: Samtools Stats/Depth, Seqkit Stats, Allele Counts, Coverage Bins, Plot Depths]
    H --> K
    I --> K
    H --> L[Assembly Evaluation: QUAST on Consensus vs. Ref, with BAM]
    L --> M[End: Outputs - Consensus FASTA, BAM, VCF, Stats JSON, QUAST Reports, Depth Plots]
    K --> M
    Z --> M
    style Z fill:#ffcccc,stroke:#ff0000
