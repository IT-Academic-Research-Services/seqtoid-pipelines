# SeqToID Documentation

High-throughput metagenomics pipelines written in **Rust**.

- [Installation](installation.md)
- **Pipelines**
    - [Consensus Genome](pipelines/consensus_genome.md)
    - [Create DB](pipelines/create_db.md)
    - …
- [Architecture](architecture/streaming.md)

---

eqToID (short for "Sequence to ID") is a collection of Rust-based bioinformatics pipelines designed for high-throughput metagenomics processing. These pipelines emphasize efficiency, leveraging asynchronous streams, in-memory processing (e.g., via RAM disks like `/dev/shm`), FIFO pipes, and parallel execution to handle large datasets quickly. A core principle is data integrity: stream functions are designed to never silently drop data, with explicit error handling, stall detection, and verbose logging for correctness and debugging. The focus is on speed without sacrificing accuracy, making it suitable for cluster environments with high-core CPUs and ample RAM.

The pipelines integrate external tools (e.g., via subprocesses) for tasks like read filtering, alignment, and variant calling, while Rust handles orchestration, data streaming, and custom utilities. Currently, the primary workflow is for generating consensus genomes from FASTQ inputs, with support for Illumina and (emerging) ONT technologies.

SeqToID has been developed primarily by Dr. Matthew Jobin at UCSF (matt.jobin@ucsf.edu), with contributions from the IT Academic Research Services team.

## Features
- **Streaming Architecture**: Processes data in chunks via Tokio async streams and Rayon thread pools, minimizing disk I/O and enabling real-time monitoring (e.g., I/O utilization via `iostat`).
- **Data Integrity**: Explicit buffer management, stall thresholds, and error propagation ensure no silent data loss. For example, T-junctions and Y-junctions for stream splitting/merging include done signals to confirm completion.
- **Resource Awareness**: Dynamically computes buffer sizes based on available RAM, CPU cores, and input size. Supports RAM temp directories for fast scratch space.
- **Modular Pipelines**: Easily extensible with modules for configuration, utilities (e.g., FASTX parsing, file manipulation), and CLI parsing.
- **Monitoring**: Background tasks for I/O utilization warnings (e.g., NVMe saturation) to prevent stalls.
- **External Tool Integration**: Seamless wrapping of tools like `fastp`, `minimap2`, `samtools`, etc., with version checks and CLI generation.
- **Output**: Generates consensus FASTA, BAM alignments, VCF variants, statistics (e.g., depth, coverage, allele counts), and QUAST evaluations.

## Workflows
SeqToID supports multiple workflows via the `--module` CLI flag. Each workflow processes FASTQ inputs and outputs results in a structured directory.



## System Requirements
SeqToID is optimized for high-performance computing environments but runs on desktops/laptops for smaller datasets.

- **OS**: POSIX-compatible (e.g., Linux, macOS). Tested on Rocky Linux (cluster), Ubuntu (EC2), and macOS (M4 MacBook Pro). Windows is not supported due to reliance on named pipes (FIFO).
- **Hardware**:
    - **Recommended (Cluster Nodes)**: Dual-socket AMD EPYC (84 cores), 1.5 TB RAM, 3.8 TB NVMe scratch space. Optional GPU nodes (NVIDIA H100 or L40S) for future accelerations.
    - **Minimum (Desktop)**: 4+ cores, 16 GB RAM, SSD storage. For large inputs (>10 GB), scale RAM proportionally (e.g., 64 GB+).
- **Rust**: Version 1.70+ (via `rustup`).
- **Dependencies**: See Software Requirements below.

## Software Requirements
SeqToID requires the following external tools, which must be in your `$PATH`. Minimum versions are based on tested compatibility; newer versions are generally fine.

| Software   | Minimum Version | Purpose                          | Installation Notes (apt or Manual) |
|:-----------|:---------------:|:---------------------------------|:----------------------------------|
| fastp      | 0.23.0         | Read trimming and QC             | `sudo apt install fastp` or download from [GitHub](https://github.com/OpenGene/fastp) |
| minimap2   | 2.24           | Alignment to references          | `sudo apt install minimap2` or download from [GitHub](https://github.com/lh3/minimap2) |
| samtools   | 1.19           | BAM processing and stats         | `sudo apt install samtools` |
| kraken2    | 2.1.2          | Taxonomic classification         | Install via [GitHub](https://github.com/DerrickWood/kraken2) or Bioconda; requires DB download |
| bcftools   | 1.19           | Variant calling and consensus    | `sudo apt install bcftools` |
| mafft      | 7.490          | Multiple sequence alignment      | `sudo apt install mafft` |
| quast      | 5.2.0          | Assembly evaluation              | Install via [GitHub](https://github.com/ablab/quast) or Bioconda; requires Python 3.8+ |
| seqkit     | 2.3.0          | Sequence stats and manipulation  | `sudo apt install seqkit` or download from [GitHub](https://github.com/shenwei356/seqkit) |
| pigz       | 2.6            | Parallel gzip compression        | `sudo apt install pigz` |
| iostat     | (sysstat 12+)  | I/O monitoring (optional)        | `sudo apt install sysstat` |
| python3    | 3.8            | Required by QUAST                | `sudo apt install python3` |

For Kraken2, download a database (e.g., MiniKraken) and set `--kraken_db` path.


## Usage
Run with `seqtoid-pipelines [OPTIONS]`. Uses Clap for parsing.

**Key Options** (run `--help` for full list):
- `--module <MODULE>`: Workflow (e.g., `consensus_genome`, `create_db`). Required.
- `--file1 <PATH>`: Primary FASTQ input (R1 or single-end). Required.
- `--file2 <PATH>`: Secondary FASTQ (R2 for paired-end). Optional.
- `--technology <TECH>`: `illumina` or `ont`. Default: illumina.
- `--out_dir <DIR>`: Output directory. Default: auto-generated `<sample_base>_YYYYMMDD`.
- `--threads <N>`: Max threads (capped by CPU cores). Default: all cores.
- `--ref_db <PATH>`: HDF5 reference database.
- `--target_accession <ID>` / `--target_sequence <PATH>` / `--target_index <PATH>`: Target reference.
- `--host_accession <ID>` / `--host_sequence <PATH>` / `--host_index <PATH>`: Host reference.
- `--ercc_sequence <PATH>` / `--ercc_index <PATH>`: ERCC controls (Illumina only).
- `--kraken_db <PATH>`: Kraken2 database path.
- `--target_taxid <ID>`: Taxonomy ID for Kraken filtering.
- `--verbose`: Enable detailed logging.
- `--max_reads <N>` / `--min_read_len <N>` / `--max_read_len <N>`: Read filters.
- `--stall_threshold <MS>`: Stream stall detection (default: varies).

**Example**:
```
seqtoid-pipelines --module consensus_genome --file1 input_R1.fastq.gz --file2 input_R2.fastq.gz --technology illumina --target_accession NC_045512 --host_sequence human_genome.fasta --kraken_db /path/to/kraken_db --target_taxid 2697049 --threads 32 --out_dir results
```

This runs consensus genome for SARS-CoV-2 (taxid 2697049), removing human reads.
