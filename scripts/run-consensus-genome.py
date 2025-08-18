import os
import subprocess
import re
import argparse
from datetime import datetime

def find_fastq_pairs(directory, sample_order):
    """
    Find paired FASTQ files in the directory and match them to the sample order.
    Returns a list of tuples: (sample_name, fastq_r1, fastq_r2)
    """
    fastq_pairs = []
    fastq_files = [f for f in os.listdir(directory) if f.endswith('.fastq') or f.endswith('.fastq.gz')]

    for sample in sample_order:
        # Match R1 and R2 files for the sample
        r1_pattern = re.compile(rf'^{sample}.*_r?1\.fastq(\.gz)?$', re.IGNORECASE)
        r2_pattern = re.compile(rf'^{sample}.*_r?2\.fastq(\.gz)?$', re.IGNORECASE)

        r1_file = None
        r2_file = None

        for fastq in fastq_files:
            if r1_pattern.match(fastq):
                r1_file = fastq
            elif r2_pattern.match(fastq):
                r2_file = fastq

        if r1_file and r2_file:
            fastq_pairs.append((sample, r1_file, r2_file))
        else:
            print(f"Warning: Could not find both R1 and R2 files for sample {sample}")

    return fastq_pairs

def run_seqtoid(fastq_dir, sample_name, r1_file, r2_file, kraken_db, adapter_fasta, quality, ref_sequence, log_file, ercc_sequence, host_sequence, ref_taxid, out):
    """
    Run the seqtoid-pipelines command for a single sample and extract runtime from console output.
    """
    command = [
        'seqtoid-pipelines',
        '--module', 'consensus_genome',
        '-k', kraken_db,
        '--adapter-fasta', adapter_fasta,
        '--quality', str(quality),
        '-i', os.path.join(fastq_dir, r1_file),
        '-I', os.path.join(fastq_dir, r2_file),
        '--ref-sequence', ref_sequence,
        '--ercc-sequences', ercc_sequence,
        '--host-sequence', host_sequence,
        '--ref-taxid', ref_taxid,
        '--out', out,
    ]

    # Format the command for logging
    command_str = ' '.join(command)
    # print(f"Running command for {sample_name}: {command_str}")
    # print(command_str)
    # print()

    # Log start time
    start_time_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Run the command and capture output
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
        status = 'SUCCESS'

        # Extract runtime from console output
        runtime_match = re.search(r'Run complete: (\d+) milliseconds\.', stdout)
        runtime_ms = int(runtime_match.group(1)) if runtime_match else None
        runtime = runtime_ms / 1000.0 if runtime_ms is not None else None
    except subprocess.CalledProcessError as e:
        stdout = e.stdout
        stderr = e.stderr
        status = 'FAILED'
        runtime = None

    # Log the results
    with open(log_file, 'a') as f:
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Start Time: {start_time_str}\n")
        f.write(f"Command: {command_str}\n")
        f.write(f"Status: {status}\n")
        f.write(f"Runtime: {runtime:.2f} seconds\n" if runtime is not None else "Runtime: Not found in output\n")
        f.write(f"Stdout:\n{stdout}\n")
        f.write(f"Stderr:\n{stderr}\n")
        f.write("-" * 80 + "\n")

    return status, runtime


def main():
    parser = argparse.ArgumentParser(description="Run seqtoid-pipelines on paired FASTQ files in specified order")
    parser.add_argument('--fastq_dir', required=True, help="Directory containing FASTQ files")
    parser.add_argument('--sample_list', required=True, help="File containing sample names in desired order")
    parser.add_argument('--kraken_db', default='/home/ubuntu/refs/kraken_db', help="Path to Kraken database")
    parser.add_argument('--adapter_fasta', default='/home/ubuntu/refs/TruSeq3-PE.fa', help="Adapter FASTA file")
    parser.add_argument('--quality', default=1, type=int, help="Quality threshold")
    parser.add_argument('--ref_sequence', default='/home/ubuntu/refs/covid-wuhan-1.fa', help="Reference sequence FASTA file")
    parser.add_argument('--log_file', default='seqtoid_run.log', help="Log file to store run information")
    # parser.add_argument('--max_reads', default='5000000000', help="Log file to store run information")
    parser.add_argument('--ercc-sequences', default='/home/ubuntu/refs/ercc_sequences.fasta')
    parser.add_argument('--host-sequence', default='/home/ubuntu/refs/hg38.fa')
    parser.add_argument('--ref-taxid', default='2697049')
    parser.add_argument('--out', default='/home/ubuntu/data/seqtoid')

    args = parser.parse_args()

    # Read sample order from file
    with open(args.sample_list, 'r') as f:
        sample_order = [line.strip() for line in f if line.strip()]

    # Find paired FASTQ files
    fastq_pairs = find_fastq_pairs(args.fastq_dir, sample_order)

    # Run seqtoid-pipelines for each pair
    for sample, r1_file, r2_file in fastq_pairs:
        print(f"Processing sample: {sample}")
        out_dir = os.path.join(args.out, sample)
        status, runtime = run_seqtoid(
            args.fastq_dir,
            sample,
            r1_file,
            r2_file,
            args.kraken_db,
            args.adapter_fasta,
            args.quality,
            args.ref_sequence,
            args.log_file,
            args.ercc_sequences,
            args.host_sequence,
            args.ref_taxid,
            out_dir
        )
        print(f"Completed {sample}: Status={status}, Runtime={runtime:.2f} seconds" if runtime is not None else f"Completed {sample}: Status={status}, Runtime=Not found")

if __name__ == "__main__":
    main()