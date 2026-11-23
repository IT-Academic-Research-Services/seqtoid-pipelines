# Gate 1 scheduling prototype

Standalone CUDA C++17 + standard-library Python. No Cargo changes. Inputs default
to `/scratch/nr-portable-v2/subset_10000/targets.faa`. Uses FASTA, not MMseqs GPU DBs.

Build on L40S: `make`. Override NVCCFLAGS for another GPU architecture.
Run: `python3 run_sweep.py`.
Optional resident control: `python3 run_sweep.py --resident`.
Results default to unique timestamp directories under `/scratch/gate1-results`.
Do not commit reference inputs, binaries, or result directories.

## Workload
BLOSUM62 ungapped local scoring over every diagonal, 32-bit accumulation, one warp
per target. No masking/composition correction, gapped alignment, E-values, or
MMseqs equivalence. Unknown residues (including U/O/J) map to X; B/Z/X/* have
explicit BLOSUM62 entries. Lowercase becomes uppercase. Targets over 100,000 aa
are rejected. Headers are not used; best-hit IDs are zero-based FASTA ordinals,
NOT original MMseqs keys. NR-derived query segments are reproducible with seed 17
within the same build/runtime. These are not actual translated-read queries.

Reference sequences reside in ordinary host RAM. Each upload stages residues,
relative offsets and lengths through reusable pinned buffers. Tiles end between
proteins; one oversized protein can exceed the nominal tile budget. GPU memory
also includes a per-tile score buffer and query/result buffers. No Q-by-N matrix.

Query-major and tile-major use the same score and reduction kernels for exactly
the same pairs. Query-major uploads each tile for each query. Tile-major uploads
each tile for each query batch. Resident mode uploads the entire reference once;
it must fit VRAM. Query upload and allocation are outside measured wall time and
reported together as initialization. Reference staging, H2D, scoring, reduction,
final D2H and synchronization are inside wall time. Warm-up is outside timing.
Device-event times are reported separately; final D2H is host wall time. Wall time
also includes event handling/launch/synchronization overhead and is not expected
to equal the sum of event times. Every query/tile call is synchronized: this is a
simple serial baseline, NOT a throughput-optimized implementation. No overlap.

Reduction returns best score, lowest ordinal on ties, and a modulo-2^64 checksum
covering all target scores. It excludes realistic top-K candidate-output costs.

The runner first compares every GPU score to the CPU oracle for two queries
against the first 64 proteins for each query length, then compares summaries
across all measured configurations. If the first 64 proteins cannot supply a
requested query length it fails explicitly. It records input/binary SHA256,
commands, stdout/stderr, machine information, JSONL and CSV. Trial order reverses
on alternating trials. Each configuration is a fresh process with warm-up;
reference parsing/packing is reported separately, not timed as scan work.

Start with defaults. Increase subset size only after correctness passes.
Forcing small tiles on a small subset is forced streaming, not true out-of-core.
The program currently accepts one FASTA file, not a shard directory; full-NR
shard iteration is deliberately not implemented in this first draft.

This kernel is intentionally simple and may be compute-bound. A speedup here is
only evidence for this implementation, not a demonstrated speedup over MMseqs.
A negative result does not rule out reuse with a faster scoring primitive.

## Status
Python syntax checked in the authoring environment. CUDA compiler/GPU unavailable
there: CUDA compilation and execution MUST be validated on the GPU instance.
