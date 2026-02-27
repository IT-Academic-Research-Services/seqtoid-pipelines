#!/usr/bin/env bash
# prepare_and_index_nt.sh
#
# Full pipeline: split large nt.fa → index chunks serially
# Uses the byte-split method that worked + the no-quotes indexing command that finally succeeded
#
# Requirements:
#   - minimap2 in PATH (or set MM2_BIN below)
#   - enough space in /scratch (local NVMe)
#   - input file: /scratch/refs/nt.fa   (843 GB in your case)
#
# Run with:
#   ./prepare_and_index_nt.sh

set -euo pipefail

# ==================== CONFIG ====================

INPUT_FA="/scratch/refs/nt.fa"                  # original full NT FASTA
OUTPUT_DIR="/scratch/refs/nt_split"

CHUNK_SIZE_GB=30                                # ~30 GB chunks → 28–29 parts
THREADS=64                                      # threads per minimap2 job (change to 96/128 later)
I_LIMIT="150G"                                  # -I value (increase to 180–220G after testing)

MM2_BIN="minimap2"                              # or full path if needed: "/usr/local/bin/minimap2"

# ================================================

echo "=== NT preparation & indexing pipeline ==="
echo "Start time     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Input FASTA    : $INPUT_FA   ($(du -sh "$INPUT_FA" 2>/dev/null || echo '?'))"
echo "Output dir     : $OUTPUT_DIR"
echo "Chunk size     : ~${CHUNK_SIZE_GB} GB"
echo "Threads/job    : $THREADS"
echo "RAM limit      : -I $I_LIMIT"
echo "minimap2       : $(which "$MM2_BIN" 2>/dev/null || echo "$MM2_BIN (not found?)")"
echo "               : $($MM2_BIN --version 2>/dev/null || echo "version check failed")"
echo ""

mkdir -p "$OUTPUT_DIR" || { echo "Cannot create $OUTPUT_DIR"; exit 1; }

# ────────────────────────────────────────────────
# 1. Splitting (byte-based, the method that worked)
# ────────────────────────────────────────────────

if ls "$OUTPUT_DIR"/nt.part_*.fa >/dev/null 2>&1; then
    echo "Parts already exist in $OUTPUT_DIR → skipping split"
    echo "Found $(ls "$OUTPUT_DIR"/nt.part_*.fa | wc -l) files"
else
    echo "Splitting $INPUT_FA into ~${CHUNK_SIZE_GB} GB chunks..."

    split -b "${CHUNK_SIZE_GB}G" \
        --numeric-suffixes=001 \
        --additional-suffix=.fa \
        "$INPUT_FA" "$OUTPUT_DIR/x"

    # Rename x001.fa → nt.part_001.fa  etc.
    echo "Renaming chunks..."
    for file in "$OUTPUT_DIR"/x*.fa; do
        if [[ -f "$file" ]]; then
            num=$(basename "$file" | cut -c2-4)   # 001, 002, ...
            dest="$OUTPUT_DIR/nt.part_${num}.fa"
            mv -v "$file" "$dest"
        fi
    done

    echo "Split finished. Files created:"
    ls -lh "$OUTPUT_DIR"/nt.part_*.fa | head -n 8
fi

# ────────────────────────────────────────────────
# 2. Serial indexing
# ────────────────────────────────────────────────

echo ""
echo "Starting serial indexing..."
echo "----------------------------------------"

cd "$OUTPUT_DIR" || { echo "Cannot cd to $OUTPUT_DIR"; exit 1; }

for fasta in nt.part_*.fa; do
    if [[ ! -f "$fasta" ]]; then continue; fi

    mmi="${fasta}.mmi"

    if [[ -s "$mmi" ]]; then
        echo "SKIP  $fasta  →  $mmi already exists ($(du -h "$mmi" | cut -f1))"
        continue
    fi

    echo ""
    echo "Indexing  $fasta"
    echo "Started   : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Command   : $MM2_BIN -d $mmi $fasta -cx sr -k 14 -w 8 -I $I_LIMIT -t $THREADS"

    # The working pattern — no quotes around filenames
    set +e   # let us capture exit code
    "$MM2_BIN" -d "$mmi" "$fasta" -cx sr -k 14 -w 8 -I "$I_LIMIT" -t "$THREADS" \
        2> >(tee "${mmi}.log.err")
    status=$?
    set -e

    if (( status == 0 )); then
        echo "Finished  : $(date '+%Y-%m-%d %H:%M:%S')"
        ls -lh "$mmi"
    else
        echo "FAILED    : exit code $status"
        echo "Error log tail:"
        tail -n 10 "${mmi}.log.err"
        echo "Removing partial index..."
        rm -f "$mmi"
        # Optional: exit 1   # stop on first error
    fi

    echo "----------------------------------------"
done

# ────────────────────────────────────────────────
# Summary
# ────────────────────────────────────────────────

echo ""
echo "=== Pipeline finished ==="
echo "End time      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Total .mmi    : $(ls -1 *.mmi 2>/dev/null | wc -l)"
echo "Total parts   : $(ls -1 nt.part_*.fa | wc -l)"
echo ""
echo "Quick checks:"
echo "  ls -lh *.mmi | head -n 6"
echo "  grep -i -l -e error -e fail *.log.err 2>/dev/null || echo 'No obvious error logs'"
echo ""
echo "Next steps:"
echo "  • If stable → try THREADS=96 or 128"
echo "  • If comfortable → add parallel=2–4 jobs (GNU parallel)"
echo "  • Then move to Rust-based alignment pipeline"
