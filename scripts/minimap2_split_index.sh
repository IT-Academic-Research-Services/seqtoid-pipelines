#!/usr/bin/env bash
# prepare_and_index_nt.sh
#
# Split the full NT FASTA once, then build persistent minimap2 shard indexes
# for either/both pipeline modes:
#   sr       : short-read NT alignment, k=14, w=8
#   map-ont  : long-read NT alignment,  k=15, w=10
#
# Usage:
#   ./prepare_and_index_nt.sh both       # build both index sets (default)
#   ./prepare_and_index_nt.sh sr         # build only short-read indexes
#   ./prepare_and_index_nt.sh map-ont    # build only long-read indexes
#
# Requirements:
#   - minimap2 in PATH (or set MM2_BIN below)
#   - enough space in /scratch
#   - input file: /scratch/refs/nt.fa

set -euo pipefail

# ==================== CONFIG ====================

INPUT_FA="/scratch/refs/nt.fa"
OUTPUT_ROOT="/scratch/refs/nt_split"
PARTS_DIR="$OUTPUT_ROOT/parts"
SR_INDEX_DIR="$OUTPUT_ROOT/sr"
ONT_INDEX_DIR="$OUTPUT_ROOT/map-ont"

CHUNK_SIZE_GB=30
THREADS=64
I_LIMIT="150G"

MM2_BIN="minimap2"
MODE="${1:-both}"

# ================================================

case "$MODE" in
    sr|map-ont|both) ;;
    *)
        echo "Usage: $0 [sr|map-ont|both]" >&2
        exit 2
        ;;
esac

echo "=== NT preparation & indexing pipeline ==="
echo "Start time     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Mode           : $MODE"
echo "Input FASTA    : $INPUT_FA   ($(du -sh "$INPUT_FA" 2>/dev/null || echo '?'))"
echo "Parts dir      : $PARTS_DIR"
echo "SR index dir   : $SR_INDEX_DIR"
echo "ONT index dir  : $ONT_INDEX_DIR"
echo "Chunk size     : ~${CHUNK_SIZE_GB} GB"
echo "Threads/job    : $THREADS"
echo "RAM limit      : -I $I_LIMIT"
echo "minimap2       : $(command -v "$MM2_BIN" 2>/dev/null || echo "$MM2_BIN (not found?)")"
echo "               : $($MM2_BIN --version 2>/dev/null || echo "version check failed")"
echo ""

mkdir -p "$PARTS_DIR"

# ----------------------------------------------------------------
# 1. Split once. Both index sets are built from these exact shards.
# ----------------------------------------------------------------

if compgen -G "$PARTS_DIR/nt.part_*.fa" >/dev/null; then
    echo "Parts already exist in $PARTS_DIR -> skipping split"
    echo "Found $(find "$PARTS_DIR" -maxdepth 1 -type f -name 'nt.part_*.fa' | wc -l) files"
else
    echo "Splitting $INPUT_FA into ~${CHUNK_SIZE_GB} GB chunks..."

    split -b "${CHUNK_SIZE_GB}G" \
        --numeric-suffixes=001 \
        --additional-suffix=.fa \
        "$INPUT_FA" "$PARTS_DIR/x"

    echo "Renaming chunks..."
    for file in "$PARTS_DIR"/x*.fa; do
        if [[ -f "$file" ]]; then
            num=$(basename "$file" | cut -c2-4)
            dest="$PARTS_DIR/nt.part_${num}.fa"
            mv -v "$file" "$dest"
        fi
    done

    echo "Split finished. Files created:"
    ls -lh "$PARTS_DIR"/nt.part_*.fa | head -n 8
fi

# ----------------------------------------------------------------
# 2. Build one index family from the shared FASTA shards.
# ----------------------------------------------------------------

build_indexes() {
    local family="$1"
    local index_dir="$2"
    local preset="$3"
    local k="$4"
    local w="$5"

    mkdir -p "$index_dir"

    echo ""
    echo "=== Building $family indexes ==="
    echo "Preset         : $preset"
    echo "Index params   : k=$k w=$w"
    echo "Output dir     : $index_dir"
    echo "----------------------------------------"

    for fasta in "$PARTS_DIR"/nt.part_*.fa; do
        [[ -f "$fasta" ]] || continue

        local base
        local mmi
        local log
        local status

        base=$(basename "$fasta" .fa)
        mmi="$index_dir/${base}.mmi"
        log="${mmi}.log.err"

        if [[ -s "$mmi" ]]; then
            echo "SKIP  $(basename "$fasta") -> $mmi already exists ($(du -h "$mmi" | cut -f1))"
            continue
        fi

        echo ""
        echo "Indexing  $fasta"
        echo "Output    : $mmi"
        echo "Started   : $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Command   : $MM2_BIN -d $mmi $fasta -cx $preset -k $k -w $w -I $I_LIMIT -t $THREADS"

        set +e
        "$MM2_BIN" -d "$mmi" "$fasta" \
            -cx "$preset" -k "$k" -w "$w" -I "$I_LIMIT" -t "$THREADS" \
            2> >(tee "$log")
        status=$?
        set -e

        if (( status == 0 )); then
            echo "Finished  : $(date '+%Y-%m-%d %H:%M:%S')"
            ls -lh "$mmi"
        else
            echo "FAILED    : exit code $status"
            echo "Error log tail:"
            tail -n 10 "$log"
            echo "Removing partial index..."
            rm -f "$mmi"
        fi

        echo "----------------------------------------"
    done
}

# Preserve the existing short-read index recipe exactly.
if [[ "$MODE" == "sr" || "$MODE" == "both" ]]; then
    build_indexes "short-read NT" "$SR_INDEX_DIR" "sr" 14 8
fi

# Long-read RunNTAlignment uses -cx map-ont; observed effective index
# parameters are k=15, w=10, HPC=false. Keep k/w explicit here.
if [[ "$MODE" == "map-ont" || "$MODE" == "both" ]]; then
    build_indexes "long-read NT" "$ONT_INDEX_DIR" "map-ont" 15 10
fi

# ----------------------------------------------------------------
# Summary
# ----------------------------------------------------------------

echo ""
echo "=== Pipeline finished ==="
echo "End time      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Total parts   : $(find "$PARTS_DIR" -maxdepth 1 -type f -name 'nt.part_*.fa' | wc -l)"

if [[ -d "$SR_INDEX_DIR" ]]; then
    echo "SR .mmi       : $(find "$SR_INDEX_DIR" -maxdepth 1 -type f -name '*.mmi' | wc -l)"
fi
if [[ -d "$ONT_INDEX_DIR" ]]; then
    echo "map-ont .mmi  : $(find "$ONT_INDEX_DIR" -maxdepth 1 -type f -name '*.mmi' | wc -l)"
fi

echo ""
echo "Index directories:"
echo "  short reads : $SR_INDEX_DIR"
echo "  long reads  : $ONT_INDEX_DIR"
