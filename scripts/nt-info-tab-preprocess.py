#!/usr/bin/env python3
"""
Build core_nt_info.tab from a local BLAST core_nt DB.

Format (one line per accession, unversioned):
  accession<TAB>title<TAB>length

Usage:
  python3 build_nt_info_tab.py \\
      --db /data/refs/core_nt/core_nt \\
      --out /data/refs/core_nt/core_nt_info.tab

Requires: BLAST+ blastdbcmd on PATH.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

VERSION_RE = re.compile(r"\.\d+$")


def strip_version(acc: str) -> str:
    return VERSION_RE.sub("", acc)


def parse_blastdbcmd_line(raw: str) -> tuple[str, str, int] | None:
    """
    blastdbcmd outfmt '%a|%t|%l' — title may contain '|' or tabs.
    Length is always the last | field; accession is the first.
    """
    raw = raw.rstrip("\n\r")
    if not raw:
        return None

    parts = raw.split("|")
    if len(parts) < 3:
        return None

    acc = strip_version(parts[0].strip())
    length_s = parts[-1].strip()
    title = "|".join(parts[1:-1]).strip()
    # Titles sometimes contain tabs/newlines from NCBI
    title = title.replace("\t", " ").replace("\r", " ").replace("\n", " ")

    if not acc or not length_s.isdigit():
        return None

    return acc, title, int(length_s)


def main() -> int:
    p = argparse.ArgumentParser(description="Build unversioned nt_info.tab from core_nt")
    p.add_argument("--db", required=True, help="BLAST DB prefix (e.g. .../core_nt)")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument(
        "--keep-first",
        action="store_true",
        help="If duplicate accession, keep first (default). Else last wins.",
    )
    args = p.parse_args()

    db = args.db
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "blastdbcmd",
        "-db", db,
        "-entry", "all",
        "-outfmt", "%a|%t|%l",
    ]

    # accession -> (title, length); first-wins unless --keep-first is false
    table: dict[str, tuple[str, int]] = {}
    n_in = n_bad = n_dup = 0

    print(f"[build_nt_info] running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    assert proc.stdout is not None

    for line in proc.stdout:
        n_in += 1
        if n_in % 5_000_000 == 0:
            print(f"[build_nt_info] read {n_in:,} lines, kept {len(table):,}", file=sys.stderr)

        parsed = parse_blastdbcmd_line(line)
        if parsed is None:
            n_bad += 1
            continue
        acc, title, length = parsed
        if acc in table:
            n_dup += 1
            if args.keep_first:
                continue
        table[acc] = (title, length)

    stderr = proc.stderr.read() if proc.stderr else ""
    rc = proc.wait()
    if rc != 0:
        print(f"[build_nt_info] blastdbcmd failed rc={rc}\n{stderr}", file=sys.stderr)
        return rc

    print(f"[build_nt_info] writing {len(table):,} rows → {out_path}", file=sys.stderr)
    with out_path.open("w", encoding="utf-8") as fh:
        for acc, (title, length) in table.items():
            fh.write(f"{acc}\t{title}\t{length}\n")

    print(
        f"[build_nt_info] done: in={n_in:,} bad={n_bad:,} dups={n_dup:,} out={len(table):,}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())