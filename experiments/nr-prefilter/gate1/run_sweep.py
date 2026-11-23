#!/usr/bin/env python3
"""Paired single-GPU scheduling benchmark. Standard library only."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('/scratch/nr-portable-v2'))
    p.add_argument('--subset', type=int, default=10000)
    p.add_argument('--fasta', type=Path)
    p.add_argument('--binary', type=Path, default=Path(__file__).resolve().parent / 'bench')
    p.add_argument('--out', type=Path, default=Path('/scratch/gate1-results'))
    p.add_argument('--queries', type=int, default=32)
    p.add_argument('--lengths', default='31')
    p.add_argument('--batches', default='1,8,32')
    p.add_argument('--tiles', default='1,16,64', help='Tile budgets in MiB')
    p.add_argument('--trials', type=int, default=3)
    p.add_argument('--device', type=int, default=0)
    p.add_argument('--timeout', type=int, default=1800)
    p.add_argument('--resident', action='store_true', help='Also attempt full-resident control; must fit VRAM')
    a = p.parse_args()
    if a.queries < 1 or a.trials < 1:
        p.error('queries and trials must be positive')
    fasta = (a.fasta or a.root / f'subset_{a.subset}' / 'targets.faa').resolve()
    binary = a.binary.resolve()
    if not fasta.is_file() or not binary.is_file():
        p.error(f'Missing input or executable: {fasta}, {binary}')
    lengths, batches, tiles = ([int(x) for x in s.split(',')] for s in (a.lengths, a.batches, a.tiles))
    if any(x < 1 for x in lengths + batches + tiles):
        p.error('sweep values must be positive')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = a.out.resolve() / stamp
    out.mkdir(parents=True, exist_ok=False)
    def sha(path):
        h = hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(8 << 20), b''):
                h.update(block)
        return h.hexdigest()
    meta = dict(vars(a))
    meta.update(fasta=str(fasta), fasta_sha256=sha(fasta), binary_sha256=sha(binary),
                python=sys.version, note='Serial event synchronization per scoring call; bounded best/checksum output; NR-derived queries.')
    (out / 'metadata.json').write_text(json.dumps(meta, default=str, indent=2) + '\n')
    for cmd, name in [(['nvidia-smi'], 'nvidia-smi.txt'), (['nvcc', '--version'], 'nvcc.txt')]:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            (out / name).write_text(r.stdout + r.stderr)
        except (OSError, subprocess.TimeoutExpired) as e:
            (out / name).write_text(str(e))
    serial = 0
    def run(schedule, length, batch, tile, validate=False, small=False):
        nonlocal serial
        serial += 1
        cmd = [str(binary), '--fasta', str(fasta), '--schedule', schedule,
               '--queries', str(2 if small else a.queries), '--qlen', str(length),
               '--batch', str(batch), '--tile-mib', str(tile), '--device', str(a.device),
               '--seed', '17', '--validate', str(int(validate))]
        if small:
            cmd += ['--max-records', '64']
        prefix = out / f'{serial:04d}'
        prefix.with_suffix('.command.json').write_text(json.dumps(cmd))
        with prefix.with_suffix('.stdout').open('w') as stdout, prefix.with_suffix('.stderr').open('w') as stderr:
            subprocess.run(cmd, stdout=stdout, stderr=stderr, check=True, timeout=a.timeout)
        result = json.loads(prefix.with_suffix('.stdout').read_text())
        result['run_id'] = serial
        return result
    # Exact CPU comparison of every query-target score on the small input.
    for length in lengths:
        expected = None
        for schedule in ('query', 'tile'):
            result = run(schedule, length, 2, min(tiles), validate=True, small=True)
            if expected is not None and result['results'] != expected:
                raise RuntimeError('Correctness schedules disagree')
            expected = result['results']
    rows = []
    for length in lengths:
        expected = None
        for tile in tiles:
            configs = [('query', 1)] + [('tile', b) for b in batches]
            if a.resident:
                configs += [('resident', a.queries)]
            for trial in range(a.trials):
                order = configs if trial % 2 == 0 else list(reversed(configs))
                for schedule, batch in order:
                    r = run(schedule, length, batch, tile)
                    if expected is not None and r['results'] != expected:
                        raise RuntimeError('Score summary mismatch across schedules/tiles/trials')
                    expected = r['results']
                    r['trial'] = trial
                    with (out / 'results.jsonl').open('a') as f:
                        f.write(json.dumps(r) + '\n')
                    rows.append({k: v for k, v in r.items() if k != 'results'})
                    with (out / 'summary.csv').open('w', newline='') as f:
                        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                        writer.writeheader()
                        writer.writerows(rows)
                    print(f"{schedule:8} L={length} B={batch} tile={tile} trial={trial} wall={r['wall_ms']:.2f} ms", flush=True)
    (out / 'COMPLETE').write_text('All measured summaries agree; small exact CPU checks passed.\n')
    print(f'Results: {out}')

if __name__ == '__main__':
    main()
