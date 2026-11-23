#!/usr/bin/env python3
"""Native-MMseqs NR preparation. Never decodes MMseqs sequence data in Python."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess

SIZES = (10000, 100000, 1000000, 10000000)

def write_json(p, obj):
    t = p.with_name(p.name + '.tmp')
    t.write_text(json.dumps(obj, indent=2) + '\n')
    os.replace(t, p)

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()

def command(args, log):
    print('RUN:', ' '.join(map(str, args)), flush=True)
    with log.open('w') as f:
        subprocess.run(list(map(str, args)), stdout=f, stderr=subprocess.STDOUT, check=True)

def fasta_stats(p):
    n = r = 0
    with p.open('rb') as f:
        for line in f:
            if line.startswith(b'>'):
                n += 1
            elif line.strip():
                if n == 0:
                    raise RuntimeError('Sequence before first FASTA header')
                r += len(line.strip())
    return n, r

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db', type=Path, default=Path('/scratch/refs/mmseqs/nrcleanDB'))
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--seed', type=int, default=20261120)
    ap.add_argument('--shard-mib', type=int, default=1024)
    a = ap.parse_args()
    if a.shard_mib < 1:
        ap.error('--shard-mib must be positive')
    db = a.db.resolve()
    idx = Path(str(db) + '.index')
    if not idx.is_file() or not Path(str(db) + '.dbtype').is_file():
        raise RuntimeError('Missing source index or dbtype')
    if not shutil.which('mmseqs'):
        raise RuntimeError('mmseqs not on PATH')
    out = a.out.resolve()
    if out == db.parent or db.parent in out.parents:
        raise RuntimeError('Output must be outside the source reference directory')
    out.mkdir(parents=True, exist_ok=True)
    with (out / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        run(a, db, idx, out)

def run(a, db, idx, out):
    version = subprocess.check_output(['mmseqs', 'version'], text=True).strip()
    cfg = {'portable_schema': 2, 'db': str(db), 'index_bytes': idx.stat().st_size,
           'index_mtime_ns': idx.stat().st_mtime_ns,
           'dbtype_hex': Path(str(db) + '.dbtype').read_bytes().hex(),
           'seed': a.seed, 'sizes': list(SIZES), 'shard_mib': a.shard_mib,
           'mmseqs_version': version, 'python_version': __import__('sys').version}
    cp = out / 'configuration.json'
    if cp.exists():
        if json.loads(cp.read_text()) != cfg:
            raise RuntimeError('Configuration changed; use a new output directory')
    else:
        if any(p.name != '.lock' for p in out.iterdir()):
            raise RuntimeError('Use an empty output directory')
        write_json(cp, cfg)
    print('Counting source index records...', flush=True)
    with idx.open('rb') as f:
        n = sum(1 for _ in f)
    if n < max(SIZES):
        raise RuntimeError('Source smaller than requested largest subset')
    keysdir = out / 'selection'
    if not (keysdir / 'complete.json').exists():
        if keysdir.exists():
            shutil.rmtree(keysdir)
        keysdir.mkdir()
        print('Selecting 10 million random row positions; this uses several GB RAM.', flush=True)
        ranks = random.Random(a.seed).sample(range(n), max(SIZES))
        tiers = {}
        start = 0
        for j, size in enumerate(SIZES):
            for row in ranks[start:size]:
                tiers[row] = j
            start = size
        del ranks
        files = [(keysdir / f'keys_{s}.txt').open('wb') for s in SIZES]
        counts = [0] * len(SIZES)
        try:
            with idx.open('rb') as f:
                for row, line in enumerate(f):
                    tier = tiers.get(row)
                    if tier is not None:
                        key = line.split()[0]
                        int(key)
                        for j in range(tier, len(SIZES)):
                            files[j].write(key + b'\n')
                            counts[j] += 1
        finally:
            for f in files:
                f.close()
        if counts != list(SIZES):
            raise RuntimeError('Key selection count mismatch')
        write_json(keysdir / 'complete.json', {'source_records': n, 'counts': counts,
                   'sha256': {f'keys_{s}.txt': digest(keysdir / f'keys_{s}.txt') for s in SIZES}})
        del tiers
    for size in SIZES:
        d = out / f'subset_{size}'
        if (d / 'complete.json').exists():
            print(f'Skipping completed subset {size}', flush=True)
            continue
        if d.exists():
            shutil.rmtree(d)
        d.mkdir()
        shutil.copyfile(keysdir / f'keys_{size}.txt', d / 'target_keys.txt')
        command(['mmseqs', 'createsubdb', d / 'target_keys.txt', db, d / 'targets',
                 '--subdb-mode', '0', '--id-mode', '0'], d / 'createsubdb.log')
        # Materialize ONLY the selected headers, never the full header DB.
        for path in list(d.iterdir()):
            if path.is_symlink():
                path.unlink()
        header_db = Path(str(db) + '_h')
        if not Path(str(header_db) + '.index').is_file():
            raise RuntimeError('Source header index is required')
        command(['mmseqs', 'createsubdb', d / 'target_keys.txt', header_db,
                 d / 'targets_h', '--subdb-mode', '0', '--id-mode', '0'],
                d / 'createheaders.log')
        # Header subsetting may itself add auxiliary links. Remove all of them.
        for path in list(d.iterdir()):
            if path.is_symlink():
                path.unlink()
        with (d / 'target_keys.txt').open('rb') as f:
            selected_keys = {int(line) for line in f}
        lookup = Path(str(db) + '.lookup')
        if not lookup.is_file():
            raise RuntimeError('Source lookup is required for portable identifier mapping')
        lookup_count = 0
        with lookup.open('rb') as src, (d / 'targets.lookup').open('wb') as dst:
            for line in src:
                if int(line.split(None, 1)[0]) in selected_keys:
                    dst.write(line)
                    lookup_count += 1
        if lookup_count != size:
            raise RuntimeError('Subset lookup count mismatch')
        del selected_keys
        source_file = Path(str(db) + '.source')
        if source_file.is_file():
            shutil.copyfile(source_file, d / 'targets.source')
        if any(path.is_symlink() for path in d.rglob('*')):
            raise RuntimeError('Subset contains a symlink; refusing completion')
        command(['mmseqs', 'convert2fasta', d / 'targets', d / 'targets.faa'], d / 'convert2fasta.log')
        count, residues = fasta_stats(d / 'targets.faa')
        with (d / 'targets.index').open('rb') as f:
            actual = [int(line.split()[0]) for line in f]
        with (d / 'target_keys.txt').open('rb') as f:
            expected = [int(line) for line in f]
        if count != size or sorted(actual) != sorted(expected):
            raise RuntimeError(f'Subset validation failed: {size}')
        del actual, expected
        write_json(d / 'complete.json', {'records': count, 'residues': residues,
                   'files': {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)}
                             for p in d.iterdir() if p.is_file()}})
    full = out / 'full'
    full.mkdir(exist_ok=True)
    exported = full / 'export_complete.json'
    fasta = full / 'nr.faa'
    if not (full / 'complete.json').exists():
        if not exported.exists():
            partial = full / 'nr.partial.faa'
            if partial.exists():
                partial.unlink()
            command(['mmseqs', 'convert2fasta', db, partial], full / 'convert2fasta.log')
            count, residues = fasta_stats(partial)
            if count != n:
                raise RuntimeError(f'Full export count {count} != source count {n}')
            os.replace(partial, fasta)
            write_json(exported, {'records': count, 'residues': residues,
                                  'bytes': fasta.stat().st_size, 'sha256': digest(fasta)})
        if not fasta.exists():
            raise RuntimeError('Export marker exists but nr.faa is missing')
        shards = full / 'shards'
        if shards.exists():
            shutil.rmtree(shards)
        shards.mkdir()
        manifests = []
        output = None
        total = 0
        limit = a.shard_mib << 20
        def close_shard():
            if output is not None:
                output.close()
                manifests.append({'file': name, 'bytes': written, 'records': records,
                                  'residues': residues, 'sha256': h.hexdigest()})
        try:
            with fasta.open('rb') as source:
                for line in source:
                    header = line.startswith(b'>')
                    if header and (output is None or written >= limit):
                        close_shard()
                        name = f'nr_{len(manifests):05d}.faa'
                        output = (shards / name).open('wb', buffering=8 << 20)
                        h = hashlib.sha256()
                        written = records = residues = 0
                    if output is None:
                        raise RuntimeError('Invalid FASTA export')
                    output.write(line)
                    h.update(line)
                    written += len(line)
                    if header:
                        records += 1
                        total += 1
                    else:
                        residues += len(line.strip())
            close_shard()
        finally:
            if output is not None and not output.closed:
                output.close()
        if total != n:
            raise RuntimeError('Shard count mismatch')
        if sum(x['residues'] for x in manifests) != json.loads(exported.read_text())['residues']:
            raise RuntimeError('Shard residue mismatch')
        write_json(full / 'complete.json', {'records': total, 'shards': manifests,
                   'headers': 'Preserved exactly as emitted by MMseqs; not assumed to be numeric keys.'})
    # Preserve full-reference identifier mapping as ordinary files.
    for suffix in ('.lookup', '.source'):
        source = Path(str(db) + suffix)
        dest = full / ('original_targets' + suffix)
        if source.is_file() and not dest.exists():
            shutil.copyfile(source, dest)
    links = [str(path) for path in out.rglob('*') if path.is_symlink()]
    if links:
        raise RuntimeError('Refusing completion: symlinks remain: ' + repr(links))
    write_json(out / 'portability.json', {
        'symlinks': 0,
        'mapping_files': {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)}
                          for p in full.glob('original_targets.*')},
        'note': 'Source paths in configuration/logs are provenance only, not runtime dependencies.'})
    write_json(out / 'complete.json', {'source_records': n, 'sizes': list(SIZES),
               'note': 'full/nr.faa retained in addition to shards; removable after validation.'})
    print(f'COMPLETE: {out}', flush=True)
    print('Full intermediate nr.faa retained; remove only after reviewing completion manifests.', flush=True)

if __name__ == '__main__':
    main()
