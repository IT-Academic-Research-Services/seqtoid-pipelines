#!/usr/bin/env python3
"""
Benchmark the current SeqToID non-host NR alignment commands.

Experimental unit
-----------------
Each trial is uniquely identified by ALL of:
    sample_id
    backend
    instance_type
    instance_id
    host_label
    trial_id

This is intentional: e.g.
    diamond + r6id.32xlarge
is a different experimental condition from
    mmseqs-cpu + r6id.32xlarge
and from another r6id.32xlarge instance.

Inputs
------
Scan only the top level of /efs/runs for directories containing:
    nonhost_R1.fastq
    nonhost_R2.fastq   (optional)

A run directory such as:
    ERR11417004_20260930-123553
becomes sample ID:
    ERR11417004

Backends
--------
    diamond
    mmseqs-cpu
    mmseqs-gpu

Reference paths are hard-coded from the SeqToID phase2 reference-preparation
scripts. The benchmark assumes the references have already been built and validated
using the canonical CPU/GPU DB preparation workflows; it does not rebuild them:
    Diamond:    /scratch/refs/diamond/diamond_07_22_2026.dmnd
    MMseqs CPU: /scratch/refs/mmseqs/nrcleanDB
    MMseqs GPU: /scratch/refs/mmseqs-gpu/nrcleanDB_gpu

Command behavior
----------------
Diamond follows the current Rust non-host implementation:
    diamond blastx R1
    diamond blastx R2, if present
    copy the R1 m8 to diamond_merged_nr.m8 (matching current Rust behavior)

MMseqs follows the current Rust implementation:
    write_combined_fastq (R1 + optional R2)
    mmseqs createdb ... --dbtype 2
    GPU only: mmseqs gpuserver + 20 second warm-up
    mmseqs search ...
    GPU only: stop gpuserver immediately after search
    GPU search uses --gpu 1 --gpu-server 1 --db-load-mode 2 --prefilter-mode 1
    mmseqs convertalis ... --format-output ...

No extractorfs step is included because it is not in the current Rust
mmseqs_fastq_to_m8_file() implementation.

Result isolation
----------------
Every trial gets its own immutable directory:
    <results-root>/trials/<backend>/<instance-type>/<instance-id>/<trial-id>/<sample>/

Each condition also gets its own directory:
    <results-root>/conditions/<backend>/<instance-type>/<instance-id>/

That condition directory contains:
    machine_metadata.json
    benchmark_results.csv
    benchmark_stages.csv

There is also a global combined results area with a lock-protected CSV, but
its correctness is not needed to identify trials: every individual trial
manifest is self-contained and contains the complete backend + machine
identity.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Sequence


RUNS_ROOT = Path("/efs/runs")
DEFAULT_RESULTS_ROOT = RUNS_ROOT / "nr_alignment_benchmarks"
DEFAULT_SCRATCH = Path("/scratch")

DIAMOND_DB = Path("/scratch/refs/diamond/diamond_07_22_2026.dmnd")
MMSEQS_CPU_DB = Path("/scratch/refs/mmseqs/nrcleanDB")
MMSEQS_GPU_DB = Path("/scratch/refs/mmseqs-gpu/nrcleanDB_gpu")

TIMESTAMP_RE = re.compile(r"^(?P<sample>.+)_(?P<stamp>\d{8}-\d{6})$")

M8_FORMAT = (
    "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,"
    "tstart,tend,evalue,bits"
)

VALID_BACKENDS = {"diamond", "mmseqs-cpu", "mmseqs-gpu"}


@dataclass(frozen=True)
class RunInput:
    run_dir: Path
    sample_id: str
    r1: Path
    r2: Optional[Path]


@dataclass(frozen=True)
class CommandResult:
    stage: str
    command: list[str]
    started_at: str
    elapsed_seconds: float
    returncode: int
    stdout_bytes: int
    stderr_bytes: int
    log_path: str


class BenchmarkError(RuntimeError):
    pass


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def safe_component(value: object) -> str:
    text = str(value)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._")
    return text or "unknown"


def quote_arg(arg: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_./:=+-]+", arg):
        return arg
    return "'" + arg.replace("'", "'\\''") + "'"


def require_command(name: str) -> None:
    if shutil.which(name) is None:
        raise BenchmarkError(f"Required executable not found on PATH: {name}")


def run_capture(
        args: Sequence[str],
        *,
        stage: str,
        log_path: Path,
        cwd: Optional[Path] = None,
        env: Optional[dict[str, str]] = None,
        timeout: Optional[float] = None,
        stdout_path: Optional[Path] = None,
) -> CommandResult:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started_at = utc_now_iso()
    t0 = time.perf_counter()

    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    stdout_fh = None
    proc: Optional[subprocess.CompletedProcess[str]] = None
    try:
        if stdout_path is not None:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_fh = stdout_path.open("wb")

        with log_path.open("w", encoding="utf-8") as log:
            log.write("COMMAND: " + " ".join(quote_arg(x) for x in args) + "\n")
            log.write("STARTED_AT: " + started_at + "\n\n")
            try:
                proc = subprocess.run(
                    list(args),
                    cwd=str(cwd) if cwd else None,
                    env=merged_env,
                    stdout=stdout_fh if stdout_fh is not None else log,
                    stderr=log,
                    check=False,
                    timeout=timeout,
                    text=False if stdout_fh is not None else True,
                )
            except subprocess.TimeoutExpired as exc:
                elapsed = time.perf_counter() - t0
                log.write(f"\nTIMEOUT after {elapsed:.3f} seconds\n")
                raise BenchmarkError(f"{stage} timed out after {elapsed:.1f}s") from exc
    finally:
        if stdout_fh is not None:
            stdout_fh.close()

    assert proc is not None
    elapsed = time.perf_counter() - t0
    stdout_bytes = (
        stdout_path.stat().st_size
        if stdout_path is not None and stdout_path.exists()
        else log_path.stat().st_size
    )
    stderr_bytes = log_path.stat().st_size

    return CommandResult(
        stage=stage,
        command=list(args),
        started_at=started_at,
        elapsed_seconds=elapsed,
        returncode=proc.returncode,
        stdout_bytes=stdout_bytes,
        stderr_bytes=stderr_bytes,
        log_path=str(log_path),
    )


def strip_run_timestamp(name: str) -> str:
    match = TIMESTAMP_RE.match(name)
    return match.group("sample") if match else name


def discover_inputs(runs_root: Path) -> list[RunInput]:
    if not runs_root.is_dir():
        raise BenchmarkError(f"Runs root does not exist or is not a directory: {runs_root}")

    discovered: list[RunInput] = []
    for entry in sorted(runs_root.iterdir(), key=lambda p: p.name):
        if not entry.is_dir():
            continue
        r1 = entry / "nonhost_R1.fastq"
        r2 = entry / "nonhost_R2.fastq"
        if not r1.is_file():
            continue
        discovered.append(
            RunInput(
                run_dir=entry,
                sample_id=strip_run_timestamp(entry.name),
                r1=r1,
                r2=r2 if r2.is_file() else None,
            )
        )
    return discovered


def read_meminfo_bytes() -> int:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError):
        pass
    return 0


def available_ram_bytes() -> int:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError):
        pass
    return read_meminfo_bytes()


def filesystem_free_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


def imds_get(path: str) -> str:
    token_cmd = [
        "curl", "-sS", "-X", "PUT",
        "-H", "X-aws-ec2-metadata-token-ttl-seconds: 60",
        "http://169.254.169.254/latest/api/token",
    ]
    try:
        token = subprocess.run(token_cmd, capture_output=True, text=True, timeout=2, check=False)
        if token.returncode != 0 or not token.stdout.strip():
            return "unknown"
        out = subprocess.run(
            [
                "curl", "-sS",
                "-H", f"X-aws-ec2-metadata-token: {token.stdout.strip()}",
                f"http://169.254.169.254/latest/meta-data/{path}",
            ],
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
        return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def detect_instance_type() -> str:
    return imds_get("instance-type")


def detect_instance_id() -> str:
    return imds_get("instance-id")


def command_version(command: Sequence[str]) -> str:
    try:
        proc = subprocess.run(command, capture_output=True, text=True, timeout=20, check=False)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    combined = (proc.stdout + "\n" + proc.stderr).strip()
    lines = [x.strip() for x in combined.splitlines() if x.strip()]
    return lines[0] if lines else "unknown"


def cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name") and ":" in line:
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def count_fastq_records(path: Path) -> int:
    lines = 0
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            lines += chunk.count(b"\n")
    if lines % 4 != 0:
        raise BenchmarkError(f"FASTQ has non-multiple-of-4 line count: {path} ({lines})")
    return lines // 4


def file_rows(path: Path) -> int:
    try:
        with path.open("rb") as fh:
            return sum(
                chunk.count(b"\n")
                for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b"")
            )
    except OSError:
        return -1


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        )
    os.replace(tmp, path)


def append_csv_locked(path: Path, fieldnames: Sequence[str], row: dict[str, object]) -> None:
    """Append one row to a shared EFS CSV with a sibling-file advisory lock."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = path.with_name(path.name + ".lock")
    with lock_path.open("a", encoding="utf-8") as lock_fh:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
        try:
            exists = path.is_file() and path.stat().st_size > 0
            with path.open("a", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(fieldnames))
                if not exists:
                    writer.writeheader()
                writer.writerow(row)
                fh.flush()
                os.fsync(fh.fileno())
        finally:
            fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)


def diamond_dbinfo(db_path: Path) -> tuple[int, int]:
    proc = subprocess.run(
        ["diamond", "dbinfo", "--db", str(db_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise BenchmarkError(f"diamond dbinfo failed for {db_path}: {proc.stderr.strip()}")
    seq_match = re.search(r"Sequences\s+(\d+)", proc.stdout)
    letters_match = re.search(r"Letters\s+(\d+)", proc.stdout)
    if not seq_match or not letters_match:
        raise BenchmarkError(f"Could not parse diamond dbinfo output for {db_path}")
    return int(seq_match.group(1)), int(letters_match.group(1))


def compute_diamond_block_size(
        db_path: Path,
        scratch: Path,
        available_ram: int,
) -> tuple[float, int, float]:
    _, letters = diamond_dbinfo(db_path)
    available_ram_gb = available_ram / 1_073_741_824.0
    total_letters_billions = letters / 1e9

    if available_ram_gb >= 1400.0:
        ram_factor = 9.0
    elif available_ram_gb >= 900.0:
        ram_factor = 10.0
    elif available_ram_gb >= 128.0:
        ram_factor = 12.0
    else:
        ram_factor = 15.0

    block_size = available_ram_gb / ram_factor
    block_size = min(block_size, total_letters_billions)
    block_size = max(block_size, 6.0)
    block_size = min(block_size, 200.0)

    scratch_free_gib = filesystem_free_bytes(scratch) / 1_073_741_824.0
    estimated_scratch_gib = block_size * 3.0 + 20.0
    if scratch_free_gib < estimated_scratch_gib * 1.2:
        block_size *= 0.5
    elif scratch_free_gib < estimated_scratch_gib * 1.5:
        block_size *= 0.7

    return block_size, letters, scratch_free_gib


def diamond_index_chunks(threads: int, available_ram: int) -> int:
    if threads >= 192:
        chunks = 32
    elif threads >= 128:
        chunks = 24
    elif threads >= 96:
        chunks = 16
    elif threads >= 64:
        chunks = 12
    else:
        chunks = 4
    return min(chunks, int(available_ram / 12_000_000_000))


def copy_fastq(r1: Path, r2: Optional[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as out:
        with r1.open("rb") as fh:
            shutil.copyfileobj(fh, out, length=8 * 1024 * 1024)
        if r2 is not None:
            with r2.open("rb") as fh:
                shutil.copyfileobj(fh, out, length=8 * 1024 * 1024)


def ensure_success(result: CommandResult) -> None:
    if result.returncode != 0:
        raise BenchmarkError(
            f"Stage {result.stage} failed with exit code {result.returncode}; see {result.log_path}"
        )


def run_diamond(
        ri: RunInput,
        scratch: Path,
        work_dir: Path,
        threads: int,
) -> list[CommandResult]:
    db = DIAMOND_DB
    if not db.is_file():
        raise BenchmarkError(f"Diamond DB does not exist: {db}")

    temp_dir = work_dir / "diamond_tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    avail_ram = available_ram_bytes()
    block_size, letters, scratch_free_gib = compute_diamond_block_size(db, scratch, avail_ram)
    chunks = diamond_index_chunks(threads, avail_ram)
    if chunks < 1:
        raise BenchmarkError("Diamond index chunk calculation produced zero chunks")

    db_prefix = db.with_suffix("")
    common = [
        "diamond", "blastx",
        "-d", str(db_prefix),
        "--threads", str(threads),
        "--mid-sensitive",
        "--block-size", f"{block_size:.1f}",
        "-c", str(chunks),
        "-f", "6",
        "--tmpdir", str(temp_dir),
        "--unal", "0",
    ]

    results: list[CommandResult] = []

    r1_out = temp_dir / "diamond_nr_r1.m8"
    r1_cmd = [*common, "--query", str(ri.r1)]
    results.append(
        run_capture(
            r1_cmd,
            stage="r1_blastx",
            log_path=work_dir / "logs" / "01_r1_blastx.log",
            stdout_path=r1_out,
        )
    )
    ensure_success(results[-1])

    r2_out: Optional[Path] = None
    if ri.r2 is not None:
        r2_out = temp_dir / "diamond_nr_r2.m8"
        r2_cmd = [*common, "--query", str(ri.r2)]
        results.append(
            run_capture(
                r2_cmd,
                stage="r2_blastx",
                log_path=work_dir / "logs" / "02_r2_blastx.log",
                stdout_path=r2_out,
            )
        )
        ensure_success(results[-1])

    merged = work_dir / "diamond_merged_nr.m8"
    merge_log = work_dir / "logs" / "03_merge_m8.log"
    started_at = utc_now_iso()
    t0 = time.perf_counter()
    with merge_log.open("w", encoding="utf-8") as log:
        log.write(f"R1 source: {r1_out}\n")
        log.write(f"R2 source: {r2_out if r2_out else 'none'}\n")
        log.write(f"OUTPUT: {merged}\n")
        log.write("NOTE: current Rust merge copies R1 only.\n")
        with r1_out.open("rb") as src, merged.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    results.append(
        CommandResult(
            stage="merge_m8",
            command=["copy", str(r1_out), str(merged)],
            started_at=started_at,
            elapsed_seconds=time.perf_counter() - t0,
            returncode=0,
            stdout_bytes=0,
            stderr_bytes=merge_log.stat().st_size,
            log_path=str(merge_log),
        )
    )

    write_json(
        work_dir / "diamond_parameters.json",
        {
            "db": str(db),
            "db_prefix": str(db_prefix),
            "db_letters": letters,
            "threads": threads,
            "index_chunks": chunks,
            "block_size_gb": round(block_size, 3),
            "scratch_free_gib_at_start": round(scratch_free_gib, 3),
        },
        )
    return results


def run_mmseqs(
        ri: RunInput,
        work_dir: Path,
        backend: str,
        threads: int,
) -> list[CommandResult]:
    db = MMSEQS_CPU_DB if backend == "mmseqs-cpu" else MMSEQS_GPU_DB
    if not db.is_file():
        raise BenchmarkError(f"MMseqs DB does not exist: {db}")

    combined = work_dir / "mmseqs_non_host_combined.fastq"
    query_db = work_dir / "nr_merged.queryDB"
    result_db = work_dir / "nr_merged.resultDB"
    tmp_dir = work_dir / "nr_merged.tmp"
    m8_path = work_dir / "nr_merged.m8"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    results: list[CommandResult] = []

    # Match the current Rust write_combined_fastq() operation.
    log = work_dir / "logs" / "00_write_combined_fastq.log"
    started_at = utc_now_iso()
    t0 = time.perf_counter()
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as fh:
        fh.write(f"R1: {ri.r1}\n")
        fh.write(f"R2: {ri.r2 if ri.r2 else 'none'}\n")
        copy_fastq(ri.r1, ri.r2, combined)
        fh.write(f"OUTPUT: {combined}\n")
    results.append(
        CommandResult(
            stage="write_combined_fastq",
            command=["copy_fastq", str(ri.r1)]
                    + ([str(ri.r2)] if ri.r2 else [])
                    + [str(combined)],
            started_at=started_at,
            elapsed_seconds=time.perf_counter() - t0,
            returncode=0,
            stdout_bytes=0,
            stderr_bytes=log.stat().st_size,
            log_path=str(log),
        )
    )

    createdb = ["mmseqs", "createdb", str(combined), str(query_db), "--dbtype", "2"]
    results.append(
        run_capture(
            createdb,
            stage="createdb",
            log_path=work_dir / "logs" / "01_createdb.log",
        )
    )
    ensure_success(results[-1])
    if not query_db.is_file() or query_db.stat().st_size < 10_000:
        raise BenchmarkError(f"MMseqs query DB looks invalid: {query_db}")

    gpuserver_proc: Optional[subprocess.Popen[bytes]] = None
    server_log_fh = None
    if backend == "mmseqs-gpu":
        server_log = work_dir / "logs" / "02_gpuserver.log"
        server_log.parent.mkdir(parents=True, exist_ok=True)
        started_at = utc_now_iso()
        t0 = time.perf_counter()
        server_log_fh = server_log.open("wb")
        server_log_fh.write(f"COMMAND: mmseqs gpuserver {db}\n".encode())
        server_log_fh.flush()
        gpuserver_proc = subprocess.Popen(
            ["mmseqs", "gpuserver", str(db)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=server_log_fh,
        )
        try:
            for _ in range(20):
                if gpuserver_proc.poll() is not None:
                    raise BenchmarkError(
                        f"MMseqs GPU server exited during warm-up with status {gpuserver_proc.returncode}"
                    )
                time.sleep(1.0)
            results.append(
                CommandResult(
                    stage="gpuserver_warmup",
                    command=["mmseqs", "gpuserver", str(db)],
                    started_at=started_at,
                    elapsed_seconds=time.perf_counter() - t0,
                    returncode=0,
                    stdout_bytes=0,
                    stderr_bytes=server_log.stat().st_size,
                    log_path=str(server_log),
                )
            )
        except Exception:
            try:
                gpuserver_proc.kill()
                gpuserver_proc.wait(timeout=10)
            except Exception:
                pass
            raise
        finally:
            server_log_fh.close()
            server_log_fh = None

    search = ["mmseqs", "search"]
    if backend == "mmseqs-gpu":
        # Match the current Rust MMseqs GPU-server client construction:
        # --gpu 1 --gpu-server 1 --db-load-mode 2 --prefilter-mode 1
        search += [
            "--gpu", "1",
            "--gpu-server", "1",
            "--db-load-mode", "2",
        ]
    search += ["--threads", str(threads)]
    if backend == "mmseqs-cpu":
        search += ["-s", "5.7"]
    search += [
        "--alignment-mode", "3",
        "--search-type", "3",
        "--max-seqs", "1000" if backend == "mmseqs-cpu" else "3000",
    ]
    # CPU benchmark preserves the established production-mode prefilter.
    # GPU gpuserver requires prefilter mode 1 for the GPU-indexed database.
    search += ["--prefilter-mode", "0" if backend == "mmseqs-cpu" else "1", str(query_db),
               str(db),
               str(result_db),
               str(tmp_dir),
               "-e", "0.001",
               "--min-seq-id", "0.25",
               ]
    results.append(
        run_capture(
            search,
            stage="search",
            log_path=work_dir / "logs" / "03_search.log",
        )
    )

    if gpuserver_proc is not None:
        stop_log = work_dir / "logs" / "04_stop_gpuserver.log"
        started_at = utc_now_iso()
        t0 = time.perf_counter()
        stop_rc = None
        with stop_log.open("w", encoding="utf-8") as fh:
            try:
                gpuserver_proc.terminate()
                gpuserver_proc.wait(timeout=30)
                stop_rc = gpuserver_proc.returncode
            except subprocess.TimeoutExpired:
                gpuserver_proc.kill()
                gpuserver_proc.wait(timeout=10)
                stop_rc = gpuserver_proc.returncode
            fh.write(f"gpuserver returncode: {stop_rc}\n")
        results.append(
            CommandResult(
                stage="stop_gpuserver",
                command=["terminate", "mmseqs", "gpuserver"],
                started_at=started_at,
                elapsed_seconds=time.perf_counter() - t0,
                returncode=0 if stop_rc is not None else 1,
                stdout_bytes=0,
                stderr_bytes=stop_log.stat().st_size,
                log_path=str(stop_log),
            )
        )

    ensure_success(results[[r.stage for r in results].index("search")])

    convert = [
        "mmseqs", "convertalis",
        str(query_db), str(db), str(result_db), str(m8_path),
        "--format-output", M8_FORMAT,
    ]
    results.append(
        run_capture(
            convert,
            stage="convertalis",
            log_path=work_dir / "logs" / "05_convertalis.log",
        )
    )
    ensure_success(results[-1])
    if not m8_path.is_file() or m8_path.stat().st_size == 0:
        raise BenchmarkError(f"MMseqs produced no m8 output: {m8_path}")

    write_json(
        work_dir / "mmseqs_parameters.json",
        {
            "backend": backend,
            "db": str(db),
            "threads": threads,
            "createdb": createdb,
            "search": search,
            "convertalis": convert,
        },
        )
    return results


def machine_metadata(
        *,
        host_label: str,
        machine_type_override: str,
) -> dict[str, object]:
    instance_id = detect_instance_id()
    instance_type = machine_type_override or detect_instance_type()
    return {
        "hostname": os.uname().nodename,
        "instance_id": instance_id,
        "instance_type": instance_type,
        "cpu_model": cpu_model(),
        "cpu_count": os.cpu_count() or 0,
        "mem_total_bytes": read_meminfo_bytes(),
        "mem_available_bytes_at_start": available_ram_bytes(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "mmseqs_version": command_version(["mmseqs", "version"]),
        "diamond_version": command_version(["diamond", "version"]),
        "benchmark_host_label": host_label,
    }


def condition_dir(results_root: Path, backend: str, machine: dict[str, object]) -> Path:
    return (
            results_root
            / "conditions"
            / safe_component(backend)
            / safe_component(machine["instance_type"])
            / safe_component(machine["instance_id"])
    )


def trial_dir(
        results_root: Path,
        backend: str,
        machine: dict[str, object],
        trial_id: str,
        sample_id: str,
) -> Path:
    return (
            results_root
            / "trials"
            / safe_component(backend)
            / safe_component(machine["instance_type"])
            / safe_component(machine["instance_id"])
            / trial_id
            / safe_component(sample_id)
    )


def run_one_backend(
        ri: RunInput,
        backend: str,
        threads: int,
        scratch: Path,
        results_root: Path,
        machine: dict[str, object],
        keep_work: bool,
) -> None:
    trial_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    work_dir = trial_dir(results_root, backend, machine, trial_id, ri.sample_id)
    work_dir.mkdir(parents=True, exist_ok=True)

    condition = condition_dir(results_root, backend, machine)
    condition.mkdir(parents=True, exist_ok=True)
    write_json(condition / "machine_metadata.json", machine)

    started_wall = utc_now_iso()
    t0 = time.perf_counter()
    stage_results: list[CommandResult] = []
    error = ""
    status = "SUCCESS"

    try:
        if backend == "diamond":
            stage_results = run_diamond(ri, scratch, work_dir, threads)
        elif backend in {"mmseqs-cpu", "mmseqs-gpu"}:
            stage_results = run_mmseqs(ri, work_dir, backend, threads)
        else:
            raise BenchmarkError(f"Unknown backend: {backend}")
    except Exception as exc:
        status = "FAILED"
        error = str(exc)

    elapsed_total = time.perf_counter() - t0
    ended_wall = utc_now_iso()

    if backend == "diamond":
        candidates = [work_dir / "diamond_merged_nr.m8"]
    else:
        candidates = [work_dir / "nr_merged.m8"]
    output_path = next((p for p in candidates if p.exists()), None)

    r1_records: object = ""
    r2_records: object = ""
    try:
        r1_records = count_fastq_records(ri.r1)
    except Exception as exc:
        r1_records = f"ERROR: {exc}"
    if ri.r2:
        try:
            r2_records = count_fastq_records(ri.r2)
        except Exception as exc:
            r2_records = f"ERROR: {exc}"

    # Immutable per-trial manifest: this is the source of truth for identity.
    manifest = {
        "trial_id": trial_id,
        "started_at": started_wall,
        "ended_at": ended_wall,
        "experimental_condition": {
            "backend": backend,
            "instance_type": machine.get("instance_type", "unknown"),
            "instance_id": machine.get("instance_id", "unknown"),
            "host_label": machine.get("benchmark_host_label", ""),
        },
        "run_input": {
            "run_directory": ri.run_dir.name,
            "sample_id": ri.sample_id,
            "r1": str(ri.r1),
            "r2": str(ri.r2) if ri.r2 else None,
            "paired_end": bool(ri.r2),
            "r1_records": r1_records,
            "r2_records": r2_records,
        },
        "threads": threads,
        "machine": machine,
        "backend": backend,
        "status": status,
        "elapsed_seconds": elapsed_total,
        "stages": [asdict(x) for x in stage_results],
        "output_m8": {
            "path": str(output_path) if output_path else None,
            "bytes": output_path.stat().st_size if output_path else None,
            "rows": file_rows(output_path) if output_path else None,
        },
        "error": error,
    }
    write_json(work_dir / "manifest.json", manifest)

    # Per-condition accounting: no two different conditions ever share these files.
    stage_fields = [
        "trial_id", "sample_id", "run_directory", "backend", "instance_type",
        "instance_id", "host_label", "stage", "elapsed_seconds", "returncode",
        "command", "log_path",
    ]
    for sr in stage_results:
        append_csv_locked(
            condition / "benchmark_stages.csv",
            stage_fields,
            {
                "trial_id": trial_id,
                "sample_id": ri.sample_id,
                "run_directory": ri.run_dir.name,
                "backend": backend,
                "instance_type": machine.get("instance_type", "unknown"),
                "instance_id": machine.get("instance_id", "unknown"),
                "host_label": machine.get("benchmark_host_label", ""),
                "stage": sr.stage,
                "elapsed_seconds": f"{sr.elapsed_seconds:.6f}",
                "returncode": sr.returncode,
                "command": json.dumps(sr.command, separators=(",", ":")),
                "log_path": sr.log_path,
            },
            )

    summary_fields = [
        "trial_id", "started_at", "ended_at", "sample_id", "run_directory", "paired_end",
        "r1_path", "r2_path", "r1_records", "r2_records", "backend", "threads",
        "instance_type", "instance_id", "host_label", "hostname", "cpu_model", "cpu_count",
        "mem_total_bytes", "mem_available_bytes_at_start", "mmseqs_version", "diamond_version",
        "status", "elapsed_seconds", "output_m8", "output_m8_bytes", "output_m8_rows",
        "error", "work_directory",
    ]
    append_csv_locked(
        condition / "benchmark_results.csv",
        summary_fields,
        {
            "trial_id": trial_id,
            "started_at": started_wall,
            "ended_at": ended_wall,
            "sample_id": ri.sample_id,
            "run_directory": ri.run_dir.name,
            "paired_end": bool(ri.r2),
            "r1_path": str(ri.r1),
            "r2_path": str(ri.r2) if ri.r2 else "",
            "r1_records": r1_records,
            "r2_records": r2_records,
            "backend": backend,
            "threads": threads,
            "instance_type": machine.get("instance_type", "unknown"),
            "instance_id": machine.get("instance_id", "unknown"),
            "host_label": machine.get("benchmark_host_label", ""),
            "hostname": machine.get("hostname", "unknown"),
            "cpu_model": machine.get("cpu_model", "unknown"),
            "cpu_count": machine.get("cpu_count", ""),
            "mem_total_bytes": machine.get("mem_total_bytes", ""),
            "mem_available_bytes_at_start": machine.get("mem_available_bytes_at_start", ""),
            "mmseqs_version": machine.get("mmseqs_version", "unknown"),
            "diamond_version": machine.get("diamond_version", "unknown"),
            "status": status,
            "elapsed_seconds": f"{elapsed_total:.6f}",
            "output_m8": str(output_path) if output_path else "",
            "output_m8_bytes": output_path.stat().st_size if output_path else "",
            "output_m8_rows": file_rows(output_path) if output_path else "",
            "error": error,
            "work_directory": str(work_dir),
        },
        )

    # Also append to one global comparison table, but still preserve the per-condition tables.
    global_root = results_root / "global"
    append_csv_locked(
        global_root / "benchmark_results.csv",
        summary_fields,
        {
            "trial_id": trial_id,
            "started_at": started_wall,
            "ended_at": ended_wall,
            "sample_id": ri.sample_id,
            "run_directory": ri.run_dir.name,
            "paired_end": bool(ri.r2),
            "r1_path": str(ri.r1),
            "r2_path": str(ri.r2) if ri.r2 else "",
            "r1_records": r1_records,
            "r2_records": r2_records,
            "backend": backend,
            "threads": threads,
            "instance_type": machine.get("instance_type", "unknown"),
            "instance_id": machine.get("instance_id", "unknown"),
            "host_label": machine.get("benchmark_host_label", ""),
            "hostname": machine.get("hostname", "unknown"),
            "cpu_model": machine.get("cpu_model", "unknown"),
            "cpu_count": machine.get("cpu_count", ""),
            "mem_total_bytes": machine.get("mem_total_bytes", ""),
            "mem_available_bytes_at_start": machine.get("mem_available_bytes_at_start", ""),
            "mmseqs_version": machine.get("mmseqs_version", "unknown"),
            "diamond_version": machine.get("diamond_version", "unknown"),
            "status": status,
            "elapsed_seconds": f"{elapsed_total:.6f}",
            "output_m8": str(output_path) if output_path else "",
            "output_m8_bytes": output_path.stat().st_size if output_path else "",
            "output_m8_rows": file_rows(output_path) if output_path else "",
            "error": error,
            "work_directory": str(work_dir),
        },
        )

    if status == "SUCCESS" and not keep_work:
        for path in work_dir.glob("*.queryDB*"):
            try:
                path.unlink()
            except OSError:
                pass
        for path in work_dir.glob("*.resultDB*"):
            try:
                path.unlink()
            except OSError:
                pass
        shutil.rmtree(work_dir / "diamond_tmp", ignore_errors=True)
        shutil.rmtree(work_dir / "nr_merged.tmp", ignore_errors=True)

    print(
        f"{ri.sample_id}\t{backend}\t{machine.get('instance_type', 'unknown')}\t"
        f"{machine.get('instance_id', 'unknown')}\t{status}\t{elapsed_total:.2f}s"
        + (f"\t{error}" if error else "")
    )


def parse_backend_list(raw: str) -> list[str]:
    backends = [x.strip() for x in raw.split(",") if x.strip()]
    invalid = [x for x in backends if x not in VALID_BACKENDS]
    if invalid:
        raise BenchmarkError(f"Unknown backend(s): {', '.join(invalid)}")
    if not backends:
        raise BenchmarkError("At least one backend must be specified")
    return backends


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Benchmark current SeqToID NR alignment commands.")
    p.add_argument("--runs-root", type=Path, default=RUNS_ROOT)
    p.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH)
    p.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    p.add_argument(
        "--backends",
        default="diamond,mmseqs-cpu,mmseqs-gpu",
        help="Comma-separated: diamond,mmseqs-cpu,mmseqs-gpu",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Common thread count; 0 means os.cpu_count().",
    )
    p.add_argument("--diamond-threads", type=int, default=0)
    p.add_argument("--mmseqs-cpu-threads", type=int, default=0)
    p.add_argument("--mmseqs-gpu-threads", type=int, default=0)
    p.add_argument(
        "--host-label",
        default="",
        help="Human-readable label for this benchmark configuration. Instance ID is also recorded.",
    )
    p.add_argument(
        "--machine-type",
        default="",
        help="Override detected EC2 instance type when needed.",
    )
    p.add_argument(
        "--sample",
        action="append",
        help="Run only this derived sample ID; may be repeated. Omit for all discovered samples.",
    )
    p.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep intermediate query/result DBs and temporary files.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Discover and print samples without executing aligners.",
    )
    return p


def main() -> int:
    args = build_arg_parser().parse_args()
    try:
        backends = parse_backend_list(args.backends)
        inputs = discover_inputs(args.runs_root)
        if args.sample:
            wanted = set(args.sample)
            inputs = [x for x in inputs if x.sample_id in wanted]
        if not inputs:
            raise BenchmarkError(f"No benchmark inputs found under {args.runs_root}")

        print("Discovered benchmark inputs:")
        for ri in inputs:
            print(
                f"  {ri.sample_id}: R1={ri.r1}"
                + (f" R2={ri.r2}" if ri.r2 else " (single-end)")
            )

        print(f"Backends: {', '.join(backends)}")
        if args.dry_run:
            return 0

        if any(x.startswith("mmseqs") for x in backends):
            require_command("mmseqs")
        if "diamond" in backends:
            require_command("diamond")
        require_command("curl")

        args.results_root.mkdir(parents=True, exist_ok=True)
        if not args.scratch.exists():
            raise BenchmarkError(f"Scratch path does not exist: {args.scratch}")

        machine = machine_metadata(
            host_label=args.host_label,
            machine_type_override=args.machine_type,
        )

        # Store a host inventory record that is independent of backend.
        host_inventory_dir = args.results_root / "hosts"
        write_json(
            host_inventory_dir
            / f"{safe_component(machine['instance_type'])}.{safe_component(machine['instance_id'])}.json",
            machine,
            )

        default_threads = args.threads or (os.cpu_count() or 1)
        thread_map = {
            "diamond": args.diamond_threads or default_threads,
            "mmseqs-cpu": args.mmseqs_cpu_threads or default_threads,
            "mmseqs-gpu": args.mmseqs_gpu_threads or default_threads,
        }
        for backend, threads in thread_map.items():
            if threads < 1:
                raise BenchmarkError(f"Threads for {backend} must be >= 1")

        print(
            "Experimental host: "
            f"instance_type={machine['instance_type']} "
            f"instance_id={machine['instance_id']} "
            f"host_label={machine['benchmark_host_label']!r} "
            f"CPUs={machine['cpu_count']} RAM={machine['mem_total_bytes']}"
        )
        print(f"Results root: {args.results_root}")
        print("")

        for ri in inputs:
            for backend in backends:
                run_one_backend(
                    ri=ri,
                    backend=backend,
                    threads=thread_map[backend],
                    scratch=args.scratch,
                    results_root=args.results_root,
                    machine=machine,
                    keep_work=args.keep_work,
                )

        return 0
    except BenchmarkError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
