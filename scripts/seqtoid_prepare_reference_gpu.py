#!/usr/bin/env python3
"""Prepare and validate the SeqToID phase2 MMseqs GPU NR reference."""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

LOG = logging.getLogger("seqtoid-prepare-reference-gpu")
DEFAULT_LOCAL_ROOT = Path("/scratch/refs")
S5CMD_NUMWORKERS = "16"
S5CMD_CONCURRENCY = "80"
S5CMD_PART_SIZE = "50"
GPU_S3_PREFIX = "s3://seqtoid-public-references/phase2/refs/nr_clean/gpu"
GPU_REQUIRED_FILES = (
    "nrcleanDB_gpu",
    "nrcleanDB_gpu.dbtype",
    "nrcleanDB_gpu.idx",
    "nrcleanDB_gpu.idx.dbtype",
    "nrcleanDB_gpu.idx.index",
    "nrcleanDB_gpu.index",
    "nrcleanDB_gpu.lookup",
    "nrcleanDB_gpu.source",
    "nrcleanDB_gpu_h",
    "nrcleanDB_gpu_h.dbtype",
    "nrcleanDB_gpu_h.index",
)

@dataclass(frozen=True)
class ReferenceSpec:
    backend: str
    reference_set: str
    s3_source: str
    local_dir: Path
    db_name: str

class ReferencePreparationError(RuntimeError):
    pass

def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    LOG.debug("Running: %s", " ".join(args))
    try:
        return subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as exc:
        raise ReferencePreparationError(f"Required executable is missing: {args[0]}") from exc
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        raise ReferencePreparationError(
            f"Command failed with exit {exc.returncode}: {' '.join(args)}" + (f"; stderr: {stderr}" if stderr else "")
        ) from exc

def require_mountpoint(path: Path) -> None:
    result = subprocess.run(["mountpoint", "-q", str(path)], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise ReferencePreparationError(f"{path} is not mounted")

def validate_tools() -> None:
    for executable in ("s5cmd", "mmseqs"):
        if shutil.which(executable) is None:
            raise ReferencePreparationError(f"Required executable not found: {executable}")

def parse_inventory(raw: Iterable[str], prefix: str) -> list[dict[str, object]]:
    prefix = prefix.rstrip("/") + "/"
    objects: list[dict[str, object]] = []
    for line_number, raw_line in enumerate(raw, start=1):
        line = raw_line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ReferencePreparationError(f"Invalid JSON from s5cmd on line {line_number}: {exc}") from exc
        if not isinstance(payload, dict):
            raise ReferencePreparationError(f"Invalid JSON record from s5cmd on line {line_number}")
        if payload.get("error"):
            raise ReferencePreparationError(f"s5cmd reported an error on line {line_number}: {payload['error']}")
        key = payload.get("key")
        if not isinstance(key, str) or not key.startswith(prefix):
            continue
        rel = key[len(prefix):]
        if not rel or rel.endswith("/"):
            continue
        if "/" in rel:
            raise ReferencePreparationError(f"Unexpected nested reference object: {rel}")
        size = payload.get("size")
        if not isinstance(size, int) or size < 0:
            raise ReferencePreparationError(f"Missing or invalid size for reference object: {rel}")
        objects.append({"key": rel, "size": size, "etag": str(payload.get("etag") or "")})
    if not objects:
        raise ReferencePreparationError(f"No reference objects found below {prefix.rstrip('/')}")
    objects.sort(key=lambda obj: str(obj["key"]))
    return objects

def validate_gpu_inventory(objects: list[dict[str, object]]) -> None:
    names = {str(item["key"]) for item in objects}
    missing = sorted(set(GPU_REQUIRED_FILES) - names)
    extra = sorted(names - set(GPU_REQUIRED_FILES))
    if missing or extra:
        pieces: list[str] = []
        if missing: pieces.append(f"missing={missing}")
        if extra: pieces.append(f"extra={extra}")
        raise ReferencePreparationError("GPU MMseqs S3 inventory does not match expected database: " + "; ".join(pieces))

def build_inventory(spec: ReferenceSpec, output_path: Path) -> list[dict[str, object]]:
    output = subprocess.run(["s5cmd", "--json", "ls", f"{spec.s3_source.rstrip('/')}/*"], check=False, text=True, capture_output=True)
    if output.returncode != 0:
        raise ReferencePreparationError(f"Failed to list reference prefix {spec.s3_source}; stderr={output.stderr.strip()}")
    objects = parse_inventory(output.stdout.splitlines(), spec.s3_source)
    validate_gpu_inventory(objects)
    payload = {"backend": spec.backend, "reference_set": spec.reference_set, "s3_source": spec.s3_source, "objects": objects}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return objects

def inventory_version(inventory_path: Path) -> str:
    return hashlib.sha256(inventory_path.read_bytes()).hexdigest()

def load_inventory(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReferencePreparationError(f"Unable to read reference inventory {path}: {exc}") from exc
    objects = payload.get("objects")
    if not isinstance(objects, list) or not objects:
        raise ReferencePreparationError(f"Invalid reference inventory: {path}")
    return payload

def inventory_object_map(inventory: dict[str, object]) -> dict[str, tuple[int, str]]:
    objects = inventory.get("objects")
    if not isinstance(objects, list):
        raise ReferencePreparationError("Reference inventory has no object list")
    result: dict[str, tuple[int, str]] = {}
    for item in objects:
        if not isinstance(item, dict):
            raise ReferencePreparationError("Invalid object entry in reference inventory")
        key = item.get("key"); size = item.get("size"); etag = item.get("etag", "")
        if not isinstance(key, str) or not isinstance(size, int) or size < 0:
            raise ReferencePreparationError("Invalid object entry in reference inventory")
        result[key] = (size, str(etag))
    return result

def validate_local_inventory(local_dir: Path, inventory: dict[str, object]) -> None:
    expected = inventory_object_map(inventory)
    if not local_dir.is_dir():
        raise ReferencePreparationError(f"Reference directory is missing: {local_dir}")
    actual: dict[str, int] = {}
    for path in local_dir.iterdir():
        if path.name.startswith("."): continue
        if not path.is_file():
            raise ReferencePreparationError(f"Unexpected non-file entry in reference directory: {path.name}")
        actual[path.name] = path.stat().st_size
    expected_names = set(expected); actual_names = set(actual)
    missing = sorted(expected_names - actual_names)
    extra = sorted(actual_names - expected_names)
    wrong_size = sorted(name for name in expected_names & actual_names if expected[name][0] != actual[name])
    if missing or extra or wrong_size:
        problems = []
        if missing: problems.append(f"missing={missing[:20]}")
        if extra: problems.append(f"extra={extra[:20]}")
        if wrong_size: problems.append(f"wrong_size={wrong_size[:20]}")
        raise ReferencePreparationError("Local GPU reference does not match S3 inventory: " + "; ".join(problems))

def validate_mmseqs(spec: ReferenceSpec) -> None:
    db = spec.local_dir / spec.db_name
    if not db.is_file():
        raise ReferencePreparationError(f"MMseqs database base file is missing: {db}")
    result = subprocess.run(["mmseqs", "dbtype", str(db)], check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise ReferencePreparationError(f"mmseqs dbtype failed for {db}; stderr={result.stderr.strip()}")
    LOG.info("MMseqs GPU database validated: %s", db)

def write_metadata(reference_dir: Path, inventory: dict[str, object], version: str) -> None:
    (reference_dir / ".reference_version").write_text(version + "\n", encoding="utf-8")
    (reference_dir / ".reference_manifest.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def set_reference_permissions(reference_dir: Path) -> None:
    reference_dir.chmod(0o755)
    for path in reference_dir.iterdir():
        if path.is_dir(): path.chmod(0o755)
        elif path.is_file(): path.chmod(0o644)

def chown_cache_root(path: Path) -> None:
    try:
        import pwd, grp
        uid = pwd.getpwnam("ec2-user").pw_uid
        gid = grp.getgrnam("ec2-user").gr_gid
        os.chown(path, uid, gid)
    except (KeyError, PermissionError, OSError) as exc:
        raise ReferencePreparationError(f"Unable to set cache-root ownership on {path}: {exc}") from exc

def local_cache_is_current(spec: ReferenceSpec, inventory: dict[str, object], version: str) -> bool:
    version_path = spec.local_dir / ".reference_version"
    manifest_path = spec.local_dir / ".reference_manifest.json"
    if not version_path.is_file() or not manifest_path.is_file():
        return False
    try:
        if version_path.read_text(encoding="utf-8").strip() != version:
            return False
        if load_inventory(manifest_path) != inventory:
            return False
        validate_local_inventory(spec.local_dir, inventory)
        validate_mmseqs(spec)
    except ReferencePreparationError:
        return False
    return True

def stage_reference(spec: ReferenceSpec, inventory: dict[str, object], version: str) -> None:
    parent = spec.local_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging_root = parent / ".staging"
    staging_root.mkdir(parents=True, exist_ok=True)
    stage_dir = Path(tempfile.mkdtemp(prefix=f"{spec.backend}-", dir=staging_root))
    previous_dir = parent / f".previous-{os.getpid()}"
    try:
        LOG.info("Staging GPU reference to %s", stage_dir)
        run_command([
            "s5cmd", "--numworkers", S5CMD_NUMWORKERS,
            "cp", "--concurrency", S5CMD_CONCURRENCY,
            "--part-size", S5CMD_PART_SIZE,
            f"{spec.s3_source.rstrip('/')}/*", f"{stage_dir}/",
        ])
        validate_local_inventory(stage_dir, inventory)
        staged_spec = ReferenceSpec(spec.backend, spec.reference_set, spec.s3_source, stage_dir, spec.db_name)
        validate_mmseqs(staged_spec)
        write_metadata(stage_dir, inventory, version)
        set_reference_permissions(stage_dir)

        if previous_dir.exists(): shutil.rmtree(previous_dir)
        if spec.local_dir.exists(): os.replace(spec.local_dir, previous_dir)
        try:
            os.replace(stage_dir, spec.local_dir)
        except OSError:
            if previous_dir.exists() and not spec.local_dir.exists():
                os.replace(previous_dir, spec.local_dir)
            raise
        if previous_dir.exists(): shutil.rmtree(previous_dir)
        chown_cache_root(spec.local_dir)
    except Exception:
        if stage_dir.exists(): shutil.rmtree(stage_dir, ignore_errors=True)
        if previous_dir.exists() and not spec.local_dir.exists():
            try: os.replace(previous_dir, spec.local_dir)
            except OSError: LOG.exception("Failed to restore previous GPU reference cache")
        raise

def prepare_reference(backend: str, reference_set: str, local_root: Path, inventory_only: bool = False) -> str:
    if backend != "mmseqs-gpu":
        raise ReferencePreparationError(f"Unsupported backend for GPU preparer: {backend}")
    if reference_set != "phase2":
        raise ReferencePreparationError(f"Unsupported reference set: {reference_set}")
    require_mountpoint(Path("/scratch"))
    validate_tools()
    spec = ReferenceSpec(backend, reference_set, GPU_S3_PREFIX, local_root / "mmseqs-gpu", "nrcleanDB_gpu")
    with tempfile.TemporaryDirectory(prefix="seqtoid-reference-") as tmp:
        inventory_path = Path(tmp) / "reference-manifest.json"
        objects = build_inventory(spec, inventory_path)
        inventory = load_inventory(inventory_path)
        if objects != inventory["objects"]:
            raise ReferencePreparationError("Reference inventory changed during build")
        version = inventory_version(inventory_path)
        LOG.info("Canonical GPU MMseqs reference version: %s", version)
        if inventory_only:
            return version
        if local_cache_is_current(spec, inventory, version):
            LOG.info("Existing GPU MMseqs reference is current and validated")
            return version
        stage_reference(spec, inventory, version)
    promoted_inventory = load_inventory(spec.local_dir / ".reference_manifest.json")
    validate_local_inventory(spec.local_dir, promoted_inventory)
    validate_mmseqs(spec)
    chown_cache_root(spec.local_dir)
    return version

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare and validate a SeqToID MMseqs GPU worker reference.")
    parser.add_argument("--backend", required=True, choices=["mmseqs-gpu"])
    parser.add_argument("--reference-set", default="phase2")
    parser.add_argument("--local-root", type=Path, default=DEFAULT_LOCAL_ROOT)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--inventory-only", action="store_true")
    return parser.parse_args()

def main() -> int:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%dT%H:%M:%S%z")
    try:
        version = prepare_reference(args.backend, args.reference_set, args.local_root, args.inventory_only)
    except ReferencePreparationError as exc:
        LOG.error("Reference preparation failed: %s", exc)
        return 1
    except Exception:
        LOG.exception("Unexpected reference preparation failure")
        return 1
    print(version)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
