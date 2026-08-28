#!/usr/bin/env python3
"""Prepare and validate a SeqToID worker reference on local scratch NVMe.

The worker bootstrap is responsible for OS setup (/scratch, EFS, s5cmd, etc.).
This utility owns reference inventory, staging, validation, and atomic promotion.
"""

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


LOG = logging.getLogger("seqtoid-prepare-reference")
DEFAULT_LOCAL_ROOT = Path("/scratch/refs")
S5CMD_NUMWORKERS = "16"
S5CMD_CONCURRENCY = "80"
S5CMD_PART_SIZE = "50"

MMSEQS_CPU_S3_PREFIX = (
    "s3://seqtoid-public-references/phase2/refs/nr_clean/cpu"
)
DIAMOND_S3_OBJECT = (
    "s3://seqtoid-public-references/phase2/refs/diamond/diamond_07_22_2026.dmnd"
)


@dataclass(frozen=True)
class ReferenceSpec:
    backend: str
    reference_set: str
    s3_source: str
    local_dir: Path
    validator: str
    db_name: str


class ReferencePreparationError(RuntimeError):
    """Raised when a reference cannot be prepared and validated."""


def run_command(
        args: list[str],
        *,
        stdout: object = subprocess.PIPE,
        capture_stderr: bool = True,
) -> subprocess.CompletedProcess[str]:
    LOG.debug("Running: %s", " ".join(args))
    try:
        return subprocess.run(
            args,
            check=True,
            text=True,
            stdout=stdout,
            stderr=subprocess.PIPE if capture_stderr else None,
        )
    except FileNotFoundError as exc:
        raise ReferencePreparationError(
            f"Required executable is missing: {args[0]}"
        ) from exc
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip() if capture_stderr else ""
        raise ReferencePreparationError(
            f"Command failed with exit {exc.returncode}: {' '.join(args)}"
            + (f"; stderr: {stderr}" if stderr else "")
        ) from exc


def require_mountpoint(path: Path) -> None:
    result = subprocess.run(
        ["mountpoint", "-q", str(path)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        raise ReferencePreparationError(f"{path} is not mounted")


def validate_tools(spec: ReferenceSpec) -> None:
    for executable in ("s5cmd", spec.validator):
        if shutil.which(executable) is None:
            raise ReferencePreparationError(
                f"Required executable not found for {spec.backend}: {executable}"
            )


def parse_s5cmd_inventory(raw: Iterable[str], prefix: str) -> list[dict[str, object]]:
    """Parse s5cmd --json ls output into a deterministic object inventory."""
    prefix = prefix.rstrip("/") + "/"
    objects: list[dict[str, object]] = []

    for line_number, raw_line in enumerate(raw, start=1):
        line = raw_line.strip()
        if not line:
            continue

        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ReferencePreparationError(
                f"Invalid JSON from s5cmd on line {line_number}: {exc}"
            ) from exc

        item = payload

        key = item.get("key")
        if isinstance(key, dict):
            key = key.get("url") or key.get("key")
        if not isinstance(key, str):
            continue

        if not key.startswith(prefix):
            continue

        relative = key[len(prefix):]
        if not relative or relative.endswith("/"):
            continue
        if "/" in relative:
            raise ReferencePreparationError(
                f"Unexpected nested reference object below {prefix}: {relative}"
            )

        size = item.get("size")
        if not isinstance(size, int) or size < 0:
            raise ReferencePreparationError(
                f"Missing or invalid size for reference object: {relative}"
            )

        etag = item.get("etag") or ""
        objects.append({"key": relative, "size": size, "etag": str(etag)})

    if not objects:
        raise ReferencePreparationError(
            f"No reference objects found below {prefix.rstrip('/') }"
        )

    objects.sort(key=lambda obj: str(obj["key"]))
    return objects


def build_inventory(spec: ReferenceSpec, output_path: Path) -> list[dict[str, object]]:
    """Query S3 and write a canonical deterministic reference inventory."""
    if spec.backend == "diamond":
        objects = build_diamond_inventory(spec)
    else:
        output = subprocess.run(
            ["s5cmd", "--json", "ls", f"{spec.s3_source.rstrip('/')}/*"],
            check=False,
            text=True,
            capture_output=True,
        )
        if output.returncode != 0:
            stderr = output.stderr.strip()
            raise ReferencePreparationError(
                f"Failed to list reference prefix {spec.s3_source}"
                + (f"; stderr: {stderr}" if stderr else "")
            )
        objects = parse_s5cmd_inventory(output.stdout.splitlines(), spec.s3_source)

    payload = {
        "backend": spec.backend,
        "reference_set": spec.reference_set,
        "s3_source": spec.s3_source,
        "objects": objects,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        )
    return objects


def build_diamond_inventory(spec: ReferenceSpec) -> list[dict[str, object]]:
    """Build an inventory for the single Diamond database object."""
    output = subprocess.run(
        ["s5cmd", "--json", "ls", spec.s3_source],
        check=False,
        text=True,
        capture_output=True,
    )

    if output.returncode != 0:
        stderr = output.stderr.strip()
        raise ReferencePreparationError(
            f"Failed to inspect Diamond reference {spec.s3_source}"
            + (f"; stderr: {stderr}" if stderr else "")
        )

    for raw_line in output.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ReferencePreparationError(
                f"Invalid JSON from s5cmd for Diamond reference: {exc}"
            ) from exc

        item = payload

        key = item.get("key")
        if isinstance(key, dict):
            key = key.get("url") or key.get("key")

        if key != spec.s3_source:
            continue

        size = item.get("size")
        if not isinstance(size, int) or size < 0:
            raise ReferencePreparationError(
                f"Missing or invalid size for Diamond reference: {spec.s3_source}"
            )

        return [{
            "key": Path(spec.s3_source).name,
            "size": size,
            "etag": str(item.get("etag") or ""),
        }]

    raise ReferencePreparationError(
        f"No usable S3 metadata found for Diamond reference: {spec.s3_source}"
    )


def inventory_version(inventory_path: Path) -> str:
    """Return the SHA256 fingerprint of the canonical inventory JSON."""
    digest = hashlib.sha256()
    digest.update(inventory_path.read_bytes())
    return digest.hexdigest()


def load_inventory(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReferencePreparationError(
            f"Unable to read reference inventory {path}: {exc}"
        ) from exc

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
        key = item.get("key")
        size = item.get("size")
        etag = item.get("etag", "")
        if not isinstance(key, str) or not isinstance(size, int) or size < 0:
            raise ReferencePreparationError("Invalid object entry in reference inventory")
        result[key] = (size, str(etag))
    return result


def validate_local_inventory(local_dir: Path, inventory: dict[str, object]) -> None:
    """Ensure local files exactly match the expected S3 names and sizes."""
    expected = inventory_object_map(inventory)
    actual: dict[str, int] = {}

    if not local_dir.is_dir():
        raise ReferencePreparationError(f"Reference directory is missing: {local_dir}")

    for path in local_dir.iterdir():
        if path.name.startswith("."):
            continue
        if not path.is_file():
            raise ReferencePreparationError(
                f"Unexpected non-file entry in reference directory: {path.name}"
            )
        actual[path.name] = path.stat().st_size

    expected_names = set(expected)
    actual_names = set(actual)

    missing = sorted(expected_names - actual_names)
    extra = sorted(actual_names - expected_names)
    wrong_size = sorted(
        name for name in expected_names & actual_names
        if expected[name][0] != actual[name]
    )

    if missing or extra or wrong_size:
        problems: list[str] = []
        if missing:
            problems.append(f"missing={missing[:10]}")
        if extra:
            problems.append(f"extra={extra[:10]}")
        if wrong_size:
            problems.append(f"wrong_size={wrong_size[:10]}")
        raise ReferencePreparationError(
            "Local reference does not match S3 inventory: " + "; ".join(problems)
        )


def validate_mmseqs(spec: ReferenceSpec) -> None:
    """Validate that MMseqs can open the expected CPU NR database."""
    db_prefix = spec.local_dir / spec.db_name
    if not db_prefix.is_file():
        raise ReferencePreparationError(
            f"MMseqs database base file is missing: {db_prefix}"
        )

    result = subprocess.run(
        ["mmseqs", "dbtype", str(db_prefix)],
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        stderr = result.stderr.strip()
        raise ReferencePreparationError(
            f"mmseqs dbtype failed for {db_prefix}"
            + (f"; stderr: {stderr}" if stderr else "")
        )

    LOG.info("MMseqs database validated: %s", db_prefix)


def validate_diamond(spec: ReferenceSpec, inventory: dict[str, object]) -> None:
    """Validate that the expected Diamond database is present and usable."""
    db_path = spec.local_dir / spec.db_name
    expected = inventory_object_map(inventory)
    expected_meta = expected.get(spec.db_name)

    if not db_path.is_file():
        raise ReferencePreparationError(
            f"Diamond database is missing: {db_path}"
        )

    if expected_meta is None:
        raise ReferencePreparationError(
            f"Diamond database {spec.db_name} is absent from reference inventory"
        )

    actual_size = db_path.stat().st_size
    if actual_size != expected_meta[0]:
        raise ReferencePreparationError(
            f"Diamond database size mismatch for {db_path}: "
            f"expected {expected_meta[0]}, found {actual_size}"
        )

    if actual_size == 0:
        raise ReferencePreparationError(f"Diamond database is empty: {db_path}")

    result = subprocess.run(
        ["diamond", "dbinfo", "--db", str(db_path)],
        check=False,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        stderr = result.stderr.strip()
        raise ReferencePreparationError(
            f"diamond dbinfo failed for {db_path}"
            + (f"; stderr: {stderr}" if stderr else "")
        )

    LOG.info("Diamond database validated: %s (%d bytes)", db_path, actual_size)


def write_metadata(reference_dir: Path, inventory: dict[str, object], version: str) -> None:
    """Write local reference metadata only after successful validation."""
    (reference_dir / ".reference_version").write_text(version + "\n", encoding="utf-8")
    (reference_dir / ".reference_manifest.json").write_text(
        json.dumps(inventory, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        )


def validate_local_reference(spec: ReferenceSpec, inventory: dict[str, object]) -> None:
    validate_local_inventory(spec.local_dir, inventory)
    if spec.validator == "mmseqs":
        validate_mmseqs(spec)
    elif spec.validator == "diamond":
        validate_diamond(spec, inventory)
    else:
        raise ReferencePreparationError(
            f"Unsupported reference validator: {spec.validator}"
        )


def local_cache_is_current(
        spec: ReferenceSpec,
        inventory: dict[str, object],
        version: str,
) -> bool:
    """Return true only when local metadata, file inventory, and DB validation all pass."""
    version_path = spec.local_dir / ".reference_version"
    manifest_path = spec.local_dir / ".reference_manifest.json"

    if not version_path.is_file() or not manifest_path.is_file():
        return False

    try:
        recorded_version = version_path.read_text(encoding="utf-8").strip()
        recorded_inventory = load_inventory(manifest_path)
    except ReferencePreparationError:
        return False

    if recorded_version != version:
        return False

    if recorded_inventory != inventory:
        return False

    try:
        validate_local_reference(spec, inventory)
    except ReferencePreparationError:
        return False

    return True



def set_reference_permissions(reference_dir: Path) -> None:
    """Make the validated reference readable by worker processes."""
    reference_dir.chmod(0o755)

    for path in reference_dir.iterdir():
        if path.is_dir():
            path.chmod(0o755)
        elif path.is_file():
            path.chmod(0o644)

def stage_reference(spec: ReferenceSpec, inventory: dict[str, object], version: str) -> None:
    """Stage, validate, and atomically promote a complete reference cache."""
    spec.local_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_root = spec.local_dir.parent / ".staging"
    staging_root.mkdir(parents=True, exist_ok=True)

    stage_dir = Path(tempfile.mkdtemp(prefix=f"{spec.backend}-", dir=staging_root))
    previous_dir = spec.local_dir.parent / f".previous-{os.getpid()}"

    try:
        LOG.info("Staging %s reference to %s", spec.backend, stage_dir)
        if spec.backend == "diamond":
            run_command(
                [
                    "s5cmd",
                    "cp",
                    spec.s3_source,
                    f"{stage_dir}/",
                ],
                stdout=subprocess.PIPE,
            )
        else:
            run_command(
                [
                    "s5cmd",
                    "--numworkers",
                    S5CMD_NUMWORKERS,
                    "cp",
                    "--concurrency",
                    S5CMD_CONCURRENCY,
                    "--part-size",
                    S5CMD_PART_SIZE,
                    f"{spec.s3_source.rstrip('/')}/*",
                    f"{stage_dir}/",
                ],
                stdout=subprocess.PIPE,
            )

        LOG.info("Validating staged file inventory")
        staged_inventory = dict(inventory)
        validate_local_inventory(stage_dir, staged_inventory)

        staged_spec = ReferenceSpec(
            backend=spec.backend,
            reference_set=spec.reference_set,
            s3_source=spec.s3_source,
            local_dir=stage_dir,
            validator=spec.validator,
            db_name=spec.db_name,
        )
        LOG.info("Validating staged %s database", spec.backend)
        validate_local_reference(staged_spec, staged_inventory)
        write_metadata(stage_dir, staged_inventory, version)
        set_reference_permissions(stage_dir)

        LOG.info("Promoting validated reference to %s", spec.local_dir)
        if previous_dir.exists():
            shutil.rmtree(previous_dir)

        if spec.local_dir.exists():
            os.replace(spec.local_dir, previous_dir)

        try:
            os.replace(stage_dir, spec.local_dir)
        except OSError:
            if previous_dir.exists() and not spec.local_dir.exists():
                os.replace(previous_dir, spec.local_dir)
            raise

        if previous_dir.exists():
            shutil.rmtree(previous_dir)

    except Exception:
        if stage_dir.exists():
            shutil.rmtree(stage_dir, ignore_errors=True)
        if previous_dir.exists() and not spec.local_dir.exists():
            try:
                os.replace(previous_dir, spec.local_dir)
            except OSError:
                LOG.exception("Failed to restore previous reference cache")
        raise


def make_spec(backend: str, reference_set: str, local_root: Path) -> ReferenceSpec:
    if reference_set != "phase2":
        raise ReferencePreparationError(
            f"Unsupported reference set: {reference_set}"
        )

    if backend == "mmseqs-cpu":
        return ReferenceSpec(
            backend=backend,
            reference_set=reference_set,
            s3_source=MMSEQS_CPU_S3_PREFIX,
            local_dir=local_root / "mmseqs",
            validator="mmseqs",
            db_name="nrcleanDB",
        )

    if backend == "diamond":
        return ReferenceSpec(
            backend=backend,
            reference_set=reference_set,
            s3_source=DIAMOND_S3_OBJECT,
            local_dir=local_root / "diamond",
            validator="diamond",
            db_name=Path(DIAMOND_S3_OBJECT).name,
        )

    raise ReferencePreparationError(f"Unsupported backend: {backend}")


def prepare_reference(
        backend: str,
        reference_set: str,
        local_root: Path,
        inventory_only: bool = False,
) -> str:
    """Prepare the requested reference and return its canonical inventory version."""
    require_mountpoint(Path("/scratch"))
    spec = make_spec(backend, reference_set, local_root)
    validate_tools(spec)

    with tempfile.TemporaryDirectory(prefix="seqtoid-reference-") as tmp:
        inventory_path = Path(tmp) / "reference-manifest.json"
        objects = build_inventory(spec, inventory_path)
        inventory = load_inventory(inventory_path)

        if objects != inventory["objects"]:
            raise ReferencePreparationError(
                "Reference inventory changed during build"
            )

        version = inventory_version(inventory_path)

        LOG.info(
            "Canonical %s reference version: %s",
            spec.backend,
            version,
        )

        if inventory_only:
            LOG.info(
                "Inventory-only mode requested; not modifying local reference cache"
            )
            return version

        if local_cache_is_current(spec, inventory, version):
            LOG.info(
                "Existing %s reference is current and validated",
                spec.backend,
            )
            return version

        LOG.info(
            "Existing %s reference is absent, stale, or invalid",
            spec.backend,
        )
        stage_reference(spec, inventory, version)

    promoted_inventory = load_inventory(
        spec.local_dir / ".reference_manifest.json"
    )
    validate_local_reference(
        spec,
        promoted_inventory,
    )

    LOG.info(
        "%s reference installed and validated",
        spec.backend,
    )
    LOG.info(
        "  S3 source: %s",
        spec.s3_source,
    )
    LOG.info(
        "  Local path: %s",
        spec.local_dir,
    )
    LOG.info(
        "  Version: %s",
        version,
    )

    return version


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare and validate a SeqToID worker reference."
    )
    parser.add_argument(
        "--backend",
        required=True,
        choices=["mmseqs-cpu", "diamond"],
        help="Reference backend to prepare.",
    )
    parser.add_argument(
        "--reference-set",
        default="phase2",
        help="Reference set to prepare.",
    )
    parser.add_argument(
        "--local-root",
        type=Path,
        default=DEFAULT_LOCAL_ROOT,
        help="Root directory for local worker reference caches.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help=(
            "Inspect and fingerprint the canonical S3 reference without "
            "downloading or modifying the local reference cache."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
    )

    try:
        version = prepare_reference(
            args.backend,
            args.reference_set,
            args.local_root,
            args.inventory_only,
        )
    except ReferencePreparationError as exc:
        LOG.error(
            "Reference preparation failed: %s",
            exc,
        )
        return 1
    except Exception:
        LOG.exception(
            "Unexpected reference preparation failure"
        )
        return 1

    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())