#!/usr/bin/env python3
"""
Prepare and validate the SeqToID launch-node reference cache on local scratch.

The launch-node bootstrap owns OS/runtime setup:
    - /scratch
    - /efs
    - /refs -> /scratch/refs
    - s5cmd
    - required runtime software

This utility owns the launch-node reference cache:
    - inventory of the canonical S3 reference tree
    - exclusion of worker-owned reference trees
    - staging
    - local inventory validation
    - atomic promotion
    - deterministic reference versioning

Canonical source:

    s3://seqtoid-public-references/phase2/refs/

Excluded prefixes:

    nr_clean/cpu/
    nr_clean/gpu/
    diamond/

Everything else under phase2/refs/ is included.

In particular, the launch node DOES include:

    nr_clean/nr_clean.fa
    nr_clean/nr_clean_index.fst

s5cmd 2.3.0 behavior:

For zero-byte S3 objects, s5cmd --json ls may omit the "size"
field entirely. Such records are accepted as zero-byte objects only
when they are explicitly files and their ETag is the standard MD5
digest of an empty object.
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
from pathlib import Path
from typing import Iterable


LOG = logging.getLogger("seqtoid-prepare-launch-references")


DEFAULT_LOCAL_ROOT = Path("/scratch/refs")

S3_REFERENCE_ROOT = (
    "s3://seqtoid-public-references/phase2/refs"
)

S5CMD_NUMWORKERS = "16"
S5CMD_CONCURRENCY = "80"
S5CMD_PART_SIZE = "50"

EMPTY_OBJECT_ETAG = (
    "d41d8cd98f00b204e9800998ecf8427e"
)


# These are the ONLY excluded subdirectories.
#
# Paths are relative to:
#
#     s3://seqtoid-public-references/phase2/refs/
#
EXCLUDED_PREFIXES = (
    "nr_clean/cpu/",
    "nr_clean/gpu/",
    "diamond/",
)


class ReferencePreparationError(RuntimeError):
    """Raised when launch-node reference preparation fails."""


def run_command(
        args: list[str],
        *,
        stdout: object = subprocess.PIPE,
        capture_stderr: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a command and convert failures into a descriptive exception."""

    LOG.debug(
        "Running: %s",
        " ".join(args),
    )

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
        stderr = (
            (exc.stderr or "").strip()
            if capture_stderr
            else ""
        )

        raise ReferencePreparationError(
            f"Command failed with exit {exc.returncode}: "
            f"{' '.join(args)}"
            + (
                f"; stderr: {stderr}"
                if stderr
                else ""
            )
        ) from exc


def require_mountpoint(path: Path) -> None:
    """Require the target path to be a mounted filesystem."""

    result = subprocess.run(
        [
            "mountpoint",
            "-q",
            str(path),
        ],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    if result.returncode != 0:
        raise ReferencePreparationError(
            f"{path} is not mounted"
        )


def validate_tools() -> None:
    """Ensure the launch-node reference preparer has its required tool."""

    if shutil.which("s5cmd") is None:
        raise ReferencePreparationError(
            "Required executable is missing: s5cmd"
        )


def is_excluded(relative_key: str) -> bool:
    """
    Return True when a relative S3 key belongs to one of the
    three explicitly excluded reference subdirectories.
    """

    normalized = relative_key.lstrip("/")

    return any(
        normalized.startswith(prefix)
        for prefix in EXCLUDED_PREFIXES
    )


def s3_key_to_relative(
        key: str,
        root: str,
) -> str | None:
    """
    Convert an absolute S3 URI into a relative key below the
    canonical launch-node reference root.

    Example:

        s3://bucket/phase2/refs/core_nt/core_nt.fa

    becomes:

        core_nt/core_nt.fa
    """

    prefix = root.rstrip("/") + "/"

    if not key.startswith(prefix):
        return None

    relative = key[len(prefix):].lstrip("/")

    if not relative:
        return None

    return relative


def get_s5cmd_object_size(
        payload: dict[str, object],
        relative: str,
) -> int:
    """
    Return the object size from an s5cmd JSON record.

    s5cmd 2.3.0 normally emits:

        "size": <integer>

    However, for a legitimate zero-byte S3 object it may omit the
    size field entirely, for example:

        {
            "key": ".../library.fna.masked",
            "etag": "d41d8cd98f00b204e9800998ecf8427e",
            "type": "file"
        }

    A missing size is interpreted as zero ONLY when:

        type == "file"

    and:

        etag == MD5(empty object)

    Any other missing or malformed size remains a hard failure.
    """

    size = payload.get("size")

    if size is None:
        object_type = payload.get("type")
        etag = payload.get("etag")

        if (
                object_type == "file"
                and etag == EMPTY_OBJECT_ETAG
        ):
            LOG.debug(
                "s5cmd omitted size for zero-byte object: %s",
                relative,
            )

            return 0

        raise ReferencePreparationError(
            "Missing size for non-empty or unrecognized "
            f"reference object: {relative}"
        )

    if not isinstance(size, int):
        raise ReferencePreparationError(
            "Invalid non-integer size for reference object: "
            f"{relative}"
        )

    if size < 0:
        raise ReferencePreparationError(
            "Invalid negative size for reference object: "
            f"{relative}"
        )

    return size


def parse_s5cmd_inventory(
        raw: Iterable[str],
) -> list[dict[str, object]]:
    """
    Parse s5cmd --json ls output.

    s5cmd 2.3.0 emits each object as a top-level JSON record:

        {
            "key": "...",
            "etag": "...",
            "size": 123,
            ...
        }

    Nested directory structure is intentional here and is preserved.
    """

    objects: list[dict[str, object]] = []

    for line_number, raw_line in enumerate(
            raw,
            start=1,
    ):
        line = raw_line.strip()

        if not line:
            continue

        try:
            payload = json.loads(line)

        except json.JSONDecodeError as exc:
            raise ReferencePreparationError(
                "Invalid JSON from s5cmd on line "
                f"{line_number}: {exc}"
            ) from exc

        if not isinstance(payload, dict):
            raise ReferencePreparationError(
                "Invalid JSON record from s5cmd on line "
                f"{line_number}"
            )

        if payload.get("error"):
            raise ReferencePreparationError(
                "s5cmd reported an error on line "
                f"{line_number}: {payload['error']}"
            )

        key = payload.get("key")

        if not isinstance(key, str):
            continue

        relative = s3_key_to_relative(
            key,
            S3_REFERENCE_ROOT,
        )

        if relative is None:
            continue

        # These are intentionally the ONLY exclusions.
        if is_excluded(relative):
            continue

        # Ignore directory markers.
        if relative.endswith("/"):
            continue

        size = get_s5cmd_object_size(
            payload,
            relative,
        )

        etag = payload.get("etag") or ""

        objects.append(
            {
                "key": relative,
                "size": size,
                "etag": str(etag),
            }
        )

    if not objects:
        raise ReferencePreparationError(
            "No launch-node reference objects found below "
            f"{S3_REFERENCE_ROOT}"
        )

    objects.sort(
        key=lambda obj: str(obj["key"])
    )

    return objects


def build_inventory(
        output_path: Path,
) -> list[dict[str, object]]:
    """
    Query the canonical S3 tree and write the deterministic
    launch-node inventory.
    """

    LOG.info(
        "Building launch-node reference inventory from %s",
        S3_REFERENCE_ROOT,
    )

    output = subprocess.run(
        [
            "s5cmd",
            "--json",
            "ls",
            f"{S3_REFERENCE_ROOT}/*",
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    if output.returncode != 0:
        stderr = output.stderr.strip()

        raise ReferencePreparationError(
            "Failed to list launch-node reference tree "
            f"{S3_REFERENCE_ROOT}"
            + (
                f"; stderr: {stderr}"
                if stderr
                else ""
            )
        )

    objects = parse_s5cmd_inventory(
        output.stdout.splitlines()
    )

    payload = {
        "backend": "launch-node",
        "reference_set": "phase2",
        "s3_source": S3_REFERENCE_ROOT,
        "excluded_prefixes": list(
            EXCLUDED_PREFIXES
        ),
        "objects": objects,
    }

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
        )

    LOG.info(
        "Launch-node S3 inventory is complete: %d included objects",
        len(objects),
    )

    LOG.info(
        "Excluded prefixes:"
    )

    for prefix in EXCLUDED_PREFIXES:
        LOG.info(
            "  %s",
            prefix,
        )

    return objects


def inventory_version(
        inventory_path: Path,
) -> str:
    """Return SHA256 fingerprint of the canonical inventory JSON."""

    digest = hashlib.sha256()

    digest.update(
        inventory_path.read_bytes()
    )

    return digest.hexdigest()


def inventory_version_from_payload(
        inventory: dict[str, object],
) -> str:
    """
    Compute the deterministic inventory version directly from its
    canonical JSON representation.
    """

    canonical = (
            json.dumps(
                inventory,
                indent=2,
                sort_keys=True,
            )
            + "\n"
    ).encode(
        "utf-8"
    )

    return hashlib.sha256(
        canonical
    ).hexdigest()


def load_inventory(
        path: Path,
) -> dict[str, object]:
    """Load and minimally validate an inventory file."""

    try:
        payload = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

    except (
            OSError,
            json.JSONDecodeError,
    ) as exc:
        raise ReferencePreparationError(
            f"Unable to read reference inventory {path}: {exc}"
        ) from exc

    objects = payload.get("objects")

    if not isinstance(objects, list) or not objects:
        raise ReferencePreparationError(
            f"Invalid launch-node reference inventory: {path}"
        )

    return payload


def inventory_object_map(
        inventory: dict[str, object],
) -> dict[str, tuple[int, str]]:
    """Return relative key -> (size, etag)."""

    objects = inventory.get("objects")

    if not isinstance(objects, list):
        raise ReferencePreparationError(
            "Reference inventory has no object list"
        )

    result: dict[str, tuple[int, str]] = {}

    for item in objects:
        if not isinstance(item, dict):
            raise ReferencePreparationError(
                "Invalid object entry in reference inventory"
            )

        key = item.get("key")
        size = item.get("size")
        etag = item.get("etag", "")

        if (
                not isinstance(key, str)
                or not isinstance(size, int)
                or size < 0
        ):
            raise ReferencePreparationError(
                "Invalid object entry in reference inventory"
            )

        if is_excluded(key):
            raise ReferencePreparationError(
                "Inventory unexpectedly contains excluded object: "
                f"{key}"
            )

        if key in result:
            raise ReferencePreparationError(
                "Duplicate object in launch-node reference inventory: "
                f"{key}"
            )

        result[key] = (
            size,
            str(etag),
        )

    return result


def validate_inventory_policy(
        inventory: dict[str, object],
) -> None:
    """
    Ensure the inventory encodes exactly the intended inclusion and
    exclusion policy.
    """

    if inventory.get("backend") != "launch-node":
        raise ReferencePreparationError(
            "Launch-node inventory has an unexpected backend"
        )

    if inventory.get("reference_set") != "phase2":
        raise ReferencePreparationError(
            "Launch-node inventory has an unexpected reference set"
        )

    if inventory.get("s3_source") != S3_REFERENCE_ROOT:
        raise ReferencePreparationError(
            "Launch-node inventory has an unexpected S3 source"
        )

    excluded = inventory.get(
        "excluded_prefixes"
    )

    expected = list(
        EXCLUDED_PREFIXES
    )

    if excluded != expected:
        raise ReferencePreparationError(
            "Launch-node reference exclusion policy does not match "
            "the expected three prefixes"
        )

    objects = inventory.get("objects")

    if not isinstance(objects, list) or not objects:
        raise ReferencePreparationError(
            "Launch-node inventory has no objects"
        )

    for item in objects:
        if not isinstance(item, dict):
            raise ReferencePreparationError(
                "Invalid object in launch-node inventory"
            )

        key = item.get("key")

        if not isinstance(key, str):
            raise ReferencePreparationError(
                "Launch-node inventory contains an invalid object key"
            )

        if is_excluded(key):
            raise ReferencePreparationError(
                "Excluded worker reference object appears in inventory: "
                f"{key}"
            )


def validate_local_inventory(
        local_dir: Path,
        inventory: dict[str, object],
) -> None:
    """
    Ensure local files exactly match the expected S3 object names
    and sizes, recursively.

    Metadata files created by this utility are ignored.
    """

    validate_inventory_policy(
        inventory
    )

    expected = inventory_object_map(
        inventory
    )

    if not local_dir.is_dir():
        raise ReferencePreparationError(
            f"Reference directory is missing: {local_dir}"
        )

    actual: dict[str, int] = {}

    for path in local_dir.rglob("*"):

        relative = path.relative_to(
            local_dir
        ).as_posix()

        # Ignore only metadata owned by this utility.
        if relative in (
                ".reference_version",
                ".reference_manifest.json",
        ):
            continue

        if path.is_symlink():
            raise ReferencePreparationError(
                "Unexpected symlink in launch-node reference tree: "
                f"{relative}"
            )

        if path.is_dir():
            continue

        if not path.is_file():
            raise ReferencePreparationError(
                "Unexpected non-file entry in launch-node reference tree: "
                f"{relative}"
            )

        actual[relative] = path.stat().st_size

    expected_names = set(
        expected
    )

    actual_names = set(
        actual
    )

    missing = sorted(
        expected_names - actual_names
    )

    extra = sorted(
        actual_names - expected_names
    )

    wrong_size = sorted(
        name
        for name in (
                expected_names
                & actual_names
        )
        if expected[name][0] != actual[name]
    )

    if missing or extra or wrong_size:
        problems: list[str] = []

        if missing:
            problems.append(
                f"missing={missing[:50]}"
            )

        if extra:
            problems.append(
                f"extra={extra[:50]}"
            )

        if wrong_size:
            problems.append(
                f"wrong_size={wrong_size[:50]}"
            )

        raise ReferencePreparationError(
            "Local launch-node reference does not match S3 inventory: "
            + "; ".join(problems)
        )


def set_reference_permissions(
        reference_dir: Path,
) -> None:
    """
    Make the complete validated reference tree readable by pipeline
    processes.

    The reference root itself remains owned/managed by the launch
    bootstrap. Reference files are read-only to ordinary users.
    """

    reference_dir.chmod(0o755)

    for path in reference_dir.rglob("*"):

        if path.is_symlink():
            continue

        if path.is_dir():
            path.chmod(0o755)

        elif path.is_file():
            path.chmod(0o644)


def write_metadata(
        reference_dir: Path,
        inventory: dict[str, object],
        version: str,
) -> None:
    """Write launch-node reference metadata after validation."""

    (
            reference_dir
            / ".reference_version"
    ).write_text(
        version + "\n",
        encoding="utf-8",
        )

    (
            reference_dir
            / ".reference_manifest.json"
    ).write_text(
        json.dumps(
            inventory,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
        )


def stage_reference(
        local_root: Path,
        inventory: dict[str, object],
        version: str,
) -> None:
    """
    Download, validate, and atomically promote the entire
    launch-node reference cache.
    """

    local_root.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    staging_root = (
            local_root.parent
            / ".staging"
    )

    staging_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    stage_dir = Path(
        tempfile.mkdtemp(
            prefix="launch-node-",
            dir=staging_root,
        )
    )

    previous_dir = (
            local_root.parent
            / f".previous-launch-{os.getpid()}"
    )

    try:
        LOG.info(
            "Staging launch-node reference to %s",
            stage_dir,
        )

        command = [
            "s5cmd",
            "--numworkers",
            S5CMD_NUMWORKERS,
            "cp",
            "--concurrency",
            S5CMD_CONCURRENCY,
            "--part-size",
            S5CMD_PART_SIZE,
        ]

        # These are deliberately the ONLY exclusions.
        for prefix in EXCLUDED_PREFIXES:
            command.extend(
                [
                    "--exclude",
                    f"{prefix}*",
                ]
            )

        command.extend(
            [
                f"{S3_REFERENCE_ROOT}/*",
                f"{stage_dir}/",
            ]
        )

        run_command(
            command,
            stdout=subprocess.PIPE,
        )

        LOG.info(
            "Validating staged launch-node reference inventory"
        )

        validate_local_inventory(
            stage_dir,
            inventory,
        )

        LOG.info(
            "Writing launch-node reference metadata"
        )

        write_metadata(
            stage_dir,
            inventory,
            version,
        )

        set_reference_permissions(
            stage_dir
        )

        LOG.info(
            "Promoting validated launch-node reference to %s",
            local_root,
        )

        if previous_dir.exists():
            shutil.rmtree(
                previous_dir
            )

        if local_root.exists():
            os.replace(
                local_root,
                previous_dir,
            )

        try:
            os.replace(
                stage_dir,
                local_root,
            )

        except OSError:
            if (
                    previous_dir.exists()
                    and not local_root.exists()
            ):
                os.replace(
                    previous_dir,
                    local_root,
                )

            raise

        if previous_dir.exists():
            shutil.rmtree(
                previous_dir
            )

    except Exception:

        if stage_dir.exists():
            shutil.rmtree(
                stage_dir,
                ignore_errors=True,
            )

        if (
                previous_dir.exists()
                and not local_root.exists()
        ):
            try:
                os.replace(
                    previous_dir,
                    local_root,
                )

            except OSError:
                LOG.exception(
                    "Failed to restore previous launch-node reference cache"
                )

        raise


def local_cache_is_current(
        local_root: Path,
        inventory: dict[str, object],
        version: str,
) -> bool:
    """
    Return True only when the existing local cache exactly matches
    the current canonical inventory.
    """

    version_path = (
            local_root
            / ".reference_version"
    )

    manifest_path = (
            local_root
            / ".reference_manifest.json"
    )

    if (
            not version_path.is_file()
            or not manifest_path.is_file()
    ):
        return False

    try:
        recorded_version = (
            version_path
            .read_text(
                encoding="utf-8"
            )
            .strip()
        )

        recorded_inventory = load_inventory(
            manifest_path
        )

    except ReferencePreparationError:
        return False

    if recorded_version != version:
        return False

    if recorded_inventory != inventory:
        return False

    try:
        validate_local_inventory(
            local_root,
            inventory,
        )

    except ReferencePreparationError:
        return False

    return True


def prepare_reference(
        reference_set: str,
        local_root: Path,
        inventory_only: bool = False,
) -> str:
    """
    Build/validate the canonical launch-node reference cache and
    return its inventory fingerprint.
    """

    if reference_set != "phase2":
        raise ReferencePreparationError(
            f"Unsupported reference set: {reference_set}"
        )

    require_mountpoint(
        Path("/scratch")
    )

    validate_tools()

    with tempfile.TemporaryDirectory(
            prefix="seqtoid-launch-reference-"
    ) as tmp:

        inventory_path = (
                Path(tmp)
                / "reference-manifest.json"
        )

        objects = build_inventory(
            inventory_path
        )

        inventory = load_inventory(
            inventory_path
        )

        validate_inventory_policy(
            inventory
        )

        if objects != inventory["objects"]:
            raise ReferencePreparationError(
                "Reference inventory changed during build"
            )

        version = inventory_version(
            inventory_path
        )

        # Cross-check that the payload itself produces the same
        # deterministic fingerprint as the serialized inventory.
        payload_version = inventory_version_from_payload(
            inventory
        )

        if payload_version != version:
            raise ReferencePreparationError(
                "Reference inventory fingerprint is not deterministic"
            )

        LOG.info(
            "Canonical launch-node reference version: %s",
            version,
        )

        LOG.info(
            "Included launch-node object count: %d",
            len(objects),
        )

        LOG.info(
            "Excluded prefixes: %s",
            ", ".join(
                EXCLUDED_PREFIXES
            ),
        )

        if inventory_only:
            LOG.info(
                "Inventory-only mode requested; "
                "not modifying local reference cache"
            )

            return version

        if local_cache_is_current(
                local_root,
                inventory,
                version,
        ):
            LOG.info(
                "Existing launch-node reference cache is "
                "current and validated"
            )

            return version

        LOG.info(
            "Existing launch-node reference cache is absent, "
            "stale, or invalid"
        )

        stage_reference(
            local_root,
            inventory,
            version,
        )

    promoted_manifest_path = (
            local_root
            / ".reference_manifest.json"
    )

    promoted_version_path = (
            local_root
            / ".reference_version"
    )

    promoted_inventory = load_inventory(
        promoted_manifest_path
    )

    validate_inventory_policy(
        promoted_inventory
    )

    validate_local_inventory(
        local_root,
        promoted_inventory,
    )

    if not promoted_version_path.is_file():
        raise ReferencePreparationError(
            "Promoted launch-node reference has no version file"
        )

    promoted_version = (
        promoted_version_path
        .read_text(
            encoding="utf-8"
        )
        .strip()
    )

    if promoted_version != version:
        raise ReferencePreparationError(
            "Promoted launch-node reference version mismatch: "
            f"expected {version}, found {promoted_version}"
        )

    if promoted_inventory != inventory:
        raise ReferencePreparationError(
            "Promoted launch-node reference manifest does not "
            "match canonical S3 inventory"
        )

    LOG.info(
        "Launch-node reference installed and validated"
    )

    LOG.info(
        "  S3 source: %s",
        S3_REFERENCE_ROOT,
    )

    LOG.info(
        "  Local path: %s",
        local_root,
    )

    LOG.info(
        "  Version: %s",
        promoted_version,
    )

    LOG.info(
        "  Objects: %d",
        len(
            inventory_object_map(
                promoted_inventory
            )
        ),
    )

    return promoted_version


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare and validate the SeqToID launch-node reference cache."
        )
    )

    parser.add_argument(
        "--reference-set",
        default="phase2",
        choices=[
            "phase2",
        ],
        help="Reference set to prepare.",
    )

    parser.add_argument(
        "--local-root",
        type=Path,
        default=DEFAULT_LOCAL_ROOT,
        help=(
            "Root directory for the launch-node reference cache."
        ),
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
        ],
    )

    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help=(
            "Inspect and fingerprint the canonical launch-node "
            "S3 reference tree without downloading or modifying "
            "the local reference cache."
        ),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    logging.basicConfig(
        level=getattr(
            logging,
            args.log_level,
        ),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
    )

    try:
        version = prepare_reference(
            args.reference_set,
            args.local_root,
            args.inventory_only,
        )

    except ReferencePreparationError as exc:
        LOG.error(
            "Launch-node reference preparation failed: %s",
            exc,
        )

        return 1

    except Exception:
        LOG.exception(
            "Unexpected launch-node reference preparation failure"
        )

        return 1

    print(version)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())