"""Machine-readable provenance manifests for dcGO analysis runs.

The manifest records exact files rather than relying on mutable upstream URLs:
almost every source this pipeline consumes is published at a ``current_release``
address that is overwritten in place. Embedded release headers (``!gaf-version``,
``data-version``, ``Release:`` …) are useful *labels*; the SHA-256 is what
actually identifies the bytes that were analysed.

A manifest is written twice: once as ``"status": "running"`` before the
expensive stages, and once as ``"status": "completed"`` with output hashes and
result counts when the run finishes. A run that crashes therefore leaves a
manifest that is visibly unfinished rather than a stale "completed" record.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

MANIFEST_SCHEMA_VERSION = "1.0"
DEFAULT_LOCK_FILE = "uv.lock"
_HASH_CHUNK_SIZE = 8 * 1024 * 1024
_MAX_HEADER_LINES = 200

#: Header fields worth recording, normalised to ``lower_snake_case``. A
#: whitelist (rather than "every ``key: value`` line") keeps incidental colons
#: in prose headers — ``Email: enzyme@expasy.org`` — out of the manifest.
_RELEASE_HEADER_FIELDS = frozenset(
    {
        "data_version",
        "date",
        "date_generated",
        "format_version",
        "gaf_version",
        "generated_by",
        "go_version",
        "release",
    }
)

#: Expasy ``enzyme.dat`` states its release without a colon:
#: ``CC   Release of 10-Jun-2026``.
_RELEASE_OF_PATTERN = re.compile(r"^Release of\s+(?P<release>.+?)\s*$")


def manifest_filename(label: str) -> str:
    """Manifest file name for one ontology's run.

    Runs for different ontologies routinely share an output directory (they
    write ``domain_<label>_associations_*.tsv`` side by side), so the manifest
    is labelled the same way instead of clobbering its predecessor.
    """
    return f"run_manifest_{label}.json"


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of ``path`` without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def embedded_release_metadata(path: Path) -> dict[str, str]:
    """Extract stable release/header fields from supported text inputs.

    Covers the header conventions of the formats this pipeline reads: GAF
    comment lines (``!gaf-version: 2.2``), OBO stanzas (``data-version:``),
    UniProt controlled-vocabulary documents (``Release:     2026_02 of …``) and
    Expasy flat files (``CC   Release of 10-Jun-2026``). Returns ``{}`` when the
    format carries no recognised release identifier — which means "not present",
    never "current".
    """
    path = Path(path)
    metadata: dict[str, str] = {}
    try:
        with _open_text(path) as handle:
            for index, raw_line in enumerate(handle):
                if index >= _MAX_HEADER_LINES:
                    break
                line = raw_line.strip()
                if line.startswith("["):
                    break
                if line.startswith("!"):
                    line = line[1:].strip()
                elif line.startswith("CC "):
                    line = line[3:].strip()
                if match := _RELEASE_OF_PATTERN.match(line):
                    metadata.setdefault("release", match.group("release"))
                    continue
                if ":" not in line:
                    continue
                key, value = (part.strip() for part in line.split(":", 1))
                normalized = key.lower().replace("-", "_").replace(" ", "_")
                if normalized in _RELEASE_HEADER_FIELDS and value:
                    metadata.setdefault(normalized, value)
    except (OSError, UnicodeError):
        return {}
    return metadata


def describe_file(
    path: Path,
    *,
    role: str,
    source_url: str | None = None,
    derived_from: str | None = None,
) -> dict[str, Any]:
    """Describe and hash one concrete pipeline input or output."""
    path = Path(path)
    stat = path.stat()
    record: dict[str, Any] = {
        "role": role,
        "path": str(path),
        "size_bytes": stat.st_size,
        "modified_at": datetime.fromtimestamp(
            stat.st_mtime, tz=timezone.utc
        ).isoformat(),
        "sha256": sha256_file(path),
    }
    if source_url:
        record["source_url"] = source_url
    if derived_from:
        record["derived_from"] = derived_from
    release = embedded_release_metadata(path)
    if release:
        record["release_metadata"] = release
    return record


def git_metadata(repository: Path) -> dict[str, Any]:
    """Return commit, branch, and working-tree state, or an unavailable marker.

    ``dirty`` covers *tracked* modifications only — the ones that mean the code
    that ran is not the code at ``commit``. Untracked files are counted
    separately: this repository accumulates scratch results and data symlinks
    that no run depends on, and folding them into ``dirty`` would make the flag
    permanently true and therefore worthless.
    """

    def run(*args: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(repository), *args],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    try:
        untracked = run("ls-files", "--others", "--exclude-standard")
        return {
            "available": True,
            "commit": run("rev-parse", "HEAD"),
            "branch": run("branch", "--show-current") or None,
            "dirty": bool(run("status", "--porcelain", "--untracked-files=no")),
            "untracked_files": len(untracked.splitlines()) if untracked else 0,
        }
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {"available": False}


def dependency_lock_metadata(
    repository: Path, *, filename: str = DEFAULT_LOCK_FILE
) -> dict[str, Any]:
    """Hash the resolved dependency set (``uv.lock``).

    The installed-version list is environment state; the lock file is the
    reproducible *instruction* for rebuilding it, so its digest is what a
    third party needs to recreate the environment (``uv sync --frozen``).
    """
    lock = Path(repository) / filename
    if not lock.is_file():
        return {"available": False, "path": filename}
    return {
        "available": True,
        "path": filename,
        "size_bytes": lock.stat().st_size,
        "sha256": sha256_file(lock),
    }


def software_metadata(repository: Path | None = None) -> dict[str, Any]:
    """Return runtime, dependency-lock, and installed dcGO package information."""
    try:
        version = importlib.metadata.version("dcgo-pipeline")
    except importlib.metadata.PackageNotFoundError:
        version = None
    metadata: dict[str, Any] = {
        "dcgo_pipeline_version": version,
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "executable": sys.executable,
    }
    if repository is not None:
        metadata["dependency_lock"] = dependency_lock_metadata(repository)
    return metadata


def json_parameters(parameters: Mapping[str, Any]) -> dict[str, Any]:
    """Convert argparse-style values into portable JSON values.

    ``vars(args)`` is full of :class:`~pathlib.Path` objects, and new options
    are added without anyone remembering this file exists, so anything that is
    not a JSON scalar or container is stringified rather than allowed to blow
    up the run it is meant to document.
    """
    return {key: _json_value(value) for key, value in parameters.items()}


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_value(item) for item in value]
    return str(value)


class RunManifest:
    """Create and atomically update the provenance record for one pipeline run."""

    def __init__(
        self,
        path: Path,
        *,
        repository: Path,
        parameters: Mapping[str, Any],
        inputs: Iterable[Mapping[str, Any]],
        analysis: Mapping[str, Any] | None = None,
        command: Iterable[str] | None = None,
    ) -> None:
        self.path = Path(path)
        self.data: dict[str, Any] = {
            "schema_version": MANIFEST_SCHEMA_VERSION,
            "status": "running",
            "started_at": utc_now(),
            "command": list(command if command is not None else sys.argv),
            "working_directory": os.getcwd(),
            "parameters": json_parameters(parameters),
            "analysis": _json_value(dict(analysis)) if analysis else {},
            "software": software_metadata(repository),
            "git": git_metadata(repository),
            "inputs": [dict(item) for item in inputs],
            "outputs": [],
        }
        self.write()

    def write(self) -> None:
        """Atomically write the current manifest state."""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(self.data, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(self.path)

    def complete(
        self,
        *,
        outputs: Iterable[Mapping[str, Any]],
        summary: Mapping[str, Any],
    ) -> None:
        """Finalize a successful run with output identities and result counts."""
        self.data.update(
            {
                "status": "completed",
                "completed_at": utc_now(),
                "outputs": [dict(item) for item in outputs],
                "summary": _json_value(dict(summary)),
            }
        )
        self.write()
