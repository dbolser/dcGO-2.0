"""Provenance sidecars for the ``data/interim/protein2ipr_<species>.dat.gz`` extracts.

Two scripts write extracts to the same path from *different selection rules*:
``extract_human_interpro.py`` selects a species' proteins via its GOA file,
``scripts/extract_species_interpro.py`` via an id-mapping accession set. The
run manifest hashes whichever file is on disk, so a past run stays
identifiable — but nothing on disk said *which rule produced the current
bytes*, and the scripts would silently clobber each other (worm has both a
GOA file and a phenotype layer).

Each extract therefore gets a sidecar —
``protein2ipr_<species>.dat.gz.provenance.json`` — recording the selection
rule, the files that defined the accession set, and the match counts. Before
overwriting, a script must call :func:`ensure_overwrite_allowed`: replacing an
extract whose marker records a *different* selection rule is refused unless
forced, because the two universes answer different questions and the caller
should have to say which one they mean. Same-rule overwrites (refreshes) pass;
a pre-marker extract passes with a warning, since refusing would strand every
legacy file.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

from loguru import logger


class ProvenanceConflictError(RuntimeError):
    """An extract exists with a different selection rule; refuse to clobber."""


@dataclass(frozen=True)
class UniverseProvenance:
    """What produced a ``protein2ipr_<species>`` extract.

    Attributes:
        selection_rule: how the accession set was chosen — ``"goa"``
            (extract_human_interpro.py), ``"idmapping"`` or
            ``"accession_list"`` (scripts/extract_species_interpro.py).
        selection_sources: the file(s) that defined the accession set.
        interpro_source: the protein2ipr file that was filtered.
        n_accessions: size of the selecting accession set.
        n_matched_lines: protein2ipr lines written to the extract.
        tool: the script that wrote it.
        created: ISO-8601 UTC timestamp.
    """

    selection_rule: str
    selection_sources: tuple[str, ...]
    interpro_source: str
    n_accessions: int
    n_matched_lines: int
    tool: str
    created: str


def marker_path(extract_path: Path) -> Path:
    """The sidecar path for *extract_path* (``<name>.provenance.json``)."""
    extract_path = Path(extract_path)
    return extract_path.with_name(extract_path.name + ".provenance.json")


def read_marker(extract_path: Path) -> Optional[UniverseProvenance]:
    """Read the sidecar for *extract_path*, or ``None`` if there is none.

    An unreadable or structurally alien marker is treated as absent (with a
    warning) rather than fatal: it cannot *authorise* anything, and the
    overwrite decision then falls back to the unknown-provenance path.
    """
    path = marker_path(extract_path)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["selection_sources"] = tuple(payload["selection_sources"])
        return UniverseProvenance(**payload)
    except (ValueError, TypeError, KeyError) as exc:
        logger.warning(f"Unreadable provenance marker {path}: {exc}")
        return None


def write_marker(
    extract_path: Path,
    *,
    selection_rule: str,
    selection_sources: Sequence[Path],
    interpro_source: Path,
    n_accessions: int,
    n_matched_lines: int,
    tool: str,
    created: Optional[str] = None,
) -> Path:
    """Write the sidecar for a freshly written extract; returns its path."""
    provenance = UniverseProvenance(
        selection_rule=selection_rule,
        selection_sources=tuple(str(source) for source in selection_sources),
        interpro_source=str(interpro_source),
        n_accessions=n_accessions,
        n_matched_lines=n_matched_lines,
        tool=tool,
        created=created or datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )
    path = marker_path(extract_path)
    path.write_text(json.dumps(asdict(provenance), indent=2) + "\n", encoding="utf-8")
    logger.info(f"  Provenance marker written: {path}")
    return path


def ensure_overwrite_allowed(
    extract_path: Path, selection_rule: str, *, force: bool = False
) -> None:
    """Refuse to replace an extract built under a different selection rule.

    Raises:
        ProvenanceConflictError: the extract exists, its marker records a
            different ``selection_rule``, and ``force`` is not set.
    """
    extract_path = Path(extract_path)
    if not extract_path.exists():
        return
    existing = read_marker(extract_path)
    if existing is None:
        logger.warning(
            f"{extract_path} exists with no provenance marker (pre-marker "
            f"extract); replacing it with a {selection_rule!r}-selected one."
        )
        return
    if existing.selection_rule != selection_rule and not force:
        raise ProvenanceConflictError(
            f"{extract_path} was built with selection rule "
            f"{existing.selection_rule!r} (from "
            f"{', '.join(existing.selection_sources)}); refusing to replace it "
            f"with a {selection_rule!r}-selected extract. Re-run with --force "
            "to overwrite, or write to a different --output path."
        )
    if existing.selection_rule != selection_rule:
        logger.warning(
            f"--force: replacing {existing.selection_rule!r}-selected "
            f"{extract_path} with a {selection_rule!r}-selected extract"
        )
