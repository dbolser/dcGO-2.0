"""Canonical readers for dcGO association tables used by validation tools."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

Pair = tuple[str, str]


def load_associations(
    path: Path,
    *,
    term_column: str = "go_term",
    required_columns: Iterable[str] = (),
) -> pd.DataFrame:
    """Load an association TSV and validate its named-column contract.

    ``term_column`` supports the ontology-specific schemas written by the runner,
    such as ``go_term``, ``ec_term``, and ``doid_term``.
    """
    frame = pd.read_csv(path, sep="\t")
    required = {"domain", term_column, *required_columns}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(
            f"Association file {path} is missing required columns: {', '.join(missing)}"
        )
    return frame


def association_pairs(
    frame: pd.DataFrame, *, term_column: str = "go_term"
) -> set[Pair]:
    """Return unique ``(domain, term)`` pairs from a validated table."""
    return set(zip(frame["domain"], frame[term_column]))
