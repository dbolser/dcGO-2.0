#!/usr/bin/env python3
"""Quantify how far one GO term dominates human experimental MF annotation.

Why this exists: the §2 read-out explains MF as the one aspect where the naive
frequency baseline beats dcGO at IC>=0, on the grounds that MF "truth" is
dominated by `GO:0005515 protein binding`. That claim carried a single number —
84.6% — that appeared only in prose, with no script behind it, which made it
impossible to check or to refresh against a newer GOA release.

It also turns out to be the *least* relevant of the three ways to measure this,
because it counts annotation lines and the benchmark scores proteins:

* **Annotation share** — what fraction of experimental MF annotation *lines*
  are protein binding. This is the ~84% figure. It over-weights proteins that
  many papers have reported the same interaction for.
* **Pair share** — what fraction of distinct (protein, MF term) pairs. Much
  lower, because the duplicate lines collapse.
* **Protein coverage** — what fraction of proteins with any experimental MF
  term have protein binding, and what fraction have *nothing else*. This is the
  one that explains the benchmark behaviour: a baseline that predicts the most
  frequent term for everything is right about that term for most proteins, and
  for the "only protein binding" proteins it is right about their whole truth.

Usage
-----
    uv run python validation/protein_binding_dominance.py
    uv run python validation/protein_binding_dominance.py \
        --gaf data/raw/goa_archive/goa_human.gaf.205.gz --aspect F
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, Set

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.goa_parser import EXPERIMENTAL_EVIDENCE  # noqa: E402

#: `protein binding` — the near-universal, near-zero-information MF term.
PROTEIN_BINDING = "GO:0005515"


def read_protein_terms(
    gaf_path: Path, aspect: str, evidence_codes: Set[str]
) -> tuple[Dict[str, Set[str]], int]:
    """Return ``{protein: {term}}`` for one aspect, plus the raw line count.

    The line count is kept separately because the annotation-share statistic
    needs the duplicates that the protein map collapses.
    """
    protein_terms: Dict[str, Set[str]] = defaultdict(set)
    n_lines = 0
    opener = gzip.open if gaf_path.suffix == ".gz" else open
    with opener(gaf_path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15:
                continue
            if fields[8] != aspect or fields[6] not in evidence_codes:
                continue
            n_lines += 1
            protein_terms[fields[1]].add(fields[4])
    return protein_terms, n_lines


def summarise(
    protein_terms: Dict[str, Set[str]], n_lines: int, term: str, gaf_path: Path
) -> list[tuple[str, str, str]]:
    """Build the ``(statistic, value, share)`` rows for the report."""
    n_proteins = len(protein_terms)
    if not n_proteins:
        return []
    n_pairs = sum(len(terms) for terms in protein_terms.values())
    with_term = sum(1 for terms in protein_terms.values() if term in terms)
    only_term = sum(1 for terms in protein_terms.values() if terms == {term})
    # The line-level share is injected by the caller: `protein_terms` has
    # already collapsed the duplicate lines it would need.
    return [
        ("proteins with any annotation", f"{n_proteins:,}", ""),
        ("distinct (protein, term) pairs", f"{n_pairs:,}", ""),
        ("annotation lines", f"{n_lines:,}", ""),
        (
            f"proteins carrying {term}",
            f"{with_term:,}",
            f"{100 * with_term / n_proteins:.1f}% of proteins",
        ),
        (
            f"proteins carrying ONLY {term}",
            f"{only_term:,}",
            f"{100 * only_term / n_proteins:.1f}% of proteins",
        ),
        (
            "median terms per protein",
            f"{median(len(t) for t in protein_terms.values()):.0f}",
            "",
        ),
        ("source", gaf_path.name, ""),
    ]


def main() -> int:  # pragma: no cover - I/O wiring
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gaf",
        type=Path,
        default=Path("data/raw/goa_annotations/goa_human.gaf.gz"),
        help="GOA GAF file (default: the current human release)",
    )
    parser.add_argument(
        "--aspect", default="F", choices=["F", "P", "C"], help="GO aspect"
    )
    parser.add_argument("--term", default=PROTEIN_BINDING, help="Term to measure")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation/protein_binding_dominance.tsv"),
    )
    args = parser.parse_args()

    if not args.gaf.exists():
        print(f"Missing GAF: {args.gaf}", file=sys.stderr)
        return 1

    protein_terms, n_lines = read_protein_terms(
        args.gaf, args.aspect, EXPERIMENTAL_EVIDENCE
    )
    if not protein_terms:
        print(
            f"No {args.aspect} annotations with experimental evidence", file=sys.stderr
        )
        return 1

    # The line-level share needs the duplicates, so count them in a second pass
    # over the same file rather than keeping every line in memory.
    n_lines_term = 0
    opener = gzip.open if args.gaf.suffix == ".gz" else open
    with opener(args.gaf, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15:
                continue
            if (
                fields[8] == args.aspect
                and fields[6] in EXPERIMENTAL_EVIDENCE
                and fields[4] == args.term
            ):
                n_lines_term += 1

    rows = summarise(protein_terms, n_lines, args.term, args.gaf)
    rows.insert(
        3,
        (
            f"annotation lines carrying {args.term}",
            f"{n_lines_term:,}",
            f"{100 * n_lines_term / n_lines:.1f}% of lines",
        ),
    )

    width = max(len(name) for name, _v, _s in rows)
    print(f"\nAspect {args.aspect}, experimental evidence only\n")
    for name, value, share in rows:
        print(f"  {name:<{width}}  {value:>12}  {share}")
    print()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("statistic\tvalue\tshare\n")
        for name, value, share in rows:
            handle.write(f"{name}\t{value}\t{share}\n")
    print(f"Wrote {args.output}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
