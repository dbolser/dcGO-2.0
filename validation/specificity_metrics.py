#!/usr/bin/env python3
"""Specificity metrics for an association table (VALIDATION_PLAN next-steps §2).

The relative-inference work is judged on whether reported terms sit at the most
specific supportable level, summarised by three numbers per configuration:

  * **significant** — association rows evaluated.
  * **mean #ancestors** — mean ancestor-closure size of the reported terms; a
    generic term near the root has few, so *lower* is more general, not better:
    read it together with the chain share.
  * **on a chain** — share of associations whose domain is *also* significant
    for an ancestor of the same term. A cascade of one signal reported at every
    level up the DAG shows up here directly.

It also reports which of the DAG roots appear at all — the vacuous associations
the ``--min-ic`` floor exists to remove (a root has no parents, so the relative
inference can never test it).

The corresponding numbers in ``VALIDATION_PLAN.md`` (28.6% / 55.2% / 82.4%)
were computed ad hoc; this module makes them reproducible, including the
``--min-ic`` sweep over the pipeline's exported ``ic`` column. The floor is
applied after BH in the pipeline, so filtering the unfloored table on that
column is *exactly* the floored run's reported set — one pipeline run yields
the whole sweep:

    uv run python validation/specificity_metrics.py \\
        --associations results/domain_go_associations_significant.tsv \\
        --obo data/raw/go_ontology/go-basic.obo \\
        --domain-type single --min-ic 0 1 2 3 5
"""

from __future__ import annotations

import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from validation.temporal_benchmark import ASPECT_ROOTS  # noqa: E402

GetAncestors = Callable[[str], set[str]]


@dataclass(frozen=True)
class SpecificityMetrics:
    """The per-configuration row of the specificity table."""

    n_associations: int
    mean_ancestors: float
    on_chain_share: float
    roots_present: tuple[str, ...]


def specificity_metrics(
    pairs: Iterable[tuple[str, str]],
    get_ancestors: GetAncestors,
    roots: frozenset[str] = ASPECT_ROOTS,
) -> SpecificityMetrics:
    """Compute the specificity summary over ``(domain, term)`` associations.

    "On a chain" counts a row when its domain is also associated with an
    *ancestor* of its term (the descendant end of an ancestor–descendant
    chain). Terms the ontology does not contain contribute an empty closure —
    zero ancestors, never on a chain — rather than failing.
    """
    by_domain: dict[str, set[str]] = defaultdict(set)
    for domain, term in pairs:
        by_domain[domain].add(term)

    n = 0
    ancestor_total = 0
    on_chain = 0
    roots_present: set[str] = set()
    for terms in by_domain.values():
        for term in terms:
            ancestors = get_ancestors(term)
            n += 1
            ancestor_total += len(ancestors)
            if ancestors & terms:
                on_chain += 1
            if term in roots:
                roots_present.add(term)
    return SpecificityMetrics(
        n_associations=n,
        mean_ancestors=ancestor_total / n if n else 0.0,
        on_chain_share=on_chain / n if n else 0.0,
        roots_present=tuple(sorted(roots_present)),
    )


def load_association_rows(
    path: Path, term_column: str = "go_term"
) -> list[Mapping[str, str]]:
    """Rows of a runner-written association TSV, keyed by header name."""
    import csv

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if rows and term_column not in rows[0]:
        raise SystemExit(f"{path}: no {term_column!r} column in header {list(rows[0])}")
    return rows


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Specificity metrics (mean #ancestors, ancestor-chain share, "
        "root presence) for a domain-association table, optionally swept over "
        "information-content floors."
    )
    parser.add_argument(
        "--associations",
        type=Path,
        required=True,
        help="A domain_<ontology>_associations_significant.tsv from the runner",
    )
    parser.add_argument(
        "--obo",
        type=Path,
        default=Path("data/raw/go_ontology/go-basic.obo"),
        help="Ontology OBO; ancestors come from the pipeline's own "
        "OntologyProcessor, so the edge policy (is_a/part_of) matches the runs "
        "being measured (default: data/raw/go_ontology/go-basic.obo)",
    )
    parser.add_argument(
        "--term-column",
        default="go_term",
        help="Term column name in the TSV (default: go_term)",
    )
    parser.add_argument(
        "--domain-type",
        choices=["any", "single", "supra"],
        default="any",
        help="Restrict to a domain family before measuring. BH corrects the "
        "families separately, so 'single' reproduces a --disable-supra-domains "
        "run's rows exactly (default: any)",
    )
    parser.add_argument(
        "--min-ic",
        type=float,
        nargs="+",
        default=[0.0],
        metavar="FLOAT",
        help="IC floors to sweep, read from the table's ic column. Because the "
        "pipeline applies --min-ic after BH, each floored view here is exactly "
        "that floored run's reported set (default: 0, no floor)",
    )
    args = parser.parse_args()

    from src.ontology_processor import OntologyProcessor

    rows = load_association_rows(args.associations, args.term_column)
    if args.domain_type != "any":
        rows = [
            r
            for r in rows
            if (r.get("domain_type") == "single") == (args.domain_type == "single")
        ]
    if any(floor > 0 for floor in args.min_ic) and (not rows or "ic" not in rows[0]):
        raise SystemExit(
            f"{args.associations}: a --min-ic sweep needs the ic column "
            "(written by runs since the --min-ic feature)"
        )

    processor = OntologyProcessor(args.obo)

    def get_ancestors(term: str) -> set[str]:
        if term not in processor.go_graph:
            return set()
        return processor.get_ancestors(term)

    print("min_ic\tsignificant\tmean_ancestors\ton_chain_share\troots_present")
    for floor in args.min_ic:
        kept = rows if floor <= 0 else [r for r in rows if float(r["ic"]) >= floor]
        m = specificity_metrics(
            ((r["domain"], r[args.term_column]) for r in kept), get_ancestors
        )
        roots = ",".join(m.roots_present) if m.roots_present else "-"
        print(
            f"{floor:g}\t{m.n_associations}\t{m.mean_ancestors:.1f}\t"
            f"{m.on_chain_share:.1%}\t{roots}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
