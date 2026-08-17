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
applied after BH in the pipeline, and the column is written at full precision
(``%.10g``) so the comparison here sees the same values the in-run floor
compared — filtering the unfloored table on that column is *exactly* the
floored run's reported set, floor boundaries included. One pipeline run
yields the whole sweep:

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
from typing import Callable, Iterable

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
) -> tuple[list[dict[str, str]], list[str]]:
    """Rows and header of a runner-written association TSV.

    The header comes back separately (``DictReader.fieldnames``) so callers
    can validate columns on an empty-but-valid table — a header-only TSV is a
    legitimate zero-association run, not a schema error.
    """
    import csv

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if term_column not in fieldnames:
        raise SystemExit(f"{path}: no {term_column!r} column in header {fieldnames}")
    return rows, fieldnames


def main(argv: list[str] | None = None) -> int:
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
    args = parser.parse_args(argv)

    from src.ontology_processor import OntologyProcessor

    rows, fieldnames = load_association_rows(args.associations, args.term_column)
    if args.domain_type != "any":
        if "domain_type" not in fieldnames:
            raise SystemExit(
                f"{args.associations}: --domain-type needs a domain_type "
                f"column; header is {fieldnames}"
            )
        # The runner writes DomainType values; anything else (a missing cell,
        # a foreign table) must fail loudly rather than land in a family by
        # accident of a truth-table trick.
        known_types = {"single", "supra_pair", "supra_triple"}
        unrecognised = {r["domain_type"] for r in rows} - known_types
        if unrecognised:
            raise SystemExit(
                f"{args.associations}: unrecognised domain_type values "
                f"{sorted(str(v) for v in unrecognised)}; "
                f"expected {sorted(known_types)}"
            )
        wanted = (
            {"single"}
            if args.domain_type == "single"
            else {"supra_pair", "supra_triple"}
        )
        rows = [r for r in rows if r["domain_type"] in wanted]
    # Validated on the header, not rows[0]: an empty-but-valid TSV (a run with
    # zero significant associations) sweeps to zero-count rows, it does not
    # crash or masquerade as a schema problem.
    if any(floor > 0 for floor in args.min_ic) and "ic" not in fieldnames:
        raise SystemExit(
            f"{args.associations}: a --min-ic sweep needs the ic column "
            "(written by runs since the --min-ic feature)"
        )

    processor = OntologyProcessor(args.obo)

    def get_ancestors(term: str) -> set[str]:
        if term not in processor.go_graph:
            return set()
        return processor.get_ancestors(term)

    # A table measured against the wrong ontology looks *perfect* — every term
    # contributes an empty closure, so 0.0% on-chain and no roots. Refuse the
    # vacuous case outright and flag partial mismatches.
    table_terms = {r[args.term_column] for r in rows}
    if table_terms:
        known_terms = sum(1 for term in table_terms if term in processor.go_graph)
        if known_terms == 0:
            raise SystemExit(
                f"{args.associations}: none of its {len(table_terms):,} distinct "
                f"terms are in {args.obo} — wrong ontology for this table; every "
                "metric would be a vacuous zero"
            )
        if known_terms / len(table_terms) < 0.5:
            print(
                f"WARNING: only {known_terms:,}/{len(table_terms):,} distinct "
                f"terms are in {args.obo}; the metrics treat the rest as having "
                "no ancestors",
                file=sys.stderr,
            )

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
