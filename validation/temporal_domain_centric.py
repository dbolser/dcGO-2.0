#!/usr/bin/env python3
"""Association-level temporal metric for the §4 ablation rungs.

The ablation (``validation/ablation.py``) judges each rung by protein-centric
F_max/AUPRC, which mostly tests the domain-to-protein *transfer* step and how
much of the cohort has a domain at all. A rung's actual product is a table of
(domain, GO) associations, so this scores that table directly against what
curators **added** to InterPro2GO between the training snapshot (t0, InterPro
85.0, 2021-04) and the current one (t1):

* **added pairs** — t1 InterPro2GO propagated up the GO DAG, minus t0
  propagated, restricted to domains present in *both* maps (new InterPro
  entries are not a prediction target) and in the rung's domain universe.
* **predicted pairs** — the rung's single-domain associations, direct plus
  ancestors, on the same domains.
* **precision** = |pred ∩ added| / |pred on domains that gained ≥1 pair|
  (domains curators did not touch cannot score a hit, so they are not
  counted against the rung). ``precision_novel`` additionally drops from the
  denominator the predictions that were already curated at t0 — those are
  not candidates for "added" either way.
* **recall** = |pred ∩ added| / |added|.
* **precision_at_budget** — the same precision over the top-``budget`` direct
  single-domain associations by q-value, so rungs of very different size are
  compared at equal cost. The driver sets the budget to the smallest rung's
  eligible association count.

Uncertainty is a **paired bootstrap over domains**: every replicate resamples
the eligible domains with replacement and re-scores every rung on that one
resample, so a rung-vs-reference delta has a legitimate interval
(``validation/resampling.py`` does the same over proteins).

Rungs keyed by a member-database signature (``--domain-key ssf`` / ``pfam``:
``SSFnnnnn`` / ``PFnnnnn`` ids) are mapped to their InterPro entry through the
protein2ipr file they were trained on (column 4 → column 2; a signature belongs
to exactly one entry) before scoring; unmapped signatures are counted.

Run on one or more rung tables::

    uv run python validation/temporal_domain_centric.py \\
        --rung supra=results/ablation-t0interpro/ipr_manual/supra/domain_go_associations_significant.tsv \\
        --rung full=results/ablation-t0interpro/ipr_manual/full/domain_go_annotations_propagated.tsv \\
        --interpro2go-t0 data/raw/interpro_archive/85.0/interpro2go \\
        --interpro2go-t1 data/raw/interpro2go/interpro2go \\
        --interpro data/interim/protein2ipr_human_t0_ipr85.dat.gz
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

Pair = tuple[str, str]
GetAncestors = Callable[[str], Iterable[str]]

#: a direct association: (domain, term, q-value)
Scored = tuple[str, str, float]

SIGNATURE_PREFIXES = {"ssf": "SSF", "pfam": "PF"}
METRICS = (
    "precision",
    "precision_novel",
    "recall",
    "precision_at_budget",
    "recall_at_budget",
)
COUNTS = (
    "n_eligible_domains",
    "n_gained_domains",
    "n_added_pairs",
    "n_direct_eligible",
    "n_pred_pairs",
    "n_pred_on_gained",
    "n_hit",
    "n_hit_at_budget",
    "n_unmapped_signatures",
)


# --------------------------------------------------------------------------- #
# Inputs                                                                       #
# --------------------------------------------------------------------------- #
def read_interpro2go(path: Path) -> set[Pair]:
    """``(IPR, GO)`` pairs of an InterPro2GO file (``InterPro:IPR… > … ; GO:…``)."""
    pairs: set[Pair] = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!") or ">" not in line:
                continue
            left, right = line.rstrip("\n").split(">", 1)
            domain = left.split()[0].removeprefix("InterPro:")
            term = right.rsplit(";", 1)[-1].strip()
            if domain and term.startswith("GO:"):
                pairs.add((domain, term))
    return pairs


def propagate_pairs(pairs: Iterable[Pair], get_ancestors: GetAncestors) -> set[Pair]:
    """Each ``(domain, term)`` plus ``(domain, ancestor)`` for every ancestor."""
    out: set[Pair] = set()
    for domain, term in pairs:
        out.add((domain, term))
        out.update((domain, ancestor) for ancestor in get_ancestors(term))
    return out


@dataclass(frozen=True)
class CuratedDelta:
    """What InterPro2GO curation added between t0 and t1 (both propagated)."""

    t0: frozenset[Pair]
    t1: frozenset[Pair]
    #: domains present in both snapshots: the only ones with an "added" notion
    shared_domains: frozenset[str]

    @property
    def added(self) -> frozenset[Pair]:
        return frozenset(p for p in self.t1 - self.t0 if p[0] in self.shared_domains)


def curated_delta(
    t0_pairs: Iterable[Pair], t1_pairs: Iterable[Pair], get_ancestors: GetAncestors
) -> CuratedDelta:
    t0 = propagate_pairs(t0_pairs, get_ancestors)
    t1 = propagate_pairs(t1_pairs, get_ancestors)
    shared = {d for d, _ in t0} & {d for d, _ in t1}
    return CuratedDelta(frozenset(t0), frozenset(t1), frozenset(shared))


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def read_protein2ipr_universe(
    path: Path, domain_key: str = "interpro"
) -> tuple[set[str], dict[str, str]]:
    """The InterPro entries a run keyed by ``domain_key`` could have tested.

    Returns ``(entries, signature -> entry)``. For ``interpro`` the universe is
    every entry in the file; for a signature key it is the entries carrying a
    signature with that prefix, and the map lets the run's signature ids be
    re-keyed onto entries. A signature seen under two entries is dropped from
    the map (counted as unmapped downstream) rather than guessed.
    """
    prefix = SIGNATURE_PREFIXES.get(domain_key)
    entries: set[str] = set()
    sig_to_entry: dict[str, str] = {}
    conflicts: set[str] = set()
    with _open_text(path) as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            entry, signature = cols[1], cols[3]
            if prefix is None:
                entries.add(entry)
                continue
            if not signature.startswith(prefix):
                continue
            entries.add(entry)
            previous = sig_to_entry.setdefault(signature, entry)
            if previous != entry:
                conflicts.add(signature)
    for signature in conflicts:
        del sig_to_entry[signature]
    return entries, sig_to_entry


def is_single_domain(domain: str) -> bool:
    """Supra-domain ids are comma-joined constituent ids."""
    return "," not in domain


def load_single_associations(path: Path) -> list[Scored]:
    """Direct single-domain ``(domain, term, q)`` rows of a rung's table.

    Accepts the significant-association table (``adj_p_value``) and the
    True-Path-propagated one (``q_value``; its ancestor rows are kept — they
    are what that rung reports).
    """
    frame = pd.read_csv(path, sep="\t")
    q_column = "adj_p_value" if "adj_p_value" in frame.columns else "q_value"
    missing = {"domain", "go_term", q_column} - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    frame = frame[frame["domain"].map(is_single_domain)]
    return [
        (str(d), str(t), float(q))
        for d, t, q in zip(frame["domain"], frame["go_term"], frame[q_column])
    ]


@dataclass(frozen=True)
class RungTable:
    """A rung's direct single-domain associations, keyed by InterPro entry."""

    name: str
    direct: tuple[Scored, ...]
    n_unmapped_signatures: int = 0


def map_signatures(
    rows: Iterable[Scored], sig_to_entry: Mapping[str, str]
) -> tuple[list[Scored], int]:
    """Re-key signature-keyed associations onto InterPro entries.

    Signatures of one entry pool: a pair predicted by two of them keeps the
    smaller q. Returns ``(rows, n_unmapped_signatures)``.
    """
    best: dict[Pair, float] = {}
    unmapped: set[str] = set()
    for signature, term, q in rows:
        entry = sig_to_entry.get(signature)
        if entry is None:
            unmapped.add(signature)
            continue
        key = (entry, term)
        if q < best.get(key, math.inf):
            best[key] = q
    return [(d, t, q) for (d, t), q in best.items()], len(unmapped)


def load_rung(
    name: str,
    path: Path,
    sig_to_entry: Mapping[str, str] | None = None,
) -> RungTable:
    rows = load_single_associations(path)
    n_unmapped = 0
    if sig_to_entry is not None:
        rows, n_unmapped = map_signatures(rows, sig_to_entry)
    return RungTable(name, tuple(rows), n_unmapped)


# --------------------------------------------------------------------------- #
# Scoring                                                                      #
# --------------------------------------------------------------------------- #
def top_by_q(rows: Iterable[Scored], budget: int | None) -> list[Scored]:
    """The ``budget`` most significant rows (ties broken by id, for determinism)."""
    ordered = sorted(rows, key=lambda r: (r[2], r[0], r[1]))
    return ordered if budget is None else ordered[:budget]


@dataclass
class DomainCounts:
    """Per-domain counts for one rung, over the eligible domains in order.

    Every metric is a ratio of sums of these, so a bootstrap replicate over
    domains is a weighted sum with the resample's multiplicities.
    """

    domains: tuple[str, ...]
    added: np.ndarray
    hit: np.ndarray
    pred_on_gained: np.ndarray
    known_on_gained: np.ndarray
    hit_b: np.ndarray
    pred_on_gained_b: np.ndarray
    known_on_gained_b: np.ndarray
    n_pred_pairs: int = 0
    n_direct_eligible: int = 0
    n_unmapped_signatures: int = 0
    extras: dict[str, float] = field(default_factory=dict)

    def metrics(self, weights: np.ndarray | None = None) -> dict[str, float]:
        w = np.ones(len(self.domains)) if weights is None else weights

        def ratio(num: np.ndarray, den: np.ndarray) -> float:
            d = float(w @ den)
            return float(w @ num) / d if d > 0 else float("nan")

        return {
            "precision": ratio(self.hit, self.pred_on_gained),
            "precision_novel": ratio(
                self.hit, self.pred_on_gained - self.known_on_gained
            ),
            "recall": ratio(self.hit, self.added),
            "precision_at_budget": ratio(self.hit_b, self.pred_on_gained_b),
            "recall_at_budget": ratio(self.hit_b, self.added),
        }

    def counts(self) -> dict[str, float]:
        return {
            "n_eligible_domains": len(self.domains),
            "n_gained_domains": int((self.added > 0).sum()),
            "n_added_pairs": int(self.added.sum()),
            "n_direct_eligible": self.n_direct_eligible,
            "n_pred_pairs": self.n_pred_pairs,
            "n_pred_on_gained": int(self.pred_on_gained.sum()),
            "n_hit": int(self.hit.sum()),
            "n_hit_at_budget": int(self.hit_b.sum()),
            "n_unmapped_signatures": self.n_unmapped_signatures,
        }


def domain_counts(
    rung: RungTable,
    delta: CuratedDelta,
    eligible: Sequence[str],
    get_ancestors: GetAncestors,
    budget: int | None,
) -> DomainCounts:
    """Score one rung's table against the curated delta, per eligible domain."""
    index = {d: i for i, d in enumerate(eligible)}
    added_pairs = {p for p in delta.added if p[0] in index}
    gained = {d for d, _ in added_pairs}
    direct = [r for r in rung.direct if r[0] in index]

    def tally(rows: Iterable[Scored]) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        pred = propagate_pairs(((d, t) for d, t, _ in rows), get_ancestors)
        hit = np.zeros(len(eligible))
        on_gained = np.zeros(len(eligible))
        known = np.zeros(len(eligible))
        for pair in pred:
            if pair[0] not in gained:
                continue
            i = index[pair[0]]
            on_gained[i] += 1
            if pair in added_pairs:
                hit[i] += 1
            elif pair in delta.t0:
                known[i] += 1
        return hit, on_gained, known, len(pred)

    hit, on_gained, known, n_pred = tally(direct)
    hit_b, on_gained_b, known_b, _ = tally(top_by_q(direct, budget))
    added = np.zeros(len(eligible))
    for d, _ in added_pairs:
        added[index[d]] += 1
    return DomainCounts(
        domains=tuple(eligible),
        added=added,
        hit=hit,
        pred_on_gained=on_gained,
        known_on_gained=known,
        hit_b=hit_b,
        pred_on_gained_b=on_gained_b,
        known_on_gained_b=known_b,
        n_pred_pairs=n_pred,
        n_direct_eligible=len(direct),
        n_unmapped_signatures=rung.n_unmapped_signatures,
    )


def eligible_domains(delta: CuratedDelta, universe: Iterable[str] | None) -> list[str]:
    """Domains in both InterPro2GO snapshots and (if given) the run's universe."""
    domains = set(delta.shared_domains)
    if universe is not None:
        domains &= set(universe)
    return sorted(domains)


def default_budget(rungs: Iterable[RungTable], eligible: Iterable[str]) -> int:
    """The smallest rung's eligible direct-association count."""
    keep = set(eligible)
    return min(sum(1 for d, _, _ in r.direct if d in keep) for r in rungs)


def _percentile_ci(samples: np.ndarray, level: float) -> tuple[float, float]:
    finite = samples[np.isfinite(samples)]
    if finite.size == 0:
        return (float("nan"), float("nan"))
    alpha = (1.0 - level) / 2.0
    lo, hi = np.quantile(finite, [alpha, 1.0 - alpha])
    return (float(lo), float(hi))


def evaluate_rungs(
    rungs: Sequence[RungTable],
    delta: CuratedDelta,
    get_ancestors: GetAncestors,
    universe: Iterable[str] | None = None,
    budget: int | None = None,
    reference: str | None = None,
    n_bootstrap: int = 1000,
    seed: int = 0,
    level: float = 0.95,
) -> list[dict]:
    """One row per rung × metric, with paired domain-bootstrap intervals.

    ``budget`` defaults to the smallest rung's eligible association count.
    ``reference`` names the rung every other rung's delta is taken against
    (default: the first rung). Count rows carry no interval.
    """
    eligible = eligible_domains(delta, universe)
    if budget is None:
        budget = default_budget(rungs, eligible)
    reference = reference or rungs[0].name
    tables = {
        r.name: domain_counts(r, delta, eligible, get_ancestors, budget) for r in rungs
    }
    observed = {name: t.metrics() for name, t in tables.items()}

    n = len(eligible)
    rng = np.random.default_rng(seed)
    reps = {(name, m): np.empty(n_bootstrap) for name in tables for m in METRICS}
    for r in range(n_bootstrap):
        weights = np.bincount(rng.integers(0, n, size=n), minlength=n).astype(float)
        for name, table in tables.items():
            for m, value in table.metrics(weights).items():
                reps[(name, m)][r] = value

    rows: list[dict] = []
    for name, table in tables.items():
        for m in METRICS:
            lo, hi = _percentile_ci(reps[(name, m)], level)
            diff = reps[(name, m)] - reps[(reference, m)]
            d_lo, d_hi = _percentile_ci(diff, level)
            rows.append(
                {
                    "rung": name,
                    "metric": m,
                    "value": observed[name][m],
                    "ci_lo": lo,
                    "ci_hi": hi,
                    "reference": reference,
                    "delta": observed[name][m] - observed[reference][m],
                    "delta_ci_lo": d_lo,
                    "delta_ci_hi": d_hi,
                    "budget": budget,
                    "n_replicates": n_bootstrap,
                }
            )
        for m, value in table.counts().items():
            rows.append(
                {
                    "rung": name,
                    "metric": m,
                    "value": value,
                    "ci_lo": float("nan"),
                    "ci_hi": float("nan"),
                    "reference": reference,
                    "delta": float("nan"),
                    "delta_ci_lo": float("nan"),
                    "delta_ci_hi": float("nan"),
                    "budget": budget,
                    "n_replicates": 0,
                }
            )
    return rows


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #
def main() -> int:  # pragma: no cover - I/O wiring
    import argparse
    import sys

    from loguru import logger

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from src.ontology_processor import OntologyProcessor

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--rung",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="Labelled association table (repeatable); the significant "
        "associations or the True-Path-propagated file",
    )
    parser.add_argument(
        "--interpro2go-t0",
        type=Path,
        default=Path("data/raw/interpro_archive/85.0/interpro2go"),
    )
    parser.add_argument(
        "--interpro2go-t1",
        type=Path,
        default=Path("data/raw/interpro2go/interpro2go"),
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--interpro",
        type=Path,
        default=None,
        help="protein2ipr the rungs were trained on: fixes the domain universe "
        "and maps signature ids to entries under --domain-key ssf/pfam",
    )
    parser.add_argument(
        "--domain-key",
        choices=["interpro", *sorted(SIGNATURE_PREFIXES)],
        default="interpro",
    )
    parser.add_argument(
        "--budget",
        type=int,
        default=None,
        help="Top-N direct associations by q for precision_at_budget "
        "(default: the smallest rung's eligible count)",
    )
    parser.add_argument(
        "--reference",
        default=None,
        help="Rung the deltas are taken against (default: the first --rung)",
    )
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ci-level", type=float, default=0.95)
    parser.add_argument(
        "--output", type=Path, default=Path("validation/temporal_domain_centric.tsv")
    )
    args = parser.parse_args()

    for path in (args.interpro2go_t0, args.interpro2go_t1, args.go_ontology):
        if not path.exists():
            logger.error(f"Missing required input: {path}")
            return 1
    if args.domain_key != "interpro" and args.interpro is None:
        logger.error("--domain-key ssf/pfam needs --interpro to map signatures")
        return 1

    logger.info("Loading GO ontology...")
    get_ancestors = OntologyProcessor(args.go_ontology).get_ancestors
    delta = curated_delta(
        read_interpro2go(args.interpro2go_t0),
        read_interpro2go(args.interpro2go_t1),
        get_ancestors,
    )
    logger.info(
        f"InterPro2GO: {len(delta.shared_domains):,} shared domains, "
        f"{len(delta.added):,} added propagated pairs"
    )

    universe = sig_map = None
    if args.interpro is not None:
        universe, sig_map = read_protein2ipr_universe(args.interpro, args.domain_key)
        if args.domain_key == "interpro":
            sig_map = None
        logger.info(f"Universe: {len(universe):,} entries in {args.interpro}")

    rungs = []
    for spec in args.rung:
        name, _, path = spec.partition("=")
        if not path:
            logger.error(f"--rung must be NAME=PATH, got {spec!r}")
            return 1
        rung = load_rung(name, Path(path), sig_map)
        logger.info(
            f"  {name}: {len(rung.direct):,} single-domain rows"
            + (
                f", {rung.n_unmapped_signatures} unmapped signatures"
                if sig_map is not None
                else ""
            )
        )
        rungs.append(rung)

    rows = evaluate_rungs(
        rungs,
        delta,
        get_ancestors,
        universe=universe,
        budget=args.budget,
        reference=args.reference,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        level=args.ci_level,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    for row in rows:
        if row["metric"] in METRICS:
            logger.info(
                f"  {row['rung']:24s} {row['metric']:20s} {row['value']:.4f} "
                f"[{row['ci_lo']:.4f}, {row['ci_hi']:.4f}]"
            )
    logger.info(f"✓ Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
