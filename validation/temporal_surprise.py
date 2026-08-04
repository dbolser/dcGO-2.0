#!/usr/bin/env python3
"""Does a high surprise score at t0 predict future curation? — VALIDATION_PLAN §5.

`SURPRISE_SCORE.md` ranks emergent domain-combination associations, but that
ranking is computed on the same proteins that produced the associations: it
measures internal consistency of the evidence, not predictive power. This is the
held-out test.

The design, in one sentence: **take the associations dcGO found in 2021, ask
which proteins each one predicts a term for that they did not have in 2021, and
count how often curators added that term by 2026.**

Concretely, for an association (supra-domain ``S``, term ``t``) scored at ``t0``:

* **Predictions** — proteins carrying ``S`` that had *no* ``t`` at ``t0``
  (propagated, over the same non-IEA evidence the pipeline trains on, so the
  gate cannot leak a label the model already saw).
* **Hits** — of those, the ones annotated ``t`` (propagated, experimental
  evidence only) at ``t1``.
* **Base rate** — the same quantity for the term at large: of *all* domain-
  carrying proteins lacking ``t`` at ``t0``, what fraction gained it by ``t1``.

The headline statistic is **enrichment = hit rate / base rate**, pooled over the
associations in a stratum. The base rate is the crucial control: popular terms
accumulate annotations regardless of any prediction, and dividing by the term's
own acquisition rate removes exactly that.

Two comparisons matter, and they answer different questions:

1. **Is the ranking predictive at all?** Enrichment by surprise stratum
   (top-K, deciles) against the unranked pool and against 1.0.
2. **Does surprise add anything over plain significance?** The same candidate
   set ranked by the dcGO q-value instead. If surprise is only a proxy for
   "small p", the two curves coincide; if it is picking out emergence, the
   surprise curve should sit above it.

Deliberate limitations, so they are not oversold:

* Domain architectures come from the *current* ``protein2ipr`` — as in §2, the
  split is purely on the annotation side, so this is an annotation-temporal
  benchmark, not a fully prospective simulation.
* An association is scored over the proteins it predicts, so associations differ
  in how many predictions they make; both the pooled rate (weighted by
  predictions) and the per-association mean (unweighted) are reported.
* "Not annotated at t0" is not the same as "known to be absent" — GO is
  open-world. A miss may be a correct prediction that no one has curated yet, so
  absolute hit rates are lower bounds; the enrichment ratio is the defensible
  number because the base rate is depressed by exactly the same effect.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Sequence, Set, Tuple

GetAncestors = Callable[[str], Iterable[str]]


@dataclass(frozen=True)
class AssociationOutcome:
    """What became of one t0 association's predictions by t1.

    Attributes:
        feature: the supra-domain (or single domain) id.
        term: the predicted term.
        n_predicted: proteins carrying the feature that lacked the term at t0.
        n_hit: of those, the ones annotated with the term at t1.
        base_rate: the term's acquisition rate across all domain-carrying
            proteins that lacked it at t0.
        rank_scores: the scores this association is ranked by, keyed by ranking
            name (e.g. ``{"surprise": 6.7, "dcgo": 56.7}``).
    """

    feature: str
    term: str
    n_predicted: int
    n_hit: int
    base_rate: float
    rank_scores: Mapping[str, float]

    @property
    def hit_rate(self) -> float:
        return self.n_hit / self.n_predicted if self.n_predicted else 0.0

    @property
    def enrichment(self) -> float:
        """Hit rate over the term's own acquisition rate (inf if base rate is 0)."""
        if self.base_rate <= 0:
            return float("inf") if self.n_hit else 0.0
        return self.hit_rate / self.base_rate


@dataclass(frozen=True)
class StratumResult:
    """Pooled outcome for a set of associations."""

    name: str
    n_associations: int
    n_predicted: int
    n_hit: int
    hit_rate: float
    expected_rate: float
    enrichment: float
    mean_per_association_enrichment: float
    ci_low: float
    ci_high: float


def propagate(terms: Iterable[str], get_ancestors: GetAncestors) -> Set[str]:
    """A term set plus all its ancestors (the True Path Rule, applied to truth)."""
    out: Set[str] = set()
    for term in terms:
        out.add(term)
        out.update(get_ancestors(term))
    return out


def acquisition_base_rates(
    terms: Iterable[str],
    t0_map: Mapping[str, Set[str]],
    t1_map: Mapping[str, Set[str]],
    universe: Set[str],
) -> Dict[str, float]:
    """For each term, the fraction of eligible proteins that gained it by t1.

    Eligible = in ``universe`` (has domains) and lacking the term at ``t0``. Both
    maps must already be propagated.

    The two snapshots are inverted to term → proteins once, rather than scanning
    the universe per term: with ~10k terms over ~19k proteins the direct loop is
    200M membership tests, the inversion is one pass over the annotations.
    """
    wanted = set(terms)
    t0_index: Dict[str, Set[str]] = {}
    t1_index: Dict[str, Set[str]] = {}
    for source, index in ((t0_map, t0_index), (t1_map, t1_index)):
        for protein, protein_terms in source.items():
            if protein not in universe:
                continue
            for term in protein_terms:
                if term in wanted:
                    index.setdefault(term, set()).add(protein)

    n_universe = len(universe)
    rates: Dict[str, float] = {}
    for term in wanted:
        had_at_t0 = t0_index.get(term, frozenset())
        eligible = n_universe - len(had_at_t0)
        gained = len(t1_index.get(term, frozenset()) - had_at_t0)
        rates[term] = gained / eligible if eligible else 0.0
    return rates


def score_association(
    feature: str,
    term: str,
    carriers: Set[str],
    t0_map: Mapping[str, Set[str]],
    t1_map: Mapping[str, Set[str]],
    base_rate: float,
    rank_scores: Mapping[str, float],
) -> AssociationOutcome:
    """Count predictions and hits for one association (both maps propagated)."""
    predicted = [p for p in carriers if term not in t0_map.get(p, ())]
    hits = sum(1 for p in predicted if term in t1_map.get(p, ()))
    return AssociationOutcome(
        feature=feature,
        term=term,
        n_predicted=len(predicted),
        n_hit=hits,
        base_rate=base_rate,
        rank_scores=rank_scores,
    )


def _bootstrap_enrichment_ci(
    outcomes: Sequence[AssociationOutcome],
    n_resamples: int,
    seed: int,
    confidence: float = 0.95,
) -> Tuple[float, float]:
    """Percentile bootstrap CI for pooled enrichment, resampling *associations*.

    Resampling associations rather than proteins is the conservative choice: the
    proteins behind one association are not independent of each other (they
    share an architecture), but different associations largely are.
    """
    import random

    if not outcomes or n_resamples <= 0:
        return (float("nan"), float("nan"))
    rng = random.Random(seed)
    indices = range(len(outcomes))
    samples: List[float] = []
    for _ in range(n_resamples):
        drawn = [outcomes[rng.choice(indices)] for _ in indices]
        hits = sum(o.n_hit for o in drawn)
        predicted = sum(o.n_predicted for o in drawn)
        expected = sum(o.base_rate * o.n_predicted for o in drawn)
        if predicted and expected > 0:
            samples.append((hits / predicted) / (expected / predicted))
    if not samples:
        return (float("nan"), float("nan"))
    samples.sort()
    tail = (1.0 - confidence) / 2
    low = samples[min(len(samples) - 1, int(tail * len(samples)))]
    high = samples[min(len(samples) - 1, int((1 - tail) * len(samples)))]
    return (low, high)


def pool(
    name: str,
    outcomes: Sequence[AssociationOutcome],
    n_resamples: int = 1000,
    seed: int = 0,
) -> StratumResult:
    """Pool a stratum: hit rate, base-rate-adjusted enrichment, bootstrap CI.

    ``expected_rate`` is the prediction-weighted mean base rate — what this exact
    set of (protein, term) predictions would have hit by chance.
    """
    n_predicted = sum(o.n_predicted for o in outcomes)
    n_hit = sum(o.n_hit for o in outcomes)
    expected_hits = sum(o.base_rate * o.n_predicted for o in outcomes)
    hit_rate = n_hit / n_predicted if n_predicted else 0.0
    expected_rate = expected_hits / n_predicted if n_predicted else 0.0
    enrichment = hit_rate / expected_rate if expected_rate > 0 else float("nan")
    finite = [o.enrichment for o in outcomes if math.isfinite(o.enrichment)]
    low, high = _bootstrap_enrichment_ci(outcomes, n_resamples, seed)
    return StratumResult(
        name=name,
        n_associations=len(outcomes),
        n_predicted=n_predicted,
        n_hit=n_hit,
        hit_rate=hit_rate,
        expected_rate=expected_rate,
        enrichment=enrichment,
        mean_per_association_enrichment=sum(finite) / len(finite) if finite else 0.0,
        ci_low=low,
        ci_high=high,
    )


def top_k(
    outcomes: Sequence[AssociationOutcome], ranking: str, k: int
) -> List[AssociationOutcome]:
    """The ``k`` highest-scoring associations under one ranking."""
    ranked = sorted(
        outcomes, key=lambda o: o.rank_scores.get(ranking, float("-inf")), reverse=True
    )
    return ranked[:k]


@dataclass(frozen=True)
class RankingComparison:
    """Paired head-to-head between two rankings at one prediction budget."""

    budget: int
    ranking_a: str
    ranking_b: str
    enrichment_a: float
    enrichment_b: float
    difference: float
    ci_low: float
    ci_high: float
    fraction_favouring_a: float

    @property
    def separated(self) -> bool:
        """True when the bootstrap CI for the difference excludes zero."""
        return self.ci_low > 0 or self.ci_high < 0


def _budget_enrichment(
    outcomes: Sequence[AssociationOutcome], ranking: str, budget: int
) -> float:
    """Pooled enrichment of the top associations under ``ranking`` up to ``budget``."""
    ranked = sorted(
        outcomes, key=lambda o: o.rank_scores.get(ranking, float("-inf")), reverse=True
    )
    hits = predicted = 0
    expected = 0.0
    for o in ranked:
        if predicted >= budget:
            break
        hits += o.n_hit
        predicted += o.n_predicted
        expected += o.base_rate * o.n_predicted
    if not predicted or expected <= 0:
        return float("nan")
    return hits / expected


def compare_rankings(
    outcomes: Sequence[AssociationOutcome],
    ranking_a: str,
    ranking_b: str,
    budget: int,
    n_resamples: int = 1000,
    seed: int = 0,
) -> RankingComparison:
    """Paired bootstrap of the enrichment *difference* between two rankings.

    Comparing two independently-computed confidence intervals is the wrong test
    here: both rankings are applied to the *same* candidate pool, so their
    estimates are correlated and overlapping intervals do not imply no
    difference. This resamples the pool once per iteration and re-ranks *both*
    ways inside each resample, so the difference is measured on identical data
    every time.

    ``fraction_favouring_a`` is the share of resamples where ``ranking_a`` came
    out ahead — a one-sided bootstrap p-value of ``1 - fraction`` for the
    hypothesis that A beats B.
    """
    import random

    point_a = _budget_enrichment(outcomes, ranking_a, budget)
    point_b = _budget_enrichment(outcomes, ranking_b, budget)

    rng = random.Random(seed)
    indices = range(len(outcomes))
    diffs: List[float] = []
    for _ in range(n_resamples):
        drawn = [outcomes[rng.choice(indices)] for _ in indices]
        a = _budget_enrichment(drawn, ranking_a, budget)
        b = _budget_enrichment(drawn, ranking_b, budget)
        if math.isfinite(a) and math.isfinite(b):
            diffs.append(a - b)
    if not diffs:
        return RankingComparison(
            budget,
            ranking_a,
            ranking_b,
            point_a,
            point_b,
            float("nan"),
            float("nan"),
            float("nan"),
            float("nan"),
        )
    diffs.sort()
    low = diffs[int(0.025 * len(diffs))]
    high = diffs[min(len(diffs) - 1, int(0.975 * len(diffs)))]
    return RankingComparison(
        budget=budget,
        ranking_a=ranking_a,
        ranking_b=ranking_b,
        enrichment_a=point_a,
        enrichment_b=point_b,
        difference=point_a - point_b,
        ci_low=low,
        ci_high=high,
        fraction_favouring_a=sum(1 for d in diffs if d > 0) / len(diffs),
    )


def strata_by_prediction_budget(
    outcomes: Sequence[AssociationOutcome],
    ranking: str,
    budgets: Sequence[int],
) -> List[Tuple[str, List[AssociationOutcome]]]:
    """Take associations in rank order until a prediction *budget* is filled.

    Top-K is not a fair comparison between rankings here: the surprise score
    favours tight architectures carried by a handful of proteins, while raw
    significance favours common domains, so the same K exposes an order of
    magnitude more predictions for one ranking than the other. Equalising the
    number of (protein, term) predictions instead asks the question a curator
    would: *given capacity to check N predictions, which ranking finds more
    that come true?*

    The association that crosses the budget is included whole, so a stratum may
    slightly exceed it.
    """
    ranked = sorted(
        outcomes, key=lambda o: o.rank_scores.get(ranking, float("-inf")), reverse=True
    )
    out: List[Tuple[str, List[AssociationOutcome]]] = []
    for budget in budgets:
        taken: List[AssociationOutcome] = []
        total = 0
        for o in ranked:
            if total >= budget:
                break
            taken.append(o)
            total += o.n_predicted
        if taken and total < sum(o.n_predicted for o in ranked):
            out.append((f"{ranking} @{budget} preds", taken))
    return out


def strata_by_rank(
    outcomes: Sequence[AssociationOutcome],
    ranking: str,
    cutoffs: Sequence[int],
) -> List[Tuple[str, List[AssociationOutcome]]]:
    """Nested top-K slices plus the remainder, for one ranking."""
    ranked = sorted(
        outcomes, key=lambda o: o.rank_scores.get(ranking, float("-inf")), reverse=True
    )
    out: List[Tuple[str, List[AssociationOutcome]]] = []
    for k in cutoffs:
        if k < len(ranked):
            out.append((f"{ranking} top-{k}", ranked[:k]))
    out.append((f"{ranking} all ({len(ranked)})", list(ranked)))
    return out


def _report(args, scored, logger) -> int:  # pragma: no cover - I/O wiring
    """Strata, paired ranking comparison and output files, shared by both paths."""
    results = []
    for ranking in ("surprise", "dcgo"):
        for name, stratum in strata_by_rank(scored, ranking, args.cutoffs):
            results.append(pool(name, stratum, args.bootstrap, args.seed))
    # The fair head-to-head: equal prediction volume, not equal association count.
    for ranking in ("surprise", "dcgo"):
        for name, stratum in strata_by_prediction_budget(scored, ranking, args.budgets):
            results.append(pool(name, stratum, args.bootstrap, args.seed))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        f.write(
            "stratum\tn_associations\tn_predictions\tn_hits\thit_rate\t"
            "expected_rate\tenrichment\tenrichment_ci_low\tenrichment_ci_high\t"
            "mean_per_association_enrichment\n"
        )
        for r in results:
            f.write(
                f"{r.name}\t{r.n_associations}\t{r.n_predicted}\t{r.n_hit}\t"
                f"{r.hit_rate:.4f}\t{r.expected_rate:.4f}\t{r.enrichment:.2f}\t"
                f"{r.ci_low:.2f}\t{r.ci_high:.2f}\t"
                f"{r.mean_per_association_enrichment:.2f}\n"
            )
    logger.info(f"✓ Wrote {len(results)} strata to {args.output}")

    # The documented replay command points --replay at the default
    # --per-association-output, so writing unconditionally would re-serialise
    # the file just read, re-rounding base_rate to 5 dp and surprise to 3 dp on
    # every replay until the stored values drifted away from the run's.
    replaying_own_input = (
        args.replay is not None
        and args.replay.resolve() == args.per_association_output.resolve()
    )
    if replaying_own_input:
        logger.info(
            f"Replay source is the per-association output "
            f"({args.per_association_output}); leaving it untouched."
        )
    else:
        with open(args.per_association_output, "w") as f:
            f.write(
                "domain\tterm\tsurprise\tdcgo_score\tn_predicted\tn_hit\tbase_rate\tenrichment\n"
            )
            for o in sorted(scored, key=lambda o: -o.rank_scores["surprise"]):
                f.write(
                    f"{o.feature}\t{o.term}\t{o.rank_scores['surprise']:.3f}\t"
                    f"{o.rank_scores['dcgo']:.2f}\t{o.n_predicted}\t{o.n_hit}\t"
                    f"{o.base_rate:.5f}\t{o.enrichment:.2f}\n"
                )
        logger.info(
            f"✓ Wrote per-association outcomes to {args.per_association_output}"
        )

    comparisons = [
        compare_rankings(scored, "surprise", "dcgo", budget, args.bootstrap, args.seed)
        for budget in args.budgets
    ]

    logger.info("")
    logger.info(
        f"{'stratum':<26} {'assoc':>6} {'preds':>7} {'hits':>6} "
        f"{'hit%':>6} {'exp%':>6} {'enrich':>7} {'95% CI':>14}"
    )
    for r in results:
        logger.info(
            f"{r.name:<26} {r.n_associations:>6,} {r.n_predicted:>7,} {r.n_hit:>6,} "
            f"{100 * r.hit_rate:>5.1f}% {100 * r.expected_rate:>5.1f}% "
            f"{r.enrichment:>7.2f} {f'[{r.ci_low:.2f}, {r.ci_high:.2f}]':>14}"
        )

    logger.info("")
    logger.info("Paired head-to-head (same resampled pool, both rankings):")
    for c in comparisons:
        verdict = "separated" if c.separated else "not separated"
        logger.info(
            f"  @{c.budget:>6,} preds: surprise {c.enrichment_a:6.2f} vs "
            f"dcgo {c.enrichment_b:6.2f}  diff {c.difference:+6.2f} "
            f"[{c.ci_low:+.2f}, {c.ci_high:+.2f}]  "
            f"favours surprise in {100 * c.fraction_favouring_a:.0f}% of resamples "
            f"({verdict})"
        )
    with open(args.output, "a") as f:
        f.write("\n# paired ranking comparison (surprise - dcgo)\n")
        f.write(
            "budget\tenrichment_surprise\tenrichment_dcgo\tdifference\t"
            "diff_ci_low\tdiff_ci_high\tfraction_favouring_surprise\n"
        )
        for c in comparisons:
            f.write(
                f"{c.budget}\t{c.enrichment_a:.2f}\t{c.enrichment_b:.2f}\t"
                f"{c.difference:.2f}\t{c.ci_low:.2f}\t{c.ci_high:.2f}\t"
                f"{c.fraction_favouring_a:.3f}\n"
            )
    return 0


def main() -> int:  # pragma: no cover - I/O wiring, exercised by running it
    import argparse
    import csv
    import sys
    from collections import defaultdict
    from pathlib import Path

    from loguru import logger

    project_root = Path(__file__).resolve().parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    from src.domain_annotation_parser import DomainAnnotationParser
    from src.goa_parser import EXPERIMENTAL_EVIDENCE, GOAParser, parse_goa_human
    from src.ontology_processor import OntologyProcessor

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--surprising",
        type=Path,
        default=Path("results_t0_2021/domain_go_surprising.tsv"),
        help="Surprise ranking computed on the t0 associations",
    )
    parser.add_argument(
        "--t0-gaf",
        type=Path,
        default=Path("data/raw/goa_archive/goa_human.gaf.205.gz"),
    )
    parser.add_argument(
        "--t1-gaf",
        type=Path,
        default=Path("data/raw/goa_annotations/goa_human.gaf.gz"),
    )
    parser.add_argument(
        "--interpro",
        type=Path,
        default=Path("data/interim/protein2ipr_human.dat.gz"),
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation/temporal_surprise_metrics.tsv"),
    )
    parser.add_argument(
        "--per-association-output",
        type=Path,
        default=Path("validation/temporal_surprise_associations.tsv"),
    )
    parser.add_argument(
        "--cutoffs",
        type=int,
        nargs="+",
        default=[25, 100, 500, 2000],
        help="Top-K slices to report for each ranking",
    )
    parser.add_argument(
        "--budgets",
        type=int,
        nargs="+",
        default=[2000, 10000, 40000],
        help="Prediction budgets for the ranking-fair comparison: how many "
        "(protein, term) predictions a curator could afford to check",
    )
    parser.add_argument(
        "--replay",
        type=Path,
        default=None,
        help="Skip all parsing and re-analyse a previous run's "
        "per-association TSV (the outcomes are already counted in it)",
    )
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.replay:
        if not args.replay.exists():
            logger.error(f"Missing replay file: {args.replay}")
            return 1
    else:
        for path in (
            args.surprising,
            args.t0_gaf,
            args.t1_gaf,
            args.interpro,
            args.go_ontology,
        ):
            if not path.exists():
                logger.error(f"Missing required input: {path}")
                return 1

    if args.replay:
        logger.info(f"Replaying counted outcomes from {args.replay}...")
        scored = []
        with open(args.replay) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                scored.append(
                    AssociationOutcome(
                        feature=row["domain"],
                        term=row["term"],
                        n_predicted=int(row["n_predicted"]),
                        n_hit=int(row["n_hit"]),
                        base_rate=float(row["base_rate"]),
                        rank_scores={
                            "surprise": float(row["surprise"]),
                            "dcgo": float(row["dcgo_score"]),
                        },
                    )
                )
        logger.info(f"  Associations: {len(scored):,}")
        return _report(args, scored, logger)

    logger.info("Loading GO ontology (for propagation)...")
    processor = OntologyProcessor(args.go_ontology)
    get_ancestors = processor.get_ancestors

    logger.info("Parsing t0 GOA (non-IEA — the evidence the pipeline trains on)...")
    t0_raw = parse_goa_human(args.t0_gaf, evidence_filter="manual")
    logger.info("Parsing t1 GOA (experimental evidence only — the truth)...")
    t1_raw = GOAParser(
        evidence_codes=EXPERIMENTAL_EVIDENCE, aspects={"P", "F", "C"}
    ).parse_gaf_file(args.t1_gaf)

    logger.info("Propagating both snapshots up the GO DAG...")
    t0_map = {p: propagate(terms, get_ancestors) for p, terms in t0_raw.items()}
    t1_map = {p: propagate(terms, get_ancestors) for p, terms in t1_raw.items()}

    logger.info("Parsing domain architectures...")
    architectures = DomainAnnotationParser(
        max_supra_domain_length=3, min_domain_length=10
    ).parse_protein2ipr_file(args.interpro)
    universe = set(architectures)
    logger.info(f"  Proteins with domains: {len(universe):,}")

    logger.info(f"Reading the t0 surprise ranking from {args.surprising}...")
    rows = []
    with open(args.surprising) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    logger.info(f"  Candidate associations: {len(rows):,}")
    if not rows:
        logger.error("No ranked associations to evaluate.")
        return 1

    # Index only the features we need, over the same universe.
    wanted = {row["domain"] for row in rows}
    carriers: Dict[str, Set[str]] = defaultdict(set)
    for protein, arch in architectures.items():
        for feature in arch.single_domains:
            if feature in wanted:
                carriers[feature].add(protein)
        for feature in arch.supra_domains:
            if feature in wanted:
                carriers[feature].add(protein)

    logger.info("Computing per-term acquisition base rates (t0 → t1)...")
    base_rates = acquisition_base_rates(
        {row["term"] for row in rows}, t0_map, t1_map, universe
    )

    logger.info("Scoring associations against the t1 snapshot...")
    outcomes = []
    for row in rows:
        term = row["term"]
        outcomes.append(
            score_association(
                feature=row["domain"],
                term=term,
                carriers=carriers.get(row["domain"], set()),
                t0_map=t0_map,
                t1_map=t1_map,
                base_rate=base_rates.get(term, 0.0),
                rank_scores={
                    "surprise": float(row["surprise"]),
                    # The incumbent ranking: plain dcGO significance.
                    "dcgo": -math.log10(max(float(row["dcgo_adj_p_value"]), 1e-320)),
                },
            )
        )
    scored = [o for o in outcomes if o.n_predicted > 0]
    logger.info(
        f"  Associations making at least one prediction: {len(scored):,} "
        f"({sum(o.n_predicted for o in scored):,} protein-term predictions)"
    )

    return _report(args, scored, logger)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
