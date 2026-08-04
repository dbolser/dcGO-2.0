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

A note on comparison 2, because it is the fragile one. Matching on prediction
budget means the *set* of associations being scored is chosen inside every
bootstrap resample, and a budget-truncated pooled ratio is a hard thing to
resample honestly: the selected set jumps discontinuously, a duplicated
high-rank association can fill the budget with copies of itself, and — on this
data — the ranking assigns the same score to 98% of the pool, so most of a
"top of the ranking" is an arbitrary tie-break. Everything from
:func:`summarise_bootstrap` downwards exists to make that visible rather than to
hide it behind one interval: percentile, basic and BCa are all reported, so are
the diagnostics that say whether any of them is meaningful, and
:func:`tie_break_spread` measures how much of a result is input order. When the
answer is "this comparison is not resolvable with this design", that is what the
output says.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from statistics import NormalDist
from typing import (
    Callable,
    ClassVar,
    Dict,
    Iterable,
    List,
    Mapping,
    Sequence,
    Set,
    Tuple,
)

GetAncestors = Callable[[str], Iterable[str]]

_NORMAL = NormalDist()


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


# --------------------------------------------------------------------------- #
# Bootstrap interval machinery
#
# The budget-matched ranking comparison below turned out to have a bootstrap
# distribution that is *not* centred on the observed statistic, which makes the
# choice of interval formula change the answer. Rather than pick one silently we
# compute all three standard intervals and a set of diagnostics that say whether
# any of them can be trusted.
# --------------------------------------------------------------------------- #


def quantile(sorted_values: Sequence[float], p: float) -> float:
    """Linear-interpolated (type-7) quantile of an already-sorted sequence.

    The previous code indexed with ``int(p * len(samples))``, which for the
    2.5% tail of 1,000 samples is the 25th order statistic rather than a
    quantile — fine at the centre, biased outward in the tails that a CI is
    made of.
    """
    if not sorted_values or not math.isfinite(p):
        return float("nan")
    if len(sorted_values) == 1:
        return sorted_values[0]
    p = min(max(p, 0.0), 1.0)
    h = (len(sorted_values) - 1) * p
    lo = math.floor(h)
    hi = math.ceil(h)
    return sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo])


def _skew(values: Sequence[float], mean: float, sd: float) -> float:
    if sd <= 0 or len(values) < 3:
        return float("nan")
    return sum((v - mean) ** 3 for v in values) / (len(values) * sd**3)


def _acceleration(jackknife: Sequence[float]) -> float:
    """Efron's acceleration from leave-one-out replicates (0.0 if degenerate).

    Leaving out the only association in a thin budget slice makes the statistic
    undefined; those replicates carry no shape information and are dropped
    rather than poisoning the whole estimate with a NaN.
    """
    values = [v for v in jackknife if math.isfinite(v)]
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    d2 = 0.0
    d3 = 0.0
    for v in values:
        d = mean - v
        d2 += d * d
        d3 += d * d * d
    if d2 <= 0:
        return 0.0
    return d3 / (6.0 * d2**1.5)


@dataclass(frozen=True)
class BootstrapSummary:
    """Percentile, basic and BCa intervals for one statistic, plus diagnostics.

    ``trustworthy`` is the honest gate. A bootstrap interval of any flavour
    assumes the resampling distribution of ``theta* - theta_hat`` approximates
    the sampling distribution of ``theta_hat - theta``. When the observed
    statistic falls outside its own resampling distribution's central mass that
    assumption has failed, percentile and basic intervals disagree violently
    (they are reflections of each other about ``theta_hat``), and BCa's bias
    correction is being asked to extrapolate. In that case no single interval
    should be quoted.
    """

    point: float
    n_resamples: int
    n_usable: int
    mean: float
    median: float
    sd: float
    skew: float
    z0: float
    acceleration: float
    fraction_positive: float
    percentile: Tuple[float, float]
    basic: Tuple[float, float]
    bca: Tuple[float, float]
    trustworthy: bool
    notes: Tuple[str, ...] = ()

    @property
    def point_inside_percentile(self) -> bool:
        low, high = self.percentile
        return low <= self.point <= high

    @property
    def note(self) -> str:
        return "; ".join(self.notes) if self.notes else "ok"

    @property
    def recommended(self) -> Tuple[float, float]:
        """The interval to quote — BCa when trustworthy, otherwise nothing."""
        if not self.trustworthy:
            return (float("nan"), float("nan"))
        return self.bca


#: |z0| above which the bias correction is extrapolating rather than correcting.
#: z0 = 0.5 means only ~31% of resamples fall below the observed statistic; Efron
#: & Tibshirani's worked examples sit an order of magnitude below that.
MAX_TRUSTWORTHY_Z0 = 0.5


def summarise_bootstrap(
    point: float,
    samples: Sequence[float],
    jackknife: Sequence[float] = (),
    confidence: float = 0.95,
    min_resamples: int = 200,
) -> BootstrapSummary:
    """Summarise a bootstrap distribution and decide whether to trust it.

    Args:
        point: the statistic computed on the observed data.
        samples: bootstrap replicates of the same statistic.
        jackknife: leave-one-out replicates, needed for BCa's acceleration.
            Without them BCa degenerates to a bias-corrected (BC) interval.
        confidence: two-sided coverage.
        min_resamples: below this the tails are not resolved at all.
    """
    nan2 = (float("nan"), float("nan"))
    usable = sorted(s for s in samples if math.isfinite(s))
    n = len(usable)
    if n == 0 or not math.isfinite(point):
        return BootstrapSummary(
            point=point,
            n_resamples=len(samples),
            n_usable=n,
            mean=float("nan"),
            median=float("nan"),
            sd=float("nan"),
            skew=float("nan"),
            z0=float("nan"),
            acceleration=float("nan"),
            fraction_positive=float("nan"),
            percentile=nan2,
            basic=nan2,
            bca=nan2,
            trustworthy=False,
            notes=("no usable bootstrap replicates",),
        )

    tail = (1.0 - confidence) / 2
    mean = sum(usable) / n
    sd = math.sqrt(sum((v - mean) ** 2 for v in usable) / (n - 1)) if n > 1 else 0.0
    pct = (quantile(usable, tail), quantile(usable, 1 - tail))
    # Basic (reverse-percentile): reflect the replicate quantiles about theta_hat.
    basic = (2 * point - pct[1], 2 * point - pct[0])

    notes: List[str] = []
    below = sum(1 for v in usable if v < point)
    equal = sum(1 for v in usable if v == point)
    prop_below = (below + 0.5 * equal) / n
    if prop_below <= 0.0 or prop_below >= 1.0:
        z0 = math.copysign(float("inf"), prop_below - 0.5)
        notes.append(
            "every bootstrap replicate falls on one side of the observed "
            "statistic (z0 undefined)"
        )
    else:
        z0 = _NORMAL.inv_cdf(prop_below)
    accel = _acceleration(jackknife) if len(jackknife) >= 2 else 0.0
    if len(jackknife) < 2:
        notes.append("no jackknife replicates: BCa reduces to BC (acceleration 0)")

    bca = nan2
    if math.isfinite(z0):
        z_lo = _NORMAL.inv_cdf(tail)
        z_hi = _NORMAL.inv_cdf(1 - tail)
        alphas: List[float] = []
        for z in (z_lo, z_hi):
            denom = 1 - accel * (z0 + z)
            if denom == 0 or not math.isfinite(denom):
                alphas = []
                notes.append("BCa acceleration is degenerate (zero denominator)")
                break
            alpha = _NORMAL.cdf(z0 + (z0 + z) / denom)
            if not math.isfinite(alpha):
                alphas = []
                notes.append("BCa tail probability is not finite")
                break
            alphas.append(alpha)
        if len(alphas) == 2:
            a_lo, a_hi = alphas
            if a_lo * n < 1 or (1 - a_hi) * n < 1:
                notes.append(
                    f"BCa tail probabilities ({a_lo:.4f}, {a_hi:.4f}) are finer "
                    f"than {n:,} resamples can resolve"
                )
            bca = (quantile(usable, a_lo), quantile(usable, a_hi))

    if n < min_resamples:
        notes.append(f"only {n:,} usable resamples")
    if len(samples) and n < 0.95 * len(samples):
        notes.append(
            f"{len(samples) - n:,}/{len(samples):,} resamples were non-finite "
            "and discarded, which is itself selection"
        )
    inside = pct[0] <= point <= pct[1]
    if not inside:
        notes.append(
            "the observed statistic lies OUTSIDE its own percentile interval: "
            "the bootstrap distribution is not centred on it, so no first-order "
            "interval is valid here"
        )
    if math.isfinite(z0) and abs(z0) > MAX_TRUSTWORTHY_Z0:
        notes.append(
            f"bias correction z0={z0:+.2f} exceeds |{MAX_TRUSTWORTHY_Z0}|: the "
            "median replicate is far from the observed statistic"
        )

    trustworthy = (
        inside
        and n >= min_resamples
        and math.isfinite(z0)
        and abs(z0) <= MAX_TRUSTWORTHY_Z0
        and math.isfinite(bca[0])
        and math.isfinite(bca[1])
    )
    return BootstrapSummary(
        point=point,
        n_resamples=len(samples),
        n_usable=n,
        mean=mean,
        median=quantile(usable, 0.5),
        sd=sd,
        skew=_skew(usable, mean, sd),
        z0=z0,
        acceleration=accel,
        fraction_positive=sum(1 for v in usable if v > 0) / n,
        percentile=pct,
        basic=basic,
        bca=bca,
        trustworthy=trustworthy,
        notes=tuple(notes),
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

    Unlike the budget-matched comparison further down, this statistic selects no
    subset inside the resample — the stratum is fixed — so its bootstrap
    distribution is well behaved and the percentile interval is fine.
    """
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
    """Paired head-to-head between two rankings at one prediction budget.

    Two resampling *designs* are reported for the same difference, because they
    answer different questions and — on this data — behave very differently:

    ``reselect``
        The original design. Resample the pool, then re-fill the budget from
        scratch inside every resample. This propagates the uncertainty in *which*
        associations a ranking would put at the top, but the budget cutoff makes
        the selected set jump discontinuously, and a duplicated high-rank
        association can fill the whole budget with copies of itself.
    ``fixed``
        Condition on the selection. The two budget slices are the ones the
        observed ranking produced; a resample only re-weights the associations
        (each gets a Binomial(n, 1/n) ≈ Poisson(1) multiplicity from the same
        single draw, so members shared by both slices move together and the
        comparison stays paired). This answers "are the sets these rankings
        actually selected different?" and drops the reselection discontinuity.

    Neither is strictly better; ``fixed`` understates selection uncertainty,
    ``reselect`` is barely a bootstrap at all when few associations fit the
    budget. Reporting both is the honest option.
    """

    budget: int
    ranking_a: str
    ranking_b: str
    enrichment_a: float
    enrichment_b: float
    difference: float
    n_associations_a: int
    n_associations_b: int
    reselect: BootstrapSummary
    fixed: BootstrapSummary
    #: Largest share of one slice's predictions contributed by a single
    #: association, per ranking — the leverage the bootstrap has to shift.
    max_leverage_a: float = float("nan")
    max_leverage_b: float = float("nan")
    #: How much of each slice is decided by an arbitrary tie-break.
    tie_a: "TieAmbiguity | None" = None
    tie_b: "TieAmbiguity | None" = None
    #: The same difference under random re-orderings of the tied score blocks.
    tie_break: "TieBreakSpread | None" = None
    notes: Tuple[str, ...] = field(default_factory=tuple)

    #: Below this many associations inside the budget, resampling the slice is
    #: not a meaningful bootstrap however the interval is computed.
    MIN_ASSOCIATIONS: ClassVar[int] = 30

    @property
    def thin_basis(self) -> bool:
        """True when one ranking fits too few associations inside the budget."""
        return min(self.n_associations_a, self.n_associations_b) < self.MIN_ASSOCIATIONS

    @property
    def resolvable(self) -> bool:
        """Whether any interval here should be quoted as a result."""
        return (
            not self.thin_basis
            and not self.mostly_tied
            and self.fixed.trustworthy
            and self.reselect.trustworthy
        )

    @property
    def separated(self) -> bool:
        """True only when a *trustworthy* interval excludes zero."""
        if not self.resolvable:
            return False
        low, high = self.fixed.recommended
        return low > 0 or high < 0

    @property
    def mostly_tied(self) -> bool:
        """True when over half of a slice is inside one tied score block."""
        return any(t is not None and t.share > 0.5 for t in (self.tie_a, self.tie_b))

    @property
    def verdict(self) -> str:
        if self.mostly_tied:
            return "unresolvable: most of a slice is an arbitrary tie-break"
        if self.thin_basis:
            return (
                f"unresolvable: only {min(self.n_associations_a, self.n_associations_b)}"
                f" associations fit the budget"
            )
        if not self.reselect.trustworthy and not self.fixed.trustworthy:
            return "unresolvable: no trustworthy interval"
        if not self.reselect.trustworthy:
            return "reselect design unreliable; fixed-selection interval only"
        if not self.fixed.trustworthy:
            return "fixed design unreliable; reselect interval only"
        return "separated" if self.separated else "not separated"


def _budget_enrichment(
    outcomes: Sequence[AssociationOutcome], ranking: str, budget: int
) -> float:
    """Pooled enrichment of the top associations under ``ranking`` up to ``budget``.

    Reference implementation: readable, O(n log n), and the thing the fast
    resampling path below is tested against.
    """
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


class _RankedPool:
    """One ranking of one pool, with prefix sums, for repeated budget queries.

    The rank *order* never changes under resampling (scores are fixed per
    association), so it is computed once; a resample only changes multiplicities.

    Ties are broken by original index unless a ``jitter`` is supplied, and that
    is deliberate. The previous implementation re-sorted the *drawn* list inside
    every resample, and Python's sort is stable in draw order, so associations
    sharing a score were ordered at random in each replicate while the observed
    statistic used one fixed order. With 98% of this pool tied at surprise 0.000
    that silently turned the bootstrap into a tie-break lottery, which is a large
    part of why the committed metrics and a replay of them disagreed. Tie-break
    variability is now measured on purpose, by :func:`tie_break_spread`.
    """

    def __init__(
        self,
        outcomes: Sequence[AssociationOutcome],
        ranking: str,
        jitter: Sequence[float] | None = None,
    ) -> None:
        self.n = len(outcomes)
        # Descending index, so that with reverse=True ties keep input order.
        secondary = jitter if jitter is not None else [-i for i in range(self.n)]
        self.order = sorted(
            range(self.n),
            key=lambda i: (
                outcomes[i].rank_scores.get(ranking, float("-inf")),
                secondary[i],
            ),
            reverse=True,
        )
        self.rank = [0] * self.n
        for position, index in enumerate(self.order):
            self.rank[index] = position
        self.pred = [outcomes[i].n_predicted for i in self.order]
        self.hit = [outcomes[i].n_hit for i in self.order]
        self.exp = [outcomes[i].base_rate * outcomes[i].n_predicted for i in self.order]
        self.cum_pred: List[int] = []
        self.cum_hit: List[int] = []
        self.cum_exp: List[float] = []
        p = h = 0
        e = 0.0
        for position in range(self.n):
            p += self.pred[position]
            h += self.hit[position]
            e += self.exp[position]
            self.cum_pred.append(p)
            self.cum_hit.append(h)
            self.cum_exp.append(e)

    def cutoff(self, budget: int) -> int:
        """Number of associations taken: every one whose *preceding* total < budget."""
        import bisect

        if self.n == 0:
            return 0
        return 1 + min(bisect.bisect_left(self.cum_pred, budget), self.n - 1)

    def slice_indices(self, budget: int) -> List[int]:
        """Original indices of the associations inside ``budget``, in rank order."""
        return self.order[: self.cutoff(budget)]

    def enrichment(self, budget: int) -> float:
        k = self.cutoff(budget)
        if k == 0 or self.cum_exp[k - 1] <= 0:
            return float("nan")
        return self.cum_hit[k - 1] / self.cum_exp[k - 1]

    def enrichment_from_counts(self, counts: Sequence[int], budget: int) -> float:
        """Budget enrichment of a resample given each association's multiplicity.

        Copies of one association are adjacent in the ranking, so the number of
        copies that fit is arithmetic rather than a loop — this is what makes a
        10,000-resample run affordable.
        """
        hits = predicted = 0
        expected = 0.0
        for position, index in enumerate(self.order):
            if predicted >= budget:
                break
            multiplicity = counts[index]
            if not multiplicity:
                continue
            per = self.pred[position]
            room = budget - predicted
            take = multiplicity
            if per > 0:
                take = min(multiplicity, -(-room // per))
            hits += take * self.hit[position]
            predicted += take * per
            expected += take * self.exp[position]
        if not predicted or expected <= 0:
            return float("nan")
        return hits / expected

    def jackknife(self, budget: int) -> List[float]:
        """Leave-one-out budget enrichment, indexed by *original* index.

        Dropping an association below the cutoff changes nothing. Dropping one
        above it shifts everything after up by ``n_predicted``, which is a single
        binary search on the prefix sums rather than a re-walk.
        """
        import bisect

        full = self.enrichment(budget)
        k = self.cutoff(budget)
        values = [full] * self.n
        for position in range(k):
            index = self.order[position]
            per = self.pred[position]
            last = min(bisect.bisect_left(self.cum_pred, budget + per), self.n - 1)
            hits = self.cum_hit[last] - self.hit[position]
            expected = self.cum_exp[last] - self.exp[position]
            values[index] = hits / expected if expected > 0 else float("nan")
        return values


@dataclass(frozen=True)
class TieAmbiguity:
    """How much of a budget slice is decided by an arbitrary tie-break.

    A ranking that assigns the same score to thousands of associations is not
    ranking them; whichever of them lands inside the budget is decided by input
    order. Because equal scores are contiguous in the sorted order, the only
    score that can straddle the cutoff is the one at the boundary — so this is
    exact, not an estimate.
    """

    ranking: str
    budget: int
    n_associations: int
    n_ambiguous: int
    n_tied_pool: int
    boundary_score: float

    @property
    def share(self) -> float:
        return self.n_ambiguous / self.n_associations if self.n_associations else 0.0


def tie_ambiguity(
    outcomes: Sequence[AssociationOutcome], ranking: str, budget: int
) -> TieAmbiguity:
    """Count the associations whose place in a budget slice is a coin toss."""
    ranked = _RankedPool(outcomes, ranking)
    k = ranked.cutoff(budget)
    if k == 0:
        return TieAmbiguity(ranking, budget, 0, 0, 0, float("nan"))

    def score(i: int) -> float:
        return outcomes[i].rank_scores.get(ranking, float("-inf"))

    boundary = score(ranked.order[k - 1])
    inside = sum(1 for i in ranked.order[:k] if score(i) == boundary)
    outside = sum(1 for i in ranked.order[k:] if score(i) == boundary)
    return TieAmbiguity(
        ranking=ranking,
        budget=budget,
        n_associations=k,
        n_ambiguous=inside if outside else 0,
        n_tied_pool=inside + outside,
        boundary_score=boundary,
    )


@dataclass(frozen=True)
class TieBreakSpread:
    """The statistic recomputed under random orderings of tied score blocks.

    This is not a bootstrap and must not be read as one: the data never changes.
    It asks a different and more basic question — *is this statistic a property
    of the ranking at all, or of the order the input file happened to be in?* If
    the observed value sits at an extreme of this spread, the ranking is not what
    produced it, and no confidence interval around it means anything.
    """

    n_shuffles: int
    observed: float
    median: float
    low: float
    high: float
    fraction_positive: float
    observed_percentile: float

    @property
    def observed_is_extreme(self) -> bool:
        """True when the input-order value sits outside the middle 90%."""
        return self.observed_percentile < 0.05 or self.observed_percentile > 0.95


def tie_break_spread(
    outcomes: Sequence[AssociationOutcome],
    ranking_a: str,
    ranking_b: str,
    budget: int,
    n_shuffles: int = 500,
    seed: int = 0,
) -> TieBreakSpread:
    """Re-break tied scores at random and watch the enrichment difference move.

    One random key per association is drawn per shuffle and used as the secondary
    sort key for *both* rankings, so the comparison stays paired: a shuffle that
    happens to favour one ranking's tied block is the same shuffle both see.
    """
    observed = _RankedPool(outcomes, ranking_a).enrichment(budget) - _RankedPool(
        outcomes, ranking_b
    ).enrichment(budget)
    rng = random.Random(seed)
    n = len(outcomes)
    values: List[float] = []
    for _ in range(n_shuffles):
        jitter = [rng.random() for _ in range(n)]
        a = _RankedPool(outcomes, ranking_a, jitter).enrichment(budget)
        b = _RankedPool(outcomes, ranking_b, jitter).enrichment(budget)
        if math.isfinite(a) and math.isfinite(b):
            values.append(a - b)
    if not values:
        nan = float("nan")
        return TieBreakSpread(n_shuffles, observed, nan, nan, nan, nan, nan)
    values.sort()
    below = sum(1 for v in values if v < observed)
    return TieBreakSpread(
        n_shuffles=len(values),
        observed=observed,
        median=quantile(values, 0.5),
        low=quantile(values, 0.025),
        high=quantile(values, 0.975),
        fraction_positive=sum(1 for v in values if v > 0) / len(values),
        observed_percentile=below / len(values),
    )


def compare_rankings(
    outcomes: Sequence[AssociationOutcome],
    ranking_a: str,
    ranking_b: str,
    budget: int,
    n_resamples: int = 1000,
    seed: int = 0,
    confidence: float = 0.95,
    n_tie_shuffles: int = 500,
) -> RankingComparison:
    """Paired bootstrap of the enrichment *difference* between two rankings.

    Comparing two independently-computed confidence intervals is the wrong test
    here: both rankings are applied to the *same* candidate pool, so their
    estimates are correlated and overlapping intervals do not imply no
    difference. Every resample below is drawn once and scored both ways, so the
    difference is measured on identical data every time.

    Both designs described on :class:`RankingComparison` are run from that same
    draw, and each is summarised with percentile, basic and BCa intervals plus
    the diagnostics needed to decide whether to believe any of them.
    """
    n = len(outcomes)
    pool_a = _RankedPool(outcomes, ranking_a)
    pool_b = _RankedPool(outcomes, ranking_b)
    point_a = pool_a.enrichment(budget)
    point_b = pool_b.enrichment(budget)
    point = point_a - point_b

    slice_a = pool_a.slice_indices(budget)
    slice_b = pool_b.slice_indices(budget)
    in_a = [False] * n
    in_b = [False] * n
    for i in slice_a:
        in_a[i] = True
    for i in slice_b:
        in_b[i] = True

    def _leverage(indices: Sequence[int]) -> float:
        total = sum(outcomes[i].n_predicted for i in indices)
        if not total:
            return float("nan")
        return max(outcomes[i].n_predicted for i in indices) / total

    # Fixed-selection sums, and the leave-one-out versions for BCa.
    hits_a = sum(outcomes[i].n_hit for i in slice_a)
    exp_a = sum(outcomes[i].base_rate * outcomes[i].n_predicted for i in slice_a)
    hits_b = sum(outcomes[i].n_hit for i in slice_b)
    exp_b = sum(outcomes[i].base_rate * outcomes[i].n_predicted for i in slice_b)

    def _fixed_ratio(h: float, e: float) -> float:
        return h / e if e > 0 else float("nan")

    fixed_point = _fixed_ratio(hits_a, exp_a) - _fixed_ratio(hits_b, exp_b)
    fixed_jack: List[float] = []
    for i in range(n):
        o = outcomes[i]
        ha, ea = hits_a, exp_a
        hb, eb = hits_b, exp_b
        if in_a[i]:
            ha -= o.n_hit
            ea -= o.base_rate * o.n_predicted
        if in_b[i]:
            hb -= o.n_hit
            eb -= o.base_rate * o.n_predicted
        fixed_jack.append(_fixed_ratio(ha, ea) - _fixed_ratio(hb, eb))

    jack_a = pool_a.jackknife(budget)
    jack_b = pool_b.jackknife(budget)
    reselect_jack = [a - b for a, b in zip(jack_a, jack_b)]

    rng = random.Random(seed)
    reselect_diffs: List[float] = []
    fixed_diffs: List[float] = []
    population = range(n)
    for _ in range(n_resamples):
        counts = [0] * n
        for index in rng.choices(population, k=n):
            counts[index] += 1
        a = pool_a.enrichment_from_counts(counts, budget)
        b = pool_b.enrichment_from_counts(counts, budget)
        if math.isfinite(a) and math.isfinite(b):
            reselect_diffs.append(a - b)
        ha = ea = hb = eb = 0.0
        for i in slice_a:
            c = counts[i]
            if c:
                ha += c * outcomes[i].n_hit
                ea += c * outcomes[i].base_rate * outcomes[i].n_predicted
        for i in slice_b:
            c = counts[i]
            if c:
                hb += c * outcomes[i].n_hit
                eb += c * outcomes[i].base_rate * outcomes[i].n_predicted
        if ea > 0 and eb > 0:
            fixed_diffs.append(ha / ea - hb / eb)

    notes: List[str] = []
    if min(len(slice_a), len(slice_b)) < RankingComparison.MIN_ASSOCIATIONS:
        thin = ranking_a if len(slice_a) <= len(slice_b) else ranking_b
        notes.append(
            f"only {min(len(slice_a), len(slice_b))} {thin} associations fit the "
            f"{budget:,}-prediction budget: too few to resample meaningfully"
        )
    spread = tie_break_spread(
        outcomes, ranking_a, ranking_b, budget, n_tie_shuffles, seed
    )
    if spread.observed_is_extreme:
        notes.append(
            f"under random tie-breaks the difference is "
            f"{spread.median:+.2f} [{spread.low:+.2f}, {spread.high:+.2f}]; the "
            f"observed {spread.observed:+.2f} sits at the "
            f"{100 * spread.observed_percentile:.0f}th percentile of that spread, "
            f"so it is largely an artefact of input order"
        )
    tie_a = tie_ambiguity(outcomes, ranking_a, budget)
    tie_b = tie_ambiguity(outcomes, ranking_b, budget)
    for tie in (tie_a, tie_b):
        if tie.n_ambiguous:
            notes.append(
                f"{tie.n_ambiguous:,}/{tie.n_associations:,} of the {tie.ranking} "
                f"slice sit in a {tie.n_tied_pool:,}-way tie at score "
                f"{tie.boundary_score:.3f}; which of them lands inside the budget "
                f"is input order, not ranking"
            )

    return RankingComparison(
        budget=budget,
        ranking_a=ranking_a,
        ranking_b=ranking_b,
        enrichment_a=point_a,
        enrichment_b=point_b,
        difference=point,
        n_associations_a=len(slice_a),
        n_associations_b=len(slice_b),
        reselect=summarise_bootstrap(
            point, reselect_diffs, reselect_jack, confidence=confidence
        ),
        fixed=summarise_bootstrap(
            fixed_point, fixed_diffs, fixed_jack, confidence=confidence
        ),
        max_leverage_a=_leverage(slice_a),
        max_leverage_b=_leverage(slice_b),
        tie_a=tie_a,
        tie_b=tie_b,
        tie_break=spread,
        notes=tuple(notes),
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
    # the file just read on every replay.
    #
    # Everything is written at 12 significant figures so that a replay is
    # *lossless*. It used not to be: base_rate was stored to 5 decimal places,
    # which is one or two significant figures for the smallest acquisition rates
    # (1/19,000 ≈ 5.3e-5), and dcgo_score to 2. That is why the committed
    # metrics file and a replay of the committed per-association file did not
    # agree — see SURPRISE_SCORE.md. Replay must reproduce the run exactly or
    # the "reproduce with" command in the docs is a fiction.
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
                    f"{o.feature}\t{o.term}\t{o.rank_scores['surprise']:.12g}\t"
                    f"{o.rank_scores['dcgo']:.12g}\t{o.n_predicted}\t{o.n_hit}\t"
                    f"{o.base_rate:.12g}\t{o.enrichment:.6g}\n"
                )
        logger.info(
            f"✓ Wrote per-association outcomes to {args.per_association_output}"
        )

    comparisons = [
        compare_rankings(
            scored,
            "surprise",
            "dcgo",
            budget,
            args.bootstrap,
            args.seed,
            n_tie_shuffles=args.tie_shuffles,
        )
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
    logger.info("Paired head-to-head (one draw per resample, scored both ways).")
    logger.info(
        "Two designs: 'reselect' re-fills the budget inside every resample; "
        "'fixed' conditions on the observed slices and only re-weights."
    )
    for c in comparisons:
        logger.info("")
        logger.info(
            f"  @{c.budget:>6,} preds: surprise {c.enrichment_a:6.2f} "
            f"({c.n_associations_a:,} assoc) vs dcgo {c.enrichment_b:6.2f} "
            f"({c.n_associations_b:,} assoc)  diff {c.difference:+6.2f}  "
            f"→ {c.verdict}"
        )
        for design, s in (("reselect", c.reselect), ("fixed", c.fixed)):
            logger.info(
                f"      {design:<8} pct [{s.percentile[0]:+7.2f},{s.percentile[1]:+7.2f}]  "
                f"basic [{s.basic[0]:+7.2f},{s.basic[1]:+7.2f}]  "
                f"BCa [{s.bca[0]:+7.2f},{s.bca[1]:+7.2f}]  "
                f"z0={s.z0:+.2f} a={s.acceleration:+.3f} skew={s.skew:+.2f} "
                f"favours-A={100 * s.fraction_positive:.0f}%"
            )
            if not s.trustworthy:
                logger.info(f"      {'':<8} ! {s.note}")
        if c.tie_break is not None:
            t = c.tie_break
            logger.info(
                f"      {'tie-break':<8} random re-orderings of tied scores give "
                f"{t.median:+.2f} [{t.low:+.2f}, {t.high:+.2f}] "
                f"(observed {t.observed:+.2f} = {100 * t.observed_percentile:.0f}th "
                f"pct, favours surprise in {100 * t.fraction_positive:.0f}%)"
            )
        for note in c.notes:
            logger.info(f"      ! {note}")

    with open(args.output, "a") as f:
        f.write("\n# paired ranking comparison (surprise - dcgo)\n")
        f.write(
            "# design=reselect: the budget is re-filled from scratch inside every\n"
            "#   resample, so the composition of the top of each ranking is\n"
            "#   resampled too. design=fixed: the two budget slices are held at the\n"
            "#   ones the observed ranking produced and the resample only re-weights\n"
            "#   associations (shared members move together, so it stays paired).\n"
            "# 'point' is the statistic each design's replicates are compared to.\n"
            "# trustworthy=False means the bootstrap on that row cannot support an\n"
            "#   interval. resolvable=False is stronger: the comparison itself is\n"
            "#   not answerable at this budget (too few associations fit, or most of\n"
            "#   a slice is an arbitrary tie-break) — see 'verdict', and do not quote\n"
            "#   any interval from that budget however trustworthy the row looks.\n"
        )
        f.write(
            "budget\tdesign\tn_assoc_surprise\tn_assoc_dcgo\tenrichment_surprise\t"
            "enrichment_dcgo\tpoint\tpct_low\tpct_high\tbasic_low\tbasic_high\t"
            "bca_low\tbca_high\tboot_mean\tboot_median\tboot_sd\tboot_skew\tz0\t"
            "acceleration\tfraction_favouring_surprise\tn_resamples\ttrustworthy\t"
            "resolvable\tverdict\tnote\n"
        )
        for c in comparisons:
            for design, s in (("reselect", c.reselect), ("fixed", c.fixed)):
                f.write(
                    f"{c.budget}\t{design}\t{c.n_associations_a}\t"
                    f"{c.n_associations_b}\t{c.enrichment_a:.2f}\t"
                    f"{c.enrichment_b:.2f}\t{s.point:.2f}\t"
                    f"{s.percentile[0]:.2f}\t{s.percentile[1]:.2f}\t"
                    f"{s.basic[0]:.2f}\t{s.basic[1]:.2f}\t"
                    f"{s.bca[0]:.2f}\t{s.bca[1]:.2f}\t"
                    f"{s.mean:.2f}\t{s.median:.2f}\t{s.sd:.2f}\t{s.skew:.2f}\t"
                    f"{s.z0:.3f}\t{s.acceleration:.4f}\t"
                    f"{s.fraction_positive:.3f}\t{s.n_usable}\t"
                    f"{str(s.trustworthy).lower()}\t{str(c.resolvable).lower()}\t"
                    f"{c.verdict}\t{s.note}\n"
                )

        f.write("\n# tie-break sensitivity of the same difference\n")
        f.write(
            "# Not a bootstrap: the data never changes, only the order of\n"
            "#   associations that share a rank score. If 'observed' sits at an\n"
            "#   extreme percentile of this spread, the statistic is a property of\n"
            "#   the input file's order rather than of the ranking.\n"
        )
        f.write(
            "budget\tn_shuffles\tobserved\tmedian\tlow\thigh\t"
            "observed_percentile\tfraction_favouring_surprise\t"
            "ambiguous_surprise\tslice_surprise\tambiguous_dcgo\tslice_dcgo\n"
        )
        for c in comparisons:
            t = c.tie_break
            if t is None:
                continue
            ta, tb = c.tie_a, c.tie_b
            f.write(
                f"{c.budget}\t{t.n_shuffles}\t{t.observed:.2f}\t{t.median:.2f}\t"
                f"{t.low:.2f}\t{t.high:.2f}\t{t.observed_percentile:.3f}\t"
                f"{t.fraction_positive:.3f}\t"
                f"{ta.n_ambiguous if ta else ''}\t{ta.n_associations if ta else ''}\t"
                f"{tb.n_ambiguous if tb else ''}\t{tb.n_associations if tb else ''}\n"
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
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Bootstrap resamples. The budget-matched comparison needs several "
        "thousand before its BCa tails settle; 1,000 is enough for the "
        "fixed-stratum CIs only",
    )
    parser.add_argument(
        "--tie-shuffles",
        type=int,
        default=500,
        help="Random re-orderings of tied rank scores, used to show how much of "
        "a budget-matched result is input order rather than ranking",
    )
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
