"""The surprise score: ranking *emergent* domain-combination predictions.

dcGO's distinctive claim is that a **combination** of domains can predict a term
that none of its constituents predicts alone — a signal homology transfer and
single-domain methods cannot see. The significant-association table contains
thousands of such combinations, but raw significance is the wrong way to rank
them: the strongest ones are usually artefacts (two InterPro signatures
describing *one* region) or restatements of curated knowledge, while the genuinely
novel ones sit on two or three supporting proteins.

The **surprise score** separates those. For a supra-domain ``S`` (the contiguous
combination ``d1,d2[,d3]``) associated with term ``t``, it is a product of three
interpretable factors::

    surprise = -log10(q_emergence) × distinctness × novelty

**1. Emergence** — the statistical core. Let ``rate(f) = P(t | protein has f)``.
The parts-only expectation combines what the pieces already tell you:

* a **noisy-OR** over the constituent single domains,
  ``1 - Π(1 - rate(d_i))`` — the rate you would predict if each domain
  contributed to ``t`` independently;
* floored by the best **proper sub-feature** rate (a triple whose contained pair
  already predicts ``t`` is not surprising);
* floored by the term's **background rate** in the protein universe (a rate no
  worse than picking proteins at random is not evidence of anything, and this
  keeps the null strictly positive).

The observed count is then tested against that expectation with a one-sided
binomial tail — "how unlikely is it that ``a`` of the ``n`` proteins carrying
``S`` are annotated ``t``, if ``S`` were no better than its parts?" —
and Benjamini–Hochberg corrected across all candidates. Support is handled
honestly by construction: two supporting proteins can only reach significance
when the parts-only expectation is genuinely tiny.

**2. Distinctness** — ``1 - overlap``, where ``overlap`` is the median (over
supporting proteins) largest pairwise overlap between the constituent domains'
matched regions. Two redundant signatures for the same region (the classic
"GPCR rhodopsin-like + 7TM" artefact) overlap almost completely and score ~0;
genuinely separate domains score ~1.

**3. Novelty** — how much of the prediction curators already record for the
*constituent* domains, from a curated domain→term reference (InterPro2GO for
GO). Terms already mapped to a constituent are discounted hard; terms more
specific than curated knowledge are discounted mildly; terms outside it keep the
full weight. Without a curated reference the factor is 1.0 and the status is
``"no-reference"``.

Every component is reported alongside the score, so a consumer can re-rank or
re-weight without recomputation. The driver is
``scripts/rank_surprising_associations.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from scipy.stats import binom

#: Weight applied to a term the curated reference already maps to a constituent.
NOVELTY_CURATED = 0.1
#: Weight for a term more *general* than curated constituent knowledge.
NOVELTY_IMPLIED = 0.3
#: Weight for a term more *specific* than curated constituent knowledge.
NOVELTY_REFINES = 0.6
#: Weight for a term the curated reference does not reach at all.
NOVELTY_NOVEL = 1.0


@dataclass(frozen=True)
class EmergenceEvidence:
    """Counts behind one candidate ``(supra-domain, term)`` pair.

    Attributes:
        feature: the supra-domain id, comma-joined InterPro entries in
            positional order (``"IPR000198,IPR001452"``).
        term: the ontology term the combination is associated with.
        n_feature: proteins carrying the combination.
        n_both: of those, how many are annotated ``term``.
        single_rates: ``rate(t | d_i)`` per constituent single domain.
        part_rates: ``rate(t | sub-feature)`` for proper sub-features that exist
            in the analysis (for a triple: its two contained pairs).
        background_rate: ``rate(t)`` across the whole protein universe.
        q_value: the FDR-corrected p-value of the original dcGO association,
            carried through for reference.
    """

    feature: str
    term: str
    n_feature: int
    n_both: int
    single_rates: Tuple[float, ...]
    part_rates: Tuple[float, ...]
    background_rate: float
    q_value: float


@dataclass(frozen=True)
class SurpriseResult:
    """A scored candidate: the surprise score and every component behind it."""

    feature: str
    term: str
    n_feature: int
    n_both: int
    observed_rate: float
    expected_rate: float
    expectation_source: str
    lift: float
    p_emergence: float
    q_emergence: float = float("nan")
    region_overlap: float = 0.0
    distinctness: float = 1.0
    novelty: float = NOVELTY_NOVEL
    novelty_status: str = "no-reference"
    #: constituents that, alone, predict the term no better than chance
    uninformative_constituents: int = 0
    q_value: float = float("nan")

    @property
    def surprise(self) -> float:
        """``-log10(q_emergence) × distinctness × novelty``, floored at 0."""
        if not math.isfinite(self.q_emergence) or self.q_emergence <= 0:
            # q == 0 only from underflow; treat as the smallest representable.
            strength = -math.log10(5e-324)
        else:
            strength = -math.log10(self.q_emergence)
        return max(0.0, strength * self.distinctness * self.novelty)


def conditional_rate(
    n_annotated: int,
    n_carriers: int,
    background_rate: float,
    pseudo_count: float = 1.0,
) -> float:
    """``P(term | feature)``, shrunk toward the term's background rate.

    The raw fraction is badly overconfident for the rare features this analysis
    is full of: a domain seen in three proteins that all carry the term gets
    ``rate = 1.0``, which would then declare *every* combination containing it
    unsurprising. Adding ``pseudo_count`` observations drawn from the background
    turns that into ``(3 + 0.001) / 4 ≈ 0.75`` — still high, but no longer a
    certainty inferred from three proteins.

    Args:
        n_annotated: carriers of the feature also annotated with the term.
        n_carriers: proteins carrying the feature.
        background_rate: the term's rate across the whole protein universe.
        pseudo_count: strength of the shrinkage, in pseudo-observations.
    """
    if n_carriers <= 0:
        return background_rate
    return (n_annotated + pseudo_count * background_rate) / (n_carriers + pseudo_count)


def noisy_or(rates: Iterable[float]) -> float:
    """Independent-contributions combination of per-part rates.

    ``1 - Π(1 - rate)``: the probability at least one part "delivers" the term,
    if the parts acted independently. Empty input → 0.0.
    """
    product = 1.0
    for rate in rates:
        product *= 1.0 - min(max(rate, 0.0), 1.0)
    return 1.0 - product


def expected_rate(evidence: EmergenceEvidence) -> Tuple[float, str]:
    """The parts-only expectation for a combination, and where it came from.

    The maximum of the noisy-OR over constituents, the best proper sub-feature
    rate, and the term's background rate — see the module docstring. Returns
    ``(rate, source)`` where ``source`` is ``"noisy_or"``, ``"best_part"`` or
    ``"background"``.
    """
    candidates = [
        (noisy_or(evidence.single_rates), "noisy_or"),
        (max(evidence.part_rates, default=0.0), "best_part"),
        (evidence.background_rate, "background"),
    ]
    rate, source = max(candidates, key=lambda pair: pair[0])
    # Clamp strictly inside (0, 1) so the binomial tail stays defined.
    return min(max(rate, 1e-12), 1.0 - 1e-12), source


def emergence_pvalue(n_both: int, n_feature: int, rate: float) -> float:
    """One-sided binomial tail ``P(X >= n_both)`` for ``X ~ Bin(n_feature, rate)``.

    The probability of seeing at least this many annotated proteins among those
    carrying the combination, if the combination were no better than its parts.
    """
    if n_feature <= 0 or n_both <= 0:
        return 1.0
    return float(binom.sf(n_both - 1, n_feature, rate))


def overlap_fraction(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    """Overlap of two closed intervals as a fraction of the shorter one."""
    (a_start, a_end), (b_start, b_end) = a, b
    overlap = min(a_end, b_end) - max(a_start, b_start) + 1
    if overlap <= 0:
        return 0.0
    shorter = min(a_end - a_start + 1, b_end - b_start + 1)
    return overlap / shorter if shorter > 0 else 0.0


def max_pairwise_overlap(intervals: Sequence[Tuple[int, int]]) -> float:
    """Largest :func:`overlap_fraction` among all pairs of matched regions.

    ~1.0 means the "combination" is really several signatures describing a
    single region — the redundant-signature artefact the score must discount.
    """
    worst = 0.0
    for i in range(len(intervals)):
        for j in range(i + 1, len(intervals)):
            worst = max(worst, overlap_fraction(intervals[i], intervals[j]))
    return worst


def median(values: Sequence[float]) -> float:
    """Median of ``values`` (0.0 when empty)."""
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def locate_feature_regions(
    domain_ids: Sequence[str],
    intervals: Sequence[Tuple[int, int]],
    parts: Sequence[str],
) -> List[Tuple[int, int]]:
    """Matched regions of the first occurrence of a contiguous domain run.

    Supra-domains are built from a protein's positionally sorted domain list, so
    the combination ``("IPR1", "IPR2")`` corresponds to *adjacent* entries in
    that list. Returns the intervals of the first match, or ``[]`` if the run
    does not occur (possible when the run came from a different protein).
    """
    if not parts or len(parts) > len(domain_ids):
        return []
    for start in range(len(domain_ids) - len(parts) + 1):
        if list(domain_ids[start : start + len(parts)]) == list(parts):
            return list(intervals[start : start + len(parts)])
    return []


def novelty_factor(
    term: str,
    curated_terms: Set[str],
    ancestors_fn: Optional[Callable[[str], Iterable[str]]] = None,
) -> Tuple[float, str]:
    """How new the prediction is relative to curated constituent annotations.

    Args:
        term: the predicted term.
        curated_terms: terms a curated reference (e.g. InterPro2GO) already maps
            to *any* constituent domain of the combination. Empty means the
            reference says nothing about the constituents.
        ancestors_fn: term → ancestors, used to tell "more general than curated"
            from "more specific than curated". ``None`` skips that distinction.

    Returns:
        ``(factor, status)``; status is ``"curated"``, ``"implied"``,
        ``"refines"``, ``"novel"`` or ``"no-reference"``.
    """
    if not curated_terms:
        return NOVELTY_NOVEL, "no-reference"
    if term in curated_terms:
        return NOVELTY_CURATED, "curated"
    if ancestors_fn is not None:
        term_ancestors = set(ancestors_fn(term))
        if term_ancestors & curated_terms:
            # The prediction sits below curated knowledge: a refinement.
            return NOVELTY_REFINES, "refines"
        for curated in curated_terms:
            if term in set(ancestors_fn(curated)):
                # The prediction is a generalisation of what curators recorded.
                return NOVELTY_IMPLIED, "implied"
    return NOVELTY_NOVEL, "novel"


def score_candidate(
    evidence: EmergenceEvidence,
    region_overlap: float = 0.0,
    novelty: float = NOVELTY_NOVEL,
    novelty_status: str = "no-reference",
) -> SurpriseResult:
    """Score one candidate. ``q_emergence`` is filled in later, by :func:`apply_fdr`."""
    rate, source = expected_rate(evidence)
    observed = evidence.n_both / evidence.n_feature if evidence.n_feature else 0.0
    return SurpriseResult(
        feature=evidence.feature,
        term=evidence.term,
        n_feature=evidence.n_feature,
        n_both=evidence.n_both,
        observed_rate=observed,
        expected_rate=rate,
        expectation_source=source,
        lift=observed / rate if rate > 0 else float("inf"),
        p_emergence=emergence_pvalue(evidence.n_both, evidence.n_feature, rate),
        region_overlap=region_overlap,
        distinctness=max(0.0, 1.0 - region_overlap),
        novelty=novelty,
        novelty_status=novelty_status,
        # After shrinkage no rate is exactly zero, so "silent" means "no better
        # than the term's background rate".
        uninformative_constituents=sum(
            1 for r in evidence.single_rates if r <= evidence.background_rate
        ),
        q_value=evidence.q_value,
    )


def apply_fdr(
    results: Sequence[SurpriseResult], alpha: float = 0.05
) -> List[SurpriseResult]:
    """Attach Benjamini–Hochberg ``q_emergence`` values to scored candidates.

    The emergence test is a second, independent family of hypotheses (one per
    candidate combination), so it needs its own multiple-testing correction —
    the dcGO q-value in the input table controls a different family.
    """
    # Imported here to keep the module importable without the Cython extension.
    from src.vectorized_fisher import benjamini_hochberg_correction
    import numpy as np

    if not results:
        return []
    adjusted, _threshold = benjamini_hochberg_correction(
        np.array([r.p_emergence for r in results]), alpha=alpha
    )
    return [
        SurpriseResult(**{**vars(result), "q_emergence": float(q)})
        for result, q in zip(results, adjusted)
    ]


def proper_subfeatures(parts: Sequence[str]) -> List[str]:
    """Contiguous proper sub-combinations of a supra-domain, longest first.

    ``("A", "B", "C")`` → ``["A,B", "B,C"]``. Single domains are handled
    separately (they are the noisy-OR inputs), so they are not returned here.
    """
    subfeatures: List[str] = []
    for length in range(len(parts) - 1, 1, -1):
        for start in range(len(parts) - length + 1):
            subfeatures.append(",".join(parts[start : start + length]))
    return subfeatures


def parse_interpro2go(lines: Iterable[str]) -> Dict[str, Set[str]]:
    """Parse the curated InterPro2GO mapping into ``{IPR id: {GO id}}``.

    Line shape (``!`` comments ignored)::

        InterPro:IPR000003 Retinoid X receptor > GO:DNA binding ; GO:0003677
    """
    mapping: Dict[str, Set[str]] = {}
    for line in lines:
        if line.startswith("!"):
            continue
        head, _, go_id = line.rstrip("\n").rpartition(";")
        go_id = go_id.strip()
        if not go_id.startswith("GO:") or not head.startswith("InterPro:"):
            continue
        interpro_id = head[len("InterPro:") :].split()[0]
        mapping.setdefault(interpro_id, set()).add(go_id)
    return mapping
