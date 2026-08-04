#!/usr/bin/env python3
"""Temporal held-out benchmark (CAFA-style) for dcGO — VALIDATION_PLAN.md §2.

This is the precision-capable test the InterPro2GO coverage comparison (§1) could
not provide. The idea, following CAFA:

  * Train domain->GO associations using only annotations available at an early
    GOA release ``t0``.
  * Build a **no-knowledge** benchmark (CAFA-style), per aspect: proteins that
    had no experimental annotation in that aspect at ``t0`` but gained one by a
    later release ``t1``. Each such protein's truth is its *full* propagated
    ``t1`` experimental term set in that aspect — so every scored term is a
    genuine prediction, never memorised from training.
  * Transfer the ``t0`` domain associations to each benchmark protein via its
    domains, then score with the CAFA **protein-centric** metrics: ``F_max``,
    ``S_min`` (information-content weighted) and ``AUPRC``, reported separately
    per GO aspect (BP / MF / CC).
  * Optionally sweep an **information-content floor** (``--min-ic``): excluding
    near-universal low-IC terms (e.g. GO:0005515 "protein binding", ~85% of human
    experimental MF annotations) tests whether the naive baseline's F_max lead is
    just base-rate recovery of generic terms.

Because the domain architectures come from ``protein2ipr`` and are not
time-varying in our data, only the GOA annotations move in time — the split is
purely on the annotation side.

The metric maths lives in small pure functions (unit-tested in
``tests/unit/test_temporal_benchmark.py``); the ``main`` at the bottom wires them
to real files.

Notes / deliberate simplifications (documented so they are not oversold):
  * Information content is the **marginal** IC, ``IC(t) = -log2(P(t))`` estimated
    from the propagated ``t0`` training frequencies. This is the common CAFA
    approximation, not the full Clark & Radivojac information-accretion (which
    conditions on parents). It is enough to weight rare, specific terms above
    generic ones for ``S_min``.
  * The three aspect roots (GO:0008150 / GO:0003674 / GO:0005575) are excluded
    from evaluation — they are trivially true for every protein.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Callable, Iterable, Mapping

# The three GO aspect roots — trivially true, excluded from scoring (CAFA does
# the same). Also used to map a term to its aspect via namespace.
BP_ROOT = "GO:0008150"
MF_ROOT = "GO:0003674"
CC_ROOT = "GO:0005575"
ASPECT_ROOTS = frozenset({BP_ROOT, MF_ROOT, CC_ROOT})

# namespace (obonet) -> short aspect code
NAMESPACE_TO_ASPECT = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

GetAncestors = Callable[[str], set[str]]
# protein -> {go_term: score}
PredScores = Mapping[str, Mapping[str, float]]
# protein -> {go_term}
TrueSets = Mapping[str, set[str]]


# --------------------------------------------------------------------------- #
# Propagation helpers                                                          #
# --------------------------------------------------------------------------- #
def propagate_terms(terms: Iterable[str], get_ancestors: GetAncestors) -> set[str]:
    """Expand a set of GO terms to include all their ancestors (True Path Rule)."""
    out: set[str] = set()
    for t in terms:
        out.add(t)
        out.update(get_ancestors(t))
    return out


def build_nk_benchmark_by_aspect(
    t0_known_map: Mapping[str, set[str]],
    t1_exp_map: Mapping[str, set[str]],
    term_aspect: Mapping[str, str],
    get_ancestors: GetAncestors,
    predictable_proteins: set[str] | None = None,
) -> dict[str, dict[str, set[str]]]:
    """CAFA **no-knowledge** benchmark, split by GO aspect.

    Returns ``{aspect: {protein: true_terms}}``. For a given aspect a protein is a
    benchmark ("no-knowledge") target only if it had **no annotation known at t0
    in that aspect** but gained experimental annotation by t1 — so every scored
    term is a genuine prediction, never available to training. Its truth set is
    then the **full** propagated t1 experimental terms in that aspect (aspect
    roots excluded), which is the correct scoring reference: a method must not be
    penalised for correctly predicting a term the protein really has. (An earlier
    delta-only design, ``t1 minus t0``, wrongly scored correct predictions of
    already-known terms as misinformation.)

    ``t0_known_map`` must be the annotations available to *training* — i.e. under
    the **same evidence filter the pipeline trained on** (``manual``/non-IEA),
    not experimental-only. Gating on a narrower set than training would leak: a
    protein with a computational (e.g. ISS/IBA) t0 label the model already saw
    could re-enter the held-out set and have that label's later experimental
    confirmation scored as a fresh prediction. ``t1_exp_map`` holds t1
    **experimental** annotations (the gold standard). ``predictable_proteins``
    (if given) keeps only proteins with at least one domain.
    """
    aspects = ("BP", "MF", "CC")
    benchmark: dict[str, dict[str, set[str]]] = {a: {} for a in aspects}
    for protein, t1_terms in t1_exp_map.items():
        if predictable_proteins is not None and protein not in predictable_proteins:
            continue
        t1_closed = propagate_terms(t1_terms, get_ancestors)
        t0_closed = propagate_terms(t0_known_map.get(protein, set()), get_ancestors)
        t0_aspects = {term_aspect.get(t) for t in t0_closed}
        for aspect in aspects:
            if aspect in t0_aspects:
                # Already had experimental knowledge here at t0 — not no-knowledge.
                continue
            targets = {
                t for t in t1_closed if term_aspect.get(t) == aspect
            } - ASPECT_ROOTS
            if targets:
                benchmark[aspect][protein] = targets
    return benchmark


# --------------------------------------------------------------------------- #
# Information content                                                          #
# --------------------------------------------------------------------------- #
def information_content(
    annotation_map: Mapping[str, set[str]], get_ancestors: GetAncestors
) -> dict[str, float]:
    """Marginal information content ``IC(t) = -log2(P(t))`` from training frequencies.

    ``P(t)`` is the fraction of annotated proteins carrying term ``t`` after
    propagating each protein's annotations to its ancestor closure. Terms never
    seen get IC 0 (they carry no information for weighting). The estimate is over
    the set of annotated proteins in ``annotation_map``.
    """
    n_proteins = len(annotation_map)
    if n_proteins == 0:
        return {}
    counts: dict[str, int] = defaultdict(int)
    for terms in annotation_map.values():
        for t in propagate_terms(terms, get_ancestors):
            counts[t] += 1
    ic: dict[str, float] = {}
    for t, c in counts.items():
        p = c / n_proteins
        ic[t] = -math.log2(p) if 0.0 < p < 1.0 else 0.0
    return ic


# --------------------------------------------------------------------------- #
# Prediction transfer                                                          #
# --------------------------------------------------------------------------- #
def transfer_predictions(
    protein_domains: Mapping[str, Iterable[str]],
    domain_go_scores: Mapping[str, Mapping[str, float]],
    get_ancestors: GetAncestors,
) -> dict[str, dict[str, float]]:
    """Transfer domain->GO associations to proteins, max-propagated up the DAG.

    A protein's predicted GO terms are the union over its domains of that
    domain's associated terms; the score for a (protein, term) is the **max**
    association score across contributing domains. Each predicted term is then
    propagated to its ancestors, an ancestor taking the max score of any
    descendant that implies it (standard CAFA score propagation). Aspect roots
    are dropped.
    """
    predictions: dict[str, dict[str, float]] = {}
    for protein, domains in protein_domains.items():
        scores: dict[str, float] = {}
        for domain in domains:
            for go, s in domain_go_scores.get(domain, {}).items():
                # propagate the raw term and its ancestors with the same score,
                # keeping the max seen for each term.
                for term in (go, *get_ancestors(go)):
                    if term in ASPECT_ROOTS:
                        continue
                    if s > scores.get(term, float("-inf")):
                        scores[term] = s
        if scores:
            predictions[protein] = scores
    return predictions


def transfer_predictions_pscore(
    protein_domains: Mapping[str, Iterable[str]],
    domain_go_scores: Mapping[str, Mapping[str, float]],
    get_ancestors: GetAncestors,
) -> dict[str, dict[str, float]]:
    """Transfer via the original dcGO Predictor **p-score** (Fang & Gough 2013).

    Two differences from :func:`transfer_predictions`: the per-(protein, term)
    score is the **sum** of contributing association scores (additive evidence
    across the protein's domains, propagated to ancestors), and each protein's
    scores are then **min-max normalised to [0, 1]** — a per-protein calibration
    so the threshold sweep ranks each protein's own predictions relative to each
    other. Aspect roots are dropped. Constant vectors normalise to all-1.
    """
    predictions: dict[str, dict[str, float]] = {}
    for protein, domains in protein_domains.items():
        sums: dict[str, float] = defaultdict(float)
        for domain in domains:
            for go, s in domain_go_scores.get(domain, {}).items():
                for term in (go, *get_ancestors(go)):
                    if term in ASPECT_ROOTS:
                        continue
                    sums[term] += s
        if not sums:
            continue
        lo, hi = min(sums.values()), max(sums.values())
        span = hi - lo
        # min-max to [0, 1]; a single/constant score maps to 1.0 (fully supported).
        predictions[protein] = {
            t: ((v - lo) / span if span > 0 else 1.0) for t, v in sums.items()
        }
    return predictions


def naive_predictions(
    eval_proteins: Iterable[str],
    term_frequency: Mapping[str, float],
) -> dict[str, dict[str, float]]:
    """CAFA ``Naive`` baseline: predict every term for every protein at its freq.

    ``term_frequency`` maps a GO term to its propagated training frequency in
    [0, 1]; every evaluation protein gets the identical prediction vector. F_max
    of this baseline is the floor any real method must clear.
    """
    base = {t: f for t, f in term_frequency.items() if t not in ASPECT_ROOTS}
    return {p: dict(base) for p in eval_proteins}


def shuffle_domain_go(
    domain_go_scores: Mapping[str, Mapping[str, float]],
    seed: int = 0,
) -> dict[str, dict[str, float]]:
    """Random-domain null: permute which GO-set each domain owns.

    Each domain is reassigned another domain's entire (term -> score) association
    set by a seeded permutation of the domain labels. Run once this is an
    anecdote; :func:`permutation_null_seeds` turns it into a distribution.

    **What is exchangeable, exactly.** The null hypothesis is *"a domain's
    identity carries no information about which functions its carrier proteins
    have"*. Permuting the labels of the domain -> GO map, and nothing else, is the
    randomisation that hypothesis licenses. Note what it preserves and what it
    destroys, because both matter for how hard the null is to beat:

    * **Preserved.** Every protein keeps its real architecture. The multiset of
      GO-sets is unchanged, so the marginal frequency of each GO term across the
      association table is exactly preserved — a permuted table still predicts
      ``protein binding`` about as often as the real one. The null is therefore
      *base-rate preserving*, which is why it scores in the same range as the
      naive baseline at IC >= 0 rather than at zero. The number of domains that
      make any prediction, and the size distribution of their term sets, are also
      preserved.
    * **Destroyed.** The pairing between a domain and its terms, and with it the
      correlation between a domain's *prevalence* and the *size/specificity* of
      its term set: a rare domain can inherit a promiscuous domain's 400-term
      set.

    **What this null does NOT test, stated plainly.** It permutes the *surviving,
    FDR-significant* association table, so it inherits the real pipeline's
    decisions about which domains are predictive at all and how many terms each
    gets. It is a null for the **transfer step**, not for the whole method: it
    cannot say whether the Fisher + BH stage itself is calibrated. The stronger
    null — permute the protein -> GO labels and re-run inference end to end — is
    still open (see VALIDATION_PLAN §2 "Baselines"), and would be a much harder
    null to beat because a permuted training set would yield far fewer
    significant associations rather than the same number of scrambled ones.
    """
    import numpy as np

    domains = list(domain_go_scores.keys())
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(domains))
    return {
        domains[i]: dict(domain_go_scores[domains[perm[i]]])
        for i in range(len(domains))
    }


def permutation_null_seeds(base_seed: int, n_permutations: int) -> list[int]:
    """The seeds of a permutation null, so it is reproducible and auditable.

    Seed ``i`` is ``base_seed + i``, which means permutation 0 is exactly the
    single shuffle the harness used to report — the old result is now the first
    draw of a distribution rather than the whole story.
    """
    if n_permutations < 1:
        raise ValueError(f"n_permutations must be >= 1, got {n_permutations}")
    return [base_seed + i for i in range(n_permutations)]


# --------------------------------------------------------------------------- #
# CAFA protein-centric metrics                                                 #
# --------------------------------------------------------------------------- #
def precision_recall_at_threshold(
    pred_scores: PredScores, true_sets: TrueSets, tau: float
) -> tuple[float, float, int]:
    """Average precision/recall at score threshold ``tau`` (CAFA convention).

    Precision is averaged only over proteins with at least one prediction at or
    above ``tau`` (count ``m``); recall is averaged over **all** benchmark
    proteins in ``true_sets``. Returns ``(precision, recall, m)``.
    """
    n_eval = len(true_sets)
    if n_eval == 0:
        return 0.0, 0.0, 0
    prec_sum = 0.0
    rec_sum = 0.0
    m = 0
    for protein, true_terms in true_sets.items():
        if not true_terms:
            continue
        predicted = {t for t, s in pred_scores.get(protein, {}).items() if s >= tau}
        if predicted:
            m += 1
            tp = len(predicted & true_terms)
            prec_sum += tp / len(predicted)
            rec_sum += tp / len(true_terms)
        # proteins with no prediction contribute 0 recall (and are excluded from
        # the precision average, per CAFA).
    precision = prec_sum / m if m else 0.0
    recall = rec_sum / n_eval
    return precision, recall, m


def _candidate_thresholds(pred_scores: PredScores, n_points: int = 51) -> list[float]:
    """Threshold sweep at *observed* score cutoffs, plus a predict-nothing endpoint.

    Cutoffs are drawn from the score distribution (quantiles of the observed
    scores), not spaced evenly over the value range: for skewed scores like
    ``-log10(p)`` a value-linspace wastes points in the empty high range and can
    miss the cutoff that separates a useful term from a slightly lower false one.
    A sentinel strictly above the maximum is always appended so ``S_min`` (and
    ``F_max``) can evaluate "predict nothing" when every prediction is a high
    false positive.
    """
    all_scores = [s for terms in pred_scores.values() for s in terms.values()]
    if not all_scores:
        return [0.0]
    hi = max(all_scores)
    sentinel = hi + (abs(hi) if hi != 0 else 1.0)  # strictly above max => empty set
    uniq = sorted(set(all_scores))
    if len(uniq) <= n_points:
        return uniq + [sentinel]
    # Quantile sample of the observed scores: each cutoff is a real score, placed
    # where the mass is. Keeps the sweep cheap while tracking the distribution.
    xs = sorted(all_scores)
    step = (len(xs) - 1) / (n_points - 1)
    sampled = sorted({xs[round(i * step)] for i in range(n_points)})
    return sampled + [sentinel]


def pr_curve(
    pred_scores: PredScores,
    true_sets: TrueSets,
    thresholds: Iterable[float] | None = None,
) -> list[tuple[float, float, float]]:
    """(threshold, precision, recall) points across the score sweep."""
    taus = (
        list(thresholds)
        if thresholds is not None
        else _candidate_thresholds(pred_scores)
    )
    curve = []
    for tau in taus:
        p, r, _ = precision_recall_at_threshold(pred_scores, true_sets, tau)
        curve.append((tau, p, r))
    return curve


def f_max(
    pred_scores: PredScores,
    true_sets: TrueSets,
    thresholds: Iterable[float] | None = None,
) -> tuple[float, float]:
    """CAFA ``F_max``: the best harmonic mean of precision/recall over the sweep.

    Returns ``(f_max, tau_star)``.
    """
    best_f = 0.0
    best_tau = 0.0
    for tau, p, r in pr_curve(pred_scores, true_sets, thresholds):
        if p + r > 0:
            f = 2 * p * r / (p + r)
            if f > best_f:
                best_f, best_tau = f, tau
    return best_f, best_tau


def s_min(
    pred_scores: PredScores,
    true_sets: TrueSets,
    ic: Mapping[str, float],
    thresholds: Iterable[float] | None = None,
) -> tuple[float, float]:
    """CAFA ``S_min``: minimum semantic distance ``sqrt(ru^2 + mi^2)`` over the sweep.

    ``ru`` (remaining uncertainty) sums the IC of missed true terms; ``mi``
    (misinformation) sums the IC of predicted-but-wrong terms; both averaged over
    all benchmark proteins. Returns ``(s_min, tau_star)``. Lower is better.
    """
    n_eval = len(true_sets)
    if n_eval == 0:
        return 0.0, 0.0
    taus = (
        list(thresholds)
        if thresholds is not None
        else _candidate_thresholds(pred_scores)
    )
    best_s = float("inf")
    best_tau = 0.0
    for tau in taus:
        ru_sum = 0.0
        mi_sum = 0.0
        for protein, true_terms in true_sets.items():
            predicted = {t for t, s in pred_scores.get(protein, {}).items() if s >= tau}
            for t in true_terms - predicted:
                ru_sum += ic.get(t, 0.0)
            for t in predicted - true_terms:
                mi_sum += ic.get(t, 0.0)
        ru = ru_sum / n_eval
        mi = mi_sum / n_eval
        s = math.sqrt(ru * ru + mi * mi)
        if s < best_s:
            best_s, best_tau = s, tau
    return best_s, best_tau


def auprc(curve: list[tuple[float, float, float]]) -> float:
    """Area under the precision-recall curve via trapezoidal integration.

    ``curve`` is a list of ``(threshold, precision, recall)`` points. Points are
    sorted by increasing recall and integrated; precision is treated as a
    function of recall. Endpoints are not extrapolated.
    """
    pts = sorted(((r, p) for _t, p, r in curve))
    # collapse duplicate recall values to their max precision (upper envelope)
    dedup: dict[float, float] = {}
    for r, p in pts:
        if r not in dedup or p > dedup[r]:
            dedup[r] = p
    xs = sorted(dedup)
    if len(xs) < 2:
        return 0.0
    area = 0.0
    for i in range(1, len(xs)):
        r0, r1 = xs[i - 1], xs[i]
        p0, p1 = dedup[r0], dedup[r1]
        area += (r1 - r0) * (p0 + p1) / 2.0
    return area


def restrict_to_aspect(
    mapping: Mapping[str, Mapping[str, float] | set[str]],
    aspect: str,
    term_aspect: Mapping[str, str],
) -> dict:
    """Restrict a predictions or true-set mapping to terms of one aspect (BP/MF/CC)."""
    out: dict = {}
    for protein, terms in mapping.items():
        if isinstance(terms, dict):
            kept = {t: s for t, s in terms.items() if term_aspect.get(t) == aspect}
        else:
            kept = {t for t in terms if term_aspect.get(t) == aspect}
        if kept:
            out[protein] = kept
    return out


def filter_by_ic(
    mapping: Mapping[str, Mapping[str, float] | set[str]],
    ic: Mapping[str, float],
    min_ic: float,
) -> dict:
    """Drop terms whose information content is below ``min_ic`` (bits).

    Applied identically to truth and to every method's predictions, this removes
    near-universal, low-information terms (e.g. GO:0005515 "protein binding",
    which covers ~85% of human experimental MF annotations) so the evaluation
    rewards *informative* predictions rather than base-rate recovery. Proteins
    left with no terms are dropped. ``min_ic <= 0`` is a no-op (returns a copy).
    """
    out: dict = {}
    for protein, terms in mapping.items():
        if isinstance(terms, dict):
            kept = {t: s for t, s in terms.items() if ic.get(t, 0.0) >= min_ic}
        else:
            kept = {t for t in terms if ic.get(t, 0.0) >= min_ic}
        if kept:
            out[protein] = kept
    return out


def evaluate_aspect(
    pred_scores: PredScores,
    true_sets: TrueSets,
    ic: Mapping[str, float],
) -> dict[str, float]:
    """All three metrics for one already-aspect-restricted (pred, true) pair."""
    fmax, fmax_tau = f_max(pred_scores, true_sets, None)
    smin, smin_tau = s_min(pred_scores, true_sets, ic, None)
    curve = pr_curve(pred_scores, true_sets, None)
    return {
        "n_eval_proteins": len(true_sets),
        "f_max": fmax,
        "f_max_tau": fmax_tau,
        "s_min": smin,
        "s_min_tau": smin_tau,
        "auprc": auprc(curve),
    }


# --------------------------------------------------------------------------- #
# I/O + orchestration                                                          #
# --------------------------------------------------------------------------- #
def load_domain_go_scores(
    predictions_file, score_column: str = "p_value", neg_log10: bool = True
) -> dict[str, dict[str, float]]:
    """Load a dcGO associations TSV into ``{domain: {go_term: score}}``.

    With ``neg_log10`` (the default, for a p-value column) the score becomes
    ``-log10(p)`` so that higher = more confident and the confidence ordering has
    real dynamic range. This matters because ``hyper_score`` saturates — ~37% of
    significant associations sit at exactly 100 — which would collapse the
    threshold sweep. Set ``neg_log10=False`` to use a column (e.g. ``hyper_score``)
    verbatim.
    """
    import pandas as pd

    df = pd.read_csv(predictions_file, sep="\t")
    required = {"domain", "go_term", score_column}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Predictions file missing columns: {missing}")
    # Drop rows with NaN in a required column: a NaN score would poison -log10 and
    # a NaN key breaks lookups (NaN != NaN). Our files are clean; this is a guard.
    df = df.dropna(subset=list(required))
    out: dict[str, dict[str, float]] = defaultdict(dict)
    for domain, go, raw in zip(df["domain"], df["go_term"], df[score_column]):
        # -log10 with a floor so underflowed p == 0.0 doesn't become +inf.
        score = -math.log10(max(float(raw), 1e-320)) if neg_log10 else float(raw)
        prev = out[domain].get(go)
        if prev is None or score > prev:
            out[domain][go] = score
    return dict(out)


def _load_resampling():  # pragma: no cover - import plumbing
    """Import the sibling ``validation/resampling.py`` (validation/ is not a package)."""
    import importlib.util
    import sys
    from pathlib import Path

    if "resampling" in sys.modules:
        return sys.modules["resampling"]
    path = Path(__file__).resolve().parent / "resampling.py"
    spec = importlib.util.spec_from_file_location("resampling", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["resampling"] = module
    spec.loader.exec_module(module)
    return module


def build_term_aspect(processor) -> dict[str, str]:
    """Map each GO term to BP/MF/CC using the ontology namespace."""
    term_aspect: dict[str, str] = {}
    for term, data in processor.go_graph.nodes(data=True):
        aspect = NAMESPACE_TO_ASPECT.get(data.get("namespace", ""))
        if aspect:
            term_aspect[term] = aspect
    return term_aspect


def main() -> int:
    import argparse
    import sys
    from pathlib import Path

    import pandas as pd
    from loguru import logger

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    project_root = Path(__file__).resolve().parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    from src.domain_annotation_parser import DomainAnnotationParser
    from src.goa_parser import EXPERIMENTAL_EVIDENCE, GOAParser, parse_goa_human
    from src.ontology_processor import OntologyProcessor

    parser = argparse.ArgumentParser(
        description="Temporal held-out CAFA-style benchmark for dcGO (VALIDATION_PLAN §2)."
    )
    parser.add_argument(
        "--t0-gaf", type=Path, required=True, help="Early (training) GOA GAF"
    )
    parser.add_argument(
        "--t1-gaf", type=Path, required=True, help="Later (test) GOA GAF"
    )
    parser.add_argument(
        "--predictions",
        type=Path,
        required=True,
        help="dcGO associations TSV trained on t0",
    )
    parser.add_argument(
        "--interpro",
        type=Path,
        default=Path("data/interim/protein2ipr_human.dat.gz"),
        help="protein2ipr human subset (domain architectures)",
    )
    parser.add_argument(
        "--go-ontology",
        type=Path,
        default=Path("data/raw/go_ontology/go-basic.obo"),
    )
    parser.add_argument(
        "--score-column",
        default="p_value",
        help="Association score column (default: p_value, ranked as -log10)",
    )
    parser.add_argument(
        "--raw-score",
        dest="neg_log10",
        action="store_false",
        default=True,
        help="Use the score column verbatim instead of -log10 (e.g. for hyper_score)",
    )
    parser.add_argument(
        "--enable-supra-domains",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--disable-supra-domains",
        dest="enable_supra_domains",
        action="store_false",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("validation"))
    parser.add_argument(
        "--seed", type=int, default=0, help="Base seed for the permutation null"
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=100,
        help="Seeded domain-label permutations forming the random-domain null "
        "(default: 100). Permutation 0 uses --seed, so the 'random_domain' row in "
        "the metrics table is the first draw of the reported distribution. The "
        "smallest attainable empirical p is 1/(n+1).",
    )
    parser.add_argument(
        "--transfer",
        choices=["max", "pscore"],
        default="pscore",
        help="Domain->protein transfer: 'pscore' (default; Fang & Gough: sum of "
        "scores, min-max normalised per protein — the calibrated choice) or 'max' "
        "(max score, propagated).",
    )
    parser.add_argument(
        "--min-ic",
        type=float,
        action="append",
        metavar="BITS",
        help="Information-content floor(s): exclude terms below this IC (bits) from "
        "scoring, for truth and all methods alike. Repeatable to sweep. Default: "
        "0 (no filter) plus 2 and 4.",
    )
    args = parser.parse_args()
    ic_floors = sorted(set(args.min_ic)) if args.min_ic else [0.0, 2.0, 4.0]

    for p in (
        args.t0_gaf,
        args.t1_gaf,
        args.predictions,
        args.interpro,
        args.go_ontology,
    ):
        if not p.exists():
            logger.error(f"Missing required input: {p}")
            return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading GO ontology (for propagation + aspects)...")
    processor = OntologyProcessor(args.go_ontology)
    get_ancestors = processor.get_ancestors
    term_aspect = build_term_aspect(processor)

    logger.info("Parsing t0 (training) GOA — all non-IEA evidence...")
    # Same evidence filter the pipeline trains on. Used for IC, the naive term
    # frequencies, AND the no-knowledge gate — gating on anything narrower than
    # training would leak already-seen labels into the held-out set.
    t0_map = parse_goa_human(args.t0_gaf, evidence_filter="manual")

    logger.info("Parsing t1 (test) GOA — experimental evidence only...")
    t1_parser = GOAParser(evidence_codes=EXPERIMENTAL_EVIDENCE, aspects={"P", "F", "C"})
    t1_exp_map = t1_parser.parse_gaf_file(args.t1_gaf)

    logger.info("Parsing domain architectures...")
    dom_parser = DomainAnnotationParser(max_supra_domain_length=3, min_domain_length=10)
    architectures = dom_parser.parse_protein2ipr_file(args.interpro)
    protein_domains: dict[str, list[str]] = {}
    for protein, arch in architectures.items():
        domains = list(arch.single_domains)
        if args.enable_supra_domains:
            domains.extend(arch.supra_domains)
        if domains:
            protein_domains[protein] = domains

    logger.info("Building no-knowledge benchmark (per aspect, full t1-exp truth)...")
    benchmark = build_nk_benchmark_by_aspect(
        t0_map,
        t1_exp_map,
        term_aspect,
        get_ancestors,
        predictable_proteins=set(protein_domains),
    )
    eval_proteins = {p for aspect in benchmark.values() for p in aspect}
    for aspect, truths in benchmark.items():
        logger.info(f"  {aspect}: {len(truths):,} no-knowledge benchmark proteins")
    if not eval_proteins:
        logger.error(
            "Empty benchmark — the two snapshots may be too close together, or no "
            "no-knowledge protein has domains. Pick a wider t0/t1 gap."
        )
        return 1

    logger.info(f"Loading t0 domain->GO associations ({args.score_column})...")
    domain_go_scores = load_domain_go_scores(
        args.predictions, args.score_column, neg_log10=args.neg_log10
    )

    logger.info("Computing information content from t0 (propagated frequencies)...")
    ic = information_content(t0_map, get_ancestors)

    # Term frequency (propagated) among t0-annotated proteins for the naive baseline.
    n_t0 = len(t0_map)
    freq_counts: dict[str, int] = defaultdict(int)
    for terms in t0_map.values():
        for t in propagate_terms(terms, get_ancestors):
            freq_counts[t] += 1
    term_freq = {t: c / n_t0 for t, c in freq_counts.items()} if n_t0 else {}

    logger.info(
        f"Transferring dcGO predictions to benchmark proteins ({args.transfer})..."
    )
    transfer = (
        transfer_predictions_pscore
        if args.transfer == "pscore"
        else transfer_predictions
    )
    eval_domains = {p: protein_domains[p] for p in eval_proteins}
    dcgo_pred = transfer(eval_domains, domain_go_scores, get_ancestors)
    naive_pred = naive_predictions(eval_proteins, term_freq)
    shuffled = shuffle_domain_go(domain_go_scores, seed=args.seed)
    shuffle_pred = transfer(eval_domains, shuffled, get_ancestors)

    methods = {"dcGO": dcgo_pred, "naive": naive_pred, "random_domain": shuffle_pred}

    logger.info(f"Scoring at IC floors (bits): {ic_floors}")
    rows = []
    for aspect in ("BP", "MF", "CC"):
        if not benchmark[aspect]:
            logger.warning(f"No benchmark proteins for aspect {aspect}; skipping.")
            continue
        # Pre-restrict predictions to this aspect once, then filter per IC floor.
        aspect_preds = {
            method: restrict_to_aspect(preds, aspect, term_aspect)
            for method, preds in methods.items()
        }
        for min_ic in ic_floors:
            true_a = filter_by_ic(benchmark[aspect], ic, min_ic)
            if not true_a:
                continue
            for method, pred_a in aspect_preds.items():
                pred_f = filter_by_ic(pred_a, ic, min_ic)
                metrics = evaluate_aspect(pred_f, true_a, ic)
                metrics.update({"aspect": aspect, "method": method, "min_ic": min_ic})
                rows.append(metrics)
                logger.info(
                    f"  [{aspect} IC≥{min_ic:g}] {method:14s} "
                    f"F_max={metrics['f_max']:.3f} S_min={metrics['s_min']:.3f} "
                    f"AUPRC={metrics['auprc']:.3f} (n={metrics['n_eval_proteins']})"
                )

    result_df = pd.DataFrame(rows)[
        [
            "aspect",
            "min_ic",
            "method",
            "n_eval_proteins",
            "f_max",
            "f_max_tau",
            "s_min",
            "s_min_tau",
            "auprc",
        ]
    ]
    out_file = args.output_dir / "temporal_benchmark_metrics.tsv"
    result_df.to_csv(out_file, sep="\t", index=False)
    logger.info(f"✓ Saved metrics: {out_file}")

    # ---------------------------------------------------------------------- #
    # Permutation null. One shuffle is an anecdote; this repeats it under N
    # seeds and reports the distribution, a percentile interval and an
    # empirical p-value for the observed dcGO statistic against it.
    # ---------------------------------------------------------------------- #
    if args.n_permutations > 1:
        rs = _load_resampling()
        seeds = permutation_null_seeds(args.seed, args.n_permutations)
        logger.info(
            f"Permutation null: {len(seeds)} seeded domain-label permutations..."
        )
        null_samples: dict[tuple[str, float, str], list[float]] = defaultdict(list)
        for n_done, s in enumerate(seeds, start=1):
            perm_pred = transfer(
                eval_domains, shuffle_domain_go(domain_go_scores, seed=s), get_ancestors
            )
            for aspect in ("BP", "MF", "CC"):
                if not benchmark[aspect]:
                    continue
                pred_a = restrict_to_aspect(perm_pred, aspect, term_aspect)
                for min_ic in ic_floors:
                    true_a = filter_by_ic(benchmark[aspect], ic, min_ic)
                    if not true_a:
                        continue
                    pred_f = filter_by_ic(pred_a, ic, min_ic)
                    panel = rs.build_panel(
                        pred_f, true_a, ic, _candidate_thresholds(pred_f)
                    )
                    vals = rs.panel_metrics(panel)
                    null_samples[(aspect, min_ic, "f_max")].append(vals["f_max"])
                    null_samples[(aspect, min_ic, "auprc")].append(vals["auprc"])
            if n_done % 10 == 0:
                logger.info(f"  permutation {n_done}/{len(seeds)}")

        observed = {(r["aspect"], r["min_ic"], r["method"]): r for r in rows}
        null_rows = []
        for (aspect, min_ic, metric), samples in sorted(null_samples.items()):
            obs = observed[(aspect, min_ic, "dcGO")][metric]
            null_rows.append(
                {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "metric": metric,
                    "n_eval_proteins": observed[(aspect, min_ic, "dcGO")][
                        "n_eval_proteins"
                    ],
                    **rs.summarise_null(obs, samples),
                }
            )
        null_file = args.output_dir / "temporal_benchmark_permutation_null.tsv"
        pd.DataFrame(null_rows).to_csv(null_file, sep="\t", index=False)
        logger.info(f"✓ Saved permutation null: {null_file}")

    # Headline: does dcGO clear the naive F_max floor, and how does that change
    # as low-information terms are filtered out?
    logger.info("=" * 70)
    logger.info("TEMPORAL BENCHMARK COMPLETE — dcGO vs naive F_max by IC floor:")
    for aspect in ("BP", "MF", "CC"):
        for min_ic in ic_floors:
            sub = result_df[
                (result_df["aspect"] == aspect) & (result_df["min_ic"] == min_ic)
            ].set_index("method")
            if "dcGO" in sub.index and "naive" in sub.index:
                d = sub.loc["dcGO", "f_max"]
                nf = sub.loc["naive", "f_max"]
                verdict = "✓ above" if d > nf else "✗ below"
                logger.info(
                    f"  {aspect} IC≥{min_ic:g}: dcGO {d:.3f} vs naive {nf:.3f} "
                    f"[{verdict} floor]  (n={int(sub.loc['dcGO', 'n_eval_proteins'])})"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
