"""Knowledge yield: how many correct new facts a predictor delivers at a fixed error.

The pre-registered criterion of ``KNOWLEDGE_PROJECTION.md`` §1.2. Every
(protein, term) prediction in an aspect's cohort is pooled into one ranked
list — a curator works down one list, not one protein at a time — and a
threshold is swept over the distinct scores. At threshold ``τ``:

- ``TP(τ)`` = predictions scoring ``≥ τ`` that are in their protein's truth;
- ``PP(τ)`` = all predictions scoring ``≥ τ``.

``KY@π`` is the largest ``TP(τ)`` over thresholds whose pooled precision
``TP(τ) / PP(τ)`` is at least ``π``: the number of true facts the method can
hand over while at most ``1 − π`` of what it hands over is wrong. A method
whose precision never reaches ``π`` has ``KY@π = 0``.

Only proteins in ``truth`` (the cohort) count; predictions for any other
protein are ignored, and cohort proteins without predictions simply add
nothing. Truth and predictions are expected on the same term vocabulary and
IC floor already (``ablation.floor_cell`` style); this module does no
propagation.

Uncertainty is a paired bootstrap over cohort proteins: each replicate draws
protein multiplicities once and applies them to every method, so the ratio of
two methods' yields is computed on the same resample.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

#: The pre-registered precision levels; 0.5 is primary.
PRECISION_LEVELS: tuple[float, ...] = (0.3, 0.5, 0.7)


@dataclass(frozen=True)
class PooledRanking:
    """One method's cohort predictions, pooled and sorted by score (descending).

    ``group_end[i]`` is the exclusive end index of the ``i``-th block of equal
    scores in the sorted arrays, so cumulative sums read at ``group_end - 1``
    give the counts at each distinct threshold.
    """

    scores: np.ndarray  # float64, sorted descending
    is_tp: np.ndarray  # bool
    ic: np.ndarray  # float64, IC of the term (0 where unknown)
    protein: np.ndarray  # int64 index into ``proteins``
    group_end: np.ndarray  # int64
    proteins: tuple[str, ...]
    n_truth: int  # total (protein, term) truth pairs in the cohort


def pool(
    truth: Mapping[str, set[str]],
    preds: Mapping[str, Mapping[str, float]],
    ic: Mapping[str, float] | None = None,
    proteins: Sequence[str] | None = None,
) -> PooledRanking:
    """Pool a method's cohort predictions into one ranking.

    ``proteins`` fixes the protein index order (pass the same sequence for
    every method so bootstrap multiplicities line up); it defaults to the
    sorted cohort.
    """
    order = tuple(sorted(truth)) if proteins is None else tuple(proteins)
    index = {p: i for i, p in enumerate(order)}
    ic = ic or {}
    scores: list[float] = []
    is_tp: list[bool] = []
    ics: list[float] = []
    prot: list[int] = []
    for p in order:
        if p not in truth:
            continue
        true_terms = truth[p]
        for term, s in preds.get(p, {}).items():
            scores.append(float(s))
            is_tp.append(term in true_terms)
            ics.append(float(ic.get(term, 0.0)))
            prot.append(index[p])
    s_arr = np.asarray(scores, dtype=np.float64)
    # Stable descending sort; ties are resolved by grouping, so their order
    # never matters.
    sort = np.argsort(-s_arr, kind="stable")
    s_sorted = s_arr[sort]
    if s_sorted.size:
        change = np.flatnonzero(np.diff(s_sorted) != 0) + 1
        group_end = np.append(change, s_sorted.size).astype(np.int64)
    else:
        group_end = np.zeros(0, dtype=np.int64)
    return PooledRanking(
        scores=s_sorted,
        is_tp=np.asarray(is_tp, dtype=bool)[sort],
        ic=np.asarray(ics, dtype=np.float64)[sort],
        protein=np.asarray(prot, dtype=np.int64)[sort],
        group_end=group_end,
        proteins=order,
        n_truth=sum(len(truth[p]) for p in order if p in truth),
    )


@dataclass(frozen=True)
class YieldAtPrecision:
    """``KY@π`` and what sits at its threshold."""

    level: float
    yield_tp: float  # KY@π (weighted counts under a bootstrap)
    n_predicted: float  # PP at the threshold
    threshold: float  # score cutoff (nan when the level is never reached)
    bits: float  # sum of IC over the true positives at the threshold
    n_proteins: int  # distinct proteins with >= 1 TP at the threshold


def _curves(
    r: PooledRanking, weights: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cumulative (TP, PP, bits) at each distinct threshold."""
    w = np.ones(r.scores.size) if weights is None else weights[r.protein]
    tp = np.cumsum(w * r.is_tp)[r.group_end - 1]
    pp = np.cumsum(w)[r.group_end - 1]
    bits = np.cumsum(w * r.is_tp * r.ic)[r.group_end - 1]
    return tp, pp, bits


def yield_at_precision(
    r: PooledRanking, level: float, weights: np.ndarray | None = None
) -> YieldAtPrecision:
    """``KY@level`` of one pooled ranking (optionally bootstrap-weighted)."""
    if r.scores.size == 0:
        return YieldAtPrecision(level, 0.0, 0.0, float("nan"), 0.0, 0)
    tp, pp, bits = _curves(r, weights)
    with np.errstate(invalid="ignore", divide="ignore"):
        prec = np.where(pp > 0, tp / pp, 0.0)
    ok = np.flatnonzero(prec >= level - 1e-12)
    if ok.size == 0:
        return YieldAtPrecision(level, 0.0, 0.0, float("nan"), 0.0, 0)
    best = ok[np.argmax(tp[ok])]
    # Earliest threshold reaching that yield: most precise list with it.
    best = ok[np.flatnonzero(tp[ok] == tp[best])[0]]
    end = int(r.group_end[best])
    hit = r.is_tp[:end]
    w = np.ones(end) if weights is None else weights[r.protein[:end]]
    n_prot = int(np.unique(r.protein[:end][hit & (w > 0)]).size)
    return YieldAtPrecision(
        level=level,
        yield_tp=float(tp[best]),
        n_predicted=float(pp[best]),
        threshold=float(r.scores[end - 1]),
        bits=float(bits[best]),
        n_proteins=n_prot,
    )


def pooled_auprc(r: PooledRanking) -> float:
    """Pooled (micro) average precision: ``Σ ΔR · P`` over distinct thresholds.

    No interpolation and no origin anchor (the CAFA-evaluator convention the
    project adopted on 2026-10-15). Recall's denominator is every truth pair
    in the cohort, including those of proteins the method never predicts for.
    """
    return _weighted_auprc(r, None, float(r.n_truth))


def _weighted_auprc(
    r: PooledRanking, weights: np.ndarray | None, n_truth: float
) -> float:
    if r.scores.size == 0 or n_truth == 0:
        return 0.0
    tp, pp, _ = _curves(r, weights)
    recall = tp / n_truth
    with np.errstate(invalid="ignore", divide="ignore"):
        prec = np.where(pp > 0, tp / pp, 0.0)
    d_recall = np.diff(np.concatenate(([0.0], recall)))
    return float(np.sum(d_recall * prec))


def summarise(r: PooledRanking, levels: Sequence[float] = PRECISION_LEVELS) -> dict:
    """Point estimates: KY, bits and proteins at each level, and pooled AUPRC."""
    out: dict[str, float] = {
        "n_truth": float(r.n_truth),
        "n_pred": float(r.scores.size),
    }
    for level in levels:
        y = yield_at_precision(r, level)
        tag = f"{int(round(level * 100))}"
        out[f"ky{tag}"] = y.yield_tp
        out[f"bits{tag}"] = y.bits
        out[f"proteins{tag}"] = float(y.n_proteins)
        out[f"tau{tag}"] = y.threshold
    out["pooled_auprc"] = pooled_auprc(r)
    return out


def bootstrap_weights(n: int, n_replicates: int, seed: int) -> np.ndarray:
    """Protein multiplicities, one row per replicate (multinomial resample)."""
    rng = np.random.default_rng(seed)
    p = np.full(n, 1.0 / n)
    return np.stack([rng.multinomial(n, p) for _ in range(n_replicates)]).astype(
        np.float64
    )


def cluster_bootstrap_weights(
    clusters: np.ndarray, n_replicates: int, seed: int
) -> np.ndarray:
    """Protein multiplicities from resampling *clusters* with replacement.

    ``clusters[i]`` labels protein ``i``; each replicate draws as many
    clusters as there are, uniformly with replacement, and every protein
    takes its cluster's multiplicity. Proteins that share evidence (for
    example, the papers that annotated them) then move together.
    """
    labels, inverse = np.unique(clusters, return_inverse=True)
    rng = np.random.default_rng(seed)
    k = labels.size
    out = np.empty((n_replicates, clusters.size))
    for i in range(n_replicates):
        mult = np.bincount(rng.integers(0, k, size=k), minlength=k)
        out[i] = mult[inverse]
    return out


def precision_at(
    r: PooledRanking, tau: float, weights: np.ndarray | None = None
) -> np.ndarray | float:
    """Pooled precision of the predictions scoring at least ``tau``.

    With a ``(replicates, proteins)`` weight matrix, one value per replicate.
    """
    end = int(np.searchsorted(-r.scores, -tau, side="right"))
    if end == 0:
        return np.full(len(weights), np.nan) if weights is not None else float("nan")
    if weights is None:
        return float(r.is_tp[:end].mean())
    w = weights[:, r.protein[:end]]
    pp = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(pp > 0, (w * r.is_tp[:end]).sum(axis=1) / pp, np.nan)


def bootstrap(
    rankings: Mapping[str, PooledRanking],
    truth_sizes: np.ndarray,
    n_replicates: int = 1000,
    seed: int = 0,
    levels: Sequence[float] = PRECISION_LEVELS,
    chunk: int = 32,
    weights: np.ndarray | None = None,
) -> dict[str, dict[str, np.ndarray]]:
    """Paired bootstrap over cohort proteins.

    Every ranking must have been pooled with the same ``proteins`` order;
    ``truth_sizes[i]`` is the number of truth terms of protein ``i``. Returns
    ``{method: {"ky50": array, ..., "pooled_auprc": array}}`` of replicate
    values, paired across methods by replicate index (the same multiplicity
    rows, :func:`bootstrap_weights` unless ``weights`` is given, for every
    method). Replicates are
    evaluated ``chunk`` at a time as matrices; each replicate's value equals
    :func:`yield_at_precision` / pooled AUPRC under its weights.
    """
    names = list(rankings)
    n = truth_sizes.size
    for name in names:
        if len(rankings[name].proteins) != n:
            raise ValueError(f"{name}: pooled on a different protein order")
    if weights is None:
        weights = bootstrap_weights(n, n_replicates, seed)
    n_replicates = weights.shape[0]
    n_truth = weights @ truth_sizes
    keys = [f"ky{int(round(lv * 100))}" for lv in levels]
    out: dict[str, dict[str, np.ndarray]] = {}
    for m in names:
        r = rankings[m]
        res = {k: np.zeros(n_replicates) for k in keys} | {
            "pooled_auprc": np.zeros(n_replicates)
        }
        if r.scores.size:
            ends = r.group_end - 1
            for start in range(0, n_replicates, chunk):
                w = weights[start : start + chunk][:, r.protein]
                tp = np.cumsum(w * r.is_tp, axis=1)[:, ends]
                pp = np.cumsum(w, axis=1)[:, ends]
                with np.errstate(invalid="ignore", divide="ignore"):
                    prec = np.where(pp > 0, tp / pp, 0.0)
                for key, level in zip(keys, levels):
                    ok = prec >= level - 1e-12
                    res[key][start : start + chunk] = np.where(ok, tp, 0.0).max(axis=1)
                nt = n_truth[start : start + chunk, None]
                with np.errstate(invalid="ignore", divide="ignore"):
                    recall = np.where(nt > 0, tp / nt, 0.0)
                d_recall = np.diff(
                    np.concatenate([np.zeros((recall.shape[0], 1)), recall], axis=1),
                    axis=1,
                )
                res["pooled_auprc"][start : start + chunk] = np.sum(
                    d_recall * prec, axis=1
                )
        out[m] = res
    return out


def ratio_interval(
    num: np.ndarray, den: np.ndarray, level: float = 0.95
) -> tuple[float, float]:
    """Percentile interval of ``num / den`` over paired replicates.

    Replicates with a zero denominator and a positive numerator are counted
    as ``+inf`` (the numerator method wins outright there). Replicates where
    both are zero say nothing about the ratio and are dropped; with none
    left the interval is ``(nan, nan)``.
    """
    keep = (num > 0) | (den > 0)
    num, den = num[keep], den[keep]
    if num.size == 0:
        return float("nan"), float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.inf)
    lo, hi = np.quantile(
        ratio, [(1 - level) / 2, (1 + level) / 2], method="inverted_cdf"
    )
    return float(lo), float(hi)
