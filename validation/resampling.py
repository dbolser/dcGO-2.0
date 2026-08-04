#!/usr/bin/env python3
"""Uncertainty for the CAFA-style temporal benchmark — bootstrap CIs + permutation nulls.

``validation/temporal_benchmark.py`` reports point estimates (``F_max``, ``S_min``,
``AUPRC``) and a *single* shuffled-domain null. The external review
(`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md`, P0 publication blockers) asked for two
things this module provides:

1. **Protein-level bootstrap confidence intervals** for every metric, and — the
   part that actually decides whether a difference is real — for the **paired
   difference** between two methods. Two independent CIs computed over the same
   candidate pool are *not* a test of a difference: this repository already
   learned that the expensive way in `SURPRISE_SCORE.md`, where independent
   intervals said the surprise score beat the q-value ranking and the paired
   bootstrap said there was no difference at all. So :func:`paired_bootstrap`
   resamples the benchmark proteins **once per replicate** and recomputes *every*
   method on that same resample.

2. **Many seeded permutations** instead of one shuffle, summarised as a null
   distribution with a percentile interval and an empirical p-value
   (:func:`percentile_ci`, :func:`empirical_p_value`).

The engineering trick that makes 1,000 bootstrap replicates affordable is the
:class:`EvaluationPanel`: for a fixed threshold grid, every CAFA metric is a
ratio of **per-protein sums**. Precompute each protein's contribution at each
threshold once, and a bootstrap replicate collapses to summing rows of a small
matrix. ``panel_metrics(panel, arange(n))`` reproduces
``temporal_benchmark.evaluate_aspect`` exactly (asserted in the unit tests).

The one approximation this makes, stated plainly: the threshold grid is computed
once on the **full** cohort and then held fixed across replicates. Recomputing
the grid per replicate would make the panel useless and would also make F_max
maximised over a resample-dependent grid, which is its own (worse) bias. Fixing
the grid is the standard choice.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


# --------------------------------------------------------------------------- #
# Panel construction                                                           #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class EvaluationPanel:
    """Per-protein, per-threshold metric contributions for one method + cohort.

    Every CAFA protein-centric metric is a ratio of sums over proteins, so the
    whole evaluation for an arbitrary (possibly resampled, possibly duplicated)
    set of proteins is a row-sum of these matrices.

    Attributes
    ----------
    proteins:
        Cohort in row order.
    thresholds:
        Ascending score cutoffs, shape ``(T,)``.
    tp:
        ``(N, T)`` true positives — predicted terms at or above the cutoff that
        are in the protein's truth set.
    n_pred:
        ``(N, T)`` number of terms predicted at or above the cutoff.
    n_true:
        ``(N,)`` size of each protein's truth set.
    ru:
        ``(N, T)`` remaining uncertainty — summed IC of true terms *not*
        predicted at the cutoff.
    mi:
        ``(N, T)`` misinformation — summed IC of predicted terms that are wrong.
    """

    proteins: tuple[str, ...]
    thresholds: np.ndarray
    tp: np.ndarray
    n_pred: np.ndarray
    n_true: np.ndarray
    ru: np.ndarray
    mi: np.ndarray

    @property
    def n_proteins(self) -> int:
        return len(self.proteins)


def build_panel(
    pred_scores: Mapping[str, Mapping[str, float]],
    true_sets: Mapping[str, set[str]],
    ic: Mapping[str, float],
    thresholds: Sequence[float],
) -> EvaluationPanel:
    """Precompute the per-protein/per-threshold contributions for one method.

    ``true_sets`` fixes the cohort and its row order (sorted for determinism).
    Proteins with an empty truth set are dropped — they are not scoreable and the
    reference implementation skips them too.
    """
    taus = np.asarray(sorted(thresholds), dtype=np.float64)
    n_tau = len(taus)
    proteins = tuple(sorted(p for p, t in true_sets.items() if t))
    n = len(proteins)

    tp = np.zeros((n, n_tau), dtype=np.float64)
    n_pred = np.zeros((n, n_tau), dtype=np.float64)
    n_true = np.zeros(n, dtype=np.float64)
    ru = np.zeros((n, n_tau), dtype=np.float64)
    mi = np.zeros((n, n_tau), dtype=np.float64)

    for i, protein in enumerate(proteins):
        truth = true_sets[protein]
        n_true[i] = len(truth)
        total_true_ic = sum(ic.get(t, 0.0) for t in truth)
        ru[i, :] = total_true_ic  # nothing predicted => all of it is missed

        preds = pred_scores.get(protein)
        if not preds:
            continue

        terms = list(preds)
        scores = np.fromiter(
            (preds[t] for t in terms), dtype=np.float64, count=len(terms)
        )
        hit = np.fromiter((t in truth for t in terms), dtype=bool, count=len(terms))
        term_ic = np.fromiter(
            (ic.get(t, 0.0) for t in terms), dtype=np.float64, count=len(terms)
        )

        # Sort by score DESCENDING: the first k entries are exactly the terms
        # predicted at a cutoff that admits k terms.
        order = np.argsort(-scores, kind="stable")
        s_desc = scores[order]
        hit_desc = hit[order]
        ic_desc = term_ic[order]

        cum_tp = np.concatenate(([0.0], np.cumsum(hit_desc.astype(np.float64))))
        cum_hit_ic = np.concatenate(([0.0], np.cumsum(ic_desc * hit_desc)))
        cum_miss_ic = np.concatenate(([0.0], np.cumsum(ic_desc * ~hit_desc)))

        # k(tau) = #{scores >= tau}. searchsorted on the ascending scores.
        s_asc = s_desc[::-1]
        k = len(s_asc) - np.searchsorted(s_asc, taus, side="left")

        n_pred[i, :] = k
        tp[i, :] = cum_tp[k]
        ru[i, :] = total_true_ic - cum_hit_ic[k]
        mi[i, :] = cum_miss_ic[k]

    return EvaluationPanel(
        proteins=proteins,
        thresholds=taus,
        tp=tp,
        n_pred=n_pred,
        n_true=n_true,
        ru=ru,
        mi=mi,
    )


# --------------------------------------------------------------------------- #
# Metrics from a panel (optionally over a resampled protein index)              #
# --------------------------------------------------------------------------- #
def _auprc_from_curve(precision: np.ndarray, recall: np.ndarray) -> float:
    """Trapezoidal AUPRC over the upper envelope of (recall, precision) points.

    Kept bit-identical in behaviour to ``temporal_benchmark.auprc`` so panel and
    reference implementations agree.
    """
    dedup: dict[float, float] = {}
    for r, p in zip(recall.tolist(), precision.tolist()):
        if r not in dedup or p > dedup[r]:
            dedup[r] = p
    xs = sorted(dedup)
    if len(xs) < 2:
        return 0.0
    area = 0.0
    for i in range(1, len(xs)):
        r0, r1 = xs[i - 1], xs[i]
        area += (r1 - r0) * (dedup[r0] + dedup[r1]) / 2.0
    return area


def panel_metrics(
    panel: EvaluationPanel, index: np.ndarray | None = None
) -> dict[str, float]:
    """CAFA metrics for the (possibly resampled) rows in ``index``.

    ``index`` is a row-index array into the panel, drawn **with replacement** for
    a bootstrap replicate; ``None`` means the full cohort in order, which
    reproduces the reference evaluator exactly.

    Returns ``f_max``, ``f_max_tau``, ``s_min``, ``s_min_tau``, ``auprc``,
    ``n_eval_proteins``, plus two coverage figures the review asked to be
    reported next to F_max: ``coverage_at_fmax`` (fraction of the cohort with at
    least one prediction at the F_max cutoff — CAFA averages precision over only
    those proteins while recall is averaged over all of them, so a high F_max on
    a thin slice of the cohort must be visible) and ``coverage_any`` (fraction
    with any prediction at all).
    """
    if panel.n_proteins == 0:
        return {
            "n_eval_proteins": 0,
            "f_max": 0.0,
            "f_max_tau": 0.0,
            "s_min": 0.0,
            "s_min_tau": 0.0,
            "auprc": 0.0,
            "coverage_at_fmax": 0.0,
            "coverage_any": 0.0,
        }
    idx = np.arange(panel.n_proteins) if index is None else np.asarray(index)
    n = len(idx)

    n_pred = panel.n_pred[idx]  # (n, T)
    has_pred = n_pred > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        prec_contrib = np.where(has_pred, panel.tp[idx] / np.maximum(n_pred, 1), 0.0)
    rec_contrib = panel.tp[idx] / panel.n_true[idx][:, None]

    m = has_pred.sum(axis=0).astype(np.float64)  # (T,)
    precision = np.where(m > 0, prec_contrib.sum(axis=0) / np.maximum(m, 1), 0.0)
    recall = rec_contrib.sum(axis=0) / n

    denom = precision + recall
    f = np.where(denom > 0, 2 * precision * recall / np.maximum(denom, 1e-300), 0.0)
    best = int(np.argmax(f))
    f_max = float(f[best])
    # Reference implementation starts at best_f = 0.0 / best_tau = 0.0 and only
    # improves on a strictly greater F, so an all-zero curve reports tau 0.0.
    f_max_tau = float(panel.thresholds[best]) if f_max > 0 else 0.0

    ru = panel.ru[idx].sum(axis=0) / n
    mi = panel.mi[idx].sum(axis=0) / n
    s = np.sqrt(ru * ru + mi * mi)
    s_best = int(np.argmin(s))

    return {
        "n_eval_proteins": n,
        "f_max": f_max,
        "f_max_tau": f_max_tau,
        "s_min": float(s[s_best]),
        "s_min_tau": float(panel.thresholds[s_best]),
        "auprc": _auprc_from_curve(precision, recall),
        "coverage_at_fmax": float(m[best] / n),
        "coverage_any": float(has_pred.any(axis=1).sum() / n),
    }


# --------------------------------------------------------------------------- #
# Bootstrap                                                                    #
# --------------------------------------------------------------------------- #
def percentile_ci(samples: Iterable[float], level: float = 0.95) -> tuple[float, float]:
    """Percentile interval of a bootstrap/null sample at the given coverage."""
    arr = np.asarray(list(samples), dtype=np.float64)
    if arr.size == 0:
        return (float("nan"), float("nan"))
    alpha = (1.0 - level) / 2.0
    lo, hi = np.quantile(arr, [alpha, 1.0 - alpha])
    return (float(lo), float(hi))


def paired_bootstrap(
    panels: Mapping[str, EvaluationPanel],
    metrics: Sequence[str] = ("f_max", "auprc"),
    n_replicates: int = 1000,
    seed: int = 0,
    level: float = 0.95,
) -> dict[str, np.ndarray]:
    """Bootstrap every method on the **same** protein resample, replicate by replicate.

    All panels must share a cohort (same proteins, same row order) — that is what
    makes the comparison paired, and it is checked rather than assumed. Each
    replicate draws ``n`` row indices with replacement and evaluates *every*
    method on those rows, so ``value[a] - value[b]`` is a paired difference and
    its percentile interval is a legitimate test of "does A beat B on these
    proteins", which two independent intervals are not.

    Returns ``{f"{method}::{metric}": array(n_replicates)}``. Differences are the
    caller's business (see :func:`summarise_paired`), because which pairs matter
    depends on the comparison.
    """
    names = list(panels)
    if not names:
        return {}
    reference = panels[names[0]].proteins
    for name in names[1:]:
        if panels[name].proteins != reference:
            raise ValueError(
                f"Panel {name!r} has a different cohort from {names[0]!r}; a paired "
                "bootstrap requires identical proteins in identical row order."
            )
    n = len(reference)
    rng = np.random.default_rng(seed)
    out = {f"{m}::{k}": np.empty(n_replicates) for m in names for k in metrics}
    for r in range(n_replicates):
        idx = rng.integers(0, n, size=n)
        for name in names:
            vals = panel_metrics(panels[name], idx)
            for k in metrics:
                out[f"{name}::{k}"][r] = vals[k]
    return out


def summarise_paired(
    replicates: Mapping[str, np.ndarray],
    method_a: str,
    method_b: str,
    metric: str,
    observed_a: float,
    observed_b: float,
    level: float = 0.95,
) -> dict[str, float]:
    """Percentile CI and two-sided bootstrap p-value for ``a - b`` on one metric.

    ``p_value`` is the standard bootstrap-percentile two-sided p: twice the
    smaller tail mass of the paired-difference distribution on the far side of
    zero, clipped at 1. It answers "how often does the resampled difference
    change sign", which is what "is this difference real" means here.
    """
    diff = replicates[f"{method_a}::{metric}"] - replicates[f"{method_b}::{metric}"]
    lo, hi = percentile_ci(diff, level)
    n = len(diff)
    if n == 0:
        p = float("nan")
    else:
        frac_le = float(np.mean(diff <= 0.0))
        frac_ge = float(np.mean(diff >= 0.0))
        p = min(1.0, 2.0 * min(frac_le, frac_ge))
    return {
        "metric": metric,
        "method_a": method_a,
        "method_b": method_b,
        "observed_a": observed_a,
        "observed_b": observed_b,
        "observed_diff": observed_a - observed_b,
        "diff_mean": float(np.mean(diff)) if n else float("nan"),
        "diff_ci_lo": lo,
        "diff_ci_hi": hi,
        "p_value": p,
        "significant": bool(n and (lo > 0.0 or hi < 0.0)),
        "n_replicates": n,
    }


# --------------------------------------------------------------------------- #
# Permutation nulls                                                            #
# --------------------------------------------------------------------------- #
def empirical_p_value(
    observed: float, null_samples: Iterable[float], alternative: str = "greater"
) -> float:
    """Empirical p-value of ``observed`` against a permutation null.

    Uses the ``(r + 1) / (n + 1)`` estimator (Phipson & Smyth 2010): a
    permutation p-value can never legitimately be 0, and the naive ``r / n`` is
    anti-conservative. With ``n`` permutations the smallest attainable p is
    ``1 / (n + 1)`` — so 100 permutations can never show more than ``p ≈ 0.0099``,
    which is a limit to report rather than to hide.
    """
    arr = np.asarray(list(null_samples), dtype=np.float64)
    n = arr.size
    if n == 0:
        return float("nan")
    if alternative == "greater":
        r = int(np.sum(arr >= observed))
    elif alternative == "less":
        r = int(np.sum(arr <= observed))
    elif alternative == "two-sided":
        centre = float(np.mean(arr))
        r = int(np.sum(np.abs(arr - centre) >= abs(observed - centre)))
    else:
        raise ValueError(
            f"alternative must be greater/less/two-sided, got {alternative!r}"
        )
    return (r + 1) / (n + 1)


def summarise_null(
    observed: float,
    null_samples: Iterable[float],
    level: float = 0.95,
    alternative: str = "greater",
) -> dict[str, float]:
    """Mean, spread, percentile interval and empirical p-value for a null sample."""
    arr = np.asarray(list(null_samples), dtype=np.float64)
    lo, hi = percentile_ci(arr, level)
    mean = float(np.mean(arr)) if arr.size else float("nan")
    sd = float(np.std(arr, ddof=1)) if arr.size > 1 else float("nan")
    z = (observed - mean) / sd if arr.size > 1 and sd > 0 else float("nan")
    return {
        "observed": observed,
        "null_n": int(arr.size),
        "null_mean": mean,
        "null_sd": sd,
        "null_min": float(np.min(arr)) if arr.size else float("nan"),
        "null_max": float(np.max(arr)) if arr.size else float("nan"),
        "null_ci_lo": lo,
        "null_ci_hi": hi,
        "ratio_observed_over_null_mean": (
            observed / mean
            if arr.size and mean not in (0.0,) and not math.isnan(mean)
            else float("nan")
        ),
        "z_score": float(z),
        "empirical_p": empirical_p_value(observed, arr, alternative),
    }
