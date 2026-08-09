"""
Vectorized Fisher's exact test implementation for dcGO pipeline.

This module provides efficient batch processing of Fisher's exact tests using
the Cython ``fisher`` package (``fisher.pvalue_npy``), which evaluates an entire
array of 2x2 tables in one compiled call instead of looping ``scipy`` per test.
"""

import numpy as np
from typing import Tuple
from fisher import pvalue_npy


def fisher_exact_vectorized_batch(
    contingency_tables: np.ndarray, alternative: str = "greater"
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Fisher's exact test for a batch of 2x2 contingency tables.

    Args:
        contingency_tables: Array of shape (n, 2, 2) containing n contingency tables
        alternative: 'greater', 'less', or 'two-sided'

    Returns:
        Tuple of (odds_ratios, pvalues) arrays of shape (n,)
    """
    # fisher.pvalue_npy requires C-contiguous uint32 input arrays.
    a = np.ascontiguousarray(contingency_tables[:, 0, 0], dtype=np.uint32)
    b = np.ascontiguousarray(contingency_tables[:, 0, 1], dtype=np.uint32)
    c = np.ascontiguousarray(contingency_tables[:, 1, 0], dtype=np.uint32)
    d = np.ascontiguousarray(contingency_tables[:, 1, 1], dtype=np.uint32)

    # Returns a (left_tail, right_tail, two_tail) tuple of p-value arrays.
    left_tail, right_tail, two_tail = pvalue_npy(a, b, c, d)
    if alternative == "greater":
        pvalues = right_tail
    elif alternative == "less":
        pvalues = left_tail
    elif alternative == "two-sided":
        pvalues = two_tail
    else:
        raise ValueError(
            f"alternative must be 'greater', 'less', or 'two-sided', got {alternative!r}"
        )

    # Sample odds ratio (a*d)/(b*c), matching scipy.stats.fisher_exact.
    # float64 avoids overflow; b*c == 0 yields inf (or nan for 0/0), as in scipy.
    numerator = a.astype(np.float64) * d.astype(np.float64)
    denominator = b.astype(np.float64) * c.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        odds_ratios = numerator / denominator

    return odds_ratios, pvalues


def fisher_exact_parallel(
    contingency_tables: np.ndarray,
    alternative: str = "greater",
    n_jobs: int = -1,
    batch_size: int = 10000,
    progress_callback=None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Fisher's exact test for a large array of contingency tables.

    The heavy lifting now runs in compiled Cython (``fisher.pvalue_npy``), which
    is fast enough that no multiprocessing is needed; tables are processed in
    in-process chunks purely so ``progress_callback`` can report incremental
    progress. ``n_jobs`` is accepted for backward compatibility but unused.

    Args:
        contingency_tables: Array of shape (n, 2, 2)
        alternative: 'greater', 'less', or 'two-sided'
        n_jobs: Unused; retained for API compatibility
        batch_size: Number of tests per progress chunk
        progress_callback: Optional callback function(completed, total) for progress updates

    Returns:
        Tuple of (odds_ratios, pvalues) arrays of shape (n,)
    """
    n_tests = contingency_tables.shape[0]
    odds_ratios = np.empty(n_tests, dtype=np.float64)
    pvalues = np.empty(n_tests, dtype=np.float64)

    for start in range(0, n_tests, batch_size):
        end = min(start + batch_size, n_tests)
        chunk_odds, chunk_pvalues = fisher_exact_vectorized_batch(
            contingency_tables[start:end], alternative
        )
        odds_ratios[start:end] = chunk_odds
        pvalues[start:end] = chunk_pvalues

        if progress_callback:
            progress_callback(end, n_tests)

    return odds_ratios, pvalues


def build_contingency_table(
    n_domain_and_go: int, n_domain_not_go: int, n_go_not_domain: int, n_neither: int
) -> np.ndarray:
    """
    Build a 2x2 contingency table from counts.

    Format:
                    Has GO term | Doesn't have GO term
    Has domain:         a       |         b
    No domain:          c       |         d

    Args:
        n_domain_and_go: Proteins with both domain and GO term (a)
        n_domain_not_go: Proteins with domain but not GO term (b)
        n_go_not_domain: Proteins with GO term but not domain (c)
        n_neither: Proteins with neither (d)

    Returns:
        2x2 numpy array
    """
    return np.array(
        [[n_domain_and_go, n_domain_not_go], [n_go_not_domain, n_neither]],
        dtype=np.int32,
    )


def benjamini_hochberg_correction(
    pvalues: np.ndarray, alpha: float = 0.05, n_hypotheses: "int | None" = None
) -> Tuple[np.ndarray, float]:
    """
    Apply Benjamini-Hochberg FDR correction to p-values.

    Vectorized. This used to be a Python loop over every p-value, which on the
    default human run (1.64e9 tests) cost ~50 of the run's ~69 minutes — three
    times the compiled Fisher stage it corrects. The arithmetic is unchanged:
    scale each sorted p-value by ``n / rank``, enforce monotonicity by taking a
    running minimum from the largest p-value down, and clip at 1.

    Clipping commutes with the running minimum (``min(1, ·)`` is monotone), so
    doing it once at the end is equivalent to the old per-element ``min(1.0, …)``.

    Memory matters more than speed here: at 1.64e9 tests the sort indices alone
    are 13 GB, so the intermediates are computed in place and released as soon
    as they are dead rather than chained into one expression.

    ``n_hypotheses`` exists for the sparse enumeration in
    ``compute_cooccurring_contingency_tables``: it tests only the domain-term
    pairs that co-occur, because every other pair has p=1 exactly under a
    one-sided ``greater`` test. BH must still divide by the size of the **whole**
    family, so the caller passes it here.

    Omitting the p=1 tail changes nothing else. Those entries occupy the last
    ranks, their scaled values ``1.0 * n / j`` for ``j > k`` bottom out at 1.0 at
    ``j = n``, and the reverse running minimum therefore contributes exactly 1.0
    — which cannot lower any adjusted value that is already at most 1. Passing
    ``n_hypotheses`` reproduces the dense result on the enumerated subset
    exactly, not approximately; leaving it at ``None`` on a subset would divide
    by the wrong denominator and inflate every rejection.

    Args:
        pvalues: Array of p-values
        alpha: FDR threshold (e.g., 0.01 for 1% FDR)
        n_hypotheses: Size of the hypothesis family, when ``pvalues`` is a
            subset whose omitted members are all exactly 1. Defaults to
            ``len(pvalues)``, the dense case.

    Returns:
        Tuple of (adjusted_pvalues, threshold) where threshold is the p-value cutoff

    Raises:
        ValueError: if ``n_hypotheses`` is smaller than the number of p-values,
            which would mean the caller mis-stated the family and silently
            deflated its own q-values.
    """
    k = len(pvalues)
    n = k if n_hypotheses is None else int(n_hypotheses)
    if n < k:
        raise ValueError(
            f"n_hypotheses ({n}) is smaller than the number of p-values "
            f"({k}); the family cannot be smaller than the tests in it"
        )
    if k == 0:
        return np.zeros(0, dtype=np.float64), 0.0

    order = np.argsort(pvalues)

    # scaled = sorted_p * n / rank, built in place to avoid a second big temporary.
    scaled = pvalues[order].astype(np.float64, copy=False)
    if scaled is pvalues:  # already float64 and unsorted-copy elided
        scaled = scaled.copy()
    # Scale by the FAMILY size n, but rank within the k tests actually present.
    # When the omitted tail is all p=1 these ranks are the true ones: p=1 sorts
    # last, so every enumerated test keeps the rank it would have had densely.
    scaled *= n
    scaled /= np.arange(1, k + 1, dtype=np.float64)

    # Monotonicity: each adjusted value is the smallest scaled value at or above
    # its rank. A reversed accumulate gives that in one pass.
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    np.clip(scaled, None, 1.0, out=scaled)

    adjusted = np.empty(k, dtype=np.float64)
    adjusted[order] = scaled
    del scaled, order

    # Find threshold: largest p-value where adjusted p-value <= alpha
    significant = adjusted <= alpha
    if np.any(significant):
        threshold = float(np.max(pvalues[significant]))
    else:
        threshold = 0.0

    return adjusted, threshold


def benjamini_hochberg_by_family(
    pvalues: np.ndarray,
    family: np.ndarray,
    alpha: float = 0.05,
    labels: "dict | None" = None,
    family_sizes: "dict | None" = None,
) -> Tuple[np.ndarray, dict]:
    """Apply BH separately within each hypothesis family.

    A supra-domain is not an exchangeable sibling of its own constituent
    domains: it is a different kind of hypothesis, tested on a nested subset of
    the same proteins. Correcting both in one family means the 5.3x larger
    supra-domain space makes the threshold stricter for single-domain
    hypotheses that nobody asked to be penalised by it, and it is part of what
    the review means by "BH is applied across a highly dependent hierarchical
    hypothesis family".

    Each family controls FDR at ``alpha`` within itself. That is a deliberate
    choice and not the same guarantee as one pooled correction: the expected
    proportion of false discoveries is controlled per family, so a reader who
    pools the two output sets is looking at a union of two separately
    controlled sets, not a set controlled at ``alpha`` overall.

    Pass ``family`` in the most compact dtype that distinguishes the families —
    a bool for two of them — and use ``labels`` to name them in the result. At
    the scale this runs (1.69e9 tests on a supra-enabled human run) a ``<U6``
    string array would be **40.6 GB**, against 1.69 GB for a bool, and would
    risk exhausting memory on the main analysis path before the correction even
    starts.

    Args:
        pvalues: Array of p-values.
        family: Array of the same length assigning each test to a family.
        alpha: FDR threshold applied within each family.
        labels: Optional map from a ``family`` value to the name it should carry
            in the returned ``thresholds``. Values absent from the map keep
            their raw form.
        family_sizes: Optional map from a ``family`` value to that family's full
            hypothesis count, for the sparse enumeration where ``pvalues`` holds
            only the co-occurring pairs and every omitted pair has p=1. Each
            family is then corrected against its own dense size — for dcGO that
            is ``(domains in the family) * n_terms``, which differs between the
            single and supra families and so cannot be inferred here. Absent
            entries fall back to the number of p-values in that family.

    Returns:
        ``(adjusted, thresholds)`` where ``thresholds`` maps each family label
        to its own p-value cutoff.
    """
    if len(pvalues) != len(family):
        raise ValueError(
            f"pvalues and family must be the same length, got "
            f"{len(pvalues)} and {len(family)}"
        )
    adjusted = np.empty(len(pvalues), dtype=np.float64)
    thresholds: dict = {}
    for label in np.unique(family):
        members = np.flatnonzero(family == label)
        key = label.item() if hasattr(label, "item") else label
        family_adjusted, family_threshold = benjamini_hochberg_correction(
            pvalues[members],
            alpha=alpha,
            n_hypotheses=(family_sizes or {}).get(key),
        )
        adjusted[members] = family_adjusted
        name = (
            labels.get(label.item() if hasattr(label, "item") else label, label)
            if labels
            else label
        )
        thresholds[name] = family_threshold
    return adjusted, thresholds
