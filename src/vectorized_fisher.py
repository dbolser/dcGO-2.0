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
    pvalues: np.ndarray, alpha: float = 0.05
) -> Tuple[np.ndarray, float]:
    """
    Apply Benjamini-Hochberg FDR correction to p-values.

    Args:
        pvalues: Array of p-values
        alpha: FDR threshold (e.g., 0.01 for 1% FDR)

    Returns:
        Tuple of (adjusted_pvalues, threshold) where threshold is the p-value cutoff
    """
    n = len(pvalues)

    # Sort p-values and track original indices
    sorted_indices = np.argsort(pvalues)
    sorted_pvalues = pvalues[sorted_indices]

    # Calculate BH-adjusted p-values
    adjusted = np.zeros(n, dtype=np.float64)

    # Work backwards from largest p-value
    for i in range(n - 1, -1, -1):
        rank = i + 1
        adjusted[sorted_indices[i]] = min(1.0, sorted_pvalues[i] * n / rank)

        # Ensure monotonicity
        if i < n - 1:
            adjusted[sorted_indices[i]] = min(
                adjusted[sorted_indices[i]], adjusted[sorted_indices[i + 1]]
            )

    # Find threshold: largest p-value where adjusted p-value <= alpha
    significant = adjusted <= alpha
    if np.any(significant):
        threshold = np.max(pvalues[significant])
    else:
        threshold = 0.0

    return adjusted, threshold
