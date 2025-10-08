"""
Vectorized Fisher's exact test implementation for dcGO pipeline.

This module provides efficient batch processing of Fisher's exact tests
using NumPy arrays and parallel processing.
"""

import numpy as np
from scipy import stats
from typing import Tuple
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp


def fisher_exact_vectorized_batch(
    contingency_tables: np.ndarray,
    alternative: str = 'greater'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Fisher's exact test for a batch of 2x2 contingency tables.

    Args:
        contingency_tables: Array of shape (n, 2, 2) containing n contingency tables
        alternative: 'greater', 'less', or 'two-sided'

    Returns:
        Tuple of (odds_ratios, pvalues) arrays of shape (n,)
    """
    n_tests = contingency_tables.shape[0]
    odds_ratios = np.zeros(n_tests, dtype=np.float64)
    pvalues = np.zeros(n_tests, dtype=np.float64)

    # Process each contingency table
    for i in range(n_tests):
        table = contingency_tables[i]
        try:
            odds_ratio, pvalue = stats.fisher_exact(table, alternative=alternative)
            odds_ratios[i] = odds_ratio
            pvalues[i] = pvalue
        except (ValueError, ZeroDivisionError):
            # Handle edge cases (e.g., all zeros)
            odds_ratios[i] = np.nan
            pvalues[i] = 1.0

    return odds_ratios, pvalues


def _fisher_batch_wrapper(args):
    """Wrapper for parallel processing (must be top-level function for pickling)."""
    batch, alternative = args
    return fisher_exact_vectorized_batch(batch, alternative)


def fisher_exact_parallel(
    contingency_tables: np.ndarray,
    alternative: str = 'greater',
    n_jobs: int = -1,
    batch_size: int = 10000,
    progress_callback=None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Fisher's exact test in parallel for large arrays of contingency tables.

    Args:
        contingency_tables: Array of shape (n, 2, 2)
        alternative: 'greater', 'less', or 'two-sided'
        n_jobs: Number of parallel jobs (-1 for all CPUs)
        batch_size: Number of tests per batch
        progress_callback: Optional callback function(completed, total) for progress updates

    Returns:
        Tuple of (odds_ratios, pvalues) arrays of shape (n,)
    """
    n_tests = contingency_tables.shape[0]

    if n_jobs == -1:
        n_jobs = mp.cpu_count()

    # Split into batches
    batches = []
    for i in range(0, n_tests, batch_size):
        end_idx = min(i + batch_size, n_tests)
        batches.append((contingency_tables[i:end_idx], alternative))

    # Process batches in parallel with progress tracking
    all_results = []
    completed = 0

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = [executor.submit(_fisher_batch_wrapper, batch) for batch in batches]

        for future in futures:
            result = future.result()
            all_results.append(result)
            completed += len(result[0])

            if progress_callback:
                progress_callback(completed, n_tests)

    # Concatenate results
    odds_ratios = np.concatenate([r[0] for r in all_results])
    pvalues = np.concatenate([r[1] for r in all_results])

    return odds_ratios, pvalues


def build_contingency_table(
    n_domain_and_go: int,
    n_domain_not_go: int,
    n_go_not_domain: int,
    n_neither: int
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
    return np.array([
        [n_domain_and_go, n_domain_not_go],
        [n_go_not_domain, n_neither]
    ], dtype=np.int32)


def build_contingency_tables_vectorized(
    domains: np.ndarray,
    go_terms: np.ndarray,
    protein_domains: dict,
    protein_go: dict
) -> np.ndarray:
    """
    Build contingency tables for all domain-GO combinations.

    Args:
        domains: Array of domain IDs
        go_terms: Array of GO term IDs
        protein_domains: Dict mapping protein_id -> set of domain IDs
        protein_go: Dict mapping protein_id -> set of GO term IDs

    Returns:
        Array of shape (n_domains * n_go_terms, 2, 2)
    """
    n_domains = len(domains)
    n_go_terms = len(go_terms)
    n_tests = n_domains * n_go_terms

    # Get all proteins
    all_proteins = set(protein_domains.keys()) | set(protein_go.keys())
    n_total_proteins = len(all_proteins)

    # Pre-allocate result array
    tables = np.zeros((n_tests, 2, 2), dtype=np.int32)

    # Build each table
    idx = 0
    for domain in domains:
        # Get proteins with this domain
        proteins_with_domain = {p for p, doms in protein_domains.items() if domain in doms}
        n_with_domain = len(proteins_with_domain)
        n_without_domain = n_total_proteins - n_with_domain

        for go_term in go_terms:
            # Get proteins with this GO term
            proteins_with_go = {p for p, gos in protein_go.items() if go_term in gos}
            n_with_go = len(proteins_with_go)

            # Count overlaps
            n_domain_and_go = len(proteins_with_domain & proteins_with_go)
            n_domain_not_go = n_with_domain - n_domain_and_go
            n_go_not_domain = n_with_go - n_domain_and_go
            n_neither = n_total_proteins - n_domain_and_go - n_domain_not_go - n_go_not_domain

            tables[idx] = build_contingency_table(
                n_domain_and_go, n_domain_not_go, n_go_not_domain, n_neither
            )
            idx += 1

    return tables


def benjamini_hochberg_correction(pvalues: np.ndarray, alpha: float = 0.05) -> Tuple[np.ndarray, float]:
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
                adjusted[sorted_indices[i]],
                adjusted[sorted_indices[i + 1]]
            )

    # Find threshold: largest p-value where adjusted p-value <= alpha
    significant = adjusted <= alpha
    if np.any(significant):
        threshold = np.max(pvalues[significant])
    else:
        threshold = 0.0

    return adjusted, threshold
