"""
Hierarchical Statistical Inference for dcGO Pipeline

This module implements empirical Bayes shrinkage to regularize supra-domain
p-values by "borrowing strength" from their constituent single domains.

The key insight: A supra-domain's association with a GO term should be informed
by how well its constituent domains associate with that same GO term.
"""

from typing import Dict, List

import numpy as np
from loguru import logger

from src.sparse_fisher import DomainMetadata


class HierarchicalInferenceEngine:
    """
    Applies hierarchical shrinkage to supra-domain statistical tests.

    This engine regularizes low-count supra-domain p-values by shrinking them
    toward the geometric mean of their constituent domains' p-values.
    """

    def __init__(self, shrinkage_strength: float = 0.5, min_observations: int = 3):
        """
        Initialize the hierarchical inference engine.

        Args:
            shrinkage_strength: Base shrinkage factor (0-1). Higher = more shrinkage.
                                0 = no shrinkage, 1 = full shrinkage to prior
            min_observations: Minimum observations before applying full shrinkage
        """
        if not 0 <= shrinkage_strength <= 1:
            raise ValueError(
                f"Shrinkage strength must be in [0, 1], got {shrinkage_strength}"
            )
        if min_observations <= 0:
            raise ValueError(f"min_observations must be > 0, got {min_observations}")

        self.shrinkage_strength = shrinkage_strength
        self.min_observations = min_observations

        logger.info("HierarchicalInferenceEngine initialized:")
        logger.info(f"  Shrinkage strength: {shrinkage_strength}")
        logger.info(f"  Minimum observations: {min_observations}")

    def shrink_pvalues(
        self,
        pvalues: np.ndarray,
        domain_list: List[str],
        go_list: List[str],
        domain_metadata: Dict[str, DomainMetadata],
    ) -> np.ndarray:
        """
        Apply hierarchical shrinkage to p-values.

        Args:
            pvalues: Array of p-values (shape: n_domains * n_go_terms)
            domain_list: Ordered list of domain IDs
            go_list: Ordered list of GO term IDs
            domain_metadata: Metadata for each domain

        Returns:
            Shrunk p-values (same shape as input)
        """
        logger.info("Applying hierarchical shrinkage to supra-domain p-values...")

        n_domains = len(domain_list)
        n_go = len(go_list)
        shrunk_pvalues = pvalues.copy()

        # Build mapping from domain ID to its p-values for each GO term
        domain_to_pvalues = {}
        for d_idx, domain_id in enumerate(domain_list):
            domain_to_pvalues[domain_id] = pvalues[d_idx * n_go : (d_idx + 1) * n_go]

        # Process only supra-domains
        supra_domains_processed = 0
        total_shrinkage = 0.0

        for d_idx, domain_id in enumerate(domain_list):
            meta = domain_metadata[domain_id]

            # Skip single domains (no shrinkage needed)
            if meta.is_single_domain:
                continue

            supra_domains_processed += 1

            # Get constituent domain p-values
            constituent_pvalues = []
            for constituent_id in meta.constituent_domains:
                if constituent_id in domain_to_pvalues:
                    constituent_pvalues.append(domain_to_pvalues[constituent_id])

            if not constituent_pvalues:
                # No constituent data available - skip shrinkage
                continue

            # Stack constituent p-values (shape: n_constituents, n_go)
            constituent_array = np.array(constituent_pvalues)

            # Compute geometric mean as prior
            # Use geometric mean because p-values are multiplicative (not additive)
            # Add small epsilon to avoid log(0)
            epsilon = 1e-300
            log_pvalues = np.log(constituent_array + epsilon)
            log_mean = np.mean(log_pvalues, axis=0)
            prior_pvalues = np.exp(log_mean)

            # Adaptive shrinkage based on observation count
            # More observations → less shrinkage (trust the data)
            # Fewer observations → more shrinkage (trust the prior)
            alpha = self._compute_adaptive_shrinkage(meta.observation_count)

            # Apply shrinkage: p_shrunk = (1-α)*p_observed + α*p_prior
            # In log space for numerical stability
            start_idx = d_idx * n_go
            end_idx = (d_idx + 1) * n_go

            observed_pvalues = pvalues[start_idx:end_idx]
            log_observed = np.log(observed_pvalues + epsilon)
            log_prior = np.log(prior_pvalues + epsilon)
            log_shrunk = (1 - alpha) * log_observed + alpha * log_prior
            shrunk_pvalues[start_idx:end_idx] = np.exp(log_shrunk)

            # Track shrinkage amount
            total_shrinkage += alpha

        avg_shrinkage = (
            total_shrinkage / supra_domains_processed
            if supra_domains_processed > 0
            else 0
        )

        logger.info(f"  Processed {supra_domains_processed:,} supra-domains")
        logger.info(f"  Average shrinkage factor: {avg_shrinkage:.3f}")

        return shrunk_pvalues

    def _compute_adaptive_shrinkage(self, n_observations: int) -> float:
        """
        Compute adaptive shrinkage factor based on observation count.

        Low observation count → high shrinkage (rely on prior)
        High observation count → low shrinkage (rely on data)

        Args:
            n_observations: Number of proteins with this domain

        Returns:
            Shrinkage factor alpha in [0, 1]
        """
        # Exponential decay: α = strength * exp(-n / threshold)
        # When n = threshold, α = strength * e^(-1) ≈ 0.37 * strength
        # When n → ∞, α → 0 (no shrinkage)
        # When n = 0, α = strength (maximum shrinkage)

        alpha = self.shrinkage_strength * np.exp(
            -n_observations / self.min_observations
        )

        # Ensure bounds
        return max(0.0, min(1.0, alpha))

    def get_shrinkage_statistics(
        self,
        original_pvalues: np.ndarray,
        shrunk_pvalues: np.ndarray,
        domain_list: List[str],
        domain_metadata: Dict[str, DomainMetadata],
        significance_threshold: float = 0.01,
    ) -> Dict:
        """
        Compute statistics about shrinkage effects.

        Args:
            original_pvalues: P-values before shrinkage
            shrunk_pvalues: P-values after shrinkage
            domain_list: List of domain IDs
            domain_metadata: Domain metadata dictionary
            significance_threshold: Threshold for significance

        Returns:
            Dictionary of shrinkage statistics
        """
        n_go = len(shrunk_pvalues) // len(domain_list)

        # Track changes for supra-domains only
        supra_original = []
        supra_shrunk = []

        for d_idx, domain_id in enumerate(domain_list):
            meta = domain_metadata[domain_id]
            if meta.is_supra_domain:
                start = d_idx * n_go
                end = (d_idx + 1) * n_go
                supra_original.extend(original_pvalues[start:end])
                supra_shrunk.extend(shrunk_pvalues[start:end])

        supra_original = np.array(supra_original)
        supra_shrunk = np.array(supra_shrunk)

        # Compute statistics
        pvalue_increases = supra_shrunk > supra_original
        n_increased = np.sum(pvalue_increases)
        pct_increased = (
            100 * n_increased / len(supra_original) if len(supra_original) > 0 else 0
        )

        # Median change
        pvalue_ratio = supra_shrunk / (supra_original + 1e-300)
        median_ratio = np.median(pvalue_ratio)

        # Significance changes
        orig_sig = supra_original < significance_threshold
        shrunk_sig = supra_shrunk < significance_threshold
        lost_significance = np.sum(orig_sig & ~shrunk_sig)
        gained_significance = np.sum(~orig_sig & shrunk_sig)

        return {
            "n_supra_tests": len(supra_original),
            "n_pvalues_increased": int(n_increased),
            "pct_pvalues_increased": float(pct_increased),
            "median_pvalue_ratio": float(median_ratio),
            "lost_significance": int(lost_significance),
            "gained_significance": int(gained_significance),
        }
