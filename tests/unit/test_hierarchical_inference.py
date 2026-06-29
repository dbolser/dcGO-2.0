"""Unit tests for hierarchical inference shrinkage statistics."""

import numpy as np

from src.hierarchical_inference import HierarchicalInferenceEngine
from src.sparse_fisher import build_domain_metadata


def _engine():
    return HierarchicalInferenceEngine(shrinkage_strength=0.5, min_observations=3)


class TestShrinkageStatistics:
    """Test suite for get_shrinkage_statistics."""

    def test_no_supra_domains_returns_finite_zeroed_stats(self):
        """With only single domains there are no supra tests to summarize.

        Regression test: np.median of an empty slice previously returned nan,
        making baseline / no-supra runs emit non-finite metrics.
        """
        domain_list = ["IPR1", "IPR2"]
        protein_domains = {"P1": {"IPR1"}, "P2": {"IPR2"}}
        metadata = build_domain_metadata(domain_list, protein_domains)

        # 2 domains x 2 GO terms
        pvalues = np.array([0.2, 0.8, 0.05, 0.5], dtype=np.float64)
        stats = _engine().get_shrinkage_statistics(
            pvalues, pvalues, domain_list, metadata
        )

        assert stats["n_supra_tests"] == 0
        assert stats["median_pvalue_ratio"] == 1.0
        assert np.isfinite(stats["median_pvalue_ratio"])
        assert stats["n_pvalues_increased"] == 0
        assert stats["pct_pvalues_increased"] == 0.0
        assert stats["lost_significance"] == 0
        assert stats["gained_significance"] == 0

    def test_with_supra_domains_reports_finite_stats(self):
        """The normal path still produces finite stats and counts supra tests."""
        domain_list = ["IPR1", "IPR2", "IPR1,IPR2"]
        protein_domains = {
            "P1": {"IPR1", "IPR2", "IPR1,IPR2"},
            "P2": {"IPR1"},
            "P3": {"IPR2"},
        }
        go_list = ["GO1", "GO2"]
        metadata = build_domain_metadata(domain_list, protein_domains)

        # 3 domains x 2 GO terms
        original = np.array([0.2, 0.8, 0.05, 0.5, 0.01, 0.4], dtype=np.float64)
        engine = _engine()
        shrunk = engine.shrink_pvalues(original, domain_list, go_list, metadata)
        stats = engine.get_shrinkage_statistics(original, shrunk, domain_list, metadata)

        # One supra domain ("IPR1,IPR2") across 2 GO terms
        assert stats["n_supra_tests"] == len(go_list)
        assert np.isfinite(stats["median_pvalue_ratio"])
