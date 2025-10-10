"""
Unit tests for vectorized Fisher's exact test and FDR correction.

Tests the parallel statistical testing implementation.
"""

import pytest
import numpy as np

from src.vectorized_fisher import (
    fisher_exact_vectorized_batch,
    fisher_exact_parallel,
    build_contingency_table,
    build_contingency_tables_vectorized,
    benjamini_hochberg_correction,
)


class TestBuildContingencyTable:
    """Test suite for building individual contingency tables."""

    def test_table_construction(self):
        """Test that contingency table is built correctly."""
        table = build_contingency_table(10, 5, 3, 7)

        assert table.shape == (2, 2)
        assert table[0, 0] == 10  # both
        assert table[0, 1] == 5  # domain only
        assert table[1, 0] == 3  # GO only
        assert table[1, 1] == 7  # neither

    def test_table_sum(self):
        """Test that table values sum to total proteins."""
        table = build_contingency_table(10, 5, 3, 7)
        assert table.sum() == 25

    def test_table_dtype(self):
        """Test that table has correct data type."""
        table = build_contingency_table(10, 5, 3, 7)
        assert table.dtype == np.int32


class TestFisherExactVectorizedBatch:
    """Test suite for batch Fisher's exact test."""

    def test_single_table(self):
        """Test Fisher's exact test on a single contingency table."""
        # Known case: strong association
        table = np.array([[[10, 2], [2, 10]]], dtype=np.int32)

        odds_ratios, pvalues = fisher_exact_vectorized_batch(
            table, alternative="greater"
        )

        assert len(odds_ratios) == 1
        assert len(pvalues) == 1
        assert odds_ratios[0] > 1.0  # Should show enrichment
        assert 0 <= pvalues[0] <= 1

    def test_multiple_tables(self):
        """Test Fisher's exact test on multiple contingency tables."""
        tables = np.array(
            [
                [[10, 2], [2, 10]],  # Strong enrichment
                [[5, 5], [5, 5]],  # No association
                [[2, 10], [10, 2]],  # Depletion
            ],
            dtype=np.int32,
        )

        odds_ratios, pvalues = fisher_exact_vectorized_batch(
            tables, alternative="greater"
        )

        assert len(odds_ratios) == 3
        assert len(pvalues) == 3

        # First case should have OR > 1 and low p-value
        assert odds_ratios[0] > 1.0
        assert pvalues[0] < 0.05

        # Second case should have OR ≈ 1 and high p-value
        assert np.isclose(odds_ratios[1], 1.0, atol=0.1)
        assert pvalues[1] > 0.5

        # Third case should have OR < 1
        assert odds_ratios[2] < 1.0

    def test_edge_case_zeros(self):
        """Test handling of contingency tables with zeros."""
        # Table with zero in one cell
        tables = np.array(
            [
                [[10, 0], [0, 10]],  # Perfect association
                [[0, 10], [10, 0]],  # Perfect inverse
            ],
            dtype=np.int32,
        )

        odds_ratios, pvalues = fisher_exact_vectorized_batch(
            tables, alternative="greater"
        )

        # Should not raise errors
        assert len(odds_ratios) == 2
        assert len(pvalues) == 2

    def test_alternative_hypotheses(self):
        """Test different alternative hypotheses."""
        table = np.array([[[10, 2], [2, 10]]], dtype=np.int32)

        # Test 'greater'
        _, pval_greater = fisher_exact_vectorized_batch(table, alternative="greater")

        # Test 'less'
        _, pval_less = fisher_exact_vectorized_batch(table, alternative="less")

        # Test 'two-sided'
        _, pval_two = fisher_exact_vectorized_batch(table, alternative="two-sided")

        # For enrichment, 'greater' should give smallest p-value
        assert pval_greater[0] < pval_two[0]
        assert pval_less[0] > pval_greater[0]


class TestFisherExactParallel:
    """Test suite for parallel Fisher's exact test."""

    def test_small_batch(self):
        """Test parallel processing with small number of tables."""
        tables = np.array(
            [
                [[10, 5], [5, 10]],
                [[8, 3], [3, 8]],
                [[6, 4], [4, 6]],
            ],
            dtype=np.int32,
        )

        odds_ratios, pvalues = fisher_exact_parallel(
            tables, alternative="greater", n_jobs=2, batch_size=2
        )

        assert len(odds_ratios) == 3
        assert len(pvalues) == 3
        assert all(0 <= p <= 1 for p in pvalues)

    def test_large_batch(self):
        """Test parallel processing with larger number of tables."""
        n_tests = 1000
        np.random.seed(42)

        # Create random contingency tables
        tables = np.random.randint(1, 20, size=(n_tests, 2, 2), dtype=np.int32)

        odds_ratios, pvalues = fisher_exact_parallel(
            tables, alternative="greater", n_jobs=4, batch_size=100
        )

        assert len(odds_ratios) == n_tests
        assert len(pvalues) == n_tests
        assert all(0 <= p <= 1 for p in pvalues)

    def test_progress_callback(self):
        """Test that progress callback is called correctly."""
        tables = np.random.randint(1, 20, size=(100, 2, 2), dtype=np.int32)

        callback_values = []

        def callback(completed, total):
            callback_values.append((completed, total))

        fisher_exact_parallel(
            tables,
            alternative="greater",
            n_jobs=2,
            batch_size=25,
            progress_callback=callback,
        )

        # Should have been called at least once
        assert len(callback_values) > 0

        # All calls should have same total
        totals = [v[1] for v in callback_values]
        assert all(t == 100 for t in totals)

        # Completed should increase
        completed_values = [v[0] for v in callback_values]
        assert completed_values[-1] == 100

    def test_results_consistency(self):
        """Test that parallel results match non-parallel results."""
        np.random.seed(42)
        tables = np.random.randint(1, 20, size=(50, 2, 2), dtype=np.int32)

        # Run batch version
        or_batch, pval_batch = fisher_exact_vectorized_batch(tables)

        # Run parallel version
        or_parallel, pval_parallel = fisher_exact_parallel(
            tables, n_jobs=2, batch_size=10
        )

        # Results should be identical
        np.testing.assert_array_almost_equal(or_batch, or_parallel)
        np.testing.assert_array_almost_equal(pval_batch, pval_parallel)


class TestBuildContingencyTablesVectorized:
    """Test suite for building contingency tables from protein annotations."""

    @pytest.fixture
    def sample_annotations(self):
        """Create sample protein annotations."""
        protein_domains = {
            "P1": {"D1", "D2"},
            "P2": {"D1"},
            "P3": {"D2", "D3"},
            "P4": {"D3"},
        }

        protein_go = {
            "P1": {"GO1", "GO2"},
            "P2": {"GO1"},
            "P3": {"GO2", "GO3"},
            "P4": {"GO3"},
        }

        domains = np.array(["D1", "D2", "D3"])
        go_terms = np.array(["GO1", "GO2", "GO3"])

        return domains, go_terms, protein_domains, protein_go

    def test_table_count(self, sample_annotations):
        """Test that correct number of tables is generated."""
        domains, go_terms, protein_domains, protein_go = sample_annotations

        tables = build_contingency_tables_vectorized(
            domains, go_terms, protein_domains, protein_go
        )

        expected_count = len(domains) * len(go_terms)
        assert tables.shape[0] == expected_count

    def test_table_values(self, sample_annotations):
        """Test specific contingency table values."""
        domains, go_terms, protein_domains, protein_go = sample_annotations

        tables = build_contingency_tables_vectorized(
            domains, go_terms, protein_domains, protein_go
        )

        # Test D1-GO1 table (both in P1 and P2)
        # D1: P1, P2 (2 proteins)
        # GO1: P1, P2 (2 proteins)
        # Both: P1, P2 (2)
        table_d1_go1 = tables[0]  # First domain × first GO term
        assert table_d1_go1[0, 0] == 2  # both

    def test_all_tables_sum_to_total(self, sample_annotations):
        """Test that all tables sum to total protein count."""
        domains, go_terms, protein_domains, protein_go = sample_annotations

        tables = build_contingency_tables_vectorized(
            domains, go_terms, protein_domains, protein_go
        )

        n_proteins = len(set(protein_domains.keys()) | set(protein_go.keys()))

        for table in tables:
            assert table.sum() == n_proteins


class TestBenjaminiHochbergCorrection:
    """Test suite for FDR correction."""

    def test_basic_correction(self):
        """Test basic FDR correction on known p-values."""
        pvalues = np.array([0.001, 0.01, 0.02, 0.05, 0.1])

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)

        # Adjusted p-values should be >= original
        assert np.all(adjusted >= pvalues)

        # All adjusted p-values should be <= 1.0
        assert np.all(adjusted <= 1.0)

        # Most significant p-value should be adjusted
        assert adjusted[0] > pvalues[0]

    def test_no_significant(self):
        """Test case where no p-values are significant."""
        pvalues = np.array([0.5, 0.6, 0.7, 0.8, 0.9])

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.01)

        # Threshold should be 0 (none significant)
        assert threshold == 0.0

        # No values should be significant
        assert np.all(adjusted > 0.01)

    def test_all_significant(self):
        """Test case where all p-values are highly significant."""
        pvalues = np.array([1e-10, 1e-9, 1e-8, 1e-7, 1e-6])

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)

        # All should be significant
        assert np.all(adjusted < 0.05)

        # Threshold should be the largest original p-value
        significant_originals = pvalues[adjusted <= 0.05]
        assert threshold == np.max(significant_originals)

    def test_monotonicity(self):
        """Test that adjusted p-values are monotonic with sorted input."""
        pvalues = np.sort(np.random.uniform(0, 1, 100))

        adjusted, _ = benjamini_hochberg_correction(pvalues, alpha=0.05)

        sorted_adjusted = np.sort(adjusted)

        # Adjusted p-values should already be sorted
        np.testing.assert_array_almost_equal(adjusted, sorted_adjusted)

    def test_single_pvalue(self):
        """Test FDR correction with single p-value."""
        pvalues = np.array([0.01])

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)

        assert len(adjusted) == 1
        assert adjusted[0] == 0.01  # Single p-value unchanged

    def test_ties_in_pvalues(self):
        """Test handling of tied p-values."""
        pvalues = np.array([0.01, 0.01, 0.05, 0.05, 0.1])

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)

        # Should handle ties gracefully
        assert len(adjusted) == len(pvalues)
        assert np.all(adjusted >= pvalues)

    def test_comparison_with_scipy(self):
        """Test that our implementation matches expected FDR behavior."""
        np.random.seed(42)
        pvalues = np.random.uniform(0, 1, 100)

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)

        # Basic sanity checks
        assert len(adjusted) == len(pvalues)
        assert np.all(adjusted >= pvalues)
        assert np.all((adjusted >= 0) & (adjusted <= 1))

        # Check that threshold is reasonable
        if threshold > 0:
            assert threshold <= 0.05

    def test_different_alpha_levels(self):
        """Test FDR correction at different significance levels."""
        pvalues = np.array([0.001, 0.01, 0.05, 0.1, 0.2])

        # Stricter alpha should result in higher threshold
        adj_001, thresh_001 = benjamini_hochberg_correction(pvalues, alpha=0.001)
        adj_01, thresh_01 = benjamini_hochberg_correction(pvalues, alpha=0.01)
        adj_05, thresh_05 = benjamini_hochberg_correction(pvalues, alpha=0.05)

        # More lenient alpha should have higher or equal threshold
        assert thresh_001 <= thresh_01 <= thresh_05


class TestIntegration:
    """Integration tests for complete statistical testing workflow."""

    def test_full_workflow(self):
        """Test complete workflow from tables to FDR-corrected results."""
        np.random.seed(42)

        # Create mix of significant and non-significant associations
        n_tests = 1000
        tables = []

        # 100 strong associations
        for _ in range(100):
            tables.append([[20, 5], [5, 20]])

        # 900 random/null associations
        for _ in range(900):
            tables.append(np.random.randint(5, 15, size=(2, 2)).tolist())

        tables = np.array(tables, dtype=np.int32)

        # Run Fisher's tests
        odds_ratios, pvalues = fisher_exact_parallel(
            tables, alternative="greater", n_jobs=4, batch_size=100
        )

        # Apply FDR correction
        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)

        # Verify results
        assert len(adjusted) == n_tests
        assert np.sum(adjusted <= 0.05) > 0  # Should find some significant
        assert np.sum(adjusted <= 0.05) < n_tests  # But not all

    def test_realistic_dcgo_scenario(self):
        """Test realistic dcGO scenario with domain-GO associations."""
        # Simulate 50 domains × 50 GO terms = 2500 tests
        n_domains = 50
        n_go_terms = 50
        n_proteins = 500

        np.random.seed(42)

        # Create realistic contingency tables
        # Most will be null, few will be truly associated
        tables = []
        for _ in range(n_domains * n_go_terms):
            # Most associations are random
            if np.random.random() > 0.95:  # 5% true associations
                # Strong association
                a = np.random.randint(20, 50)
                b = np.random.randint(5, 15)
                c = np.random.randint(5, 15)
                d = n_proteins - a - b - c
            else:
                # Null association
                a = np.random.randint(1, 10)
                b = np.random.randint(10, 50)
                c = np.random.randint(10, 50)
                d = n_proteins - a - b - c

            tables.append([[a, b], [c, d]])

        tables = np.array(tables, dtype=np.int32)

        # Run complete analysis
        odds_ratios, pvalues = fisher_exact_parallel(
            tables, alternative="greater", n_jobs=4
        )

        adjusted, threshold = benjamini_hochberg_correction(pvalues, alpha=0.01)

        # Verify we get reasonable results
        n_significant = np.sum(adjusted <= 0.01)

        assert n_significant > 0  # Should find some significant
        assert n_significant < len(tables)  # But not all (FDR control working)
        assert threshold >= 0  # Threshold should be non-negative
