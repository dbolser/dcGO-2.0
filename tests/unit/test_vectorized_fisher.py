"""
Unit tests for vectorized Fisher's exact test and FDR correction.

Tests the parallel statistical testing implementation.
"""

import numpy as np
import pytest

from src.vectorized_fisher import (
    benjamini_hochberg_by_family,
    fisher_exact_vectorized_batch,
    fisher_exact_parallel,
    build_contingency_table,
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


class TestBenjaminiHochbergVectorizationEquivalence:
    """The vectorized BH must be bit-identical to the loop it replaced.

    That loop cost ~50 of the ~69 minutes of a default human run (1.64e9
    tests) — three times the compiled Fisher stage it corrects. Replacing it is
    only safe if the numbers do not move, so this pins the reference
    implementation rather than trusting that the rewrite "looks equivalent".
    """

    @staticmethod
    def _reference(pvalues: np.ndarray, alpha: float = 0.05):
        """The pre-vectorization implementation, verbatim."""
        n = len(pvalues)
        sorted_indices = np.argsort(pvalues)
        sorted_pvalues = pvalues[sorted_indices]
        adjusted = np.zeros(n, dtype=np.float64)
        for i in range(n - 1, -1, -1):
            rank = i + 1
            adjusted[sorted_indices[i]] = min(1.0, sorted_pvalues[i] * n / rank)
            if i < n - 1:
                adjusted[sorted_indices[i]] = min(
                    adjusted[sorted_indices[i]], adjusted[sorted_indices[i + 1]]
                )
        significant = adjusted <= alpha
        threshold = np.max(pvalues[significant]) if np.any(significant) else 0.0
        return adjusted, threshold

    @pytest.mark.parametrize("alpha", [0.01, 0.05])
    @pytest.mark.parametrize(
        "name",
        [
            "uniform",
            "underflow",
            "ties",
            "all_ones",
            "all_zeros",
            "single",
            "mixed",
        ],
    )
    def test_matches_the_reference_exactly(self, name: str, alpha: float) -> None:
        rng = np.random.default_rng(0)
        cases = {
            # Ordinary case.
            "uniform": rng.random(5000),
            # p-values below 1e-300 are common and meaningful here.
            "underflow": rng.random(2000) ** 40,
            # Ties are where a naive rewrite diverges: the running minimum and
            # the per-element loop must agree on every member of a tied block.
            "ties": np.repeat(rng.random(50), 40),
            "all_ones": np.ones(200),
            "all_zeros": np.zeros(200),
            "single": np.array([0.5]),
            "mixed": np.concatenate([np.zeros(50), np.ones(50), rng.random(200)]),
        }
        pvalues = cases[name]

        want_adjusted, want_threshold = self._reference(pvalues.copy(), alpha)
        got_adjusted, got_threshold = benjamini_hochberg_correction(
            pvalues.copy(), alpha
        )

        assert np.array_equal(want_adjusted, got_adjusted)
        assert want_threshold == got_threshold

    def test_empty_input_returns_empty(self) -> None:
        """The reference indexed into an empty sort; this must not raise."""
        adjusted, threshold = benjamini_hochberg_correction(np.array([]), alpha=0.05)
        assert adjusted.shape == (0,)
        assert threshold == 0.0

    def test_adjusted_values_are_monotone_in_the_p_value_order(self) -> None:
        """The property BH's step-up procedure exists to guarantee."""
        rng = np.random.default_rng(7)
        pvalues = rng.random(5000)
        adjusted, _ = benjamini_hochberg_correction(pvalues, alpha=0.05)
        in_p_order = adjusted[np.argsort(pvalues)]
        assert np.all(np.diff(in_p_order) >= 0)


class TestBenjaminiHochbergByFamily:
    """Single domains and supra-domains are corrected as separate families.

    A supra-domain is not an exchangeable sibling of its own constituents, and
    pooling them let the 5.3x larger supra space tighten the threshold for
    single-domain hypotheses that gain nothing from it.
    """

    def test_each_family_matches_correcting_it_alone(self) -> None:
        """The whole contract: a family's result must not depend on the other."""
        rng = np.random.default_rng(0)
        single_p = rng.random(400) ** 3
        supra_p = rng.random(1600) ** 2
        pvalues = np.concatenate([single_p, supra_p])
        family = np.array(["single"] * 400 + ["supra"] * 1600)

        adjusted, thresholds = benjamini_hochberg_by_family(pvalues, family, alpha=0.01)

        single_alone, single_threshold = benjamini_hochberg_correction(
            single_p, alpha=0.01
        )
        supra_alone, supra_threshold = benjamini_hochberg_correction(
            supra_p, alpha=0.01
        )

        assert np.array_equal(adjusted[:400], single_alone)
        assert np.array_equal(adjusted[400:], supra_alone)
        assert thresholds["single"] == single_threshold
        assert thresholds["supra"] == supra_threshold

    def test_a_large_second_family_no_longer_penalises_the_first(self) -> None:
        """The reason for the change, stated as a test.

        Pooling makes every single-domain q-value depend on how many
        supra-domain hypotheses happen to be in the run. Splitting removes that.
        """
        rng = np.random.default_rng(1)
        single_p = rng.random(200) ** 3

        small = np.concatenate([single_p, rng.random(200)])
        large = np.concatenate([single_p, rng.random(20000)])
        family_small = np.array(["single"] * 200 + ["supra"] * 200)
        family_large = np.array(["single"] * 200 + ["supra"] * 20000)

        adj_small, _ = benjamini_hochberg_by_family(small, family_small, alpha=0.01)
        adj_large, _ = benjamini_hochberg_by_family(large, family_large, alpha=0.01)
        assert np.array_equal(adj_small[:200], adj_large[:200])

        # Pooled, the same single-domain p-values would have been penalised.
        pooled_small, _ = benjamini_hochberg_correction(small, alpha=0.01)
        pooled_large, _ = benjamini_hochberg_correction(large, alpha=0.01)
        assert not np.array_equal(pooled_small[:200], pooled_large[:200])

    def test_single_family_reduces_to_plain_bh(self) -> None:
        rng = np.random.default_rng(2)
        pvalues = rng.random(500)
        family = np.array(["single"] * 500)
        adjusted, thresholds = benjamini_hochberg_by_family(pvalues, family, alpha=0.05)
        want, want_threshold = benjamini_hochberg_correction(pvalues, alpha=0.05)
        assert np.array_equal(adjusted, want)
        assert thresholds == {"single": want_threshold}

    def test_mismatched_lengths_raise(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            benjamini_hochberg_by_family(
                np.array([0.1, 0.2]), np.array(["single"]), alpha=0.05
            )
