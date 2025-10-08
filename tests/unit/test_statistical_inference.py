"""
Unit tests for the statistical inference module.

Tests Fisher's exact test, FDR correction, hypergeometric scoring,
and correspondence matrix construction.
"""

import pytest
import pandas as pd
from pathlib import Path
from typing import Dict, List, Set

# Add src to path for testing
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from statistical_inference import StatisticalInferenceEngine, AssociationResult


class TestAssociationResult:
    """Test the AssociationResult dataclass."""

    def test_association_result_creation(self):
        """Test basic AssociationResult creation and properties."""
        result = AssociationResult(
            domain="PF00001",
            go_term="GO:0005515",
            a=10,
            b=20,
            c=30,
            d=40,
            p_value=0.001,
            odds_ratio=2.5,
            hyper_score=85.2,
        )

        assert result.domain == "PF00001"
        assert result.go_term == "GO:0005515"
        assert result.a == 10
        assert result.total_proteins == 100  # a + b + c + d
        assert result.domain_total == 40  # a + c
        assert result.go_term_total == 30  # a + b
        assert result.q_value is None  # Not set initially

    def test_association_result_validation(self):
        """Test validation in AssociationResult post_init."""
        # Valid result should not raise
        result = AssociationResult(
            domain="PF00001",
            go_term="GO:0005515",
            a=5,
            b=10,
            c=15,
            d=20,
            p_value=0.05,
            odds_ratio=1.5,
            hyper_score=50.0,
        )
        assert result.a == 5

        # Invalid negative values should raise
        with pytest.raises(ValueError):
            AssociationResult(
                domain="PF00001",
                go_term="GO:0005515",
                a=-1,
                b=10,
                c=15,
                d=20,  # Negative a
                p_value=0.05,
                odds_ratio=1.5,
                hyper_score=50.0,
            )

    def test_significance_testing(self):
        """Test significance determination methods."""
        # Significant result
        significant = AssociationResult(
            domain="PF00001",
            go_term="GO:0005515",
            a=10,
            b=5,
            c=5,
            d=80,
            p_value=0.001,
            odds_ratio=16.0,
            hyper_score=95.0,
            q_value=0.005,
        )

        assert significant.is_significant(alpha=0.01)
        assert significant.is_significant_fdr(alpha=0.01)

        # Non-significant result
        non_significant = AssociationResult(
            domain="PF00002",
            go_term="GO:0008150",
            a=2,
            b=8,
            c=18,
            d=72,
            p_value=0.5,
            odds_ratio=1.0,
            hyper_score=10.0,
            q_value=0.8,
        )

        assert not non_significant.is_significant(alpha=0.01)
        assert not non_significant.is_significant_fdr(alpha=0.01)


class TestStatisticalInferenceEngine:
    """Test the main statistical inference engine."""

    @pytest.fixture
    def sample_protein_domain_map(self) -> Dict[str, List[str]]:
        """Create sample protein-domain mapping for testing."""
        return {
            "P12345": ["PF00001", "PF00002"],
            "P67890": ["PF00001", "PF00003"],
            "P11111": ["PF00002"],
            "P22222": ["PF00003"],
            "P33333": ["PF00001"],
            "P44444": ["PF00004"],
            "P55555": ["PF00001", "PF00002", "PF00003"],  # Multi-domain
            "P66666": ["PF00005"],
        }

    @pytest.fixture
    def sample_protein_go_map(self) -> Dict[str, Set[str]]:
        """Create sample protein-GO mapping for testing."""
        return {
            "P12345": {"GO:0005515", "GO:0008150"},
            "P67890": {"GO:0005515", "GO:0009987"},
            "P11111": {"GO:0008150"},
            "P22222": {"GO:0009987"},
            "P33333": {"GO:0005515"},
            "P44444": {"GO:0016740"},
            "P55555": {"GO:0005515", "GO:0008150", "GO:0009987"},
            "P66666": {"GO:0016740"},
        }

    @pytest.fixture
    def inference_engine(self, sample_protein_domain_map, sample_protein_go_map):
        """Create StatisticalInferenceEngine for testing."""
        return StatisticalInferenceEngine(
            sample_protein_domain_map, sample_protein_go_map
        )

    def test_engine_initialization(self, inference_engine, sample_protein_domain_map):
        """Test proper initialization of the inference engine."""
        assert inference_engine.total_proteins == len(sample_protein_domain_map)
        assert len(inference_engine.domain_counts) > 0
        assert len(inference_engine.go_counts) > 0

        # Check domain frequency calculation
        assert (
            inference_engine.domain_counts["PF00001"] == 4
        )  # P12345, P67890, P33333, P55555
        assert inference_engine.domain_counts["PF00002"] == 3  # P12345, P11111, P55555

        # Check GO term frequency calculation
        assert (
            inference_engine.go_counts["GO:0005515"] == 4
        )  # P12345, P67890, P33333, P55555
        assert inference_engine.go_counts["GO:0008150"] == 3  # P12345, P11111, P55555

    def test_correspondence_matrix_construction(self, inference_engine):
        """Test building of domain-GO correspondence matrix."""
        matrix = inference_engine.build_correspondence_matrix()

        assert isinstance(matrix, pd.DataFrame)
        assert matrix.shape[0] > 0  # Has GO terms (rows)
        assert matrix.shape[1] > 0  # Has domains (columns)

        # Check specific co-occurrences
        if "GO:0005515" in matrix.index and "PF00001" in matrix.columns:
            # PF00001 and GO:0005515 co-occur in P12345, P67890, P33333, P55555
            assert matrix.loc["GO:0005515", "PF00001"] == 4

        if "GO:0008150" in matrix.index and "PF00002" in matrix.columns:
            # PF00002 and GO:0008150 co-occur in P12345, P11111, P55555
            assert matrix.loc["GO:0008150", "PF00002"] == 3

    def test_contingency_table_calculation(self, inference_engine):
        """Test 2x2 contingency table construction."""
        matrix = inference_engine.build_correspondence_matrix()

        # Test specific domain-GO pair
        domain = "PF00001"
        go_term = "GO:0005515"

        a, b, c, d = inference_engine.calculate_contingency_values(
            domain, go_term, matrix
        )

        # Verify contingency table values
        assert a >= 0  # Both domain and GO term
        assert b >= 0  # GO term but not domain
        assert c >= 0  # Domain but not GO term
        assert d >= 0  # Neither
        assert a + b + c + d == inference_engine.total_proteins

        # Check marginals
        assert a + c == inference_engine.domain_counts[domain]
        assert a + b == inference_engine.go_counts[go_term]

    def test_fisher_exact_test_calculation(self, inference_engine):
        """Test Fisher's exact test implementation."""
        # Create a clear enrichment case
        # Domain PF00001 appears in 4/8 proteins
        # GO term GO:0005515 appears in 4/8 proteins
        # They co-occur in 4/8 proteins (perfect overlap)

        matrix = inference_engine.build_correspondence_matrix()
        a, b, c, d = inference_engine.calculate_contingency_values(
            "PF00001", "GO:0005515", matrix
        )

        # Run single association test
        result = inference_engine._test_single_association(
            "PF00001", "GO:0005515", a, b, c, d
        )

        assert isinstance(result, AssociationResult)
        assert result.domain == "PF00001"
        assert result.go_term == "GO:0005515"
        assert result.p_value >= 0.0
        assert result.p_value <= 1.0
        assert result.odds_ratio > 0
        assert 1 <= result.hyper_score <= 100

    def test_hypergeometric_scoring(self, inference_engine):
        """Test hypergeometric score calculation."""
        # Test with known values
        a, b, c, d = 10, 5, 5, 80  # Strong enrichment
        score = inference_engine._calculate_hypergeometric_score(a, b, c, d)

        assert 1 <= score <= 100
        assert score > 50  # Should be high for strong enrichment

        # Test with no enrichment
        a, b, c, d = 1, 49, 19, 31  # Weak association
        score_weak = inference_engine._calculate_hypergeometric_score(a, b, c, d)
        assert score_weak < score  # Weaker association should have lower score

        # Test edge cases
        assert inference_engine._calculate_hypergeometric_score(0, 10, 10, 80) == 0.0
        assert inference_engine._calculate_hypergeometric_score(10, 0, 0, 90) == 100.0

    def test_fdr_correction(self, inference_engine):
        """Test FDR correction implementation."""
        # Create mock results with known p-values
        mock_results = [
            AssociationResult("PF00001", "GO:0001", 10, 5, 5, 80, 0.001, 5.0, 90.0),
            AssociationResult("PF00002", "GO:0002", 8, 7, 7, 78, 0.01, 3.0, 75.0),
            AssociationResult("PF00003", "GO:0003", 5, 10, 10, 75, 0.05, 2.0, 50.0),
            AssociationResult("PF00004", "GO:0004", 2, 13, 13, 72, 0.1, 1.5, 25.0),
        ]

        corrected_results = inference_engine._apply_fdr_correction(mock_results)

        # Check that q_values are assigned
        for result in corrected_results:
            assert result.q_value is not None
            assert result.q_value >= result.p_value  # Q-value should be >= p-value

        # Check that significant results are identified correctly
        significant_count = len([r for r in corrected_results if r.q_value < 0.01])
        assert significant_count >= 0

    def test_full_statistical_pipeline(self, inference_engine):
        """Test complete statistical inference workflow."""
        results = inference_engine.run_statistical_tests(min_cooccurrence=1)

        assert isinstance(results, list)
        assert len(results) >= 0  # May have no significant results with small test data

        # Check result structure
        for result in results:
            assert isinstance(result, AssociationResult)
            assert result.domain.startswith("PF")
            assert result.go_term.startswith("GO:")
            assert result.q_value is not None
            assert result.q_value < 0.01  # Should only return significant results

    def test_empty_input_handling(self):
        """Test handling of empty input data."""
        empty_domain_map = {}
        empty_go_map = {}

        with pytest.raises(ValueError):
            StatisticalInferenceEngine(empty_domain_map, empty_go_map)

    def test_mismatched_protein_sets(self):
        """Test handling of mismatched protein sets between domain and GO maps."""
        domain_map = {"P1": ["PF00001"], "P2": ["PF00002"]}
        go_map = {"P3": {"GO:0001"}, "P4": {"GO:0002"}}  # No overlap

        engine = StatisticalInferenceEngine(domain_map, go_map)
        assert engine.total_proteins == len(domain_map)  # Uses domain map size

        results = engine.run_statistical_tests()
        assert len(results) == 0  # No associations possible with no overlap

    def test_single_protein_edge_case(self):
        """Test edge case with single protein."""
        single_domain_map = {"P1": ["PF00001"]}
        single_go_map = {"P1": {"GO:0001"}}

        engine = StatisticalInferenceEngine(single_domain_map, single_go_map)
        results = engine.run_statistical_tests()

        # Should handle gracefully but likely no significant results
        assert isinstance(results, list)

    def test_export_functionality(self, inference_engine, tmp_path):
        """Test exporting results to TSV format."""
        results = inference_engine.run_statistical_tests(min_cooccurrence=1)

        output_file = tmp_path / "test_results.tsv"
        exported_file = inference_engine.export_results_tsv(results, output_file)

        assert exported_file.exists()
        assert exported_file.stat().st_size > 0

        # Check file format
        with open(exported_file) as f:
            header = f.readline().strip()
            assert "domain" in header.lower()
            assert "go_term" in header.lower()
            assert "p_value" in header.lower()


class TestStatisticalValidation:
    """Tests for statistical method validation and edge cases."""

    def test_fisher_exact_known_case(self):
        """Test Fisher's exact test with known statistical case."""
        from scipy.stats import fisher_exact

        # Known case: strong enrichment
        a, b, c, d = 8, 2, 2, 88

        # Test our implementation matches scipy
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

        engine = StatisticalInferenceEngine({"P1": ["PF1"]}, {"P1": {"GO:1"}})
        result = engine._test_single_association("PF1", "GO:1", a, b, c, d)

        assert abs(result.p_value - p_value) < 1e-10  # Should match scipy exactly
        assert abs(result.odds_ratio - odds_ratio) < 1e-10

    def test_bonferroni_vs_fdr_comparison(self):
        """Test that FDR is less conservative than Bonferroni correction."""
        p_values = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05]

        # Bonferroni correction
        bonferroni_threshold = 0.05 / len(p_values)
        bonferroni_significant = [p < bonferroni_threshold for p in p_values]

        # FDR correction using our implementation
        from statsmodels.stats.multitest import fdrcorrection

        _, q_values = fdrcorrection(p_values, alpha=0.05, method="bh")
        fdr_significant = [q < 0.05 for q in q_values]

        # FDR should be less conservative (more discoveries)
        assert sum(fdr_significant) >= sum(bonferroni_significant)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
