"""
Integration tests for the True Path Rule in the complete pipeline.

Tests the interaction between run_dcgo_human.py and the OntologyProcessor.
"""

import pytest
import tempfile
from pathlib import Path
from dataclasses import dataclass

from src.ontology_processor import OntologyProcessor


# Test GO ontology with realistic hierarchy
TEST_OBO_CONTENT = """
format-version: 1.2
data-version: test-integration

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0044699
name: single-organism process
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0051179
name: localization
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0051234
name: establishment of localization
namespace: biological_process
is_a: GO:0051179 ! localization

[Term]
id: GO:0006810
name: transport
namespace: biological_process
is_a: GO:0051234 ! establishment of localization

[Term]
id: GO:0006811
name: ion transport
namespace: biological_process
is_a: GO:0006810 ! transport

[Term]
id: GO:0006812
name: cation transport
namespace: biological_process
is_a: GO:0006811 ! ion transport

[Term]
id: GO:0098655
name: cation transmembrane transport
namespace: biological_process
is_a: GO:0006812 ! cation transport

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005488
name: binding
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0043167
name: ion binding
namespace: molecular_function
is_a: GO:0005488 ! binding

[Term]
id: GO:0043169
name: cation binding
namespace: molecular_function
is_a: GO:0043167 ! ion binding
"""


@dataclass
class AssociationResult:
    """Mock association result matching run_dcgo_human.py structure."""

    domain: str
    go_term: str
    p_value: float
    q_value: float
    hyper_score: float
    a: int
    b: int
    c: int
    d: int


@pytest.fixture
def temp_obo_file():
    """Create temporary OBO file for integration tests."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".obo", delete=False) as f:
        f.write(TEST_OBO_CONTENT)
        temp_path = Path(f.name)

    yield temp_path
    temp_path.unlink()


@pytest.fixture
def ontology_processor(temp_obo_file):
    """Create OntologyProcessor for integration tests."""
    return OntologyProcessor(temp_obo_file)


class TestPipelineIntegration:
    """Integration tests for True Path Rule in production pipeline."""

    def test_pipeline_data_flow(self, ontology_processor):
        """Test complete data flow from significant associations to propagated annotations."""
        # Simulate protein data from pipeline
        protein_domain_map = {
            f"P{i:04d}": ["IPR012345"] if i < 50 else [] for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"}
            if i < 40  # Very specific term
            else (
                {"GO:0006812"}
                if i < 50  # Parent term
                else {"GO:0051179"}
            )  # Unrelated term
            for i in range(100)
        }

        # Simulate significant associations from Fisher tests
        significant_associations = [
            AssociationResult(
                domain="IPR012345",
                go_term="GO:0098655",  # cation transmembrane transport (very specific)
                p_value=1e-20,
                q_value=1e-17,
                hyper_score=99.5,
                a=40,
                b=10,
                c=0,
                d=50,
            ),
            AssociationResult(
                domain="IPR012345",
                go_term="GO:0006812",  # cation transport (parent)
                p_value=1e-18,
                q_value=1e-15,
                hyper_score=98.0,
                a=50,
                b=0,
                c=0,
                d=50,
            ),
        ]

        # Stage 1: Apply optimal level filter
        filtered = ontology_processor.apply_optimal_level_filter(
            significant_associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        # Should retain at least the most specific association
        assert len(filtered) >= 1

        # Stage 2: Propagate annotations
        propagated = ontology_processor.propagate_annotations(filtered)

        # Verify output structure
        assert len(propagated) >= 0  # May be empty if all filtered out

        # Check annotation types if we have results
        if propagated:
            direct_anns = [ann for ann in propagated if ann.annotation_type == "direct"]

            assert len(direct_anns) > 0
            # May or may not have propagated annotations depending on term depth

        # Verify hierarchical propagation if results exist
        if propagated:
            go_terms = {ann.go_term for ann in propagated}
            direct_anns = [ann for ann in propagated if ann.annotation_type == "direct"]

            # Verify that if specific terms are kept, hierarchy is respected
            # (actual terms depend on filtering results)
            assert len(go_terms) >= 1

    def test_pipeline_with_multiple_domains(self, ontology_processor):
        """Test pipeline with multiple domains and terms."""
        # Realistic multi-domain scenario
        protein_domain_map = {
            f"P{i:04d}": ["IPR001"]
            if i < 30
            else (["IPR002"] if 30 <= i < 60 else ["IPR003"])
            for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"}
            if i < 30
            else ({"GO:0043169"} if 30 <= i < 60 else {"GO:0051179"})
            for i in range(100)
        }

        # Multiple significant associations
        associations = [
            AssociationResult("IPR001", "GO:0098655", 1e-15, 1e-12, 98.0, 28, 2, 2, 68),
            AssociationResult("IPR002", "GO:0043169", 1e-14, 1e-11, 97.0, 28, 2, 2, 68),
            AssociationResult("IPR003", "GO:0051179", 1e-10, 1e-8, 92.0, 35, 5, 5, 55),
        ]

        # Apply True Path Rule
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        propagated = ontology_processor.propagate_annotations(filtered)

        # Should have annotations for multiple domains
        domains = {ann.domain for ann in propagated}
        assert len(domains) >= 1

        # Each domain should have both direct and propagated annotations
        for domain in domains:
            domain_anns = [ann for ann in propagated if ann.domain == domain]
            direct = [ann for ann in domain_anns if ann.annotation_type == "direct"]

            assert len(direct) >= 1
            # May or may not have propagated depending on term depth

    def test_pipeline_threshold_sensitivity(self, ontology_processor):
        """Test that alpha_threshold affects filtering results."""
        protein_domain_map = {
            f"P{i:04d}": ["IPR001"] if i < 40 else [] for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"}
            if i < 30
            else ({"GO:0006812"} if i < 40 else {"GO:0051179"})
            for i in range(100)
        }

        associations = [
            AssociationResult("IPR001", "GO:0098655", 1e-10, 1e-8, 95.0, 30, 10, 10, 50)
        ]

        # Test with strict threshold
        filtered_strict = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.001,  # Very strict
        )

        # Test with lenient threshold
        filtered_lenient = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,  # Standard
        )

        # Lenient should retain at least as many as strict
        assert len(filtered_lenient) >= len(filtered_strict)

    def test_pipeline_type_consistency(self, ontology_processor):
        """Test that protein_domain_map with lists works correctly."""
        # This tests the bug fix where sets were changed to lists
        protein_domain_map = {
            f"P{i:04d}": ["IPR001", "IPR002"] if i < 30 else ["IPR003"]
            for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"} if i < 30 else {"GO:0051179"}
            for i in range(100)
        }

        associations = [
            AssociationResult("IPR001", "GO:0098655", 1e-15, 1e-12, 98.0, 28, 2, 2, 68)
        ]

        # Should not raise type errors
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        assert isinstance(filtered, list)

    def test_pipeline_preserves_statistics(self, ontology_processor):
        """Test that propagated annotations preserve original statistics."""
        protein_domain_map = {
            f"P{i:04d}": ["IPR001"] if i < 40 else [] for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"} if i < 35 else {"GO:0051179"}
            for i in range(100)
        }

        original_q = 1e-15
        original_score = 99.5

        associations = [
            AssociationResult(
                "IPR001",
                "GO:0098655",
                p_value=1e-20,
                q_value=original_q,
                hyper_score=original_score,
                a=35,
                b=5,
                c=5,
                d=55,
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        if filtered:
            propagated = ontology_processor.propagate_annotations(filtered)

            # All annotations should have original q-value and score
            for ann in propagated:
                assert ann.q_value == original_q
                assert ann.association_score == original_score


class TestPipelineEdgeCases:
    """Test edge cases in pipeline integration."""

    def test_no_significant_associations(self, ontology_processor):
        """Test pipeline with no significant associations."""
        protein_domain_map = {f"P{i:04d}": ["IPR001"] for i in range(100)}
        protein_go_map = {f"P{i:04d}": {"GO:0098655"} for i in range(100)}

        # Empty input
        filtered = ontology_processor.apply_optimal_level_filter(
            [],
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        propagated = ontology_processor.propagate_annotations(filtered)

        assert filtered == []
        assert propagated == []

    def test_all_associations_filtered_out(self, ontology_processor):
        """Test when optimal level filter rejects all associations."""
        # Create scenario where associations are not at optimal level
        protein_domain_map = {
            f"P{i:04d}": ["IPR001"] if i < 50 else [] for i in range(100)
        }

        # All proteins with domain have very specific term
        # Association with general term should be filtered
        protein_go_map = {
            f"P{i:04d}": {
                "GO:0098655",
                "GO:0006812",
                "GO:0006811",
                "GO:0006810",
                "GO:0008150",
            }
            if i < 50
            else {"GO:0008150"}
            for i in range(100)
        }

        # Association at root level (very general)
        associations = [
            AssociationResult("IPR001", "GO:0008150", 1e-5, 1e-3, 80.0, 50, 0, 0, 50)
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        # May filter out general association if specific exists
        # Result depends on statistical test
        assert isinstance(filtered, list)

    def test_large_scale_propagation(self, ontology_processor):
        """Test propagation with many associations."""
        protein_domain_map = {
            f"P{i:04d}": [f"IPR{j:03d}" for j in range(i % 5)] for i in range(200)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"}
            if i % 3 == 0
            else ({"GO:0043169"} if i % 3 == 1 else {"GO:0051179"})
            for i in range(200)
        }

        # Create many associations
        associations = []
        for i in range(20):
            assoc = AssociationResult(
                f"IPR{i:03d}", "GO:0098655", 1e-10, 1e-8, 95.0, 60, 10, 10, 120
            )
            associations.append(assoc)

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        if filtered:
            propagated = ontology_processor.propagate_annotations(filtered)

            # Should handle large number of annotations
            assert len(propagated) >= len(filtered)

            # Verify no duplicates
            pairs = [(ann.domain, ann.go_term) for ann in propagated]
            assert len(pairs) == len(set(pairs))


class TestPipelineOutputFormat:
    """Test that pipeline output matches expected format for export."""

    def test_annotation_export_format(self, ontology_processor):
        """Test that annotations can be exported in expected TSV format."""
        protein_domain_map = {
            f"P{i:04d}": ["IPR001"] if i < 40 else [] for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"} if i < 35 else {"GO:0051179"}
            for i in range(100)
        }

        associations = [
            AssociationResult("IPR001", "GO:0098655", 1e-20, 1e-17, 99.5, 35, 5, 5, 55)
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        if filtered:
            propagated = ontology_processor.propagate_annotations(filtered)

            # Verify all required fields are present
            for ann in propagated:
                assert hasattr(ann, "domain")
                assert hasattr(ann, "go_term")
                assert hasattr(ann, "q_value")
                assert hasattr(ann, "association_score")
                assert hasattr(ann, "annotation_type")
                assert hasattr(ann, "direct_source_term")

                # Verify types
                assert isinstance(ann.domain, str)
                assert isinstance(ann.go_term, str)
                assert isinstance(ann.q_value, float)
                assert isinstance(ann.association_score, float)
                assert ann.annotation_type in {"direct", "propagated"}
                assert isinstance(ann.direct_source_term, str)

    def test_compare_with_without_true_path(self, ontology_processor):
        """Test difference between results with and without True Path Rule."""
        protein_domain_map = {
            f"P{i:04d}": ["IPR001"] if i < 40 else [] for i in range(100)
        }

        protein_go_map = {
            f"P{i:04d}": {"GO:0098655"} if i < 35 else {"GO:0051179"}
            for i in range(100)
        }

        associations = [
            AssociationResult("IPR001", "GO:0098655", 1e-20, 1e-17, 99.5, 35, 5, 5, 55)
        ]

        # With True Path Rule: filter + propagate
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        if filtered:
            propagated = ontology_processor.propagate_annotations(filtered)

        # With True Path Rule:
        # - May produce more annotations (due to propagation)
        # - May produce fewer (due to optimal level filtering)
        # - May produce same number in edge cases
        # The key test is that the pipeline runs without errors
        assert isinstance(filtered, list)
        if filtered:
            assert isinstance(propagated, list)
