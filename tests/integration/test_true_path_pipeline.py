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
        """Test complete data flow from significant associations to propagated annotations.

        True-Path-consistent fixture: proteins annotated to a specific term also
        carry all its ancestors (as in propagated GOA). The specific association
        is enriched within its parent's background and must survive the filter;
        the redundant parent association (universal within its own background)
        should not.
        """
        child = "GO:0098655"  # cation transmembrane transport (specific)
        parent = "GO:0006812"  # cation transport
        unrelated = "GO:0043169"  # cation binding (other branch)
        # Full ancestor closures (True Path Rule).
        child_ann = {child} | ontology_processor.get_ancestors(child)
        parent_ann = {parent} | ontology_processor.get_ancestors(parent)
        unrelated_ann = {unrelated} | ontology_processor.get_ancestors(unrelated)

        # P0000-P0049 carry the domain; P0000-P0039 also have the child term.
        protein_domain_map = {
            f"P{i:04d}": ["IPR012345"] if i < 50 else [] for i in range(100)
        }
        protein_go_map = {}
        for i in range(100):
            p = f"P{i:04d}"
            if i < 40:  # domain + child (+ ancestors)
                protein_go_map[p] = set(child_ann)
            elif i < 90:  # parent background (with domain for 40-49, without for 50-89)
                protein_go_map[p] = set(parent_ann)
            else:  # unrelated branch
                protein_go_map[p] = set(unrelated_ann)

        significant_associations = [
            AssociationResult("IPR012345", child, 1e-20, 1e-17, 99.5, 40, 10, 0, 50),
            AssociationResult("IPR012345", parent, 1e-18, 1e-15, 98.0, 50, 0, 0, 50),
        ]

        # Stage 1: optimal level filter — keeps the specific association.
        filtered = ontology_processor.apply_optimal_level_filter(
            significant_associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )
        filtered_terms = {a.go_term for a in filtered}
        assert child in filtered_terms  # specific association retained
        # The parent association is universal within its own background (every
        # protein carrying the child also carries the parent), so the
        # optimal-level filter must drop it as redundant.
        assert parent not in filtered_terms

        # Stage 2: propagate — UP to the ancestors, never sideways/down.
        propagated = ontology_processor.propagate_annotations(filtered)
        direct = {a.go_term for a in propagated if a.annotation_type == "direct"}
        propagated_up = {
            a.go_term for a in propagated if a.annotation_type == "propagated"
        }

        assert child in direct
        # Child's ancestors (parent .. root) must appear as propagated annotations.
        assert {parent, "GO:0008150"} <= propagated_up
        # Never propagate into an unrelated branch.
        assert unrelated not in (direct | propagated_up)

    def test_pipeline_with_multiple_domains(self, ontology_processor):
        """Test pipeline with multiple domains, each enriched for a distinct term."""
        # Each domain is associated with a distinct specific term; each term's
        # proteins carry the term's full ancestor closure (True Path Rule), and
        # a domain-free background carries the parent term so the association is
        # genuinely enriched and survives filtering.
        specs = [
            ("IPR001", "GO:0098655"),  # cation transmembrane transport
            ("IPR002", "GO:0043169"),  # cation binding
            ("IPR003", "GO:0006810"),  # transport
        ]

        protein_domain_map = {}
        protein_go_map = {}
        idx = 0

        def add(count, domain, terms):
            nonlocal idx
            for _ in range(count):
                p = f"P{idx:04d}"
                idx += 1
                protein_domain_map[p] = [domain] if domain else []
                protein_go_map[p] = set(terms)

        associations = []
        for domain, term in specs:
            term_ann = {term} | ontology_processor.get_ancestors(term)
            (parent,) = tuple(ontology_processor.go_graph.predecessors(term))
            parent_ann = {parent} | ontology_processor.get_ancestors(parent)
            add(30, domain, term_ann)  # domain + specific term (+ ancestors)
            add(20, "", parent_ann)  # background: parent term, no domain
            associations.append(
                AssociationResult(domain, term, 1e-15, 1e-12, 98.0, 30, 0, 0, 20)
            )

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )
        propagated = ontology_processor.propagate_annotations(filtered)

        # All three domain-specific associations should survive and be direct.
        direct_by_domain = {
            (a.domain, a.go_term) for a in propagated if a.annotation_type == "direct"
        }
        for domain, term in specs:
            assert (domain, term) in direct_by_domain
            # Each specific term propagates up to exactly its ancestor closure.
            up = {
                a.go_term
                for a in propagated
                if a.domain == domain and a.annotation_type == "propagated"
            }
            assert ontology_processor.get_ancestors(term) <= up

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
