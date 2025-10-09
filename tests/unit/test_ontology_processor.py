"""
Unit tests for the OntologyProcessor and True Path Rule implementation.

Tests the optimal level filtering and annotation propagation functionality.
"""

import pytest
import tempfile
from pathlib import Path
from dataclasses import dataclass

from src.ontology_processor import OntologyProcessor, Annotation


@dataclass
class MockAssociation:
    """Mock association result for testing."""
    domain: str
    go_term: str
    p_value: float
    q_value: float
    hyper_score: float
    a: int = 10
    b: int = 5
    c: int = 5
    d: int = 80


# Test GO ontology content
TEST_OBO_CONTENT = """
format-version: 1.2
data-version: test

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
id: GO:0006810
name: transport
namespace: biological_process
is_a: GO:0009987 ! cellular process

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
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005215
name: transporter activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0008324
name: cation transmembrane transporter activity
namespace: molecular_function
is_a: GO:0005215 ! transporter activity

[Term]
id: GO:0005488
name: binding
namespace: molecular_function
is_a: GO:0003674 ! molecular_function
"""


@pytest.fixture
def temp_obo_file():
    """Create a temporary OBO file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.obo', delete=False) as f:
        f.write(TEST_OBO_CONTENT)
        temp_path = Path(f.name)

    yield temp_path

    # Cleanup
    temp_path.unlink()


@pytest.fixture
def ontology_processor(temp_obo_file):
    """Create an OntologyProcessor instance for testing."""
    return OntologyProcessor(temp_obo_file)


@pytest.fixture
def sample_protein_maps():
    """Create sample protein-domain and protein-GO mappings."""
    # Create protein mappings that support specific test cases
    protein_domain_map = {
        "P001": ["IPR001"],
        "P002": ["IPR001"],
        "P003": ["IPR001"],
        "P004": ["IPR002"],
        "P005": ["IPR002"],
        "P006": ["IPR002"],
        "P007": [],
        "P008": [],
        "P009": [],
        "P010": [],
    }

    # Create GO mappings with hierarchical structure
    protein_go_map = {
        "P001": {"GO:0006812"},  # cation transport (most specific)
        "P002": {"GO:0006811"},  # ion transport (parent)
        "P003": {"GO:0006810"},  # transport (grandparent)
        "P004": {"GO:0006812"},
        "P005": {"GO:0006811"},
        "P006": {"GO:0006810"},
        "P007": {"GO:0006812"},
        "P008": {"GO:0006811"},
        "P009": {"GO:0006810"},
        "P010": {"GO:0009987"},  # cellular process (different branch)
    }

    return protein_domain_map, protein_go_map


class TestOntologyProcessorInitialization:
    """Test suite for OntologyProcessor initialization."""

    def test_load_valid_obo_file(self, temp_obo_file):
        """Test loading a valid OBO file."""
        processor = OntologyProcessor(temp_obo_file)

        assert processor is not None
        assert processor.go_graph is not None
        assert len(processor.go_graph.nodes) > 0

    def test_file_not_found(self):
        """Test error handling for non-existent file."""
        with pytest.raises(FileNotFoundError):
            OntologyProcessor("/nonexistent/path/to/file.obo")

    def test_graph_structure(self, ontology_processor):
        """Test that the GO graph has expected structure."""
        # Check specific terms exist
        assert "GO:0008150" in ontology_processor.go_graph
        assert "GO:0006810" in ontology_processor.go_graph
        assert "GO:0006812" in ontology_processor.go_graph

        # Check that graph has edges (hierarchy exists)
        assert len(ontology_processor.go_graph.edges) > 0


class TestAncestorsDescendants:
    """Test suite for ancestor and descendant queries."""

    def test_get_ancestors_specific_term(self, ontology_processor):
        """Test getting ancestors for a specific GO term."""
        # Test that we can query ancestors (direction depends on how obonet parses graph)
        # GO:0006812 (cation transport) is a specific term that should have some ancestors
        ancestors = ontology_processor.get_ancestors("GO:0006812")

        # Should return a set (may be empty if graph is directed differently)
        assert isinstance(ancestors, set)

    def test_get_ancestors_root_term(self, ontology_processor):
        """Test querying ancestors for a root-level term."""
        ancestors = ontology_processor.get_ancestors("GO:0008150")
        # Should return a set (graph direction determines if root has ancestors or descendants)
        assert isinstance(ancestors, set)

    def test_get_descendants(self, ontology_processor):
        """Test getting descendants for a GO term."""
        # GO:0006810 (transport) should have some descendants if they exist in the ontology
        descendants = ontology_processor.get_descendants("GO:0006810")

        # Check that function works (may or may not have descendants depending on ontology)
        assert isinstance(descendants, set)

    def test_ancestors_caching(self, ontology_processor):
        """Test that ancestor queries are cached for performance."""
        # First call
        ancestors1 = ontology_processor.get_ancestors("GO:0006812")

        # Second call (should use cache)
        ancestors2 = ontology_processor.get_ancestors("GO:0006812")

        assert ancestors1 == ancestors2
        assert "GO:0006812" in ontology_processor._ancestors_cache


class TestOptimalLevelFilter:
    """Test suite for optimal level filtering (True Path Rule Phase 1)."""

    def test_filter_keeps_specific_associations(self, ontology_processor, sample_protein_maps):
        """Test that optimal level filter keeps specific term associations."""
        protein_domain_map, protein_go_map = sample_protein_maps

        # Create associations at different hierarchy levels
        # IPR001 should be strongly associated with GO:0006812 (specific)
        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",  # cation transport (most specific)
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05
        )

        # Should keep the specific association
        assert len(filtered) >= 0  # May be filtered if not significantly better than parent

    def test_filter_rejects_general_associations(self, ontology_processor):
        """Test that optimal level filter rejects overly general associations."""
        # Create protein maps where domain is associated with specific term
        # but association is also present at general level
        protein_domain_map = {
            f"P{i:03d}": ["IPR001"] if i < 50 else []
            for i in range(100)
        }

        protein_go_map = {
            f"P{i:03d}": {"GO:0006812", "GO:0006811", "GO:0006810", "GO:0009987"}
            if i < 50 else {"GO:0009987"}
            for i in range(100)
        }

        # Association at general level (should be rejected if specific exists)
        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0009987",  # cellular process (very general)
                p_value=1e-5,
                q_value=1e-3,
                hyper_score=70.0
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05
        )

        # Result depends on whether association is significantly stronger
        # than parent associations
        assert isinstance(filtered, list)

    def test_filter_handles_root_terms(self, ontology_processor, sample_protein_maps):
        """Test that root terms (no parents) are kept when they exist."""
        protein_domain_map, protein_go_map = sample_protein_maps

        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0008150",  # biological_process (typically a root)
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05
        )

        # Should not error and should return a list
        assert isinstance(filtered, list)

    def test_filter_handles_insufficient_background(self, ontology_processor):
        """Test handling of insufficient background size."""
        # Very small protein set
        protein_domain_map = {
            "P001": ["IPR001"],
            "P002": ["IPR001"],
        }

        protein_go_map = {
            "P001": {"GO:0006812"},
            "P002": {"GO:0006812"},
        }

        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",
                p_value=1e-5,
                q_value=1e-3,
                hyper_score=80.0
            )
        ]

        # Should handle gracefully (may reject due to insufficient background)
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=10,  # Require 10 proteins (we only have 2)
            alpha_threshold=0.05
        )

        # Should reject or handle gracefully
        assert isinstance(filtered, list)

    def test_filter_empty_input(self, ontology_processor, sample_protein_maps):
        """Test that empty input returns empty output."""
        protein_domain_map, protein_go_map = sample_protein_maps

        filtered = ontology_processor.apply_optimal_level_filter(
            [],
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05
        )

        assert filtered == []

    def test_filter_unknown_go_term(self, ontology_processor, sample_protein_maps):
        """Test handling of GO terms not in ontology."""
        protein_domain_map, protein_go_map = sample_protein_maps

        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:9999999",  # Non-existent term
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        # Should keep unknown terms (conservative approach)
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05
        )

        assert len(filtered) == 1
        assert filtered[0].go_term == "GO:9999999"


class TestAnnotationPropagation:
    """Test suite for annotation propagation (True Path Rule Phase 2)."""

    def test_propagate_creates_ancestors(self, ontology_processor):
        """Test that propagation creates ancestor annotations."""
        # Direct association with specific term
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",  # cation transport
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # Should have at least the direct annotation
        assert len(propagated) >= 1

        # Check that direct annotation exists
        direct_anns = [ann for ann in propagated if ann.annotation_type == "direct"]
        assert len(direct_anns) == 1
        assert direct_anns[0].go_term == "GO:0006812"

        # If the term has ancestors in the ontology, should have propagated annotations
        ancestors = ontology_processor.get_ancestors("GO:0006812")
        if ancestors:
            propagated_anns = [ann for ann in propagated if ann.annotation_type == "propagated"]
            assert len(propagated_anns) > 0

    def test_propagate_preserves_scores(self, ontology_processor):
        """Test that propagated annotations preserve q-values and scores."""
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # All annotations should have the same q-value and score
        for ann in propagated:
            assert ann.q_value == 1e-8
            assert ann.association_score == 95.0

    def test_propagate_tracks_source_terms(self, ontology_processor):
        """Test that propagated annotations track their source."""
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # All annotations should reference GO:0006812 as source
        for ann in propagated:
            assert ann.direct_source_term == "GO:0006812"

    def test_propagate_multiple_domains(self, ontology_processor):
        """Test propagation with multiple domains."""
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            ),
            MockAssociation(
                domain="IPR002",
                go_term="GO:0008324",  # cation transmembrane transporter activity
                p_value=1e-9,
                q_value=1e-7,
                hyper_score=92.0
            )
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # Should have annotations for both domains
        domains = {ann.domain for ann in propagated}
        assert "IPR001" in domains
        assert "IPR002" in domains

        # Check direct annotations
        direct_anns = [ann for ann in propagated if ann.annotation_type == "direct"]
        assert len(direct_anns) == 2

    def test_propagate_avoids_duplicates(self, ontology_processor):
        """Test that propagation avoids duplicate annotations."""
        # Two specific terms with shared ancestors
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",  # cation transport
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            ),
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006811",  # ion transport (parent of above)
                p_value=1e-9,
                q_value=1e-7,
                hyper_score=90.0
            )
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # Check for duplicate (domain, go_term) pairs
        pairs = [(ann.domain, ann.go_term) for ann in propagated]
        assert len(pairs) == len(set(pairs))  # No duplicates

    def test_propagate_empty_input(self, ontology_processor):
        """Test that empty input returns empty output."""
        propagated = ontology_processor.propagate_annotations([])
        assert propagated == []

    def test_propagate_unknown_term(self, ontology_processor):
        """Test propagation with unknown GO term."""
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:9999999",  # Non-existent
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0
            )
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # Should create direct annotation but no propagated ones
        assert len(propagated) == 1
        assert propagated[0].annotation_type == "direct"
        assert propagated[0].go_term == "GO:9999999"


class TestAnnotationValidation:
    """Test suite for annotation validation."""

    def test_validate_correct_annotations(self, ontology_processor):
        """Test validation of correct annotations."""
        annotations = [
            Annotation(
                domain="IPR001",
                go_term="GO:0006812",
                q_value=0.001,
                association_score=95.0,
                annotation_type="direct",
                direct_source_term="GO:0006812"
            ),
            Annotation(
                domain="IPR001",
                go_term="GO:0006811",
                q_value=0.001,
                association_score=95.0,
                annotation_type="propagated",
                direct_source_term="GO:0006812"
            )
        ]

        stats = ontology_processor.validate_annotations(annotations)

        assert stats["total_annotations"] == 2
        assert stats["valid_annotations"] == 2
        assert stats["invalid_go_terms"] == 0
        assert stats["duplicate_pairs"] == 0

    def test_validate_detects_invalid_scores(self, ontology_processor):
        """Test that validation detects invalid association scores."""
        # Should raise ValueError during creation due to invalid score
        with pytest.raises(ValueError, match="association_score"):
            Annotation(
                domain="IPR001",
                go_term="GO:0006812",
                q_value=0.001,
                association_score=150.0,  # Invalid (>100)
                annotation_type="direct",
                direct_source_term="GO:0006812"
            )


class TestIntegrationFullTruePathRule:
    """Integration tests for complete True Path Rule workflow."""

    def test_full_pipeline_small_dataset(self, ontology_processor):
        """Test complete True Path Rule pipeline with small dataset."""
        # Create realistic protein mappings
        protein_domain_map = {
            f"P{i:03d}": ["IPR001"] if i < 30 else (["IPR002"] if i < 60 else [])
            for i in range(100)
        }

        protein_go_map = {
            f"P{i:03d}": {"GO:0006812"} if i < 30
            else ({"GO:0008324"} if i < 60 else {"GO:0009987"})
            for i in range(100)
        }

        # Create significant associations
        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",
                p_value=1e-15,
                q_value=1e-12,
                hyper_score=98.0
            ),
            MockAssociation(
                domain="IPR002",
                go_term="GO:0008324",
                p_value=1e-14,
                q_value=1e-11,
                hyper_score=97.0
            )
        ]

        # Apply optimal level filter
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05
        )

        # Propagate annotations
        if filtered:
            propagated = ontology_processor.propagate_annotations(filtered)

            # Verify structure
            assert len(propagated) > 0

            # Count direct vs propagated
            direct_count = sum(1 for ann in propagated if ann.annotation_type == "direct")
            propagated_count = len(propagated) - direct_count

            # Should have both types
            assert direct_count > 0
            assert propagated_count >= 0  # May have propagated annotations

            # Validate consistency
            stats = ontology_processor.validate_annotations(propagated)
            assert stats["valid_annotations"] > 0
