"""
Unit tests for the ontology processor module.

Tests True Path Rule implementation, optimal level filtering,
and annotation propagation using mock GO structures.
"""

import pytest
import tempfile
import networkx as nx
from pathlib import Path
from unittest.mock import patch

# Add src to path for testing
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from ontology_processor import OntologyProcessor, Annotation
from statistical_inference import AssociationResult


class TestOntologyProcessor:
    """Test suite for OntologyProcessor class."""

    @pytest.fixture
    def mock_go_graph(self):
        """Create a mock GO graph for testing."""
        # Create a simple DAG representing GO structure
        graph = nx.DiGraph()

        # Add nodes (GO terms)
        go_terms = {
            "GO:0008150": "biological_process",  # root
            "GO:0065007": "biological regulation",
            "GO:0050789": "regulation of biological process",
            "GO:0006355": "regulation of transcription",
            "GO:0045893": "positive regulation of transcription",
            "GO:0005515": "protein binding",  # molecular function branch
            "GO:0003674": "molecular_function",  # root
            "GO:0003824": "catalytic activity",
            "GO:0016740": "transferase activity",
        }

        for term_id, name in go_terms.items():
            graph.add_node(term_id, name=name)

        # Add edges (is_a relationships, child -> parent)
        edges = [
            (
                "GO:0065007",
                "GO:0008150",
            ),  # biological regulation is_a biological process
            (
                "GO:0050789",
                "GO:0065007",
            ),  # regulation of biological process is_a biological regulation
            (
                "GO:0006355",
                "GO:0050789",
            ),  # regulation of transcription is_a regulation of biological process
            (
                "GO:0045893",
                "GO:0006355",
            ),  # positive regulation of transcription is_a regulation of transcription
            ("GO:0005515", "GO:0003674"),  # protein binding is_a molecular function
            ("GO:0003824", "GO:0003674"),  # catalytic activity is_a molecular function
            (
                "GO:0016740",
                "GO:0003824",
            ),  # transferase activity is_a catalytic activity
        ]

        graph.add_edges_from(edges)
        return graph

    @pytest.fixture
    def mock_obo_file(self):
        """Create a mock OBO file for testing."""
        obo_content = """format-version: 1.2
data-version: releases/2023-01-01
saved-by: dcgo_test
default-namespace: gene_ontology

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
def: "A biological process represents a specific objective that the organism is genetically programmed to achieve."

[Term]
id: GO:0065007
name: biological regulation
namespace: biological_process
def: "Any process that modulates a measurable attribute of any biological process, quality or function."
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0006355
name: regulation of transcription
namespace: biological_process
def: "Any process that modulates the frequency, rate or extent of cellular DNA-templated transcription."
is_a: GO:0065007 ! biological regulation

[Term]
id: GO:0005515
name: protein binding
namespace: molecular_function
def: "Interacting selectively and non-covalently with any protein or protein complex."
is_a: GO:0003674 ! molecular_function
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".obo", delete=False) as f:
            f.write(obo_content)
            return Path(f.name)

    @pytest.fixture
    def ontology_processor(self, mock_obo_file):
        """Create OntologyProcessor instance for testing."""
        with patch("obonet.read_obo") as mock_read:
            # Return our mock graph instead of reading actual file
            mock_graph = self.create_mock_graph()
            mock_read.return_value = mock_graph

            processor = OntologyProcessor(mock_obo_file)
            processor.go_graph = mock_graph  # Ensure we use our mock
            return processor

    def create_mock_graph(self):
        """Helper to create consistent mock graph."""
        graph = nx.DiGraph()

        # Simple hierarchy for testing
        terms = {
            "GO:0008150": {"name": "biological_process"},
            "GO:0065007": {"name": "biological regulation"},
            "GO:0006355": {"name": "regulation of transcription"},
            "GO:0045893": {"name": "positive regulation of transcription"},
            "GO:0003674": {"name": "molecular_function"},
            "GO:0005515": {"name": "protein binding"},
        }

        for term_id, data in terms.items():
            graph.add_node(term_id, **data)

        # Add hierarchical relationships (child -> parent)
        edges = [
            ("GO:0065007", "GO:0008150"),
            ("GO:0006355", "GO:0065007"),
            ("GO:0045893", "GO:0006355"),
            ("GO:0005515", "GO:0003674"),
        ]

        graph.add_edges_from(edges)
        return graph

    def test_initialization(self, mock_obo_file):
        """Test OntologyProcessor initialization."""
        with patch("obonet.read_obo") as mock_read:
            mock_read.return_value = self.create_mock_graph()

            processor = OntologyProcessor(mock_obo_file)

            assert processor.go_graph is not None
            assert len(processor.go_graph.nodes()) > 0
            mock_read.assert_called_once()

    def test_graph_preparation(self, ontology_processor):
        """Test graph preparation and obsolete term removal."""
        # Add an obsolete term to test removal
        ontology_processor.go_graph.add_node("GO:9999999", is_obsolete=True)
        initial_count = len(ontology_processor.go_graph.nodes())

        ontology_processor._prepare_graph()

        # Obsolete term should be removed
        final_count = len(ontology_processor.go_graph.nodes())
        assert final_count == initial_count - 1
        assert "GO:9999999" not in ontology_processor.go_graph.nodes()

    def test_optimal_level_filtering_basic(self, ontology_processor):
        """Test basic optimal level filtering logic."""
        # Create mock associations for testing
        associations = [
            AssociationResult(
                domain="PF00001",
                go_term="GO:0045893",  # most specific
                a=10,
                b=5,
                c=5,
                d=80,
                p_value=0.001,
                odds_ratio=8.0,
                hyper_score=90.0,
                q_value=0.005,
            ),
            AssociationResult(
                domain="PF00001",
                go_term="GO:0006355",  # parent of above
                a=12,
                b=8,
                c=8,
                d=72,
                p_value=0.01,
                odds_ratio=4.0,
                hyper_score=70.0,
                q_value=0.02,
            ),
        ]

        # Mock protein maps for background testing
        protein_domain_map = {f"P{i}": ["PF00001"] for i in range(20)}
        protein_go_map = {f"P{i}": {"GO:0045893", "GO:0006355"} for i in range(15)}
        protein_go_map.update({f"P{i}": {"GO:0006355"} for i in range(15, 20)})

        with patch.object(
            ontology_processor, "_test_against_parent_background"
        ) as mock_test:
            # Mock that child is significantly more specific than parent
            mock_test.return_value = 0.01  # Significant difference

            filtered = ontology_processor.apply_optimal_level_filter(
                associations, protein_domain_map, protein_go_map
            )

            # Should keep the more specific term
            assert len(filtered) <= len(associations)
            if len(filtered) > 0:
                assert any(assoc.go_term == "GO:0045893" for assoc in filtered)

    def test_parent_background_testing(self, ontology_processor):
        """Test testing associations against parent term background."""
        # Create test data where domain is enriched in child but not parent
        protein_domain_map = {
            "P1": ["PF00001"],
            "P2": ["PF00001"],
            "P3": ["PF00002"],
            "P4": ["PF00002"],
            "P5": ["PF00003"],
            "P6": ["PF00003"],
        }

        protein_go_map = {
            "P1": {"GO:0045893", "GO:0006355"},  # Has both child and parent
            "P2": {"GO:0045893", "GO:0006355"},  # Has both child and parent
            "P3": {"GO:0006355"},  # Has parent only
            "P4": {"GO:0006355"},  # Has parent only
            "P5": {"GO:0006355"},  # Has parent only
            "P6": {"GO:0006355"},  # Has parent only
        }

        p_value = ontology_processor._test_against_parent_background(
            "PF00001", "GO:0045893", "GO:0006355", protein_domain_map, protein_go_map
        )

        assert 0.0 <= p_value <= 1.0
        # PF00001 should be significantly enriched in GO:0045893 within GO:0006355 background
        assert p_value < 0.05  # Should be significant

    def test_annotation_propagation(self, ontology_processor):
        """Test annotation propagation up the ontology hierarchy."""
        # Create direct associations
        direct_associations = [
            AssociationResult(
                domain="PF00001",
                go_term="GO:0045893",  # Leaf term
                a=10,
                b=5,
                c=5,
                d=80,
                p_value=0.001,
                odds_ratio=8.0,
                hyper_score=90.0,
                q_value=0.005,
            ),
            AssociationResult(
                domain="PF00002",
                go_term="GO:0005515",  # Another leaf term
                a=8,
                b=7,
                c=7,
                d=78,
                p_value=0.002,
                odds_ratio=4.0,
                hyper_score=80.0,
                q_value=0.008,
            ),
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # Should have more annotations than input (due to propagation)
        assert len(propagated) >= len(direct_associations)

        # Check that direct annotations are preserved
        direct_annotations = [
            ann for ann in propagated if ann.annotation_type == "direct"
        ]
        assert len(direct_annotations) == len(direct_associations)

        # Check that propagated annotations exist
        propagated_annotations = [
            ann for ann in propagated if ann.annotation_type == "propagated"
        ]
        assert len(propagated_annotations) > 0

        # Verify propagation follows hierarchy
        pf00001_terms = {ann.go_term for ann in propagated if ann.domain == "PF00001"}
        assert "GO:0045893" in pf00001_terms  # Direct term
        assert "GO:0006355" in pf00001_terms  # Parent term
        assert "GO:0008150" in pf00001_terms  # Root term

    def test_annotation_dataclass(self):
        """Test Annotation dataclass functionality."""
        annotation = Annotation(
            domain="PF00001",
            go_term="GO:0005515",
            q_value=0.001,
            association_score=85.0,
            annotation_type="direct",
            direct_source_term="GO:0005515",
        )

        assert annotation.domain == "PF00001"
        assert annotation.go_term == "GO:0005515"
        assert annotation.is_direct()
        assert not annotation.is_propagated()

        # Test propagated annotation
        prop_annotation = Annotation(
            domain="PF00001",
            go_term="GO:0003674",
            q_value=0.001,
            association_score=85.0,
            annotation_type="propagated",
            direct_source_term="GO:0005515",
        )

        assert not prop_annotation.is_direct()
        assert prop_annotation.is_propagated()

    def test_edge_cases(self, ontology_processor):
        """Test edge cases and error handling."""
        # Test with empty associations
        empty_result = ontology_processor.propagate_annotations([])
        assert len(empty_result) == 0

        # Test with term not in graph
        unknown_term_association = AssociationResult(
            domain="PF00001",
            go_term="GO:9999999",  # Unknown term
            a=5,
            b=5,
            c=5,
            d=85,
            p_value=0.001,
            odds_ratio=2.0,
            hyper_score=60.0,
            q_value=0.005,
        )

        # Should handle gracefully
        result = ontology_processor.propagate_annotations([unknown_term_association])
        assert len(result) >= 1  # Should include the original annotation

    def test_root_term_handling(self, ontology_processor):
        """Test handling of root terms in the ontology."""
        # Test with root term association
        root_association = AssociationResult(
            domain="PF00001",
            go_term="GO:0008150",  # Root term
            a=10,
            b=10,
            c=10,
            d=70,
            p_value=0.001,
            odds_ratio=2.0,
            hyper_score=50.0,
            q_value=0.005,
        )

        # Root terms should pass optimal level test by default
        protein_domain_map = {"P1": ["PF00001"]}
        protein_go_map = {"P1": {"GO:0008150"}}

        passes = ontology_processor._passes_optimal_level_test(
            "PF00001", "GO:0008150", protein_domain_map, protein_go_map
        )
        assert passes  # Root terms should always pass

    def test_insufficient_background_handling(self, ontology_processor):
        """Test handling when parent background is too small."""
        # Create scenario with very small parent background
        protein_domain_map = {"P1": ["PF00001"], "P2": ["PF00002"]}
        protein_go_map = {"P1": {"GO:0045893"}, "P2": {"GO:0006355"}}

        p_value = ontology_processor._test_against_parent_background(
            "PF00001", "GO:0045893", "GO:0006355", protein_domain_map, protein_go_map
        )

        # Should handle small background gracefully
        assert 0.0 <= p_value <= 1.0


class TestOntologyIntegration:
    """Integration tests for ontology processing workflows."""

    @pytest.fixture
    def realistic_test_data(self):
        """Create realistic test data for integration testing."""
        # Protein domain mappings
        protein_domain_map = {
            f"P{i:05d}": [f"PF{j:05d}" for j in range(1, 4)]
            for i in range(1, 101)  # 100 proteins with 3 domains each
        }

        # GO annotations with hierarchical structure
        protein_go_map = {}
        for i in range(1, 101):
            protein_id = f"P{i:05d}"
            go_terms = set()

            # Add leaf terms
            if i <= 30:
                go_terms.add("GO:0045893")  # positive regulation of transcription
            if i <= 50:
                go_terms.add("GO:0005515")  # protein binding

            # Add parent terms (following True Path Rule)
            if "GO:0045893" in go_terms:
                go_terms.update(["GO:0006355", "GO:0065007", "GO:0008150"])
            if "GO:0005515" in go_terms:
                go_terms.add("GO:0003674")

            protein_go_map[protein_id] = go_terms

        return protein_domain_map, protein_go_map

    def test_complete_ontology_workflow(self, realistic_test_data):
        """Test complete ontology processing workflow with realistic data."""
        protein_domain_map, protein_go_map = realistic_test_data

        # Create mock statistical results
        mock_associations = [
            AssociationResult(
                "PF00001", "GO:0045893", 20, 10, 10, 60, 0.001, 8.0, 90.0, 0.005
            ),
            AssociationResult(
                "PF00001", "GO:0006355", 25, 15, 15, 45, 0.01, 4.0, 70.0, 0.02
            ),
            AssociationResult(
                "PF00002", "GO:0005515", 30, 20, 20, 30, 0.001, 6.0, 85.0, 0.008
            ),
        ]

        # Mock OBO file and processor
        with tempfile.NamedTemporaryFile(suffix=".obo") as temp_obo:
            with patch("obonet.read_obo") as mock_read:
                # Create hierarchical graph
                graph = nx.DiGraph()
                graph.add_nodes_from(
                    [
                        "GO:0008150",
                        "GO:0065007",
                        "GO:0006355",
                        "GO:0045893",
                        "GO:0003674",
                        "GO:0005515",
                    ]
                )
                graph.add_edges_from(
                    [
                        ("GO:0065007", "GO:0008150"),
                        ("GO:0006355", "GO:0065007"),
                        ("GO:0045893", "GO:0006355"),
                        ("GO:0005515", "GO:0003674"),
                    ]
                )
                mock_read.return_value = graph

                processor = OntologyProcessor(Path(temp_obo.name))

                # Apply optimal level filtering
                filtered = processor.apply_optimal_level_filter(
                    mock_associations, protein_domain_map, protein_go_map
                )

                # Propagate annotations
                final_annotations = processor.propagate_annotations(filtered)

                # Verify results
                assert len(final_annotations) > 0

                # Check annotation types
                direct_count = len(
                    [a for a in final_annotations if a.annotation_type == "direct"]
                )
                propagated_count = len(
                    [a for a in final_annotations if a.annotation_type == "propagated"]
                )

                assert direct_count > 0
                assert propagated_count > 0
                assert (
                    propagated_count > direct_count
                )  # Should have more propagated than direct


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
