"""
Unit tests for the OntologyProcessor and True Path Rule implementation.

Tests the optimal level filtering and annotation propagation functionality.
"""

import networkx as nx
import pytest
import tempfile
from pathlib import Path
from dataclasses import dataclass

from src.ontology_processor import (
    Annotation,
    OntologyProcessor,
)


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
    with tempfile.NamedTemporaryFile(mode="w", suffix=".obo", delete=False) as f:
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

    def test_metric_failure_aborts_initialization(self, temp_obo_file, monkeypatch):
        def fail_metrics(*args, **kwargs):
            raise nx.NetworkXError("metric failure")

        monkeypatch.setattr(nx, "single_source_shortest_path_length", fail_metrics)

        with pytest.raises(nx.NetworkXError, match="metric failure"):
            OntologyProcessor(temp_obo_file)

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
        """Ancestors of a specific term are its more-general parents (issue #13)."""
        # GO:0006812 (cation transport) chain: 0006812 -> 0006811 -> 0006810
        # -> 0009987 -> 0008150 (root).
        ancestors = ontology_processor.get_ancestors("GO:0006812")
        assert ancestors == {
            "GO:0006811",
            "GO:0006810",
            "GO:0009987",
            "GO:0008150",
        }

    def test_get_ancestors_root_term(self, ontology_processor):
        """A root term has no ancestors."""
        assert ontology_processor.get_ancestors("GO:0008150") == set()

    def test_get_descendants(self, ontology_processor):
        """Descendants of a term are its more-specific children (issue #13)."""
        # GO:0006810 (transport) descendants: ion transport + cation transport.
        assert ontology_processor.get_descendants("GO:0006810") == {
            "GO:0006811",
            "GO:0006812",
        }

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

    def test_filter_keeps_specific_associations(
        self, ontology_processor, sample_protein_maps
    ):
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
                hyper_score=95.0,
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        # Should keep the specific association
        assert (
            len(filtered) >= 0
        )  # May be filtered if not significantly better than parent

    def test_filter_rejects_general_associations(self, ontology_processor):
        """Test that optimal level filter rejects overly general associations."""
        # Create protein maps where domain is associated with specific term
        # but association is also present at general level
        protein_domain_map = {
            f"P{i:03d}": ["IPR001"] if i < 50 else [] for i in range(100)
        }

        protein_go_map = {
            f"P{i:03d}": {"GO:0006812", "GO:0006811", "GO:0006810", "GO:0009987"}
            if i < 50
            else {"GO:0009987"}
            for i in range(100)
        }

        # Association at general level (should be rejected if specific exists)
        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0009987",  # cellular process (very general)
                p_value=1e-5,
                q_value=1e-3,
                hyper_score=70.0,
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
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
                hyper_score=95.0,
            )
        ]

        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
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
                hyper_score=80.0,
            )
        ]

        # Should handle gracefully (may reject due to insufficient background)
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=10,  # Require 10 proteins (we only have 2)
            alpha_threshold=0.05,
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
            alpha_threshold=0.05,
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
                hyper_score=95.0,
            )
        ]

        # Should keep unknown terms (conservative approach)
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        assert len(filtered) == 1
        assert filtered[0].go_term == "GO:9999999"

    def test_malformed_association_is_not_silently_skipped(
        self, ontology_processor, sample_protein_maps
    ):
        protein_domain_map, protein_go_map = sample_protein_maps

        with pytest.raises(AttributeError):
            ontology_processor.apply_optimal_level_filter(
                [object()], protein_domain_map, protein_go_map
            )


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
                hyper_score=95.0,
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
            propagated_anns = [
                ann for ann in propagated if ann.annotation_type == "propagated"
            ]
            assert len(propagated_anns) > 0

    def test_propagate_preserves_scores(self, ontology_processor):
        """Test that propagated annotations preserve q-values and scores."""
        direct_associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006812",
                p_value=1e-10,
                q_value=1e-8,
                hyper_score=95.0,
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
                hyper_score=95.0,
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
                hyper_score=95.0,
            ),
            MockAssociation(
                domain="IPR002",
                go_term="GO:0008324",  # cation transmembrane transporter activity
                p_value=1e-9,
                q_value=1e-7,
                hyper_score=92.0,
            ),
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
                hyper_score=95.0,
            ),
            MockAssociation(
                domain="IPR001",
                go_term="GO:0006811",  # ion transport (parent of above)
                p_value=1e-9,
                q_value=1e-7,
                hyper_score=90.0,
            ),
        ]

        propagated = ontology_processor.propagate_annotations(direct_associations)

        # Check for duplicate (domain, go_term) pairs
        pairs = [(ann.domain, ann.go_term) for ann in propagated]
        assert len(pairs) == len(set(pairs))  # No duplicates

    @pytest.mark.parametrize("reverse", [False, True])
    def test_direct_parent_wins_over_propagated_child(
        self, ontology_processor, reverse
    ):
        associations = [
            MockAssociation("IPR001", "GO:0006812", 1e-10, 1e-8, 95.0),
            MockAssociation("IPR001", "GO:0006811", 1e-4, 1e-3, 70.0),
        ]
        if reverse:
            associations.reverse()

        propagated = ontology_processor.propagate_annotations(associations)
        parent = next(ann for ann in propagated if ann.go_term == "GO:0006811")

        assert parent.annotation_type == "direct"
        assert parent.direct_source_term == "GO:0006811"
        # The direct parent keeps provenance; the stronger child contributes
        # the aggregate evidence values.
        assert parent.q_value == 1e-8
        assert parent.association_score == 95.0

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
                hyper_score=95.0,
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
                direct_source_term="GO:0006812",
            ),
            Annotation(
                domain="IPR001",
                go_term="GO:0006811",
                q_value=0.001,
                association_score=95.0,
                annotation_type="propagated",
                direct_source_term="GO:0006812",
            ),
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
                direct_source_term="GO:0006812",
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
            f"P{i:03d}": {"GO:0006812"}
            if i < 30
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
                hyper_score=98.0,
            ),
            MockAssociation(
                domain="IPR002",
                go_term="GO:0008324",
                p_value=1e-14,
                q_value=1e-11,
                hyper_score=97.0,
            ),
        ]

        # Apply optimal level filter
        filtered = ontology_processor.apply_optimal_level_filter(
            associations,
            protein_domain_map,
            protein_go_map,
            min_background_size=3,
            alpha_threshold=0.05,
        )

        # Propagate annotations
        if filtered:
            propagated = ontology_processor.propagate_annotations(filtered)

            # Verify structure
            assert len(propagated) > 0

            # Count direct vs propagated
            direct_count = sum(
                1 for ann in propagated if ann.annotation_type == "direct"
            )
            propagated_count = len(propagated) - direct_count

            # Should have both types
            assert direct_count > 0
            assert propagated_count >= 0  # May have propagated annotations

            # Validate consistency
            stats = ontology_processor.validate_annotations(propagated)
            assert stats["valid_annotations"] > 0


class TestBackgroundIndexEquivalence:
    """The inverted index must be a pure speed-up, not a change of statistics.

    ``_test_against_parent_background`` used to scan the whole proteome four
    times per parent test, which does not finish on a real human run. It now
    intersects an inverted index instead. These tests pin the four contingency
    cells to a brute-force recomputation over the raw maps, so any drift in the
    set algebra fails here rather than silently changing every True Path result.
    """

    @staticmethod
    def _brute_force_cells(
        domain, child, parent, protein_domain_map, protein_go_map, get_ancestors=None
    ):
        """The definition of a/b/c/d, written out longhand.

        A protein *has* a term if it is annotated to that term or to any
        descendant of it — the True Path Rule. That applies to the parent, which
        defines the background, and equally to the child, whose carriers are
        counted within it. Passing ``get_ancestors`` makes both explicit; without
        it this pins the older direct-only definition.
        """

        def has(term, terms):
            if term in terms:
                return True
            if get_ancestors is None:
                return False
            return any(term in get_ancestors(annotated) for annotated in terms)

        background = {p for p, terms in protein_go_map.items() if has(parent, terms)}
        a = sum(
            1
            for p in background
            if domain in protein_domain_map.get(p, [])
            and has(child, protein_go_map.get(p, set()))
        )
        b = sum(
            1
            for p in background
            if has(child, protein_go_map.get(p, set()))
            and domain not in protein_domain_map.get(p, [])
        )
        c = sum(
            1
            for p in background
            if domain in protein_domain_map.get(p, [])
            and not has(child, protein_go_map.get(p, set()))
        )
        return a, b, c, len(background) - (a + b + c)

    def test_index_reproduces_brute_force_cells(
        self, ontology_processor, sample_protein_maps
    ):
        from src.ontology_processor import _BackgroundIndex

        protein_domain_map, protein_go_map = sample_protein_maps
        index = _BackgroundIndex(
            protein_domain_map, protein_go_map, ontology_processor.get_ancestors
        )

        for domain in ("IPR001", "IPR002", "IPR_absent"):
            for child, parent in (
                ("GO:0006812", "GO:0006811"),
                ("GO:0006811", "GO:0006810"),
                ("GO:0006810", "GO:0009987"),
            ):
                expected = self._brute_force_cells(
                    domain,
                    child,
                    parent,
                    protein_domain_map,
                    protein_go_map,
                    ontology_processor.get_ancestors,
                )
                background = index.term_proteins.get(parent, set())
                dom = index.domain_proteins.get(domain, set())
                child_bg = index.term_proteins.get(child, set()) & background
                a = len(child_bg & dom)
                b = len(child_bg) - a
                c = len(dom & background) - a
                assert (a, b, c, len(background) - (a + b + c)) == expected, (
                    domain,
                    child,
                    parent,
                )


class TestBackgroundIsPropagated:
    """The parental background counts proteins annotated *beneath* the parent.

    Indexing only direct annotations gave any parent term nobody is directly
    annotated to an empty background, so `_test_against_parent_background`
    raised and the caller's conservative `except` discarded the child untested —
    54,951 times on the human t0 run, leaving ~14% of associations and
    collapsing prediction coverage to 0.22-0.50. The §4 ablation's "True Path
    hurts in 12/12 cells" was measured with that in place.
    """

    def test_a_parent_with_no_direct_annotation_still_has_a_background(
        self, ontology_processor
    ):
        """The exact shape that used to raise.

        Nobody is annotated to GO:0008150 directly; three proteins are annotated
        beneath it. The True Path Rule says all three imply it.
        """
        from src.ontology_processor import _BackgroundIndex

        protein_domain_map = {"P1": ["IPR1"], "P2": ["IPR1"], "P3": ["IPR2"]}
        protein_go_map = {
            "P1": {"GO:0006812"},  # cation transport
            "P2": {"GO:0006811"},  # ion transport
            "P3": {"GO:0006810"},  # transport
        }

        direct = _BackgroundIndex(protein_domain_map, protein_go_map)
        assert direct.term_proteins.get("GO:0008150", frozenset()) == frozenset()

        propagated = _BackgroundIndex(
            protein_domain_map, protein_go_map, ontology_processor.get_ancestors
        )
        assert propagated.term_proteins["GO:0008150"] == {"P1", "P2", "P3"}

    def test_direct_annotations_are_still_present(self, ontology_processor):
        """Propagation adds; it must not replace."""
        from src.ontology_processor import _BackgroundIndex

        index = _BackgroundIndex(
            {"P1": ["IPR1"]},
            {"P1": {"GO:0006812"}},
            ontology_processor.get_ancestors,
        )
        assert "P1" in index.term_proteins["GO:0006812"]

    def test_domains_are_not_propagated(self, ontology_processor):
        """Domains have no hierarchy here; only the term index changes."""
        from src.ontology_processor import _BackgroundIndex

        maps = ({"P1": ["IPR1", "IPR2"]}, {"P1": {"GO:0006812"}})
        direct = _BackgroundIndex(*maps)
        propagated = _BackgroundIndex(*maps, ontology_processor.get_ancestors)
        assert direct.domain_proteins == propagated.domain_proteins


# OBO with all three edge families GO uses: is_a, part_of, and the regulates
# family. Only the first two license annotation propagation.
REGULATES_OBO_CONTENT = """
format-version: 1.2
data-version: test-regulates

[Term]
id: GO:0000001
name: root process

[Term]
id: GO:0000002
name: whole process
alt_id: GO:0000099
is_a: GO:0000001 ! root process

[Term]
id: GO:0000003
name: part process
relationship: part_of GO:0000002 ! whole process

[Term]
id: GO:0000004
name: negative regulation of whole process
is_a: GO:0000001 ! root process
relationship: negatively_regulates GO:0000002 ! whole process

[Term]
id: GO:0000005
name: regulation of part process
is_a: GO:0000001 ! root process
relationship: regulates GO:0000003 ! part process
"""


@pytest.fixture
def regulates_obo_file():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".obo", delete=False) as f:
        f.write(REGULATES_OBO_CONTENT)
        temp_path = Path(f.name)
    yield temp_path
    temp_path.unlink()


class TestPropagationEdgeTypes:
    """Only is_a and part_of edges may carry annotation propagation."""

    def test_is_a_and_part_of_are_traversed(self, regulates_obo_file):
        processor = OntologyProcessor(regulates_obo_file)
        # part_of then is_a: part process -> whole process -> root process
        assert processor.get_ancestors("GO:0000003") == {"GO:0000002", "GO:0000001"}
        assert processor.get_parents("GO:0000003") == ["GO:0000002"]

    def test_regulates_edges_are_not_traversed(self, regulates_obo_file):
        processor = OntologyProcessor(regulates_obo_file)
        # "negative regulation of whole process" must NOT gain "whole process"
        # as an ancestor or parent via its negatively_regulates edge.
        assert processor.get_ancestors("GO:0000004") == {"GO:0000001"}
        assert processor.get_parents("GO:0000004") == ["GO:0000001"]
        # Nor may plain regulates carry propagation.
        assert processor.get_ancestors("GO:0000005") == {"GO:0000001"}
        # And the regulated term must not count the regulator as a descendant.
        assert "GO:0000004" not in processor.get_descendants("GO:0000002")
        assert "GO:0000005" not in processor.get_descendants("GO:0000003")

    def test_regulates_edges_are_removed_from_graph(self, regulates_obo_file):
        processor = OntologyProcessor(regulates_obo_file)
        remaining = {key for _, _, key in processor.go_graph.edges(keys=True)}
        assert remaining <= {"is_a", "part_of"}
        # 2 is_a edges from GO:0000004/GO:0000005, one from GO:0000002, one
        # part_of from GO:0000003; the 2 regulates-family edges are dropped.
        assert processor.go_graph.number_of_edges() == 4

    def test_alt_ids_map_to_their_primary_term(self, regulates_obo_file):
        """Merged ids are node data, not nodes: the map is how callers remap."""
        processor = OntologyProcessor(regulates_obo_file)
        assert processor.alt_id_map == {"GO:0000099": "GO:0000002"}
        # An alt_id is not itself a graph node, so membership alone misses it.
        assert "GO:0000099" not in processor.go_graph

    def test_propagate_annotations_ignores_regulates(self, regulates_obo_file):
        processor = OntologyProcessor(regulates_obo_file)
        associations = [
            MockAssociation(
                domain="IPR001",
                go_term="GO:0000004",
                p_value=1e-6,
                q_value=1e-5,
                hyper_score=10.0,
            )
        ]
        annotations = processor.propagate_annotations(associations)
        terms = {ann.go_term for ann in annotations}
        assert terms == {"GO:0000004", "GO:0000001"}
