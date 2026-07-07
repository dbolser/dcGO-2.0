"""Regression tests for GO DAG traversal direction (issue #13).

obonet emits edges child -> parent; OntologyProcessor reverses the graph so that
ancestors are the more-general terms and predecessors are parents. Before the
fix, get_ancestors returned descendants and the True Path Rule filter tested
associations against children instead of parents.

Fixture ontology is a single chain:
    GO:root -> GO:mid -> GO:leaf   (leaf is_a mid is_a root)
plus a sibling GO:leaf2 is_a GO:mid.
"""

import tempfile
from pathlib import Path

import pytest

from src.ontology_processor import OntologyProcessor

pytestmark = pytest.mark.unit

MINI_OBO = """
format-version: 1.2
data-version: test

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: mid
namespace: biological_process
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: leaf
namespace: biological_process
is_a: GO:0000002 ! mid

[Term]
id: GO:0000004
name: leaf2
namespace: biological_process
is_a: GO:0000002 ! mid
"""

ROOT, MID, LEAF, LEAF2 = "GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004"


@pytest.fixture
def processor():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".obo", delete=False) as f:
        f.write(MINI_OBO)
        path = Path(f.name)
    try:
        yield OntologyProcessor(path)
    finally:
        path.unlink()


class TestAncestorDirection:
    def test_leaf_ancestors_are_general_terms(self, processor):
        assert processor.get_ancestors(LEAF) == {MID, ROOT}

    def test_root_has_no_ancestors(self, processor):
        assert processor.get_ancestors(ROOT) == set()

    def test_mid_ancestors(self, processor):
        assert processor.get_ancestors(MID) == {ROOT}


class TestDescendantDirection:
    def test_root_descendants_are_specific_terms(self, processor):
        assert processor.get_descendants(ROOT) == {MID, LEAF, LEAF2}

    def test_leaf_has_no_descendants(self, processor):
        assert processor.get_descendants(LEAF) == set()


class TestParentLookup:
    def test_predecessors_are_parents(self, processor):
        # The optimal-level filter relies on predecessors() being parents.
        assert set(processor.go_graph.predecessors(LEAF)) == {MID}
        assert set(processor.go_graph.predecessors(MID)) == {ROOT}
        assert set(processor.go_graph.predecessors(ROOT)) == set()

    def test_roots_have_zero_in_degree(self, processor):
        roots = [
            n
            for n in processor.go_graph.nodes()
            if processor.go_graph.in_degree(n) == 0
        ]
        assert roots == [ROOT]


class TestPropagationDirection:
    def test_propagation_goes_up_to_general_terms(self, processor):
        from dataclasses import dataclass

        @dataclass
        class Assoc:
            domain: str
            go_term: str
            q_value: float
            hyper_score: float

        anns = processor.propagate_annotations(
            [Assoc(domain="IPR001", go_term=LEAF, q_value=1e-5, hyper_score=90.0)]
        )
        propagated = {a.go_term for a in anns if a.annotation_type == "propagated"}
        direct = {a.go_term for a in anns if a.annotation_type == "direct"}
        assert direct == {LEAF}
        # Propagates UP to the general ancestors, never to the sibling leaf.
        assert propagated == {MID, ROOT}
        assert LEAF2 not in propagated
