"""Unit tests for the generic term-hierarchy propagation engine."""

from types import SimpleNamespace

from src.hierarchy import closure_ancestors, propagate_via_ancestors


def _assoc(domain, term, q_value, hyper_score=50.0):
    return SimpleNamespace(
        domain=domain, go_term=term, q_value=q_value, hyper_score=hyper_score
    )


class TestClosureAncestors:
    def test_transitive_chain(self):
        anc = closure_ancestors({"c": {"b"}, "b": {"a"}})
        assert anc("c") == {"a", "b"}
        assert anc("b") == {"a"}
        assert anc("a") == set()

    def test_multiple_parents_dag(self):
        anc = closure_ancestors({"x": {"p1", "p2"}, "p1": {"root"}, "p2": {"root"}})
        assert anc("x") == {"p1", "p2", "root"}

    def test_unknown_term(self):
        assert closure_ancestors({"c": {"b"}})("zzz") == set()

    def test_cycle_does_not_hang(self):
        # Not expected in real DAGs, but must terminate with a finite result.
        anc = closure_ancestors({"a": {"b"}, "b": {"a"}})
        assert isinstance(anc("a"), set)


class TestPropagateViaAncestors:
    def test_direct_plus_ancestors(self):
        anc = closure_ancestors({"c": {"b"}, "b": {"a"}})
        anns = propagate_via_ancestors([_assoc("D1", "c", 0.001)], anc)
        by_pair = {(a.domain, a.go_term): a for a in anns}
        assert by_pair[("D1", "c")].annotation_type == "direct"
        assert by_pair[("D1", "b")].annotation_type == "propagated"
        assert by_pair[("D1", "a")].annotation_type == "propagated"
        assert all(a.direct_source_term == "c" for a in anns)

    def test_no_duplicate_pairs(self):
        anc = closure_ancestors({"c": {"b"}, "b": {"a"}})
        anns = propagate_via_ancestors([_assoc("D1", "c", 0.01)], anc)
        keys = [(a.domain, a.go_term) for a in anns]
        assert len(keys) == len(set(keys))

    def test_shared_ancestor_attributed_to_most_significant(self):
        # b is both a direct term and an ancestor of c; the lower-q source wins.
        anc = closure_ancestors({"c": {"b"}, "b": {"a"}})
        anns = propagate_via_ancestors(
            [_assoc("D1", "c", 0.05), _assoc("D1", "b", 0.001)], anc
        )
        shared = [a for a in anns if a.go_term == "a"]
        assert len(shared) == 1
        assert shared[0].direct_source_term == "b"  # from the more significant chain
