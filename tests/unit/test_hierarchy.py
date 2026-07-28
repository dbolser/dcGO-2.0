"""Unit tests for the generic term-hierarchy propagation engine."""

import gzip
from types import SimpleNamespace

import pytest

from src.hierarchy import (
    alpha_prefix_ancestors,
    closure_ancestors,
    dotted_ancestors,
    parse_obo_child_parents,
    propagate_via_ancestors,
)


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


class TestDottedAncestors:
    def test_tcdb_number_truncates_one_level_at_a_time(self):
        assert dotted_ancestors("8.A.98.1.10") == ["8.A.98.1", "8.A.98", "8.A", "8"]

    def test_two_level_id(self):
        assert dotted_ancestors("3.6") == ["3"]

    def test_top_level_id_has_no_ancestors(self):
        assert dotted_ancestors("8") == []

    def test_custom_separator(self):
        assert dotted_ancestors("a:b:c", separator=":") == ["a:b", "a"]

    def test_empty_components_ignored(self):
        assert dotted_ancestors("1..2") == ["1"]


class TestAlphaPrefixAncestors:
    def test_merops_id_yields_family_then_catalytic_type(self):
        assert alpha_prefix_ancestors("S01.151") == ["S01", "S"]

    def test_cazy_family_yields_its_class(self):
        assert alpha_prefix_ancestors("GT32") == ["GT"]

    def test_class_only_id_has_no_ancestors_of_its_own(self):
        assert alpha_prefix_ancestors("GT") == []

    def test_non_matching_id(self):
        assert alpha_prefix_ancestors("R-HSA-71384") == []


SAMPLE_OBO = """format-version: 1.2

[Term]
id: CHEBI:29105
name: zinc(2+)
is_a: CHEBI:25213 ! metal cation
relationship: part_of CHEBI:33697

[Term]
id: CHEBI:25213
name: metal cation
is_a: CHEBI:36916

[Term]
id: CHEBI:99999
name: obsolete thing
is_a: CHEBI:25213
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


class TestParseObo:
    def _write(self, tmp_path, text=SAMPLE_OBO, name="test.obo"):
        path = tmp_path / name
        path.write_text(text)
        return path

    def test_is_a_edges(self, tmp_path):
        parents = parse_obo_child_parents(self._write(tmp_path))
        assert "CHEBI:25213" in parents["CHEBI:29105"]
        assert parents["CHEBI:25213"] == {"CHEBI:36916"}

    def test_part_of_included_by_default(self, tmp_path):
        parents = parse_obo_child_parents(self._write(tmp_path))
        assert parents["CHEBI:29105"] == {"CHEBI:25213", "CHEBI:33697"}

    def test_relations_can_be_restricted(self, tmp_path):
        parents = parse_obo_child_parents(self._write(tmp_path), relations=())
        assert parents["CHEBI:29105"] == {"CHEBI:25213"}

    def test_obsolete_terms_dropped_by_default(self, tmp_path):
        parents = parse_obo_child_parents(self._write(tmp_path))
        assert "CHEBI:99999" not in parents

    def test_obsolete_terms_can_be_kept(self, tmp_path):
        parents = parse_obo_child_parents(self._write(tmp_path), include_obsolete=True)
        assert parents["CHEBI:99999"] == {"CHEBI:25213"}

    def test_typedef_stanzas_ignored(self, tmp_path):
        assert "part_of" not in parse_obo_child_parents(self._write(tmp_path))

    def test_gzipped(self, tmp_path):
        path = tmp_path / "test.obo.gz"
        with gzip.open(path, "wt") as f:
            f.write(SAMPLE_OBO)
        assert parse_obo_child_parents(path)["CHEBI:25213"] == {"CHEBI:36916"}

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_obo_child_parents(tmp_path / "nope.obo")

    def test_feeds_closure_ancestors(self, tmp_path):
        parents = parse_obo_child_parents(self._write(tmp_path))
        assert closure_ancestors(parents)("CHEBI:29105") == {
            "CHEBI:25213",
            "CHEBI:33697",
            "CHEBI:36916",
        }
