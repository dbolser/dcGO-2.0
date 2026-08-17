"""The specificity summary behind VALIDATION_PLAN's chain/ancestor numbers.

Those numbers (28.6% / 55.2% / 82.4% "on a chain") were originally computed ad
hoc; ``validation/specificity_metrics.py`` makes them reproducible. The fixture
here is a three-level chain — every expected value is computable by hand.
"""

from __future__ import annotations

from validation.specificity_metrics import SpecificityMetrics, specificity_metrics

ANCESTORS = {
    "LEAF": {"MID", "ROOT"},
    "MID": {"ROOT"},
    "ROOT": set(),
}


def get_ancestors(term: str) -> set[str]:
    return ANCESTORS.get(term, set())


class TestSpecificityMetrics:
    def test_hand_computable_chain(self) -> None:
        """Domain X reports a full chain; domain Y only the leaf."""
        pairs = [
            ("X", "LEAF"),
            ("X", "MID"),
            ("X", "ROOT"),
            ("Y", "LEAF"),
        ]
        m = specificity_metrics(pairs, get_ancestors, roots=frozenset({"ROOT"}))
        assert m == SpecificityMetrics(
            n_associations=4,
            # (2 + 1 + 0 + 2) / 4
            mean_ancestors=1.25,
            # X-LEAF and X-MID have an ancestor in X's set; X-ROOT (nothing
            # above it) and Y-LEAF (Y has no ancestor rows) do not.
            on_chain_share=0.5,
            roots_present=("ROOT",),
        )

    def test_chains_never_cross_domains(self) -> None:
        """Y holding MID does not put X's LEAF on a chain."""
        m = specificity_metrics(
            [("X", "LEAF"), ("Y", "MID")], get_ancestors, roots=frozenset({"ROOT"})
        )
        assert m.on_chain_share == 0.0
        assert m.roots_present == ()

    def test_unknown_terms_count_but_cannot_chain(self) -> None:
        m = specificity_metrics(
            [("X", "NOT-IN-ONTOLOGY")], get_ancestors, roots=frozenset({"ROOT"})
        )
        assert m.n_associations == 1
        assert m.mean_ancestors == 0.0
        assert m.on_chain_share == 0.0

    def test_empty_input_is_all_zero(self) -> None:
        m = specificity_metrics([], get_ancestors)
        assert m.n_associations == 0
        assert m.mean_ancestors == 0.0
        assert m.on_chain_share == 0.0
        assert m.roots_present == ()
