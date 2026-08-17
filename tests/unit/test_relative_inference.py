"""The parental-background test, now ontology-agnostic.

Relative inference is Step 2's second statistical inference in the dcGO paper,
not the True Path Rule (Step 3). It used to live inside `OntologyProcessor` and
could therefore only serve GO. These tests cover it directly, including with a
hierarchy that is a plain dict rather than an OBO graph — which is the point of
the extraction.

The rejection-counting and background-index cases were moved here from
`test_ontology_processor.py` when the implementation moved.
"""

from __future__ import annotations

import pytest

from src.relative_inference import (
    BackgroundIndex,
    union_parental_p_value,
    InsufficientBackgroundError,
    InvalidContingencyTableError,
    filter_by_parental_background,
    parental_p_value,
    relative_p_value,
    relative_p_values,
)


class Assoc:
    """Minimal stand-in for AssociationResult."""

    def __init__(self, domain: str, term: str) -> None:
        self.domain = domain
        self.go_term = term

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Assoc({self.domain}, {self.go_term})"


# A hierarchy with no OBO file behind it at all: child -> direct parents.
CHILD_TO_PARENTS = {
    # T:leaf has TWO parents, so the union background differs from any single
    # one of them. Without a multi-parent term the union-vs-per-parent
    # distinction is invisible, and 56.1% of real GO terms have several parents.
    "T:leaf": {"T:mid", "T:other"},
    "T:mid": {"T:root"},
    "T:other": {"T:root"},
}


def parents_of(term):
    return CHILD_TO_PARENTS.get(term, set())


def ancestors_of(term):
    seen = set()
    stack = list(parents_of(term))
    while stack:
        node = stack.pop()
        if node not in seen:
            seen.add(node)
            stack.extend(parents_of(node))
    return seen


@pytest.fixture
def maps():
    """11 proteins; IPR_A tracks T:leaf exactly, IPR_B is spread under T:root.

    Sized so the leaf test clears the 0.05 threshold rather than landing on it:
    with 3 T:leaf proteins in an 8-protein T:mid background the one-tailed
    Fisher p is 1/C(8,3) = 0.018. Three proteins per term also keeps every
    background at or above the default `min_background_size` of 3.
    """
    protein_domains = {
        "P1": ["IPR_A"],
        "P2": ["IPR_A"],
        "P3": ["IPR_A"],
        "P4": ["IPR_B"],
        "P5": ["IPR_B"],
        "P6": ["IPR_B"],
        "P7": ["IPR_B"],
        "P8": ["IPR_B"],
        "P9": ["IPR_B"],
        "P10": ["IPR_B"],
        "P11": ["IPR_B"],
    }
    protein_terms = {
        "P1": {"T:leaf"},
        "P2": {"T:leaf"},
        "P3": {"T:leaf"},
        "P4": {"T:mid"},
        "P5": {"T:mid"},
        "P6": {"T:mid"},
        "P7": {"T:mid"},
        "P8": {"T:mid"},
        "P9": {"T:other"},
        "P10": {"T:other"},
        "P11": {"T:other"},
    }
    return protein_domains, protein_terms


class TestWorksWithoutAnOboGraph:
    """The extraction's reason for existing."""

    def test_the_test_runs_on_a_dict_hierarchy(self, maps):
        """No OBO graph anywhere: parents and ancestors are plain functions."""
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        kept = [
            term
            for domain, term in [("IPR_A", "T:leaf"), ("IPR_B", "T:other")]
            if union_parental_p_value(domain, term, index, parents_of, 3) < 0.05
        ]
        # IPR_A is exactly T:leaf within its parents' union; IPR_B is spread
        # across T:root and so is not specific to T:other.
        assert kept == ["T:leaf"]

    def test_background_index_propagates_term_membership(self, maps):
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        # Nobody is annotated to T:root directly, but everything sits beneath it.
        assert index.term_proteins["T:root"] == set(protein_terms)
        assert index.term_proteins["T:mid"] == {f"P{i}" for i in range(1, 9)}

    def test_domains_are_never_propagated(self, maps):
        protein_domains, protein_terms = maps
        direct = BackgroundIndex(protein_domains, protein_terms)
        propagated = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        assert direct.domain_proteins == propagated.domain_proteins


class TestRelativePValueSemantics:
    def test_a_term_with_no_parents_scores_zero(self, maps):
        """0.0 is the identity for max(overall_p, relative_p).

        A root has no parental background, so the overall inference alone
        decides it — rather than the association being discarded untested.
        """
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        assert relative_p_value("IPR_A", "T:root", index, parents_of, 3) == 0.0

    def test_the_weakest_parent_governs(self, maps):
        """An association must survive every direct parent, so max wins."""
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)

        two_parents = {"T:leaf": {"T:mid", "T:other"}}
        singles = [
            parental_p_value("IPR_A", "T:leaf", parent, index, 3)
            for parent in ("T:mid", "T:other")
        ]
        combined = relative_p_value(
            "IPR_A", "T:leaf", index, lambda t: two_parents.get(t, set()), 3
        )
        assert combined == max(singles)

    def test_short_circuit_does_not_change_the_decision(self, maps):
        """`reject_at` skips parents; it must not flip a keep into a drop."""
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        two_parents = {"T:leaf": {"T:mid", "T:other"}}
        fn = lambda t: two_parents.get(t, set())  # noqa: E731

        full = relative_p_value("IPR_A", "T:leaf", index, fn, 3)
        short = relative_p_value("IPR_A", "T:leaf", index, fn, 3, reject_at=0.05)
        assert (full < 0.05) == (short < 0.05)

    def test_threshold_is_strict(self, maps):
        """`p >= alpha` rejects, matching the original implementation.

        A p-value of exactly 1.0 (the "no association possible" return) must be
        dropped by an alpha of 1.0 rather than kept.
        """
        protein_domains, protein_terms = maps
        kept = filter_by_parental_background(
            [Assoc("IPR_B", "T:leaf")],  # IPR_B is in no T:leaf protein: p = 1.0
            protein_domains,
            protein_terms,
            parents_fn=parents_of,
            ancestors_fn=ancestors_of,
            min_background_size=3,
            alpha_threshold=1.0,
        )
        assert kept == []


class TestRejectionAccounting:
    @pytest.mark.parametrize(
        "error, counter",
        [
            (InsufficientBackgroundError("too small"), "InsufficientBackgroundError"),
            (
                InvalidContingencyTableError("invalid cells"),
                "InvalidContingencyTableError",
            ),
        ],
    )
    def test_expected_failures_are_counted_by_type(
        self, maps, monkeypatch, error, counter
    ):
        protein_domains, protein_terms = maps

        def reject(*args, **kwargs):
            raise error

        monkeypatch.setattr("src.relative_inference.parental_p_value", reject)

        p_values, rejections = relative_p_values(
            [Assoc("IPR_A", "T:leaf")],
            protein_domains,
            protein_terms,
            parents_of,
            ancestors_of,
        )
        assert p_values == [None]
        assert rejections == {counter: 1}

    def test_an_untestable_association_is_dropped_not_kept(self, maps, monkeypatch):
        """Conservative: what cannot be shown specific does not survive."""
        protein_domains, protein_terms = maps

        def reject(*args, **kwargs):
            raise InsufficientBackgroundError("too small")

        monkeypatch.setattr("src.relative_inference.parental_p_value", reject)

        assert (
            filter_by_parental_background(
                [Assoc("IPR_A", "T:leaf")],
                protein_domains,
                protein_terms,
                parents_fn=parents_of,
                ancestors_fn=ancestors_of,
            )
            == []
        )

    def test_unexpected_errors_propagate(self, maps, monkeypatch):
        protein_domains, protein_terms = maps

        def fail_unexpectedly(*args, **kwargs):
            raise RuntimeError("implementation defect")

        monkeypatch.setattr(
            "src.relative_inference.parental_p_value", fail_unexpectedly
        )

        with pytest.raises(RuntimeError, match="implementation defect"):
            relative_p_values(
                [Assoc("IPR_A", "T:leaf")],
                protein_domains,
                protein_terms,
                parents_of,
                ancestors_of,
            )

    def test_malformed_association_is_not_silently_skipped(self, maps):
        protein_domains, protein_terms = maps
        with pytest.raises(AttributeError):
            filter_by_parental_background(
                [object()],
                protein_domains,
                protein_terms,
                parents_fn=parents_of,
                ancestors_fn=ancestors_of,
            )


class TestGuards:
    def test_small_background_raises(self, maps):
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        with pytest.raises(InsufficientBackgroundError):
            parental_p_value("IPR_A", "T:leaf", "T:leaf", index, min_background_size=99)

    def test_empty_association_list_returns_empty(self, maps):
        protein_domains, protein_terms = maps
        assert (
            filter_by_parental_background(
                [], protein_domains, protein_terms, parents_of, ancestors_of
            )
            == []
        )

    def test_empty_maps_are_rejected(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            filter_by_parental_background(
                [Assoc("IPR_A", "T:leaf")], {}, {}, parents_of, ancestors_of
            )


class TestVectorisedMatchesTheLoop:
    """The pre-BH path must compute the same p-values as the post-hoc one.

    The loop tests a handful of already-significant associations with scipy and
    set algebra; the vectorised path tests every co-occurring pair with sparse
    matmuls. They must agree exactly, or "combine before BH" is not the same
    statistic the filter was applying.
    """

    @staticmethod
    def _vectorised(protein_domains, protein_terms, min_background_size=3):
        from src.relative_inference import compute_relative_p_values
        from src.sparse_fisher import (
            build_sparse_matrices,
            compute_cooccurring_contingency_tables,
        )

        domain_list = sorted({d for ds in protein_domains.values() for d in ds})
        term_list = sorted({t for ts in protein_terms.values() for t in ts})
        pdm, ptm, _ = build_sparse_matrices(
            {p: set(ds) for p, ds in protein_domains.items()},
            protein_terms,
            domain_list,
            term_list,
        )
        _, pair_index, _ = compute_cooccurring_contingency_tables(pdm, ptm)
        relative_p, _tables, rejections = compute_relative_p_values(
            pdm,
            ptm,
            pair_index,
            term_list,
            parents_of,
            ancestors_of,
            min_background_size=min_background_size,
        )
        n_terms = len(term_list)
        by_pair = {
            (domain_list[int(idx) // n_terms], term_list[int(idx) % n_terms]): p
            for idx, p in zip(pair_index, relative_p)
        }
        return by_pair, rejections

    @staticmethod
    def _loop(protein_domains, protein_terms, pairs, min_background_size=3):
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        return {
            (domain, term): union_parental_p_value(
                domain, term, index, parents_of, min_background_size
            )
            for domain, term in pairs
        }

    def test_every_cooccurring_pair_agrees(self, maps):
        protein_domains, protein_terms = maps
        vectorised, _ = self._vectorised(protein_domains, protein_terms)
        looped = self._loop(protein_domains, protein_terms, vectorised)

        assert vectorised.keys() == looped.keys()
        for pair, p in sorted(vectorised.items()):
            assert p == pytest.approx(looped[pair], abs=1e-12), pair

    def test_untestable_parents_become_p_one(self, maps):
        """The vectorised spelling of the loop's conservative rejection."""
        protein_domains, protein_terms = maps
        vectorised, rejections = self._vectorised(
            protein_domains, protein_terms, min_background_size=99
        )
        assert rejections["InsufficientBackgroundError"] > 0
        # Every term here has a parent, so nothing escapes the guard.
        assert set(vectorised.values()) == {1.0}

    def test_a_root_term_scores_zero(self, maps):
        """max(overall, 0.0) leaves a parentless term to the overall inference."""
        protein_domains, protein_terms = maps
        # Annotate directly to the root so it becomes a column of its own.
        terms = dict(protein_terms)
        terms["P1"] = {"T:root"}
        vectorised, _ = self._vectorised(protein_domains, terms)
        assert vectorised[("IPR_A", "T:root")] == 0.0


class TestPropagatedTermMatrix:
    def test_ancestors_absent_from_the_annotation_map_gain_a_column(self, maps):
        """The #46 defect in matrix form: a parent nobody is annotated to.

        Nothing is annotated to T:root directly, so it is not a column of the
        annotation matrix. Without extending the axis it would have an empty
        background and every child of it would be rejected untested.
        """
        from src.relative_inference import build_propagated_term_matrix
        from src.sparse_fisher import build_sparse_matrices

        protein_domains, protein_terms = maps
        term_list = sorted({t for ts in protein_terms.values() for t in ts})
        assert "T:root" not in term_list

        _, ptm, _ = build_sparse_matrices(
            {p: set(ds) for p, ds in protein_domains.items()},
            protein_terms,
            sorted({d for ds in protein_domains.values() for d in ds}),
            term_list,
        )
        propagated, extended, index_of = build_propagated_term_matrix(
            ptm, term_list, ancestors_of
        )
        assert "T:root" in index_of
        counts = propagated.sum(axis=0).A1
        assert counts[index_of["T:root"]] == len(protein_terms)
        assert counts[index_of["T:mid"]] == 8

    def test_multiple_routes_to_one_ancestor_count_the_protein_once(self):
        """A binary matrix, not a multiplicity count."""
        from src.relative_inference import build_propagated_term_matrix
        from src.sparse_fisher import build_sparse_matrices

        # Both children of T:root, one protein annotated to both.
        protein_terms = {"P1": {"T:mid", "T:other"}}
        term_list = ["T:mid", "T:other"]
        _, ptm, _ = build_sparse_matrices(
            {"P1": {"IPR_A"}}, protein_terms, ["IPR_A"], term_list
        )
        propagated, _, index_of = build_propagated_term_matrix(
            ptm, term_list, ancestors_of
        )
        assert propagated.sum(axis=0).A1[index_of["T:root"]] == 1


class TestUnionBackgroundIsThePaperDefinition:
    """`N_pa` is the count over proteins annotated by *any* direct parent.

    The paper's Figure 1 caption defines the parental background as the union
    over direct parents, tested once. An earlier implementation ran one test per
    parent and took the maximum, which is uniformly more conservative: the union
    is a superset of every single parent, and the child count is unchanged, so
    the same association scores at least as well against it.
    """

    def test_union_is_never_more_conservative_than_per_parent_max(self, maps):
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)

        for domain in ("IPR_A", "IPR_B"):
            for term in ("T:leaf", "T:mid", "T:other"):
                if not parents_of(term):
                    continue
                union = union_parental_p_value(domain, term, index, parents_of, 3)
                per_parent = relative_p_value(domain, term, index, parents_of, 3)
                assert union <= per_parent + 1e-12, (domain, term, union, per_parent)

    def test_a_multi_parent_term_actually_differs(self, maps):
        """Guard against the fixture silently losing its multi-parent term."""
        protein_domains, protein_terms = maps
        assert len(parents_of("T:leaf")) > 1
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        union = union_parental_p_value("IPR_A", "T:leaf", index, parents_of, 3)
        per_parent = relative_p_value("IPR_A", "T:leaf", index, parents_of, 3)
        assert union < per_parent

    def test_the_background_is_the_union_not_the_intersection(self, maps):
        """The body text says "all", the caption says "any"; we follow "any"."""
        protein_domains, protein_terms = maps
        index = BackgroundIndex(protein_domains, protein_terms, ancestors_of)
        parents = list(parents_of("T:leaf"))
        union = set().union(*(index.term_proteins[p] for p in parents))
        intersection = set.intersection(*(index.term_proteins[p] for p in parents))
        assert len(union) > len(intersection)  # the fixture distinguishes them


class TestNegativeCellsNeverReachFisher:
    """Inconsistent parents_fn/ancestors_fn must not wrap through uint32.

    If the two hierarchy views are ever derived from different relations, the
    child's propagated proteins stop being a subset of the parental union and
    c = dp - a goes negative. fisher_exact casts tables to uint32, so an
    unguarded negative cell becomes ~4.29e9 — silent garbage entering BH. The
    guard must count the pair as an invalid table and hand back the documented
    conservative verdict: p = 1 and an all-zero table row.
    """

    @staticmethod
    def _run(parents_fn, ancestors_fn):
        import numpy as np

        from src.relative_inference import compute_relative_p_values
        from src.sparse_fisher import (
            build_sparse_matrices,
            compute_cooccurring_contingency_tables,
        )

        # D co-occurs with T:child on P1-P3 and with T:parent only on P4, so
        # under an ancestors_fn that refuses to propagate, a = 3 while dp = 1:
        # c = dp - a = -2 with a + c != 0, which the trivial-margin mask cannot
        # catch — only the negative-cell guard stands between this and uint32.
        protein_domains = {
            "P1": {"D"},
            "P2": {"D"},
            "P3": {"D"},
            "P4": {"D"},
            "P5": {"D2"},
        }
        protein_terms = {
            "P1": {"T:child"},
            "P2": {"T:child"},
            "P3": {"T:child"},
            "P4": {"T:parent"},
            "P5": {"T:parent"},
        }
        domain_list = sorted({d for ds in protein_domains.values() for d in ds})
        term_list = sorted({t for ts in protein_terms.values() for t in ts})
        pdm, ptm, _ = build_sparse_matrices(
            protein_domains, protein_terms, domain_list, term_list
        )
        _, pair_index, _ = compute_cooccurring_contingency_tables(pdm, ptm)
        relative_p, tables, rejections = compute_relative_p_values(
            pdm,
            ptm,
            pair_index,
            term_list,
            parents_fn,
            ancestors_fn,
            min_background_size=1,
        )
        n_terms = len(term_list)
        by_pair = {
            (domain_list[int(idx) // n_terms], term_list[int(idx) % n_terms]): (
                p,
                tables[i],
            )
            for i, (idx, p) in enumerate(zip(pair_index, relative_p))
        }
        assert np.all(tables >= 0)
        return by_pair, rejections

    def test_inconsistent_hierarchies_yield_p_one_and_a_zero_table(self):
        # parents_fn knows about T:parent; ancestors_fn refuses to propagate —
        # the inconsistency that makes c = dp - a = -2 for (D, T:child).
        by_pair, rejections = self._run(
            parents_fn=lambda t: {"T:parent"} if t == "T:child" else set(),
            ancestors_fn=lambda t: set(),
        )
        assert rejections.get("InvalidContingencyTableError", 0) >= 1
        p, table = by_pair[("D", "T:child")]
        assert p == 1.0
        assert not table.any()

    def test_consistent_hierarchies_are_untouched_by_the_guard(self):
        # The same data with an ancestors_fn that agrees with parents_fn: the
        # child's proteins land inside the parental union, every cell is
        # non-negative, and the pair is actually tested.
        by_pair, rejections = self._run(
            parents_fn=lambda t: {"T:parent"} if t == "T:child" else set(),
            ancestors_fn=lambda t: {"T:parent"} if t == "T:child" else set(),
        )
        assert rejections == {}
        p, table = by_pair[("D", "T:child")]
        assert 0.0 < p <= 1.0
        assert table.any()
