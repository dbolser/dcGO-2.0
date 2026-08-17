"""The shared annotation-frequency IC — one convention for pipeline and benchmark.

``src.information_content.information_content`` is the counting convention
behind the pipeline's exported ``ic`` column and the ``--min-ic`` reporting
floor, and (via a propagating wrapper) the temporal benchmark's IC cells. The
fixtures are hand-computable: powers of two over a 4-protein universe, so every
expected value is exact in floating point.
"""

from __future__ import annotations

import pytest

from src.information_content import (
    information_content,
    information_content_from_term_sets,
)
from validation.temporal_benchmark import information_content as benchmark_ic


class TestHandComputableFixture:
    """4 proteins; term frequencies 4/4, 2/4, 1/4 → IC 0, 1, 2 exactly."""

    MAP = {
        "P1": {"ROOT", "MID", "LEAF"},
        "P2": {"ROOT", "MID"},
        "P3": {"ROOT"},
        "P4": {"ROOT"},
    }

    def test_frequencies_become_exact_ics(self) -> None:
        ic = information_content(self.MAP)
        assert ic == {"ROOT": 0.0, "MID": 1.0, "LEAF": 2.0}

    def test_a_universal_term_has_ic_zero_by_construction(self) -> None:
        """P(root) = 1, so any positive --min-ic floor removes DAG roots."""
        assert information_content(self.MAP)["ROOT"] == 0.0

    def test_terms_never_annotated_are_absent_not_zero(self) -> None:
        """Callers read absent terms as IC 0 via .get(term, 0.0)."""
        assert "UNSEEN" not in information_content(self.MAP)


class TestEdgeCases:
    def test_empty_map_gives_empty_ic(self) -> None:
        assert information_content({}) == {}

    def test_proteins_with_empty_term_sets_still_widen_the_universe(self) -> None:
        """A protein whose annotations were all dropped (dead ids) stays in the
        analysable universe, so it belongs in the denominator."""
        ic = information_content({"P1": {"T"}, "P2": set()})
        assert ic["T"] == 1.0  # 1/2, not 1/1

    def test_single_protein_universe_makes_every_term_ic_zero(self) -> None:
        assert information_content({"P1": {"T"}}) == {"T": 0.0}


class TestStreamingMode:
    """information_content_from_term_sets is the same estimate, streamed."""

    def test_matches_the_mapping_form(self) -> None:
        streamed = information_content_from_term_sets(
            iter(TestHandComputableFixture.MAP.values())
        )
        assert streamed == information_content(TestHandComputableFixture.MAP)

    def test_generator_input_counts_the_universe_itself(self) -> None:
        """The denominator is the number of sets yielded — no len() needed."""
        ic = information_content_from_term_sets(s for s in ({"T"}, set()))
        assert ic == {"T": 1.0}

    def test_empty_stream_gives_empty_ic(self) -> None:
        assert information_content_from_term_sets(iter(())) == {}


class TestBenchmarkWrapperSharesTheConvention:
    """The temporal benchmark's IC = shared IC over the propagated map."""

    def test_wrapper_propagates_then_counts(self) -> None:
        parents = {"LEAF": {"MID"}, "MID": {"ROOT"}}

        def ancestors(term: str) -> set[str]:
            out: set[str] = set()
            frontier = [term]
            while frontier:
                for parent in parents.get(frontier.pop(), ()):
                    if parent not in out:
                        out.add(parent)
                        frontier.append(parent)
            return out

        unpropagated = {
            "P1": {"LEAF"},
            "P2": {"MID"},
            "P3": {"ROOT"},
            "P4": {"ROOT"},
        }
        ic = benchmark_ic(unpropagated, ancestors)
        # Propagated frequencies: ROOT 4/4, MID 2/4, LEAF 1/4.
        assert ic == {"ROOT": 0.0, "MID": 1.0, "LEAF": 2.0}

        # Same map, already propagated, through the shared function directly.
        propagated = {
            "P1": {"LEAF", "MID", "ROOT"},
            "P2": {"MID", "ROOT"},
            "P3": {"ROOT"},
            "P4": {"ROOT"},
        }
        assert information_content(propagated) == ic

    def test_unpropagated_counting_would_inflate_mid_level_ic(self) -> None:
        """The defect the propagated-map requirement exists to prevent."""
        unpropagated = {
            "P1": {"LEAF"},
            "P2": {"MID"},
            "P3": {"ROOT"},
            "P4": {"ROOT"},
        }
        ic = information_content(unpropagated)  # deliberately wrong input
        assert ic["MID"] == pytest.approx(2.0)  # inflated: truth is 1.0
