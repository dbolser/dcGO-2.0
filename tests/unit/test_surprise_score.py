"""Unit tests for the surprise score (emergent domain-combination ranking)."""

import math

import pytest

from src.surprise_score import (
    NOVELTY_CURATED,
    NOVELTY_IMPLIED,
    NOVELTY_NOVEL,
    NOVELTY_REFINES,
    EmergenceEvidence,
    apply_fdr,
    emergence_pvalue,
    expected_rate,
    locate_feature_regions,
    max_pairwise_overlap,
    median,
    noisy_or,
    novelty_factor,
    overlap_fraction,
    parse_interpro2go,
    proper_subfeatures,
    score_candidate,
)


def evidence(**overrides) -> EmergenceEvidence:
    """An emergent pair: 3/3 carriers annotated, parts and background near-silent."""
    defaults = dict(
        feature="IPR000198,IPR001452",
        term="GO:0051963",
        n_feature=3,
        n_both=3,
        single_rates=(0.0, 0.0),
        part_rates=(),
        background_rate=0.001,
        q_value=1e-8,
    )
    defaults.update(overrides)
    return EmergenceEvidence(**defaults)


class TestNoisyOr:
    def test_silent_parts_predict_nothing(self):
        assert noisy_or([0.0, 0.0]) == 0.0

    def test_single_part_passes_through(self):
        assert noisy_or([0.25]) == pytest.approx(0.25)

    def test_independent_parts_combine(self):
        # 1 - 0.5*0.8 = 0.6
        assert noisy_or([0.5, 0.2]) == pytest.approx(0.6)

    def test_certain_part_saturates(self):
        assert noisy_or([1.0, 0.3]) == pytest.approx(1.0)

    def test_empty_is_zero(self):
        assert noisy_or([]) == 0.0


class TestExpectedRate:
    def test_background_floors_silent_constituents(self):
        rate, source = expected_rate(evidence())
        assert source == "background"
        assert rate == pytest.approx(0.001)

    def test_noisy_or_wins_when_constituents_informative(self):
        rate, source = expected_rate(evidence(single_rates=(0.4, 0.2)))
        assert source == "noisy_or"
        assert rate == pytest.approx(0.52)

    def test_contained_subcombination_raises_the_bar(self):
        # A triple whose contained pair already predicts the term at 0.9 is not
        # surprising, even though its single domains are silent.
        rate, source = expected_rate(evidence(part_rates=(0.9,)))
        assert source == "best_part"
        assert rate == pytest.approx(0.9)

    def test_rate_stays_inside_open_unit_interval(self):
        rate, _ = expected_rate(evidence(single_rates=(1.0,), background_rate=1.0))
        assert 0.0 < rate < 1.0


class TestEmergencePvalue:
    def test_all_carriers_annotated_against_rare_term_is_significant(self):
        assert emergence_pvalue(3, 3, 0.001) == pytest.approx(1e-9)

    def test_expected_outcome_is_not_significant(self):
        assert emergence_pvalue(1, 2, 0.5) == pytest.approx(0.75)

    def test_no_support_is_not_evidence(self):
        assert emergence_pvalue(0, 5, 0.01) == 1.0
        assert emergence_pvalue(3, 0, 0.01) == 1.0

    def test_more_support_is_more_significant(self):
        assert emergence_pvalue(5, 5, 0.01) < emergence_pvalue(2, 2, 0.01)


class TestOverlap:
    def test_disjoint_regions_do_not_overlap(self):
        assert overlap_fraction((1, 100), (200, 300)) == 0.0

    def test_abutting_regions_do_not_overlap(self):
        assert overlap_fraction((1, 100), (101, 200)) == 0.0

    def test_nested_region_overlaps_completely(self):
        # The shorter region is entirely inside the longer one.
        assert overlap_fraction((1, 300), (50, 100)) == pytest.approx(1.0)

    def test_partial_overlap_is_relative_to_the_shorter_region(self):
        assert overlap_fraction((1, 100), (51, 250)) == pytest.approx(0.5)

    def test_max_pairwise_finds_the_worst_pair(self):
        intervals = [(1, 100), (300, 400), (310, 390)]
        assert max_pairwise_overlap(intervals) == pytest.approx(1.0)

    def test_max_pairwise_of_separate_domains_is_zero(self):
        assert max_pairwise_overlap([(1, 100), (200, 300)]) == 0.0


class TestMedian:
    def test_empty(self):
        assert median([]) == 0.0

    def test_odd_length(self):
        assert median([0.1, 0.9, 0.5]) == 0.5

    def test_even_length_averages_the_middle(self):
        assert median([0.0, 0.2, 0.4, 1.0]) == pytest.approx(0.3)


class TestLocateFeatureRegions:
    domains = ["IPR001", "IPR002", "IPR003"]
    intervals = [(1, 50), (60, 120), (130, 200)]

    def test_adjacent_pair_is_located(self):
        assert locate_feature_regions(
            self.domains, self.intervals, ["IPR002", "IPR003"]
        ) == [(60, 120), (130, 200)]

    def test_non_adjacent_pair_is_not_a_supra_domain(self):
        assert (
            locate_feature_regions(self.domains, self.intervals, ["IPR001", "IPR003"])
            == []
        )

    def test_absent_domain_yields_nothing(self):
        assert locate_feature_regions(self.domains, self.intervals, ["IPR999"]) == []

    def test_longer_than_architecture(self):
        assert locate_feature_regions(["IPR001"], [(1, 50)], ["IPR001", "IPR002"]) == []


class TestProperSubfeatures:
    def test_pair_has_no_sub_combinations(self):
        assert proper_subfeatures(("A", "B")) == []

    def test_triple_yields_its_contained_pairs(self):
        assert proper_subfeatures(("A", "B", "C")) == ["A,B", "B,C"]


class TestNoveltyFactor:
    # A → B → C, so ancestors(C) = {B, A}.
    ancestors = {"C": {"B", "A"}, "B": {"A"}, "A": set()}

    def ancestors_fn(self, term):
        return self.ancestors.get(term, set())

    def test_no_curated_reference(self):
        factor, status = novelty_factor("C", set(), self.ancestors_fn)
        assert (factor, status) == (NOVELTY_NOVEL, "no-reference")

    def test_exactly_curated_is_discounted(self):
        factor, status = novelty_factor("B", {"B"}, self.ancestors_fn)
        assert (factor, status) == (NOVELTY_CURATED, "curated")

    def test_more_specific_than_curated_is_a_refinement(self):
        factor, status = novelty_factor("C", {"A"}, self.ancestors_fn)
        assert (factor, status) == (NOVELTY_REFINES, "refines")

    def test_more_general_than_curated_is_implied(self):
        factor, status = novelty_factor("A", {"C"}, self.ancestors_fn)
        assert (factor, status) == (NOVELTY_IMPLIED, "implied")

    def test_unrelated_term_is_novel(self):
        factor, status = novelty_factor("Z", {"C"}, self.ancestors_fn)
        assert (factor, status) == (NOVELTY_NOVEL, "novel")

    def test_without_ancestors_only_exact_matches_are_discounted(self):
        assert novelty_factor("C", {"A"}, None) == (NOVELTY_NOVEL, "novel")


class TestScoreCandidate:
    def test_components_are_reported(self):
        result = score_candidate(evidence(), region_overlap=0.0)
        assert result.n_both == 3
        assert result.observed_rate == pytest.approx(1.0)
        assert result.expected_rate == pytest.approx(0.001)
        assert result.lift == pytest.approx(1000.0)
        assert result.distinctness == pytest.approx(1.0)
        assert result.uninformative_constituents == 2

    def test_redundant_signatures_score_zero(self):
        result = score_candidate(evidence(), region_overlap=1.0)
        scored = apply_fdr([result])[0]
        assert scored.distinctness == 0.0
        assert scored.surprise == 0.0

    def test_curated_prediction_scores_below_an_equally_strong_novel_one(self):
        novel, curated = apply_fdr(
            [
                score_candidate(evidence(), 0.0, NOVELTY_NOVEL, "novel"),
                score_candidate(evidence(), 0.0, NOVELTY_CURATED, "curated"),
            ]
        )
        assert novel.surprise > curated.surprise

    def test_surprise_is_finite_when_q_underflows(self):
        result = score_candidate(evidence(n_feature=40, n_both=40))
        scored = apply_fdr([result])[0]
        assert math.isfinite(scored.surprise)
        assert scored.surprise > 0


class TestApplyFdr:
    def test_q_values_are_attached_and_never_below_p(self):
        results = apply_fdr(
            [
                score_candidate(evidence(n_both=3, n_feature=3)),
                score_candidate(evidence(n_both=1, n_feature=3)),
            ]
        )
        assert all(r.q_emergence >= r.p_emergence for r in results)

    def test_empty_input(self):
        assert apply_fdr([]) == []

    def test_weak_candidate_is_not_significant(self):
        results = apply_fdr(
            [score_candidate(evidence(n_both=1, n_feature=2, background_rate=0.5))]
        )
        assert results[0].q_emergence > 0.05


class TestParseInterpro2Go:
    def test_parses_mapping_lines(self):
        lines = [
            "!version date: 2026/06/01",
            "InterPro:IPR000003 Retinoid X receptor > GO:DNA binding ; GO:0003677",
            "InterPro:IPR000003 Retinoid X receptor > GO:nucleus ; GO:0005634",
            "InterPro:IPR001452 SH3 domain > GO:protein binding ; GO:0005515",
        ]
        mapping = parse_interpro2go(lines)
        assert mapping["IPR000003"] == {"GO:0003677", "GO:0005634"}
        assert mapping["IPR001452"] == {"GO:0005515"}

    def test_comments_and_junk_ignored(self):
        assert parse_interpro2go(["! comment", "", "garbage line"]) == {}
