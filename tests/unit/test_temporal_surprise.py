"""Unit tests for the temporal validation of the surprise score."""

import math

import pytest

from validation.temporal_surprise import (
    AssociationOutcome,
    acquisition_base_rates,
    pool,
    propagate,
    score_association,
    strata_by_prediction_budget,
    strata_by_rank,
    top_k,
)

# A → B → C, so ancestors(C) = {B, A}.
ANCESTORS = {"C": {"B", "A"}, "B": {"A"}, "A": set()}


def get_ancestors(term):
    return ANCESTORS.get(term, set())


def outcome(n_predicted, n_hit, base_rate=0.01, surprise=1.0, dcgo=1.0, feature="F"):
    return AssociationOutcome(
        feature=feature,
        term="T",
        n_predicted=n_predicted,
        n_hit=n_hit,
        base_rate=base_rate,
        rank_scores={"surprise": surprise, "dcgo": dcgo},
    )


class TestPropagate:
    def test_adds_ancestors(self):
        assert propagate(["C"], get_ancestors) == {"A", "B", "C"}

    def test_keeps_terms_without_ancestors(self):
        assert propagate(["A"], get_ancestors) == {"A"}

    def test_empty(self):
        assert propagate([], get_ancestors) == set()


class TestAcquisitionBaseRates:
    # P1 already had the term at t0; P2 gains it; P3 and P4 do not.
    t0 = {"P1": {"T"}}
    t1 = {"P1": {"T"}, "P2": {"T"}}
    universe = {"P1", "P2", "P3", "P4"}

    def test_excludes_proteins_that_already_had_the_term(self):
        # Eligible = P2, P3, P4; one of them gained it.
        rates = acquisition_base_rates(["T"], self.t0, self.t1, self.universe)
        assert rates["T"] == pytest.approx(1 / 3)

    def test_term_nobody_gains_has_zero_rate(self):
        rates = acquisition_base_rates(["OTHER"], self.t0, self.t1, self.universe)
        assert rates["OTHER"] == 0.0

    def test_empty_universe(self):
        assert acquisition_base_rates(["T"], {}, {}, set())["T"] == 0.0


class TestScoreAssociation:
    t0 = {"P1": {"T"}, "P2": set(), "P3": set()}
    t1 = {"P1": {"T"}, "P2": {"T"}, "P3": set()}

    def test_proteins_already_annotated_are_not_predictions(self):
        result = score_association(
            "F", "T", {"P1", "P2", "P3"}, self.t0, self.t1, 0.01, {}
        )
        # P1 knew it at t0, so only P2 and P3 are predictions.
        assert result.n_predicted == 2

    def test_hits_are_the_predictions_that_came_true(self):
        result = score_association(
            "F", "T", {"P1", "P2", "P3"}, self.t0, self.t1, 0.01, {}
        )
        assert result.n_hit == 1
        assert result.hit_rate == pytest.approx(0.5)

    def test_no_carriers_makes_no_predictions(self):
        result = score_association("F", "T", set(), self.t0, self.t1, 0.01, {})
        assert (result.n_predicted, result.n_hit, result.hit_rate) == (0, 0, 0.0)

    def test_carrier_unknown_to_both_snapshots_is_a_prediction_and_a_miss(self):
        result = score_association("F", "T", {"P9"}, self.t0, self.t1, 0.01, {})
        assert (result.n_predicted, result.n_hit) == (1, 0)


class TestEnrichment:
    def test_beating_the_base_rate(self):
        # 50% hit rate against a 10% base rate.
        assert outcome(10, 5, base_rate=0.1).enrichment == pytest.approx(5.0)

    def test_matching_the_base_rate_is_no_enrichment(self):
        assert outcome(10, 1, base_rate=0.1).enrichment == pytest.approx(1.0)

    def test_zero_base_rate_with_hits_is_infinite(self):
        assert outcome(10, 1, base_rate=0.0).enrichment == math.inf

    def test_zero_base_rate_without_hits_is_zero(self):
        assert outcome(10, 0, base_rate=0.0).enrichment == 0.0


class TestPool:
    def test_pooled_rates_weight_by_predictions(self):
        # 8/10 and 0/90: the pooled rate must not be the mean of 80% and 0%.
        result = pool("s", [outcome(10, 8), outcome(90, 0)], n_resamples=0)
        assert result.n_predicted == 100
        assert result.hit_rate == pytest.approx(0.08)

    def test_expected_rate_is_the_prediction_weighted_base_rate(self):
        result = pool(
            "s",
            [outcome(10, 1, base_rate=0.5), outcome(90, 1, base_rate=0.1)],
            n_resamples=0,
        )
        # (10*0.5 + 90*0.1) / 100
        assert result.expected_rate == pytest.approx(0.14)

    def test_enrichment_is_observed_over_expected(self):
        result = pool("s", [outcome(100, 20, base_rate=0.05)], n_resamples=0)
        assert result.enrichment == pytest.approx(4.0)

    def test_per_association_mean_ignores_infinities(self):
        result = pool(
            "s",
            [outcome(10, 5, base_rate=0.1), outcome(10, 1, base_rate=0.0)],
            n_resamples=0,
        )
        # The infinite one is dropped, leaving the 5x association.
        assert result.mean_per_association_enrichment == pytest.approx(5.0)

    def test_empty_stratum_does_not_divide_by_zero(self):
        result = pool("s", [], n_resamples=0)
        assert (result.n_associations, result.n_predicted, result.hit_rate) == (
            0,
            0,
            0.0,
        )
        assert math.isnan(result.enrichment)

    def test_bootstrap_ci_brackets_the_point_estimate(self):
        outcomes = [outcome(20, 4, base_rate=0.05) for _ in range(30)]
        result = pool("s", outcomes, n_resamples=200, seed=1)
        assert result.ci_low <= result.enrichment <= result.ci_high

    def test_bootstrap_is_deterministic_for_a_seed(self):
        outcomes = [outcome(20, i % 5, base_rate=0.05) for i in range(30)]
        first = pool("s", outcomes, n_resamples=100, seed=7)
        second = pool("s", outcomes, n_resamples=100, seed=7)
        assert (first.ci_low, first.ci_high) == (second.ci_low, second.ci_high)


class TestRanking:
    outcomes = [
        outcome(10, 1, surprise=5.0, dcgo=1.0, feature="high-surprise"),
        outcome(10, 1, surprise=1.0, dcgo=9.0, feature="high-dcgo"),
        outcome(10, 1, surprise=3.0, dcgo=5.0, feature="middling"),
    ]

    def test_top_k_uses_the_named_ranking(self):
        assert top_k(self.outcomes, "surprise", 1)[0].feature == "high-surprise"
        assert top_k(self.outcomes, "dcgo", 1)[0].feature == "high-dcgo"

    def test_top_k_beyond_the_end_returns_everything(self):
        assert len(top_k(self.outcomes, "surprise", 99)) == 3

    def test_strata_are_nested_slices_plus_the_whole_set(self):
        strata = strata_by_rank(self.outcomes, "surprise", [1, 2])
        assert [name for name, _ in strata] == [
            "surprise top-1",
            "surprise top-2",
            "surprise all (3)",
        ]
        assert [len(s) for _, s in strata] == [1, 2, 3]

    def test_cutoffs_at_or_beyond_the_set_size_are_skipped(self):
        # A "top-3" of three associations is just the whole set; no duplicate row.
        strata = strata_by_rank(self.outcomes, "surprise", [3, 10])
        assert [name for name, _ in strata] == ["surprise all (3)"]


class TestPredictionBudgetStrata:
    """Equal prediction volume, not equal association count — the fair head-to-head."""

    # Surprise favours tight architectures (few predictions each); dcgo favours
    # common domains (many). Same three associations, opposite orders.
    outcomes = [
        outcome(5, 1, surprise=9.0, dcgo=1.0, feature="tight"),
        outcome(50, 1, surprise=5.0, dcgo=5.0, feature="medium"),
        outcome(500, 1, surprise=1.0, dcgo=9.0, feature="broad"),
    ]

    def test_budget_takes_associations_until_it_is_filled(self):
        strata = strata_by_prediction_budget(self.outcomes, "surprise", [50])
        (_name, taken) = strata[0]
        # "tight" (5) then "medium" (50) crosses the budget of 50.
        assert [o.feature for o in taken] == ["tight", "medium"]

    def test_the_crossing_association_is_included_whole(self):
        (_name, taken) = strata_by_prediction_budget(self.outcomes, "surprise", [1])[0]
        assert sum(o.n_predicted for o in taken) == 5

    def test_rankings_select_different_associations_at_the_same_budget(self):
        by_surprise = strata_by_prediction_budget(self.outcomes, "surprise", [10])[0][1]
        by_dcgo = strata_by_prediction_budget(self.outcomes, "dcgo", [10])[0][1]
        assert [o.feature for o in by_surprise] == ["tight", "medium"]
        assert [o.feature for o in by_dcgo] == ["broad"]

    def test_budget_covering_everything_is_not_reported(self):
        # A stratum equal to the whole pool adds no information over "all".
        assert strata_by_prediction_budget(self.outcomes, "surprise", [10_000]) == []

    def test_names_identify_ranking_and_budget(self):
        (name, _) = strata_by_prediction_budget(self.outcomes, "dcgo", [100])[0]
        assert name == "dcgo @100 preds"
