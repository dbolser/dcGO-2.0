"""Unit tests for the temporal validation of the surprise score."""

import math
import random

import pytest

from validation.temporal_surprise import (
    MAX_TRUSTWORTHY_Z0,
    AssociationOutcome,
    RankingComparison,
    _budget_enrichment,
    _RankedPool,
    acquisition_base_rates,
    compare_rankings,
    pool,
    propagate,
    quantile,
    score_association,
    strata_by_prediction_budget,
    strata_by_rank,
    summarise_bootstrap,
    tie_ambiguity,
    tie_break_spread,
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


class TestQuantile:
    """Linear-interpolated quantiles — the old code used a truncated index."""

    values = [float(v) for v in range(11)]  # 0..10

    def test_endpoints(self):
        assert quantile(self.values, 0.0) == 0.0
        assert quantile(self.values, 1.0) == 10.0

    def test_median(self):
        assert quantile(self.values, 0.5) == pytest.approx(5.0)

    def test_interpolates_between_order_statistics(self):
        # h = 10 * 0.25 = 2.5, so halfway between the 3rd and 4th values.
        assert quantile(self.values, 0.25) == pytest.approx(2.5)

    def test_single_value(self):
        assert quantile([3.0], 0.9) == 3.0

    def test_empty_is_nan(self):
        assert math.isnan(quantile([], 0.5))


class TestSummariseBootstrap:
    """Percentile / basic / BCa side by side, and when to refuse all three."""

    @staticmethod
    def _symmetric(point=0.0, n=2000):
        rng = random.Random(11)
        return [rng.gauss(point, 1.0) for _ in range(n)]

    def test_symmetric_distribution_makes_the_three_intervals_agree(self):
        samples = self._symmetric()
        s = summarise_bootstrap(0.0, samples, jackknife=[1.0, 1.0, 1.0])
        assert s.trustworthy
        for got, want in zip(s.percentile, s.bca):
            assert got == pytest.approx(want, abs=0.15)
        for got, want in zip(s.percentile, s.basic):
            assert got == pytest.approx(want, abs=0.15)

    def test_basic_is_the_reflection_of_percentile_about_the_point(self):
        s = summarise_bootstrap(4.0, self._symmetric(point=1.0))
        assert s.basic[0] == pytest.approx(2 * 4.0 - s.percentile[1])
        assert s.basic[1] == pytest.approx(2 * 4.0 - s.percentile[0])

    def test_biased_distribution_moves_bca_away_from_percentile(self):
        # 65% of replicates below the observed statistic: z0 = +0.385, so BCa
        # shifts both ends upward relative to the percentile interval. This is
        # the case the old single-percentile-interval reporting hid.
        samples = [float(i) for i in range(1000)]  # uniform 0..999
        point = 650.0
        s = summarise_bootstrap(point, samples, jackknife=[0.0] * 5)
        assert s.z0 == pytest.approx(0.385, abs=0.01)
        assert s.bca[0] > s.percentile[0]
        assert s.bca[1] > s.percentile[1]
        assert s.trustworthy  # biased, but not so biased as to be useless

    def test_acceleration_is_zero_without_jackknife_and_is_noted(self):
        s = summarise_bootstrap(0.0, self._symmetric())
        assert s.acceleration == 0.0
        assert "BCa reduces to BC" in s.note

    def test_skewed_jackknife_gives_nonzero_acceleration(self):
        # One association with all the influence: strongly skewed jackknife.
        jack = [0.0] * 50 + [-10.0]
        s = summarise_bootstrap(0.0, self._symmetric(), jackknife=jack)
        assert s.acceleration != 0.0
        # A positive acceleration pushes BCa's tails apart from percentile's.
        assert s.bca != s.percentile

    def test_point_outside_percentile_interval_is_refused(self):
        # Every replicate below the observed statistic by a wide margin.
        samples = [float(i) / 100 for i in range(1000)]
        s = summarise_bootstrap(50.0, samples)
        assert not s.trustworthy
        assert "OUTSIDE its own percentile interval" in s.note
        assert math.isnan(s.recommended[0]) and math.isnan(s.recommended[1])

    def test_large_bias_correction_is_refused_even_when_the_point_is_inside(self):
        # 80% of replicates below the point but still inside the interval.
        samples = [float(i) for i in range(1000)]
        s = summarise_bootstrap(800.0, samples)
        assert s.point_inside_percentile
        assert abs(s.z0) > MAX_TRUSTWORTHY_Z0
        assert not s.trustworthy

    def test_all_replicates_on_one_side_leaves_z0_undefined(self):
        s = summarise_bootstrap(-1.0, [float(i) for i in range(100)])
        assert math.isinf(s.z0)
        assert math.isnan(s.bca[0])
        assert not s.trustworthy

    def test_no_usable_replicates(self):
        s = summarise_bootstrap(1.0, [float("nan"), float("inf")])
        assert not s.trustworthy
        assert s.n_usable == 0

    def test_too_few_resamples_is_refused(self):
        s = summarise_bootstrap(0.0, self._symmetric(n=50))
        assert not s.trustworthy
        assert "usable resamples" in s.note


class TestRankedPool:
    """The fast budget machinery must match the readable reference exactly."""

    @staticmethod
    def _pool(n=60, seed=3):
        rng = random.Random(seed)
        return [
            outcome(
                n_predicted=rng.randint(1, 40),
                n_hit=rng.randint(0, 3),
                base_rate=rng.uniform(0.001, 0.05),
                surprise=rng.choice([0.0, 0.0, 0.0, 1.0, 2.0, 3.0]),
                dcgo=rng.uniform(0, 50),
                feature=f"F{i}",
            )
            for i in range(n)
        ]

    @pytest.mark.parametrize("budget", [10, 100, 400, 1000])
    def test_matches_the_reference_implementation(self, budget):
        outcomes = self._pool()
        for ranking in ("surprise", "dcgo"):
            ranked = _RankedPool(outcomes, ranking)
            assert ranked.enrichment(budget) == pytest.approx(
                _budget_enrichment(outcomes, ranking, budget)
            )

    def test_counts_of_one_reproduce_the_observed_slice(self):
        outcomes = self._pool()
        ranked = _RankedPool(outcomes, "dcgo")
        ones = [1] * len(outcomes)
        assert ranked.enrichment_from_counts(ones, 300) == pytest.approx(
            ranked.enrichment(300)
        )

    def test_a_duplicated_association_fills_the_budget_with_copies(self):
        # The top association makes 10 predictions; at multiplicity 5 it alone
        # fills a 40-prediction budget, so nothing else is ever reached.
        outcomes = [
            outcome(10, 5, base_rate=0.1, surprise=9.0, feature="top"),
            outcome(10, 0, base_rate=0.1, surprise=1.0, feature="rest"),
        ]
        ranked = _RankedPool(outcomes, "surprise")
        assert ranked.enrichment_from_counts([5, 1], 40) == pytest.approx(
            25 / 5.0  # 5 copies: 25 hits against 5 * 10 * 0.1 expected
        )

    @pytest.mark.parametrize("budget", [10, 100, 400])
    def test_jackknife_matches_brute_force_leave_one_out(self, budget):
        outcomes = self._pool()
        for ranking in ("surprise", "dcgo"):
            fast = _RankedPool(outcomes, ranking).jackknife(budget)
            for i in range(len(outcomes)):
                brute = _budget_enrichment(
                    [o for j, o in enumerate(outcomes) if j != i], ranking, budget
                )
                assert fast[i] == pytest.approx(brute), f"{ranking}@{budget} i={i}"

    def test_ties_break_by_input_order_not_by_draw_order(self):
        # Two associations with identical scores: the first in the input wins.
        outcomes = [
            outcome(5, 1, surprise=1.0, feature="first"),
            outcome(5, 0, surprise=1.0, feature="second"),
        ]
        assert _RankedPool(outcomes, "surprise").slice_indices(1) == [0]


class TestTieAmbiguity:
    def test_a_tied_block_straddling_the_cutoff_is_flagged(self):
        outcomes = [outcome(10, 1, surprise=0.0, feature=f"F{i}") for i in range(5)]
        amb = tie_ambiguity(outcomes, "surprise", 15)
        # Two associations fit 15 predictions; all five share the score.
        assert amb.n_associations == 2
        assert amb.n_ambiguous == 2
        assert amb.n_tied_pool == 5
        assert amb.share == pytest.approx(1.0)

    def test_distinct_scores_are_never_ambiguous(self):
        outcomes = [
            outcome(10, 1, surprise=float(5 - i), feature=f"F{i}") for i in range(5)
        ]
        assert tie_ambiguity(outcomes, "surprise", 15).n_ambiguous == 0

    def test_a_slice_covering_everything_has_nothing_outside_to_tie_with(self):
        outcomes = [outcome(10, 1, surprise=0.0) for _ in range(3)]
        assert tie_ambiguity(outcomes, "surprise", 10_000).n_ambiguous == 0


class TestTieBreakSpread:
    def test_an_unambiguous_ranking_does_not_move(self):
        outcomes = [
            outcome(10, 1, surprise=float(9 - i), dcgo=float(i), feature=f"F{i}")
            for i in range(10)
        ]
        spread = tie_break_spread(outcomes, "surprise", "dcgo", 30, n_shuffles=50)
        assert spread.low == pytest.approx(spread.high)
        assert spread.median == pytest.approx(spread.observed)

    def test_a_wholly_tied_ranking_moves_a_lot(self):
        # Every surprise score is 0, so the "surprise slice" is whichever
        # associations the input order happened to put first.
        rng = random.Random(5)
        outcomes = [
            outcome(
                10,
                1 if i < 5 else 0,
                base_rate=0.01,
                surprise=0.0,
                dcgo=rng.uniform(0, 1),
                feature=f"F{i}",
            )
            for i in range(40)
        ]
        spread = tie_break_spread(outcomes, "surprise", "dcgo", 100, n_shuffles=200)
        assert spread.high > spread.low
        # The input order puts all five hits first, which is the best case.
        assert spread.observed >= spread.high


class TestCompareRankings:
    @staticmethod
    def _separable(n=200, seed=2):
        """A pool where surprise genuinely ranks and dcgo does not."""
        rng = random.Random(seed)
        outcomes = []
        for i in range(n):
            good = i < n // 4
            outcomes.append(
                outcome(
                    n_predicted=10,
                    n_hit=4 if good else 0,
                    base_rate=0.01,
                    surprise=float(n - i),
                    dcgo=rng.uniform(0, 1),
                    feature=f"F{i}",
                )
            )
        return outcomes

    def test_deterministic_for_a_seed(self):
        outcomes = self._separable()
        first = compare_rankings(
            outcomes,
            "surprise",
            "dcgo",
            300,
            n_resamples=100,
            seed=4,
            n_tie_shuffles=20,
        )
        second = compare_rankings(
            outcomes,
            "surprise",
            "dcgo",
            300,
            n_resamples=100,
            seed=4,
            n_tie_shuffles=20,
        )
        assert first.fixed.percentile == second.fixed.percentile
        assert first.reselect.bca == second.reselect.bca

    def test_a_real_difference_is_detected(self):
        outcomes = self._separable()
        c = compare_rankings(
            outcomes,
            "surprise",
            "dcgo",
            300,
            n_resamples=600,
            seed=1,
            n_tie_shuffles=50,
        )
        assert c.difference > 0
        assert c.separated
        assert c.verdict == "separated"

    def test_reports_both_designs_and_all_three_intervals(self):
        c = compare_rankings(
            self._separable(), "surprise", "dcgo", 300, n_resamples=300, seed=1
        )
        for s in (c.reselect, c.fixed):
            assert len(s.percentile) == 2
            assert len(s.basic) == 2
            assert len(s.bca) == 2

    def test_a_thin_slice_is_declared_unresolvable(self):
        # One enormous dcgo association swallows the whole budget on its own.
        outcomes = [
            outcome(1000, 5, base_rate=0.01, surprise=0.0, dcgo=99.0, feature="huge"),
        ] + [
            outcome(
                5, 1, base_rate=0.01, surprise=float(50 - i), dcgo=1.0, feature=f"F{i}"
            )
            for i in range(50)
        ]
        c = compare_rankings(
            outcomes,
            "surprise",
            "dcgo",
            100,
            n_resamples=100,
            seed=0,
            n_tie_shuffles=20,
        )
        assert c.n_associations_b == 1
        assert c.thin_basis
        assert not c.resolvable
        assert not c.separated
        assert "unresolvable" in c.verdict

    def test_a_mostly_tied_slice_is_declared_unresolvable(self):
        outcomes = [
            outcome(
                10, 1 if i % 3 == 0 else 0, base_rate=0.01, surprise=0.0, dcgo=float(i)
            )
            for i in range(100)
        ]
        c = compare_rankings(
            outcomes,
            "surprise",
            "dcgo",
            200,
            n_resamples=100,
            seed=0,
            n_tie_shuffles=20,
        )
        assert c.mostly_tied
        assert not c.resolvable
        assert "arbitrary tie-break" in c.verdict

    def test_min_associations_is_a_class_level_constant(self):
        assert RankingComparison.MIN_ASSOCIATIONS == 30
