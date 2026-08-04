"""Unit tests for bootstrap CIs and permutation nulls (VALIDATION_PLAN §4).

Two things need pinning here.

1. **The panel must be the same evaluator.** ``resampling.build_panel`` +
   ``panel_metrics`` is a fast re-expression of
   ``temporal_benchmark.evaluate_aspect`` for resampling; if it ever drifts from
   the reference implementation, every bootstrap interval in the paper is
   measuring something other than the reported point estimate. The equivalence
   test below is therefore the most important test in this file.
2. **The paired machinery must actually be paired**, and the empirical p-value
   must use the ``(r+1)/(n+1)`` estimator rather than the anti-conservative
   ``r/n``.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, _ROOT / "validation" / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


tb = _load("temporal_benchmark")
rs = _load("resampling")


# --------------------------------------------------------------------------- #
# Fixtures                                                                     #
# --------------------------------------------------------------------------- #
TRUE_SETS = {
    "P1": {"GO:a", "GO:b"},
    "P2": {"GO:b", "GO:c"},
    "P3": {"GO:a"},
    "P4": {"GO:d"},  # nothing is ever predicted for P4
}

PRED_A = {
    "P1": {"GO:a": 0.9, "GO:b": 0.4, "GO:z": 0.2},
    "P2": {"GO:b": 0.7, "GO:x": 0.6},
    "P3": {"GO:a": 0.5, "GO:c": 0.1},
}

PRED_B = {
    "P1": {"GO:z": 0.8, "GO:a": 0.3},
    "P2": {"GO:c": 0.55},
    "P3": {"GO:q": 0.95},
    "P4": {"GO:d": 0.05},
}

IC = {
    "GO:a": 1.0,
    "GO:b": 2.0,
    "GO:c": 3.0,
    "GO:d": 4.0,
    "GO:x": 5.0,
    "GO:z": 6.0,
    "GO:q": 0.5,
}


def _panel(pred):
    return rs.build_panel(pred, TRUE_SETS, IC, tb._candidate_thresholds(pred))


# --------------------------------------------------------------------------- #
# The equivalence that everything else rests on                                #
# --------------------------------------------------------------------------- #
class TestPanelMatchesReferenceEvaluator:
    @pytest.mark.parametrize("pred", [PRED_A, PRED_B])
    def test_full_cohort_reproduces_evaluate_aspect(self, pred):
        reference = tb.evaluate_aspect(pred, TRUE_SETS, IC)
        panel = _panel(pred)
        got = rs.panel_metrics(panel)
        assert got["n_eval_proteins"] == reference["n_eval_proteins"]
        assert got["f_max"] == pytest.approx(reference["f_max"])
        assert got["f_max_tau"] == pytest.approx(reference["f_max_tau"])
        assert got["s_min"] == pytest.approx(reference["s_min"])
        assert got["s_min_tau"] == pytest.approx(reference["s_min_tau"])
        assert got["auprc"] == pytest.approx(reference["auprc"])

    def test_matches_reference_when_a_method_predicts_nothing(self):
        panel = rs.build_panel({}, TRUE_SETS, IC, tb._candidate_thresholds({}))
        got = rs.panel_metrics(panel)
        reference = tb.evaluate_aspect({}, TRUE_SETS, IC)
        assert got["f_max"] == pytest.approx(reference["f_max"]) == 0.0
        assert got["s_min"] == pytest.approx(reference["s_min"])
        assert got["coverage_any"] == 0.0

    def test_empty_cohort_is_all_zeros(self):
        panel = rs.build_panel(PRED_A, {}, IC, [0.0])
        got = rs.panel_metrics(panel)
        assert got["n_eval_proteins"] == 0
        assert got["f_max"] == 0.0


class TestPanelStructure:
    def test_row_order_is_sorted_and_drops_empty_truth(self):
        panel = rs.build_panel(PRED_A, {**TRUE_SETS, "P0": set()}, IC, [0.0, 1.0])
        assert panel.proteins == ("P1", "P2", "P3", "P4")

    def test_counts_at_a_known_threshold(self):
        # tau = 0.5 keeps P1:{a(0.9)}, P2:{b(0.7), x(0.6)}, P3:{a(0.5)}
        panel = rs.build_panel(PRED_A, TRUE_SETS, IC, [0.5])
        by_name = dict(zip(panel.proteins, range(len(panel.proteins))))
        assert panel.n_pred[by_name["P1"], 0] == 1
        assert panel.tp[by_name["P1"], 0] == 1
        assert panel.n_pred[by_name["P2"], 0] == 2
        assert panel.tp[by_name["P2"], 0] == 1
        # P1 missed GO:b (IC 2.0); P2 wrongly predicted GO:x (IC 5.0)
        assert panel.ru[by_name["P1"], 0] == pytest.approx(2.0)
        assert panel.mi[by_name["P2"], 0] == pytest.approx(5.0)
        # P4 has no predictions at all: all its truth IC is remaining uncertainty
        assert panel.n_pred[by_name["P4"], 0] == 0
        assert panel.ru[by_name["P4"], 0] == pytest.approx(4.0)


class TestCoverage:
    def test_coverage_any_counts_proteins_with_predictions(self):
        # 3 of the 4 benchmark proteins get any prediction from PRED_A.
        got = rs.panel_metrics(_panel(PRED_A))
        assert got["coverage_any"] == pytest.approx(0.75)

    def test_coverage_at_fmax_never_exceeds_coverage_any(self):
        got = rs.panel_metrics(_panel(PRED_A))
        assert got["coverage_at_fmax"] <= got["coverage_any"] + 1e-12


# --------------------------------------------------------------------------- #
# Bootstrap                                                                    #
# --------------------------------------------------------------------------- #
class TestPairedBootstrap:
    def test_same_resample_is_used_for_every_method(self):
        """A paired bootstrap must reuse one protein draw across methods.

        If the draws were independent, ``A - A`` (the same panel entered twice
        under two names) would have non-zero spread. Paired, it is exactly zero.
        """
        panels = {"a1": _panel(PRED_A), "a2": _panel(PRED_A), "b": _panel(PRED_B)}
        reps = rs.paired_bootstrap(panels, ("f_max",), n_replicates=50, seed=7)
        assert np.allclose(reps["a1::f_max"], reps["a2::f_max"])
        # ...and the two genuinely different methods do vary against each other.
        assert not np.allclose(reps["a1::f_max"], reps["b::f_max"])

    def test_rejects_mismatched_cohorts(self):
        good = _panel(PRED_A)
        other = rs.build_panel(PRED_A, {"P1": {"GO:a"}}, IC, [0.0, 1.0])
        with pytest.raises(ValueError, match="paired bootstrap"):
            rs.paired_bootstrap({"a": good, "b": other}, ("f_max",), n_replicates=2)

    def test_is_reproducible_from_the_seed(self):
        panels = {"a": _panel(PRED_A), "b": _panel(PRED_B)}
        one = rs.paired_bootstrap(panels, ("f_max",), n_replicates=30, seed=11)
        two = rs.paired_bootstrap(panels, ("f_max",), n_replicates=30, seed=11)
        assert np.array_equal(one["a::f_max"], two["a::f_max"])

    def test_empty_panel_mapping(self):
        assert rs.paired_bootstrap({}, ("f_max",), n_replicates=5) == {}


class TestPercentileCI:
    def test_brackets_the_middle_of_the_sample(self):
        lo, hi = rs.percentile_ci(np.arange(101), level=0.90)
        assert lo == pytest.approx(5.0)
        assert hi == pytest.approx(95.0)

    def test_empty_sample_is_nan(self):
        lo, hi = rs.percentile_ci([])
        assert np.isnan(lo) and np.isnan(hi)


class TestSummarisePaired:
    def test_identical_methods_give_a_zero_width_interval_and_p_of_one(self):
        panels = {"a1": _panel(PRED_A), "a2": _panel(PRED_A)}
        reps = rs.paired_bootstrap(panels, ("f_max",), n_replicates=40, seed=3)
        out = rs.summarise_paired(reps, "a1", "a2", "f_max", 0.5, 0.5)
        assert out["diff_ci_lo"] == 0.0 and out["diff_ci_hi"] == 0.0
        assert out["p_value"] == pytest.approx(1.0)
        assert out["significant"] is False

    def test_flags_a_separated_difference(self):
        reps = {
            "a::f_max": np.full(200, 0.5) + np.linspace(0, 0.01, 200),
            "b::f_max": np.full(200, 0.2),
        }
        out = rs.summarise_paired(reps, "a", "b", "f_max", 0.5, 0.2)
        assert out["diff_ci_lo"] > 0
        assert out["significant"] is True
        assert out["p_value"] == pytest.approx(0.0)
        assert out["observed_diff"] == pytest.approx(0.3)


# --------------------------------------------------------------------------- #
# Permutation null                                                             #
# --------------------------------------------------------------------------- #
class TestEmpiricalPValue:
    def test_uses_the_r_plus_one_over_n_plus_one_estimator(self):
        # No null sample reaches the observed value: p is 1/(n+1), never 0.
        null = np.zeros(99)
        assert rs.empirical_p_value(1.0, null) == pytest.approx(1 / 100)

    def test_counts_ties_as_exceedances(self):
        null = np.array([1.0, 1.0, 0.0, 0.0])
        assert rs.empirical_p_value(1.0, null) == pytest.approx(3 / 5)

    def test_less_alternative(self):
        null = np.array([1.0, 2.0, 3.0])
        assert rs.empirical_p_value(0.0, null, "less") == pytest.approx(1 / 4)

    def test_two_sided_is_symmetric_about_the_null_mean(self):
        null = np.array([-1.0, 0.0, 1.0])
        assert rs.empirical_p_value(5.0, null, "two-sided") == pytest.approx(
            rs.empirical_p_value(-5.0, null, "two-sided")
        )

    def test_rejects_unknown_alternative(self):
        with pytest.raises(ValueError, match="greater/less/two-sided"):
            rs.empirical_p_value(1.0, [0.0], "sideways")

    def test_empty_null_is_nan(self):
        assert np.isnan(rs.empirical_p_value(1.0, []))


class TestSummariseNull:
    def test_reports_mean_spread_interval_and_p(self):
        rng = np.random.default_rng(0)
        null = rng.normal(0.10, 0.01, size=500)
        out = rs.summarise_null(0.25, null)
        assert out["null_n"] == 500
        assert out["null_mean"] == pytest.approx(0.10, abs=0.005)
        assert out["null_ci_lo"] < out["null_mean"] < out["null_ci_hi"]
        assert out["empirical_p"] == pytest.approx(1 / 501)
        assert out["z_score"] > 5
        assert out["ratio_observed_over_null_mean"] == pytest.approx(2.5, abs=0.2)

    def test_single_sample_has_no_sd(self):
        out = rs.summarise_null(1.0, [0.5])
        assert np.isnan(out["null_sd"])
        assert np.isnan(out["z_score"])
