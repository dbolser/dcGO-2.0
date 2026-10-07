"""Unit tests for the evaluator-convention check (VALIDATION_PLAN §4).

The committed numbers are only comparable with the alternatives if the
alternatives are computed on the same curves and the same bootstrap
resamples, so that is what these pin: the committed metrics come back
unchanged, the AUPRC anchor is exactly the triangle the predict-nothing point
adds, the CAFA floor only removes thresholds, and the chunked exact curve is
the one-panel curve.
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
cm = _load("check_metric_conventions")

TRUTH = {"P1": {"GO:a", "GO:b"}, "P2": {"GO:b", "GO:c"}, "P3": {"GO:a"}}
# p-score-like: each protein's weakest term scaled to exactly 0.
PREDS = {
    "P1": {"GO:a": 1.0, "GO:b": 0.4, "GO:z": 0.0},
    "P2": {"GO:x": 1.0, "GO:b": 0.0005, "GO:c": 0.0},
    "P3": {"GO:a": 0.7, "GO:c": 0.0},
}
IC = {"GO:a": 1.0, "GO:b": 2.0, "GO:c": 3.0, "GO:x": 5.0, "GO:z": 6.0}


def _panel(pred=PREDS):
    return rs.build_panel(pred, TRUTH, IC, tb._candidate_thresholds(pred))


class TestConventionMetrics:
    def test_committed_metrics_come_back_unchanged(self):
        panel = _panel()
        got, ref = cm.convention_metrics(panel), rs.panel_metrics(panel)
        assert got["f_max"] == ref["f_max"]
        assert got["auprc"] == ref["auprc"]

    def test_anchor_is_the_triangle_to_the_strictest_cutoff(self):
        panel = _panel()
        curve = rs.panel_curve(panel)
        # The sentinel predicts nothing: recall 0, precision 0.
        assert curve.coverage[-1] == 0.0 and curve.recall[-1] == 0.0
        top = np.nonzero(curve.coverage > 0)[0][-1]
        got = cm.convention_metrics(panel)
        assert got["auprc_anchor"] == pytest.approx(
            curve.recall[top] * curve.precision[top] / 2
        )
        assert got["auprc_no_anchor"] == pytest.approx(
            got["auprc"] - got["auprc_anchor"]
        )

    def test_cafa_floor_drops_only_the_thresholds_below_th_step(self):
        panel = _panel()
        curve = rs.panel_curve(panel)
        got = cm.convention_metrics(panel, th_step=0.001)
        reachable = panel.thresholds >= 0.001
        assert got["f_max_cafa_floor"] == pytest.approx(curve.f[reachable].max())
        assert got["f_max_tau_cafa_floor"] >= 0.001
        # Here the best cutoff is 0 (predict every scored term), below the grid.
        assert got["f_max"] > got["f_max_cafa_floor"]
        assert (
            cm.convention_metrics(panel, th_step=0.0)["f_max_cafa_floor"]
            == (got["f_max"])
        )

    def test_resampled_rows(self):
        panel = _panel()
        idx = np.array([0, 0, 2])
        got, ref = cm.convention_metrics(panel, idx), rs.panel_metrics(panel, idx)
        assert got["f_max"] == ref["f_max"] and got["auprc"] == ref["auprc"]


class TestSweeps:
    def test_exact_thresholds_are_every_score_plus_the_sentinel(self):
        taus = cm.exact_thresholds(PREDS)
        scores = sorted({s for t in PREDS.values() for s in t.values()})
        assert taus[:-1] == scores
        assert taus[-1] == tb._candidate_thresholds(PREDS)[-1] > scores[-1]

    def test_quantile_sweep_has_no_even_grid(self):
        many = {"P1": {f"GO:t{i}": (i / 999) ** 4 for i in range(1000)}}
        taus = cm.quantile_sweep(many)
        assert len(taus) == 52  # 51 quantiles + the sentinel
        assert set(taus[:-1]) <= set(many["P1"].values())

    def test_chunked_curve_equals_one_panel(self):
        taus = cm.exact_thresholds(PREDS)
        ref = rs.panel_metrics(rs.build_panel(PREDS, TRUTH, IC, taus))
        got = cm.sweep_metrics(PREDS, TRUTH, IC, taus, chunk=2)
        for key in ("f_max", "auprc", "s_min"):
            assert got[key] == pytest.approx(ref[key], abs=1e-15), key


class TestSameResamples:
    def test_bootstrap_seed_matches_the_former_inline_derivation(self):
        rng = np.random.default_rng(0)
        per_aspect = {a: int(rng.integers(0, 2**31 - 1)) for a in ("BP", "MF", "CC")}
        for aspect in ("BP", "MF", "CC"):
            for min_ic in (0.0, 2.0, 4.0):
                assert cm.abl.bootstrap_seed(0, aspect, min_ic) == (
                    per_aspect[aspect] + int(min_ic * 1000)
                )

    def test_another_evaluate_scores_the_same_resamples(self):
        panels = {"a": _panel(), "b": _panel({"P1": {"GO:a": 0.9}})}
        ref = rs.paired_bootstrap(panels, ("f_max", "auprc"), n_replicates=25, seed=5)
        got = rs.paired_bootstrap(
            panels,
            ("f_max", "auprc", "auprc_no_anchor"),
            n_replicates=25,
            seed=5,
            evaluate=cm.convention_metrics,
        )
        for key, values in ref.items():
            np.testing.assert_array_equal(got[key], values)
