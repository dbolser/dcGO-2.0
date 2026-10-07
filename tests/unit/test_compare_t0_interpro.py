"""Unit tests for the cross-arm comparison of the §4 cell (VALIDATION_PLAN §4).

The arms are only comparable as a paired difference if every arm x method is
bootstrapped on the same resamples of the same proteins, so that is what these
pin: an arm against itself differs by exactly zero in every replicate, a real
difference keeps its observed sign, and arms scored on different cohorts are
refused rather than compared.
"""

import importlib.util
import sys
from pathlib import Path

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
ct = _load("compare_t0_interpro")

TRUTH = {"P1": {"GO:a", "GO:b"}, "P2": {"GO:b", "GO:c"}, "P3": {"GO:a"}}
IC = {"GO:a": 1.0, "GO:b": 2.0, "GO:c": 3.0, "GO:x": 5.0}
GOOD = {
    "P1": {"GO:a": 1.0, "GO:b": 0.8},
    "P2": {"GO:b": 0.9, "GO:c": 0.7},
    "P3": {"GO:a": 0.9, "GO:x": 0.1},
}
POOR = {"P1": {"GO:x": 1.0, "GO:a": 0.2}, "P2": {"GO:x": 0.9}, "P3": {}}


def _panel(pred, truth=TRUTH):
    return rs.build_panel(pred, truth, IC, tb._candidate_thresholds(pred))


def _rows(panels, contrast):
    return ct.compare_cell(panels, [contrast], n_replicates=50, seed=3)


def test_an_arm_against_itself_differs_by_exactly_zero():
    panels = {"x": {"m": _panel(GOOD)}, "y": {"m": _panel(GOOD)}}
    for row in _rows(panels, ("x", "y")):
        assert row["observed_diff"] == 0.0
        assert (row["diff_ci_lo"], row["diff_ci_hi"]) == (0.0, 0.0)
        assert not row["significant"]


def test_rows_name_the_contrast_and_keep_the_observed_values():
    panels = {"t0": {"m": _panel(GOOD)}, "ctrl": {"m": _panel(POOR)}}
    rows = _rows(panels, ("t0", "ctrl"))
    assert [(r["arm_a"], r["arm_b"], r["method"], r["metric"]) for r in rows] == [
        ("t0", "ctrl", "m", "f_max"),
        ("t0", "ctrl", "m", "auprc"),
    ]
    good, poor = rs.panel_metrics(_panel(GOOD)), rs.panel_metrics(_panel(POOR))
    for row in rows:
        assert row["observed_a"] == good[row["metric"]]
        assert row["observed_b"] == poor[row["metric"]]
        assert row["observed_diff"] > 0
        assert "method_a" not in row


def test_arms_scored_on_different_cohorts_are_refused():
    other = {p: TRUTH[p] for p in ("P1", "P2")}
    panels = {"x": {"m": _panel(GOOD)}, "y": {"m": _panel(GOOD, other)}}
    with pytest.raises(ValueError, match="different cohort"):
        _rows(panels, ("x", "y"))
