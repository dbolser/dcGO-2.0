"""Unit tests for the F_max operating point in validation/ablation_pr_decomposition.py.

With no cutoff where F > 0 there is no F_max operating point, so the
decomposition must not read a sensitivity or a prediction budget off an
arbitrary cutoff (resampling.panel_metrics reports no F_max cutoff either).
"""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(
        name, _ROOT / "validation" / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


pd_ = _load("ablation_pr_decomposition")
rs, tb = pd_.rs, pd_.tb

TRUTH = {"P1": {"GO:a"}, "P2": {"GO:b"}}
IC = {"GO:a": 1.0, "GO:b": 1.0, "GO:z": 1.0}


def _panel(pred):
    return rs.build_panel(pred, TRUTH, IC, tb._candidate_thresholds(pred))


def test_all_wrong_panel_has_no_operating_point():
    panel = _panel({"P1": {"GO:z": 0.4}, "P2": {"GO:z": 0.9}})
    assert pd_.f_max_index(panel) is None
    sens, n_pred = pd_.reference_operating_point(panel)
    assert math.isnan(sens) and math.isnan(n_pred)
    assert pd_.pr_at_fmax(panel)["coverage_at_fmax"] == 0.0


def test_a_hit_gives_the_panel_metrics_cutoff():
    panel = _panel({"P1": {"GO:a": 0.9}, "P2": {"GO:z": 0.5}})
    best = pd_.f_max_index(panel)
    assert best is not None
    assert panel.thresholds[best] == rs.panel_metrics(panel)["f_max_tau"]
