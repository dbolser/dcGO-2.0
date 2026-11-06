"""Tests for the pooled knowledge-yield criterion (validation/knowledge_yield.py)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_PATH = Path(__file__).resolve().parents[2] / "validation" / "knowledge_yield.py"
_spec = importlib.util.spec_from_file_location("knowledge_yield", _PATH)
ky = importlib.util.module_from_spec(_spec)
sys.modules["knowledge_yield"] = ky
_spec.loader.exec_module(ky)


TRUTH = {"A": {"t1", "t2"}, "B": {"t3"}, "C": {"t4"}}


def test_yield_is_largest_tp_at_or_above_precision():
    preds = {
        "A": {"t1": 0.9, "t2": 0.5, "x": 0.4},
        "B": {"t3": 0.8, "y": 0.7},
        "C": {"z": 0.3},
    }
    r = ky.pool(TRUTH, preds)
    # Ranked: t1 .9 (T), t3 .8 (T), y .7 (F), t2 .5 (T), x .4 (F), z .3 (F)
    # precision: 1, 1, 2/3, 3/4, 3/5, 1/2
    assert ky.yield_at_precision(r, 0.5).yield_tp == 3
    assert ky.yield_at_precision(r, 0.7).yield_tp == 3  # 3/4 >= 0.7
    assert ky.yield_at_precision(r, 0.8).yield_tp == 2
    y = ky.yield_at_precision(r, 0.7)
    assert y.threshold == 0.5 and y.n_predicted == 4 and y.n_proteins == 2


def test_ties_are_one_threshold():
    preds = {"A": {"t1": 1.0, "x": 1.0}, "B": {"y": 1.0}}
    r = ky.pool(TRUTH, preds)
    # One block of three at 1.0 with one TP: precision 1/3 only.
    assert ky.yield_at_precision(r, 0.5).yield_tp == 0
    assert ky.yield_at_precision(r, 0.3).yield_tp == 1


def test_non_cohort_predictions_ignored_and_unreached_level_is_zero():
    preds = {"Z": {"t1": 1.0}, "A": {"x": 0.9}}
    r = ky.pool(TRUTH, preds)
    assert r.scores.size == 1
    y = ky.yield_at_precision(r, 0.5)
    assert y.yield_tp == 0 and np.isnan(y.threshold)


def test_pooled_auprc_perfect_and_recall_denominator():
    perfect = {"A": {"t1": 0.9, "t2": 0.8}, "B": {"t3": 0.7}, "C": {"t4": 0.6}}
    assert ky.pooled_auprc(ky.pool(TRUTH, perfect)) == pytest.approx(1.0)
    half = {"A": {"t1": 0.9, "t2": 0.8}}  # half the truth, all correct
    assert ky.pooled_auprc(ky.pool(TRUTH, half)) == pytest.approx(0.5)


def test_bits_weight_by_ic():
    preds = {"A": {"t1": 0.9, "t2": 0.8}}
    r = ky.pool(TRUTH, preds, ic={"t1": 2.0, "t2": 3.0})
    assert ky.yield_at_precision(r, 0.5).bits == pytest.approx(5.0)


def test_bootstrap_unit_weights_reproduce_point_estimate():
    preds = {"A": {"t1": 0.9, "x": 0.4}, "B": {"t3": 0.8}}
    order = sorted(TRUTH)
    r = ky.pool(TRUTH, preds, proteins=order)
    w = np.ones(len(order))
    assert (
        ky.yield_at_precision(r, 0.5, w).yield_tp
        == ky.yield_at_precision(r, 0.5).yield_tp
    )
    sizes = np.array([len(TRUTH[p]) for p in order], dtype=float)
    reps = ky.bootstrap({"m": r}, sizes, n_replicates=20, seed=1)
    assert reps["m"]["ky50"].shape == (20,)
    assert np.all(reps["m"]["pooled_auprc"] >= 0)


def test_ratio_interval_handles_zero_denominator():
    lo, hi = ky.ratio_interval(np.array([2.0, 4.0, 6.0]), np.array([1.0, 2.0, 0.0]))
    assert lo == pytest.approx(2.0, abs=0.1) and hi == np.inf


def test_ratio_interval_drops_double_zero_replicates():
    lo, hi = ky.ratio_interval(np.array([0.0, 2.0, 4.0]), np.array([0.0, 1.0, 1.0]))
    assert (lo, hi) == (2.0, 4.0)
    assert all(np.isnan(ky.ratio_interval(np.zeros(3), np.zeros(3))))


def test_vectorised_bootstrap_matches_per_replicate_definition():
    rng = np.random.default_rng(3)
    truth = {f"P{i}": {f"t{j}" for j in range(rng.integers(1, 5))} for i in range(40)}
    preds = {
        p: {f"t{j}": float(rng.choice([0.2, 0.5, 0.9, rng.random()])) for j in range(6)}
        for p in truth
    }
    order = sorted(truth)
    r = ky.pool(truth, preds, proteins=order)
    sizes = np.array([len(truth[p]) for p in order], dtype=float)
    reps = ky.bootstrap({"m": r}, sizes, n_replicates=25, seed=7, chunk=4)
    weights = ky.bootstrap_weights(len(order), 25, 7)
    for i, w in enumerate(weights):
        for lv, key in ((0.3, "ky30"), (0.5, "ky50"), (0.7, "ky70")):
            assert reps["m"][key][i] == ky.yield_at_precision(r, lv, w).yield_tp
        assert reps["m"]["pooled_auprc"][i] == pytest.approx(
            ky._weighted_auprc(r, w, float(w @ sizes))
        )


def test_cluster_weights_move_clusters_together():
    clusters = np.array([0, 0, 1, 2, 2, 2])
    w = ky.cluster_bootstrap_weights(clusters, n_replicates=50, seed=3)
    assert w.shape == (50, 6)
    assert np.all(w[:, 0] == w[:, 1]) and np.all(w[:, 3] == w[:, 5])
    assert np.allclose(w.sum(axis=1) % 1, 0)


def test_precision_at_threshold_with_and_without_weights():
    preds = {"A": {"t1": 0.9, "x": 0.8}, "B": {"t3": 0.7}}
    order = sorted(TRUTH)
    r = ky.pool(TRUTH, preds, proteins=order)
    assert ky.precision_at(r, 0.8) == pytest.approx(0.5)
    assert ky.precision_at(r, 0.7) == pytest.approx(2 / 3)
    w = np.array([[1.0, 2.0, 1.0], [0.0, 1.0, 0.0]])  # order A, B, C
    np.testing.assert_allclose(ky.precision_at(r, 0.7, w), [3 / 4, 1.0])
