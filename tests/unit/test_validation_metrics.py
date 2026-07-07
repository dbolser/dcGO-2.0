"""Unit tests for the reframed InterPro2GO validation logic (VALIDATION_PLAN.md §1)."""

import importlib.util
import sys
from pathlib import Path

import pytest

# validation/ is not a package; load validate_results.py by path.
_VALIDATION = Path(__file__).resolve().parents[2] / "validation" / "validate_results.py"
_spec = importlib.util.spec_from_file_location("validate_results", _VALIDATION)
vr = importlib.util.module_from_spec(_spec)
sys.modules["validate_results"] = vr
_spec.loader.exec_module(vr)


# A tiny 3-level ontology: child -> parent -> root
ANCESTORS = {
    "GO:child": {"GO:parent", "GO:root"},
    "GO:parent": {"GO:root"},
    "GO:root": set(),
    "GO:other": set(),
}


def get_ancestors(go):
    return ANCESTORS.get(go, set())


class TestPropagatePairs:
    def test_expands_to_ancestor_closure(self):
        result = vr.propagate_pairs({("D1", "GO:child")}, get_ancestors)
        assert result == {
            ("D1", "GO:child"),
            ("D1", "GO:parent"),
            ("D1", "GO:root"),
        }

    def test_original_term_retained_when_no_ancestors(self):
        assert vr.propagate_pairs({("D1", "GO:root")}, get_ancestors) == {
            ("D1", "GO:root")
        }

    def test_no_propagation_callable(self):
        pairs = {("D1", "GO:child"), ("D2", "GO:parent")}
        assert vr.propagate_pairs(pairs, lambda _g: set()) == pairs

    def test_domains_do_not_cross_contaminate(self):
        result = vr.propagate_pairs(
            {("D1", "GO:child"), ("D2", "GO:other")}, get_ancestors
        )
        assert ("D2", "GO:parent") not in result
        assert ("D1", "GO:parent") in result


class TestRestrictToSharedDomains:
    def test_keeps_only_shared_domains(self):
        pred = {("D1", "GO:a"), ("D2", "GO:b")}
        ref = {("D1", "GO:c"), ("D3", "GO:d")}
        p, r, shared = vr.restrict_to_shared_domains(pred, ref)
        assert shared == {"D1"}
        assert p == {("D1", "GO:a")}
        assert r == {("D1", "GO:c")}

    def test_no_overlap_yields_empty(self):
        p, r, shared = vr.restrict_to_shared_domains({("D1", "GO:a")}, {("D2", "GO:b")})
        assert shared == set() and p == set() and r == set()


class TestComputeMetrics:
    def test_propagation_credits_specific_prediction_against_general_reference(self):
        # Prediction is the specific child; reference curates the parent.
        pred = vr.propagate_pairs({("D1", "GO:child")}, get_ancestors)
        ref = vr.propagate_pairs({("D1", "GO:parent")}, get_ancestors)
        pred, ref, _ = vr.restrict_to_shared_domains(pred, ref)
        m = vr.compute_metrics(pred, ref, "t")
        # (D1, parent) and (D1, root) are recovered against the propagated reference.
        assert m["recovered"] == 2
        assert m["reference_coverage"] == pytest.approx(1.0)

    def test_candidates_are_not_false_positives_naming(self):
        pred = {("D1", "GO:a"), ("D1", "GO:b")}
        ref = {("D1", "GO:a")}
        m = vr.compute_metrics(pred, ref, "t")
        assert m["candidate_predictions"] == 1
        assert "false_positives" not in m
        assert m["reference_coverage"] == pytest.approx(1.0)
        assert m["precision_lower_bound"] == pytest.approx(0.5)

    def test_empty_sets_are_safe(self):
        m = vr.compute_metrics(set(), set(), "t")
        assert m["reference_coverage"] == 0.0
        assert m["precision_lower_bound"] == 0.0


class TestThresholdSweepDedup:
    def test_sweep_deduplicates_thresholds(self):
        import pandas as pd

        preds = pd.DataFrame(
            {
                "domain": ["D1", "D1"],
                "go_term": ["GO:child", "GO:other"],
                "p_value": [1e-20, 1e-3],
                "adj_p_value": [1e-15, 1e-2],
                "hyper_score": [100.0, 50.0],
            }
        )
        reference = {("D1", "GO:parent")}
        # 1e-2 appears twice; the sweep must not emit duplicate rows.
        thresholds = {
            "p_value": [1e-2, 1e-2],
            "adj_p_value": [1e-2, 0.01],
            "hyper_score": [50, 50],
        }
        df = vr.calculate_overlap_at_thresholds(
            preds, reference, thresholds, get_ancestors
        )
        assert df["threshold"].is_unique
        # one p_value row, one adj_p row, one score row
        assert len(df) == 3
