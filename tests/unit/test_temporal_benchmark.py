"""Unit tests for the temporal CAFA-style benchmark metrics (VALIDATION_PLAN §2).

The metric maths is exercised on tiny synthetic fixtures — a 3-level ontology and
a handful of proteins — so the CAFA conventions (precision averaged over
predicted proteins, recall over all benchmark proteins, IC-weighted S_min,
max-propagated prediction scores) are pinned exactly.
"""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

# validation/ is not a package; load temporal_benchmark.py by path (mirrors the
# §1 test bootstrap in test_validation_metrics.py).
_TB = Path(__file__).resolve().parents[2] / "validation" / "temporal_benchmark.py"
_spec = importlib.util.spec_from_file_location("temporal_benchmark", _TB)
tb = importlib.util.module_from_spec(_spec)
sys.modules["temporal_benchmark"] = tb
_spec.loader.exec_module(tb)


# A tiny biological-process ontology: leaf -> mid -> BP_ROOT
ANCESTORS = {
    "GO:leaf": {"GO:mid", tb.BP_ROOT},
    "GO:mid": {tb.BP_ROOT},
    tb.BP_ROOT: set(),
    "GO:other": {tb.BP_ROOT},
}


def anc(go):
    return ANCESTORS.get(go, set())


class TestPropagateTerms:
    def test_expands_to_ancestor_closure(self):
        assert tb.propagate_terms({"GO:leaf"}, anc) == {"GO:leaf", "GO:mid", tb.BP_ROOT}

    def test_empty(self):
        assert tb.propagate_terms(set(), anc) == set()


# All the BP-ontology terms map to BP; add an MF term for the aspect-split test.
TERM_ASPECT = {
    "GO:leaf": "BP",
    "GO:mid": "BP",
    tb.BP_ROOT: "BP",
    "GO:other": "BP",
    "GO:mf1": "MF",
    tb.MF_ROOT: "MF",
}


class TestBuildNKBenchmark:
    def test_truth_is_full_t1_not_delta(self):
        # t0 experimental has nothing in BP; t1 adds leaf. Truth is the FULL
        # propagated t1 BP set (leaf + mid), NOT just the delta — a correct
        # prediction of an already-general term must count, not be penalised.
        t0_exp = {"P1": set()}
        t1_exp = {"P1": {"GO:leaf"}}
        bench = tb.build_nk_benchmark_by_aspect(
            t0_exp, t1_exp, TERM_ASPECT, anc, predictable_proteins={"P1"}
        )
        assert bench["BP"] == {"P1": {"GO:leaf", "GO:mid"}}  # BP_ROOT excluded

    def test_no_knowledge_gate_excludes_prior_experimental(self):
        # P1 already had experimental BP annotation at t0 -> excluded from BP.
        t0_exp = {"P1": {"GO:mid"}}
        t1_exp = {"P1": {"GO:leaf"}}
        bench = tb.build_nk_benchmark_by_aspect(
            t0_exp, t1_exp, TERM_ASPECT, anc, predictable_proteins={"P1"}
        )
        assert "P1" not in bench["BP"]

    def test_roots_are_dropped(self):
        t0_exp = {"P1": set()}
        t1_exp = {"P1": {"GO:mid"}}
        bench = tb.build_nk_benchmark_by_aspect(
            t0_exp, t1_exp, TERM_ASPECT, anc, predictable_proteins={"P1"}
        )
        assert bench["BP"] == {"P1": {"GO:mid"}}

    def test_unpredictable_proteins_excluded(self):
        t0_exp = {"P1": set(), "P2": set()}
        t1_exp = {"P1": {"GO:leaf"}, "P2": {"GO:leaf"}}
        bench = tb.build_nk_benchmark_by_aspect(
            t0_exp, t1_exp, TERM_ASPECT, anc, predictable_proteins={"P1"}
        )
        assert set(bench["BP"]) == {"P1"}

    def test_aspects_are_independent(self):
        # P1 had BP experimental at t0 (blocks BP) but nothing in MF: gains an MF
        # term at t1 -> appears only in the MF benchmark.
        t0_exp = {"P1": {"GO:mid"}}
        t1_exp = {"P1": {"GO:leaf", "GO:mf1"}}
        bench = tb.build_nk_benchmark_by_aspect(
            t0_exp, t1_exp, TERM_ASPECT, anc, predictable_proteins={"P1"}
        )
        assert "P1" not in bench["BP"]
        assert bench["MF"] == {"P1": {"GO:mf1"}}


class TestInformationContent:
    def test_marginal_ic_from_frequencies(self):
        # 4 proteins; leaf appears in 1 (=> mid,root also via propagation in all
        # that have leaf). Build so P(mid)=1/2, P(leaf)=1/4.
        ann = {
            "P1": {"GO:leaf"},  # -> leaf, mid, root
            "P2": {"GO:mid"},  # -> mid, root
            "P3": {"GO:other"},  # -> other, root
            "P4": {"GO:other"},  # -> other, root
        }
        ic = tb.information_content(ann, anc)
        assert ic["GO:leaf"] == pytest.approx(-math.log2(1 / 4))
        assert ic["GO:mid"] == pytest.approx(-math.log2(2 / 4))
        # root present in all 4 => P=1 => IC 0 (carries no information)
        assert ic[tb.BP_ROOT] == pytest.approx(0.0)

    def test_empty(self):
        assert tb.information_content({}, anc) == {}


class TestTransferPredictions:
    def test_max_propagation_over_domains(self):
        # P1 has domains D1 (leaf@80) and D2 (other@50).
        protein_domains = {"P1": ["D1", "D2"]}
        domain_go = {"D1": {"GO:leaf": 80.0}, "D2": {"GO:other": 50.0}}
        pred = tb.transfer_predictions(protein_domains, domain_go, anc)
        # leaf propagates to mid+root@80; other@50 (root excluded).
        assert pred["P1"]["GO:leaf"] == 80.0
        assert pred["P1"]["GO:mid"] == 80.0
        assert pred["P1"]["GO:other"] == 50.0
        assert tb.BP_ROOT not in pred["P1"]

    def test_conflicting_scores_take_max(self):
        protein_domains = {"P1": ["D1", "D2"]}
        # Both domains imply mid (D1 via leaf@30, D2 directly@70) -> max 70.
        domain_go = {"D1": {"GO:leaf": 30.0}, "D2": {"GO:mid": 70.0}}
        pred = tb.transfer_predictions(protein_domains, domain_go, anc)
        assert pred["P1"]["GO:mid"] == 70.0
        assert pred["P1"]["GO:leaf"] == 30.0


class TestTransferPscore:
    def test_sums_scores_and_minmax_normalises_per_protein(self):
        # P1: D1 gives leaf@2 (=> leaf,mid each +2), D2 gives mid@3 (=> mid +3).
        # sums: leaf=2, mid=5. min-max over {2,5}: leaf->0.0, mid->1.0.
        protein_domains = {"P1": ["D1", "D2"]}
        domain_go = {"D1": {"GO:leaf": 2.0}, "D2": {"GO:mid": 3.0}}
        pred = tb.transfer_predictions_pscore(protein_domains, domain_go, anc)
        assert pred["P1"]["GO:mid"] == pytest.approx(1.0)
        assert pred["P1"]["GO:leaf"] == pytest.approx(0.0)

    def test_single_term_normalises_to_one(self):
        protein_domains = {"P1": ["D1"]}
        domain_go = {"D1": {"GO:leaf": 7.0}}
        pred = tb.transfer_predictions_pscore(protein_domains, domain_go, anc)
        # leaf and mid both get sum 7 (constant) -> all map to 1.0
        assert pred["P1"]["GO:leaf"] == pytest.approx(1.0)
        assert pred["P1"]["GO:mid"] == pytest.approx(1.0)


class TestPrecisionRecall:
    def test_cafa_averaging_conventions(self):
        # Two eval proteins. P1 predicts {leaf(90), mid(90)}; true {leaf, mid}.
        # P2 predicts nothing above tau; true {other}.
        pred = {"P1": {"GO:leaf": 90.0, "GO:mid": 90.0}}
        true = {"P1": {"GO:leaf", "GO:mid"}, "P2": {"GO:other"}}
        p, r, m = tb.precision_recall_at_threshold(pred, true, tau=50.0)
        # precision averaged over m=1 predicted protein: 2/2 = 1.0
        assert p == pytest.approx(1.0)
        assert m == 1
        # recall averaged over BOTH benchmark proteins: (1.0 + 0.0)/2 = 0.5
        assert r == pytest.approx(0.5)

    def test_threshold_filters_predictions(self):
        pred = {"P1": {"GO:leaf": 40.0}}
        true = {"P1": {"GO:leaf"}}
        p, r, m = tb.precision_recall_at_threshold(pred, true, tau=50.0)
        assert m == 0 and p == 0.0 and r == 0.0


class TestFmax:
    def test_perfect_prediction_reaches_one(self):
        pred = {"P1": {"GO:leaf": 100.0}}
        true = {"P1": {"GO:leaf"}}
        fmax, tau = tb.f_max(pred, true)
        assert fmax == pytest.approx(1.0)

    def test_partial_prediction(self):
        # Predict leaf correctly but also a wrong term at same score: precision 1/2.
        pred = {"P1": {"GO:leaf": 90.0, "GO:other": 90.0}}
        true = {"P1": {"GO:leaf"}}
        fmax, _ = tb.f_max(pred, true)
        # best F: precision 0.5, recall 1.0 -> F = 2*.5*1/(1.5) = 0.667
        assert fmax == pytest.approx(2 / 3, abs=1e-6)


class TestSmin:
    def test_perfect_prediction_zero_distance(self):
        pred = {"P1": {"GO:leaf": 100.0}}
        true = {"P1": {"GO:leaf"}}
        ic = {"GO:leaf": 2.0}
        smin, _ = tb.s_min(pred, true, ic)
        assert smin == pytest.approx(0.0)

    def test_missed_term_contributes_remaining_uncertainty(self):
        # Predict nothing above any tau -> ru = IC(leaf), mi = 0.
        pred = {"P1": {}}
        true = {"P1": {"GO:leaf"}}
        ic = {"GO:leaf": 3.0}
        smin, _ = tb.s_min(pred, true, ic)
        assert smin == pytest.approx(3.0)


class TestAUPRC:
    def test_constant_precision(self):
        # precision 1.0 across recall 0->1 => area 1.0
        curve = [(0.0, 1.0, 1.0), (1.0, 1.0, 0.0)]
        assert tb.auprc(curve) == pytest.approx(1.0)

    def test_single_point_is_zero(self):
        assert tb.auprc([(0.0, 0.5, 0.5)]) == 0.0


class TestShuffleDeterminism:
    def test_seed_is_deterministic(self):
        dg = {"D1": {"GO:leaf": 1.0}, "D2": {"GO:mid": 2.0}, "D3": {"GO:other": 3.0}}
        a = tb.shuffle_domain_go(dg, seed=42)
        b = tb.shuffle_domain_go(dg, seed=42)
        assert a == b
        # keys preserved, association sets permuted among them
        assert set(a) == set(dg)


class TestFilterByIC:
    def test_drops_low_information_terms_from_truth_and_preds(self):
        ic = {"GO:leaf": 5.0, "GO:mid": 1.0}  # mid is near-universal, low IC
        true = {"P1": {"GO:leaf", "GO:mid"}}
        preds = {"P1": {"GO:leaf": 9.0, "GO:mid": 9.0}}
        assert tb.filter_by_ic(true, ic, min_ic=2.0) == {"P1": {"GO:leaf"}}
        assert tb.filter_by_ic(preds, ic, min_ic=2.0) == {"P1": {"GO:leaf": 9.0}}

    def test_zero_floor_is_noop(self):
        ic = {"GO:leaf": 5.0}
        true = {"P1": {"GO:leaf"}}
        assert tb.filter_by_ic(true, ic, min_ic=0.0) == true

    def test_protein_emptied_by_filter_is_dropped(self):
        ic = {"GO:mid": 1.0}
        true = {"P1": {"GO:mid"}}
        assert tb.filter_by_ic(true, ic, min_ic=2.0) == {}


class TestRestrictToAspect:
    def test_splits_predictions_and_truesets(self):
        term_aspect = {"GO:leaf": "BP", "GO:mf1": "MF"}
        preds = {"P1": {"GO:leaf": 5.0, "GO:mf1": 9.0}}
        true = {"P1": {"GO:leaf", "GO:mf1"}}
        assert tb.restrict_to_aspect(preds, "BP", term_aspect) == {
            "P1": {"GO:leaf": 5.0}
        }
        assert tb.restrict_to_aspect(true, "MF", term_aspect) == {"P1": {"GO:mf1"}}
