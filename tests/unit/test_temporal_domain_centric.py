"""Unit tests for the association-level temporal metric (ablation add-on)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
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


tdc = _load("temporal_domain_centric")
prd = _load("ablation_pr_decomposition")
rs = prd.rs
ts = _load("temporal_surprise")

# GO: leaf -> mid -> root ; other -> root
ANCESTORS = {"leaf": {"mid", "root"}, "mid": {"root"}, "other": {"root"}, "root": set()}


def get_ancestors(term):
    return ANCESTORS[term]


# Three domains. t0 curation: D1 has mid (so root too); D2 has root; D3 has root.
# t1 adds leaf to D1 (propagated: only (D1, leaf) is new), other to D2
# (new: (D2, other)); D3 is untouched; D4 is a new entry (not shared).
T0 = {("D1", "mid"), ("D2", "root"), ("D3", "root")}
T1 = T0 | {("D1", "leaf"), ("D2", "other"), ("D4", "leaf")}


@pytest.fixture
def delta():
    return tdc.curated_delta(T0, T1, get_ancestors)


def test_added_pairs_are_propagated_and_restricted_to_shared_domains(delta):
    assert delta.shared_domains == {"D1", "D2", "D3"}
    assert delta.added == {("D1", "leaf"), ("D2", "other")}


def _rung(name, rows):
    return tdc.RungTable(name, tuple(rows))


def test_precision_recall_and_budget_on_a_known_example(delta):
    # good: predicts the leaf on D1 (hit) and D3's root (D3 gained nothing, so
    # it is not counted against the rung); one wrong guess on D2.
    good = _rung(
        "good",
        [("D1", "leaf", 1e-9), ("D3", "root", 1e-5), ("D2", "mid", 1e-3)],
    )
    # weak: only re-states t0 knowledge on D1.
    weak = _rung("weak", [("D1", "mid", 1e-4)])
    rows = tdc.evaluate_rungs(
        [good, weak], delta, get_ancestors, budget=1, n_bootstrap=50, seed=1
    )
    by = {(r["rung"], r["metric"]): r for r in rows}

    # good, all: propagated pred on gained domains = (D1,leaf),(D1,mid),(D1,root),
    # (D2,mid),(D2,root) = 5; hits = 1; known-at-t0 among them = (D1,mid),
    # (D1,root),(D2,root) = 3.
    assert by[("good", "precision")]["value"] == pytest.approx(1 / 5)
    assert by[("good", "precision_novel")]["value"] == pytest.approx(1 / 2)
    assert by[("good", "recall")]["value"] == pytest.approx(1 / 2)
    # budget 1 keeps only (D1, leaf): pred = leaf, mid, root -> 1 hit of 3.
    assert by[("good", "precision_at_budget")]["value"] == pytest.approx(1 / 3)
    assert by[("good", "recall_at_budget")]["value"] == pytest.approx(1 / 2)
    assert by[("good", "n_hit")]["value"] == 1
    assert by[("good", "n_direct_eligible")]["value"] == 3

    # weak: (D1,mid),(D1,root) both known at t0 -> precision 0, novel NaN.
    assert by[("weak", "precision")]["value"] == 0.0
    assert np.isnan(by[("weak", "precision_novel")]["value"])
    assert by[("weak", "recall")]["value"] == 0.0

    # Deltas are against the first rung by default and paired per replicate.
    assert by[("weak", "recall")]["reference"] == "good"
    assert by[("weak", "recall")]["delta"] == pytest.approx(-0.5)
    assert by[("good", "recall")]["delta"] == 0.0
    assert (
        by[("good", "recall")]["delta_ci_lo"]
        == 0.0
        == by[("good", "recall")]["delta_ci_hi"]
    )
    assert by[("good", "recall")]["ci_lo"] <= 0.5 <= by[("good", "recall")]["ci_hi"]
    assert by[("good", "recall")]["budget"] == 1


def test_default_budget_is_the_smallest_rungs_eligible_count(delta):
    big = _rung("big", [("D1", "leaf", 1e-9), ("D2", "mid", 1e-3), ("D9", "x", 0.1)])
    small = _rung("small", [("D1", "leaf", 1e-9), ("D2", "mid", 1e-3)])
    rows = tdc.evaluate_rungs(
        [big, small], delta, get_ancestors, universe={"D1", "D2"}, n_bootstrap=5
    )
    assert {r["budget"] for r in rows} == {2}
    # D3 is outside the universe, D9 outside InterPro2GO: neither is eligible.
    assert next(r for r in rows if r["metric"] == "n_eligible_domains")["value"] == 2


def test_bootstrap_is_seeded(delta):
    rung = _rung("r", [("D1", "leaf", 1e-9), ("D2", "mid", 1e-3)])
    a = tdc.evaluate_rungs([rung], delta, get_ancestors, n_bootstrap=30, seed=3)
    b = tdc.evaluate_rungs([rung], delta, get_ancestors, n_bootstrap=30, seed=3)
    assert [r["ci_lo"] for r in a if r["metric"] == "precision"] == [
        r["ci_lo"] for r in b if r["metric"] == "precision"
    ]


def test_signature_mapping_pools_and_counts_unmapped():
    sig_map = {"SSF1": "D1", "SSF2": "D1", "SSF3": "D2"}
    rows, n_unmapped = tdc.map_signatures(
        [
            ("SSF1", "leaf", 1e-3),
            ("SSF2", "leaf", 1e-9),
            ("SSF3", "mid", 1e-2),
            ("SSF9", "mid", 1e-2),
        ],
        sig_map,
    )
    assert n_unmapped == 1
    assert sorted(rows) == [("D1", "leaf", 1e-9), ("D2", "mid", 1e-2)]


def test_protein2ipr_universe_maps_signatures_to_entries(tmp_path):
    lines = [
        "P1\tIPR1\tname\tSSF1\t1\t50",
        "P1\tIPR1\tname\tPF1\t1\t50",
        "P2\tIPR2\tname\tSSF2\t1\t50",
        "P3\tIPR3\tname\tPF2\t1\t50",
        "P4\tIPR2\tname\tSSF3\t1\t50",
        "P4\tIPR3\tname\tSSF3\t1\t50",  # SSF3 under two entries: dropped
    ]
    path = tmp_path / "protein2ipr.dat"
    path.write_text("\n".join(lines) + "\n")
    entries, sig_map = tdc.read_protein2ipr_universe(path, "interpro")
    assert entries == {"IPR1", "IPR2", "IPR3"} and sig_map == {}
    entries, sig_map = tdc.read_protein2ipr_universe(path, "ssf")
    assert entries == {"IPR1", "IPR2", "IPR3"}
    assert sig_map == {"SSF1": "IPR1", "SSF2": "IPR2"}
    entries, sig_map = tdc.read_protein2ipr_universe(path, "pfam")
    assert entries == {"IPR1", "IPR3"} and sig_map == {"PF1": "IPR1", "PF2": "IPR3"}


def test_read_interpro2go_and_single_rows(tmp_path):
    i2g = tmp_path / "interpro2go"
    i2g.write_text(
        "!date: x\n"
        "InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677\n"
        "InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:a ; b ; GO:0003707\n"
    )
    assert tdc.read_interpro2go(i2g) == {
        ("IPR000003", "GO:0003677"),
        ("IPR000003", "GO:0003707"),
    }
    table = tmp_path / "assoc.tsv"
    table.write_text(
        "domain\tgo_term\tq_value\nIPR1\tGO:1\t0.001\nIPR1,IPR2\tGO:1\t0.002\n"
    )
    assert tdc.load_single_associations(table) == [("IPR1", "GO:1", 0.001)]


# --------------------------------------------------------------------------- #
# ablation_pr_decomposition fixed-operating-point columns                      #
# --------------------------------------------------------------------------- #
def _panel(preds, truth):
    taus = sorted({s for p in preds.values() for s in p.values()})
    return rs.build_panel(preds, truth, {}, taus)


def test_ppv_at_sensitivity_and_f_at_budget_on_synthetic_curves():
    truth = {"P1": {"a", "b"}, "P2": {"a"}}
    # Reference: at tau=2 predicts a for both (PPV 1, sens 0.75, 2 preds);
    # at tau=1 adds a wrong term to P2 (PPV 0.75, sens 0.75, 3 preds).
    ref = _panel({"P1": {"a": 2.0}, "P2": {"a": 2.0, "x": 1.0}}, truth)
    sens, n_pred = prd.reference_operating_point(ref)
    assert sens == pytest.approx(0.75) and n_pred == 2.0

    # Rung: tau=3 -> P1:a only (PPV 1, sens 0.25, 1 pred);
    #       tau=1 -> P1:a,b and P2:a,x (PPV 0.75, sens 1.0, 4 preds).
    rung = _panel({"P1": {"a": 3.0, "b": 1.0}, "P2": {"a": 1.0, "x": 1.0}}, truth)
    # Interpolated between (0.25, 1.0) and (1.0, 0.75) at 0.75.
    assert prd.ppv_at_sensitivity(rung, 0.75) == pytest.approx(1 - 0.25 * (0.5 / 0.75))
    curve = rs.panel_curve(rung)
    f1, f4 = float(curve.f[1]), float(curve.f[0])  # thresholds ascend: 1.0, 3.0
    assert prd.f_at_budget(rung, 2.0) == pytest.approx(f1 + (f4 - f1) * (1 / 3))
    # Outside the rung's reach -> NaN, never extrapolated.
    assert np.isnan(prd.ppv_at_sensitivity(ref, 0.9))
    assert np.isnan(prd.f_at_budget(ref, 10.0))


# --------------------------------------------------------------------------- #
# temporal_surprise.emergence_enrichment (the per-rung wrapper)                #
# --------------------------------------------------------------------------- #
def test_emergence_enrichment_pools_all_associations_without_a_ranking():
    # Universe of 4 carriers of A,B plus 4 non-carriers; t0: nobody has "t".
    # By t1 both A,B carriers P1,P2 and one non-carrier gained it.
    features = {f"P{i}": ["A", "B", "A,B"] if i <= 4 else ["C"] for i in range(1, 9)}
    t0 = {p: set() for p in features}
    t1 = {"P1": {"t"}, "P2": {"t"}, "P5": {"t"}}
    result = ts.emergence_enrichment(
        [("A,B", "t"), ("A,B", "t"), ("Z,Z", "t")], features, t0, t1, n_resamples=20
    )
    # Duplicates collapse; Z,Z has no carriers and is dropped.
    assert result.n_associations == 1
    assert result.n_predicted == 4 and result.n_hit == 2
    assert result.expected_rate == pytest.approx(3 / 8)
    assert result.enrichment == pytest.approx(0.5 / (3 / 8))
    assert result.ci_low <= result.enrichment <= result.ci_high
