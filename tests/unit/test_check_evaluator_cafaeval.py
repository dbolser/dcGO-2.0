"""Unit tests for the CAFA-evaluator cross-check (TODO P0, VALIDATION_PLAN §4).

The export must hand cafaeval exactly what our evaluator scores: zero scores
are "not predicted" in CAFA format, scores must survive the text round trip
bit for bit, and the IC floor must keep the terms ``filter_by_ic`` keeps. The
end-to-end test scores a tiny fixture with cafaeval itself.
"""

import importlib.util
import json
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
ck = _load("check_evaluator_cafaeval")


# A five-term BP chain-and-branch plus the other two roots, so cafaeval builds
# all three namespaces. GO:0000003 hangs off GO:0000001 by part_of only.
OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0000001
name: bp a
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000002
name: bp b
namespace: biological_process
is_a: GO:0000001 ! bp a

[Term]
id: GO:0000003
name: bp c
namespace: biological_process
relationship: part_of GO:0000001 ! bp a

[Term]
id: GO:0000004
name: bp d
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component
"""

PARENTS = {
    "GO:0000001": {"GO:0008150"},
    "GO:0000002": {"GO:0000001"},
    "GO:0000003": {"GO:0000001"},
    "GO:0000004": {"GO:0008150"},
}
TERM_ASPECT = {t: "BP" for t in (*PARENTS, "GO:0008150")}
IC = {
    "GO:0008150": 0.1,
    "GO:0000001": 1.0,
    "GO:0000002": 3.0,
    "GO:0000003": 2.5,
    "GO:0000004": 0.5,
}


def get_ancestors(term):
    out, stack = set(), list(PARENTS.get(term, ()))
    while stack:
        parent = stack.pop()
        if parent not in out:
            out.add(parent)
            stack.extend(PARENTS.get(parent, ()))
    return out


# Raw (unpropagated) t1 annotations, and the propagated truth our benchmark
# builds from them (aspect root excluded).
T1 = {"P1": {"GO:0000002"}, "P2": {"GO:0000003", "GO:0000004"}, "P3": {"GO:0000004"}}
TRUTH = {
    p: tb.propagate_terms(terms, get_ancestors) - tb.ASPECT_ROOTS
    for p, terms in T1.items()
}
# Consistent (ancestor >= descendant) per-protein scores, as the p-score makes.
PREDS = {
    "P1": {"GO:0000001": 1.0, "GO:0000002": 0.7, "GO:0000003": 0.3},
    "P2": {"GO:0000001": 0.9, "GO:0000003": 0.9, "GO:0000004": 0.0},
    "P3": {"GO:0000001": 1.0, "GO:0000002": 1 / 3},
}


class TestExport:
    def test_ground_truth_is_sorted_protein_term_pairs(self, tmp_path):
        path = tmp_path / "gt.tsv"
        n = ck.write_ground_truth(path, {"P2": {"GO:2", "GO:1"}, "P1": {"GO:3"}})
        assert n == 3
        assert path.read_text().splitlines() == [
            "P1\tGO:3",
            "P2\tGO:1",
            "P2\tGO:2",
        ]

    def test_zero_scores_are_left_out(self, tmp_path):
        path = tmp_path / "pred.tsv"
        n = ck.write_predictions(path, {"P1": {"GO:1": 0.5, "GO:2": 0.0}})
        assert n == 1
        assert path.read_text() == "P1\tGO:1\t0.5\n"

    def test_scores_round_trip_bit_for_bit(self, tmp_path):
        path = tmp_path / "pred.tsv"
        scores = {"GO:1": 1 / 3, "GO:2": 0.1 + 0.2, "GO:3": np.nextafter(1.0, 0.0)}
        ck.write_predictions(path, {"P1": scores})
        back = {
            term: float(score)
            for _, term, score in (
                line.split("\t") for line in path.read_text().splitlines()
            )
        }
        assert back == scores

    @pytest.mark.parametrize("score", [-0.1, 1.5])
    def test_out_of_range_score_is_refused(self, tmp_path, score):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            ck.write_predictions(tmp_path / "p.tsv", {"P1": {"GO:1": score}})

    def test_information_accretion_file(self, tmp_path):
        path = tmp_path / "ia.tsv"
        assert ck.write_information_accretion(path, {"GO:2": 1.5, "GO:1": 0.25}) == 2
        assert path.read_text() == "GO:1\t0.25\nGO:2\t1.5\n"


class TestRawTruth:
    def test_keeps_cohort_and_aspect_only(self):
        t1 = {"P1": {"GO:bp", "GO:mf"}, "P2": {"GO:bp"}, "P3": {"GO:mf"}}
        aspects = {"GO:bp": "BP", "GO:mf": "MF"}
        assert ck.raw_truth(t1, ["P1", "P3"], "BP", aspects) == {"P1": {"GO:bp"}}


class TestFloorTerms:
    def test_matches_filter_by_ic_and_drops_roots(self):
        terms = set(IC) | {"GO:unseen"}
        assert ck.floor_terms(terms, IC, 0.0) == terms - {"GO:0008150"}
        assert ck.floor_terms(terms, IC, 2.0) == {"GO:0000002", "GO:0000003"}

    def test_unseen_terms_have_ic_zero(self):
        # filter_by_ic gives an IC-less term 0 bits, so any positive floor drops it.
        assert "GO:unseen" not in ck.floor_terms({"GO:unseen"}, IC, 0.5)


class TestCurveSummary:
    def test_reproduces_panel_metrics_on_our_curve(self):
        taus = tb._candidate_thresholds(PREDS)
        panel = rs.build_panel(PREDS, TRUTH, IC, taus)
        ours = rs.panel_metrics(panel)
        got = ck.curve_summary(ck.our_curve(PREDS, TRUTH, IC, taus))
        for key in (
            "f_max",
            "f_max_tau",
            "coverage_at_fmax",
            "s_min",
            "s_min_tau",
            "auprc",
        ):
            assert got[key] == pytest.approx(ours[key], abs=1e-15), key

    def test_curves_on_different_thresholds_are_refused(self):
        a = ck.our_curve(PREDS, TRUTH, IC, [0.5])
        b = ck.our_curve(PREDS, TRUTH, IC, [0.6])
        with pytest.raises(ValueError, match="different thresholds"):
            ck.max_abs_diff(a, b)


class TestManifestHashes:
    def test_reports_only_disagreeing_roles(self, tmp_path):
        gaf = tmp_path / "t0.gaf"
        gaf.write_text("x")
        good = ck.sha256(gaf)
        (tmp_path / "m").mkdir()
        (tmp_path / "m" / "a.json").write_text(
            json.dumps(
                {
                    "inputs": [
                        {"role": "gaf", "sha256": good},
                        {"role": "domain_annotations", "sha256": "ignored"},
                    ]
                }
            )
        )
        (tmp_path / "m" / "b.json").write_text(
            json.dumps({"inputs": [{"role": "gaf", "sha256": "0" * 64}]})
        )
        problems = ck.manifest_hash_mismatches(tmp_path / "m", {"gaf": gaf})
        assert len(problems) == 1 and problems[0].startswith("b.json: gaf")

    def test_a_role_no_manifest_records_is_reported(self, tmp_path):
        gaf = tmp_path / "t0.gaf"
        gaf.write_text("x")
        (tmp_path / "m").mkdir()
        (tmp_path / "m" / "a.json").write_text(
            json.dumps({"inputs": [{"role": "gaf", "sha256": ck.sha256(gaf)}]})
        )
        expected = {"gaf": gaf, "go_obo": gaf}
        problems = ck.manifest_hash_mismatches(tmp_path / "m", expected)
        assert problems == [f"no manifest under {tmp_path / 'm'} records role go_obo"]

    def test_a_wrong_manifest_dir_does_not_pass(self, tmp_path):
        gaf = tmp_path / "t0.gaf"
        gaf.write_text("x")
        problems = ck.manifest_hash_mismatches(tmp_path / "missing", {"gaf": gaf})
        assert problems == [
            f"no manifest under {tmp_path / 'missing'} records role gaf"
        ]


class TestAgainstCafaeval:
    """Score the fixture with cafaeval and require our numbers back."""

    @pytest.fixture
    def job(self, tmp_path, monkeypatch):
        pytest.importorskip("cafaeval")
        obo = tmp_path / "go.obo"
        obo.write_text(OBO)
        ia = tmp_path / "ia.tsv"
        ck.write_information_accretion(ia, IC)
        pred_file = tmp_path / "pred.tsv"
        ck.write_predictions(pred_file, PREDS)
        gt_files, cells = {}, {}
        for min_ic in (0.0, 2.0):
            cell = ck.abl.floor_cell(TRUTH, PREDS, IC, min_ic)
            cells[("m", "BP", min_ic)] = cell
            gt_files[min_ic] = tmp_path / f"gt{min_ic:g}.tsv"
            ck.write_ground_truth(
                gt_files[min_ic], ck.raw_truth(T1, cell[0], "BP", TERM_ASPECT)
            )
        ontologies, gts = ck.parse_cafa(obo, ia, gt_files)
        job = ck.Job(
            ontologies=ontologies,
            gts=gts,
            pred_files={"m": pred_file},
            cells=cells,
            ic=IC,
            grid=np.arange(0.01, 1, 0.01),
            n_threads=1,
        )
        monkeypatch.setattr(ck, "_JOB", job)
        return job

    def test_truth_predictions_and_curves_agree(self, job):
        rows = {
            row["min_ic"]: row
            for row in [*ck._score_cell(("m", 0.0)), *ck._score_cell(("m", 2.0))]
        }
        assert set(rows) == {0.0, 2.0}
        for min_ic, row in rows.items():
            truth_f, pred_f, _ = job.cells[("m", "BP", min_ic)]
            assert row["n_eval_cafaeval"] == len(truth_f)
            assert row["truth_mismatch"] == 0
            assert row["prediction_mismatch"] == 0
            assert row["curve_max_abs_diff"] < 1e-12
            # As reported: the grid, zero-coverage rows dropped.
            grid = ck.our_curve(pred_f, truth_f, IC, job.grid)
            keep = grid.coverage > 0
            ours = ck.curve_summary(
                ck.Curve(
                    *(
                        x[keep]
                        for x in (
                            grid.tau,
                            grid.precision,
                            grid.recall,
                            grid.coverage,
                            grid.s,
                        )
                    )
                )
            )
            for key in (
                "f_max",
                "f_max_tau",
                "coverage_at_fmax",
                "s_min",
                "s_min_tau",
                "auprc",
            ):
                assert row[f"{key}_cafaeval"] == pytest.approx(ours[key], abs=1e-12)

    def test_the_floor_shrinks_the_cohort(self, job):
        rows = {
            row["min_ic"]: row
            for row in [*ck._score_cell(("m", 0.0)), *ck._score_cell(("m", 2.0))]
        }
        # P3's truth is {GO:0000004} (0.5 bits): gone at IC >= 2.
        assert rows[0.0]["n_eval_cafaeval"] == 3
        assert rows[2.0]["n_eval_cafaeval"] == 2
