"""Unit tests for the ablation ladder's pure parts (VALIDATION_PLAN §4).

The driver's I/O is deliberately untested (it needs a 90-minute pipeline run and
an 8 GB GAF); what is testable — and what a reviewer would actually challenge —
is the selection-stage accounting, because the review's specific objection was
that IC thresholds silently change the evaluation cohort.
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
ab = _load("ablation")


IC = {"GO:generic": 0.5, "GO:mid": 2.5, "GO:rare": 6.0}

T0 = {"P1": {"GO:generic"}, "P2": {"GO:mid"}}
T1 = {"P3": {"GO:generic"}, "P4": {"GO:rare"}, "P5": {"GO:mid"}}

BENCH_ALL = {
    "BP": {"P3": {"GO:generic"}, "P4": {"GO:rare"}, "P5": {"GO:mid"}},
    "MF": {},
    "CC": {},
}
# P5 has no InterPro domain, so it never reaches the scored cohort.
BENCH_DOM = {
    "BP": {"P3": {"GO:generic"}, "P4": {"GO:rare"}},
    "MF": {},
    "CC": {},
}


def _rows(**kwargs):
    return ab.selection_stage_counts(
        T0, T1, ["P1", "P2", "P3", "P4"], BENCH_ALL, BENCH_DOM, IC, [0.0, 2.0, 4.0]
    )


class TestLadder:
    def test_every_rung_after_the_first_declares_its_parent(self):
        names = {r.name for r in ab.LADDER}
        assert ab.LADDER[0].parent is None
        for rung in ab.LADDER[1:]:
            assert rung.parent in names, rung.name

    def test_rung_names_are_unique(self):
        names = [r.name for r in ab.LADDER]
        assert len(names) == len(set(names))

    def test_true_path_rungs_read_the_propagated_file(self, tmp_path):
        by_name = {r.name: r for r in ab.LADDER}
        assert (
            ab.rung_prediction_file(tmp_path, by_name["supra_tpr"]).name
            == "domain_go_annotations_propagated.tsv"
        )
        assert (
            ab.rung_prediction_file(tmp_path, by_name["supra"]).name
            == "domain_go_associations_significant.tsv"
        )

    def test_true_path_rungs_reuse_their_base_run_directory(self):
        by_name = {r.name: r for r in ab.LADDER}
        # supra_tpr is post-processing of the supra run, not a separate one.
        assert by_name["supra_tpr"].run_dir == by_name["supra"].run_dir
        assert by_name["full"].run_dir == by_name["supra_shrink"].run_dir


class TestSelectionStageCounts:
    def test_reports_the_input_stages(self):
        rows = {r["stage"]: r for r in _rows() if r["aspect"] == "-"}
        assert rows["t0_annotated_proteins"]["n_proteins"] == 2
        assert rows["t1_experimental_proteins"]["n_proteins"] == 3
        assert rows["proteins_with_domains"]["n_proteins"] == 4

    def test_shows_the_proteins_lost_to_the_has_a_domain_filter(self):
        rows = _rows()
        candidates = next(
            r
            for r in rows
            if r["stage"] == "no_knowledge_candidates" and r["aspect"] == "BP"
        )
        scored = next(
            r
            for r in rows
            if r["stage"] == "no_knowledge_with_domains" and r["aspect"] == "BP"
        )
        assert candidates["n_proteins"] == 3
        assert scored["n_proteins"] == 2  # P5 has no domains

    def test_ic_floors_change_the_cohort_and_that_is_reported(self):
        """The whole point of this ledger: IC>=0 and IC>=4 are different proteins."""
        rows = [
            r
            for r in _rows()
            if r["stage"] == "cohort_at_ic_floor" and r["aspect"] == "BP"
        ]
        by_floor = {r["min_ic"]: r for r in rows}
        assert by_floor["0"]["n_proteins"] == 2  # P3 (generic) + P4 (rare)
        assert by_floor["2"]["n_proteins"] == 1  # P3's only term drops out
        assert by_floor["4"]["n_proteins"] == 1  # only P4's GO:rare survives
        # ...and the ledger says how much of the IC>=0 cohort each floor loses.
        assert by_floor["0"]["n_dropped_vs_ic0"] == 0
        assert by_floor["4"]["n_dropped_vs_ic0"] == 1
        assert by_floor["4"]["pct_of_ic0"] == pytest.approx(50.0)

    def test_one_row_per_aspect_and_floor(self):
        rows = [r for r in _rows() if r["stage"] == "cohort_at_ic_floor"]
        assert len(rows) == 3 * 3  # BP/MF/CC x three floors

    def test_empty_aspects_are_still_reported_as_zero(self):
        rows = [
            r
            for r in _rows()
            if r["stage"] == "no_knowledge_with_domains" and r["aspect"] == "MF"
        ]
        assert rows[0]["n_proteins"] == 0


class TestAssociationCount:
    def test_excludes_the_header(self, tmp_path):
        path = tmp_path / "assoc.tsv"
        path.write_text("domain\tgo_term\n" + "IPR1\tGO:1\n" * 5)
        assert ab.association_count(path) == 5

    def test_header_only_file_is_zero(self, tmp_path):
        path = tmp_path / "assoc.tsv"
        path.write_text("domain\tgo_term\n")
        assert ab.association_count(path) == 0

    def test_empty_file_is_zero(self, tmp_path):
        path = tmp_path / "assoc.tsv"
        path.write_text("")
        assert ab.association_count(path) == 0
