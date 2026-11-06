"""Tests for pure helpers of validation/projection_benchmark.py."""

import gzip
import importlib.util
import sys
from pathlib import Path

import pytest

pytest.importorskip("sklearn")
_PATH = Path(__file__).resolve().parents[2] / "validation" / "projection_benchmark.py"
_spec = importlib.util.spec_from_file_location("projection_benchmark", _PATH)
pb = importlib.util.module_from_spec(_spec)
sys.modules["projection_benchmark"] = pb
_spec.loader.exec_module(pb)

PARENTS = {"leaf": {"mid", "top"}, "mid": {"top"}, "other": {"top"}}


def anc(t):
    return PARENTS.get(t, set())


def test_leaf_yield_counts_each_correct_chain_once():
    truth = {"P": {"leaf", "mid", "top", "other"}}
    preds = {"P": {"leaf": 0.9, "mid": 0.9, "top": 0.95, "other": 0.4}}
    # At tau 0.5: leaf, mid, top are correct; mid and top are implied by leaf.
    assert pb.leaf_yield(truth, preds, 0.5, anc) == 1.0
    assert pb.leaf_yield(truth, preds, 0.3, anc) == 2.0  # + other
    assert pb.leaf_yield(truth, preds, float("nan"), anc) == 0.0


def test_gaf_rows_parses_symbol_taxon_and_not(tmp_path):
    path = tmp_path / "x.gaf.gz"
    cols = [
        "UniProtKB",
        "P1",
        "GENE1",
        "NOT|enables",
        "GO:1",
        "PMID:1",
        "IDA",
        "",
        "F",
        "",
        "",
        "protein",
        "taxon:9606",
        "20200101",
        "UniProt",
    ]
    with gzip.open(path, "wt") as fh:
        fh.write("!gaf-version: 2.2\n")
        fh.write("\t".join(cols) + "\n")
        cols[3] = "enables"
        cols[12] = "taxon:10090|taxon:9606"
        fh.write("\t".join(cols) + "\n")
    rows = list(pb.gaf_rows(path))
    assert rows[0] == ("P1", "GENE1", "taxon:9606", "IDA", "MF", True)
    assert rows[1][2] == "taxon:10090" and rows[1][5] is False
    assert pb.gaf_accessions(path) == {"P1"}
