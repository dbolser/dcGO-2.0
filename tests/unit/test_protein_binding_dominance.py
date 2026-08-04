"""Unit tests for the protein-binding dominance measurement (§2 read-out)."""

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.protein_binding_dominance import (  # noqa: E402
    PROTEIN_BINDING,
    read_protein_terms,
    summarise,
)

EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}


def _gaf_line(protein: str, term: str, evidence: str, aspect: str) -> str:
    """One GAF 2.2 row; only columns 1, 4, 6 and 8 are read."""
    fields = ["UniProtKB", protein, protein, "", term, "PMID:1"]
    fields += [
        evidence,
        "",
        aspect,
        "",
        "",
        "protein",
        "taxon:9606",
        "20200101",
        "UniProt",
    ]
    return "\t".join(fields) + "\n"


@pytest.fixture
def gaf_file(tmp_path: Path) -> Path:
    """Three proteins: one binding-only, one mixed, one with no binding.

    P1 carries `protein binding` twice — two papers, one pair — so the
    line-level and protein-level shares must come out different.
    """
    path = tmp_path / "goa_test.gaf.gz"
    rows = [
        "!gaf-version: 2.2\n",
        _gaf_line("P1", PROTEIN_BINDING, "IPI", "F"),
        _gaf_line("P1", PROTEIN_BINDING, "IDA", "F"),
        _gaf_line("P2", PROTEIN_BINDING, "IPI", "F"),
        _gaf_line("P2", "GO:0004672", "IDA", "F"),
        _gaf_line("P3", "GO:0003700", "IDA", "F"),
        _gaf_line("P4", PROTEIN_BINDING, "IEA", "F"),  # dropped: not experimental
        _gaf_line("P5", PROTEIN_BINDING, "IDA", "P"),  # dropped: wrong aspect
    ]
    with gzip.open(path, "wt") as handle:
        handle.writelines(rows)
    return path


class TestReadProteinTerms:
    def test_filters_evidence_and_aspect(self, gaf_file: Path) -> None:
        protein_terms, n_lines = read_protein_terms(gaf_file, "F", EXPERIMENTAL)
        assert set(protein_terms) == {"P1", "P2", "P3"}
        assert n_lines == 5

    def test_collapses_duplicate_lines_into_one_pair(self, gaf_file: Path) -> None:
        protein_terms, _n = read_protein_terms(gaf_file, "F", EXPERIMENTAL)
        assert protein_terms["P1"] == {PROTEIN_BINDING}


class TestSummarise:
    def test_line_and_protein_measures_differ(self, gaf_file: Path) -> None:
        """The point of the script: the two shares are not the same number."""
        protein_terms, n_lines = read_protein_terms(gaf_file, "F", EXPERIMENTAL)
        rows = dict(
            (name, value)
            for name, value, _share in summarise(
                protein_terms, n_lines, PROTEIN_BINDING, gaf_file
            )
        )
        # 5 experimental F lines, but only 4 distinct pairs (P1 duplicated).
        assert rows["annotation lines"] == "5"
        assert rows["distinct (protein, term) pairs"] == "4"

    def test_counts_carriers_and_exclusive_carriers(self, gaf_file: Path) -> None:
        protein_terms, n_lines = read_protein_terms(gaf_file, "F", EXPERIMENTAL)
        shares = {
            name: share
            for name, _value, share in summarise(
                protein_terms, n_lines, PROTEIN_BINDING, gaf_file
            )
        }
        # P1 and P2 carry it (2/3); only P1 carries nothing else (1/3).
        assert shares[f"proteins carrying {PROTEIN_BINDING}"] == "66.7% of proteins"
        assert (
            shares[f"proteins carrying ONLY {PROTEIN_BINDING}"] == "33.3% of proteins"
        )

    def test_empty_input_returns_no_rows(self, tmp_path: Path) -> None:
        assert summarise({}, 0, PROTEIN_BINDING, tmp_path / "x.gaf") == []
