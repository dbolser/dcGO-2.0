"""The specificity summary behind VALIDATION_PLAN's chain/ancestor numbers.

Those numbers (28.6% / 55.2% / 82.4% "on a chain") were originally computed ad
hoc; ``validation/specificity_metrics.py`` makes them reproducible. The fixture
here is a three-level chain — every expected value is computable by hand.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from validation.specificity_metrics import (
    SpecificityMetrics,
    load_association_rows,
    main,
    specificity_metrics,
)

ANCESTORS = {
    "LEAF": {"MID", "ROOT"},
    "MID": {"ROOT"},
    "ROOT": set(),
}


def get_ancestors(term: str) -> set[str]:
    return ANCESTORS.get(term, set())


class TestSpecificityMetrics:
    def test_hand_computable_chain(self) -> None:
        """Domain X reports a full chain; domain Y only the leaf."""
        pairs = [
            ("X", "LEAF"),
            ("X", "MID"),
            ("X", "ROOT"),
            ("Y", "LEAF"),
        ]
        m = specificity_metrics(pairs, get_ancestors, roots=frozenset({"ROOT"}))
        assert m == SpecificityMetrics(
            n_associations=4,
            # (2 + 1 + 0 + 2) / 4
            mean_ancestors=1.25,
            # X-LEAF and X-MID have an ancestor in X's set; X-ROOT (nothing
            # above it) and Y-LEAF (Y has no ancestor rows) do not.
            on_chain_share=0.5,
            roots_present=("ROOT",),
        )

    def test_chains_never_cross_domains(self) -> None:
        """Y holding MID does not put X's LEAF on a chain."""
        m = specificity_metrics(
            [("X", "LEAF"), ("Y", "MID")], get_ancestors, roots=frozenset({"ROOT"})
        )
        assert m.on_chain_share == 0.0
        assert m.roots_present == ()

    def test_unknown_terms_count_but_cannot_chain(self) -> None:
        m = specificity_metrics(
            [("X", "NOT-IN-ONTOLOGY")], get_ancestors, roots=frozenset({"ROOT"})
        )
        assert m.n_associations == 1
        assert m.mean_ancestors == 0.0
        assert m.on_chain_share == 0.0

    def test_empty_input_is_all_zero(self) -> None:
        m = specificity_metrics([], get_ancestors)
        assert m.n_associations == 0
        assert m.mean_ancestors == 0.0
        assert m.on_chain_share == 0.0
        assert m.roots_present == ()


HEADER = "domain\tgo_term\tdomain_type\tic\n"

OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0006811
name: ion transport
namespace: biological_process
is_a: GO:0008150 ! biological_process
"""


@pytest.fixture
def obo(tmp_path: Path) -> Path:
    path = tmp_path / "mini.obo"
    path.write_text(OBO, encoding="utf-8")
    return path


def write_table(tmp_path: Path, *lines: str) -> Path:
    path = tmp_path / "assoc.tsv"
    path.write_text(HEADER + "".join(f"{line}\n" for line in lines))
    return path


class TestLoadAssociationRows:
    def test_header_only_table_is_valid_and_empty(self, tmp_path: Path) -> None:
        rows, fieldnames = load_association_rows(write_table(tmp_path))
        assert rows == []
        assert fieldnames == ["domain", "go_term", "domain_type", "ic"]

    def test_missing_term_column_fails_even_on_an_empty_table(
        self, tmp_path: Path
    ) -> None:
        with pytest.raises(SystemExit, match="ec_term"):
            load_association_rows(write_table(tmp_path), term_column="ec_term")


class TestCliGuards:
    def test_empty_but_valid_table_sweeps_to_zero_count_rows(
        self, tmp_path: Path, obo: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """The ic column lives in the header, so a zero-association run is a
        legitimate all-zero sweep, not a schema error."""
        table = write_table(tmp_path)
        assert (
            main(["--associations", str(table), "--obo", str(obo), "--min-ic", "1"])
            == 0
        )
        out = capsys.readouterr().out.splitlines()
        assert out[1].startswith("1\t0\t")

    def test_wrong_ontology_is_refused_not_reported_as_perfect(
        self, tmp_path: Path, obo: Path
    ) -> None:
        """Zero overlap yields empty closures everywhere — a 0.0%-chain result
        indistinguishable from a perfect one. Refuse it."""
        table = write_table(tmp_path, "IPR000001\tHP:0000001\tsingle\t2.0")
        with pytest.raises(SystemExit, match="wrong ontology"):
            main(["--associations", str(table), "--obo", str(obo)])

    def test_partial_ontology_mismatch_warns(
        self, tmp_path: Path, obo: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        table = write_table(
            tmp_path,
            "IPR000001\tGO:0006811\tsingle\t2.0",
            "IPR000001\tHP:0000001\tsingle\t2.0",
            "IPR000001\tHP:0000002\tsingle\t2.0",
        )
        assert main(["--associations", str(table), "--obo", str(obo)]) == 0
        assert "WARNING" in capsys.readouterr().err

    def test_unrecognised_domain_type_value_is_an_error(
        self, tmp_path: Path, obo: Path
    ) -> None:
        table = write_table(tmp_path, "IPR000001\tGO:0006811\tsupra\t2.0")
        with pytest.raises(SystemExit, match="unrecognised domain_type"):
            main(
                [
                    "--associations",
                    str(table),
                    "--obo",
                    str(obo),
                    "--domain-type",
                    "supra",
                ]
            )

    def test_missing_domain_type_column_is_an_error(
        self, tmp_path: Path, obo: Path
    ) -> None:
        path = tmp_path / "assoc.tsv"
        path.write_text("domain\tgo_term\nIPR000001\tGO:0006811\n")
        with pytest.raises(SystemExit, match="domain_type"):
            main(
                [
                    "--associations",
                    str(path),
                    "--obo",
                    str(obo),
                    "--domain-type",
                    "single",
                ]
            )

    def test_supra_selects_both_supra_families(
        self, tmp_path: Path, obo: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        table = write_table(
            tmp_path,
            "IPR000001\tGO:0006811\tsingle\t2.0",
            "IPR000001,IPR000002\tGO:0006811\tsupra_pair\t2.0",
            "IPR000001,IPR000002,IPR000003\tGO:0006811\tsupra_triple\t2.0",
        )
        assert (
            main(
                [
                    "--associations",
                    str(table),
                    "--obo",
                    str(obo),
                    "--domain-type",
                    "supra",
                ]
            )
            == 0
        )
        out = capsys.readouterr().out.splitlines()
        assert out[1].split("\t")[1] == "2"
