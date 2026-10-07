"""``extract_human_interpro.py`` cuts a subset of any protein2ipr by GAF selection.

A temporal benchmark needs the t0 architectures of two protein sets: the
proteins training saw (t0 GAF) and the held-out ones annotated only later (t1
GAF). So the selection is the union of several GAFs, and the source can be an
archived protein2ipr rather than the current one.
"""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

#: P1 only in the t0 GAF, P2 only in the t1 GAF (IEA), P3 in both; P9 in neither.
T0_ROWS = [("P1", "IDA"), ("P3", "IDA")]
T1_ROWS = [("P2", "IEA"), ("P3", "IDA")]
SOURCE_PROTEINS = ["P1", "P2", "P3", "P9"]


def _gaf(path: Path, rows: list[tuple[str, str]]) -> Path:
    lines = ["!gaf-version: 2.2"] + [
        f"UniProtKB\t{protein}\t{protein}\t\tGO:0006811\tPMID:1\t{evidence}\t\tP"
        f"\t{protein}\t\tprotein\ttaxon:9606\t20210101\tGOA"
        for protein, evidence in rows
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    return path


def _protein2ipr(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for protein in SOURCE_PROTEINS:
            handle.write(f"{protein}\tIPR000001\tDomain\tPF00001\t10\t110\n")
    return path


def _extract(cwd: Path, *args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "extract_human_interpro.py"), *map(str, args)],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=60,
    )


def _proteins(path: Path) -> set[str]:
    with gzip.open(path, "rt") as handle:
        return {line.split("\t", 1)[0] for line in handle}


@pytest.mark.parametrize(
    "evidence,expected", [("all", {"P1", "P2", "P3"}), ("manual", {"P1", "P3"})]
)
def test_union_of_gafs_selects_from_an_archived_source(
    tmp_path: Path, evidence: str, expected: set[str]
) -> None:
    source = _protein2ipr(tmp_path / "archive" / "protein2ipr.dat.gz")
    t0 = _gaf(tmp_path / "t0.gaf.gz", T0_ROWS)
    t1 = _gaf(tmp_path / "t1.gaf.gz", T1_ROWS)
    output = tmp_path / "subset" / "protein2ipr_human_t0.dat.gz"

    result = _extract(
        tmp_path,
        "--source",
        source,
        "--gaf",
        t0,
        "--gaf",
        t1,
        "--evidence-filter",
        evidence,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert _proteins(output) == expected
    # The selecting list sits beside the custom output, not over the species'.
    listed = (output.parent / "protein2ipr_human_t0_proteins.txt").read_text().split()
    assert set(listed) == expected
    assert not (tmp_path / "data").exists()
    marker = json.loads(Path(f"{output}.provenance.json").read_text())
    assert marker["selection_sources"] == [str(t0), str(t1)]
    assert marker["interpro_source"] == str(source)
    assert marker["evidence_filter"] == evidence


def test_defaults_are_the_species_files(tmp_path: Path) -> None:
    _protein2ipr(tmp_path / "data/raw/interpro_mappings/protein2ipr.dat.gz")
    _gaf(tmp_path / "data/raw/goa_annotations/goa_toy.gaf.gz", T0_ROWS)

    result = _extract(tmp_path, "--species", "toy")

    assert result.returncode == 0, result.stderr
    assert _proteins(tmp_path / "data/interim/protein2ipr_toy.dat.gz") == {"P1", "P3"}
    assert (tmp_path / "data/interim/toy_proteins.txt").read_text().split() == [
        "P1",
        "P3",
    ]


def test_a_missing_gaf_fails_before_scanning(tmp_path: Path) -> None:
    source = _protein2ipr(tmp_path / "protein2ipr.dat.gz")
    output = tmp_path / "subset.dat.gz"

    result = _extract(
        tmp_path,
        "--source",
        source,
        "--gaf",
        _gaf(tmp_path / "t0.gaf.gz", T0_ROWS),
        "--gaf",
        tmp_path / "absent.gaf.gz",
        "--output",
        output,
    )

    assert result.returncode == 1
    assert "absent.gaf.gz" in result.stderr
    assert not output.exists()


@pytest.mark.parametrize("flag", ["--source", "--gaf"])
def test_a_non_default_selection_needs_an_explicit_output(
    tmp_path: Path, flag: str
) -> None:
    # Without --output it would replace the species' default extract, which
    # every run and evaluator reads by default.
    default = _protein2ipr(tmp_path / "data/interim/protein2ipr_toy.dat.gz")
    before = default.read_bytes()
    value = (
        _protein2ipr(tmp_path / "archive.dat.gz")
        if flag == "--source"
        else _gaf(tmp_path / "t0.gaf.gz", T0_ROWS)
    )

    result = _extract(tmp_path, "--species", "toy", flag, value)

    assert result.returncode == 2
    assert "--output" in result.stderr
    assert default.read_bytes() == before
