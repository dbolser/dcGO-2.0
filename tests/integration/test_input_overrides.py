"""``--gaf``/``--interpro`` train on files other than the species' current ones.

The temporal benchmark trains on an archived GOA release and a subset of an
archived InterPro release. These run the real pipeline from a directory that
holds *no* default inputs, so a run can only succeed by reading the overrides,
and check that the manifest identifies them exactly as it does the defaults.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from run_dcgo_human import main

DATE_GENERATED = "2021-04-21 13:14"


def gaf_line(protein: str, term: str) -> str:
    return (
        f"UniProtKB\t{protein}\t{protein}\t\t{term}\tPMID:1\tIDA\t\tP"
        f"\t{protein}\t\tprotein\ttaxon:9606\t20210101\tGOA"
    )


@pytest.fixture
def archived(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """An archived GAF and protein2ipr subset outside the default layout."""
    monkeypatch.chdir(tmp_path)
    gaf_lines = ["!gaf-version: 2.2", f"!date-generated: {DATE_GENERATED}"]
    ipr_lines = []
    for i in range(40):
        protein = f"P{i:05d}"
        term, domain = (
            ("GO:0006811", "IPR000001") if i < 20 else ("GO:0051179", "IPR000002")
        )
        gaf_lines.append(gaf_line(protein, term))
        ipr_lines.append(f"{protein}\t{domain}\tDomain\tPF00001\t10\t110")

    archive = tmp_path / "archive"
    archive.mkdir()
    gaf = archive / "goa_human.gaf.205.gz"
    gaf.write_bytes(gzip.compress(("\n".join(gaf_lines) + "\n").encode()))
    interpro = archive / "protein2ipr_human_t0.dat.gz"
    interpro.write_bytes(gzip.compress(("\n".join(ipr_lines) + "\n").encode()))
    return {"gaf": gaf, "interpro": interpro}


def run(archived: dict[str, Path], out: Path, *extra: str) -> int:
    return main(
        [
            "--gaf",
            str(archived["gaf"]),
            "--interpro",
            str(archived["interpro"]),
            "--output-dir",
            str(out),
            *extra,
        ]
    )


def test_a_run_reads_and_records_the_overrides(
    archived: dict[str, Path], tmp_path: Path
) -> None:
    assert run(archived, tmp_path / "out") == 0

    rows = (tmp_path / "out" / "domain_go_associations_significant.tsv").read_text()
    assert "IPR000001\tGO:0006811" in rows

    manifest = json.loads((tmp_path / "out" / "run_manifest_go.json").read_text())
    by_role = {record["role"]: record for record in manifest["inputs"]}
    for role, key in (("gaf", "gaf"), ("domain_annotations", "interpro")):
        record = by_role[role]
        assert record["path"] == str(archived[key])
        assert (
            record["sha256"] == hashlib.sha256(archived[key].read_bytes()).hexdigest()
        )
        # Not the current-release URL the species name would have implied.
        assert "source_url" not in record and "derived_from" not in record
    assert by_role["gaf"]["release_metadata"]["date_generated"] == DATE_GENERATED
    assert manifest["status"] == "completed"


@pytest.mark.parametrize("key", ["gaf", "interpro"])
def test_a_missing_override_fails_before_the_run(
    archived: dict[str, Path], tmp_path: Path, key: str
) -> None:
    archived[key] = tmp_path / "absent.gz"

    assert run(archived, tmp_path / "out") == 1
    assert not (tmp_path / "out" / "run_manifest_go.json").exists()
