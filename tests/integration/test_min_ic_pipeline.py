"""End-to-end tests for the --min-ic reporting floor and the exported ic column.

The fixture reproduces, in miniature, the defect the floor exists for: a
curator annotates a protein set directly to a DAG root, so the root looks
Fisher-enriched for their domain even though, once annotations are True-Path
propagated, every protein in the universe carries it (P = 1, IC = 0). The
relative inference can never remove that association — a root has no parents to
test against — so the floor is the mechanism that does.

Design (hand-computable):
  * 40 proteins. P0000–P0019 carry IPR000001 and are annotated to the root
    GO:0008150 *directly* plus the specific GO:0006811. P0020–P0039 carry
    IPR000002 and are annotated to GO:0051179 (another child of the root).
  * Propagated frequencies: root 40/40 → IC 0; each child 20/40 → IC 1.
  * All three co-occurring (domain, term) pairs have the maximally enriched
    20/0/0/20 table, p = 1/C(40,20) ≈ 7.3e-12 — all significant under BH.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from run_dcgo_human import main

ROOT = "GO:0008150"
SPECIFIC = "GO:0006811"
OTHER = "GO:0051179"

OBO = """format-version: 1.2
data-version: test-min-ic

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0006811
name: ion transport
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0051179
name: localization
namespace: biological_process
is_a: GO:0008150 ! biological_process
"""


def gaf_line(protein: str, term: str) -> str:
    return (
        f"UniProtKB\t{protein}\t{protein}\t\t{term}\tPMID:1\tIDA\t\tP"
        f"\t{protein}\t\tprotein\ttaxon:9606\t20260101\tGOA"
    )


@pytest.fixture
def pipeline_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A working directory holding the miniature species 'ictest'."""
    monkeypatch.chdir(tmp_path)

    gaf_lines = ["!gaf-version: 2.2"]
    ipr_lines = []
    for i in range(40):
        protein = f"P{i:05d}"
        if i < 20:
            gaf_lines += [gaf_line(protein, ROOT), gaf_line(protein, SPECIFIC)]
            ipr_lines.append(f"{protein}\tIPR000001\tDomain A\tPF00001\t10\t110")
        else:
            gaf_lines.append(gaf_line(protein, OTHER))
            ipr_lines.append(f"{protein}\tIPR000002\tDomain B\tPF00002\t10\t110")

    gaf = tmp_path / "data/raw/goa_annotations/goa_ictest.gaf.gz"
    gaf.parent.mkdir(parents=True)
    gaf.write_bytes(gzip.compress(("\n".join(gaf_lines) + "\n").encode()))

    ipr = tmp_path / "data/interim/protein2ipr_ictest.dat.gz"
    ipr.parent.mkdir(parents=True)
    ipr.write_bytes(gzip.compress(("\n".join(ipr_lines) + "\n").encode()))

    (tmp_path / "go.obo").write_text(OBO, encoding="utf-8")
    return tmp_path


def run_pipeline(pipeline_dir: Path, out: str, *extra: str) -> tuple[list[dict], dict]:
    """Run main() and return (significant rows as dicts, manifest json)."""
    code = main(
        [
            "--species",
            "ictest",
            "--go-ontology",
            str(pipeline_dir / "go.obo"),
            "--output-dir",
            str(pipeline_dir / out),
            *extra,
        ]
    )
    assert code == 0
    tsv = (pipeline_dir / out / "domain_go_associations_significant.tsv").read_text()
    header, *lines = tsv.splitlines()
    columns = header.split("\t")
    rows = [dict(zip(columns, line.split("\t"))) for line in lines]
    manifest = json.loads(
        (pipeline_dir / out / "run_manifest_go.json").read_text(encoding="utf-8")
    )
    return rows, manifest


class TestIcColumn:
    def test_ic_is_exported_from_propagated_frequencies(self, pipeline_dir: Path):
        """With a hierarchy engaged, the root's directly-annotated rows show IC 0."""
        rows, manifest = run_pipeline(pipeline_dir, "results", "--enable-true-path")

        by_pair = {(r["domain"], r["go_term"]): r for r in rows}
        assert set(by_pair) == {
            ("IPR000001", ROOT),
            ("IPR000001", SPECIFIC),
            ("IPR000002", OTHER),
        }
        # Propagated frequencies: root 40/40, children 20/40.
        assert float(by_pair[("IPR000001", ROOT)]["ic"]) == 0.0
        assert float(by_pair[("IPR000001", SPECIFIC)]["ic"]) == 1.0
        assert float(by_pair[("IPR000002", OTHER)]["ic"]) == 1.0

        thresholds = manifest["analysis"]["thresholds"]
        assert thresholds["min_ic"] is None
        assert thresholds["ic_source"] == "propagated"

        # The top100 and propagated-annotations exports carry the column too.
        top = pipeline_dir / "results" / "domain_go_associations_top100.tsv"
        assert "\tic" in top.read_text().splitlines()[0]
        propagated = (
            pipeline_dir / "results" / "domain_go_annotations_propagated.tsv"
        ).read_text()
        prop_header, *prop_lines = propagated.splitlines()
        ic_idx = prop_header.split("\t").index("ic")
        prop_ic = {
            (f[0], f[1]): float(f[ic_idx])
            for f in (line.split("\t") for line in prop_lines)
        }
        # The ancestor row rolled up from the specific term reports the
        # *ancestor's* IC — 0 at the root — not its source's.
        assert prop_ic[("IPR000002", ROOT)] == 0.0
        assert prop_ic[("IPR000002", OTHER)] == 1.0

    def test_bare_run_falls_back_to_direct_frequencies(self, pipeline_dir: Path):
        """No hierarchy flags, no floor: the column is the direct-map estimate,
        and the manifest says so (the OBO is not an input of a bare run)."""
        rows, manifest = run_pipeline(pipeline_dir, "results_bare")

        by_pair = {(r["domain"], r["go_term"]): r for r in rows}
        # Direct frequency of the root is 20/40 → IC 1: exactly the mid-level
        # inflation the propagated estimate exists to prevent.
        assert float(by_pair[("IPR000001", ROOT)]["ic"]) == 1.0
        assert manifest["analysis"]["thresholds"]["ic_source"] == "direct"


class TestMinIcFloor:
    def test_floor_removes_the_root_and_only_the_root(self, pipeline_dir: Path):
        rows, manifest = run_pipeline(pipeline_dir, "results_floor", "--min-ic", "0.5")

        pairs = {(r["domain"], r["go_term"]) for r in rows}
        assert pairs == {("IPR000001", SPECIFIC), ("IPR000002", OTHER)}
        assert all(float(r["ic"]) >= 0.5 for r in rows)
        thresholds = manifest["analysis"]["thresholds"]
        assert thresholds["min_ic"] == 0.5
        # The floor alone engages the hierarchy: IC must be propagated.
        assert thresholds["ic_source"] == "propagated"
        assert manifest["summary"]["significant_associations"] == 2

    def test_floor_is_reporting_only_q_values_are_unchanged(self, pipeline_dir: Path):
        """Post-BH placement: surviving rows keep the exact q-values of an
        unfloored run, so the floor never alters the hypothesis family."""
        unfloored, _ = run_pipeline(pipeline_dir, "results_all", "--enable-true-path")
        floored, _ = run_pipeline(pipeline_dir, "results_ic", "--min-ic", "0.5")

        unfloored_q = {(r["domain"], r["go_term"]): r["adj_p_value"] for r in unfloored}
        floored_q = {(r["domain"], r["go_term"]): r["adj_p_value"] for r in floored}
        # The floored run is exactly the unfloored run minus low-IC rows.
        assert set(floored_q) == {
            pair for pair in unfloored_q if pair != ("IPR000001", ROOT)
        }
        for pair, q in floored_q.items():
            assert q == unfloored_q[pair]
