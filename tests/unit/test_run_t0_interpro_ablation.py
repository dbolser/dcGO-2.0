"""The t0-InterPro ablation runbook re-runs exactly the committed rungs.

Its contract is narrow: every rung command is the one the committed manifest
records, with only the inputs and the output directory changed, and the cell
covers the whole ablation ladder.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


runbook = _load(
    "run_t0_interpro_ablation", _ROOT / "scripts" / "run_t0_interpro_ablation.py"
)
ab = _load("ablation", _ROOT / "validation" / "ablation.py")

SUBSET = Path("data/interim/protein2ipr_human_t0_ipr85.dat.gz")


def test_the_cell_covers_the_ablation_ladder():
    assert set(runbook.RUNGS) == {rung.run_dir for rung in ab.LADDER}


@pytest.mark.parametrize("rung", runbook.RUNGS)
def test_rung_commands_change_only_inputs_and_output(rung):
    manifest = runbook.MANIFESTS / f"{rung}.json"
    recorded = json.loads(manifest.read_text())["command"][1:]

    command = runbook.rung_command(manifest, Path("out") / rung, SUBSET)

    assert command[0] == sys.executable
    assert command[-4:] == ["--gaf", str(runbook.T0_GAF), "--interpro", str(SUBSET)]
    out = command.index("--output-dir")
    assert command[out + 1] == str(Path("out") / rung)
    # Everything else, in order, is what the committed run recorded.
    assert command[1:out] + command[out + 2 : -4] == (
        recorded[: out - 1] + recorded[out + 1 :]
    )


def test_the_archive_is_verified_and_the_subset_extracted(tmp_path):
    steps = dict(
        runbook.build_steps(
            runbook.ARCHIVE, tmp_path / "subset.dat.gz", tmp_path, tmp_path
        )
    )
    assert steps["verify"][-2:] == ["--group", "interpro-85"]
    extract = steps["extract"]
    assert extract[extract.index("--source") + 1] == str(runbook.ARCHIVE)
    assert [extract[i + 1] for i, a in enumerate(extract) if a == "--gaf"] == [
        str(runbook.T0_GAF),
        str(runbook.T1_GAF),
    ]
    assert extract[extract.index("--evidence-filter") + 1] == "all"
    evaluate = steps["eval"]
    assert evaluate[evaluate.index("--interpro") + 1] == str(tmp_path / "subset.dat.gz")


def test_a_subset_already_cut_from_the_same_source_is_reused(tmp_path):
    from src.universe_provenance import write_marker

    source = tmp_path / "protein2ipr.dat.gz"
    subset = tmp_path / "subset.dat.gz"
    subset.write_bytes(b"")
    write_marker(
        subset,
        selection_rule="goa",
        selection_sources=[runbook.T0_GAF, runbook.T1_GAF],
        interpro_source=source,
        n_accessions=0,
        n_matched_lines=0,
        tool="extract_human_interpro.py",
    )

    names = [
        name for name, _ in runbook.build_steps(source, subset, tmp_path, tmp_path)
    ]
    assert "verify" not in names  # not the archive: nothing to verify
    assert "extract" not in names
    other = [
        name
        for name, _ in runbook.build_steps(
            tmp_path / "other.dat.gz", subset, tmp_path, tmp_path
        )
    ]
    assert "extract" in other


def test_dry_run_prints_and_runs_nothing(tmp_path, capsys):
    run_dir = tmp_path / "runs"
    assert runbook.main(["--dry-run", "--run-dir", str(run_dir)]) == 0
    printed = capsys.readouterr().out
    assert "# full" in printed and "validation/ablation.py" in printed
    assert not run_dir.exists()
