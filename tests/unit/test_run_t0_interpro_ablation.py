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
    # Every arm is scored on the committed cell's cohort.
    assert evaluate[evaluate.index("--cohort-interpro") + 1] == str(runbook.COHORT)


#: Stand-ins for the two selecting GAFs' SHA-256s (the real files are not in CI).
GAF_HASHES = ("t0-sha", "t1-sha")


@pytest.fixture(autouse=True)
def _gaf_hashes(monkeypatch):
    monkeypatch.setattr(runbook, "selection_hashes", lambda: GAF_HASHES)


def _cut(
    subset,
    source,
    gafs=(runbook.T0_GAF, runbook.T1_GAF),
    evidence="all",
    hashes=GAF_HASHES,
):
    """A subset on disk whose marker records how it was extracted."""
    from src.universe_provenance import write_marker

    subset.write_bytes(b"")
    write_marker(
        subset,
        selection_rule="goa",
        selection_sources=list(gafs),
        interpro_source=source,
        n_accessions=0,
        n_matched_lines=0,
        tool="extract_human_interpro.py",
        evidence_filter=evidence,
        selection_sha256=hashes,
    )


def _names(source, subset, run_dir):
    return [name for name, _ in runbook.build_steps(source, subset, run_dir, run_dir)]


def test_a_subset_cut_by_this_exact_extraction_is_reused(tmp_path):
    source = tmp_path / "protein2ipr.dat.gz"
    subset = tmp_path / "subset.dat.gz"
    _cut(subset, source)

    names = _names(source, subset, tmp_path)
    assert "verify" not in names  # not the archive: nothing to verify
    assert "extract" not in names
    assert "extract" in _names(tmp_path / "other.dat.gz", subset, tmp_path)


@pytest.mark.parametrize(
    "gafs,evidence,hashes",
    [
        ((runbook.T0_GAF,), "all", GAF_HASHES),  # t1-only no-knowledge missing
        ((runbook.T0_GAF, runbook.T1_GAF), "manual", GAF_HASHES),  # extract default
        ((runbook.T0_GAF, runbook.T1_GAF), None, GAF_HASHES),  # predates the field
        ((runbook.T0_GAF, runbook.T1_GAF), "all", ("t0-sha", "old-t1")),  # new t1
        ((runbook.T0_GAF, runbook.T1_GAF), "all", None),  # predates the hashes
    ],
)
def test_a_subset_selected_differently_is_cut_again(tmp_path, gafs, evidence, hashes):
    source = tmp_path / "protein2ipr.dat.gz"
    subset = tmp_path / "subset.dat.gz"
    _cut(subset, source, gafs, evidence, hashes)

    assert "extract" in _names(source, subset, tmp_path)


def _finish(run_dir, rung, status="completed"):
    (run_dir / rung).mkdir()
    (run_dir / rung / "run_manifest_go.json").write_text(json.dumps({"status": status}))


def test_completed_rungs_are_not_rerun(tmp_path):
    source = tmp_path / "src.dat.gz"
    subset = tmp_path / "subset.dat.gz"
    _cut(subset, source)
    _finish(tmp_path, "single")
    _finish(tmp_path, "supra", "running")  # interrupted: never finalised

    names = _names(source, subset, tmp_path)
    assert "single" not in names
    assert names[names.index("supra") :] == [*runbook.RUNGS[1:], "eval"]


def test_a_new_subset_reruns_completed_rungs(tmp_path):
    _finish(tmp_path, "single")  # trained on whatever subset was there before

    names = _names(tmp_path / "src.dat.gz", tmp_path / "subset.dat.gz", tmp_path)
    assert names == ["extract", *runbook.RUNGS, "eval"]


def test_dry_run_prints_and_runs_nothing(tmp_path, capsys):
    run_dir = tmp_path / "runs"
    assert runbook.main(["--dry-run", "--run-dir", str(run_dir)]) == 0
    printed = capsys.readouterr().out
    assert "# full" in printed and "validation/ablation.py" in printed
    assert not run_dir.exists()


def test_the_subset_sidecar_is_kept_with_the_rung_manifests(tmp_path, monkeypatch):
    # The subset sits in gitignored data/interim; its sidecar (archive, GAFs,
    # evidence filter) must travel with the committed evaluation.
    source = tmp_path / "src.dat.gz"
    subset = tmp_path / "subset.dat.gz"
    _cut(subset, source)
    run_dir = tmp_path / "runs"
    run_dir.mkdir()
    for rung in runbook.RUNGS:
        _finish(run_dir, rung)
    monkeypatch.setattr(
        runbook.subprocess,
        "run",
        lambda *a, **k: runbook.subprocess.CompletedProcess(a, 0),
    )

    eval_dir = tmp_path / "eval"
    argv = ["--source", str(source), "--subset", str(subset)]
    argv += ["--run-dir", str(run_dir), "--eval-dir", str(eval_dir)]
    assert runbook.main(argv) == 0

    copied = eval_dir / "manifests" / "subset.dat.gz.provenance.json"
    assert json.loads(copied.read_text())["evidence_filter"] == "all"
    assert {p.name for p in (eval_dir / "manifests").glob("*.json")} == {
        copied.name,
        *(f"{rung}.json" for rung in runbook.RUNGS),
    }
