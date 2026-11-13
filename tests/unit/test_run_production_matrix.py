"""A finished production cell is reused only if it was made with today's flags.

The matrix driver is resumable; a cell whose manifest records another command
(e.g. a changed ``--paperparity-flags``) must be neither silently kept nor
overwritten.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "run_production_matrix", _ROOT / "scripts" / "run_production_matrix.py"
)
matrix = importlib.util.module_from_spec(_spec)
sys.modules["run_production_matrix"] = matrix
_spec.loader.exec_module(matrix)

CELL = "go_experimental_paperparity"


def _finish(root, command):
    (root / CELL).mkdir(parents=True)
    (root / CELL / "run_manifest_go.json").write_text(
        json.dumps({"status": "completed", "command": command})
    )


def _command(root, flags):
    return [
        sys.executable,
        "run_dcgo_human.py",
        "--ontology",
        "go",
        "--evidence-filter",
        "experimental",
        *flags,
        "--num-cores",
        "8",
        "--output-dir",
        str(root / CELL),
        "--domain-key",
        "interpro",
    ]


@pytest.fixture
def no_runs(monkeypatch):
    calls = []
    monkeypatch.setattr(matrix.subprocess, "run", lambda *a, **k: calls.append(a))
    flags = list(matrix.PAPER_PARITY_FLAGS)
    yield calls
    matrix.PAPER_PARITY_FLAGS[:] = flags


def test_a_cell_made_with_the_same_flags_is_skipped(tmp_path, no_runs):
    _finish(tmp_path, _command(tmp_path, matrix.PAPER_PARITY_FLAGS))
    argv = ["--output-root", str(tmp_path), "--only", CELL, "--num-cores", "2"]
    assert matrix.main(argv) == 0
    assert no_runs == []


def test_a_cell_made_with_other_flags_is_reported_not_reused(tmp_path, no_runs):
    _finish(tmp_path, _command(tmp_path, matrix.PAPER_PARITY_FLAGS))
    argv = [
        "--output-root",
        str(tmp_path),
        "--only",
        CELL,
        "--paperparity-flags",
        "--propagate-annotations --enable-true-path --min-ic 1",
    ]
    assert matrix.main(argv) == 1
    assert no_runs == []
