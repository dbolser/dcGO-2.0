"""Tests for machine-readable run provenance."""

import gzip
import hashlib
import json
import subprocess
from pathlib import Path

from src.run_manifest import (
    RunManifest,
    dependency_lock_metadata,
    describe_file,
    embedded_release_metadata,
    git_metadata,
    json_parameters,
    manifest_filename,
    sha256_file,
)


def test_sha256_and_file_description(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("dcGO\n", encoding="utf-8")

    expected = hashlib.sha256(b"dcGO\n").hexdigest()
    assert sha256_file(path) == expected
    record = describe_file(path, role="test", source_url="https://example.test/x")
    assert record["role"] == "test"
    assert record["size_bytes"] == 5
    assert record["sha256"] == expected
    assert record["source_url"] == "https://example.test/x"


def test_extracts_gaf_release_metadata_from_gzip(tmp_path):
    path = tmp_path / "goa.gaf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("!gaf-version: 2.2\n")
        handle.write("!date-generated: 2026-07-01\n")
        handle.write("!generated-by: UniProt\n")
        handle.write("UniProtKB\tP1\n")

    assert embedded_release_metadata(path) == {
        "gaf_version": "2.2",
        "date_generated": "2026-07-01",
        "generated_by": "UniProt",
    }


def test_extracts_obo_data_version(tmp_path):
    path = tmp_path / "ontology.obo"
    path.write_text(
        "format-version: 1.2\ndata-version: releases/2026-07-01\n\n[Term]\n",
        encoding="utf-8",
    )
    assert embedded_release_metadata(path) == {
        "format_version": "1.2",
        "data_version": "releases/2026-07-01",
    }


def test_extracts_uniprot_vocabulary_release(tmp_path):
    """UniProt's keywlist.txt / subcell.txt state their release in a header."""
    path = tmp_path / "subcell.txt"
    path.write_text(
        "Description: Controlled vocabulary of subcellular locations\n"
        "Name:        subcell.txt\n"
        "Release:     2026_02 of 10-Jun-2026\n",
        encoding="utf-8",
    )
    assert embedded_release_metadata(path) == {"release": "2026_02 of 10-Jun-2026"}


def test_extracts_expasy_colonless_release(tmp_path):
    """Expasy enzyme.dat writes 'CC   Release of <date>' with no colon."""
    path = tmp_path / "enzyme.dat"
    path.write_text(
        "CC   -------------------------------------------------------------\n"
        "CC   ENZYME nomenclature database\n"
        "CC   Release of 10-Jun-2026\n"
        "CC   Email: enzyme@expasy.org\n"
        "ID   1.1.1.1\n",
        encoding="utf-8",
    )
    # The e-mail line has a colon but is not a release field, so it is dropped.
    assert embedded_release_metadata(path) == {"release": "10-Jun-2026"}


def test_release_metadata_is_empty_for_unrecognized_formats(tmp_path):
    path = tmp_path / "relations.tsv"
    path.write_text("R-HSA-1\tR-HSA-2\n", encoding="utf-8")
    assert embedded_release_metadata(path) == {}


def test_json_parameters_stringifies_non_json_values(tmp_path):
    parameters = {
        "output_dir": Path("results"),
        "fdr_threshold": 0.01,
        "enable_true_path": False,
        "aspects": {"P"},
        "hierarchy_needs": ("go_obo",),
        "nested": {"path": Path("a/b")},
        "missing": None,
    }
    converted = json_parameters(parameters)
    assert converted["output_dir"] == "results"
    assert converted["aspects"] == ["P"]
    assert converted["hierarchy_needs"] == ["go_obo"]
    assert converted["nested"] == {"path": "a/b"}
    assert converted["missing"] is None
    # The whole point is that the manifest always serialises.
    json.dumps(converted)


def test_dependency_lock_is_hashed_when_present(tmp_path):
    lock = tmp_path / "uv.lock"
    lock.write_text("version = 1\n", encoding="utf-8")
    metadata = dependency_lock_metadata(tmp_path)
    assert metadata["available"] is True
    assert metadata["sha256"] == sha256_file(lock)

    assert dependency_lock_metadata(tmp_path / "elsewhere") == {
        "available": False,
        "path": "uv.lock",
    }


def test_manifest_filename_is_per_ontology():
    assert manifest_filename("go") == "run_manifest_go.json"
    assert manifest_filename("ec") == "run_manifest_ec.json"


def test_manifest_transitions_from_running_to_completed(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    input_path.write_text("input", encoding="utf-8")
    output_path.write_text("output", encoding="utf-8")
    (tmp_path / "uv.lock").write_text("version = 1\n", encoding="utf-8")
    manifest_path = tmp_path / "results" / "run_manifest_go.json"

    manifest = RunManifest(
        manifest_path,
        repository=tmp_path,
        parameters={"output_dir": Path("results"), "fdr_threshold": 0.01},
        inputs=[describe_file(input_path, role="annotations")],
        analysis={"ontology": {"key": "go"}, "thresholds": {"fdr_threshold": 0.01}},
        command=["dcgo", "--fdr-threshold", "0.01"],
    )
    running = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert running["status"] == "running"
    assert running["parameters"]["output_dir"] == "results"
    assert running["analysis"]["ontology"]["key"] == "go"
    assert running["software"]["dependency_lock"]["sha256"] == sha256_file(
        tmp_path / "uv.lock"
    )
    assert running["outputs"] == []

    manifest.complete(
        outputs=[describe_file(output_path, role="associations")],
        summary={"significant_associations": 4},
    )
    completed = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert completed["status"] == "completed"
    assert completed["summary"]["significant_associations"] == 4
    assert completed["outputs"][0]["sha256"] == sha256_file(output_path)
    assert completed["started_at"] <= completed["completed_at"]


def test_manifest_overwrites_a_previous_completed_manifest(tmp_path):
    """A new run must not inherit the previous run's 'completed' status."""
    manifest_path = tmp_path / "run_manifest_ec.json"
    first = RunManifest(
        manifest_path, repository=tmp_path, parameters={}, inputs=[], command=["dcgo"]
    )
    first.complete(outputs=[], summary={"significant_associations": 1})
    assert json.loads(manifest_path.read_text())["status"] == "completed"

    RunManifest(
        manifest_path, repository=tmp_path, parameters={}, inputs=[], command=["dcgo"]
    )
    assert json.loads(manifest_path.read_text())["status"] == "running"


def test_git_metadata_is_explicit_when_not_a_repository(tmp_path):
    assert git_metadata(tmp_path) == {"available": False}


def test_git_metadata_separates_tracked_changes_from_untracked_files(tmp_path):
    """`dirty` must mean "the code differs from `commit`", nothing else."""

    def git(*args):
        subprocess.run(
            ["git", "-C", str(tmp_path), *args], check=True, capture_output=True
        )

    git("init", "--initial-branch=main")
    git("config", "user.email", "test@example.test")
    git("config", "user.name", "Test")
    (tmp_path / "tracked.py").write_text("x = 1\n", encoding="utf-8")
    git("add", "tracked.py")
    git("commit", "-m", "initial")

    clean = git_metadata(tmp_path)
    assert clean["dirty"] is False
    assert clean["untracked_files"] == 0
    assert clean["branch"] == "main"

    # Scratch output next to the checkout is not a change to the code.
    (tmp_path / "scratch.tsv").write_text("noise\n", encoding="utf-8")
    with_scratch = git_metadata(tmp_path)
    assert with_scratch["dirty"] is False
    assert with_scratch["untracked_files"] == 1

    # Editing a tracked file is.
    (tmp_path / "tracked.py").write_text("x = 2\n", encoding="utf-8")
    assert git_metadata(tmp_path)["dirty"] is True
