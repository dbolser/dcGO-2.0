"""Tests for the live species/ontology-generic input-resolution stage."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from src.runner import parse_run_request, resolve_inputs


def test_run_request_captures_generic_identity_and_input_fields(tmp_path: Path) -> None:
    enzyme_dat = tmp_path / "enzyme.dat"
    enzyme_dat.write_text("//\n")
    request = parse_run_request(
        [
            "--species",
            "mouse",
            "--ontology",
            "ec",
            "--domain-key",
            "ssf",
            "--enzyme-dat",
            str(enzyme_dat),
            "--output-dir",
            "mouse-ec",
            "--enable-true-path",
            "--disable-supra-domains",
        ]
    )

    assert request.species == "mouse"
    assert request.ontology == "ec"
    assert request.domain_key == "ssf"
    assert request.output_dir == Path("mouse-ec")
    assert request.enzyme_dat == enzyme_dat
    assert request.enable_true_path is True
    assert request.enable_supra_domains is False


def test_input_resolution_consumes_request_and_registry(tmp_path: Path) -> None:
    enzyme_dat = tmp_path / "enzyme.dat"
    enzyme_dat.write_text("//\n")
    request = parse_run_request(
        ["--ontology", "ec", "--enzyme-dat", str(enzyme_dat)]
    )

    resolved = resolve_inputs(request)

    assert resolved.ontology_entry.key == "ec"
    assert resolved.ontology_label == "ec"
    assert resolved.ontology_paths["enzyme_dat"] == enzyme_dat
    assert resolved.missing_inputs == ()
    assert resolved.true_path_unsupported is False


def test_xref_identity_uses_selected_database_label() -> None:
    request = parse_run_request(
        ["--ontology", "xref", "--xref-db", "KEGG"]
    )

    assert resolve_inputs(request).ontology_label == "kegg"


def test_run_request_is_immutable() -> None:
    request = parse_run_request([])

    with pytest.raises(FrozenInstanceError):
        request.species = "mouse"  # type: ignore[misc]
