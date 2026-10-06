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
    request = parse_run_request(["--ontology", "ec", "--enzyme-dat", str(enzyme_dat)])

    resolved = resolve_inputs(request)

    assert resolved.ontology_entry.key == "ec"
    assert resolved.ontology_label == "ec"
    assert resolved.ontology_paths["enzyme_dat"] == enzyme_dat
    assert resolved.missing_inputs == ()
    assert resolved.true_path_unsupported is False


def test_xref_identity_uses_selected_database_label() -> None:
    request = parse_run_request(["--ontology", "xref", "--xref-db", "KEGG"])

    assert resolve_inputs(request).ontology_label == "kegg"


def test_run_request_is_immutable() -> None:
    request = parse_run_request([])

    with pytest.raises(FrozenInstanceError):
        request.species = "mouse"  # type: ignore[misc]


def test_gaf_and_interpro_default_to_the_species_files() -> None:
    resolved = resolve_inputs(parse_run_request(["--species", "mouse"]))

    assert resolved.ontology_paths["gaf"] == Path(
        "data/raw/goa_annotations/goa_mouse.gaf.gz"
    )
    assert resolved.interpro_file == Path("data/interim/protein2ipr_mouse.dat.gz")


def test_gaf_and_interpro_overrides_replace_the_species_files(tmp_path: Path) -> None:
    """The temporal benchmark trains on archived files, not the species' current ones."""
    gaf = tmp_path / "goa_human.gaf.205.gz"
    gaf.write_text("!gaf-version: 2.2\n")
    interpro = tmp_path / "protein2ipr_human_t0.dat.gz"

    resolved = resolve_inputs(
        parse_run_request(["--gaf", str(gaf), "--interpro", str(interpro)])
    )

    assert resolved.ontology_paths["gaf"] == gaf
    assert resolved.interpro_file == interpro
    assert resolved.missing_inputs == ()


def test_a_missing_gaf_override_is_caught_by_the_early_input_check(
    tmp_path: Path,
) -> None:
    resolved = resolve_inputs(
        parse_run_request(["--gaf", str(tmp_path / "absent.gaf.gz")])
    )

    assert any("absent.gaf.gz" in missing for missing in resolved.missing_inputs)
