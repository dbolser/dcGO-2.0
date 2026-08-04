"""The runner records the inputs the selected ontology actually declares.

These tests exist because the manifest's value is entirely in it being complete:
a per-ontology ``if/elif`` would silently omit inputs for the ontologies nobody
remembered to add. ``start_run_manifest`` instead reads ``needs`` /
``hierarchy_needs`` off the registry entry, so this checks that contract for a
hierarchy-carrying ontology (GO), a hierarchy-free one (EC), and one that lists
the same file in both roles (subcellular).
"""

import argparse
import json

import pytest

from run_dcgo_human import start_run_manifest
from src.ontology_registry import ONTOLOGIES, get_ontology


@pytest.fixture
def fake_inputs(tmp_path):
    """A stand-in file for every path the runner can resolve."""
    paths = {}
    for name, filename in {
        "gaf": "goa_human.gaf.gz",
        "go_obo": "go-basic.obo",
        "enzyme_dat": "enzyme.dat",
        "uniprot_dat": "uniprot_sprot.dat.gz",
        "reactome_relations": "ReactomePathwaysRelation.txt",
        "keywlist": "keywlist.txt",
        "subcell": "subcell.txt",
        "chebi_obo": "chebi_lite.obo",
    }.items():
        path = tmp_path / filename
        path.write_text(f"contents of {name}\n", encoding="utf-8")
        paths[name] = path
    return paths


def make_args(tmp_path, **overrides):
    defaults = dict(
        species="human",
        ontology="go",
        evidence_filter="manual",
        fdr_threshold=0.01,
        enable_true_path=False,
        enable_supra_domains=True,
        enable_shrinkage=False,
        shrinkage_strength=0.5,
        num_cores=4,
        batch_size=50000,
        output_dir=tmp_path / "out",
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def run_manifest_json(tmp_path, fake_inputs, **overrides):
    interpro = tmp_path / "protein2ipr_human.dat.gz"
    interpro.write_text("domains\n", encoding="utf-8")
    args = make_args(tmp_path, **overrides)
    entry = get_ontology(args.ontology)
    manifest = start_run_manifest(
        args,
        ontology_entry=entry,
        ontology_label=args.ontology,
        ontology_paths=fake_inputs,
        interpro_file=interpro,
    )
    return manifest, json.loads(manifest.path.read_text(encoding="utf-8"))


def test_go_run_records_the_gaf_and_nothing_else_without_true_path(
    tmp_path, fake_inputs
):
    manifest, data = run_manifest_json(tmp_path, fake_inputs, ontology="go")

    assert manifest.path.name == "run_manifest_go.json"
    assert data["status"] == "running"
    assert [record["role"] for record in data["inputs"]] == [
        "domain_annotations",
        "gaf",
    ]
    assert all(record["sha256"] for record in data["inputs"])
    assert data["analysis"]["ontology"]["hierarchy_inputs"] == []
    assert data["analysis"]["ontology"]["supports_true_path"] is True


def test_go_run_with_true_path_records_the_ontology_file(tmp_path, fake_inputs):
    _, data = run_manifest_json(
        tmp_path, fake_inputs, ontology="go", enable_true_path=True
    )

    roles = [record["role"] for record in data["inputs"]]
    assert roles == ["domain_annotations", "gaf", "go_obo"]
    obo = next(record for record in data["inputs"] if record["role"] == "go_obo")
    assert obo["path"] == str(fake_inputs["go_obo"])
    assert obo["sha256"]
    assert data["analysis"]["ontology"]["hierarchy_inputs"] == ["go_obo"]
    assert data["analysis"]["ontology"]["true_path_enabled"] is True
    assert (
        data["analysis"]["ontology"]["propagation"]
        == "go_dag_with_parental_background_filter"
    )


def test_ec_run_records_enzyme_dat_and_an_implicit_hierarchy(tmp_path, fake_inputs):
    _, data = run_manifest_json(
        tmp_path, fake_inputs, ontology="ec", enable_true_path=True
    )

    assert [record["role"] for record in data["inputs"]] == [
        "domain_annotations",
        "enzyme_dat",
    ]
    # EC's hierarchy is implicit in the numbering, so there is no extra input.
    assert data["analysis"]["ontology"]["hierarchy_inputs"] == []
    assert data["analysis"]["ontology"]["propagation"] == "ancestor_closure"


def test_an_input_used_for_both_roles_is_recorded_once(tmp_path, fake_inputs):
    """`subcellular` declares subcell.txt as annotation *and* hierarchy input."""
    _, data = run_manifest_json(
        tmp_path, fake_inputs, ontology="subcellular", enable_true_path=True
    )

    assert [record["role"] for record in data["inputs"]] == [
        "domain_annotations",
        "uniprot_dat",
        "subcell",
    ]
    assert data["analysis"]["ontology"]["hierarchy_inputs"] == ["subcell"]


def test_every_registered_ontology_can_be_recorded(tmp_path, fake_inputs):
    """No registered ontology may declare an input the runner cannot resolve."""
    for key, entry in ONTOLOGIES.items():
        declared = set(entry.needs) | set(entry.hierarchy_needs)
        assert declared <= set(fake_inputs), (
            f"--ontology {key} declares input(s) {declared - set(fake_inputs)} "
            "that run_dcgo_human.py does not resolve into ontology_paths"
        )


def test_thresholds_are_recorded(tmp_path, fake_inputs):
    _, data = run_manifest_json(
        tmp_path, fake_inputs, fdr_threshold=0.05, enable_shrinkage=True
    )

    thresholds = data["analysis"]["thresholds"]
    assert thresholds["fdr_threshold"] == 0.05
    assert thresholds["evidence_filter"] == "manual"
    assert thresholds["enable_shrinkage"] is True
    assert thresholds["shrinkage_strength"] == 0.5
    assert thresholds["max_supra_domain_length"] == 3
    assert thresholds["fisher_alternative"] == "greater"
    # Every parameter of the invocation is preserved verbatim as well.
    assert data["parameters"]["fdr_threshold"] == 0.05
    assert data["parameters"]["output_dir"] == str(tmp_path / "out")
