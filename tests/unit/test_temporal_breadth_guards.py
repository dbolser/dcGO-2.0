"""Preflight guards for the breadth benchmark (validation/temporal_breadth.py).

These exist because each failure they catch used to surface only *after* the
architecture pass — the expensive part of the run, and unrelated to any of them.
The third guard is the one that matters most, because its failure is silent
rather than loud: comparing a snapshot against itself yields zero hits, which
reads as "this layer has no predictive signal" rather than "this was never
tested".
"""

from pathlib import Path

import pytest

from validation.temporal_breadth import snapshot_problems


ALL_KEYS = (
    "gaf",
    "go_obo",
    "uniprot_dat",
    "subcell",
    "chebi_obo",
    "reactome_relations",
    "keywlist",
    "doid_obo",
    "enzyme_dat",
)


@pytest.fixture
def snapshots(tmp_path: Path):
    """Resolvable, on-disk path maps whose annotation inputs differ by snapshot."""

    def touch(name: str) -> Path:
        path = tmp_path / name
        path.write_text("x")
        return path

    t0 = {key: touch(f"t0_{key}") for key in ALL_KEYS}
    # Only the two genuinely time-varying annotation sources are re-pointed,
    # mirroring the script: hierarchies are shared between snapshots on purpose.
    t1 = dict(t0, uniprot_dat=touch("t1_uniprot_dat"), gaf=touch("t1_gaf"))
    return t0, t1


class TestIdenticalSnapshots:
    def test_ec_with_one_enzyme_release_is_refused(self, snapshots):
        # The live case: EC annotations come from enzyme_dat alone, so a single
        # ENZYME release makes t0_raw == t1_raw and every hit count zero by
        # construction.
        t0, t1 = snapshots
        problems = snapshot_problems(["ec"], t0, t1)
        assert len(problems) == 1
        assert "snapshot against itself" in problems[0]
        assert "enzyme_dat" in problems[0]

    def test_ec_with_a_distinct_t0_release_is_allowed(self, snapshots, tmp_path):
        t0, t1 = snapshots
        archived = tmp_path / "enzyme_2021.dat"
        archived.write_text("x")
        assert snapshot_problems(["ec"], dict(t0, enzyme_dat=archived), t1) == []

    def test_ontologies_with_a_time_varying_source_are_allowed(self, snapshots):
        # reactome's annotations come from uniprot_dat, which does differ. Its
        # hierarchy input (reactome_relations) is shared, and that is fine —
        # propagating both snapshots up today's hierarchy is a stated
        # approximation, not an empty comparison.
        t0, t1 = snapshots
        assert snapshot_problems(["reactome", "go", "subcellular"], t0, t1) == []

    def test_a_shared_hierarchy_alone_does_not_trigger_the_guard(self, snapshots):
        t0, t1 = snapshots
        assert t0["reactome_relations"] == t1["reactome_relations"]
        assert snapshot_problems(["reactome"], t0, t1) == []


class TestUnresolvedAndMissingInputs:
    def test_input_the_script_cannot_resolve(self, snapshots):
        # The doid regression: paths_t0 duplicates build_ontology_paths, so a
        # registry entry added later is simply absent from it.
        t0, t1 = snapshots
        t0.pop("doid_obo")
        t1.pop("doid_obo")
        problems = snapshot_problems(["doid"], t0, t1)
        assert len(problems) == 1
        assert "doid_obo" in problems[0] and "does not resolve" in problems[0]

    def test_input_named_but_absent_from_disk(self, snapshots):
        t0, t1 = snapshots
        t0["uniprot_dat"] = Path("/nonexistent/uniprot_sprot.dat.gz")
        problems = snapshot_problems(["reactome"], t0, t1)
        assert any("missing t0 input" in p and "uniprot_dat" in p for p in problems)

    def test_absence_is_reported_per_snapshot(self, snapshots):
        t0, t1 = snapshots
        t1["uniprot_dat"] = Path("/nonexistent/uniprot_sprot.dat.gz")
        problems = snapshot_problems(["reactome"], t0, t1)
        assert any("missing t1 input" in p for p in problems)

    def test_hierarchy_inputs_are_only_required_when_propagating(self, snapshots):
        t0, t1 = snapshots
        t0["reactome_relations"] = Path("/nonexistent/rel.txt")
        t1["reactome_relations"] = Path("/nonexistent/rel.txt")
        assert snapshot_problems(["reactome"], t0, t1, check_hierarchy=True)
        assert snapshot_problems(["reactome"], t0, t1, check_hierarchy=False) == []

    def test_an_unresolvable_input_short_circuits_its_ontology(self, snapshots):
        """One ontology reporting twice would bury the actionable line."""
        t0, t1 = snapshots
        t0.pop("enzyme_dat")
        t1.pop("enzyme_dat")
        assert len(snapshot_problems(["ec"], t0, t1)) == 1

    def test_clean_input_reports_nothing(self, snapshots):
        t0, t1 = snapshots
        assert snapshot_problems([], t0, t1) == []
