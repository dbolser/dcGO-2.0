"""Unit tests for the ontology registry (the --ontology dispatch table)."""

import pytest

from src.annotation_source import AnnotationSource
from src.ontology_registry import (
    ONTOLOGIES,
    describe_ontologies,
    get_ontology,
    missing_inputs,
    ontology_keys,
)


@pytest.fixture
def paths(tmp_path):
    """Every registered input, present and non-empty."""
    made = {}
    for name in (
        "gaf",
        "go_obo",
        "enzyme_dat",
        "uniprot_dat",
        "reactome_relations",
        "keywlist",
        "subcell",
        "chebi_obo",
        "doid_obo",
        "hpo_g2p",
        "hpo_obo",
        "syngo_zip",
    ):
        path = tmp_path / name
        path.write_text("")
        made[name] = path
    return made


class TestRegistryShape:
    def test_keys_are_self_consistent(self):
        assert all(key == entry.key for key, entry in ONTOLOGIES.items())

    def test_go_and_ec_are_still_registered(self):
        assert {"go", "ec"} <= set(ontology_keys())

    def test_every_entry_has_a_description(self):
        assert all(entry.description for entry in ONTOLOGIES.values())

    def test_uniprot_native_ontologies_declare_the_flat_file(self):
        for key in ("reactome", "keyword", "disease", "subcellular", "ligand", "rhea"):
            assert "uniprot_dat" in ONTOLOGIES[key].needs

    def test_describe_lists_every_ontology(self):
        described = describe_ontologies()
        assert all(key in described for key in ONTOLOGIES)

    def test_unknown_key_names_the_alternatives(self):
        with pytest.raises(KeyError, match="reactome"):
            get_ontology("not-an-ontology")


class TestTruePathSupport:
    @pytest.mark.parametrize(
        "key",
        [
            "go",
            "ec",
            "reactome",
            "keyword",
            "tcdb",
            "merops",
            "cazy",
            "subcellular",
            "ligand",
            "cofactor",
            "doid",
            "orphanet_doid",
            "hpo",
            "syngo",
        ],
    )
    def test_hierarchical_ontologies_can_propagate(self, key):
        assert get_ontology(key).supports_true_path

    @pytest.mark.parametrize("key", ["disease", "rhea", "drugbank", "xref", "pharos"])
    def test_flat_ontologies_cannot(self, key):
        assert not get_ontology(key).supports_true_path

    def test_only_go_propagates_outside_the_shared_engine(self):
        external = [k for k, e in ONTOLOGIES.items() if e.external_propagation]
        assert external == ["go"]


class TestImplicitHierarchies:
    """Ontologies whose hierarchy is encoded in the term id itself."""

    def test_tcdb_walks_up_the_tc_number(self, paths):
        ancestors = get_ontology("tcdb").build_ancestors(paths)
        assert list(ancestors("8.A.98.1.10")) == ["8.A.98.1", "8.A.98", "8.A", "8"]

    def test_merops_walks_family_then_catalytic_type(self, paths):
        ancestors = get_ontology("merops").build_ancestors(paths)
        assert list(ancestors("S01.151")) == ["S01", "S"]

    def test_cazy_walks_to_its_class(self, paths):
        assert list(get_ontology("cazy").build_ancestors(paths)("GT32")) == ["GT"]

    def test_ec_keeps_its_dash_padded_form(self, paths):
        ancestors = get_ontology("ec").build_ancestors(paths)
        assert list(ancestors("1.1.1.1")) == ["1.1.1.-", "1.1.-.-", "1.-.-.-"]


class TestMissingInputs:
    def test_nothing_missing_when_all_present(self, paths):
        assert missing_inputs(get_ontology("reactome"), paths) == []

    def test_absent_file_is_reported(self, paths, tmp_path):
        paths["uniprot_dat"] = tmp_path / "gone.dat"
        missing = missing_inputs(get_ontology("reactome"), paths)
        assert len(missing) == 1 and "uniprot_dat" in missing[0]

    def test_hierarchy_inputs_only_checked_when_propagating(self, paths, tmp_path):
        paths["reactome_relations"] = tmp_path / "gone.txt"
        entry = get_ontology("reactome")
        assert missing_inputs(entry, paths) == []
        assert missing_inputs(entry, paths, for_hierarchy=True)

    def test_unset_path_key_is_reported(self, paths):
        del paths["chebi_obo"]
        assert missing_inputs(get_ontology("ligand"), paths, for_hierarchy=True) == [
            "chebi_obo"
        ]


class TestSourceConstruction:
    def test_dr_vocabulary_targets_the_right_database(self, paths):
        source = get_ontology("orphanet").build_source(paths, {})
        assert isinstance(source, AnnotationSource)
        assert source.database == "Orphanet"

    def test_disease_filters_to_phenotype_mim_entries(self, paths):
        source = get_ontology("disease").build_source(paths, {})
        assert (source.database, source.id_type) == ("MIM", "phenotype")

    def test_pharos_reads_the_term_from_the_id_type_field(self, paths):
        assert get_ontology("pharos").build_source(paths, {}).term_from_id_type

    def test_go_source_honours_the_evidence_filter(self, paths):
        source = get_ontology("go").build_source(
            paths, {"evidence_filter": "experimental"}
        )
        assert source.evidence_filter == "experimental"

    def test_xref_takes_the_database_from_options(self, paths):
        source = get_ontology("xref").build_source(
            paths, {"xref_db": "BRENDA", "xref_type": None}
        )
        assert source.database == "BRENDA"
        assert source.spec.ontology_id == "BRENDA"

    def test_xref_passes_the_type_filter_through(self, paths):
        source = get_ontology("xref").build_source(
            paths, {"xref_db": "MIM", "xref_type": "gene"}
        )
        assert source.id_type == "gene"

    def test_subcellular_source_needs_both_inputs(self, paths):
        source = get_ontology("subcellular").build_source(paths, {})
        assert source.dat_path == paths["uniprot_dat"]
        assert source.subcell_path == paths["subcell"]

    def test_doid_rekeys_the_phenotype_mim_layer(self, paths):
        source = get_ontology("doid").build_source(paths, {})
        assert (source.database, source.id_type) == ("MIM", "phenotype")
        assert source.xref_prefix == "MIM"
        assert source.doid_obo_path == paths["doid_obo"]

    def test_orphanet_doid_rekeys_the_orphanet_layer(self, paths):
        source = get_ontology("orphanet_doid").build_source(paths, {})
        assert (source.database, source.id_type) == ("Orphanet", None)
        assert source.xref_prefix == "ORDO"

    def test_doid_needs_the_obo_even_without_true_path(self, paths, tmp_path):
        # The OBO supplies the mapping table, not just the hierarchy, so a run
        # without --enable-true-path must still fail loudly when it is absent.
        paths["doid_obo"] = tmp_path / "gone.obo"
        assert missing_inputs(get_ontology("doid"), paths)

    def test_hpo_reads_the_g2p_file_and_the_flat_file(self, paths):
        source = get_ontology("hpo").build_source(paths, {})
        assert source.genes_to_phenotype_path == paths["hpo_g2p"]
        assert source.dat_path == paths["uniprot_dat"]

    def test_syngo_reads_the_zip_and_the_flat_file(self, paths):
        source = get_ontology("syngo").build_source(paths, {})
        assert source.zip_path == paths["syngo_zip"]
        assert source.dat_path == paths["uniprot_dat"]


class TestGeneKeyedPreflight:
    """The gene-keyed layers must demand every input they read."""

    def test_hpo_needs_the_annotations_and_the_flat_file(self, paths, tmp_path):
        entry = get_ontology("hpo")
        assert set(entry.needs) == {"hpo_g2p", "uniprot_dat"}
        paths["hpo_g2p"] = tmp_path / "gone.txt"
        assert missing_inputs(entry, paths)

    def test_hpo_hierarchy_needs_the_obo_only_when_propagating(self, paths, tmp_path):
        entry = get_ontology("hpo")
        paths["hpo_obo"] = tmp_path / "gone.obo"
        assert missing_inputs(entry, paths) == []
        assert missing_inputs(entry, paths, for_hierarchy=True)

    def test_syngo_zip_is_both_annotation_and_hierarchy_input(self, paths, tmp_path):
        entry = get_ontology("syngo")
        assert "syngo_zip" in entry.needs
        assert entry.hierarchy_needs == ("syngo_zip",)
        paths["syngo_zip"] = tmp_path / "gone.zip"
        assert missing_inputs(entry, paths)
