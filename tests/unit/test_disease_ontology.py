"""Unit tests for the OMIM/Orphanet → Disease Ontology re-keying.

The mapping is not one-to-one, so the tests are organised around the three ways
it can go wrong — a source id with no DO term, one with several, and one whose
DO term is obsolete — plus the coverage bookkeeping that has to make each of
those visible rather than silent.
"""

import pytest

from src.disease_ontology import (
    DiseaseOntologyAnnotationSource,
    build_doid_xref_map,
    remap_protein_terms,
)
from src.hierarchy import closure_ancestors, parse_obo_child_parents

# A miniature doid.obo exercising every case the real file contains:
#   100  → one live term (the ordinary case)
#   200  → two live terms (one-to-many)
#   300  → an obsolete term with replaced_by (rescued)
#   400  → an obsolete term with no replacement (dropped)
#   500  → an obsolete term replaced by another obsolete term (chained, rescued)
#   PS99 → a phenotypic series, which UniProt DR MIM lines never carry
#   (600 appears in the annotations but nowhere here: unmapped)
MINI_OBO = """format-version: 1.2
data-version: releases/2026-07-31/doid.obo

[Term]
id: DOID:1
name: disease

[Term]
id: DOID:10
name: nervous system disease
is_a: DOID:1 ! disease
xref: MIM:PS99

[Term]
id: DOID:100
name: ataxia
is_a: DOID:10 ! nervous system disease
xref: MIM:100
xref: ORDO:900

[Term]
id: DOID:101
name: myopathy
is_a: DOID:10 ! nervous system disease
xref: MIM:200

[Term]
id: DOID:102
name: neuropathy
is_a: DOID:10 ! nervous system disease
xref: MIM:200

[Term]
id: DOID:103
name: obsolete ataxia 2
xref: MIM:300
is_obsolete: true
replaced_by: DOID:100

[Term]
id: DOID:104
name: obsolete myopathy 7
xref: MIM:400
is_obsolete: true

[Term]
id: DOID:105
name: obsolete neuropathy 3
xref: MIM:500
is_obsolete: true
replaced_by: DOID:106

[Term]
id: DOID:106
name: obsolete neuropathy 3, again
is_obsolete: true
replaced_by: DOID:102

[Typedef]
id: part_of
name: part of
"""


@pytest.fixture
def obo(tmp_path):
    path = tmp_path / "doid.obo"
    path.write_text(MINI_OBO)
    return path


@pytest.fixture
def mapping(obo):
    return build_doid_xref_map(obo, "MIM")


class TestBuildXrefMap:
    def test_simple_one_to_one(self, mapping):
        assert mapping.targets("100") == {"DOID:100"}

    def test_one_to_many_keeps_every_target(self, mapping):
        # A MIM entry that DO split into two disease classes is a curation
        # statement, not noise: the protein gets both, and the DAG pools them.
        assert mapping.targets("200") == {"DOID:101", "DOID:102"}
        assert mapping.n_one_to_many == 1

    def test_obsolete_target_resolved_via_replaced_by(self, mapping):
        assert mapping.targets("300") == {"DOID:100"}
        assert mapping.n_obsolete_resolved == 2  # MIM:300 and the MIM:500 chain

    def test_obsolete_chain_followed_to_a_live_term(self, mapping):
        assert mapping.targets("500") == {"DOID:102"}

    def test_obsolete_without_replacement_is_dropped_not_kept(self, mapping):
        assert mapping.targets("400") == set()
        assert "400" in mapping.source_ids_without_target
        assert mapping.n_obsolete_dropped == 1

    def test_unknown_source_id_maps_to_nothing(self, mapping):
        assert mapping.targets("600") == set()

    def test_phenotypic_series_counted_separately(self, mapping):
        # "PS99" is a legitimate DO xref but can never match a UniProt DR MIM id.
        assert mapping.n_non_numeric == 1
        assert mapping.targets("PS99") == {"DOID:10"}

    def test_audit_counts_add_up(self, mapping):
        assert mapping.n_terms == 9
        assert mapping.n_obsolete_terms == 4
        assert mapping.n_xrefs == 7  # every "xref: MIM:" line, obsolete included
        assert mapping.n_source_ids == 6  # 100, 200, 300, 400, 500, PS99
        assert len(mapping) == 5  # 400 has no live target

    def test_obsolete_terms_stay_out_of_the_hierarchy(self, obo):
        # Which is exactly why replaced_by has to be resolved at mapping time:
        # an obsolete target would be an orphan with no ancestors to propagate to.
        child_to_parents = parse_obo_child_parents(obo, relations=())
        assert "DOID:103" not in child_to_parents
        assert closure_ancestors(child_to_parents)("DOID:100") == {
            "DOID:10",
            "DOID:1",
        }

    def test_orphanet_prefix_uses_the_same_machinery(self, obo):
        assert build_doid_xref_map(obo, "ORDO").targets("900") == {"DOID:100"}

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            build_doid_xref_map(tmp_path / "nope.obo")


class TestRemapProteinTerms:
    def test_terms_are_replaced_and_expanded(self, mapping):
        remapped, _ = remap_protein_terms(
            {"P1": {"100"}, "P2": {"200"}}, mapping, label="test"
        )
        assert remapped == {
            "P1": {"DOID:100"},
            "P2": {"DOID:101", "DOID:102"},
        }

    def test_distinct_omim_ids_pool_onto_one_do_term(self, mapping):
        # The whole point of parse-time re-keying: two sparse OMIM phenotypes
        # become one better-supported DO class before any Fisher test runs.
        remapped, coverage = remap_protein_terms(
            {"P1": {"100"}, "P2": {"300"}}, mapping, label="test"
        )
        assert remapped == {"P1": {"DOID:100"}, "P2": {"DOID:100"}}
        assert coverage.n_source_values == 2
        assert coverage.n_result_values == 1

    def test_unmapped_terms_are_dropped_but_counted(self, mapping):
        remapped, coverage = remap_protein_terms(
            {"P1": {"100", "600"}}, mapping, label="test"
        )
        assert remapped == {"P1": {"DOID:100"}}
        assert coverage.unmapped_values == ["600"]
        assert coverage.n_mapped_values == 1
        assert coverage.n_source_values == 2

    def test_protein_with_no_mappable_term_leaves_the_layer(self, mapping):
        remapped, coverage = remap_protein_terms(
            {"P1": {"100"}, "P2": {"600"}}, mapping, label="test"
        )
        assert set(remapped) == {"P1"}
        assert (coverage.n_source_keys, coverage.n_result_keys) == (2, 1)

    def test_coverage_is_reported_over_annotations_not_just_terms(self, mapping):
        # A rare unmapped id costs less than a widely used one; coverage has to
        # be weighted by use, or it overstates the damage (or hides it).
        protein_terms = {"P1": {"100"}, "P2": {"100"}, "P3": {"100"}, "P4": {"600"}}
        _, coverage = remap_protein_terms(protein_terms, mapping, label="test")
        assert coverage.value_coverage == pytest.approx(0.5)
        assert coverage.annotation_coverage == pytest.approx(0.75)

    def test_expansion_is_counted(self, mapping):
        _, coverage = remap_protein_terms({"P1": {"200"}}, mapping, label="test")
        assert coverage.n_expanded_annotations == 1
        assert coverage.n_source_annotations == 1
        assert coverage.n_result_annotations == 2

    def test_unmapped_terms_are_ordered_by_use(self, mapping):
        protein_terms = {"P1": {"600", "700"}, "P2": {"700"}}
        _, coverage = remap_protein_terms(protein_terms, mapping, label="test")
        assert coverage.unmapped_values == ["700", "600"]

    def test_empty_input_is_not_a_division_by_zero(self, mapping):
        remapped, coverage = remap_protein_terms({}, mapping, label="test")
        assert remapped == {}
        assert coverage.value_coverage == 0.0
        assert coverage.annotation_coverage == 0.0


class TestDiseaseOntologyAnnotationSource:
    """End-to-end over a two-entry flat file: DR MIM in, DOID out."""

    DAT = """AC   P00001;
DR   MIM; 100; phenotype.
DR   MIM; 999; gene.
//
AC   P00002;
DR   MIM; 200; phenotype.
DR   MIM; 600; phenotype.
//
AC   P00003;
DR   MIM; 600; phenotype.
//
"""

    @pytest.fixture
    def dat(self, tmp_path):
        path = tmp_path / "uniprot_sprot.dat"
        path.write_text(self.DAT)
        return path

    def test_parse_produces_doid_terms(self, dat, obo):
        source = DiseaseOntologyAnnotationSource(dat, obo)
        assert source.parse() == {
            "P00001": {"DOID:100"},
            "P00002": {"DOID:101", "DOID:102"},
        }

    def test_gene_typed_mim_links_are_still_excluded(self, dat, obo):
        # DR MIM carries both the gene and the phenotype entry; only the latter
        # is a disease. (999 is a gene id and has no DO term anyway.)
        assert "999" not in DiseaseOntologyAnnotationSource(dat, obo).parse()

    def test_coverage_is_exposed_after_parsing(self, dat, obo):
        source = DiseaseOntologyAnnotationSource(dat, obo)
        assert source.coverage is None
        source.parse()
        assert source.coverage.n_source_keys == 3
        assert source.coverage.n_result_keys == 2
        assert source.coverage.unmapped_values == ["600"]

    def test_spec_declares_the_doid_prefix(self, dat, obo):
        spec = DiseaseOntologyAnnotationSource(dat, obo).spec
        assert spec.term_prefix == "DOID:"
