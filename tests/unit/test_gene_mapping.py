"""Unit tests for the gene id → UniProt accession mapping layer.

Organised like the DOID re-key tests: around the ways a mapping can go wrong —
a gene with no accession, a gene carried by several accessions — plus the
coverage bookkeeping that must make each visible rather than silent.
"""

import pytest

from src.gene_mapping import parse_gene_accession_index, remap_gene_annotations

# A miniature Swiss-Prot flat file:
#   GeneID 10       → one accession (the ordinary case)
#   GeneID 20       → two accessions (one-to-many)
#   HGNC:4 / SYMB4  → an entry reachable by HGNC id and symbol
#   (GeneID 99 appears in annotations but nowhere here: unmapped)
MINI_DAT = """AC   P10001; Q99999;
DR   GeneID; 10; -.
DR   HGNC; HGNC:1; NAT2.
//
AC   P10002;
DR   GeneID; 20; -.
DR   HGNC; HGNC:2; AAK1.
//
AC   P10003;
DR   GeneID; 20; -.
//
AC   P10004;
DR   HGNC; HGNC:4; SYMB4.
//
"""


@pytest.fixture
def dat(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(MINI_DAT)
    return path


@pytest.fixture
def index(dat):
    return parse_gene_accession_index(dat)


class TestParseGeneAccessionIndex:
    def test_geneid_maps_to_the_primary_accession(self, index):
        assert index.geneid.targets("10") == {"P10001"}

    def test_one_gene_on_two_entries_keeps_both(self, index):
        assert index.geneid.targets("20") == {"P10002", "P10003"}
        assert index.geneid.n_one_to_many == 1

    def test_hgnc_id_and_symbol_are_both_indexed(self, index):
        assert index.hgnc.targets("HGNC:4") == {"P10004"}
        assert index.symbol.targets("SYMB4") == {"P10004"}

    def test_unknown_gene_maps_to_nothing(self, index):
        assert index.geneid.targets("99") == set()

    def test_every_index_counts_the_same_scan(self, index):
        assert index.geneid.n_entries == index.hgnc.n_entries == 4
        assert (len(index.geneid), len(index.hgnc), len(index.symbol)) == (2, 3, 3)

    def test_missing_flat_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_gene_accession_index(tmp_path / "gone.dat")


class TestRemapGeneAnnotations:
    GENE_TERMS = {
        "10": {"T:1"},
        "20": {"T:1", "T:2"},
        "99": {"T:3"},
    }

    def test_genes_are_replaced_and_expanded(self, index):
        remapped, _ = remap_gene_annotations(self.GENE_TERMS, index.geneid, "test")
        assert remapped == {
            "P10001": {"T:1"},
            "P10002": {"T:1", "T:2"},
            "P10003": {"T:1", "T:2"},
        }

    def test_term_carried_only_by_an_unmapped_gene_leaves_the_layer(self, index):
        remapped, _ = remap_gene_annotations(self.GENE_TERMS, index.geneid, "test")
        assert not any("T:3" in terms for terms in remapped.values())

    def test_coverage_counts_the_gene_axis(self, index):
        # The coverage's "values" are the gene ids being remapped, its "keys"
        # are the ontology terms (see remap_gene_annotations).
        _, coverage = remap_gene_annotations(self.GENE_TERMS, index.geneid, "test")
        assert coverage.n_source_values == 3  # distinct gene ids
        assert coverage.n_mapped_values == 2
        assert coverage.unmapped_values == ["99"]
        assert coverage.n_source_keys == 3  # distinct ontology terms
        assert coverage.n_result_keys == 2  # T:3 left the layer

    def test_annotation_coverage_is_over_pairs_not_genes(self, index):
        _, coverage = remap_gene_annotations(self.GENE_TERMS, index.geneid, "test")
        assert coverage.n_source_annotations == 4  # gene-term pairs
        assert coverage.n_mapped_annotations == 3
        assert coverage.annotation_coverage == pytest.approx(3 / 4)

    def test_one_to_many_expansion_is_counted(self, index):
        _, coverage = remap_gene_annotations(self.GENE_TERMS, index.geneid, "test")
        # Both of gene 20's annotations landed on two accessions each.
        assert coverage.n_expanded_annotations == 2

    def test_empty_input_is_not_a_division_by_zero(self, index):
        remapped, coverage = remap_gene_annotations({}, index.geneid, "test")
        assert remapped == {}
        assert coverage.annotation_coverage == 0.0
