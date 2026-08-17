"""Unit tests for the HPO gene → phenotype layer (NCBI GeneID re-keyed to UniProt)."""

import pytest

from src.hpo_annotation_source import HPOAnnotationSource, parse_genes_to_phenotype

# The real header, columns read by name so the frequency/disease tail is inert.
G2P = """ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id
10\tNAT2\tHP:0000007\tAutosomal recessive inheritance\t-\tOMIM:243400
10\tNAT2\tHP:0001939\tAbnormality of metabolism/homeostasis\t-\tOMIM:243400
20\tAAK1\tHP:0000007\tAutosomal recessive inheritance\t1/2\tOMIM:111111
99\tGONE\tHP:0009999\tNot mappable\t-\tOMIM:222222
"""

DAT = """AC   P10001;
DR   GeneID; 10; -.
//
AC   P10002;
DR   GeneID; 20; -.
//
"""


@pytest.fixture
def g2p(tmp_path):
    path = tmp_path / "genes_to_phenotype.txt"
    path.write_text(G2P)
    return path


@pytest.fixture
def dat(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(DAT)
    return path


class TestParseGenesToPhenotype:
    def test_terms_are_grouped_by_gene(self, g2p):
        assert parse_genes_to_phenotype(g2p) == {
            "10": {"HP:0000007", "HP:0001939"},
            "20": {"HP:0000007"},
            "99": {"HP:0009999"},
        }

    def test_columns_are_found_by_name_not_position(self, tmp_path):
        # A release that reorders or extends the columns must still parse.
        path = tmp_path / "reordered.txt"
        path.write_text("hpo_id\textra\tncbi_gene_id\nHP:0000001\tx\t7\n")
        assert parse_genes_to_phenotype(path) == {"7": {"HP:0000001"}}

    def test_missing_column_names_the_gap(self, tmp_path):
        path = tmp_path / "bad.txt"
        path.write_text("gene\tterm\n1\tHP:1\n")
        with pytest.raises(ValueError, match="ncbi_gene_id"):
            parse_genes_to_phenotype(path)

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_genes_to_phenotype(tmp_path / "gone.txt")


class TestHPOAnnotationSource:
    """End-to-end: GeneID-keyed TSV in, accession-keyed HP terms out."""

    def test_parse_produces_accession_keyed_hp_terms(self, g2p, dat):
        source = HPOAnnotationSource(g2p, dat)
        assert source.parse() == {
            "P10001": {"HP:0000007", "HP:0001939"},
            "P10002": {"HP:0000007"},
        }

    def test_unmapped_gene_is_counted_not_silent(self, g2p, dat):
        source = HPOAnnotationSource(g2p, dat)
        assert source.coverage is None
        source.parse()
        assert source.coverage.unmapped_terms == ["99"]
        assert source.coverage.n_mapped_terms == 2  # of 3 gene ids

    def test_spec_declares_the_hp_prefix(self, g2p, dat):
        spec = HPOAnnotationSource(g2p, dat).spec
        assert spec.term_prefix == "HP:"
        assert spec.ontology_id == "HP"
