"""Unit tests for the MGI mouse phenotype layer (MGI markers re-keyed to UniProt)."""

import pytest

from src.mgi_annotation_source import (
    MGIAnnotationSource,
    parse_mgi_genepheno,
    parse_mrk_swissprot,
)

# Real MGI_GenePheno.rpt shape: 8 unnamed columns, marker in column 7.
# The Atm|Rad50 row is a genuine multi-gene genotype from the live file.
GENEPHENO = (
    "Rb1<tm1Tyj>/Rb1<tm1Tyj>\tRb1<tm1Tyj>\tMGI:1857242\tinvolves: 129S2/SvPas\t"
    "MP:0000600\t12529408\tMGI:97874\tMGI:2166359\n"
    "Rb1<tm1Tyj>/Rb1<tm1Tyj>\tRb1<tm1Tyj>\tMGI:1857242\tinvolves: 129S2/SvPas\t"
    "MP:0001716\t16449662\tMGI:97874\tMGI:2166359\n"
    "Rbpj<tm1Kyo>/Rbpj<tm1Kyo>\tRbpj<tm1Kyo>\tMGI:1857411\tinvolves: 129S2\t"
    "MP:0001614\t15466160\tMGI:96522\tMGI:2166381\n"
    "Atm<tm1Awb>/Atm<+>\tAtm<+>|Atm<tm1Awb>\tMGI:1857132|MGI:5614069\t"
    "involves: 129/Sv\tMP:0002216\t24532689\tMGI:107202|MGI:109292\tMGI:5614077\n"
)

MRK = (
    "MGI:97874\tRb1\tO\tRB transcriptional corepressor 1\t89.61\t14\tP13405\n"
    "MGI:96522\tRbpj\tO\trecombination signal binding protein\t50.0\t5\t"
    "P31266 Q3TLR6\n"
    "MGI:1914088\t0610009L18Rik\tO\tRIKEN cDNA\t84.07\t11\tQ9CVY3\n"
)


@pytest.fixture
def genepheno(tmp_path):
    path = tmp_path / "MGI_GenePheno.rpt"
    path.write_text(GENEPHENO)
    return path


@pytest.fixture
def mrk(tmp_path):
    path = tmp_path / "MRK_SwissProt_TrEMBL.rpt"
    path.write_text(MRK)
    return path


class TestParseMGIGenePheno:
    def test_terms_are_grouped_by_marker(self, genepheno):
        assert parse_mgi_genepheno(genepheno) == {
            "MGI:97874": {"MP:0000600", "MP:0001716"},
            "MGI:96522": {"MP:0001614"},
        }

    def test_multi_gene_genotypes_are_dropped(self, genepheno):
        # The Atm|Rad50 double mutant must not credit either gene.
        assert "MGI:107202" not in parse_mgi_genepheno(genepheno)
        assert "MGI:107202|MGI:109292" not in parse_mgi_genepheno(genepheno)

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_mgi_genepheno(tmp_path / "gone.rpt")


class TestParseMrkSwissprot:
    def test_space_separated_accessions_stay_one_to_many(self, mrk):
        gene_map = parse_mrk_swissprot(mrk)
        assert gene_map.targets("MGI:96522") == {"P31266", "Q3TLR6"}
        assert gene_map.n_one_to_many == 1

    def test_unknown_marker_is_empty(self, mrk):
        assert parse_mrk_swissprot(mrk).targets("MGI:404") == set()

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_mrk_swissprot(tmp_path / "gone.rpt")


class TestMGIAnnotationSource:
    """End-to-end: genotype rows in, accession-keyed MP terms out."""

    def test_parse_produces_accession_keyed_mp_terms(self, genepheno, mrk):
        source = MGIAnnotationSource(genepheno, mrk)
        assert source.parse() == {
            "P13405": {"MP:0000600", "MP:0001716"},
            "P31266": {"MP:0001614"},
            "Q3TLR6": {"MP:0001614"},
        }

    def test_coverage_counts_the_expansion(self, genepheno, mrk):
        source = MGIAnnotationSource(genepheno, mrk)
        source.parse()
        assert source.coverage.n_mapped_values == 2  # both markers mapped
        assert source.coverage.n_expanded_annotations == 1  # Rbpj → 2 accessions

    def test_unmapped_marker_is_counted_not_silent(self, tmp_path, genepheno):
        mrk = tmp_path / "empty_mrk.rpt"
        mrk.write_text("MGI:97874\tRb1\tO\tdesc\t89.61\t14\tP13405\n")
        source = MGIAnnotationSource(genepheno, mrk)
        source.parse()
        assert source.coverage.unmapped_values == ["MGI:96522"]

    def test_spec_declares_the_mp_prefix(self, genepheno, mrk):
        spec = MGIAnnotationSource(genepheno, mrk).spec
        assert spec.term_prefix == "MP:"
        assert spec.ontology_id == "MP"
