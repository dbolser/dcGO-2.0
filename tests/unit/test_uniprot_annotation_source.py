"""Unit tests for UniProt-native annotation sources (DR cross-refs + keywords)."""

import gzip

import pytest

from src.annotation_source import AnnotationSource
from src.uniprot_annotation_source import (
    KEYWORD_SPEC,
    REACTOME_SPEC,
    UniProtCrossRefAnnotationSource,
    UniProtKeywordAnnotationSource,
    parse_uniprot_cross_refs,
    parse_uniprot_keywords,
    reactome_source,
)

# Two entries. P07327 has a secondary accession, Reactome/KEGG/GO cross-refs, and
# a multi-line KW block. Q00000 has one Reactome id and no keywords.
SAMPLE_UNIPROT_DAT = """ID   ADH1A_HUMAN             Reviewed;         375 AA.
AC   P07327; B2R5V5;
DR   Reactome; R-HSA-71384; Ethanol oxidation.
DR   Reactome; R-HSA-9033241; Peroxisomal protein import.
DR   KEGG; hsa:124; .
DR   GO; GO:0004022; F:alcohol dehydrogenase (NAD+) activity; IDA:UniProtKB.
KW   Cytoplasm; Metal-binding; NAD;
KW   Oxidoreductase; Zinc.
//
ID   TEST2_HUMAN             Reviewed;         100 AA.
AC   Q00000;
DR   Reactome; R-HSA-71384; Ethanol oxidation.
//
"""


@pytest.fixture
def uniprot_dat_file(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(SAMPLE_UNIPROT_DAT)
    return path


@pytest.fixture
def uniprot_dat_gz(tmp_path):
    path = tmp_path / "uniprot_sprot.dat.gz"
    with gzip.open(path, "wt") as f:
        f.write(SAMPLE_UNIPROT_DAT)
    return path


class TestParseCrossRefs:
    def test_reactome_ids_by_primary_accession(self, uniprot_dat_file):
        result = parse_uniprot_cross_refs(uniprot_dat_file, "Reactome")
        assert result["P07327"] == {"R-HSA-71384", "R-HSA-9033241"}
        assert result["Q00000"] == {"R-HSA-71384"}
        # Keyed by primary accession only, not the secondary B2R5V5.
        assert "B2R5V5" not in result

    def test_other_database_selects_only_that_db(self, uniprot_dat_file):
        assert parse_uniprot_cross_refs(uniprot_dat_file, "KEGG") == {
            "P07327": {"hsa:124"}
        }

    def test_go_cross_refs_available_too(self, uniprot_dat_file):
        # The same parser can harvest GO from DR lines (another UniProt-native vocab).
        assert parse_uniprot_cross_refs(uniprot_dat_file, "GO")["P07327"] == {
            "GO:0004022"
        }

    def test_gzip_supported(self, uniprot_dat_gz):
        assert parse_uniprot_cross_refs(uniprot_dat_gz, "Reactome")["P07327"] == {
            "R-HSA-71384",
            "R-HSA-9033241",
        }

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_uniprot_cross_refs(tmp_path / "nope.dat", "Reactome")


class TestParseKeywords:
    def test_multiline_keywords(self, uniprot_dat_file):
        result = parse_uniprot_keywords(uniprot_dat_file)
        assert result["P07327"] == {
            "Cytoplasm",
            "Metal-binding",
            "NAD",
            "Oxidoreductase",
            "Zinc",
        }

    def test_entry_without_keywords_absent(self, uniprot_dat_file):
        assert "Q00000" not in parse_uniprot_keywords(uniprot_dat_file)


class TestSources:
    def test_reactome_source_is_annotation_source(self, uniprot_dat_file):
        source = reactome_source(uniprot_dat_file)
        assert isinstance(source, AnnotationSource)
        assert source.spec is REACTOME_SPEC
        assert source.parse()["P07327"] == {"R-HSA-71384", "R-HSA-9033241"}

    def test_cross_ref_source_arbitrary_db(self, uniprot_dat_file):
        source = UniProtCrossRefAnnotationSource(
            uniprot_dat_file, "KEGG", REACTOME_SPEC
        )
        assert source.parse() == {"P07327": {"hsa:124"}}

    def test_keyword_source(self, uniprot_dat_file):
        source = UniProtKeywordAnnotationSource(uniprot_dat_file)
        assert isinstance(source, AnnotationSource)
        assert source.spec is KEYWORD_SPEC
        assert "Zinc" in source.parse()["P07327"]
