"""Unit tests for UniProt-native annotation sources (DR cross-refs + keywords)."""

import gzip

import pytest

from src.annotation_source import AnnotationSource
from src.uniprot_annotation_source import (
    DISEASE_SPEC,
    KEYWORD_SPEC,
    REACTOME_SPEC,
    UniProtCrossRefAnnotationSource,
    UniProtKeywordAnnotationSource,
    disease_source,
    parse_keyword_hierarchy,
    parse_reactome_relations,
    parse_uniprot_cross_refs,
    parse_uniprot_keywords,
    reactome_source,
)

# Two entries. P07327 has a secondary accession, Reactome/KEGG/GO cross-refs, a
# multi-line KW block, and MIM links (one gene, one phenotype). Q00000 has one
# Reactome id, no keywords, and only a gene-typed MIM link.
SAMPLE_UNIPROT_DAT = """ID   ADH1A_HUMAN             Reviewed;         375 AA.
AC   P07327; B2R5V5;
DR   Reactome; R-HSA-71384; Ethanol oxidation.
DR   Reactome; R-HSA-9033241; Peroxisomal protein import.
DR   KEGG; hsa:124; .
DR   GO; GO:0004022; F:alcohol dehydrogenase (NAD+) activity; IDA:UniProtKB.
DR   MIM; 103700; gene.
DR   MIM; 300100; phenotype.
KW   Cytoplasm; Metal-binding; NAD;
KW   Oxidoreductase; Zinc.
//
ID   TEST2_HUMAN             Reviewed;         100 AA.
AC   Q00000;
DR   Reactome; R-HSA-71384; Ethanol oxidation.
DR   MIM; 999999; gene.
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


class TestIdTypeFilterAndDisease:
    def test_no_filter_keeps_all_mim(self, uniprot_dat_file):
        result = parse_uniprot_cross_refs(uniprot_dat_file, "MIM")
        assert result["P07327"] == {"103700", "300100"}
        assert result["Q00000"] == {"999999"}

    def test_phenotype_filter_drops_gene_entries(self, uniprot_dat_file):
        result = parse_uniprot_cross_refs(uniprot_dat_file, "MIM", id_type="phenotype")
        assert result["P07327"] == {"300100"}  # gene 103700 dropped
        assert "Q00000" not in result  # only a gene MIM link

    def test_disease_source(self, uniprot_dat_file):
        source = disease_source(uniprot_dat_file)
        assert isinstance(source, AnnotationSource)
        assert source.spec is DISEASE_SPEC
        assert source.parse() == {"P07327": {"300100"}}

    def test_cross_ref_source_passes_id_type(self, uniprot_dat_file):
        source = UniProtCrossRefAnnotationSource(
            uniprot_dat_file, "MIM", DISEASE_SPEC, id_type="phenotype"
        )
        assert source.parse() == {"P07327": {"300100"}}


class TestReactomeHierarchy:
    def test_parse_relations(self, tmp_path):
        p = tmp_path / "rel.txt"
        p.write_text("R-HSA-1\tR-HSA-2\nR-HSA-2\tR-HSA-3\n")
        assert parse_reactome_relations(p) == {
            "R-HSA-2": {"R-HSA-1"},
            "R-HSA-3": {"R-HSA-2"},
        }

    def test_species_prefix_filter(self, tmp_path):
        p = tmp_path / "rel.txt"
        p.write_text("R-HSA-1\tR-HSA-2\nR-MMU-1\tR-MMU-2\n")
        assert parse_reactome_relations(p, species_prefix="R-HSA-") == {
            "R-HSA-2": {"R-HSA-1"}
        }


class TestKeywordHierarchy:
    def test_parse_keyword_paths(self, tmp_path):
        p = tmp_path / "keywlist.txt"
        p.write_text(
            "ID   2Fe-2S.\n"
            "HI   Ligand: Iron; Iron-sulfur; 2Fe-2S.\n"
            "HI   Ligand: Metal-binding; 2Fe-2S.\n"
            "CA   Ligand.\n"
            "//\n"
            "ID   Kinase.\n"
            "HI   Molecular function: Transferase; Kinase.\n"
            "CA   Molecular function.\n"
            "//\n"
        )
        result = parse_keyword_hierarchy(p)
        assert result["2Fe-2S"] == {"Iron-sulfur", "Metal-binding"}
        assert result["Kinase"] == {"Transferase"}
