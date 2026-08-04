"""Unit tests for UniProt-native annotation sources (DR cross-refs + keywords)."""

import gzip

import pytest

from src.annotation_source import AnnotationSource
from src.uniprot_annotation_source import (
    DISEASE_SPEC,
    KEYWORD_SPEC,
    REACTOME_SPEC,
    UniProtCofactorAnnotationSource,
    UniProtCrossRefAnnotationSource,
    UniProtKeywordAnnotationSource,
    UniProtLigandAnnotationSource,
    UniProtRheaAnnotationSource,
    UniProtSubcellularAnnotationSource,
    disease_source,
    parse_binding_ligands,
    parse_catalysed_reactions,
    parse_cofactors,
    parse_subcell_vocabulary,
    parse_subcellular_locations,
    parse_uniprot_accessions,
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


# An entry exercising the layers curated into the entry *body* rather than DR
# lines: a wrapped SUBCELLULAR LOCATION comment with an isoform tag, a topology
# qualifier and a Note tail; COFACTOR and CATALYTIC ACTIVITY ChEBI/Rhea
# cross-references; and FT binding sites naming their ligand.
SAMPLE_BODY_DAT = """ID   ADH1A_HUMAN             Reviewed;         375 AA.
AC   P07327;
CC   -!- CATALYTIC ACTIVITY:
CC       Reaction=a primary alcohol + NAD(+) = an aldehyde + NADH + H(+);
CC         Xref=Rhea:RHEA:10736, ChEBI:CHEBI:15378, ChEBI:CHEBI:57540;
CC         EC=1.1.1.1; Evidence={ECO:0000269|PubMed:2738060};
CC   -!- COFACTOR:
CC       Name=Zn(2+); Xref=ChEBI:CHEBI:29105;
CC         Evidence={ECO:0000269|PubMed:11274460};
CC   -!- SUBCELLULAR LOCATION: Cytoplasm {ECO:0000269|PubMed:1}. Cell membrane;
CC       Single-pass type II membrane protein. Note=Also seen in the nucleus.
CC   -!- SIMILARITY: Belongs to the zinc-containing alcohol dehydrogenase family.
CC   ---------------------------------------------------------------------------
CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
CC   ---------------------------------------------------------------------------
FT   BINDING         47
FT                   /ligand="Zn(2+)"
FT                   /ligand_id="ChEBI:CHEBI:29105"
FT                   /ligand_label="1"
FT   BINDING         48..52
FT                   /ligand="NAD(+)"
FT                   /ligand_id="ChEBI:CHEBI:57540"
//
ID   ISO_HUMAN               Reviewed;         100 AA.
AC   Q00001;
CC   -!- SUBCELLULAR LOCATION: [Isoform 2]: Nucleus. Cytoplasm.
CC   ---------------------------------------------------------------------------
CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
CC   ---------------------------------------------------------------------------
//
ID   BARE_HUMAN              Reviewed;         100 AA.
AC   Q00002;
CC   -!- FUNCTION: Nothing to harvest here.
//
"""

# A miniature subcell.txt: a location (ID), a topology (IT), and the parents
# they hang from, with the SL content lines the CC comments actually carry.
SAMPLE_SUBCELL = """ID   Cytoplasm.
AC   SL-0086
SL   Cytoplasm.
//
ID   Nucleus.
AC   SL-0191
SY   Cell nucleus.
SL   Nucleus.
//
ID   Cell membrane.
AC   SL-0039
SL   Cell membrane.
HI   Membrane.
//
ID   Membrane.
AC   SL-0162
SL   Membrane.
//
IT   Single-pass type II membrane protein.
AC   SL-9903
SL   Single-pass type II membrane protein.
HI   Single-pass membrane protein.
//
IT   Single-pass membrane protein.
AC   SL-9904
SL   Single-pass membrane protein.
//
"""


@pytest.fixture
def body_dat_file(tmp_path):
    path = tmp_path / "body.dat"
    path.write_text(SAMPLE_BODY_DAT)
    return path


@pytest.fixture
def subcell_file(tmp_path):
    path = tmp_path / "subcell.txt"
    path.write_text(SAMPLE_SUBCELL)
    return path


class TestSubcellVocabulary:
    def test_names_and_accessions(self, subcell_file):
        vocab = parse_subcell_vocabulary(subcell_file)
        assert vocab.name_of["SL-0086"] == "Cytoplasm"
        assert vocab.accession_of["cytoplasm"] == "SL-0086"

    def test_synonyms_and_topologies_are_indexed(self, subcell_file):
        vocab = parse_subcell_vocabulary(subcell_file)
        assert vocab.accession_of["cell nucleus"] == "SL-0191"
        assert vocab.accession_of["single-pass type ii membrane protein"] == "SL-9903"

    def test_hierarchy_is_rekeyed_to_accessions(self, subcell_file):
        vocab = parse_subcell_vocabulary(subcell_file)
        assert vocab.child_to_parents["SL-0039"] == {"SL-0162"}
        assert vocab.child_to_parents["SL-9903"] == {"SL-9904"}

    def test_terms_without_parents_absent_from_hierarchy(self, subcell_file):
        assert "SL-0086" not in parse_subcell_vocabulary(subcell_file).child_to_parents


class TestSubcellularLocations:
    def test_locations_mapped_to_sl_accessions(self, body_dat_file, subcell_file):
        result = parse_subcellular_locations(body_dat_file, subcell_file)
        assert result["P07327"] == {"SL-0086", "SL-0039", "SL-9903"}

    def test_isoform_prefix_stripped(self, body_dat_file, subcell_file):
        result = parse_subcellular_locations(body_dat_file, subcell_file)
        assert result["Q00001"] == {"SL-0191", "SL-0086"}

    def test_note_text_is_not_harvested(self, body_dat_file, subcell_file):
        # "Note=Also seen in the nucleus" must not add the Nucleus term.
        assert (
            "SL-0191"
            not in parse_subcellular_locations(body_dat_file, subcell_file)["P07327"]
        )

    def test_entry_without_locations_absent(self, body_dat_file, subcell_file):
        assert "Q00002" not in parse_subcellular_locations(body_dat_file, subcell_file)

    def test_trailing_copyright_block_is_not_read_as_location_text(
        self, body_dat_file, subcell_file
    ):
        # Every entry ends with a CC-tagged copyright block; it must not be
        # appended to whichever comment topic happened to come last.
        result = parse_subcellular_locations(body_dat_file, subcell_file)
        assert result["Q00001"] == {"SL-0191", "SL-0086"}


class TestLigandsCofactorsReactions:
    def test_binding_site_ligands(self, body_dat_file):
        assert parse_binding_ligands(body_dat_file) == {
            "P07327": {"CHEBI:29105", "CHEBI:57540"}
        }

    def test_cofactors_only_from_the_cofactor_block(self, body_dat_file):
        # CHEBI:57540 appears in CATALYTIC ACTIVITY, not COFACTOR.
        assert parse_cofactors(body_dat_file) == {"P07327": {"CHEBI:29105"}}

    def test_catalysed_reactions(self, body_dat_file):
        assert parse_catalysed_reactions(body_dat_file) == {"P07327": {"RHEA:10736"}}

    def test_entries_without_the_layer_are_absent(self, body_dat_file):
        for parsed in (
            parse_binding_ligands(body_dat_file),
            parse_cofactors(body_dat_file),
            parse_catalysed_reactions(body_dat_file),
        ):
            assert "Q00002" not in parsed


class TestBodyLayerSources:
    def test_subcellular_source(self, body_dat_file, subcell_file):
        source = UniProtSubcellularAnnotationSource(body_dat_file, subcell_file)
        assert isinstance(source, AnnotationSource)
        assert source.spec.term_prefix == "SL-"
        assert "SL-0086" in source.parse()["P07327"]

    def test_ligand_source(self, body_dat_file):
        assert UniProtLigandAnnotationSource(body_dat_file).parse() == {
            "P07327": {"CHEBI:29105", "CHEBI:57540"}
        }

    def test_cofactor_source(self, body_dat_file):
        assert UniProtCofactorAnnotationSource(body_dat_file).parse() == {
            "P07327": {"CHEBI:29105"}
        }

    def test_rhea_source(self, body_dat_file):
        assert UniProtRheaAnnotationSource(body_dat_file).parse() == {
            "P07327": {"RHEA:10736"}
        }


class TestTermFromIdType:
    """Databases whose DR id is the accession and whose term is the third field."""

    SAMPLE = (
        "ID   T_HUMAN                 Reviewed;         100 AA.\n"
        "AC   P00001;\n"
        "DR   Pharos; P00001; Tbio.\n"
        "DR   CD-CODE; 91857CE7; Nucleolus.\n"
        "//\n"
        "ID   U_HUMAN                 Reviewed;         100 AA.\n"
        "AC   P00002;\n"
        "DR   Pharos; P00002; Tdark.\n"
        "DR   KEGG; hsa:124; -.\n"
        "//\n"
    )

    @pytest.fixture
    def dat(self, tmp_path):
        path = tmp_path / "typed.dat"
        path.write_text(self.SAMPLE)
        return path

    def test_third_field_becomes_the_term(self, dat):
        assert parse_uniprot_cross_refs(dat, "Pharos", term_from_id_type=True) == {
            "P00001": {"Tbio"},
            "P00002": {"Tdark"},
        }

    def test_default_still_uses_the_id(self, dat):
        assert parse_uniprot_cross_refs(dat, "Pharos") == {
            "P00001": {"P00001"},
            "P00002": {"P00002"},
        }

    def test_placeholder_third_field_is_not_a_term(self, dat):
        assert parse_uniprot_cross_refs(dat, "KEGG", term_from_id_type=True) == {}


class TestAccessionSet:
    """The t0 membership a temporal evaluation needs but a source cannot give."""

    def test_every_primary_accession(self, uniprot_dat_file):
        assert parse_uniprot_accessions(uniprot_dat_file) == {"P07327", "Q00000"}

    def test_secondary_accessions_excluded(self, uniprot_dat_file):
        # P07327's entry also lists B2R5V5; only the primary keys protein2ipr.
        assert "B2R5V5" not in parse_uniprot_accessions(uniprot_dat_file)

    def test_entries_without_annotations_still_count_as_present(self, body_dat_file):
        # BARE_HUMAN carries no harvestable term, but it existed in the release.
        assert "Q00002" in parse_uniprot_accessions(body_dat_file)

    def test_gzip_supported(self, uniprot_dat_gz):
        assert parse_uniprot_accessions(uniprot_dat_gz) == {"P07327", "Q00000"}
