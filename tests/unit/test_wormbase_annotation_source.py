"""Unit tests for the WormBase layers (WBGene re-keyed to UniProt).

Covers the phenotype GAF (NOT policy) and the expression-based anatomy GAF
(WBbt, whose qualifier column grades the call — Uncertain dropped).
"""

import gzip

import pytest

from src.gene_mapping import parse_idmapping_accession_map, parse_idmapping_accessions
from src.wormbase_annotation_source import (
    WormBaseAnatomyAnnotationSource,
    WormBasePhenotypeAnnotationSource,
    parse_wb_anatomy_association,
    parse_wb_phenotype_association,
)

# GAF 2.0 shape: header comments, WBGene in column 2, NOT qualifier in
# column 4, WBPhenotype term in column 5.
GAF = (
    "!gaf-version: 2.0\n"
    "!generated-by: WormBase\n"
    "WB\tWBGene00000001\taap-1\t\tWBPhenotype:0000061\tWB_REF:x\tIMP\t"
    "WB:WBVar1\tP\t\tY110A7A.10\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000001\taap-1\t\tWBPhenotype:0000295\tWB_REF:y\tIMP\t"
    "WB:WBVar2\tP\t\tY110A7A.10\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000001\taap-1\tNOT\tWBPhenotype:0000062\tWB:WBVar3\tIMP\t"
    "WB:Person\tP\t\tY110A7A.10\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000002\taat-1\t\tWBPhenotype:0000062\tWB:WBVar4\tIMP\t"
    "WB:Person\tP\t\tF27C8.1\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000099\tgone\t\tWBPhenotype:0000001\tWB:WBVar5\tIMP\t"
    "WB:Person\tP\t\tZ99.9\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
)

IDMAPPING = (
    "P41932\tUniProtKB-ID\t14331_CAEEL\n"
    "P41932\tWormBase\tWBGene00000001\n"
    "Q20655\tWormBase\tWBGene00000002\n"
    "Q20655\tGene_Name\taat-1\n"
    "A0A000\tWormBase\tWBGene00000001\n"
    # Isoform rows: first column carries a -N suffix.
    "Q20655-2\tWormBase_TRS\tF27C8.1b.1\n"
)


@pytest.fixture
def gaf(tmp_path):
    path = tmp_path / "phenotype_association.wb.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(GAF)
    return path


@pytest.fixture
def idmapping(tmp_path):
    path = tmp_path / "CAEEL_idmapping.dat"
    path.write_text(IDMAPPING)
    return path


class TestParseWBPhenotypeAssociation:
    def test_terms_are_grouped_by_gene(self, gaf):
        parsed = parse_wb_phenotype_association(gaf)
        assert parsed["WBGene00000001"] == {
            "WBPhenotype:0000061",
            "WBPhenotype:0000295",
        }
        assert parsed["WBGene00000002"] == {"WBPhenotype:0000062"}

    def test_not_qualified_rows_are_dropped(self, gaf):
        parsed = parse_wb_phenotype_association(gaf)
        assert "WBPhenotype:0000062" not in parsed["WBGene00000001"]

    def test_plain_text_file_also_parses(self, tmp_path):
        path = tmp_path / "assoc.wb"
        path.write_text(GAF)
        assert "WBGene00000001" in parse_wb_phenotype_association(path)

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_wb_phenotype_association(tmp_path / "gone.wb.gz")


class TestIdmappingAccessionMap:
    def test_only_the_requested_id_type_is_kept(self, idmapping):
        gene_map = parse_idmapping_accession_map(idmapping, "WormBase")
        assert gene_map.targets("aat-1") == set()  # Gene_Name row ignored
        assert gene_map.targets("WBGene00000002") == {"Q20655"}

    def test_shared_gene_is_one_to_many(self, idmapping):
        gene_map = parse_idmapping_accession_map(idmapping, "WormBase")
        assert gene_map.targets("WBGene00000001") == {"P41932", "A0A000"}
        assert gene_map.n_one_to_many == 1

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_idmapping_accession_map(tmp_path / "gone.dat.gz", "WormBase")


class TestIdmappingAccessions:
    """The accession-universe reader used by extract_species_interpro."""

    def test_isoform_suffixes_collapse_to_canonical(self, idmapping):
        # Q20655-2 must not survive as its own accession: protein2ipr is
        # keyed by canonical accession, so isoform ids can never match.
        assert parse_idmapping_accessions(idmapping) == {
            "P41932",
            "Q20655",
            "A0A000",
        }

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_idmapping_accessions(tmp_path / "gone.dat.gz")


# Anatomy GAF: same shape, WBbt terms, expression qualifiers in column 4
# (Certain / Enriched / Partial / Uncertain / blank).
ANATOMY_GAF = (
    "!gaf-version: 2.0\n"
    "WB\tWBGene00000001\taap-1\tCertain\tWBbt:0003679\tWB_REF:x\tIDA\t"
    "WB:Expr1\tA\t\tY110A7A.10\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000001\taap-1\tEnriched\tWBbt:0005772\tWB_REF:x\tIDA\t"
    "WB:Expr2\tA\t\tY110A7A.10\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000001\taap-1\tUncertain\tWBbt:0005812\tWB_REF:x\tIDA\t"
    "WB:Expr3\tA\t\tY110A7A.10\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000002\taat-1\t\tWBbt:0005772\tWB_REF:y\tIDA\t"
    "WB:Expr4\tA\t\tF27C8.1\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
    "WB\tWBGene00000002\taat-1\tPartial\tWBbt:0003679\tWB_REF:y\tIDA\t"
    "WB:Expr5\tA\t\tF27C8.1\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
)


@pytest.fixture
def anatomy_gaf(tmp_path):
    path = tmp_path / "anatomy_association.wb.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(ANATOMY_GAF)
    return path


class TestParseWBAnatomyAssociation:
    def test_terms_are_grouped_by_gene(self, anatomy_gaf):
        parsed = parse_wb_anatomy_association(anatomy_gaf)
        assert parsed["WBGene00000001"] == {"WBbt:0003679", "WBbt:0005772"}

    def test_uncertain_rows_are_dropped(self, anatomy_gaf):
        # Curators flagged the call as doubtful; it must not become evidence.
        parsed = parse_wb_anatomy_association(anatomy_gaf)
        assert "WBbt:0005812" not in parsed["WBGene00000001"]

    def test_uncertain_matches_pipe_separated_tokens(self, tmp_path):
        # Qualifiers are |-lists; "Enriched|Uncertain" must drop like a bare
        # "Uncertain", and a token merely containing the word must not.
        import gzip as _gzip

        gaf = (
            "!gaf-version: 2.0\n"
            "WB\tWBGene00000003\tabc-1\tEnriched|Uncertain\tWBbt:0005813\t"
            "WB_REF:z\tIDA\tWB:Expr6\tA\t\tX.1\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
            "WB\tWBGene00000003\tabc-1\tUncertainty\tWBbt:0005814\t"
            "WB_REF:z\tIDA\tWB:Expr7\tA\t\tX.1\tgene\ttaxon:6239\t20251105\tWB\t\t\n"
        )
        path = tmp_path / "anatomy_pipe.wb.gz"
        with _gzip.open(path, "wt") as handle:
            handle.write(gaf)
        parsed = parse_wb_anatomy_association(path)
        assert parsed["WBGene00000003"] == {"WBbt:0005814"}

    def test_blank_and_partial_qualifiers_are_kept(self, anatomy_gaf):
        parsed = parse_wb_anatomy_association(anatomy_gaf)
        assert parsed["WBGene00000002"] == {"WBbt:0005772", "WBbt:0003679"}

    def test_phenotype_terms_do_not_leak_in(self, anatomy_gaf):
        # The parser is prefix-guarded, so a WBPhenotype row in the wrong file
        # would be counted malformed rather than annotated.
        parsed = parse_wb_anatomy_association(anatomy_gaf)
        all_terms = {term for terms in parsed.values() for term in terms}
        assert all(term.startswith("WBbt:") for term in all_terms)


class TestWormBaseAnatomyAnnotationSource:
    """End-to-end: anatomy GAF rows in, accession-keyed WBbt terms out."""

    def test_parse_produces_accession_keyed_terms(self, anatomy_gaf, idmapping):
        source = WormBaseAnatomyAnnotationSource(anatomy_gaf, idmapping)
        parsed = source.parse()
        assert parsed["P41932"] == {"WBbt:0003679", "WBbt:0005772"}
        assert parsed["A0A000"] == parsed["P41932"]  # one-to-many expansion
        assert parsed["Q20655"] == {"WBbt:0005772", "WBbt:0003679"}

    def test_spec_declares_the_wbbt_prefix(self, anatomy_gaf, idmapping):
        spec = WormBaseAnatomyAnnotationSource(anatomy_gaf, idmapping).spec
        assert spec.term_prefix == "WBbt:"
        assert spec.ontology_id == "WBbt"


class TestWormBasePhenotypeAnnotationSource:
    """End-to-end: GAF rows in, accession-keyed WBPhenotype terms out."""

    def test_parse_produces_accession_keyed_terms(self, gaf, idmapping):
        source = WormBasePhenotypeAnnotationSource(gaf, idmapping)
        parsed = source.parse()
        assert parsed["P41932"] == {"WBPhenotype:0000061", "WBPhenotype:0000295"}
        assert parsed["A0A000"] == parsed["P41932"]  # one-to-many expansion
        assert parsed["Q20655"] == {"WBPhenotype:0000062"}

    def test_unmapped_gene_is_counted_not_silent(self, gaf, idmapping):
        source = WormBasePhenotypeAnnotationSource(gaf, idmapping)
        source.parse()
        assert source.coverage.unmapped_values == ["WBGene00000099"]
        assert source.coverage.n_mapped_values == 2

    def test_spec_declares_the_wbphenotype_prefix(self, gaf, idmapping):
        spec = WormBasePhenotypeAnnotationSource(gaf, idmapping).spec
        assert spec.term_prefix == "WBPhenotype:"
        assert spec.ontology_id == "WBPhenotype"
