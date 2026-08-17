"""Unit tests for the ZFIN affected-anatomy layer (ZDB-GENE re-keyed to UniProt)."""

import pytest

from src.zfin_annotation_source import (
    ZFINAnatomyAnnotationSource,
    parse_pheno_gene_clean_data,
    parse_zfin_uniprot,
)


def row(
    gene="ZDB-GENE-1",
    e1_sub="",
    e1_super="ZFA:0001559",
    quality="PATO:0001592",
    tag="abnormal",
    e2_sub="",
    e2_super="",
):
    """One 25-column phenoGeneCleanData_fish row with the fields under test."""
    fields = [""] * 25
    fields[0], fields[1], fields[2] = "1", "sym", gene
    fields[3], fields[7] = e1_sub, e1_super
    fields[9], fields[11] = quality, tag
    fields[12], fields[16] = e2_sub, e2_super
    return "\t".join(fields)


@pytest.fixture
def pheno(tmp_path):
    path = tmp_path / "phenoGeneCleanData_fish.txt"
    path.write_text(
        "\n".join(
            [
                row(gene="ZDB-GENE-1", e1_super="ZFA:0001559"),
                row(gene="ZDB-GENE-1", e1_sub="ZFA:0009035", e1_super="ZFA:0001620"),
                # GO process entity: not anatomy, must be skipped.
                row(gene="ZDB-GENE-1", e1_super="GO:0048747"),
                # E2 columns are affected structures too.
                row(gene="ZDB-GENE-2", e1_super="ZFA:0000092", e2_super="ZFA:0000051"),
                # A hypothetical normal-tagged row is not an abnormality.
                row(gene="ZDB-GENE-2", e1_super="ZFA:0001559", tag="normal"),
                # Unmappable gene, for the coverage test.
                row(gene="ZDB-GENE-9", e1_super="ZFA:0001559"),
            ]
        )
        + "\n"
    )
    return path


@pytest.fixture
def uniprot_map(tmp_path):
    path = tmp_path / "uniprot.txt"
    path.write_text(
        "ZDB-GENE-1\tSO:0001217\tppardb\tA9C4A5\n"
        "ZDB-GENE-2\tSO:0001217\tdlc\tQ9IAT6\n"
        "ZDB-GENE-2\tSO:0001217\tdlc\tA4JYS0\n"
        # Non-gene ZDB objects must be ignored.
        "ZDB-BAC-111024-49\tSO:0000153\tCH73-34H11\tA0A8M2B4B3\n"
    )
    return path


class TestParsePhenoGeneCleanData:
    def test_sub_and_superterms_from_both_entities_are_kept(self, pheno):
        parsed = parse_pheno_gene_clean_data(pheno)
        assert parsed["ZDB-GENE-1"] == {
            "ZFA:0001559",
            "ZFA:0009035",
            "ZFA:0001620",
        }
        assert parsed["ZDB-GENE-2"] == {"ZFA:0000092", "ZFA:0000051"}

    def test_non_zfa_entities_are_skipped(self, pheno):
        terms = set().union(*parse_pheno_gene_clean_data(pheno).values())
        assert not any(term.startswith("GO:") for term in terms)

    def test_non_abnormal_rows_are_dropped(self, pheno):
        # The normal-tagged ZFA:0001559 row must not credit ZDB-GENE-2.
        assert "ZFA:0001559" not in parse_pheno_gene_clean_data(pheno)["ZDB-GENE-2"]

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_pheno_gene_clean_data(tmp_path / "gone.txt")


class TestParseZfinUniprot:
    def test_gene_rows_map_and_stay_one_to_many(self, uniprot_map):
        gene_map = parse_zfin_uniprot(uniprot_map)
        assert gene_map.targets("ZDB-GENE-1") == {"A9C4A5"}
        assert gene_map.targets("ZDB-GENE-2") == {"Q9IAT6", "A4JYS0"}
        assert gene_map.n_one_to_many == 1

    def test_non_gene_objects_are_ignored(self, uniprot_map):
        assert parse_zfin_uniprot(uniprot_map).targets("ZDB-BAC-111024-49") == set()

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_zfin_uniprot(tmp_path / "gone.txt")


class TestZFINAnatomyAnnotationSource:
    """End-to-end: EQ rows in, accession-keyed ZFA terms out."""

    def test_parse_produces_accession_keyed_zfa_terms(self, pheno, uniprot_map):
        source = ZFINAnatomyAnnotationSource(pheno, uniprot_map)
        parsed = source.parse()
        assert parsed["A9C4A5"] == {"ZFA:0001559", "ZFA:0009035", "ZFA:0001620"}
        assert parsed["Q9IAT6"] == {"ZFA:0000092", "ZFA:0000051"}
        assert parsed["A4JYS0"] == parsed["Q9IAT6"]

    def test_unmapped_gene_is_counted_not_silent(self, pheno, uniprot_map):
        source = ZFINAnatomyAnnotationSource(pheno, uniprot_map)
        source.parse()
        assert source.coverage.unmapped_values == ["ZDB-GENE-9"]
        assert source.coverage.n_mapped_values == 2

    def test_spec_declares_the_zfa_prefix(self, pheno, uniprot_map):
        spec = ZFINAnatomyAnnotationSource(pheno, uniprot_map).spec
        assert spec.term_prefix == "ZFA:"
        assert spec.ontology_id == "ZFA"
