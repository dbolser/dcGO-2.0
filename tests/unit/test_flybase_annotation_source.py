"""Unit tests for the FlyBase phenotype layers (FBal→FBgn→UniProt)."""

import pytest

from src.flybase_annotation_source import (
    FBBT_SPEC,
    FBCV_SPEC,
    FlyBasePhenotypeAnnotationSource,
    parse_fbal_to_fbgn,
    parse_fbgn_uniprot,
    parse_genotype_phenotype,
)

GENOTYPES = (
    "## FlyBase genotype-phenotype report\n"
    "#genotype_symbols\tgenotype_FBids\tphenotype_name\tphenotype_id\t"
    "qualifier_names\tqualifier_ids\treference\n"
    # Single allele, FBcv phenotype class.
    "a[1]\tFBal0000001\tviable\tFBcv:0000349\t\t\tFBrf1\n"
    "a[1]\tFBal0000001\tabnormal learning\tFBcv:0000397\tadult\tFBdv:1\tFBrf1\n"
    # Same allele, FBbt anatomy row.
    "a[1]\tFBal0000001\tpigment cell\tFBbt:00004230\t\t\tFBrf1\n"
    # Homozygote: two tokens, one distinct allele — still single-allele.
    "b[1]/b[1]\tFBal0000002/FBal0000002\tfertile\tFBcv:0000374\t\t\tFBrf2\n"
    # Multi-allele genotype: dropped.\n
    "a[1]/b[1]\tFBal0000001/FBal0000002\tlethal\tFBcv:0000351\t\t\tFBrf3\n"
    # Space-separated multi-locus genotype: dropped.\n
    "a[1] c[1]\tFBal0000001 FBal0000003\tviable\tFBcv:0000349\t\t\tFBrf4\n"
    # Non-allele id (aberration): dropped.\n
    "Df(1)x\tFBab0000010\tlethal\tFBcv:0000351\t\t\tFBrf5\n"
    # Allele with no FBgn mapping, for the counting test.\n
    "z[1]\tFBal0000099\tviable\tFBcv:0000349\t\t\tFBrf6\n"
)

FBAL_TO_FBGN = (
    "# Generated from FlyBase release FB2026_02\n"
    "#\tAlleleID\tAlleleSymbol\tGeneID\tGeneSymbol\n"
    "FBal0000001\ta[1]\tFBgn0000010\ta\n"
    "FBal0000002\tb[1]\tFBgn0000020\tb\n"
    "FBal0000003\tc[1]\tFBgn0000030\tc\n"
)

FBGN_UNIPROT = (
    "## FlyBase FBgn-Major Accessions Table\n"
    "## gene_symbol\torganism\tprimary_FBgn#\tna_acc\tna_prot\tuniprot\tentrez\n"
    "a\tDmel\tFBgn0000010\tAQ1\t\tP10001\t1\n"
    "a\tDmel\tFBgn0000010\tAQ2\t\t\t1\n"
    "b\tDmel\tFBgn0000020\tAQ3\t\tP10002\t2\n"
    "b\tDmel\tFBgn0000020\tAQ4\t\tQ90002\t2\n"
)


@pytest.fixture
def genotypes(tmp_path):
    path = tmp_path / "genotype_phenotype_data.tsv"
    path.write_text(GENOTYPES)
    return path


@pytest.fixture
def fbal_map(tmp_path):
    path = tmp_path / "fbal_to_fbgn.tsv"
    path.write_text(FBAL_TO_FBGN)
    return path


@pytest.fixture
def fbgn_map(tmp_path):
    path = tmp_path / "fbgn_NAseq_Uniprot.tsv"
    path.write_text(FBGN_UNIPROT)
    return path


class TestParseGenotypePhenotype:
    def test_single_allele_fbcv_rows_are_kept(self, genotypes):
        parsed = parse_genotype_phenotype(genotypes, "FBcv:")
        assert parsed["FBal0000001"] == {"FBcv:0000349", "FBcv:0000397"}
        assert parsed["FBal0000002"] == {"FBcv:0000374"}

    def test_homozygote_counts_as_single_allele(self, genotypes):
        assert "FBal0000002" in parse_genotype_phenotype(genotypes, "FBcv:")

    def test_multi_allele_genotypes_are_dropped(self, genotypes):
        parsed = parse_genotype_phenotype(genotypes, "FBcv:")
        assert "FBcv:0000351" not in parsed.get("FBal0000001", set())
        assert "FBal0000003" not in parsed

    def test_prefix_selects_the_vocabulary(self, genotypes):
        assert parse_genotype_phenotype(genotypes, "FBbt:") == {
            "FBal0000001": {"FBbt:00004230"}
        }

    def test_non_allele_ids_are_dropped(self, genotypes):
        assert "FBab0000010" not in parse_genotype_phenotype(genotypes, "FBcv:")

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_genotype_phenotype(tmp_path / "gone.tsv", "FBcv:")


class TestFlyBaseMappings:
    def test_fbal_to_fbgn_skips_comments(self, fbal_map):
        assert parse_fbal_to_fbgn(fbal_map) == {
            "FBal0000001": "FBgn0000010",
            "FBal0000002": "FBgn0000020",
            "FBal0000003": "FBgn0000030",
        }

    def test_fbgn_uniprot_skips_blank_accessions(self, fbgn_map):
        gene_map = parse_fbgn_uniprot(fbgn_map)
        assert gene_map.targets("FBgn0000010") == {"P10001"}
        assert gene_map.targets("FBgn0000020") == {"P10002", "Q90002"}
        assert gene_map.n_one_to_many == 1

    def test_missing_files_fail_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_fbal_to_fbgn(tmp_path / "gone.tsv")
        with pytest.raises(FileNotFoundError):
            parse_fbgn_uniprot(tmp_path / "gone.tsv")


class TestFlyBasePhenotypeAnnotationSource:
    """End-to-end: genotype rows in, accession-keyed FBcv/FBbt terms out."""

    def test_parse_produces_accession_keyed_fbcv_terms(
        self, genotypes, fbal_map, fbgn_map
    ):
        source = FlyBasePhenotypeAnnotationSource(
            genotypes, fbal_map, fbgn_map, spec=FBCV_SPEC
        )
        parsed = source.parse()
        assert parsed["P10001"] == {"FBcv:0000349", "FBcv:0000397"}
        assert parsed["P10002"] == {"FBcv:0000374"}
        assert parsed["Q90002"] == parsed["P10002"]  # one-to-many expansion

    def test_fbbt_spec_selects_the_anatomy_terms(self, genotypes, fbal_map, fbgn_map):
        source = FlyBasePhenotypeAnnotationSource(
            genotypes, fbal_map, fbgn_map, spec=FBBT_SPEC
        )
        assert source.parse() == {"P10001": {"FBbt:00004230"}}
        assert source.spec.term_prefix == "FBbt:"

    def test_unmapped_allele_and_gene_are_counted(self, genotypes, fbal_map, fbgn_map):
        source = FlyBasePhenotypeAnnotationSource(
            genotypes, fbal_map, fbgn_map, spec=FBCV_SPEC
        )
        source.parse()
        assert source.n_unmapped_alleles == 1  # FBal0000099 has no FBgn
        assert source.coverage.n_mapped_values == 2  # both known genes mapped
