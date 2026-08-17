"""Unit tests for the HPA cell-type layer (elevated expression, flat vocabulary).

Organised around the elevation policy — detection cutoff, the 4×-mean
enhanced criterion — and the Ensembl-gene re-keying via the flat file's
``DR Ensembl`` fourth field.
"""

import zipfile

import pytest

from src.gene_mapping import parse_ensembl_gene_accession_map
from src.hpa_annotation_source import (
    HPACellTypeAnnotationSource,
    parse_hpa_single_cell,
)

# ENSG1: 10.0 in adipocytes vs mean(1.0, 1.0) elsewhere → elevated (10 ≥ 4×1).
#        The 1.0 values are expressed but not elevated (1 < 4×5.5).
# ENSG2: expressed everywhere at the same level → nothing elevated.
# ENSG3: below the detection cutoff everywhere → nothing expressed at all.
# ENSG4: elevated in two cell types at once — needs enough near-zero types
#        that neither high value drags the other's "mean of others" over the
#        bar (24 ≥ 4 × 24/7).
TSV = """Gene\tGene name\tCell type\tnCPM
ENSG1\tG1\tadipocytes\t10.0
ENSG1\tG1\thepatocytes\t1.0
ENSG1\tG1\tmicroglia\t1.0
ENSG2\tG2\tadipocytes\t5.0
ENSG2\tG2\thepatocytes\t5.0
ENSG2\tG2\tmicroglia\t5.0
ENSG3\tG3\tadipocytes\t0.4
ENSG3\tG3\thepatocytes\t0.0
ENSG4\tG4\tadipocytes\t24.0
ENSG4\tG4\thepatocytes\t24.0
ENSG4\tG4\tmicroglia\t0.0
ENSG4\tG4\tb-cells\t0.0
ENSG4\tG4\tt-cells\t0.0
ENSG4\tG4\tcone photoreceptor cells\t0.0
ENSG4\tG4\tmelanocytes\t0.0
ENSG4\tG4\tenterocytes\t0.0
"""

DAT = """AC   P00001;
DR   Ensembl; ENST0001.1; ENSP0001.1; ENSG1.14. [P00001-1]
DR   Ensembl; ENST0002.1; ENSP0002.1; ENSG1.14. [P00001-2]
//
AC   P00002;
DR   Ensembl; ENST0003.2; ENSP0003.2; ENSG2.3.
//
AC   P00004;
DR   Ensembl; ENST0004.2; ENSP0004.2; ENSG4.9.
//
AC   P00005;
DR   Ensembl; ENST0005.2; ENSP0005.2; ENSG4.9.
//
"""


@pytest.fixture
def tsv(tmp_path):
    path = tmp_path / "rna_single_cell_type.tsv"
    path.write_text(TSV)
    return path


@pytest.fixture
def dat(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(DAT)
    return path


class TestElevationPolicy:
    def test_only_elevated_pairs_are_annotations(self, tsv):
        gene_terms, _ = parse_hpa_single_cell(tsv)
        assert gene_terms["ENSG1"] == {"adipocytes"}

    def test_uniform_expression_is_not_elevated(self, tsv):
        # Expressed in every cell type at the same level: real expression,
        # zero specificity — the exact case the 4x-mean criterion excludes.
        gene_terms, _ = parse_hpa_single_cell(tsv)
        assert "ENSG2" not in gene_terms

    def test_below_detection_cutoff_is_nothing(self, tsv):
        gene_terms, _ = parse_hpa_single_cell(tsv)
        assert "ENSG3" not in gene_terms

    def test_a_gene_may_be_elevated_in_several_cell_types(self, tsv):
        gene_terms, _ = parse_hpa_single_cell(tsv)
        assert gene_terms["ENSG4"] == {"adipocytes", "hepatocytes"}

    def test_counts_audit_the_filter(self, tsv):
        _, counts = parse_hpa_single_cell(tsv)
        assert counts.n_pairs == 16
        assert counts.n_expressed == 8  # nCPM >= 1
        assert counts.n_elevated == 3
        assert (counts.n_genes, counts.n_genes_elevated) == (4, 2)

    def test_zip_distribution_is_read_in_place(self, tmp_path):
        path = tmp_path / "rna_single_cell_type.tsv.zip"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("rna_single_cell_type.tsv", TSV)
        gene_terms, _ = parse_hpa_single_cell(path)
        assert "ENSG1" in gene_terms

    def test_ntpm_header_from_older_releases_is_accepted(self, tmp_path):
        path = tmp_path / "old.tsv"
        path.write_text(TSV.replace("nCPM", "nTPM"))
        gene_terms, _ = parse_hpa_single_cell(path)
        assert gene_terms["ENSG1"] == {"adipocytes"}

    def test_missing_expression_column_fails_loudly(self, tmp_path):
        path = tmp_path / "bad.tsv"
        path.write_text("Gene\tGene name\tCell type\tTPM\nENSG1\tG1\tx\t1.0\n")
        with pytest.raises(ValueError, match="expression column"):
            parse_hpa_single_cell(path)

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_hpa_single_cell(tmp_path / "gone.tsv")


class TestEnsemblGeneMap:
    def test_fourth_dr_field_supplies_the_gene(self, dat):
        gene_map = parse_ensembl_gene_accession_map(dat)
        # Version suffix and isoform bracket both stripped; two transcripts of
        # the same gene collapse onto one accession.
        assert gene_map.targets("ENSG1") == {"P00001"}
        assert gene_map.targets("ENSG2") == {"P00002"}

    def test_shared_gene_is_one_to_many(self, dat):
        gene_map = parse_ensembl_gene_accession_map(dat)
        assert gene_map.targets("ENSG4") == {"P00004", "P00005"}
        assert gene_map.n_one_to_many == 1

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_ensembl_gene_accession_map(tmp_path / "gone.dat")


class TestHPACellTypeAnnotationSource:
    """End-to-end: expression rows in, accession-keyed cell-type names out."""

    def test_parse_produces_accession_keyed_cell_types(self, tsv, dat):
        source = HPACellTypeAnnotationSource(tsv, dat)
        parsed = source.parse()
        assert parsed["P00001"] == {"adipocytes"}
        assert parsed["P00004"] == {"adipocytes", "hepatocytes"}
        assert parsed["P00005"] == parsed["P00004"]  # one-to-many expansion

    def test_unmapped_gene_is_counted_not_silent(self, tsv, dat):
        source = HPACellTypeAnnotationSource(tsv, dat)
        parsed = source.parse()
        # ENSG2/ENSG3 never reach the remap (not elevated); every elevated
        # gene here maps, so the unmapped list is empty — assert the audit
        # trail exists and agrees.
        assert source.coverage.unmapped_values == []
        assert source.filter_counts.n_genes_elevated == 2
        assert set(parsed) == {"P00001", "P00004", "P00005"}

    def test_spec_is_flat_hpa_vocabulary(self, tsv, dat):
        spec = HPACellTypeAnnotationSource(tsv, dat).spec
        assert spec.ontology_id == "HPA-CellType"
        assert spec.term_prefix is None
