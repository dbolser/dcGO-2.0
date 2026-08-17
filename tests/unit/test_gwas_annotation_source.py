"""Unit tests for the GWAS Catalog → EFO layer (gene symbols re-keyed to UniProt).

Organised around the row policy: the genome-wide-significance filter, the two
intergenic forms (flag and flanking pair), the mapped-gene column choice, and
the URI → CURIE normalisation that has to agree with the efo.obo id forms.
"""

import zipfile

import pytest

from src.gwas_annotation_source import (
    GWASCatalogAnnotationSource,
    normalise_efo_id,
    parse_efo_child_parents,
    parse_gwas_associations,
    trait_uri_to_curie,
)

HEADER = "\t".join(
    [
        "PUBMEDID",
        "DISEASE/TRAIT",
        "REPORTED GENE(S)",
        "MAPPED_GENE",
        "INTERGENIC",
        "P-VALUE",
        "MAPPED_TRAIT",
        "MAPPED_TRAIT_URI",
    ]
)


def row(
    mapped_gene,
    p_value="1E-12",
    intergenic="0",
    trait_uri="http://www.ebi.ac.uk/efo/EFO_0004574",
    reported="SOME",
):
    return "\t".join(
        ["1", "trait", reported, mapped_gene, intergenic, p_value, "t", trait_uri]
    )


TSV = "\n".join(
    [
        HEADER,
        # Ordinary row: one gene, one EFO trait.
        row("CDK2"),
        # Same gene, second trait as an OBO PURL (MONDO namespace kept).
        row("CDK2", trait_uri="http://purl.obolibrary.org/obo/MONDO_0005148"),
        # Below genome-wide significance: dropped.
        row("CDK2", p_value="2E-6", trait_uri="http://www.ebi.ac.uk/efo/EFO_0000001"),
        # Exactly at the threshold: dropped (the filter is strict <).
        row("CDK2", p_value="5E-8", trait_uri="http://www.ebi.ac.uk/efo/EFO_0000001"),
        # Intergenic flag: dropped even though a gene is mapped.
        row("CDK2", intergenic="1"),
        # Flanking-pair notation: intergenic in effect, dropped. The spaced
        # dash is the marker; hyphenated symbols (HLA-B below) must survive.
        row("MIR3143 - RPL10P2"),
        # Multiple mapped genes, both separators, duplicates collapsed.
        row("HLA-B, BRCA1; BRCA1"),
        # Multi-trait row: both CURIEs land on the gene.
        row(
            "TP53",
            trait_uri="http://www.ebi.ac.uk/efo/EFO_0000616, "
            "http://purl.obolibrary.org/obo/HP_0000001",
        ),
        # No trait URI: dropped and counted.
        row("TP53", trait_uri=""),
        # Unparsable p-value: counted with the significance drops.
        row("TP53", p_value="NR"),
        "",
    ]
)

DAT = """AC   P24941;
DR   HGNC; HGNC:1771; CDK2.
//
AC   P38398;
DR   HGNC; HGNC:1100; BRCA1.
//
AC   P04637;
DR   HGNC; HGNC:11998; TP53.
//
"""

MINI_EFO_OBO = """format-version: 1.2

[Term]
id: efo:EFO_0000001
name: experimental factor

[Term]
id: efo:EFO_0004574
name: total cholesterol measurement
is_a: efo:EFO_0000001 ! experimental factor

[Term]
id: MONDO:0005148
name: type 2 diabetes mellitus
is_a: efo:EFO_0000001 ! experimental factor
"""


@pytest.fixture
def tsv(tmp_path):
    path = tmp_path / "gwas-catalog-associations.tsv"
    path.write_text(TSV)
    return path


@pytest.fixture
def dat(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(DAT)
    return path


class TestIdNormalisation:
    def test_uri_forms_become_curies(self):
        assert (
            trait_uri_to_curie("http://www.ebi.ac.uk/efo/EFO_0004574") == "EFO:0004574"
        )
        assert (
            trait_uri_to_curie("http://purl.obolibrary.org/obo/MONDO_0005148")
            == "MONDO:0005148"
        )
        assert (
            trait_uri_to_curie("http://www.orpha.net/ORDO/Orphanet_140162")
            == "Orphanet:140162"
        )

    def test_efo_idspace_prefix_is_stripped(self):
        assert normalise_efo_id("efo:EFO_0000001") == "EFO:0000001"

    def test_plain_curies_pass_through(self):
        assert normalise_efo_id("MONDO:0000001") == "MONDO:0000001"
        assert normalise_efo_id("OBA:0000015") == "OBA:0000015"


class TestParseGWASAssociations:
    def test_gene_trait_pairs_with_policy_counts(self, tsv):
        gene_terms, counts = parse_gwas_associations(tsv)
        assert gene_terms["CDK2"] == {"EFO:0004574", "MONDO:0005148"}
        assert counts.n_rows == 10
        assert counts.n_below_significance == 3  # 2E-6, the 5E-8 boundary, "NR"
        assert counts.n_intergenic == 1
        assert counts.n_flanking == 1
        assert counts.n_no_trait == 1
        assert counts.n_kept == 4

    def test_multi_gene_rows_credit_each_symbol(self, tsv):
        gene_terms, _ = parse_gwas_associations(tsv)
        assert gene_terms["HLA-B"] == {"EFO:0004574"}
        assert gene_terms["BRCA1"] == {"EFO:0004574"}

    def test_multi_trait_rows_credit_each_curie(self, tsv):
        gene_terms, _ = parse_gwas_associations(tsv)
        assert gene_terms["TP53"] == {"EFO:0000616", "HP:0000001"}

    def test_zip_distribution_is_read_in_place(self, tmp_path):
        path = tmp_path / "associations.zip"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("gwas-catalog-download.tsv", TSV)
        gene_terms, _ = parse_gwas_associations(path)
        assert "CDK2" in gene_terms

    def test_missing_columns_fail_loudly(self, tmp_path):
        path = tmp_path / "bad.tsv"
        path.write_text("A\tB\n1\t2\n")
        with pytest.raises(ValueError, match="MAPPED_GENE"):
            parse_gwas_associations(path)

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_gwas_associations(tmp_path / "gone.tsv")


class TestEFOHierarchy:
    def test_edges_are_normalised_on_both_sides(self, tmp_path):
        path = tmp_path / "efo.obo"
        path.write_text(MINI_EFO_OBO)
        child_parents = parse_efo_child_parents(path)
        # The annotation side produces "EFO:0004574"; the hierarchy must agree.
        assert child_parents["EFO:0004574"] == {"EFO:0000001"}
        assert child_parents["MONDO:0005148"] == {"EFO:0000001"}


class TestGWASCatalogAnnotationSource:
    """End-to-end: catalog rows in, accession-keyed trait CURIEs out."""

    def test_parse_produces_accession_keyed_traits(self, tsv, dat):
        source = GWASCatalogAnnotationSource(tsv, dat)
        parsed = source.parse()
        assert parsed["P24941"] == {"EFO:0004574", "MONDO:0005148"}
        assert parsed["P38398"] == {"EFO:0004574"}
        assert parsed["P04637"] == {"EFO:0000616", "HP:0000001"}

    def test_unmapped_symbol_is_counted_not_silent(self, tsv, dat):
        source = GWASCatalogAnnotationSource(tsv, dat)
        source.parse()
        # HLA-B has no HGNC line in the mini flat file.
        assert source.coverage.unmapped_values == ["HLA-B"]
        assert source.filter_counts.n_kept == 4

    def test_spec_declares_the_efo_prefix(self, tsv, dat):
        spec = GWASCatalogAnnotationSource(tsv, dat).spec
        assert spec.ontology_id == "EFO"
