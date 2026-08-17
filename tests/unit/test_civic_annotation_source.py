"""Unit tests for the CIViC-derived cancer layers (NCIt and OncoTree).

Organised around the chain's three stages — evidence-row policy and
molecular-profile gene extraction, DOID → target cross-reference hops, and
the final symbol → UniProt re-key — each of which must count what it drops.
"""

import json

import pytest

from src.civic_annotation_source import (
    CIViCNCItAnnotationSource,
    CIViCOncoTreeAnnotationSource,
    build_doid_to_ncit_map,
    build_doid_to_oncotree_map,
    genes_of_molecular_profile,
    parse_civic_evidence,
    parse_oncotree,
)

HEADER = "\t".join(
    [
        "molecular_profile",
        "molecular_profile_id",
        "disease",
        "doid",
        "evidence_direction",
        "evidence_status",
    ]
)


def row(profile, doid, direction="Supports", status="accepted"):
    return "\t".join([profile, "1", "disease", doid, direction, status])


CIVIC_TSV = "\n".join(
    [
        HEADER,
        row("JAK2 V617F", "1037"),
        # Fusion: both partners credited.
        row("BCR::ABL1 Fusion", "8552"),
        # Complex profile: every AND component contributes its gene.
        row("BRAF V600E AND TP53 Mutation", "1037"),
        # Negative evidence: dropped and counted.
        row("JAK2 V617F", "8552", direction="Does Not Support"),
        # No direction at all.
        row("JAK2 V617F", "8552", direction="N/A"),
        # No disease id: cannot enter either target vocabulary.
        row("JAK2 V617F", ""),
        # Not accepted (never in the nightly; counted defensively).
        row("JAK2 V617F", "1037", status="submitted"),
        # Leading zeros are part of the DOID id and must survive.
        row("ERBB2 Amplification", "0060079"),
        "",
    ]
)

# Miniature doid.obo carrying NCI and UMLS_CUI xrefs (plus one obsolete term
# whose replacement carries the xref chain forward).
DOID_OBO = """format-version: 1.2

[Term]
id: DOID:1037
name: lymphoblastic leukemia
xref: NCI:C3167
xref: UMLS_CUI:C0023448

[Term]
id: DOID:8552
name: chronic myeloid leukemia
xref: NCI:C3174

[Term]
id: DOID:0060079
name: her2-receptor positive breast cancer

[Typedef]
id: part_of
"""

ONCOTREE_JSON = [
    {
        "code": "ALL",
        "name": "Acute Lymphoid Leukemia",
        "externalReferences": {"NCI": ["C3167"], "UMLS": ["C0023448"]},
        "parent": "LYMPH",
        "level": 2,
    },
    {
        "code": "LYMPH",
        "name": "Lymphoid",
        "externalReferences": {},
        "parent": "TISSUE",
        "level": 1,
    },
    # UMLS-only node: reachable only through the UMLS route.
    {
        "code": "UMLSONLY",
        "name": "Umls Only Type",
        "externalReferences": {"UMLS": ["C0023448"]},
        "parent": "TISSUE",
        "level": 1,
    },
]

DAT = """AC   O60674;
DR   HGNC; HGNC:6192; JAK2.
//
AC   P00519;
DR   HGNC; HGNC:76; ABL1.
//
AC   P15056;
DR   HGNC; HGNC:1097; BRAF.
//
AC   P04637;
DR   HGNC; HGNC:11998; TP53.
//
"""


@pytest.fixture
def civic(tmp_path):
    path = tmp_path / "nightly-ClinicalEvidenceSummaries.tsv"
    path.write_text(CIVIC_TSV)
    return path


@pytest.fixture
def doid_obo(tmp_path):
    path = tmp_path / "doid.obo"
    path.write_text(DOID_OBO)
    return path


@pytest.fixture
def oncotree_json(tmp_path):
    path = tmp_path / "oncotree_tumortypes.json"
    path.write_text(json.dumps(ONCOTREE_JSON))
    return path


@pytest.fixture
def dat(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(DAT)
    return path


class TestGeneExtraction:
    def test_first_token_is_the_gene(self):
        assert genes_of_molecular_profile("JAK2 V617F") == {"JAK2"}
        assert genes_of_molecular_profile("ACVR1 Gain-of-Function") == {"ACVR1"}

    def test_fusions_credit_both_partners(self):
        assert genes_of_molecular_profile("BCR::ABL1 Fusion") == {"BCR", "ABL1"}

    def test_and_or_components_each_contribute(self):
        assert genes_of_molecular_profile("BRAF V600E AND TP53 Mutation") == {
            "BRAF",
            "TP53",
        }
        assert genes_of_molecular_profile("EGFR L858R OR EGFR T790M") == {"EGFR"}


class TestParseCIViCEvidence:
    def test_gene_doid_pairs_with_policy_counts(self, civic):
        gene_doids, counts = parse_civic_evidence(civic)
        assert gene_doids["JAK2"] == {"DOID:1037"}
        assert gene_doids["ABL1"] == {"DOID:8552"}
        assert counts.n_rows == 8
        assert counts.n_not_supports == 2  # Does Not Support + N/A
        assert counts.n_not_accepted == 1
        assert counts.n_no_doid == 1
        assert counts.n_kept == 4

    def test_leading_zeros_survive_the_doid_prefix(self, civic):
        gene_doids, _ = parse_civic_evidence(civic)
        assert gene_doids["ERBB2"] == {"DOID:0060079"}

    def test_missing_file_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_civic_evidence(tmp_path / "gone.tsv")


class TestDiseaseMappings:
    def test_doid_to_ncit_inverts_the_do_xrefs(self, doid_obo):
        mapping = build_doid_to_ncit_map(doid_obo)
        assert mapping.targets("DOID:1037") == {"NCIT:C3167"}
        assert mapping.targets("DOID:8552") == {"NCIT:C3174"}
        # No NCI xref: unmapped, and the remap stage will count it.
        assert mapping.targets("DOID:0060079") == set()

    def test_doid_to_oncotree_unions_nci_and_umls_routes(self, doid_obo, oncotree_json):
        mapping = build_doid_to_oncotree_map(doid_obo, oncotree_json)
        # NCI:C3167 → ALL, plus UMLS_CUI:C0023448 → ALL and UMLSONLY.
        assert mapping.targets("DOID:1037") == {"ALL", "UMLSONLY"}
        # C3174 appears in no OncoTree node: OncoTree simply lacks the type.
        assert mapping.targets("DOID:8552") == set()

    def test_oncotree_hierarchy_chains_to_the_tissue_root(self, oncotree_json):
        vocabulary = parse_oncotree(oncotree_json)
        assert vocabulary.child_to_parents["ALL"] == {"LYMPH"}
        assert vocabulary.child_to_parents["LYMPH"] == {"TISSUE"}


class TestCIViCNCItAnnotationSource:
    """End-to-end: evidence rows in, accession-keyed NCIT terms out."""

    def test_parse_produces_accession_keyed_ncit_terms(self, civic, doid_obo, dat):
        source = CIViCNCItAnnotationSource(civic, doid_obo, dat)
        parsed = source.parse()
        assert parsed["O60674"] == {"NCIT:C3167"}  # JAK2
        assert parsed["P00519"] == {"NCIT:C3174"}  # ABL1
        assert parsed["P15056"] == {"NCIT:C3167"}  # BRAF (via DOID:1037)

    def test_every_stage_is_audited(self, civic, doid_obo, dat):
        source = CIViCNCItAnnotationSource(civic, doid_obo, dat)
        source.parse()
        assert source.filter_counts.n_kept == 4
        # DOID:0060079 has no NCI xref; ERBB2 leaves at the disease stage.
        assert source.disease_coverage.unmapped_values == ["DOID:0060079"]
        # BCR has no HGNC line in the mini flat file.
        assert source.coverage.unmapped_values == ["BCR"]

    def test_spec_declares_the_ncit_prefix(self, civic, doid_obo, dat):
        assert CIViCNCItAnnotationSource(civic, doid_obo, dat).spec.term_prefix == (
            "NCIT:"
        )


class TestCIViCOncoTreeAnnotationSource:
    def test_parse_produces_accession_keyed_codes(
        self, civic, doid_obo, oncotree_json, dat
    ):
        source = CIViCOncoTreeAnnotationSource(civic, doid_obo, oncotree_json, dat)
        parsed = source.parse()
        assert parsed["O60674"] == {"ALL", "UMLSONLY"}  # JAK2 via DOID:1037
        # ABL1's only disease (DOID:8552) has no OncoTree node, so the
        # protein leaves the layer entirely — and is counted, not silent.
        assert "P00519" not in parsed
        assert "DOID:8552" in source.disease_coverage.unmapped_values
