"""
Unit tests for GOA (Gene Ontology Annotation) parser.

Tests parsing of GAF files with evidence code filtering.
"""

import pytest
import gzip
from pathlib import Path

from src.goa_parser import (
    GOAAnnotation,
    GOAParser,
    parse_goa_human,
    EXPERIMENTAL_EVIDENCE,
    MANUAL_EVIDENCE,
    ELECTRONIC_EVIDENCE,
    ALL_EVIDENCE,
)


class TestGOAAnnotation:
    """Test suite for GOAAnnotation dataclass."""

    def test_creation(self):
        """Test basic annotation creation."""
        annotation = GOAAnnotation(
            protein_id="P12345",
            go_term="GO:0008150",
            evidence_code="IDA",
            aspect="P",
            db_reference="PMID:12345678",
            taxon_id="taxon:9606",
        )

        assert annotation.protein_id == "P12345"
        assert annotation.go_term == "GO:0008150"
        assert annotation.evidence_code == "IDA"
        assert annotation.aspect == "P"

    def test_invalid_go_term(self):
        """Test that invalid GO term format raises error."""
        with pytest.raises(ValueError, match="Invalid GO term format"):
            GOAAnnotation(
                protein_id="P12345",
                go_term="INVALID:0008150",  # Should start with GO:
                evidence_code="IDA",
                aspect="P",
                db_reference="PMID:12345678",
                taxon_id="taxon:9606",
            )


class TestGOAParser:
    """Test suite for GOAParser class."""

    @pytest.fixture
    def sample_gaf_file(self, tmp_path):
        """Create a sample GAF 2.2 file for testing."""
        file_path = tmp_path / "sample.gaf.gz"

        # GAF 2.2 format (15+ columns)
        # Fields: DB, DB_Object_ID, DB_Object_Symbol, Qualifier, GO_ID, DB:Reference,
        #         Evidence_Code, With/From, Aspect, DB_Object_Name, DB_Object_Synonym,
        #         DB_Object_Type, Taxon, Date, Assigned_By
        data = [
            "!gaf-version: 2.2",
            "!header line",
            # Experimental evidence (IDA)
            "UniProtKB\tP12345\tPROT1\t\tGO:0008150\tPMID:12345678\tIDA\t\tP\tProtein 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            # Electronic evidence (IEA)
            "UniProtKB\tP12345\tPROT1\t\tGO:0005515\tPMID:11111111\tIEA\t\tF\tProtein 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            # Another protein with experimental evidence
            "UniProtKB\tP23456\tPROT2\t\tGO:0008150\tPMID:22222222\tIMP\t\tP\tProtein 2\t\tprotein\ttaxon:9606\t20210101\tGOA",
            # Computational evidence (ISS)
            "UniProtKB\tP23456\tPROT2\t\tGO:0005737\tPMID:33333333\tISS\t\tC\tProtein 2\t\tprotein\ttaxon:9606\t20210101\tGOA",
            # NOT qualifier (should be excluded)
            "UniProtKB\tP34567\tPROT3\tNOT\tGO:0008150\tPMID:44444444\tIDA\t\tP\tProtein 3\t\tprotein\ttaxon:9606\t20210101\tGOA",
            # Valid annotation for PROT3
            "UniProtKB\tP34567\tPROT3\t\tGO:0003674\tPMID:55555555\tTAS\t\tF\tProtein 3\t\tprotein\ttaxon:9606\t20210101\tGOA",
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        return file_path

    def test_parser_initialization(self):
        """Test parser initialization with default parameters."""
        parser = GOAParser()

        assert parser.evidence_codes is None  # No filtering
        assert parser.aspects == {"P", "F", "C"}
        assert parser.exclude_qualifiers is True

    def test_parser_with_evidence_filter(self):
        """Test parser with evidence code filtering."""
        parser = GOAParser(evidence_codes=EXPERIMENTAL_EVIDENCE)

        assert parser.evidence_codes == EXPERIMENTAL_EVIDENCE

    def test_parse_all_evidence(self, sample_gaf_file):
        """Test parsing with all evidence codes."""
        parser = GOAParser(evidence_codes=None)  # No filtering
        protein_go_map = parser.parse_gaf_file(sample_gaf_file)

        # Should have 3 proteins (NOT qualifier excluded)
        assert len(protein_go_map) == 3

        # P12345 should have 2 GO terms (IDA and IEA)
        assert len(protein_go_map["P12345"]) == 2
        assert "GO:0008150" in protein_go_map["P12345"]
        assert "GO:0005515" in protein_go_map["P12345"]

    def test_parse_experimental_only(self, sample_gaf_file):
        """Test parsing with experimental evidence only."""
        parser = GOAParser(evidence_codes=EXPERIMENTAL_EVIDENCE)
        protein_go_map = parser.parse_gaf_file(sample_gaf_file)

        # P12345 should only have IDA annotation, not IEA
        assert "GO:0008150" in protein_go_map["P12345"]
        assert "GO:0005515" not in protein_go_map["P12345"]

        # P23456 should have IMP but not ISS
        assert "GO:0008150" in protein_go_map["P23456"]
        assert "GO:0005737" not in protein_go_map["P23456"]

    def test_parse_manual_evidence(self, sample_gaf_file):
        """Test parsing with manual evidence (non-IEA)."""
        parser = GOAParser(evidence_codes=MANUAL_EVIDENCE)
        protein_go_map = parser.parse_gaf_file(sample_gaf_file)

        # Should exclude IEA but include IDA, IMP, ISS, TAS
        assert "GO:0008150" in protein_go_map["P12345"]  # IDA
        assert "GO:0005515" not in protein_go_map["P12345"]  # IEA excluded

    def test_aspect_filtering(self, sample_gaf_file):
        """Test filtering by GO aspect."""
        # Only biological process (P)
        parser = GOAParser(aspects={"P"})
        protein_go_map = parser.parse_gaf_file(sample_gaf_file)

        # Should only have P aspect terms
        for protein, go_terms in protein_go_map.items():
            for go_term in go_terms:
                # Find corresponding annotation
                ann = [a for a in parser.annotations if a.go_term == go_term][0]
                assert ann.aspect == "P"

    def test_qualifier_exclusion(self, sample_gaf_file):
        """Test exclusion of NOT qualifiers."""
        parser = GOAParser(exclude_qualifiers=True)
        protein_go_map = parser.parse_gaf_file(sample_gaf_file)

        # P34567 should not have the NOT-qualified GO:0008150
        if "P34567" in protein_go_map:
            assert "GO:0008150" not in protein_go_map["P34567"]

        # But should have the TAS annotation
        assert "GO:0003674" in protein_go_map["P34567"]

    def test_qualifier_inclusion(self, sample_gaf_file):
        """Test including NOT qualifiers when requested."""
        parser = GOAParser(exclude_qualifiers=False)
        protein_go_map = parser.parse_gaf_file(sample_gaf_file)

        # P34567 should now have the NOT-qualified annotation
        assert "GO:0008150" in protein_go_map["P34567"]

    def test_get_annotations_by_protein(self, sample_gaf_file):
        """Test retrieving annotations for specific protein."""
        parser = GOAParser()
        parser.parse_gaf_file(sample_gaf_file)

        annotations = parser.get_annotations_by_protein("P12345")

        assert len(annotations) == 2
        go_terms = {ann.go_term for ann in annotations}
        assert "GO:0008150" in go_terms
        assert "GO:0005515" in go_terms

    def test_get_proteins_with_go_term(self, sample_gaf_file):
        """Test retrieving proteins with specific GO term."""
        parser = GOAParser()
        parser.parse_gaf_file(sample_gaf_file)

        proteins = parser.get_proteins_with_go_term("GO:0008150")

        # Both P12345 and P23456 have GO:0008150
        assert "P12345" in proteins
        assert "P23456" in proteins

    def test_statistics(self, sample_gaf_file):
        """Test that parser collects statistics."""
        parser = GOAParser(evidence_codes=EXPERIMENTAL_EVIDENCE)
        parser.parse_gaf_file(sample_gaf_file)

        stats = parser.statistics

        assert "accepted_annotations" in stats
        assert "filtered_evidence" in stats
        assert "excluded_qualifiers" in stats
        assert stats["accepted_annotations"] > 0

    def test_malformed_lines(self, tmp_path):
        """Test handling of malformed GAF lines."""
        file_path = tmp_path / "malformed.gaf.gz"

        data = [
            "!gaf-version: 2.2",
            # Valid line
            "UniProtKB\tP12345\tPROT1\t\tGO:0008150\tPMID:12345678\tIDA\t\tP\tProtein 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            # Too few fields
            "UniProtKB\tP23456\tPROT2",
            # Empty line
            "",
            # Another valid line
            "UniProtKB\tP34567\tPROT3\t\tGO:0005515\tPMID:22222222\tIMP\t\tF\tProtein 3\t\tprotein\ttaxon:9606\t20210101\tGOA",
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        parser = GOAParser()
        protein_go_map = parser.parse_gaf_file(file_path)

        # Should successfully parse valid lines
        assert len(protein_go_map) >= 2
        assert parser.statistics["malformed_lines"] >= 1

    def test_non_gzipped_file(self, tmp_path):
        """Test parsing of non-gzipped GAF file."""
        file_path = tmp_path / "sample.gaf"

        data = [
            "!gaf-version: 2.2",
            "UniProtKB\tP12345\tPROT1\t\tGO:0008150\tPMID:12345678\tIDA\t\tP\tProtein 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
        ]

        with open(file_path, "w") as f:
            f.write("\n".join(data))

        parser = GOAParser()
        protein_go_map = parser.parse_gaf_file(file_path)

        assert len(protein_go_map) == 1
        assert "P12345" in protein_go_map

    def test_file_not_found(self):
        """Test handling of missing file."""
        parser = GOAParser()

        with pytest.raises(FileNotFoundError):
            parser.parse_gaf_file(Path("/nonexistent/file.gaf.gz"))


class TestParseGoaHuman:
    """Test suite for parse_goa_human convenience function."""

    @pytest.fixture
    def sample_gaf_file(self, tmp_path):
        """Create sample GAF file."""
        file_path = tmp_path / "goa_human.gaf.gz"

        data = [
            "!gaf-version: 2.2",
            "UniProtKB\tP12345\tPROT1\t\tGO:0008150\tPMID:12345678\tIDA\t\tP\tProtein 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tP12345\tPROT1\t\tGO:0005515\tPMID:11111111\tIEA\t\tF\tProtein 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tP23456\tPROT2\t\tGO:0008150\tPMID:22222222\tIMP\t\tP\tProtein 2\t\tprotein\ttaxon:9606\t20210101\tGOA",
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        return file_path

    def test_parse_all(self, sample_gaf_file):
        """Test parsing with 'all' evidence filter."""
        protein_go_map = parse_goa_human(sample_gaf_file, evidence_filter="all")

        # Should include IEA
        assert "GO:0005515" in protein_go_map["P12345"]

    def test_parse_manual(self, sample_gaf_file):
        """Test parsing with 'manual' evidence filter (default)."""
        protein_go_map = parse_goa_human(sample_gaf_file, evidence_filter="manual")

        # Should exclude IEA
        assert "GO:0008150" in protein_go_map["P12345"]  # IDA included
        assert "GO:0005515" not in protein_go_map["P12345"]  # IEA excluded

    def test_parse_experimental(self, sample_gaf_file):
        """Test parsing with 'experimental' evidence filter."""
        protein_go_map = parse_goa_human(
            sample_gaf_file, evidence_filter="experimental"
        )

        # Should only have experimental evidence
        assert "GO:0008150" in protein_go_map["P12345"]  # IDA
        assert "GO:0008150" in protein_go_map["P23456"]  # IMP

    def test_aspect_filtering(self, sample_gaf_file):
        """Test filtering by GO aspects."""
        # Only molecular function (F)
        protein_go_map = parse_goa_human(sample_gaf_file, aspects={"F"})

        # Should only have F aspect terms
        for protein, go_terms in protein_go_map.items():
            # GO:0005515 is F aspect, GO:0008150 is P aspect
            if "GO:0005515" in go_terms:
                assert "GO:0008150" not in go_terms


class TestEvidenceCodeConstants:
    """Test evidence code constant definitions."""

    def test_experimental_evidence(self):
        """Test experimental evidence codes."""
        assert "IDA" in EXPERIMENTAL_EVIDENCE
        assert "IMP" in EXPERIMENTAL_EVIDENCE
        assert "IEA" not in EXPERIMENTAL_EVIDENCE

    def test_electronic_evidence(self):
        """Test electronic evidence codes."""
        assert "IEA" in ELECTRONIC_EVIDENCE
        assert len(ELECTRONIC_EVIDENCE) == 1

    def test_manual_evidence(self):
        """Test manual evidence (non-IEA)."""
        assert "IDA" in MANUAL_EVIDENCE
        assert "IMP" in MANUAL_EVIDENCE
        assert "ISS" in MANUAL_EVIDENCE
        assert "TAS" in MANUAL_EVIDENCE
        assert "IEA" not in MANUAL_EVIDENCE

    def test_all_evidence(self):
        """Test all evidence codes."""
        # ALL should include everything
        assert "IDA" in ALL_EVIDENCE
        assert "IEA" in ALL_EVIDENCE
        assert "ISS" in ALL_EVIDENCE


class TestIntegration:
    """Integration tests for complete GOA parsing workflow."""

    def test_realistic_workflow(self, tmp_path):
        """Test realistic parsing workflow with various evidence codes."""
        file_path = tmp_path / "test.gaf.gz"

        # Create realistic GAF data
        data = [
            "!gaf-version: 2.2",
            "!Header lines",
            # Mix of evidence codes and aspects
            "UniProtKB\tPROT01\tGENE1\t\tGO:0008150\tPMID:1\tIDA\t\tP\tGene 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tPROT01\tGENE1\t\tGO:0003674\tPMID:2\tIEA\t\tF\tGene 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tPROT01\tGENE1\t\tGO:0005575\tPMID:3\tISS\t\tC\tGene 1\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tPROT02\tGENE2\t\tGO:0008150\tPMID:4\tIMP\t\tP\tGene 2\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tPROT02\tGENE2\tNOT\tGO:0003674\tPMID:5\tIDA\t\tF\tGene 2\t\tprotein\ttaxon:9606\t20210101\tGOA",
            "UniProtKB\tPROT03\tGENE3\t\tGO:0005575\tPMID:6\tTAS\t\tC\tGene 3\t\tprotein\ttaxon:9606\t20210101\tGOA",
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        # Parse with manual evidence (default)
        protein_go_map = parse_goa_human(file_path, evidence_filter="manual")

        # Verify results
        assert len(protein_go_map) == 3

        # PROT01 should have IDA and ISS, but not IEA
        assert "GO:0008150" in protein_go_map["PROT01"]
        assert "GO:0005575" in protein_go_map["PROT01"]
        assert "GO:0003674" not in protein_go_map["PROT01"]

        # PROT02 should have IMP, NOT qualifier excluded
        assert "GO:0008150" in protein_go_map["PROT02"]
        assert "GO:0003674" not in protein_go_map["PROT02"]

        # PROT03 should have TAS
        assert "GO:0005575" in protein_go_map["PROT03"]
