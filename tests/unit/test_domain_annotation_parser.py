"""
Unit tests for domain annotation parser.

Tests parsing of InterPro protein2ipr.dat files and supra-domain generation.
"""

import pytest
import gzip
from pathlib import Path

from src.domain_annotation_parser import (
    DomainAnnotation,
    ProteinDomainArchitecture,
    DomainAnnotationParser,
)


class TestDomainAnnotation:
    """Test suite for DomainAnnotation dataclass."""

    def test_creation(self):
        """Test basic domain annotation creation."""
        annotation = DomainAnnotation(
            protein_id="P12345",
            interpro_id="IPR001234",
            interpro_name="Test Domain",
            signature_id="PF001234",
            start=100,
            end=200,
        )

        assert annotation.protein_id == "P12345"
        assert annotation.interpro_id == "IPR001234"
        assert annotation.start == 100
        assert annotation.end == 200

    def test_length_property(self):
        """Test domain length calculation."""
        annotation = DomainAnnotation(
            protein_id="P12345",
            interpro_id="IPR001234",
            interpro_name="Test Domain",
            signature_id="PF001234",
            start=100,
            end=200,
        )

        assert annotation.length == 101  # 200 - 100 + 1


class TestProteinDomainArchitecture:
    """Test suite for ProteinDomainArchitecture."""

    def test_creation(self):
        """Test architecture creation."""
        arch = ProteinDomainArchitecture(
            protein_id="P12345",
            single_domains=["IPR001", "IPR002", "IPR003"],
            supra_domains=["IPR001,IPR002", "IPR002,IPR003"],
            domain_annotations=[],
        )

        assert arch.protein_id == "P12345"
        assert len(arch.single_domains) == 3
        assert len(arch.supra_domains) == 2

    def test_all_domains_property(self):
        """Test all_domains property combines single and supra-domains."""
        arch = ProteinDomainArchitecture(
            protein_id="P12345",
            single_domains=["IPR001", "IPR002"],
            supra_domains=["IPR001,IPR002"],
            domain_annotations=[],
        )

        all_domains = arch.all_domains
        assert len(all_domains) == 3
        assert "IPR001" in all_domains
        assert "IPR002" in all_domains
        assert "IPR001,IPR002" in all_domains


class TestDomainAnnotationParser:
    """Test suite for DomainAnnotationParser class."""

    @pytest.fixture
    def sample_protein2ipr_file(self, tmp_path):
        """Create a sample protein2ipr.dat file for testing."""
        file_path = tmp_path / "protein2ipr.dat.gz"

        # Sample data in protein2ipr.dat format
        # Format: protein_id, interpro_id, interpro_name, signature_id, start, end
        data = [
            "P00001\tIPR001234\tDomain A\tPF001234\t10\t110",
            "P00001\tIPR005678\tDomain B\tPF005678\t150\t250",
            "P00002\tIPR001234\tDomain A\tPF001234\t20\t120",
            "P00002\tIPR009876\tDomain C\tPF009876\t200\t300",
            "P00002\tIPR005678\tDomain B\tPF005678\t350\t450",
            "P00003\tIPR001234\tDomain A\tPF001234\t5\t105",
            # Small domain that should be filtered
            "P00003\tIPR999999\tSmall\tPF999999\t110\t115",
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        return file_path

    def test_parser_initialization(self):
        """Test parser initialization with default parameters."""
        parser = DomainAnnotationParser()

        assert parser.max_supra_domain_length == 3
        assert parser.min_domain_length == 10
        assert parser.species_filter is None

    def test_parser_custom_parameters(self):
        """Test parser with custom parameters."""
        parser = DomainAnnotationParser(
            max_supra_domain_length=5, min_domain_length=20, species_filter={"P0"}
        )

        assert parser.max_supra_domain_length == 5
        assert parser.min_domain_length == 20
        assert parser.species_filter == {"P0"}

    def test_parse_file(self, sample_protein2ipr_file):
        """Test parsing of protein2ipr file."""
        parser = DomainAnnotationParser(min_domain_length=10)
        architectures = parser.parse_protein2ipr_file(sample_protein2ipr_file)

        # Should have 3 proteins
        assert len(architectures) == 3
        assert "P00001" in architectures
        assert "P00002" in architectures
        assert "P00003" in architectures

    def test_domain_filtering_by_length(self, sample_protein2ipr_file):
        """Test that small domains are filtered out."""
        parser = DomainAnnotationParser(min_domain_length=10)
        architectures = parser.parse_protein2ipr_file(sample_protein2ipr_file)

        # P00003 should have IPR001234 but not IPR999999 (too small)
        arch_p00003 = architectures["P00003"]
        assert "IPR001234" in arch_p00003.single_domains
        assert "IPR999999" not in arch_p00003.single_domains

    def test_protein_filter(self, sample_protein2ipr_file):
        """Test filtering to specific proteins."""
        parser = DomainAnnotationParser()
        protein_filter = {"P00001", "P00002"}

        architectures = parser.parse_protein2ipr_file(
            sample_protein2ipr_file, protein_filter=protein_filter
        )

        # Should only have filtered proteins
        assert len(architectures) == 2
        assert "P00001" in architectures
        assert "P00002" in architectures
        assert "P00003" not in architectures

    def test_max_proteins_limit(self, sample_protein2ipr_file):
        """Test limiting maximum number of proteins."""
        parser = DomainAnnotationParser()
        architectures = parser.parse_protein2ipr_file(
            sample_protein2ipr_file, max_proteins=2
        )

        # Should stop after 2 proteins
        assert len(architectures) == 2

    def test_supra_domain_generation(self):
        """Test generation of supra-domains."""
        parser = DomainAnnotationParser(max_supra_domain_length=3)

        # Test with 3 domains
        domains = ["IPR001", "IPR002", "IPR003"]
        supra_domains = parser._generate_supra_domains(domains)

        # Should generate:
        # - 2-domain: IPR001,IPR002 | IPR002,IPR003
        # - 3-domain: IPR001,IPR002,IPR003
        assert len(supra_domains) == 3
        assert "IPR001,IPR002" in supra_domains
        assert "IPR002,IPR003" in supra_domains
        assert "IPR001,IPR002,IPR003" in supra_domains

    def test_supra_domain_max_length(self):
        """Test that supra-domain length is limited correctly."""
        parser = DomainAnnotationParser(max_supra_domain_length=2)

        domains = ["IPR001", "IPR002", "IPR003", "IPR004"]
        supra_domains = parser._generate_supra_domains(domains)

        # Should only generate 2-domain combinations
        assert all(supra.count(",") == 1 for supra in supra_domains)
        assert len(supra_domains) == 3  # Adjacent pairs

    def test_supra_domain_single_domain(self):
        """Test supra-domain generation with single domain."""
        parser = DomainAnnotationParser(max_supra_domain_length=3)

        domains = ["IPR001"]
        supra_domains = parser._generate_supra_domains(domains)

        # Should not generate any supra-domains (need at least 2)
        assert len(supra_domains) == 0

    def test_get_protein_domain_map(self, sample_protein2ipr_file):
        """Test getting simple protein-domain mapping."""
        parser = DomainAnnotationParser(max_supra_domain_length=2)
        parser.parse_protein2ipr_file(sample_protein2ipr_file)

        protein_domain_map = parser.get_protein_domain_map()

        # Check that mapping includes both single and supra-domains
        assert "P00001" in protein_domain_map
        domains_p00001 = protein_domain_map["P00001"]

        # Should have single domains
        assert "IPR001234" in domains_p00001
        assert "IPR005678" in domains_p00001

        # Should have supra-domain
        assert "IPR001234,IPR005678" in domains_p00001

    def test_get_domain_statistics(self, sample_protein2ipr_file):
        """Test getting domain statistics."""
        parser = DomainAnnotationParser()
        parser.parse_protein2ipr_file(sample_protein2ipr_file)

        stats = parser.get_domain_statistics()

        assert stats["total_proteins"] == 3
        assert stats["total_unique_domains"] > 0
        assert "domain_counts" in stats

        # IPR001234 appears in all 3 proteins
        assert stats["domain_counts"]["IPR001234"] == 3

    def test_file_not_found(self):
        """Test handling of missing file."""
        parser = DomainAnnotationParser()

        with pytest.raises(FileNotFoundError):
            parser.parse_protein2ipr_file(Path("/nonexistent/file.dat.gz"))

    def test_malformed_lines(self, tmp_path):
        """Test handling of malformed input lines."""
        file_path = tmp_path / "malformed.dat.gz"

        data = [
            "P00001\tIPR001234\tDomain A\tPF001234\t10\t110",  # Valid
            "P00002\tIPR005678",  # Missing fields
            "",  # Empty line
            "P00003\tIPR009876\tDomain C\tPF009876\t200\t300",  # Valid
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        parser = DomainAnnotationParser()
        architectures = parser.parse_protein2ipr_file(file_path)

        # Should successfully parse valid lines
        assert len(architectures) >= 2

    def test_non_gzipped_file(self, tmp_path):
        """Test parsing of non-gzipped file."""
        file_path = tmp_path / "protein2ipr.dat"

        data = [
            "P00001\tIPR001234\tDomain A\tPF001234\t10\t110",
            "P00002\tIPR005678\tDomain B\tPF005678\t20\t120",
        ]

        with open(file_path, "w") as f:
            f.write("\n".join(data))

        parser = DomainAnnotationParser()
        architectures = parser.parse_protein2ipr_file(file_path)

        assert len(architectures) == 2


class TestIntegration:
    """Integration tests for domain annotation parsing workflow."""

    def test_complete_workflow(self, tmp_path):
        """Test complete workflow from file parsing to domain mapping."""
        # Create test file
        file_path = tmp_path / "test.dat.gz"

        data = [
            "PROT1\tIPR001\tDomain 1\tPF001\t10\t100",
            "PROT1\tIPR002\tDomain 2\tPF002\t150\t250",
            "PROT1\tIPR003\tDomain 3\tPF003\t300\t400",
            "PROT2\tIPR001\tDomain 1\tPF001\t10\t100",
            "PROT2\tIPR004\tDomain 4\tPF004\t150\t250",
        ]

        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        # Parse with supra-domains
        parser = DomainAnnotationParser(max_supra_domain_length=3)
        architectures = parser.parse_protein2ipr_file(file_path)

        # Verify architectures
        assert len(architectures) == 2

        # PROT1 should have 3 single domains
        assert len(architectures["PROT1"].single_domains) == 3

        # PROT1 should have supra-domains
        supra = architectures["PROT1"].supra_domains
        assert len(supra) > 0
        assert "IPR001,IPR002" in supra
        assert "IPR002,IPR003" in supra
        assert "IPR001,IPR002,IPR003" in supra

        # Test protein-domain map
        protein_domain_map = parser.get_protein_domain_map()
        assert len(protein_domain_map["PROT1"]) == 6  # 3 single + 3 supra

        # Test statistics
        stats = parser.get_domain_statistics()
        assert stats["total_proteins"] == 2
        assert stats["total_unique_domains"] == 4
