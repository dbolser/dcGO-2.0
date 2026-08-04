"""
Unit tests for domain annotation parser.

Tests parsing of InterPro protein2ipr.dat files and supra-domain generation.
"""

import pytest
import gzip
from pathlib import Path

from src.domain_annotation_parser import (
    DOMAIN_KEYS,
    DomainAnnotation,
    ProteinDomainArchitecture,
    DomainAnnotationParser,
    superfamily_sunid,
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


class TestSuperfamilySunid:
    """SSFnnnnn ↔ SCOP sunid, the join key for the published dcGO tables."""

    def test_extracts_sunid(self):
        assert superfamily_sunid("SSF53649") == 53649
        assert superfamily_sunid("SSF52540") == 52540

    def test_leading_zeros_and_short_ids(self):
        # dcGO's tables store bare integers, so the round trip has to be numeric
        # rather than string-prefix based.
        assert superfamily_sunid("SSF0053649") == 53649

    def test_rejects_other_member_databases(self):
        for signature in ("PF00001", "IPR000001", "G3DSA:1.10.10.10", "cd00001"):
            assert superfamily_sunid(signature) is None

    def test_rejects_structure_function_linkage_db(self):
        # SFLD signatures start with "SF" too — mistaking them for SUPERFAMILY
        # would inject non-SCOP domains into the §3 comparison.
        for signature in ("SFLDF00001", "SFLDS00001", "SFLDG01135"):
            assert superfamily_sunid(signature) is None

    def test_rejects_non_numeric_suffix(self):
        assert superfamily_sunid("SSFabcde") is None
        assert superfamily_sunid("SSF") is None


class TestSuperfamilyDomainKey:
    """The ``--domain-key ssf`` parse path (VALIDATION_PLAN §3)."""

    @pytest.fixture
    def mixed_member_db_file(self, tmp_path):
        """protein2ipr rows mixing SUPERFAMILY with other member databases.

        The SSF rows are deliberately *not* the positionally-adjacent ones: for
        PROT1 a Pfam-only entry sits between the two SUPERFAMILY matches, and the
        rows are listed out of positional order (as the real file often is).
        """
        file_path = tmp_path / "protein2ipr.dat.gz"
        data = [
            # protein, InterPro entry, name, signature, start, end
            "PROT1\tIPR000100\tSF domain A\tSSF50001\t10\t110",
            "PROT1\tIPR000200\tPfam-only domain\tPF00200\t150\t250",
            "PROT1\tIPR000300\tSF domain B\tSSF50003\t300\t400",
            "PROT1\tIPR000400\tSF domain C\tSSF50004\t450\t550",
            # A second signature for the same InterPro entry (the real file
            # lists one row per member signature).
            "PROT1\tIPR000100\tSF domain A\tSM00100\t10\t110",
            "PROT2\tIPR000100\tSF domain A\tSSF50001\t20\t120",
            "PROT2\tIPR000500\tCDD-only domain\tcd00500\t200\t300",
            # SFLD looks like SUPERFAMILY on a sloppy prefix test.
            "PROT2\tIPR000600\tSFLD domain\tSFLDG01135\t400\t500",
        ]
        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))
        return file_path

    def test_domain_key_is_validated(self):
        assert DOMAIN_KEYS == ("interpro", "ssf")
        with pytest.raises(ValueError, match="Unknown domain_key"):
            DomainAnnotationParser(domain_key="pfam")

    def test_default_key_is_interpro(self, mixed_member_db_file):
        """The default keying must be untouched by the new option."""
        parser = DomainAnnotationParser()
        architectures = parser.parse_protein2ipr_file(mixed_member_db_file)

        assert architectures["PROT1"].single_domains == [
            "IPR000100",
            "IPR000100",
            "IPR000200",
            "IPR000300",
            "IPR000400",
        ]
        assert parser.key_to_interpro == {}

    def test_ssf_key_uses_the_signature(self, mixed_member_db_file):
        parser = DomainAnnotationParser(domain_key="ssf")
        architectures = parser.parse_protein2ipr_file(mixed_member_db_file)

        assert architectures["PROT1"].single_domains == [
            "SSF50001",
            "SSF50003",
            "SSF50004",
        ]
        assert set(parser.domain_counts) == {"SSF50001", "SSF50003", "SSF50004"}

    def test_ssf_key_drops_other_member_databases(self, mixed_member_db_file):
        parser = DomainAnnotationParser(domain_key="ssf")
        architectures = parser.parse_protein2ipr_file(mixed_member_db_file)

        every_domain = {
            domain for arch in architectures.values() for domain in arch.single_domains
        }
        assert all(domain.startswith("SSF") for domain in every_domain)
        # SFLD starts with "SF" but is not SUPERFAMILY.
        assert "SFLDG01135" not in every_domain
        # PROT2's only SSF row is SSF50001; its CDD and SFLD rows are gone.
        assert architectures["PROT2"].single_domains == ["SSF50001"]

    def test_dropped_rows_are_not_in_domain_annotations(self, mixed_member_db_file):
        """The filter must apply to the stored annotations, not just the keys.

        ``domain_annotations`` is what the surprise score re-reads to locate a
        feature's matched regions; a non-SSF row surviving there would shift
        every positional lookup.
        """
        parser = DomainAnnotationParser(domain_key="ssf")
        architectures = parser.parse_protein2ipr_file(mixed_member_db_file)

        signatures = [
            ann.signature_id for ann in architectures["PROT1"].domain_annotations
        ]
        assert signatures == ["SSF50001", "SSF50003", "SSF50004"]

    def test_supra_domains_ignore_dropped_rows(self, mixed_member_db_file):
        """The drop-before-sort trap.

        PROT1's architecture is SSF50001 · (Pfam-only) · SSF50003 · SSF50004.
        Keyed by SSF, SSF50001 and SSF50003 *are* adjacent, because the Pfam row
        is not a domain in this universe. If the non-SSF rows were dropped after
        the sort-by-start — or left in the positional list and filtered later —
        the contiguous windows would be computed over a mixed ordering and the
        supra-domains would be wrong (e.g. no ``SSF50001,SSF50003`` pair, or a
        3-mer spanning the discarded Pfam row).
        """
        parser = DomainAnnotationParser(domain_key="ssf", max_supra_domain_length=3)
        architectures = parser.parse_protein2ipr_file(mixed_member_db_file)

        supra = architectures["PROT1"].supra_domains
        assert set(supra) == {
            "SSF50001,SSF50003",
            "SSF50003,SSF50004",
            "SSF50001,SSF50003,SSF50004",
        }
        # Nothing from another member database may leak into a combination.
        assert all("PF" not in feature and "cd" not in feature for feature in supra)

    def test_supra_domains_follow_position_not_file_order(self, tmp_path):
        """Contiguity is by start coordinate, even when the file is unordered."""
        file_path = tmp_path / "unordered.dat.gz"
        data = [
            "PROT1\tIPR000300\tC\tSSF50003\t300\t400",
            "PROT1\tIPR000200\tPfam-only\tPF00200\t150\t250",
            "PROT1\tIPR000100\tA\tSSF50001\t10\t110",
        ]
        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        parser = DomainAnnotationParser(domain_key="ssf", max_supra_domain_length=3)
        architectures = parser.parse_protein2ipr_file(file_path)

        assert architectures["PROT1"].single_domains == ["SSF50001", "SSF50003"]
        assert architectures["PROT1"].supra_domains == ["SSF50001,SSF50003"]

    def test_proteins_without_any_ssf_hit_are_absent(self, tmp_path):
        """Coverage cost of SSF keying: SSF-less proteins leave the universe."""
        file_path = tmp_path / "no_ssf.dat.gz"
        data = [
            "PROT1\tIPR000100\tA\tSSF50001\t10\t110",
            "PROT2\tIPR000200\tPfam-only\tPF00200\t10\t110",
        ]
        with gzip.open(file_path, "wt") as f:
            f.write("\n".join(data))

        parser = DomainAnnotationParser(domain_key="ssf")
        architectures = parser.parse_protein2ipr_file(file_path)

        assert set(architectures) == {"PROT1"}

    def test_interpro_cross_reference(self, mixed_member_db_file):
        """SSF→InterPro is emitted alongside, for joining the two keyings."""
        parser = DomainAnnotationParser(domain_key="ssf")
        parser.parse_protein2ipr_file(mixed_member_db_file)

        assert parser.interpro_for("SSF50001") == "IPR000100"
        assert parser.interpro_for("SSF50001,SSF50003") == "IPR000100,IPR000300"
        assert parser.interpro_for("SSF99999") == "-"

    def test_interpro_cross_reference_is_identity_for_interpro_key(
        self, mixed_member_db_file
    ):
        parser = DomainAnnotationParser()
        parser.parse_protein2ipr_file(mixed_member_db_file)

        assert parser.interpro_for("IPR000100") == "IPR000100"

    def test_get_protein_domain_map_uses_the_key(self, mixed_member_db_file):
        parser = DomainAnnotationParser(domain_key="ssf", max_supra_domain_length=2)
        parser.parse_protein2ipr_file(mixed_member_db_file)

        mapping = parser.get_protein_domain_map()
        assert "SSF50001" in mapping["PROT1"]
        assert "SSF50001,SSF50003" in mapping["PROT1"]
        assert not any(feature.startswith("IPR") for feature in mapping["PROT1"])


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
