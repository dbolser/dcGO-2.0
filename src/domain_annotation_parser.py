"""
Domain Annotation Parser for dcGO Pipeline

This module parses pre-computed domain annotations from InterPro's protein2ipr.dat file
instead of running InterProScan locally. This is the recommended approach for the dcGO
methodology when using UniProt proteins with existing domain annotations.

The protein2ipr.dat file format:
- Tab-separated values
- Fields: UniProt_accession, InterPro_accession, InterPro_name, signature_accession, start_location, end_location
- One line per domain annotation
- Gzipped file (~20GB compressed)

Author: dcGO Pipeline
"""

import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from loguru import logger


@dataclass
class DomainAnnotation:
    """Represents a single domain annotation from InterPro."""

    protein_id: str
    interpro_id: str
    interpro_name: str
    signature_id: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Calculate domain length."""
        return self.end - self.start + 1


@dataclass
class ProteinDomainArchitecture:
    """Represents the complete domain architecture of a protein."""

    protein_id: str
    single_domains: List[str]  # Individual domain IDs
    supra_domains: List[str]   # Contiguous domain combinations
    domain_annotations: List[DomainAnnotation]

    @property
    def all_domains(self) -> List[str]:
        """Get all domains including single and supra-domains."""
        return self.single_domains + self.supra_domains


class DomainAnnotationParser:
    """
    Parser for pre-computed InterPro domain annotations.

    This class handles parsing of protein2ipr.dat files and generates both
    single domain and supra-domain annotations for the dcGO analysis.
    """

    def __init__(
        self,
        max_supra_domain_length: int = 3,
        min_domain_length: int = 10,
        species_filter: Optional[Set[str]] = None
    ):
        """
        Initialize the domain annotation parser.

        Args:
            max_supra_domain_length: Maximum number of domains in a supra-domain
            min_domain_length: Minimum domain length to consider
            species_filter: Set of UniProt accession prefixes to filter by species
                           (e.g., human proteins typically start with specific patterns)
        """
        self.max_supra_domain_length = max_supra_domain_length
        self.min_domain_length = min_domain_length
        self.species_filter = species_filter

        # Storage for parsed annotations
        self.protein_domains: Dict[str, List[DomainAnnotation]] = defaultdict(list)
        self.domain_counts: Dict[str, int] = defaultdict(int)

        logger.info(f"DomainAnnotationParser initialized:")
        logger.info(f"  Max supra-domain length: {max_supra_domain_length}")
        logger.info(f"  Min domain length: {min_domain_length}")
        if species_filter:
            logger.info(f"  Species filter: {len(species_filter)} patterns")

    def parse_protein2ipr_file(
        self,
        protein2ipr_path: Path,
        max_proteins: Optional[int] = None,
        protein_filter: Optional[Set[str]] = None
    ) -> Dict[str, ProteinDomainArchitecture]:
        """
        Parse the protein2ipr.dat file to extract domain annotations.

        Args:
            protein2ipr_path: Path to protein2ipr.dat.gz file
            max_proteins: Maximum number of proteins to process (for testing)
            protein_filter: Set of protein IDs to include (all others ignored for memory efficiency)

        Returns:
            Dictionary mapping protein IDs to their domain architectures
        """
        logger.info(f"Parsing domain annotations from {protein2ipr_path}")
        if protein_filter:
            logger.info(f"  Filtering to {len(protein_filter):,} specific proteins")

        if not protein2ipr_path.exists():
            raise FileNotFoundError(f"protein2ipr file not found: {protein2ipr_path}")

        # Parse annotations
        protein_count = 0
        annotation_count = 0
        filtered_count = 0
        skipped_count = 0

        # Determine if file is gzipped
        open_func = gzip.open if protein2ipr_path.suffix == '.gz' else open

        with open_func(protein2ipr_path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 1000000 == 0:
                    logger.info(f"Processed {line_num:,} lines, {protein_count:,} proteins, "
                              f"{annotation_count:,} annotations, {skipped_count:,} skipped")

                # Skip empty lines
                line = line.strip()
                if not line:
                    continue

                # Parse tab-separated fields
                try:
                    fields = line.split('\t')
                    if len(fields) < 6:
                        logger.debug(f"Skipping malformed line {line_num}: insufficient fields")
                        continue

                    protein_id = fields[0]

                    # Skip proteins not in filter set (for memory efficiency)
                    if protein_filter and protein_id not in protein_filter:
                        skipped_count += 1
                        continue

                    interpro_id = fields[1]
                    interpro_name = fields[2]
                    signature_id = fields[3]
                    start = int(fields[4])
                    end = int(fields[5])

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

                # Apply species filter if provided
                if self.species_filter:
                    if not any(protein_id.startswith(prefix) for prefix in self.species_filter):
                        filtered_count += 1
                        continue

                # Create domain annotation
                annotation = DomainAnnotation(
                    protein_id=protein_id,
                    interpro_id=interpro_id,
                    interpro_name=interpro_name,
                    signature_id=signature_id,
                    start=start,
                    end=end
                )

                # Filter by minimum domain length
                if annotation.length < self.min_domain_length:
                    continue

                # Store annotation
                if protein_id not in self.protein_domains:
                    protein_count += 1

                    # Check max proteins limit
                    if max_proteins and protein_count > max_proteins:
                        logger.info(f"Reached maximum protein limit: {max_proteins}")
                        break

                self.protein_domains[protein_id].append(annotation)
                self.domain_counts[interpro_id] += 1
                annotation_count += 1

        logger.info(f"Parsing complete:")
        logger.info(f"  Total proteins: {protein_count:,}")
        logger.info(f"  Total annotations: {annotation_count:,}")
        logger.info(f"  Unique domains: {len(self.domain_counts):,}")
        if self.species_filter:
            logger.info(f"  Filtered out: {filtered_count:,} annotations")

        # Generate domain architectures
        return self._generate_domain_architectures()

    def _generate_domain_architectures(self) -> Dict[str, ProteinDomainArchitecture]:
        """
        Generate domain architectures including supra-domains.

        Returns:
            Dictionary mapping protein IDs to their complete domain architectures
        """
        logger.info("Generating domain architectures with supra-domains...")

        architectures = {}

        for protein_id, annotations in self.protein_domains.items():
            # Sort annotations by start position
            sorted_annotations = sorted(annotations, key=lambda x: x.start)

            # Extract single domain IDs
            single_domains = [ann.interpro_id for ann in sorted_annotations]

            # Generate supra-domains (contiguous domain combinations)
            supra_domains = self._generate_supra_domains(single_domains)

            # Create architecture
            architecture = ProteinDomainArchitecture(
                protein_id=protein_id,
                single_domains=single_domains,
                supra_domains=supra_domains,
                domain_annotations=sorted_annotations
            )

            architectures[protein_id] = architecture

        # Calculate statistics
        total_supra_domains = sum(len(arch.supra_domains) for arch in architectures.values())
        logger.info(f"Generated {len(architectures):,} domain architectures")
        logger.info(f"  Total supra-domains: {total_supra_domains:,}")

        return architectures

    def _generate_supra_domains(self, domain_ids: List[str]) -> List[str]:
        """
        Generate supra-domains from a list of domain IDs.

        Supra-domains are contiguous combinations of domains, representing
        domain architectures that occur together in proteins.

        Args:
            domain_ids: List of domain IDs in positional order

        Returns:
            List of supra-domain strings (comma-separated domain IDs)
        """
        supra_domains = []

        # Generate all contiguous combinations up to max length
        for length in range(2, min(len(domain_ids) + 1, self.max_supra_domain_length + 1)):
            for i in range(len(domain_ids) - length + 1):
                # Create supra-domain from contiguous domains
                supra_domain = ','.join(domain_ids[i:i+length])
                supra_domains.append(supra_domain)

        return supra_domains

    def get_protein_domain_map(self) -> Dict[str, List[str]]:
        """
        Get a simple mapping from protein IDs to all their domains.

        Returns:
            Dictionary mapping protein IDs to lists of domain IDs (including supra-domains)
        """
        protein_domain_map = {}

        for protein_id, annotations in self.protein_domains.items():
            single_domains = [ann.interpro_id for ann in annotations]
            supra_domains = self._generate_supra_domains(single_domains)
            protein_domain_map[protein_id] = single_domains + supra_domains

        return protein_domain_map

    def get_domain_statistics(self) -> Dict[str, Dict[str, int]]:
        """
        Get statistics about parsed domains.

        Returns:
            Dictionary with domain statistics including counts and coverage
        """
        return {
            'total_proteins': len(self.protein_domains),
            'total_unique_domains': len(self.domain_counts),
            'domain_counts': dict(sorted(
                self.domain_counts.items(),
                key=lambda x: x[1],
                reverse=True
            )[:100])  # Top 100 domains
        }

    def filter_by_human_uniprot(self) -> None:
        """
        Filter annotations to keep only human UniProt proteins.

        Human reviewed (Swiss-Prot) proteins can be identified by their
        accession patterns, though the most reliable way is to use the
        human-specific GOA file which already contains only human proteins.
        """
        # This is a simplified filter - in practice, we'll rely on the
        # intersection with the human GOA file
        logger.info("Note: For human-only analysis, use human-specific GOA file")
        logger.info("The pipeline will automatically intersect domain and GO annotations")


def parse_human_domains(
    protein2ipr_path: Path,
    max_supra_domain_length: int = 3,
    max_proteins: Optional[int] = None
) -> Dict[str, ProteinDomainArchitecture]:
    """
    Convenience function to parse human domain annotations.

    Args:
        protein2ipr_path: Path to protein2ipr.dat.gz file
        max_supra_domain_length: Maximum supra-domain length
        max_proteins: Maximum proteins to process (for testing)

    Returns:
        Dictionary mapping protein IDs to their domain architectures
    """
    parser = DomainAnnotationParser(
        max_supra_domain_length=max_supra_domain_length
    )

    return parser.parse_protein2ipr_file(
        protein2ipr_path,
        max_proteins=max_proteins
    )
