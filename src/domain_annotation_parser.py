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

A "domain" can be keyed either by the integrated InterPro entry (field 2, the
default) or by a member-database signature (field 4) via ``domain_key``; the
``ssf`` key reproduces the SCOP superfamily universe the published dcGO used.

Author: dcGO Pipeline
"""

import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set

from loguru import logger

# Which column of ``protein2ipr`` names a "domain".
#
# ``interpro`` (the default, and the only behaviour before this option existed)
# keys on the integrated InterPro entry in column 2. ``ssf`` keys on the
# SUPERFAMILY member-database signature in column 4 — ``SSFnnnnn``, whose numeric
# part is the SCOP sunid. That is the domain universe the published dcGO (Fang &
# Gough 2013) used, so it is what makes VALIDATION_PLAN §3 an apples-to-apples
# comparison rather than a cross-namespace one.
#
# Note the prefix test must not confuse SUPERFAMILY with the Structure-Function
# Linkage Database, whose signatures are ``SFLDF``/``SFLDS``/``SFLDG`` — they
# start with ``SF`` but not with ``SSF``, so an exact ``SSF`` prefix is safe.
SIGNATURE_PREFIXES: Dict[str, str] = {"ssf": "SSF"}

DOMAIN_KEYS = ("interpro", *sorted(SIGNATURE_PREFIXES))


def superfamily_sunid(signature_id: str) -> Optional[int]:
    """SCOP sunid behind a SUPERFAMILY signature (``SSF53649`` → ``53649``).

    InterPro's SUPERFAMILY member database is built on the SCOP 1.75 HMM library
    and its accessions are the SCOP sunid with an ``SSF`` prefix, which is what
    lets our domains be joined directly to the published dcGO tables (keyed by
    bare sunid). Returns ``None`` for anything that is not such an accession, so
    callers can filter rather than guess.
    """
    if not signature_id.startswith("SSF"):
        return None
    suffix = signature_id[3:]
    if not suffix.isdigit():
        return None
    return int(suffix)


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
    supra_domains: List[str]  # Contiguous domain combinations
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
        species_filter: Optional[Set[str]] = None,
        domain_key: str = "interpro",
    ):
        """
        Initialize the domain annotation parser.

        Args:
            max_supra_domain_length: Maximum number of domains in a supra-domain
            min_domain_length: Minimum domain length to consider
            species_filter: Set of UniProt accession prefixes to filter by species
                           (e.g., human proteins typically start with specific patterns)
            domain_key: Which ``protein2ipr`` column defines a domain —
                        ``interpro`` (column 2, the integrated entry; the
                        default and historical behaviour) or ``ssf`` (column 4's
                        SUPERFAMILY/SCOP signature). See ``DOMAIN_KEYS``.
        """
        if domain_key not in DOMAIN_KEYS:
            raise ValueError(
                f"Unknown domain_key {domain_key!r}; expected one of "
                + ", ".join(DOMAIN_KEYS)
            )

        self.max_supra_domain_length = max_supra_domain_length
        self.min_domain_length = min_domain_length
        self.species_filter = species_filter
        self.domain_key = domain_key
        self.signature_prefix = SIGNATURE_PREFIXES.get(domain_key)

        # Storage for parsed annotations
        self.protein_domains: Dict[str, List[DomainAnnotation]] = defaultdict(list)
        self.domain_counts: Dict[str, int] = defaultdict(int)
        # Signature → InterPro entry, recorded while parsing so a signature-keyed
        # run can still be cross-referenced against an InterPro-keyed one without
        # a second pass over the 20 GB source. Empty for ``interpro`` keying.
        self.key_to_interpro: Dict[str, Set[str]] = defaultdict(set)

        logger.info("DomainAnnotationParser initialized:")
        logger.info(f"  Domain key: {domain_key}")
        logger.info(f"  Max supra-domain length: {max_supra_domain_length}")
        logger.info(f"  Min domain length: {min_domain_length}")
        if species_filter:
            logger.info(f"  Species filter: {len(species_filter)} patterns")

    def domain_key_of(self, annotation: DomainAnnotation) -> str:
        """The domain identifier this parser keys ``annotation`` by."""
        if self.signature_prefix is None:
            return annotation.interpro_id
        return annotation.signature_id

    def interpro_for(self, domain_id: str) -> str:
        """InterPro entry (or comma-joined entries, for a supra-domain) for a key.

        For ``interpro`` keying this is the identity. For ``ssf`` keying it is
        the entry InterPro integrates that signature into, which over human data
        is a 1:1 bijection — so it is a free cross-reference column rather than a
        lossy mapping. Unmappable parts are reported as ``-``; genuinely
        ambiguous ones are joined with ``|``.
        """
        if self.signature_prefix is None:
            return domain_id
        parts = []
        for part in domain_id.split(","):
            entries = self.key_to_interpro.get(part)
            parts.append("|".join(sorted(entries)) if entries else "-")
        return ",".join(parts)

    def parse_protein2ipr_file(
        self,
        protein2ipr_path: Path,
        max_proteins: Optional[int] = None,
        protein_filter: Optional[Set[str]] = None,
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
        open_func = gzip.open if protein2ipr_path.suffix == ".gz" else open

        with open_func(protein2ipr_path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 1000000 == 0:
                    logger.info(
                        f"Processed {line_num:,} lines, {protein_count:,} proteins, "
                        f"{annotation_count:,} annotations, {skipped_count:,} skipped"
                    )

                # Skip empty lines
                line = line.strip()
                if not line:
                    continue

                # Parse tab-separated fields
                try:
                    fields = line.split("\t")
                    if len(fields) < 6:
                        logger.debug(
                            f"Skipping malformed line {line_num}: insufficient fields"
                        )
                        continue

                    protein_id = fields[0]

                    # Skip proteins not in filter set (for memory efficiency)
                    if protein_filter and protein_id not in protein_filter:
                        skipped_count += 1
                        continue

                    interpro_id = fields[1]
                    interpro_name = fields[2]
                    signature_id = fields[3]

                    # Signature-keyed runs drop the other member databases HERE,
                    # at parse time — before the row is appended to
                    # self.protein_domains and therefore before
                    # _generate_domain_architectures sorts by start position.
                    # Filtering later would leave the discarded rows interleaved
                    # in the positional ordering, so _generate_supra_domains
                    # would treat domains separated by a dropped row as
                    # "contiguous" and emit combinations that do not exist in the
                    # architecture. See tests/unit/test_domain_annotation_parser.py
                    # ::TestSuperfamilyDomainKey::test_supra_domains_ignore_dropped_rows.
                    if (
                        self.signature_prefix is not None
                        and not signature_id.startswith(self.signature_prefix)
                    ):
                        skipped_count += 1
                        continue

                    start = int(fields[4])
                    end = int(fields[5])

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

                # Apply species filter if provided
                if self.species_filter:
                    if not any(
                        protein_id.startswith(prefix) for prefix in self.species_filter
                    ):
                        filtered_count += 1
                        continue

                # Create domain annotation
                annotation = DomainAnnotation(
                    protein_id=protein_id,
                    interpro_id=interpro_id,
                    interpro_name=interpro_name,
                    signature_id=signature_id,
                    start=start,
                    end=end,
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
                self.domain_counts[self.domain_key_of(annotation)] += 1
                if self.signature_prefix is not None:
                    self.key_to_interpro[signature_id].add(interpro_id)
                annotation_count += 1

        logger.info("Parsing complete:")
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
            single_domains = [self.domain_key_of(ann) for ann in sorted_annotations]

            # Generate supra-domains (contiguous domain combinations)
            supra_domains = self._generate_supra_domains(single_domains)

            # Create architecture
            architecture = ProteinDomainArchitecture(
                protein_id=protein_id,
                single_domains=single_domains,
                supra_domains=supra_domains,
                domain_annotations=sorted_annotations,
            )

            architectures[protein_id] = architecture

        # Calculate statistics
        total_supra_domains = sum(
            len(arch.supra_domains) for arch in architectures.values()
        )
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
        for length in range(
            2, min(len(domain_ids) + 1, self.max_supra_domain_length + 1)
        ):
            for i in range(len(domain_ids) - length + 1):
                # Create supra-domain from contiguous domains
                supra_domain = ",".join(domain_ids[i : i + length])
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
            single_domains = [self.domain_key_of(ann) for ann in annotations]
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
            "total_proteins": len(self.protein_domains),
            "total_unique_domains": len(self.domain_counts),
            "domain_counts": dict(
                sorted(self.domain_counts.items(), key=lambda x: x[1], reverse=True)[
                    :100
                ]
            ),  # Top 100 domains
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
        logger.info(
            "The pipeline will automatically intersect domain and GO annotations"
        )


def parse_human_domains(
    protein2ipr_path: Path,
    max_supra_domain_length: int = 3,
    max_proteins: Optional[int] = None,
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
    parser = DomainAnnotationParser(max_supra_domain_length=max_supra_domain_length)

    return parser.parse_protein2ipr_file(protein2ipr_path, max_proteins=max_proteins)
