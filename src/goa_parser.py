"""
GOA (Gene Ontology Annotation) file parser with evidence code filtering.

This module parses GAF (Gene Association Format) files and provides
filtering based on evidence codes to control annotation quality.

Evidence Code Categories:
- Experimental (high confidence): EXP, IDA, IPI, IMP, IGI, IEP
- Computational analysis: ISS, ISO, ISA, ISM, IGC, IBA, IBD, IKR, IRD, RCA
- Author statements: TAS, NAS
- Curator statements: IC, ND
- Electronic annotation (automated): IEA

Author: dcGO Pipeline
"""

import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set

from loguru import logger


# Standard evidence code sets for filtering
EXPERIMENTAL_EVIDENCE = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}
COMPUTATIONAL_EVIDENCE = {
    "ISS",
    "ISO",
    "ISA",
    "ISM",
    "IGC",
    "IBA",
    "IBD",
    "IKR",
    "IRD",
    "RCA",
}
AUTHOR_EVIDENCE = {"TAS", "NAS"}
CURATOR_EVIDENCE = {"IC", "ND"}
ELECTRONIC_EVIDENCE = {"IEA"}

# Common evidence code presets
MANUAL_EVIDENCE = (
    EXPERIMENTAL_EVIDENCE | COMPUTATIONAL_EVIDENCE | AUTHOR_EVIDENCE | CURATOR_EVIDENCE
)
NON_IEA_EVIDENCE = (
    EXPERIMENTAL_EVIDENCE | COMPUTATIONAL_EVIDENCE | AUTHOR_EVIDENCE | CURATOR_EVIDENCE
)
ALL_EVIDENCE = (
    EXPERIMENTAL_EVIDENCE
    | COMPUTATIONAL_EVIDENCE
    | AUTHOR_EVIDENCE
    | CURATOR_EVIDENCE
    | ELECTRONIC_EVIDENCE
)


@dataclass
class GOAAnnotation:
    """Represents a single GO annotation from a GAF file."""

    protein_id: str
    go_term: str
    evidence_code: str
    aspect: (
        str  # P (biological process), F (molecular function), C (cellular component)
    )
    db_reference: str
    taxon_id: str

    def __post_init__(self):
        """Validate GO term format."""
        if not self.go_term.startswith("GO:"):
            raise ValueError(f"Invalid GO term format: {self.go_term}")


class GOAParser:
    """
    Parser for GOA (Gene Ontology Annotation) files in GAF format.

    Supports evidence code filtering to control annotation quality and
    aspect filtering to focus on specific GO categories.
    """

    def __init__(
        self,
        evidence_codes: Optional[Set[str]] = None,
        aspects: Optional[Set[str]] = None,
        exclude_qualifiers: bool = True,
    ):
        """
        Initialize GOA parser with filtering options.

        Args:
            evidence_codes: Set of evidence codes to include (None = all)
                           Use EXPERIMENTAL_EVIDENCE, MANUAL_EVIDENCE, etc.
            aspects: Set of GO aspects to include (None = all)
                    'P' = biological process
                    'F' = molecular function
                    'C' = cellular component
            exclude_qualifiers: Exclude annotations with qualifiers like 'NOT'
        """
        self.evidence_codes = evidence_codes
        self.aspects = aspects or {"P", "F", "C"}
        self.exclude_qualifiers = exclude_qualifiers

        self.annotations: List[GOAAnnotation] = []
        self.protein_go_map: Dict[str, Set[str]] = defaultdict(set)
        self.statistics: Dict[str, int] = defaultdict(int)

        logger.info("GOAParser initialized:")
        if evidence_codes:
            logger.info(f"  Evidence codes: {sorted(evidence_codes)}")
        else:
            logger.info("  Evidence codes: ALL")
        logger.info(f"  Aspects: {sorted(self.aspects)}")
        logger.info(f"  Exclude qualifiers: {exclude_qualifiers}")

    def parse_gaf_file(self, gaf_path: Path) -> Dict[str, Set[str]]:
        """
        Parse a GAF file and return protein-GO mappings.

        Args:
            gaf_path: Path to GAF file (can be gzipped)

        Returns:
            Dictionary mapping protein IDs to sets of GO terms
        """
        logger.info(f"Parsing GAF file: {gaf_path}")

        if not gaf_path.exists():
            raise FileNotFoundError(f"GAF file not found: {gaf_path}")

        # Determine if file is gzipped
        open_func = gzip.open if gaf_path.suffix == ".gz" else open

        line_count = 0
        with open_func(gaf_path, "rt") as f:
            for line in f:
                line_count += 1

                # Skip comments and header lines
                if line.startswith("!"):
                    continue

                # Parse GAF 2.2 format
                try:
                    fields = line.strip().split("\t")
                    if len(fields) < 15:
                        self.statistics["malformed_lines"] += 1
                        continue

                    # Extract relevant fields (GAF 2.2 format)
                    db = fields[0]
                    db_object_id = fields[1]
                    qualifier = fields[3]
                    go_id = fields[4]
                    db_reference = fields[5]
                    evidence_code = fields[6]
                    taxon = fields[12]
                    aspect = fields[8]

                    # Create protein ID (typically UniProt accession)
                    protein_id = db_object_id

                    # Apply filters
                    # 1. Exclude qualifiers like 'NOT' if requested
                    if self.exclude_qualifiers and qualifier:
                        if "NOT" in qualifier.upper():
                            self.statistics["excluded_qualifiers"] += 1
                            continue

                    # 2. Filter by evidence code
                    if self.evidence_codes and evidence_code not in self.evidence_codes:
                        self.statistics["filtered_evidence"] += 1
                        continue

                    # 3. Filter by aspect
                    if aspect not in self.aspects:
                        self.statistics["filtered_aspect"] += 1
                        continue

                    # Create annotation
                    annotation = GOAAnnotation(
                        protein_id=protein_id,
                        go_term=go_id,
                        evidence_code=evidence_code,
                        aspect=aspect,
                        db_reference=db_reference,
                        taxon_id=taxon,
                    )

                    self.annotations.append(annotation)
                    self.protein_go_map[protein_id].add(go_id)
                    self.statistics["accepted_annotations"] += 1
                    self.statistics[f"evidence_{evidence_code}"] += 1
                    self.statistics[f"aspect_{aspect}"] += 1

                except (IndexError, ValueError) as e:
                    self.statistics["parse_errors"] += 1
                    logger.debug(f"Error parsing line {line_count}: {e}")
                    continue

        # Log statistics
        logger.info("GAF parsing complete:")
        logger.info(f"  Total lines: {line_count:,}")
        logger.info(
            f"  Accepted annotations: {self.statistics['accepted_annotations']:,}"
        )
        logger.info(f"  Unique proteins: {len(self.protein_go_map):,}")
        logger.info(f"  Filtered by evidence: {self.statistics['filtered_evidence']:,}")
        logger.info(f"  Filtered by aspect: {self.statistics['filtered_aspect']:,}")
        logger.info(
            f"  Excluded qualifiers: {self.statistics['excluded_qualifiers']:,}"
        )

        # Show evidence code distribution
        evidence_stats = {
            k: v for k, v in self.statistics.items() if k.startswith("evidence_")
        }
        if evidence_stats:
            logger.info("  Evidence code distribution:")
            for code, count in sorted(
                evidence_stats.items(), key=lambda x: x[1], reverse=True
            )[:10]:
                logger.info(f"    {code.replace('evidence_', '')}: {count:,}")

        return dict(self.protein_go_map)

    def get_annotations_by_protein(self, protein_id: str) -> List[GOAAnnotation]:
        """Get all annotations for a specific protein."""
        return [ann for ann in self.annotations if ann.protein_id == protein_id]

    def get_proteins_with_go_term(self, go_term: str) -> Set[str]:
        """Get all proteins annotated with a specific GO term."""
        return {ann.protein_id for ann in self.annotations if ann.go_term == go_term}


def parse_goa_human(
    gaf_path: Path, evidence_filter: str = "manual", aspects: Optional[Set[str]] = None
) -> Dict[str, Set[str]]:
    """
    Convenience function to parse human GOA file with common presets.

    Args:
        gaf_path: Path to goa_human.gaf.gz file
        evidence_filter: Evidence filter preset:
            'all' - All evidence codes (including IEA)
            'manual' - Exclude electronic (IEA) annotations
            'experimental' - Only experimental evidence
        aspects: GO aspects to include (default: all)

    Returns:
        Dictionary mapping protein IDs to GO term sets
    """
    # Select evidence codes based on preset
    evidence_codes = {
        "all": None,  # No filtering
        "manual": NON_IEA_EVIDENCE,
        "experimental": EXPERIMENTAL_EVIDENCE,
    }.get(evidence_filter, NON_IEA_EVIDENCE)

    parser = GOAParser(
        evidence_codes=evidence_codes, aspects=aspects or {"P", "F", "C"}
    )

    return parser.parse_gaf_file(gaf_path)
