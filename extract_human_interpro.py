#!/usr/bin/env python3
"""
Extract human protein annotations from the full protein2ipr.dat.gz file.
This creates a much smaller human-specific file for faster processing.
"""

import gzip
from pathlib import Path
from loguru import logger
import sys

logger.remove()
logger.add(sys.stderr, level="INFO")

def extract_human_proteins(
    protein_set_file: Path,
    interpro_file: Path,
    output_file: Path
):
    """
    Extract lines from protein2ipr.dat for specific proteins.

    Args:
        protein_set_file: Text file with one protein ID per line
        interpro_file: Full protein2ipr.dat.gz file
        output_file: Output file (will be gzipped)
    """
    # Load protein set
    logger.info(f"Loading protein IDs from {protein_set_file}")
    with open(protein_set_file) as f:
        protein_ids = set(line.strip() for line in f if line.strip())
    logger.info(f"  Loaded {len(protein_ids):,} protein IDs")

    # Extract matching lines
    logger.info(f"Extracting annotations from {interpro_file}")
    logger.info(f"  This will take ~10 minutes to scan the 20GB file")

    matched_lines = 0
    total_lines = 0

    with gzip.open(interpro_file, 'rt') as fin, \
         gzip.open(output_file, 'wt') as fout:
        for line in fin:
            total_lines += 1

            if total_lines % 10000000 == 0:
                logger.info(f"  Processed {total_lines:,} lines, found {matched_lines:,} matches")

            # Quick check: does the line start with any of our protein IDs?
            protein_id = line.split('\t')[0] if '\t' in line else ''

            if protein_id in protein_ids:
                fout.write(line)
                matched_lines += 1

    logger.info(f"✓ Extraction complete!")
    logger.info(f"  Total lines scanned: {total_lines:,}")
    logger.info(f"  Matching lines found: {matched_lines:,}")
    logger.info(f"  Output file: {output_file}")
    logger.info(f"  Output size: {output_file.stat().st_size / 1e6:.1f} MB")

def main():
    # Parse GOA to get human protein IDs
    from src.goa_parser import parse_goa_human

    logger.info("Step 1: Parsing GOA to get human protein IDs...")
    goa_file = Path("data/raw/goa_annotations/goa_human.gaf.gz")
    protein_go_map = parse_goa_human(
        goa_file,
        evidence_filter='manual',
        aspects={'P', 'F', 'C'}
    )

    # Write protein IDs to temp file
    protein_list_file = Path("data/interim/human_proteins.txt")
    protein_list_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing {len(protein_go_map):,} protein IDs to {protein_list_file}")
    with open(protein_list_file, 'w') as f:
        for protein_id in sorted(protein_go_map.keys()):
            f.write(f"{protein_id}\n")

    # Extract human annotations
    logger.info("")
    logger.info("Step 2: Extracting human protein annotations from InterPro...")
    interpro_file = Path("data/raw/interpro_mappings/protein2ipr.dat.gz")
    output_file = Path("data/interim/protein2ipr_human.dat.gz")

    extract_human_proteins(protein_list_file, interpro_file, output_file)

    logger.info("")
    logger.info("=" * 60)
    logger.info("SUCCESS! Human-specific InterPro file created")
    logger.info("=" * 60)
    logger.info(f"You can now use: {output_file}")
    logger.info(f"This file contains only the {len(protein_go_map):,} human proteins from GOA")

if __name__ == "__main__":
    sys.exit(main() or 0)
