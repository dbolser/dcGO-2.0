#!/usr/bin/env python3
"""
Extract a species' protein annotations from the full protein2ipr.dat.gz file.

This creates a much smaller species-specific file for faster processing. The
name is kept for backward compatibility, but the script works for any organism
via ``--species`` (default: human), mirroring ``run_dcgo_human.py``.

    uv run python extract_human_interpro.py                 # human (default)
    uv run python extract_human_interpro.py --species mouse # any organism
"""

import argparse
import gzip
from pathlib import Path
from loguru import logger
import sys

logger.remove()
logger.add(sys.stderr, level="INFO")


def extract_species_proteins(
    protein_set_file: Path, interpro_file: Path, output_file: Path
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
    logger.info("  This will take ~10 minutes to scan the 20GB file")

    matched_lines = 0
    total_lines = 0

    with gzip.open(interpro_file, "rt") as fin, gzip.open(output_file, "wt") as fout:
        for line in fin:
            total_lines += 1

            if total_lines % 10000000 == 0:
                logger.info(
                    f"  Processed {total_lines:,} lines, found {matched_lines:,} matches"
                )

            # Quick check: does the line start with any of our protein IDs?
            protein_id = line.split("\t")[0] if "\t" in line else ""

            if protein_id in protein_ids:
                fout.write(line)
                matched_lines += 1

    logger.info("✓ Extraction complete!")
    logger.info(f"  Total lines scanned: {total_lines:,}")
    logger.info(f"  Matching lines found: {matched_lines:,}")
    logger.info(f"  Output file: {output_file}")
    logger.info(f"  Output size: {output_file.stat().st_size / 1e6:.1f} MB")


def main():
    parser = argparse.ArgumentParser(
        description="Extract a species' subset of protein2ipr.dat.gz for a fast dcGO run."
    )
    parser.add_argument(
        "--species",
        default="human",
        help="Species to extract: 'human', 'mouse', etc. Must match the "
        "goa_<species>.gaf.gz file name (default: human)",
    )
    parser.add_argument(
        "--evidence-filter",
        default="manual",
        choices=["all", "manual", "experimental"],
        help="Evidence filter used to select the protein set (default: manual)",
    )
    args = parser.parse_args()

    # Parse GOA to get this species' protein IDs
    from src.goa_parser import parse_goa

    logger.info(f"Step 1: Parsing GOA to get {args.species} protein IDs...")
    goa_file = Path(f"data/raw/goa_annotations/goa_{args.species}.gaf.gz")
    if not goa_file.exists():
        logger.error(f"GOA file not found: {goa_file}")
        logger.error(
            f"Download it first, e.g. "
            f"uv run python scripts/download_data.py --species {args.species}"
        )
        return 1

    protein_go_map = parse_goa(
        goa_file, evidence_filter=args.evidence_filter, aspects={"P", "F", "C"}
    )

    # Write protein IDs to temp file
    protein_list_file = Path(f"data/interim/{args.species}_proteins.txt")
    protein_list_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing {len(protein_go_map):,} protein IDs to {protein_list_file}")
    with open(protein_list_file, "w") as f:
        for protein_id in sorted(protein_go_map.keys()):
            f.write(f"{protein_id}\n")

    # Extract species annotations
    logger.info("")
    logger.info(
        f"Step 2: Extracting {args.species} protein annotations from InterPro..."
    )
    interpro_file = Path("data/raw/interpro_mappings/protein2ipr.dat.gz")
    if not interpro_file.exists():
        logger.error(f"InterPro mappings file not found: {interpro_file}")
        logger.error(
            "Download it first: uv run python scripts/download_data.py "
            "--datasets interpro_mappings"
        )
        return 1
    output_file = Path(f"data/interim/protein2ipr_{args.species}.dat.gz")

    extract_species_proteins(protein_list_file, interpro_file, output_file)

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"SUCCESS! {args.species.title()}-specific InterPro file created")
    logger.info("=" * 60)
    logger.info(f"You can now use: {output_file}")
    logger.info(
        f"This file contains only the {len(protein_go_map):,} "
        f"{args.species} proteins from GOA"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
