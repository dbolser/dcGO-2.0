#!/usr/bin/env python3
"""Build ``data/interim/protein2ipr_<species>.dat.gz`` from an accession set.

``extract_human_interpro.py`` selects a species' proteins *via its GOA file*,
which the model-organism phenotype layers (``--ontology mp / wbphenotype /
zfa / fbcv / fbbt``) do not need and may not have. This script selects them
via an accession set instead — by default the first column of a UniProt
per-organism idmapping file (``data/raw/uniprot_idmapping/``), which lists
every Swiss-Prot *and* TrEMBL accession of the taxon.

By default it filters the existing all-species extract
(``data/interim/protein2ipr_allspecies.dat.gz``, seconds) rather than the raw
~20 GB ``protein2ipr.dat.gz``. The caveat is inherited honestly: the
all-species extract only contains proteins that were in the all-species GOA
universe, so accessions outside it yield no domains. Measure what fraction of
your *annotated* accessions the output covers (the run logs it as the
domain-universe intersection); if it is poor, re-run with ``--full-scan`` to
stream the raw file (~10 min).

    uv run python scripts/extract_species_interpro.py --species worm \
        --idmapping data/raw/uniprot_idmapping/CAEEL_6239_idmapping.dat.gz

    uv run python scripts/extract_species_interpro.py --species fly \
        --idmapping data/raw/uniprot_idmapping/DROME_7227_idmapping.dat.gz \
        --full-scan
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

from loguru import logger

# Make repo-root modules importable when run as a plain script.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from extract_human_interpro import extract_species_proteins  # noqa: E402

logger.remove()
logger.add(sys.stderr, level="INFO")

ALLSPECIES_EXTRACT = Path("data/interim/protein2ipr_allspecies.dat.gz")
RAW_PROTEIN2IPR = Path("data/raw/interpro_mappings/protein2ipr.dat.gz")


def accessions_from_idmapping(path: Path) -> set[str]:
    """Distinct accessions (first column) of a per-organism idmapping file."""
    open_func = gzip.open if path.suffix == ".gz" else open
    accessions: set[str] = set()
    with open_func(path, "rt") as handle:
        for line in handle:
            accession = line.split("\t", 1)[0].strip()
            if accession:
                accessions.add(accession)
    return accessions


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract a species' protein2ipr subset from an accession set."
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species name used in the output filename, matching the "
        "run_dcgo_human.py --species value (e.g. worm, zebrafish, fly)",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--idmapping",
        type=Path,
        help="UniProt per-organism idmapping file; its first column defines "
        "the species' accession set",
    )
    group.add_argument(
        "--accession-list",
        type=Path,
        help="Plain text file with one UniProt accession per line",
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help=f"protein2ipr file to filter (default: {ALLSPECIES_EXTRACT}; "
        "see --full-scan)",
    )
    parser.add_argument(
        "--full-scan",
        action="store_true",
        help=f"Filter the raw {RAW_PROTEIN2IPR} instead of the all-species "
        "extract — slower (~10 min) but not limited to the all-species GOA "
        "universe",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: data/interim/protein2ipr_<species>.dat.gz)",
    )
    args = parser.parse_args()

    source = args.source or (RAW_PROTEIN2IPR if args.full_scan else ALLSPECIES_EXTRACT)
    if not source.exists():
        logger.error(f"Source file not found: {source}")
        return 1

    if args.idmapping:
        logger.info(f"Reading accession set from idmapping file {args.idmapping}")
        accessions = accessions_from_idmapping(args.idmapping)
    else:
        logger.info(f"Reading accession list {args.accession_list}")
        accessions = {
            line.strip() for line in open(args.accession_list) if line.strip()
        }
    logger.info(f"  {len(accessions):,} accessions define the {args.species} universe")

    accession_file = Path(f"data/interim/{args.species}_proteins.txt")
    accession_file.parent.mkdir(parents=True, exist_ok=True)
    with open(accession_file, "w") as handle:
        for accession in sorted(accessions):
            handle.write(f"{accession}\n")
    logger.info(f"  Accession list written to {accession_file}")

    output = args.output or Path(f"data/interim/protein2ipr_{args.species}.dat.gz")
    extract_species_proteins(accession_file, source, output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
