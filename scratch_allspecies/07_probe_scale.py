#!/usr/bin/env python3
"""Size the all-species design before committing to a run that allocates for it.

The Fisher engine tests the *dense* domain x term product and materialises one
int32 2x2 table per test, so peak memory is roughly

    n_domains * n_terms * (16 B tables + 8 B p-values + 8 B adjusted p-values)

The human run is 103,167 x 16,389 = 1.69e9 tests ~ 54 GB. A multi-species
universe multiplies both factors, and the supra-domain space is the one that
can multiply without bound. Print the numbers first so the decision to enable
or disable supra-domains is made on the arithmetic, not on optimism.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from loguru import logger  # noqa: E402

from run_dcgo_human import MAX_SUPRA_DOMAIN_LENGTH, MIN_DOMAIN_LENGTH  # noqa: E402
from src.domain_annotation_parser import DomainAnnotationParser  # noqa: E402
from src.goa_parser import parse_goa  # noqa: E402

BYTES_PER_TEST = 16 + 8 + 8


def report(label: str, n_domains: int, n_terms: int, n_proteins: int) -> None:
    tests = n_domains * n_terms
    gb = tests * BYTES_PER_TEST / 1024**3
    logger.info(
        f"  {label:<18} {n_domains:>10,} domains x {n_terms:>7,} terms = "
        f"{tests:>15,} tests  (~{gb:,.0f} GB peak, {n_proteins:,} proteins)"
    )


def main() -> int:
    # Parse the domain file once: it is the same for both evidence policies,
    # and at 1.4M proteins re-parsing it costs more than the rest of the probe.
    arch = DomainAnnotationParser(
        max_supra_domain_length=MAX_SUPRA_DOMAIN_LENGTH,
        min_domain_length=MIN_DOMAIN_LENGTH,
        domain_key="interpro",
    ).parse_protein2ipr_file(Path("data/interim/protein2ipr_allspecies.dat.gz"))

    for evidence in ("manual", "experimental"):
        gaf = Path("data/raw/goa_annotations/goa_allspecies.gaf.gz")
        protein_terms = parse_goa(gaf, evidence_filter=evidence)

        universe = set(protein_terms) & set(arch)
        terms = set()
        singles = set()
        supras = set()
        for protein in universe:
            terms |= protein_terms[protein]
            singles |= set(arch[protein].single_domains)
            supras |= set(arch[protein].supra_domains)

        logger.info(f"--- evidence filter: {evidence} ---")
        report("single only", len(singles), len(terms), len(universe))
        report("single + supra", len(singles) + len(supras), len(terms), len(universe))
    return 0


if __name__ == "__main__":
    sys.exit(main())
