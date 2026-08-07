#!/usr/bin/env python3
"""Measure phylogenetic non-independence in a multi-species association set.

Fisher's exact test counts proteins and assumes each is an independent
observation. In a multi-species universe that assumption fails in a specific,
directional way: fifty mammalian orthologs of one protein are nearer to *one*
observation than to fifty, so every cell of the contingency table is inflated
and the p-value is optimistic by an unknown factor.

This does not try to correct the test. It reports, for every significant
association, the pooled support next to two collapsed supports, so the
inflation is visible rather than assumed away:

* ``n_proteins``  — the support the p-value was actually computed from
* ``n_uniref50``  — distinct UniRef50 clusters among those proteins, the
                    ortholog-group proxy (50% identity groups orthologs and
                    close paralogs across all species)
* ``n_uniref90``  — a tighter clustering, closer to "same protein, near species"
* ``n_taxa``      — distinct organisms

The ratio ``n_proteins / n_uniref50`` is the per-association inflation factor.
An association at ratio 1.0 is carried by genuinely distinct sequences; one at
ratio 30 is one observation counted thirty times.

Usage:
    uv run python scratch_allspecies/05_non_independence.py \
        --associations results_allspecies_manual/domain_go_associations_significant.tsv \
        --evidence-filter manual \
        --out scratch_allspecies/out/non_independence_manual.tsv
"""

import argparse
import csv
import gzip
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from loguru import logger  # noqa: E402

# Take the domain-parsing constants from the pipeline entry point, not a copy:
# a divergence here would silently measure a different domain universe than the
# run whose associations we are auditing.
from run_dcgo_human import MAX_SUPRA_DOMAIN_LENGTH, MIN_DOMAIN_LENGTH  # noqa: E402
from src.domain_annotation_parser import DomainAnnotationParser  # noqa: E402
from src.goa_parser import parse_goa  # noqa: E402

UNIREF = Path("data/interim/uniref_taxon.tsv.gz")


def load_uniref(accessions: set) -> dict:
    """accession -> (uniref90, uniref50, taxon), restricted to what we need."""
    mapping = {}
    with gzip.open(UNIREF, "rt") as fh:
        for line in fh:
            acc, sep, rest = line.partition("\t")
            if acc not in accessions:
                continue
            parts = rest.rstrip("\n").split("\t")
            u90 = parts[0] or ""
            u50 = parts[1] or ""
            taxon = parts[2] or ""
            mapping[acc] = (u90, u50, taxon)
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--associations", type=Path, required=True)
    ap.add_argument("--species", default="allspecies")
    ap.add_argument("--evidence-filter", default="manual")
    ap.add_argument("--domain-key", default="interpro")
    ap.add_argument("--supra-domains", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    gaf = Path(f"data/raw/goa_annotations/goa_{args.species}.gaf.gz")
    ipr = Path(f"data/interim/protein2ipr_{args.species}.dat.gz")

    logger.info("Parsing annotations...")
    protein_terms = parse_goa(gaf, evidence_filter=args.evidence_filter)
    logger.info("Parsing domains...")
    arch = DomainAnnotationParser(
        max_supra_domain_length=MAX_SUPRA_DOMAIN_LENGTH,
        min_domain_length=MIN_DOMAIN_LENGTH,
        domain_key=args.domain_key,
    ).parse_protein2ipr_file(ipr)

    universe = set(protein_terms) & set(arch)
    logger.info(f"Universe: {len(universe):,} proteins")

    # Inverted indices over the same universe the run used.
    domain_proteins = defaultdict(set)
    term_proteins = defaultdict(set)
    for protein in universe:
        a = arch[protein]
        domains = set(a.single_domains)
        if args.supra_domains:
            domains |= set(a.supra_domains)
        for d in domains:
            domain_proteins[d].add(protein)
        for t in protein_terms[protein]:
            term_proteins[t].add(protein)

    logger.info("Loading UniRef clusters and taxa...")
    uniref = load_uniref(universe)
    logger.info(f"  mapped {len(uniref):,} / {len(universe):,} universe proteins")

    rows = []
    unmapped_associations = 0
    with args.associations.open() as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            domain, term = rec["domain"], rec["go_term"]
            both = domain_proteins.get(domain, set()) & term_proteins.get(term, set())
            if not both:
                continue
            u90 = {uniref[p][0] for p in both if p in uniref and uniref[p][0]}
            u50 = {uniref[p][1] for p in both if p in uniref and uniref[p][1]}
            taxa = {uniref[p][2] for p in both if p in uniref and uniref[p][2]}
            if not u50:
                unmapped_associations += 1
                continue
            rows.append(
                {
                    "domain": domain,
                    "term": term,
                    "n_proteins": len(both),
                    "n_uniref90": len(u90),
                    "n_uniref50": len(u50),
                    "n_taxa": len(taxa),
                    "inflation_uniref50": round(len(both) / len(u50), 4),
                    "adj_p_value": rec["adj_p_value"],
                }
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    ratios = [r["inflation_uniref50"] for r in rows]
    pooled = sum(r["n_proteins"] for r in rows)
    collapsed = sum(r["n_uniref50"] for r in rows)
    logger.info(f"Wrote {args.out} ({len(rows):,} associations)")
    logger.info(f"  associations with no UniRef50 mapping: {unmapped_associations:,}")
    logger.info(f"  total pooled support     : {pooled:,}")
    logger.info(f"  total UniRef50-collapsed : {collapsed:,}")
    logger.info(f"  overall inflation        : {pooled / collapsed:.2f}x")
    logger.info(f"  median per-association   : {statistics.median(ratios):.2f}x")
    logger.info(
        f"  associations whose support is a single UniRef50 cluster: "
        f"{sum(1 for r in rows if r['n_uniref50'] == 1):,} "
        f"({sum(1 for r in rows if r['n_uniref50'] == 1) / len(rows):.1%})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
