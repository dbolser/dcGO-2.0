#!/usr/bin/env python3
"""Characterise the all-species universe before believing anything it produces.

Two of the three traps recorded in TODO.md are properties of the *input*, so
they can be measured before a single Fisher test runs:

* **Annotation-transfer circularity.** How much of the multi-species annotation
  is projected rather than observed, and how much of the projected part cites a
  human protein as its source (GAF column 8, With/From). An association learned
  from annotations that were themselves inferred from human function cannot be
  used as evidence about human function.
* **Taxonomic composition.** How the universe divides across organisms, which
  decides whether "all species" means a broad sweep or a handful of well-curated
  model organisms wearing a broad label.

Writes one tidy TSV per artefact so the numbers can be cited without rerunning.
"""

import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

GAF = Path("data/raw/goa_annotations/goa_allspecies.gaf.gz")
HUMAN_GAF = Path("data/raw/goa_annotations/goa_human.gaf.gz")
OUT = Path("scratch_allspecies/out")

# Codes that assert function was transferred from another protein rather than
# observed in this one. Sequence/orthology/phylogeny based, per the GO docs.
PROJECTED = {"ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA"}
EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}


def human_accessions() -> set:
    """Accessions annotated in the human GOA file — the circularity reference."""
    acc = set()
    with gzip.open(HUMAN_GAF, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.split("\t")
            if len(f) > 1:
                acc.add(f[1])
    return acc


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    human = human_accessions()
    print(f"human reference accessions: {len(human):,}", file=sys.stderr)

    evidence = Counter()
    taxon_annotations = Counter()
    taxon_proteins = defaultdict(set)
    # Circularity: projected annotations on NON-human proteins whose With/From
    # names a human accession. Those are the ones that would feed human function
    # back into a "background" we then use to predict human function.
    projected_nonhuman = 0
    projected_nonhuman_from_human = 0
    aspect = Counter()
    lines = 0

    with gzip.open(GAF, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            lines += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            acc, qualifier, code, with_from = f[1], f[3], f[6], f[7]
            taxon, asp = f[12], f[8]
            if "NOT" in qualifier.upper():
                continue
            evidence[code] += 1
            aspect[asp] += 1
            # taxon field is like "taxon:9606" or "taxon:9606|taxon:1234"
            tax = taxon.split("|")[0].replace("taxon:", "")
            taxon_annotations[tax] += 1
            taxon_proteins[tax].add(acc)
            if code in PROJECTED and acc not in human:
                projected_nonhuman += 1
                if any(part.split(":")[-1] in human for part in with_from.split("|")):
                    projected_nonhuman_from_human += 1

    with (OUT / "evidence_composition.tsv").open("w") as fh:
        fh.write("evidence_code\tclass\tannotations\tfraction\n")
        total = sum(evidence.values())
        for code, n in evidence.most_common():
            cls = (
                "experimental"
                if code in EXPERIMENTAL
                else "projected"
                if code in PROJECTED
                else "other"
            )
            fh.write(f"{code}\t{cls}\t{n}\t{n / total:.6f}\n")

    with (OUT / "taxon_composition.tsv").open("w") as fh:
        fh.write("taxon\tannotations\tproteins\n")
        for tax, n in taxon_annotations.most_common():
            fh.write(f"{tax}\t{n}\t{len(taxon_proteins[tax])}\n")

    total = sum(evidence.values())
    exp = sum(n for c, n in evidence.items() if c in EXPERIMENTAL)
    proj = sum(n for c, n in evidence.items() if c in PROJECTED)
    with (OUT / "universe_summary.tsv").open("w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"gaf_data_lines\t{lines}\n")
        fh.write(f"annotations_kept\t{total}\n")
        fh.write(f"distinct_taxa\t{len(taxon_annotations)}\n")
        fh.write(f"distinct_proteins\t{len(set().union(*taxon_proteins.values()))}\n")
        fh.write(f"experimental_annotations\t{exp}\n")
        fh.write(f"projected_annotations\t{proj}\n")
        fh.write(f"projected_fraction\t{proj / total:.6f}\n")
        fh.write(f"projected_nonhuman\t{projected_nonhuman}\n")
        fh.write(f"projected_nonhuman_from_human\t{projected_nonhuman_from_human}\n")
        fh.write(
            "projected_nonhuman_from_human_fraction\t"
            f"{projected_nonhuman_from_human / max(projected_nonhuman, 1):.6f}\n"
        )
        for asp, n in aspect.most_common():
            fh.write(f"aspect_{asp}\t{n}\n")

    print(f"kept {total:,} annotations over {len(taxon_annotations):,} taxa")
    print(f"projected: {proj:,} ({proj / total:.1%}); experimental: {exp:,}")
    print(
        f"projected non-human citing a human protein: "
        f"{projected_nonhuman_from_human:,} / {projected_nonhuman:,} "
        f"({projected_nonhuman_from_human / max(projected_nonhuman, 1):.1%})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
