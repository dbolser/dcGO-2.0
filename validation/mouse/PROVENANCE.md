# Mouse held-out evaluation — provenance

Everything needed to reproduce `temporal_benchmark_metrics.tsv` and
`temporal_benchmark_permutation_null.tsv` in this directory. Recorded because
the result turns on a choice that is easy to get wrong (see "Release matching").

## Release matching — read this first

GOA release numbers are **per species**. Mouse release 205 is dated
**2023-09-21**, not 2021-04 like human release 205. Using it would compare a
33-month window against human's five-year one and make mouse look weak for a
reason that has nothing to do with the method.

| Snapshot | File | `date-generated` |
| --- | --- | --- |
| human t0 | `goa_human.gaf.205.gz` | 2021-04-21 13:14 |
| **mouse t0** | **`goa_mouse.gaf.191.gz`** | **2021-04-08 09:07** |
| mouse t1 | `goa_mouse.gaf.gz` (current) | 2026-06 |

Mouse 191 is thirteen days from human 205. Verified by reading the header of
each file rather than by assuming the numbering aligns.

## Inputs

| Role | Path | Source |
| --- | --- | --- |
| t0 annotations | `data/raw/goa_archive/goa_mouse.gaf.191.gz` | `https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/MOUSE/goa_mouse.gaf.191.gz` (note: uppercase `MOUSE`; the lowercase path 404s) |
| t1 annotations | `data/raw/goa_annotations/goa_mouse.gaf.gz` | `scripts/download_data.py --species mouse --datasets goa_annotations` |
| domains | `data/interim/protein2ipr_mouse.dat.gz` | `extract_human_interpro.py --species mouse` — 20,990 mouse proteins from the 13 GB `protein2ipr.dat.gz` |
| GO ontology | `data/raw/go_ontology/go-basic.obo` | as for the human runs |

## Commands

```bash
# 1. Inputs
uv run python scripts/download_data.py --species mouse --datasets goa_annotations
uv run python extract_human_interpro.py --species mouse
curl -o data/raw/goa_archive/goa_mouse.gaf.191.gz \
  https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/MOUSE/goa_mouse.gaf.191.gz

# 2. Train on t0. `--species mouse_t0_2021` resolves to the symlinks
#    goa_mouse_t0_2021.gaf.gz -> goa_mouse.gaf.191.gz and
#    protein2ipr_mouse_t0_2021.dat.gz -> protein2ipr_mouse.dat.gz
uv run python run_dcgo_human.py --species mouse_t0_2021 --num-cores 8 \
    --output-dir results_mouse_t0_2021

# 3. Score against t1
uv run python validation/temporal_benchmark.py \
    --t0-gaf data/raw/goa_archive/goa_mouse.gaf.191.gz \
    --t1-gaf data/raw/goa_annotations/goa_mouse.gaf.gz \
    --predictions results_mouse_t0_2021/domain_go_associations_significant.tsv \
    --interpro data/interim/protein2ipr_mouse.dat.gz \
    --min-ic 0 --min-ic 2 --min-ic 4 --min-ic 6 \
    --output-dir validation/mouse
```

Defaults in force, none overridden: `--transfer pscore`, `--evidence-filter
manual`, `--fdr-threshold 0.01`, `--seed 0`, 100 permutations, supra-domains on,
True Path off, no minimum support.

## The trained model

`results_mouse_t0_2021/domain_go_associations_significant.tsv` — 20,732
proteins, 103,229 domain features, 18,588 GO terms, **158,378 associations** at
FDR<0.01 (single-family cutoff 1.21e-06, supra-family 7.79e-07). The table
itself is gitignored for size; its run manifest records the input SHA-256s and
the git commit.

## Code version

Produced on the code after #44–#47: shrinkage removed, single and supra-domain
BH families corrected separately, True Path parental background propagated,
`--min-support` available but unused here. The human numbers this is compared
against were regenerated on the same commit, so the two are directly comparable.
