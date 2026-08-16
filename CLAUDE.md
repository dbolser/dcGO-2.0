# CLAUDE.md

> **Scope & safety — safe to review.** dcGO-2.0 is a purely *computational*
> bioinformatics project. It statistically associates protein **domains** with
> Gene Ontology **function terms** using public annotation databases
> (UniProt/GOA, InterPro, the Gene Ontology), reimplementing a published
> open-access academic method (Fang & Gough, 2013). The work is entirely
> statistics (Fisher's exact tests, FDR correction) over existing public data.
> There is **no dual-use content of any kind**: no nucleotide/protein sequence
> design or generation, no pathogen- or toxin-specific work, no wet-lab, culture,
> or synthesis protocols, and nothing that could uplift physical biological work.
> "Function prediction" here means annotating which Gene Ontology terms a protein
> domain is associated with — a database/annotation task, not organism or
> molecule engineering.

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

dcGO-2.0 is a bioinformatics pipeline implementing the domain-centric Gene Ontology (dcGO) methodology for protein function prediction. It transforms protein-level GO annotations into statistically validated domain-level associations using Fisher's exact tests, FDR correction, and hypergeometric scoring.

The supported entry point today is `run_dcgo_human.py`. It consumes **pre-computed InterPro domain annotations** (`protein2ipr.dat`) — it does not run InterProScan or scan sequences itself. There is no general `main_pipeline` orchestrator yet.

The path is **species-parameterised**: download, extract, and run all take a `--species` flag (default `human`) and share `goa_<species>` / `protein2ipr_<species>` file naming, so any organism with a GOA file works (the script name keeps `human` for historical reasons). Annotations enter through the **`AnnotationSource` seam** (`src/annotation_source.py`): the Fisher/FDR engine only ever sees `{protein → {term}}` dicts, so associating domains with a non-GO ontology (Disease Ontology, HPO, EC, Reactome …) means adding an `AnnotationSource` subclass, not touching the statistics. Broader multi-ontology work is tracked in `FUTURE_WORK.md`.

## Core Architecture

The human path runs as a sequence of scripts backed by modules in `src/`:

1. **Download** (`scripts/download_data.py`) - Fetches required datasets (GOA, GO ontology, InterPro `protein2ipr`) from source into `data/raw/<source>/`. Reads URLs from `config/settings.py`.
2. **Extract** (`extract_human_interpro.py`) - Filters the ~20 GB `protein2ipr.dat.gz` down to human proteins found in GOA, writing `data/interim/protein2ipr_human.dat.gz`.
3. **Run** (`run_dcgo_human.py`) - Orchestrates the analysis using the `src/` modules below.

Key `src/` modules:
- `annotation_source.py` - `AnnotationSource` abstraction (protein → ontology term). `GAFAnnotationSource` is the GO reference implementation; the seam for adding non-GO ontologies.
- `ontology_registry.py` - the `--ontology` dispatch table: one `OntologyEntry` per ontology, holding its `AnnotationSource` factory, its ancestors factory (or `None` when it has no hierarchy), and the input files each needs so a run fails early on a missing input. Adding an ontology is one entry here.
- `disease_ontology.py` - Disease Ontology adapter (`--ontology doid` / `orphanet_doid`). UniProt's disease layer is OMIM/Orphanet ids with no DAG; DO cross-references both (`xref: MIM:`, `xref: ORDO:`), so this re-keys the annotations onto DOID terms **at parse time** — sparse per-locus phenotypes pool into a DO class before the Fisher tests, then propagate up DO's `is_a` DAG. The module docstring states the policy for unmapped / one-to-many / obsolete ids, and every case is counted and logged (`XrefMapping`, `RemapCoverage`).
- `ec_annotation_source.py` - Enzyme Commission adapter (`run_dcgo_human.py --ontology ec`): parses Expasy `enzyme.dat` (already UniProt-keyed, so no id mapping) and provides `propagate_ec_annotations`/`ec_ancestors` for EC True Path propagation (hierarchy is implicit in the numbering — no OBO). First non-GO ontology on the seam.
- `uniprot_annotation_source.py` - UniProt-native adapters. Three layers of the Swiss-Prot flat file: `DR` cross-references (`reactome`, `disease`, `orphanet`, `tcdb`, `merops`, `cazy`, `unipathway`, `complex`, `drugbank`, `pharos`, `condensate`, plus `xref --xref-db NAME` for anything else), `KW` keywords, and — new — the layers curated into the entry *body*: `CC SUBCELLULAR LOCATION` (mapped to `SL-` terms via `subcell.txt`), `FT /ligand_id` ChEBI ligands, `CC COFACTOR` ChEBI cofactors and `CC CATALYTIC ACTIVITY` Rhea reactions. UniProt is the protein universe, so all of these are already accession-keyed — no id mapping. Also parses the Reactome/keyword/subcellular hierarchies for True Path propagation.
- `hierarchy.py` - Shared, ontology-agnostic True Path engine: `closure_ancestors` (child→parents map → transitive-ancestors fn) and `propagate_via_ancestors`, plus the hierarchy *loaders* — `dotted_ancestors` (TCDB), `alpha_prefix_ancestors` (MEROPS/CAZy) and `parse_obo_child_parents` (a light OBO reader used for ChEBI). Everything except GO propagates through this engine (GO keeps its obonet `OntologyProcessor` path, which also does parental-background filtering).
- `run_manifest.py` - Machine-readable run provenance (`run_manifest_<ontology>.json`): input/output SHA-256s and release headers, Git state, `uv.lock` hash, command line, every parameter and threshold, timestamps and summary counts. Written by `start_run_manifest`/`manifest.complete` in `run_dcgo_human.py`; checklist in `REPRODUCIBILITY.md`.
- `surprise_score.py` - Ranks *emergent* supra-domain associations: a binomial test of the combination against what its parts already predict, times a redundant-signature penalty, times a novelty discount vs InterPro2GO. Driver: `scripts/rank_surprising_associations.py`; method and results in `SURPRISE_SCORE.md`.
- `goa_parser.py` - Parses GOA GAF files (protein → GO), with evidence-code filtering. `parse_goa` is the species-agnostic API (`parse_goa_human` is a kept alias).
- `domain_annotation_parser.py` - Parses `protein2ipr` domain mappings and builds domain architectures. `--domain-key` chooses *which column is a domain*: `interpro` (the integrated entry, default) or `ssf` (the SUPERFAMILY signature, whose numeric part is the SCOP sunid — the published dcGO's domain universe, used by VALIDATION_PLAN §3). Non-matching member-database rows are dropped **at parse time**, before the sort-by-start, so supra-domain contiguity is computed over the chosen universe only.
- `sparse_fisher.py` - Sparse contingency-table construction for domain × GO.
- `vectorized_fisher.py` - Vectorized Fisher's exact tests (Cython `fisher`) + Benjamini–Hochberg FDR.
- `ontology_processor.py` - True Path Rule / GO DAG propagation (opt-in).

## Development Commands

### Environment Setup
```bash
uv sync              # runtime deps
uv sync --group dev  # + dev deps (pytest, ruff, mypy)
```

### Code Quality
CI (`.github/workflows/ci.yml`) uses ruff:
```bash
uv run ruff check src/ tests/ config/ scripts/ validation/ \
    run_dcgo_human.py extract_human_interpro.py
uv run ruff format --check
```
`mypy` is available via the dev group. `uv run mypy` reproduces CI's configured
strict fence from `pyproject.toml`; expand its `files` list as legacy annotations
are fixed. `uv run mypy src/` reports the full migration backlog.

### Testing
```bash
uv run pytest                          # all tests (632, ~8 s)
uv run pytest tests/unit -v            # unit tests
uv run pytest tests/integration -v     # integration tests
uv run pytest --cov=src --cov-report=html
uv run pytest tests/unit/test_vectorized_fisher.py -v   # single file
```

### Pipeline Execution (human path)
```bash
uv run python scripts/download_data.py        # download required datasets
uv run python extract_human_interpro.py       # one-time human subset extraction
uv run python run_dcgo_human.py --num-cores 8 # statistical inference

# With True Path Rule propagation (paper Step 3 — adds ancestor annotations)
uv run python run_dcgo_human.py --num-cores 8 \
    --enable-true-path --go-ontology data/raw/go_ontology/go-basic.obo

# With the parental-background filter (paper Step 2 "relative inference" — GO
# only, and it *removes* associations). Independent of --enable-true-path;
# passing both reproduces what --enable-true-path alone used to do for GO.
uv run python run_dcgo_human.py --enable-relative-inference

# Other ontologies (see src/ontology_registry.py or --help for all 21)
uv run python run_dcgo_human.py --ontology subcellular --enable-true-path
uv run python run_dcgo_human.py --ontology ligand      # FT /ligand_id (ChEBI)
uv run python run_dcgo_human.py --ontology doid --enable-true-path  # OMIM re-keyed to DO

# Calibration control: shuffle protein↔term-set assignment; a well-behaved layer
# returns ~0 significant associations. Writes domain_<ontology>_permuted<seed>_*.
uv run python run_dcgo_human.py --ontology doid --permute-annotations 7

# Key domains by SCOP superfamily instead of InterPro entry (VALIDATION_PLAN §3)
uv run python run_dcgo_human.py --domain-key ssf --output-dir results_ssf

# Rank the emergent domain-combination predictions
uv run python scripts/rank_surprising_associations.py --ontology go

# Compare against the published dcGO (needs --group dcgo-reference downloaded)
uv run python scripts/download_data.py --group dcgo-reference
uv run python validation/compare_original_dcgo.py \
    --associations results_ssf/domain_ssf_go_associations_significant.tsv

# Component ablation + permutation null + paired bootstrap CIs (VALIDATION_PLAN §4).
# Needs one pipeline run per rung under --run-dir; see validation/ablation.py.
uv run python validation/ablation.py \
    --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \
    --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
    --run-dir results/ablation --n-bootstrap 1000 --n-permutations 200

# HPC batch script
sbatch scripts/run_dcgo_hpc.sh
```

## Data Flow

**External inputs** (downloaded by `scripts/download_data.py`):
- GOA annotations — protein → GO mappings (GAF 2.2)
- GO ontology — `go-basic.obo` (needed for `--enable-true-path` or `--enable-relative-inference`)
- InterPro mappings — pre-computed `protein2ipr.dat.gz` domain annotations

**Internal flow:**
```
download_data.py → extract_human_interpro.py → run_dcgo_human.py
    goa_parser + domain_annotation_parser
        → sparse_fisher → vectorized_fisher (+ hierarchical_inference)
        → ontology_processor (optional)
        → results TSV
```

## Configuration System

`config/settings.py` centralizes dataset URLs and parameters. Note that the
supported `run_dcgo_human.py` path takes most parameters via CLI flags rather
than from `Config`. Key defaults:
- `FDR_THRESHOLD = 0.01` - False discovery rate cutoff
- `MIN_PROTEINS_PER_ASSOCIATION = 3` - **not applied by the supported path.**
  `run_dcgo_human.py` keeps an association on FDR significance alone unless
  `--min-support N` is passed. The default is deliberately no filter: the
  emergent domain combinations this method exists to find sit at n = 2-8
  proteins, so a non-zero default would delete them. `--min-support` is applied
  *after* the BH correction, so it never alters the hypothesis family.
- `MAX_SUPRA_DOMAIN_LENGTH = 3` - Maximum contiguous domain combinations
- `NUM_CORES = 8` - Parallel processing default

## Performance Considerations

For the human path (single machine):
- **Extraction**: streams the ~20 GB `protein2ipr.dat.gz` once (~10 min). This is
  now the dominant cost of a first run.
- **Inference**: seconds, not minutes. The Fisher stage enumerates only the
  domain–term pairs that actually co-occur; every other pair has p=1 exactly
  under the one-sided `greater` test, so it is skipped and BH corrects against
  the full hypothesis count. Human single-domain is 318,749,661 hypotheses of
  which 655,659 are evaluated (**1.6 s**, was 323 s dense); human supra-domain is
  1,690,803,963 / 3,008,670 (**5.4 s**, was 2,123 s).
- **Bottleneck**: the sparse product and the BH sort, both proportional to the
  *co-occurring* pair count rather than the dense product.

Multi-organism runs are implemented — see `MULTISPECIES_BACKGROUND.md`. The
all-species supra design is 13.1e9 hypotheses, unrunnable densely (~389 GB of
tables); enumerating co-occurring pairs makes it 9.5M tables and 268 s.

## Known Limitations

- No local domain scanning — only pre-computed InterPro annotations are consumed.
- **The component ablation (§4, 2026-08-04) is negative for two of three
  components — read `VALIDATION_PLAN.md` §4 before claiming any of them helps.**
  Over 12 aspect × IC cells with a paired protein-level bootstrap: supra-domains
  improve 0/12 and the True Path Rule is *significantly
  worse* in 12/12. On the §2 benchmark the best configuration is single domains
  only. The supra-domain machinery's demonstrated value is the emergent
  combinations in `SURPRISE_SCORE.md`, not protein-centric F_max.
- **`--enable-shrinkage` was removed** (2026-08-05). It geometrically
  interpolated each supra-domain p-value toward the geometric mean of its
  constituents', which pulled thin evidence *toward* its well-supported parts —
  a 3-protein combination at p=0.01 became p=1e-24 — taking FDR<0.01 rejections
  from 163,277 to 463,924 (+184%). The output was not a valid p-value under any
  null, so BH did not control FDR on it, and the ablation found no effect on
  prediction quality in any of 12 cells. A genuine version would shrink the
  observed *rate* and recompute Fisher, which is a different method.
- True Path Rule is opt-in (`--enable-true-path`), not part of the default run;
  it now errors out for ontologies with no hierarchy instead of silently
  skipping. **It is propagation and nothing else.** The parental-background
  filter it used to run alongside for GO is the paper's separate *relative
  inference* (Step 2, vs. the true-path rule at Step 3) and is now
  `--enable-relative-inference`, available for all 12 ontologies with a
  hierarchy, and it *removes* associations rather than adding them. Our version is still a post-hoc `alpha < 0.05` filter
  applied after BH, where the paper combines the overall and relative p-values
  *before* correcting; that gap is VALIDATION_PLAN next-steps item 2.
- **The §4 ablation cannot attribute its True Path result to either stage.** It
  was run when one flag drove both, so "the True Path Rule is significantly
  worse in 12/12 cells" is a statement about filter-plus-propagation, measured
  additionally with the unpropagated-background defect in place (54,951
  rejections; fixed in #46, now 237 on the same run). Propagation only adds
  annotations and cannot by itself lower recall. Re-run the ablation against the
  split flags before citing that number.
- The surprise score re-ranks associations that already passed the dcGO FDR
  filter, using the same proteins — it measures internal consistency of the
  evidence, not out-of-sample performance. It is also **not a total order**: on
  the GO run 9,923 of 10,136 evaluated associations score exactly 0.000, so any
  comparison reaching deeper than its ~213 scored associations is comparing an
  arbitrary tie-break. `validation/temporal_surprise.py` detects this and refuses
  to quote an interval rather than reporting one — see `SURPRISE_SCORE.md`.
- Validation covers InterPro2GO coverage (§1, ~65%), a temporal CAFA-style
  benchmark (§2, `validation/temporal_benchmark.py`), the multi-ontology breadth
  test (§2, `validation/temporal_breadth.py`), the published-dcGO comparison
  (§3, `validation/compare_original_dcgo.py`: precision 0.54–0.63 on the shared
  SCOP superfamily space; **recall against them is not interpretable** — they
  are all-species and 2016, we are human-only and 2026), and the §4 ablation
  with bootstrap CIs and a permutation null (`validation/ablation.py`,
  `validation/resampling.py`). dcGO clears a 100–200-permutation random-domain
  null in every cell at the attainable p-floor, and beats the naive **F_max**
  baseline on informative terms — but **loses to naive on AUPRC** at IC≥0 in all
  aspects. Still open: score calibration (§5), an untouched evaluation interval.
  See `VALIDATION_PLAN.md`.
- **The background is human-only by default, and widening it helps** — see
  `MULTISPECIES_BACKGROUND.md`. `--species allspecies` runs a 1,464,355-protein,
  9,074-taxon universe. On the held-out 2021→2026 split, with the evaluation
  fixed and only the training universe changed, it wins 8/9 F_max and 9/9 AUPRC
  cells; under `--evidence-filter experimental` it wins 9/9 and 9/9. Two caveats
  belong with any number taken from it: the `manual` universe is **75.8%
  projected** annotation, of which 55.2% of the non-human part cites a human
  protein, and support is **inflated ~2.44×** by orthology with half the
  associations resting on ≤3 UniRef50 clusters. Neither touches the held-out
  result; both bear on the significant counts and their FDR.

## Conventions

The codebase targets Python 3.12 and uses dataclasses, type hints, pathlib, and context managers throughout.
