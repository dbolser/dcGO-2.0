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
- `ec_annotation_source.py` - Enzyme Commission adapter (`run_dcgo_human.py --ontology ec`): parses Expasy `enzyme.dat` (already UniProt-keyed, so no id mapping) and provides `propagate_ec_annotations`/`ec_ancestors` for EC True Path propagation (hierarchy is implicit in the numbering — no OBO). First non-GO ontology on the seam.
- `uniprot_annotation_source.py` - UniProt-native adapters (`--ontology reactome|keyword|disease|xref`): harvest DR cross-references and KW keywords straight from the Swiss-Prot flat file. `disease` = `DR MIM` phenotype-only; `xref --xref-db NAME` opens any DR database (KEGG/Orphanet/DisGeNET/…), with an optional `--xref-type` filter. UniProt is the protein universe, so these terms are already accession-keyed — no id mapping. Also parses the Reactome/keyword hierarchies (`parse_reactome_relations`, `parse_keyword_hierarchy`) for True Path propagation.
- `hierarchy.py` - Shared, ontology-agnostic True Path engine: `closure_ancestors` (child→parents map → transitive-ancestors fn) and `propagate_via_ancestors`. EC, Reactome, and keywords all propagate through it (GO keeps its obonet `OntologyProcessor` path).
- `goa_parser.py` - Parses GOA GAF files (protein → GO), with evidence-code filtering. `parse_goa` is the species-agnostic API (`parse_goa_human` is a kept alias).
- `domain_annotation_parser.py` - Parses `protein2ipr` domain mappings and builds domain architectures.
- `sparse_fisher.py` - Sparse contingency-table construction for domain × GO.
- `vectorized_fisher.py` - Vectorized Fisher's exact tests (Cython `fisher`) + Benjamini–Hochberg FDR.
- `hierarchical_inference.py` - Supra-domain generation and optional empirical-Bayes shrinkage.
- `ontology_processor.py` - True Path Rule / GO DAG propagation (opt-in).
- `database_manager.py` - SQLite storage/export helpers.
- `data_acquisition.py` - Older async (aiohttp) downloader library; **not** on the supported path (use `scripts/download_data.py`).

## Development Commands

### Environment Setup
```bash
uv sync              # runtime deps
uv sync --group dev  # + dev deps (pytest, ruff, mypy)
```

### Code Quality
CI (`.github/workflows/ci.yml`) uses ruff:
```bash
uv run ruff check src/ tests/
uv run ruff format --check
```
`mypy` is available via the dev group (`uv run mypy src/`) but is not enforced in CI.

### Testing
```bash
uv run pytest                          # all tests (155, ~5 s)
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

# With True Path Rule propagation
uv run python run_dcgo_human.py --num-cores 8 \
    --enable-true-path --go-ontology data/raw/go_ontology/go-basic.obo

# HPC batch script
sbatch scripts/run_dcgo_hpc.sh
```

## Data Flow

**External inputs** (downloaded by `scripts/download_data.py`):
- GOA annotations — protein → GO mappings (GAF 2.2)
- GO ontology — `go-basic.obo` (only needed for `--enable-true-path`)
- InterPro mappings — pre-computed `protein2ipr.dat.gz` domain annotations

**Internal flow:**
```
download_data.py → extract_human_interpro.py → run_dcgo_human.py
    goa_parser + domain_annotation_parser
        → sparse_fisher → vectorized_fisher (+ hierarchical_inference)
        → ontology_processor (optional)
        → results TSV / database_manager
```

## Configuration System

`config/settings.py` centralizes dataset URLs and parameters. Note that the
supported `run_dcgo_human.py` path takes most parameters via CLI flags rather
than from `Config`. Key defaults:
- `FDR_THRESHOLD = 0.01` - False discovery rate cutoff
- `MIN_PROTEINS_PER_ASSOCIATION = 3` - Minimum evidence requirement
- `MAX_SUPRA_DOMAIN_LENGTH = 3` - Maximum contiguous domain combinations
- `NUM_CORES = 8` - Parallel processing default

## Testing Architecture

- `tests/unit/` - Individual component tests (Fisher, parsers, ontology, supra-domains)
- `tests/integration/` - Multi-component workflow tests (e.g. True Path pipeline)
- `tests/e2e/` - Reserved for full-pipeline tests

## Performance Considerations

For the human path (single machine):
- **Extraction**: streams the ~20 GB `protein2ipr.dat.gz` once (~10 min).
- **Inference**: ~300M Fisher tests in ~50 min on 8 cores (~100k tests/s).
- **Bottlenecks**: (1) contingency-table construction (memory), (2) Fisher tests
  on millions of domain–GO pairs (CPU).

A full multi-organism run would be substantially heavier and is not yet implemented.

## Known Limitations

- No local domain scanning — only pre-computed InterPro annotations are consumed.
- True Path Rule is opt-in (`--enable-true-path`), not part of the default run.
- Validation covers InterPro2GO coverage (§1, ~65%) and a temporal CAFA-style
  benchmark (§2, `validation/temporal_benchmark.py`: dcGO beats the random-domain
  null 1.7–3.2× but trails the naive F_max baseline). Still open: ablation (§4),
  score calibration (§5), original-dcGO comparison (§3). See `VALIDATION_PLAN.md`.

## Package Structure

```
run_dcgo_human.py            # Main entry point
extract_human_interpro.py    # Human subset extraction
scripts/download_data.py     # Dataset downloader
src/
├── annotation_source.py     # AnnotationSource seam (protein → ontology term)
├── ec_annotation_source.py  # Enzyme Commission adapter (Expasy enzyme.dat)
├── uniprot_annotation_source.py # UniProt-native adapters (Reactome/keywords/disease/xref)
├── hierarchy.py             # Shared True Path engine (ancestors + propagation)
├── goa_parser.py            # GOA GAF parsing
├── domain_annotation_parser.py  # protein2ipr parsing
├── sparse_fisher.py         # Sparse contingency tables
├── vectorized_fisher.py     # Fisher's exact tests + BH-FDR
├── hierarchical_inference.py    # Supra-domains + shrinkage
├── ontology_processor.py    # True Path Rule
├── database_manager.py      # SQLite operations
└── data_acquisition.py      # Legacy async downloader (unused by supported path)
config/settings.py           # Dataset URLs + configuration
```

The codebase targets Python 3.12 and uses dataclasses, type hints, pathlib, and context managers throughout.
