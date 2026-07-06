# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

dcGO-2.0 is a bioinformatics pipeline implementing the domain-centric Gene Ontology (dcGO) methodology for protein function prediction. It transforms protein-level GO annotations into statistically validated domain-level associations using Fisher's exact tests, FDR correction, and hypergeometric scoring.

The supported entry point today is the **human analysis path** (`run_dcgo_human.py`). It consumes **pre-computed InterPro domain annotations** (`protein2ipr.dat`) — it does not run InterProScan or scan sequences itself. There is no general `main_pipeline` orchestrator yet; a multi-organism version is future work (see `FUTURE_WORK.md`).

## Core Architecture

The human path runs as a sequence of scripts backed by modules in `src/`:

1. **Download** (`scripts/download_data.py`) - Fetches required datasets (GOA, GO ontology, InterPro `protein2ipr`) from source into `data/raw/<source>/`. Reads URLs from `config/settings.py`.
2. **Extract** (`extract_human_interpro.py`) - Filters the ~20 GB `protein2ipr.dat.gz` down to human proteins found in GOA, writing `data/interim/protein2ipr_human.dat.gz`.
3. **Run** (`run_dcgo_human.py`) - Orchestrates the analysis using the `src/` modules below.

Key `src/` modules:
- `goa_parser.py` - Parses GOA GAF files (protein → GO), with evidence-code filtering.
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
uv run pytest                          # all tests (113, ~3 s)
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
- Validation is against InterPro2GO only; no CAFA / temporal benchmark or
  comparison to the original dcGO results yet. See `VALIDATION_PLAN.md`.

## Package Structure

```
run_dcgo_human.py            # Main entry point
extract_human_interpro.py    # Human subset extraction
scripts/download_data.py     # Dataset downloader
src/
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
