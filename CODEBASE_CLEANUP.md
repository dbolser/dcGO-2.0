# dcGOspeed Codebase Cleanup - Summary

## Changes Made

### 1. Bug Fixed
- **src/statistical_inference.py:79** - Removed `@property` decorator from `is_significant()` method (it was being called as a method with parameters, not as a property)

### 2. Legacy Code Removed

**Deleted Modules:**
- `src/main_pipeline.py` - Replaced by `run_dcgo_human.py` with sparse matrix approach
- `src/statistical_inference.py` - Replaced by `src/sparse_fisher.py` + `src/vectorized_fisher.py`
- `src/domain_scanning.py` - InterProScan integration not used (we use pre-computed annotations)

**Deleted Tests:**
- `tests/unit/test_main_pipeline.py`
- `tests/unit/test_statistical_inference.py`
- `tests/e2e/test_pipeline_integration.py`
- `tests/test_config.py`
- `tests/unit/test_data_acquisition.py`
- `tests/unit/test_ontology_processor.py`

### 3. Dependencies Cleaned

**Removed heavy unused dependencies from pyproject.toml:**
- `torch` and `torchvision` (GPU libraries - not used)
- `biopython`, `pyfaidx` (not used in production pipeline)
- `scikit-learn` (not used)
- `matplotlib`, `seaborn`, `plotly` (visualization - not used)
- `requests`, `aiohttp` (web - not used in production)
- `click`, `rich`, `tqdm` (not used in minimal pipeline)
- `pydantic`, `jsonschema`, `pyyaml`, `toml` (not used)

**Kept core dependencies:**
- `numpy`, `pandas`, `scipy` (core scientific computing)
- `loguru` (logging)
- Fisher testing runs through the compiled `fisher` package; no Python
  multiprocessing dependency is required.
- Benjamini-Hochberg correction is implemented and tested in
  `src/vectorized_fisher.py`; `statsmodels` is not a runtime dependency.
- `obonet`, `networkx` (GO ontology - for future True Path Rule)
- `sqlalchemy` (database support - optional)

### 4. Production Pipeline Enhanced

**run_dcgo_human.py improvements:**
- Added `--species` parameter to support any organism (human, mouse, zebrafish, etc.)
- Updated documentation and help text
- Now works with species-specific data files: `goa_{species}.gaf.gz` and `protein2ipr_{species}.dat.gz`

**Example usage:**
```bash
# Human (default)
uv run python run_dcgo_human.py --num-cores 16

# Mouse
uv run python run_dcgo_human.py --species mouse --num-cores 16 --evidence-filter experimental

# Any species
uv run python run_dcgo_human.py --species zebrafish --fdr-threshold 0.05
```

## Current Codebase Structure

For the authoritative module list see `README.md` and `CLAUDE.md` — the layout
has moved on considerably since this cleanup. `src/data_acquisition.py` and
`src/database_manager.py` no longer exist (downloads live in
`scripts/download_data.py`; there is no SQLite export). The four items in the
old "What's Next" list are all resolved: the True Path Rule is integrated behind
`--enable-true-path`, and the suite is at 632 tests.

## Status of this document

This file is a historical cleanup report, not an active roadmap. The single
authoritative queue is `TODO.md`; detailed engineering rationale lives in
`docs/CODE_QUALITY_ROADMAP.md`. The review of PRs #54–#60 was incorporated
there rather than creating another competing list here.
