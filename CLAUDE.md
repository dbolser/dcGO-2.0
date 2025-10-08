# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

dcGOspeed is a production-ready bioinformatics pipeline implementing the domain-centric Gene Ontology (dcGO) methodology for protein function prediction. The pipeline transforms protein-level functional annotations into statistically validated domain-level associations using rigorous statistical inference.

## Core Architecture

The codebase follows a modular pipeline architecture with five main processing stages:

1. **Data Acquisition** (`src/data_acquisition.py`) - Downloads and manages UniProt, GOA, Pfam, and GO ontology datasets
2. **Domain Scanning** (`src/domain_scanning.py`) - Integrates with InterProScan to identify protein domain architectures and generates supra-domains
3. **Statistical Inference** (`src/statistical_inference.py`) - Performs Fisher's exact tests with FDR correction and hypergeometric scoring
4. **Ontology Processing** (`src/ontology_processor.py`) - Implements the True Path Rule with optimal level filtering and annotation propagation
5. **Database Management** (`src/database_manager.py`) - SQLite backend with performance optimization and export capabilities

The main orchestrator (`src/main_pipeline.py`) coordinates all modules with checkpoint/resume functionality for long-running HPC jobs.

## Development Commands

### Environment Setup
```bash
# Setup with uv (recommended)
uv sync
uv sync --group dev  # Include development dependencies

# Alternative with pip
pip install -e .
pip install -e ".[dev]"
```

### Code Quality
```bash
# Format code
black src/ tests/
isort src/ tests/

# Lint
flake8 src/ tests/

# Type check
mypy src/

# Run all quality checks
black src/ tests/ && isort src/ tests/ && flake8 src/ tests/ && mypy src/
```

### Testing
```bash
# Run all tests
pytest

# Run specific test categories
pytest -m unit          # Unit tests only
pytest -m integration   # Integration tests only
pytest -m e2e           # End-to-end tests only

# Run with coverage
pytest --cov=src/dcgo_pipeline --cov-report=html

# Run single test file
pytest tests/test_statistical_inference.py -v

# Run specific test
pytest tests/test_domain_scanning.py::TestDomainScanner::test_supra_domain_generation -v
```

### Pipeline Execution
```bash
# Run complete pipeline (requires large datasets and InterProScan)
python -m src.main_pipeline --num-cores 8 --output-dir results/

# Run with custom configuration
python -m src.main_pipeline --config config/custom.yaml

# HPC execution
sbatch scripts/run_dcgo_hpc.sh

# Install InterProScan and dependencies
bash scripts/install.sh
```

## Module Dependencies and Data Flow

**External Dependencies:**
- InterProScan 5.67+ (domain scanning)
- UniProt datasets (protein sequences)
- GOA annotations (protein-GO mappings)
- Pfam HMM profiles (domain definitions)
- GO ontology (hierarchical structure)

**Internal Data Flow:**
```
DataAcquisition → DomainScanner → StatisticalInference
                      ↓               ↓
        OntologyProcessor ← DatabaseManager
```

**Key Data Structures:**
- `AssociationResult` - Statistical test results with contingency table data
- `Annotation` - Final domain-GO associations with propagation info
- `DomainArchitecture` - Protein domain arrangements including supra-domains
- `PipelineCheckpoint` - State management for long-running processes

## Configuration System

The pipeline uses a centralized configuration system (`config/settings.py`) with:
- Data source URLs and paths
- Statistical thresholds (FDR < 0.01)
- Processing parameters (cores, batch sizes)
- Database settings and paths

Key configuration points:
- `FDR_THRESHOLD = 0.01` - False discovery rate cutoff
- `MIN_PROTEINS_PER_ASSOCIATION = 3` - Minimum evidence requirement
- `MAX_SUPRA_DOMAIN_LENGTH = 3` - Maximum contiguous domain combinations
- `NUM_CORES = 8` - Parallel processing setting

## Testing Architecture

The test suite is organized by scope:
- `tests/unit/` - Individual component tests
- `tests/integration/` - Multi-component workflow tests  
- `tests/e2e/` - Full pipeline tests with mock data

Critical test markers:
- `@pytest.mark.slow` - Computationally intensive tests
- `@pytest.mark.integration` - Cross-module tests
- `@pytest.mark.unit` - Isolated component tests

## Performance Considerations

The pipeline is designed for HPC environments processing millions of proteins:
- **Memory**: Peak ~200GB during InterProScan execution
- **Runtime**: 3-7 days on HPC cluster depending on dataset size
- **Storage**: ~500GB for intermediate files, ~10GB for final database
- **Parallelization**: Scales linearly with available cores for domain scanning

Critical performance bottlenecks:
1. InterProScan domain scanning (most compute-intensive)
2. Correspondence matrix construction (memory-intensive)
3. Fisher's exact tests on millions of domain-GO pairs

## Production Deployment

The pipeline includes production-ready features:
- Checkpoint/resume functionality for HPC job interruptions
- Comprehensive logging with structured output
- Error recovery and validation at each stage
- SQLite database with proper indexing for performance
- Docker containerization support

For HPC deployment, use `scripts/run_dcgo_hpc.sh` with appropriate SLURM parameters adjusted for your cluster specifications.

## Package Structure

```
src/
├── dcgo_pipeline/           # Main package (future modular components)
│   ├── core/               # Core pipeline functionality
│   ├── utils/              # Utility functions
│   ├── validation/         # Validation utilities
│   ├── integration/        # External tool integrations
│   └── visualization/      # Plotting and visualization
├── data_acquisition.py     # Dataset download and management
├── domain_scanning.py      # InterProScan integration
├── statistical_inference.py # Fisher's tests and FDR correction
├── ontology_processor.py   # True Path Rule implementation
├── database_manager.py     # SQLite operations
└── main_pipeline.py        # Pipeline orchestration
```

The codebase uses modern Python 3.12 features including dataclasses, type hints, pathlib, and context managers throughout. All modules include comprehensive error handling, logging, and documentation.