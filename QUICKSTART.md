# dcGOspeed Quick Start Guide

This guide shows how to run the dcGO pipeline for human proteins from start to finish.

## Prerequisites

```bash
# Install dependencies
uv sync

# Or with pip
pip install -e .
```

## Quick Start: Human Proteins (Recommended)

The pipeline has been optimized for human protein analysis using pre-computed InterPro annotations.

### Step 1: Download Required Data

```bash
# Download GOA, GO ontology, and InterPro mappings
uv run python -m src.data_acquisition \
    --datasets goa_annotations \
    --datasets go_ontology \
    --datasets interpro_mappings
```

This downloads:
- `goa_human.gaf.gz` (~11 MB) - Human GO annotations
- `go-basic.obo` (~30 MB) - GO ontology structure
- `protein2ipr.dat.gz` (~20 GB) - InterPro domain annotations for all organisms

### Step 2: Extract Human-Specific Annotations

Since `protein2ipr.dat.gz` contains ALL organisms, we extract only human proteins:

```bash
uv run python extract_human_interpro.py
```

This creates:
- `data/interim/protein2ipr_human.dat.gz` (~3 MB) - Human proteins only
- `data/interim/human_proteins.txt` - List of 18,966 human protein IDs

**Time**: ~10 minutes (one-time operation)

### Step 3: Run Statistical Inference

Run the complete statistical inference pipeline:

```bash
uv run python test_sparse_fisher.py
```

This performs:
1. Parse GOA annotations (2s) - 18,966 proteins, 595K annotations, evidence filtering
2. Parse InterPro domains (<1s) - 18,783 proteins, 232K annotations
3. Build sparse matrices (<1s) - Efficient protein×domain and protein×GO matrices
4. Compute contingency tables (~5 min) - All 303M domain-GO combinations
5. Run Fisher's exact tests (~20 min) - Parallel processing at 30K+ tests/second
6. Apply FDR correction (<1s) - Benjamini-Hochberg with α=0.01

**Total time**: ~25 minutes
**Output**: Significant domain-GO associations with FDR < 0.01

## Configuration Options

### Evidence Code Filtering

Control GO annotation quality in `config/settings.py` or via command line:

```python
# config/settings.py
evidence_filter: str = 'manual'  # Default: exclude electronic (IEA) annotations

# Options:
# 'all'          - Include ALL annotations (including IEA)
# 'manual'       - Exclude IEA (only manually curated) - RECOMMENDED
# 'experimental' - Only experimental evidence (highest confidence)
```

See `ONTOLOGY_OPTIONS.md` for detailed evidence code information.

### Statistical Parameters

In `config/settings.py`:

```python
fdr_threshold: float = 0.01              # FDR cutoff (1%)
min_proteins_per_association: int = 3    # Minimum evidence
max_supra_domain_length: int = 3         # Max contiguous domains
alpha_threshold: float = 0.05            # Significance threshold
```

## Output Files

After running the pipeline:

```
results/
├── domain_go_associations.tsv       # Significant associations
├── domain_go_associations_all.tsv   # All tested associations
├── statistics_summary.json          # Performance metrics
└── dcgo_database.db                 # SQLite database (if using full pipeline)

data/interim/
├── protein2ipr_human.dat.gz         # Human domain annotations
└── human_proteins.txt               # Human protein IDs
```

## Advanced Usage

### Running with Different Parameters

```bash
# Use experimental evidence only
uv run python test_sparse_fisher.py --evidence-filter experimental

# Adjust FDR threshold
uv run python test_sparse_fisher.py --fdr-threshold 0.05

# Use more/fewer CPU cores
uv run python test_sparse_fisher.py --num-cores 16
```

### Testing with Subsets

To test with a smaller dataset (faster):

```bash
# Test with top 100 domains and GO terms only
uv run python test_vectorized_fisher.py
```

## Performance Notes

**Human dataset** (18,783 proteins):
- 18,705 unique domains
- 16,241 unique GO terms
- 303,787,905 total tests
- ~25 minutes total time on 8 cores

**Memory usage**:
- Peak: ~15 GB (during Fisher tests)
- Typical: ~5 GB

**Scaling**:
- Fisher tests: Linear with CPU cores
- Memory: Linear with number of tests
- For 1B tests: ~50 GB RAM recommended

## Troubleshooting

### Out of Memory
- Reduce batch size in Fisher tests
- Process domains in chunks
- Use fewer GO terms (filter by annotation count)

### Slow Performance
- Ensure `data/interim/protein2ipr_human.dat.gz` exists (don't parse full 20GB file)
- Use `evidence_filter='manual'` (excludes ~30% low-quality annotations)
- Increase `num_cores` for parallel processing

### Missing Files
```bash
# Check file structure
ls -lh data/raw/goa_annotations/
ls -lh data/raw/go_ontology/
ls -lh data/raw/interpro_mappings/
ls -lh data/interim/

# Re-extract human proteins if needed
uv run python extract_human_interpro.py
```

## Next Steps

1. **Explore Results**: Open `results/domain_go_associations.tsv` in Excel/R
2. **Integrate with Ontology**: Use `src/ontology_processor.py` for GO term propagation
3. **Database Queries**: Query `dcgo_database.db` with SQL for custom analyses
4. **Visualizations**: Plot enrichment heatmaps, domain-GO networks

## File Organization

```
dcGOspeed/
├── data/
│   ├── raw/                          # Downloaded data (original filenames preserved)
│   │   ├── goa_annotations/
│   │   │   └── goa_human.gaf.gz
│   │   ├── go_ontology/
│   │   │   └── go-basic.obo
│   │   └── interpro_mappings/
│   │       └── protein2ipr.dat.gz    # 20GB - all organisms
│   ├── interim/                      # Processed data
│   │   ├── protein2ipr_human.dat.gz  # 3MB - human only
│   │   └── human_proteins.txt
│   └── processed/                    # Final outputs
├── results/                          # Analysis results
├── logs/                             # Log files
├── config/                           # Configuration
│   └── settings.py
└── src/                              # Source code
    ├── goa_parser.py                 # GO annotation parser
    ├── domain_annotation_parser.py   # InterPro parser
    ├── sparse_fisher.py              # Sparse matrix operations
    ├── vectorized_fisher.py          # Parallel Fisher tests
    └── ...
```

## References

- **dcGO Methodology**: Fang & Gough (2013) Nucleic Acids Research
- **InterPro**: https://www.ebi.ac.uk/interpro/
- **Gene Ontology**: http://geneontology.org/
- **GOA**: https://www.ebi.ac.uk/GOA

For more details, see:
- `README.md` - Project overview
- `ONTOLOGY_OPTIONS.md` - Evidence code filtering
- `DATA_FLOW_CLARIFICATION.md` - Data sources and workflow
