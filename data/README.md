# Data Directory Structure

This directory contains all raw and processed data files for the dcGO pipeline.


## Organization

The data acquisition module preserves original filenames and organizes files by
source for better provenance tracking:

```
data/
├── raw/                        # Raw downloaded files organized by source
│   ├── goa_annotations/
│   │   └── goa_human.gaf.gz    # GOA annotations (from EBI)
│   ├── go_ontology/
│   │   └── go-basic.obo        # GO ontology (from GO Consortium)
│   ├── interpro_mappings/
│   │   └── protein2ipr.dat.gz  # Pre-computed InterPro mappings (from EBI)
│   └── [other sources]/
├── processed/                  # Intermediate processed files
│   ├── domains/                # Parsed domain architectures
│   ├── ontology/               # Processed GO data
│   └── annotations/            # Parsed annotations
├── interim/                    # Temporary processing files
└── cache/                      # Cached results for faster reruns
```


## Data Sources

### Required Files

| File | Source | Size | URL Pattern |
|------|--------|------|-------------|
| `goa_human.gaf.gz` | GOA @ EBI | ~11MB | https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz |
| `go-basic.obo` | GO Consortium | ~30MB | http://current.geneontology.org/ontology/go-basic.obo |
| `protein2ipr.dat.gz` | InterPro @ EBI | ~21GB | https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz |


### Optional Files

| File | Source | Purpose |
|------|--------|---------|
| `uniprot_sprot.fasta.gz` | UniProt | Only for local HMM scanning |
| `Pfam-A.hmm.gz` | Pfam | Only for local HMM scanning |


## File Formats

- **GAF (Gene Association Format)**: Tab-separated GO annotations
- **OBO (Open Biomedical Ontologies)**: GO ontology hierarchy
- **protein2ipr.dat**: Tab-separated InterPro domain mappings


## Provenance

All files are downloaded directly from authoritative sources and organized in
subdirectories matching the source name. Original filenames are preserved to
facilitate:

- Verification against upstream sources
- Understanding data versioning
- Reproducibility of analyses
- Troubleshooting and debugging


## Updating Data

To download the latest versions:

```bash
# Download specific dataset
python -m src.data_acquisition download --dataset goa_annotations

# Force redownload
python -m src.data_acquisition download --force-redownload

# Download all required datasets
uv run python -m src.main_pipeline --datasets goa_annotations --datasets go_ontology --datasets interpro_mappings
```
