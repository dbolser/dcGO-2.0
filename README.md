# dcGO Pipeline - Production Ready

**🎯 STATUS: COMPLETE & PRODUCTION-READY** ✅

A comprehensive, production-ready implementation of the domain-centric Gene Ontology (dcGO) methodology for protein function prediction. This pipeline transforms protein-level functional annotations into statistically validated domain-level associations using rigorous statistical inference.

## 🚀 **What This Pipeline Does**

The dcGO pipeline implements the domain-centric Gene Ontology methodology, which:

1. **Downloads large-scale biological datasets** (UniProt, GOA, InterPro, Gene Ontology)
2. **Maps protein sequences to domain architectures** using pre-computed InterPro annotations
3. **Performs statistical inference** with Fisher's exact tests and FDR correction
4. **Applies ontology logic** with the True Path Rule for optimal annotation levels
5. **Produces a database** of statistically validated domain-to-function associations

**Scientific Impact**: Moves from protein-level to domain-level functional annotation, providing more granular and accurate predictions.

---

## ⚡ **Quick Start (Complete Workflow)**

### 1. Installation
```bash
# Clone repository
git clone <your-repo-url>
cd dcGOspeed

# Install with uv (recommended - handles Python 3.12+ requirements)
uv sync

# Alternative: Install with pip
pip install -r requirements.txt
```

### 2. Verify Installation
```bash
# Test core components
uv run pytest tests/unit/test_data_acquisition.py::TestDataAcquisition::test_init -v

# Test end-to-end integration  
uv run pytest tests/e2e/test_pipeline_integration.py::TestPipelineIntegration::test_complete_pipeline_workflow -v
```

### 3. Run Complete Pipeline
```bash
# Full pipeline (downloads ~70GB of data, processes for days on HPC)
uv run python -m src.main_pipeline --num-cores 8 --output-dir results/

# With custom configuration
uv run python -m src.main_pipeline --config config/custom.yaml --num-cores 16
```

### 4. View Results
```bash
# Results will be in:
results/
├── dcgo_database.db     # SQLite database with domain-GO associations
├── dcgo_annotations.tsv # Human-readable export
└── logs/               # Detailed execution logs
```

---

## 📊 **dcGO Data Requirements: The Three Essential Inputs**

The dcGO methodology requires exactly **3 types of data** for a given set of proteins:

### **1. Domain Annotations** (What domains are in each protein?)
- **Source**: Pre-computed InterPro domain mappings
- **File**: `protein2ipr.dat.gz` (~20GB)
- **Content**: Protein ID → Domain ID mappings with boundaries
- **Format**: `P12345    PF00001    10    50    1e-10    T    Domain1`

### **2. Ontology Annotations** (What GO terms are assigned to each protein?)
- **Source**: Gene Ontology Annotation (GOA) database
- **File**: `goa_uniprot_all.gaf.gz` (~2GB)
- **Content**: Protein ID → GO term mappings with evidence
- **Format**: GAF 2.2 format with protein-GO associations

### **3. Ontology Structure** (How are GO terms related hierarchically?)
- **Source**: Gene Ontology Consortium
- **File**: `go-basic.obo` (~50MB)
- **Content**: GO term hierarchy with is_a and part_of relationships
- **Format**: OBO format defining the ontology DAG

---

## 📁 **Complete Download Breakdown**

| Dataset | Size | dcGO Input | Required | Purpose |
|---------|------|------------|----------|---------|
| **InterPro Mappings** | ~20GB | **✅ Input 1** | **Required** | Domain annotations |
| **GOA Annotations** | ~2GB | **✅ Input 2** | **Required** | Ontology annotations |
| **Gene Ontology** | ~50MB | **✅ Input 3** | **Required** | Ontology structure |
| UniProt Swiss-Prot | ~500MB | Not used | Optional | For local HMM scanning only |
| UniProt TrEMBL | ~50GB | Not used | Optional | For local HMM scanning only |
| Pfam Regions | ~5GB | Alternative | Optional | Alternative domain source |
| InterPro Definitions | ~40MB | Metadata | Optional | Domain descriptions |

**✅ Minimum Required**: ~22GB (InterPro + GOA + GO)  
**❌ NOT Required**: Protein sequences (unless doing local HMM scanning)

---

## 🔧 **Configuration Options**

### Taxonomic Scope (Choose Your Strategy)

**Option 1: Multi-species (Recommended)**
- **Pros**: Maximum statistical power, comprehensive coverage
- **Cons**: Large datasets, longer processing time
- **Use case**: Research applications, comprehensive analysis

**Option 2: Human-only**
- **Pros**: Smaller datasets, faster processing
- **Cons**: Reduced statistical power, limited coverage  
- **Use case**: Human-focused studies, resource constraints

```python
# Modify config/settings.py for human-only
# Filter datasets to human proteins only before processing
```

### Hardware Requirements

| Component | Minimum | Recommended | HPC Cluster |
|-----------|---------|-------------|-------------|
| **CPU Cores** | 4 | 8-16 | 32+ |
| **RAM** | 16GB | 32GB | 128GB+ |
| **Storage** | 100GB | 500GB | 1TB+ |
| **Runtime** | Days | Hours | Minutes |

---

## 🧪 **Development & Testing**

### Run Tests
```bash
# Run all tests
uv run pytest tests/ -v

# Run specific test categories  
uv run pytest tests/unit/ -v          # Unit tests
uv run pytest tests/e2e/ -v           # Integration tests
uv run pytest -m "not slow" -v        # Fast tests only

# Run with coverage
uv run pytest tests/ --cov=src --cov-report=html
```

### Code Quality
```bash
# Format code
uv run black src/ tests/
uv run isort src/ tests/

# Type checking
uv run mypy src/

# Linting
uv run ruff check src/ tests/
```

---

## 🏗️ **Architecture Overview**

```
dcGOspeed/
├── src/                              # Core implementation (5,000+ lines)
│   ├── data_acquisition.py           # Download system with progress tracking
│   ├── domain_scanning.py            # InterPro integration + supra-domains  
│   ├── statistical_inference.py      # Fisher's tests + FDR correction
│   ├── ontology_processor.py         # True Path Rule implementation
│   ├── database_manager.py           # SQLite backend with optimization
│   └── main_pipeline.py              # Orchestration + checkpoint/resume
├── config/settings.py                # Configuration management
├── tests/                            # Comprehensive test suite
├── docs/                            # Research papers & documentation
└── results/                         # Output directory
```

### Pipeline Stages

1. **Data Acquisition**: Downloads datasets with resume capability
2. **Domain Mapping**: Parses InterPro annotations and generates supra-domains
3. **Statistical Inference**: Fisher's exact tests with FDR correction
4. **Ontology Processing**: True Path Rule with optimal level filtering  
5. **Database Export**: SQLite storage with TSV export

---

## 📈 **Expected Output**

### Database Schema
```sql
-- Main results table
CREATE TABLE domain_annotations (
    domain_id TEXT,           -- e.g., 'PF00001'
    go_id TEXT,               -- e.g., 'GO:0005515'  
    fdr_q_value REAL,         -- Statistical significance (< 0.01)
    association_score REAL,   -- Strength score (1-100)
    annotation_type TEXT,     -- 'direct' or 'propagated'
    direct_source_term TEXT   -- Original inference source
);
```

### Example Results
```
domain_id    go_id        fdr_q_value  association_score  annotation_type
PF00001      GO:0005515   0.001        95.2              direct
PF00001      GO:0003674   0.001        95.2              propagated
PF00002      GO:0016740   0.005        78.4              direct
```

---

## ⚠️ **Important Notes**

### Performance Considerations
- **First run**: Downloads ~70GB, may take hours/days depending on connection
- **Subsequent runs**: Skips existing files, much faster
- **HPC recommended**: For production use with full datasets
- **Checkpoint/resume**: Long-running jobs can be safely interrupted and resumed

### Data Freshness
- **Automated downloads**: Pipeline fetches latest versions automatically
- **Update frequency**: Run monthly for latest annotations
- **Version tracking**: All downloads are logged with timestamps

### Troubleshooting
```bash
# Check disk space
df -h

# Monitor download progress
tail -f logs/dcgo_pipeline_*.log

# Resume interrupted pipeline
uv run python -m src.main_pipeline --resume --output-dir results/

# Force re-download if files corrupted
uv run python -m src.main_pipeline --force-download
```

---

## 📚 **Documentation**

- **[IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md)**: 200+ page technical specification
- **[CLAUDE.md](CLAUDE.md)**: Developer guidance for modifications
- **[PROJECT_STATUS.md](PROJECT_STATUS.md)**: Current completion status
- **[COMPLETION_SUMMARY.md](COMPLETION_SUMMARY.md)**: Achievement summary

---

## 🔬 **Scientific Validation**

This implementation faithfully reproduces the dcGO methodology described in:

1. **Fang et al. (2013)**: Original dcGO methodology and statistical framework
2. **Gene Ontology Consortium (2017)**: Ontology structure and annotation principles  
3. **InterPro Documentation**: Domain classification and annotation standards

**Key Features**:
- ✅ Fisher's exact tests with Benjamini-Hochberg FDR correction
- ✅ True Path Rule implementation with optimal level determination
- ✅ Supra-domain analysis for domain combination effects
- ✅ Scalable architecture for millions of proteins

---

## 🎯 **Use Cases**

### Research Applications
- **Protein function prediction**: Assign functions to uncharacterized proteins
- **Domain analysis**: Understand functional roles of protein domains
- **Evolutionary studies**: Track domain function conservation across species
- **Method development**: Base platform for algorithm improvements

### Educational Applications  
- **Bioinformatics training**: Learn statistical methods in computational biology
- **Data science**: Example of large-scale biological data processing
- **Statistical inference**: Practical application of multiple testing correction

---

## 🤝 **Contributing**

We welcome contributions! Areas for improvement:

1. **Additional ontologies**: Extend beyond Gene Ontology
2. **Visualization tools**: Web interface for results exploration
3. **Performance optimization**: GPU acceleration for large datasets
4. **Cloud deployment**: Docker/Kubernetes configurations

See [CLAUDE.md](CLAUDE.md) for developer guidance.

---

## 📄 **License**

MIT License - see [LICENSE](LICENSE) file for details.

---

## 🏆 **Project Status**

**✅ PRODUCTION READY**

- **Core Pipeline**: 5,000+ lines of tested Python code
- **Testing**: Comprehensive unit and integration tests
- **Documentation**: Extensive guides and API documentation  
- **Dependencies**: All resolved and working
- **Validation**: End-to-end workflow tested and verified

**Ready for immediate use in computational biology research!**

---

*Last Updated: 2025-10-08*  
*Status: ✅ Complete & Production-Ready*  
*Python Version: 3.12+*