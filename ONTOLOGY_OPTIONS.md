# Ontology and Annotation Quality Options

## 1. GOA Evidence Code Filtering ✅ IMPLEMENTED

### What are Evidence Codes?

GOA (Gene Ontology Annotation) files include **evidence codes** that indicate how an annotation was derived. These range from high-confidence experimental evidence to computationally inferred or automatically assigned annotations.

### Evidence Code Categories

| Category | Codes | Description | Reliability |
|----------|-------|-------------|-------------|
| **Experimental** | EXP, IDA, IPI, IMP, IGI, IEP | Direct experimental evidence | ⭐⭐⭐ Highest |
| **Computational** | ISS, ISO, ISA, ISM, IGC, IBA, IBD, IKR, IRD, RCA | Sequence/structure similarity | ⭐⭐ High |
| **Author Statement** | TAS, NAS | Traceable/Non-traceable author statements | ⭐⭐ High |
| **Curator Statement** | IC, ND | Inferred by curator, No biological data | ⭐ Medium |
| **Electronic** | IEA | Inferred from electronic annotation (automated) | ⚠️ Lower confidence |

### Evidence Code Details

#### Experimental Evidence (Highest Confidence)
- **EXP** - Inferred from Experiment
- **IDA** - Inferred from Direct Assay (e.g., enzyme assay, Western blot)
- **IPI** - Inferred from Physical Interaction (e.g., co-immunoprecipitation)
- **IMP** - Inferred from Mutant Phenotype (e.g., gene knockout)
- **IGI** - Inferred from Genetic Interaction (e.g., synthetic lethality)
- **IEP** - Inferred from Expression Pattern (e.g., RNA-seq, microarray)

#### Computational Analysis
- **ISS** - Inferred from Sequence or structural Similarity
- **ISO** - Inferred from Sequence Orthology
- **ISA** - Inferred from Sequence Alignment
- **ISM** - Inferred from Sequence Model (e.g., HMM)
- **IGC** - Inferred from Genomic Context
- **IBA** - Inferred from Biological aspect of Ancestor
- **IBD** - Inferred from Biological aspect of Descendant
- **IKR** - Inferred from Key Residues
- **IRD** - Inferred from Rapid Divergence
- **RCA** - Inferred from Reviewed Computational Analysis

#### Author/Curator Statements
- **TAS** - Traceable Author Statement (cited in publication)
- **NAS** - Non-traceable Author Statement
- **IC** - Inferred by Curator
- **ND** - No biological Data available

#### Electronic Annotation (Lowest Confidence)
- **IEA** - Inferred from Electronic Annotation (automated, no human review)

### Configuration Options

The pipeline now supports three evidence filter presets:

```python
# In config/settings.py
evidence_filter: str = 'manual'  # Default

# Options:
# 'all'          - Include ALL evidence codes (including IEA)
# 'manual'       - Exclude IEA (only manually curated annotations)
# 'experimental' - Only experimental evidence (highest confidence)
```

### Usage

#### Via Configuration File
```python
# config/settings.py
evidence_filter: str = 'experimental'  # Use only experimental evidence
```

#### Via Command Line
```bash
# Use manual curation only (default)
python -m src.main_pipeline --evidence-filter manual

# Include all annotations (including IEA)
python -m src.main_pipeline --evidence-filter all

# Use only experimental evidence
python -m src.main_pipeline --evidence-filter experimental
```

### Impact on Results

| Filter | Annotations | Proteins | Quality | Coverage |
|--------|------------|----------|---------|----------|
| `all` | ~100% | ~100% | Mixed | Maximum |
| `manual` | ~60-70% | ~95% | High | Good |
| `experimental` | ~20-30% | ~70% | Highest | Limited |

**Recommendation**: Start with `'manual'` (default) to exclude IEA annotations while maintaining good coverage. Use `'experimental'` for high-confidence analyses.

### Implementation

The pipeline uses the new `GOAParser` class in `src/goa_parser.py`:

```python
from src.goa_parser import parse_goa_human

protein_go_map = parse_goa_human(
    gaf_path="data/raw/goa_annotations/goa_human.gaf.gz",
    evidence_filter='manual',  # or 'all', 'experimental'
    aspects={'P', 'F', 'C'}    # GO aspects to include
)
```

The parser automatically:
- Filters by evidence codes
- Excludes NOT qualifiers (negative annotations)
- Reports statistics on evidence code distribution
- Logs filtering decisions

---

## 2. Other Ontologies (HPO, DO, MP) ⚠️ NOT CURRENTLY SUPPORTED

### Available from UniProt/GOA?

**Short answer**: No, UniProt/GOA only provides Gene Ontology (GO) annotations.

### Other Phenotype/Disease Ontologies

| Ontology | Name | Focus | Source |
|----------|------|-------|--------|
| **HPO** | Human Phenotype Ontology | Human phenotypes | Monarch Initiative |
| **DO** | Disease Ontology | Human diseases | Disease Ontology Project |
| **MP** | Mammalian Phenotype Ontology | Mouse phenotypes | MGI |
| **OMIM** | Online Mendelian Inheritance in Man | Human genetic disorders | OMIM |

### Why Stick with GO for Now?

1. **Data Availability**: GOA provides comprehensive, well-curated GO annotations for all UniProt proteins
2. **InterPro Integration**: InterPro domains already have GO term associations
3. **Established Methodology**: dcGO was originally designed for GO terms
4. **Statistical Power**: GO has the most annotations, providing better statistical inference

### Future Extension Possibilities

While not currently implemented, the pipeline could theoretically be extended to support other ontologies:

#### HPO (Human Phenotype Ontology)
- **Source**: https://hpo.jax.org/
- **Data**: Gene-phenotype associations
- **Format**: Similar to GAF, could use same methodology
- **Use case**: Link domains to human phenotypes

#### Disease Ontology (DO)
- **Source**: https://disease-ontology.org/
- **Data**: Gene-disease associations
- **Challenge**: Less comprehensive domain-disease data

#### Mammalian Phenotype (MP)
- **Source**: http://www.informatics.jax.org/
- **Data**: Mouse gene-phenotype associations
- **Challenge**: Need mouse protein-domain mappings

### Potential Implementation Path

If you wanted to add HPO support in the future:

1. Download HPO annotations (e.g., from Monarch Initiative)
2. Parse HPO phenotype-gene associations
3. Use the same dcGO methodology (Fisher's exact test)
4. Generate domain-phenotype associations

However, this would require:
- HPO annotation files in a parseable format
- HPO ontology structure (OBO format)
- Validation that domain-phenotype associations make biological sense

### Recommendation

**Stick with GO for now** because:
- ✅ Comprehensive data availability
- ✅ Well-established methodology
- ✅ Strong statistical power
- ✅ Direct InterPro-GO associations already exist
- ✅ Widely used and validated

Once you have a working dcGO pipeline for GO terms, you could consider extending to other ontologies as a separate project.

---

## Summary

| Feature | Status | Configuration |
|---------|--------|---------------|
| GO evidence filtering | ✅ Implemented | `evidence_filter = 'manual'` |
| GO aspect filtering | ✅ Implemented | `aspects = {'P', 'F', 'C'}` |
| HPO support | ❌ Not available | - |
| DO support | ❌ Not available | - |
| MP support | ❌ Not available | - |

The pipeline now provides robust GO annotation quality control through evidence code filtering, which is the most important quality metric for functional annotations. Other ontologies would require separate data sources and validation.
