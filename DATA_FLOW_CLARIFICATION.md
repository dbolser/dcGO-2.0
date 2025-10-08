# dcGO Data Flow Clarification

## 🚨 **IMPORTANT: Current Pipeline Architecture Issue**

The current implementation has a **mismatch** between the configuration (which correctly identifies pre-computed domain annotations) and the main pipeline (which still tries to run InterProScan on protein sequences).

---

## ✅ **CORRECT dcGO Data Flow**

### **What dcGO Actually Needs:**

```
Input 1: Domain Annotations    Input 2: GO Annotations    Input 3: Ontology DAG
     ↓                              ↓                           ↓
   protein2ipr.dat.gz         goa_uniprot_all.gaf.gz      go-basic.obo
   (P12345 → PF00001)         (P12345 → GO:0005515)      (GO hierarchy)
     ↓                              ↓                           ↓
        ↘                          ↙                           ↓
         Correspondence Matrix                            True Path Rule
         (Domain × GO co-occurrence)                      (Ontology logic)
                     ↓                                          ↓
              Fisher's Exact Tests                    Optimal Level Filter
              + FDR Correction                        + Propagation
                     ↓                                          ↓
                dcGO Database (Domain → GO associations)
```

### **What We DON'T Need:**
- ❌ Protein sequences (FASTA files)
- ❌ HMM profiles (unless doing local scanning)
- ❌ InterProScan execution (unless custom proteins)

---

## 🔧 **Required Configuration Updates**

### **1. Data Sources (config/settings.py)**

```python
# REQUIRED (The 3 Essential Inputs)
'interpro_mappings': DataSource(
    name='interpro_mappings',
    url='https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz',
    description='Pre-computed InterPro domain mappings for UniProt proteins',
    required=True  # ✅ INPUT 1: Domain annotations
),
'goa_annotations': DataSource(
    name='goa_annotations', 
    url='ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz',
    description='Gene Ontology Annotation (GOA) database',
    required=True  # ✅ INPUT 2: GO annotations
),
'go_ontology': DataSource(
    name='go_ontology',
    url='http://current.geneontology.org/ontology/go-basic.obo', 
    description='Gene Ontology basic ontology file',
    required=True  # ✅ INPUT 3: Ontology DAG
),

# OPTIONAL (Only for local HMM scanning or alternatives)
'uniprot_sprot': DataSource(
    name='uniprot_sprot',
    url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
    description='UniProt Swiss-Prot protein sequences (only needed for local HMM scanning)',
    required=False  # ❌ NOT needed for standard dcGO
),
'uniprot_trembl': DataSource(
    name='uniprot_trembl',
    url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz', 
    description='UniProt TrEMBL protein sequences (only needed for local HMM scanning)',
    required=False  # ❌ NOT needed for standard dcGO
),
```

### **2. Pipeline Stages (src/main_pipeline.py)**

The domain scanning stage should be replaced with domain annotation parsing:

```python
# CURRENT (INCORRECT)
async def run_domain_scanning(self, parameters, data_files):
    # Tries to run InterProScan on sequences - WRONG!
    sequence_file = data_files.get('uniprot_sprot')
    # ... InterProScan execution

# CORRECT (UPDATED)  
async def run_domain_annotation_parsing(self, parameters, data_files):
    # Parse pre-computed domain annotations - RIGHT!
    interpro_file = data_files.get('interpro_mappings')
    # ... Parse protein2ipr.dat.gz format
```

---

## 📊 **Corrected Data Requirements**

| File | Size | dcGO Role | Source | Required |
|------|------|-----------|--------|----------|
| `protein2ipr.dat.gz` | ~20GB | **Input 1**: Domain annotations | InterPro | ✅ Essential |
| `goa_uniprot_all.gaf.gz` | ~2GB | **Input 2**: GO annotations | GOA | ✅ Essential |
| `go-basic.obo` | ~50MB | **Input 3**: Ontology DAG | GO Consortium | ✅ Essential |
| `uniprot_sprot.fasta.gz` | ~500MB | Not used in dcGO | UniProt | ❌ Optional |
| `uniprot_trembl.fasta.gz` | ~50GB | Not used in dcGO | UniProt | ❌ Optional |

**Total Essential Data**: ~22GB (not ~70GB!)

---

## 🔄 **Implementation Updates Needed**

### **1. Domain Annotation Parser**
Create `src/domain_annotation_parser.py`:
- Parse `protein2ipr.dat.gz` format
- Generate supra-domains from domain architectures
- Output: `protein_domain_map = {'P12345': ['PF00001', 'PF00002', 'PF00001,PF00002']}`

### **2. GO Annotation Parser** 
Update `src/go_annotation_parser.py`:
- Parse GAF 2.2 format from `goa_uniprot_all.gaf.gz`
- Output: `protein_go_map = {'P12345': {'GO:0005515', 'GO:0008150'}}`

### **3. Main Pipeline**
Update `src/main_pipeline.py`:
- Replace `run_domain_scanning()` with `run_domain_annotation_parsing()`
- Remove InterProScan execution logic
- Focus on parsing pre-computed annotations

---

## ✅ **Simplified Workflow**

```bash
# 1. Download the 3 essential files (~22GB)
uv run python -c "
from src.data_acquisition import DataAcquisition
from config.settings import Config
da = DataAcquisition(Config())
da.download_specific_dataset('interpro_mappings')  # ~20GB
da.download_specific_dataset('goa_annotations')    # ~2GB  
da.download_specific_dataset('go_ontology')        # ~50MB
"

# 2. Run dcGO analysis (no HMM scanning needed)
uv run python -m src.main_pipeline --num-cores 8 --output-dir results/
```

---

## 🎯 **Bottom Line**

**The dcGO methodology does NOT require protein sequences or HMM scanning when using pre-computed domain annotations from InterPro. The current pipeline architecture needs to be updated to reflect this reality.**

**Essential Downloads**: 22GB  
**Optional Downloads**: 50GB+  
**Current Implementation**: Mixed (correct config, incorrect pipeline execution)**