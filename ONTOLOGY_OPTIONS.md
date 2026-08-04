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

## 2. Other Ontologies ✅ IMPLEMENTED (18 beyond GO)

> This section previously said UniProt/GOA supplies only GO. That is wrong: GOA
> supplies only GO, but the **UniProt flat file** carries a dozen more
> vocabularies per accession, and the entry body carries several more.

`run_dcgo_human.py --ontology <name>` selects the annotation layer; the
Fisher/FDR engine is unchanged, because everything enters through the
`AnnotationSource` seam as `{protein → {term}}`. The dispatch table lives in
`src/ontology_registry.py`; `--help` prints the current list.

| `--ontology` | Terms | Source | True Path hierarchy |
| --- | --- | --- | --- |
| `go` | Gene Ontology | GOA GAF (evidence-filtered) | `go-basic.obo` |
| `ec` | Enzyme Commission | Expasy `enzyme.dat` | implicit in the number |
| `reactome` | pathways | `DR Reactome` | `ReactomePathwaysRelation.txt` |
| `keyword` | UniProt keywords | `KW` lines | `keywlist.txt` |
| `subcellular` | `SL-` locations | `CC SUBCELLULAR LOCATION` | `subcell.txt` |
| `ligand` | ChEBI ligands | `FT …/ligand_id` | `chebi_lite.obo` |
| `cofactor` | ChEBI cofactors | `CC COFACTOR` | `chebi_lite.obo` |
| `rhea` | Rhea reactions | `CC CATALYTIC ACTIVITY` | — |
| `tcdb` | transporter classes | `DR TCDB` | implicit |
| `merops` | peptidase families | `DR MEROPS` | implicit |
| `cazy` | CAZy families | `DR CAZy` | implicit |
| `disease` | OMIM phenotypes | `DR MIM` (phenotype) | — |
| `doid` | Disease Ontology | `DR MIM`, re-keyed via `doid.obo` xrefs | `doid.obo` |
| `orphanet` | rare diseases | `DR Orphanet` | — |
| `orphanet_doid` | Disease Ontology | `DR Orphanet`, re-keyed via `doid.obo` xrefs | `doid.obo` |
| `unipathway` | metabolic pathways | `DR UniPathway` | — |
| `complex` | protein complexes | `DR ComplexPortal` | — |
| `drugbank` | drugs | `DR DrugBank` | — |
| `pharos` | target development level | `DR Pharos` (3rd field) | — |
| `condensate` | condensates | `DR CD-CODE` (3rd field) | — |
| `xref` | any other DR database | `--xref-db NAME` | — |

### Evidence codes

Only GO has them. The UniProt-native layers are curated cross-references and
comments, so there is no IEA/experimental distinction to filter — which is also
why GO is taken from GOA (with `--evidence-filter`) rather than from UniProt's
own `DR GO` lines, which include IEA.

### What is still missing, and why

HPO and MP are **not** UniProt-keyed: HPO annotates HGNC genes and MGI annotates
mouse genes. Those need the identifier mapping backbone described in
`FUTURE_WORK.md` §3 — the only remaining reason to build it.

The Disease Ontology *is* now reachable (`doid`, `orphanet_doid`), but by a
different route: DO carries the cross-ontology mapping itself, as `xref: MIM:`
and `xref: ORDO:` lines, so re-keying UniProt's own disease cross-references at
parse time needs no protein-identifier mapping at all — only a term-space
translation (`src/disease_ontology.py`). Mondo, which dcGO 2023 switched to,
would work the same way. Note that this is a mapping *through* OMIM/Orphanet, so
it inherits their coverage: a disease UniProt does not cross-reference, or that
DO does not cross-reference back, is invisible to it.

Deliberately excluded, though reachable through `xref`:

* **1:1 accession mirrors** (AlphaFoldDB, STRING, GeneCards, DisGeNET, KEGG gene
  ids, …). Every "term" has one protein, so no contingency table has signal.
* **Domain databases** (Pfam, PANTHER, SUPFAM, CDD, PROSITE, …). Associating
  InterPro domains with domain signatures is circular.

The measurements behind those calls are in `docs/uniprot_ontology_survey.md`.

---

## Summary

| Feature | Status | Configuration |
|---------|--------|---------------|
| GO evidence filtering | ✅ Implemented | `evidence_filter = 'manual'` |
| GO aspect filtering | ✅ Implemented | `aspects = {'P', 'F', 'C'}` |
| 18 non-GO ontologies | ✅ Implemented | `--ontology reactome\|ligand\|subcellular\|…` |
| HPO / MP / DO-proper support | ❌ Needs the id-mapping backbone | see `FUTURE_WORK.md` §3 |

Evidence-code filtering remains the main quality control for GO. For the
UniProt-native layers the equivalent control is *which layer you choose*: the
`CC`/`FT` layers are manually curated with ECO evidence in the flat file, while
`DR` cross-references inherit whatever the source database asserts.
