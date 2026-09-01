# Future Work: Expanding Ontology Coverage for Human Proteins

> **Status (2026-08-31): the expansion described below is built.** The registry
> holds **35 `--ontology` keys**, and the full production matrix (63 cells:
> every layer × baseline, plus paper-parity for every hierarchical layer) ran
> clean on 2026-08-18. See [STATUS.md](STATUS.md) for the snapshot and
> [TODO.md](TODO.md) for what remains — chiefly evaluation and reproducibility,
> not new adapters. The sections below are kept as the design record.
>
> **Status (foundations landed).** Two enabling seams are now in place:
> (1) the whole path is **species-parameterised** (`--species` on download /
> extract / run), so non-human annotations no longer need code changes; and
> (2) annotations enter through the **`AnnotationSource` abstraction**
> (`src/annotation_source.py`) — the Fisher/FDR engine only sees
> `{protein → {term}}`, so a new ontology is a new `AnnotationSource` subclass
> plus, to propagate, an ontology-specific hierarchy — an OBO DAG for graph-based
> ontologies (GO), or none at all when the hierarchy is implicit (EC's numbering).
> `GAFAnnotationSource` is the reference (GO) implementation, and
> `ECAnnotationSource` (`src/ec_annotation_source.py`, `run_dcgo_human.py
> --ontology ec`) is the first non-GO ontology — Enzyme Commission, parsed from
> Expasy `enzyme.dat`, which is already UniProt-keyed so it needs no id mapping.
> EC also supports **True Path Rule propagation** (`--enable-true-path`) via
> `propagate_ec_annotations` / `ec_ancestors` — the EC hierarchy is implicit in
> the numbering, so it needs no OBO.
>
> **Strategy: UniProt is the protein universe → prefer UniProt-native terms.**
> Because the domain side (`protein2ipr`), GOA, and Expasy ENZYME are all keyed
> by UniProt accession, the cheapest annotations are the ones UniProt already
> carries per accession — they need *no* identifier mapping.
>
> **Status (the UniProt-native sweep is done).** Every vocabulary reachable that
> way has now been enumerated and, where usable, wired up.
> `scripts/survey_uniprot_ontologies.py` measured all ~150 `DR` databases in the
> human subset and sorted them by proteins-per-term
> (`docs/uniprot_ontology_survey.md`); the usable ones are registered in
> **`src/ontology_registry.py`**, which is now the single dispatch table for
> `--ontology` (source factory + hierarchy factory + required inputs). 35 keys:
> the UniProt-native set (`go`, `ec`, `reactome`, `keyword`, `disease`, `doid`,
> `orphanet`, `orphanet_doid`, `tcdb`, `merops`, `cazy`, `unipathway`,
> `complex`, `drugbank`, `pharos`, `condensate`, `subcellular`, `ligand`,
> `cofactor`, `rhea`, `xref`), the gene-keyed layers (`hpo`, `syngo`, #66), the
> model-organism phenotype layers (`mp`, `wbphenotype`, `wbbt`, `zfa`, `fbcv`,
> `fbbt`, #68) and wave 3 (`mondo`, `orphanet_mondo`, `efo`, `celltype`,
> `ncit`, `oncotree`, #71). Beyond `DR` lines, three
> layers curated into the entry *body* are now harvested too: subcellular
> location (`CC` prose matched against `subcell.txt`), ChEBI ligands (`FT
> /ligand_id`) and cofactors, and Rhea reactions. Deliberately **not** registered:
> ~1:1 accession mirrors (AlphaFoldDB/STRING/GeneCards/DisGeNET/KEGG-gene …,
> where every term has one protein and no test has power) and domain databases
> (Pfam/PANTHER/SUPFAM …, circular) — both still reachable via `xref` for
> deliberate experiments. Hierarchies landed for all of GO, EC, Reactome,
> keywords, subcellular locations, ChEBI (ligand/cofactor), TCDB, MEROPS and
> CAZy; `--enable-true-path` now *fails* rather than silently skipping for the
> flat ones.
>
> The Disease Ontology has since landed too, and *without* the backbone: DO
> carries the cross-ontology mapping itself (`xref: MIM:`, `xref: ORDO:`), so
> `--ontology doid` / `orphanet_doid` re-key UniProt's own disease
> cross-references onto DOID terms at parse time and then propagate up DO's
> `is_a` DAG (`src/disease_ontology.py`). Mondo now works identically
> (`--ontology mondo` / `orphanet_mondo`, #71).
>
> **The central identifier-mapping backbone (§3) turned out not to be needed.**
> HPO, SynGO and the model-organism layers re-key gene → UniProt accession at
> parse time using each source's own id-mapping data (Swiss-Prot `DR GeneID` /
> `DR HGNC` lines; MGI/WormBase/ZFIN/FlyBase mapping files, TrEMBL included),
> with the DOID layer's counted policy for unmapped and one-to-many ids. Also
> notable: the model-organism layers need **no orthology projection** — the
> dcGO trick is that the association is learned on the model organism's own
> proteins and domains are species-agnostic.
>
> **What is left is blocked on inputs, not code:** MAxO (ontology on disk, no
> annotation source — needs Monarch's `maxo-annotations`) and the gated
> sources (SNOMED CT and MedDRA licences, OMIM `genemap2.txt` registration).
> The full acquisition ledger is `data/ACQUISITION_MATRIX.md` (untracked).

## Objectives
- Extend dcGO beyond Gene Ontology (GO) annotations to cover a broader ontology landscape relevant to human protein function, disease, phenotypes, and enzymatic activity.
- Establish a unified integration framework so new ontology layers can be ingested, transformed, and queried through existing dcGO interfaces without bespoke code.

## Target Ontologies and Primary Data Sources

> Progress markers below: **[done]** = selectable via `--ontology` today;
> **[partial]** = a UniProt-native proxy is available but not the ontology
> itself; **[open]** = still needs the identifier-mapping backbone (§3).

1. **[done]** **Disease Ontology (DO)** — `--ontology doid` and
   `--ontology orphanet_doid` re-key UniProt's `DR MIM` (phenotype) and
   `DR Orphanet` cross-references onto DOID terms at parse time, using DO's own
   `xref: MIM:` / `xref: ORDO:` lines, and propagate up DO's `is_a` DAG. The
   raw-id layers (`disease`, `orphanet`) remain for comparison.
   - *Ontology*: Disease Ontology (OBO Foundry), pinned release + checksum in
     `config/settings.py` (`disease_ontology`).
   - *Implementation*: `src/disease_ontology.py`; the mapping policy for
     unmapped, one-to-many and obsolete ids is documented there and counted at
     parse time.
   - *Landed since*: Mondo as the unifying layer (`--ontology mondo` /
     `orphanet_mondo`, #71 — same mechanism, different OBO). *Still open*:
     annotation sources beyond UniProt's own cross-references (Monarch,
     DisGeNET).
2. **[done]** **Human Phenotype Ontology (HPO)** — `--ontology hpo` (#66):
   `genes_to_phenotype.txt` (NCBI-GeneID-keyed) re-keyed to UniProt via
   Swiss-Prot's own `DR GeneID` lines, propagated over `hp.obo`
   (`src/hpo_annotation_source.py`). SynGO landed in the same PR
   (`--ontology syngo`, HGNC-keyed).
   - *Open question*: the paper-parity configuration collapses HPO 996 → 38
     associations, driven by relative inference — see TODO.md P0.
3. **[done]** **Mammalian Phenotype (MP)** — `--ontology mp --species mouse`
   (#68): MGI gene-to-MP associations restricted to single-gene genotypes,
   re-keyed via MGI's own id-mapping files (TrEMBL included,
   `src/mgi_annotation_source.py`). No orthology projection — the association
   is learned on mouse proteins; domains are species-agnostic. The same PR
   added `wbphenotype`/`wbbt` (worm), `zfa` (zebrafish affected anatomy — ZP is
   not derivable from the on-disk files) and `fbcv`/`fbbt` (fly).
4. **[done]** **Enzyme Commission (EC)** — `--ontology ec`, plus the finer
   `--ontology rhea` (individual catalysed reactions from `CC CATALYTIC
   ACTIVITY`).
   - *Ontology*: IUBMB EC hierarchy.
   - *Annotations*:
     - UniProtKB enzyme annotations (EC numbers per protein, canonical source).
     - BRENDA (for supplemental enzyme function data; licensing check required).
5. **[done]** **Pathway/Process Ontologies** — `--ontology reactome` (with
   hierarchy) and `--ontology unipathway`. Pathway Commons and KEGG *pathways*
   remain open: UniProt's `DR KEGG` line is a gene id, not a pathway id.
   - Reactome pathways (hierarchical, stable IDs) with UniProt mappings.
   - Pathway Ontology (PW) via Pathway Commons (covers KEGG, Reactome, WikiPathways).
   - KEGG Orthology (KO) to pathway relationships (check licensing for redistribution).
6. **[done]** **Chemical/Drug Ontologies** — `--ontology ligand` and
   `--ontology cofactor` (ChEBI, with the ChEBI hierarchy) and `--ontology
   drugbank`. TTD remains open.
   - ChEBI for molecular functions/binding.
   - DrugBank or Therapeutic Target Database (TTD) for drug-protein associations (licensing dependent).

## Acquisition & Harmonization Strategy
1. **Ontology Harvesting Layer**
   - Extend existing ontology ingestion module to accept a registry of sources (URL, format, update cadence).
   - Support OBO, OWL, and JSON formats via `obonet`, `owlready`, and `rdflib` parsers.
   - Normalize ontology metadata (ID, label, synonyms, definition) into a shared graph schema.
2. **Annotation Harvesting Layer**
   - Define standard annotation schema: `{protein_id, subject_type, ontology_id, relation_type, evidence_code, source_db, reference, taxon, mapping_method, inferred_from}`.
     - `subject_type` / `relation_type`: keep direct, cross-referenced, and cross-species annotations distinguishable (e.g. a direct gene–phenotype assertion vs. a protein–disease cross-reference) rather than collapsing them into one shape.
     - `mapping_method` / `inferred_from`: record how each annotation was derived (e.g. ortholog projection from a mouse MP term), so inferred annotations stay auditable and separable from direct evidence — see §4 Versioning & Provenance.
   - Build source-specific adapters:
     - UniProt REST API & FTP for DO/HPO/EC/Reactome cross-references.
     - HPOA TSV parser mapping HGNC → UniProt (via precomputed mapping table from UniProt or Ensembl BioMart).
     - MGI annotations filtered through orthology table (Ensembl Compara, Alliance of Genome Resources).
     - Monarch API exports for unified disease/phenotype data.
     - Reactome Neo4j export or `Reactome ContentService` for pathways.
   - Capture evidence codes where available (e.g., HPO uses ECO terms).
3. **Identifier Mapping Backbone**
   - Maintain central mapping tables:
     - UniProt accession ↔ HGNC symbol/ID ↔ Ensembl Gene/Protein ↔ NCBI Gene.
     - For phenotypes/diseases, use cross-ontology mappings (Mondo ↔ DO, HPO ↔ MP bridging) to avoid duplication.
   - Automate updates using UniProt ID mapping API; store as versioned artifact.
4. **Versioning & Provenance**
   - Record source release versions and timestamps.
   - Store raw downloads under `data/raw/{source}/{version}`; processed outputs under `data/processed/{ontology}/{version}`.
   - Embed provenance metadata in final annotation tables for reproducibility.

## Integration into dcGO Framework
1. **Schema Extensions**
   - Generalize current GO-specific data model to `Ontology` and `Annotation` tables.
   - Introduce ontology type metadata (GO_BP, DO, HPO, etc.) to allow filtering and downstream analytics.
   - Support hierarchical propagation rules per ontology (e.g., propagate through `is_a`, `part_of`, custom relations defined in each ontology).
2. **Pipeline Refactoring**
   - Parameterize ETL pipeline so new ontologies can be registered via configuration (YAML/JSON) specifying:
     - download URL(s)
     - parser class
     - propagation rules
     - evidence filters
   - Implement incremental update mechanism to avoid reprocessing stable releases.
3. **APIs and Data Access**
   - Update Python API to expose ontology-agnostic query functions (`get_annotations(protein, ontology_type=None)`).
   - Ensure existing GO workflows continue to function (regression tests).
   - Add utilities for cross-ontology queries (e.g., fetch HPO terms for a disease mapped to protein via DO).
4. **Testing & Validation**
   - Unit tests for each ingestion adapter using fixture downloads.
   - Integration tests verifying sample proteins produce expected combined annotation sets.
   - Consistency checks: verify ontology IDs resolve to labels; ensure no orphan annotations.

## Operational Considerations
- **Licensing**: Review redistribution terms for KEGG, BRENDA, DrugBank; may need user-provided credentials.
- **Update Cadence**: Implement monthly cron to refresh ontologies/annotations, with change logs.
- **Performance**: Cache intermediate graph representations; consider using graph database (e.g., Neo4j) if relational model becomes limiting.
- **Community Alignment**: Track OBO Foundry conventions; align evidence codes with ECO to facilitate interoperability.

## Milestones
1. **[done]** Prototype ingestion for DO and HPO (MVP multi-ontology support).
2. **[done]** Integrate MP — via the model organism's own proteins rather than
   orthology mapping; cross-species pipeline validated by the production
   matrix.
3. **[done]** Add EC and Reactome pathway layers.
4. **[partial]** Configuration-driven onboarding: the registry
   (`src/ontology_registry.py`) makes a new ontology one declarative entry;
   YAML/JSON-driven onboarding without code remains open.
5. **[open]** Monitoring dashboards summarizing ontology coverage per release.

## Open Questions
- Best approach for harmonizing conflicting disease ontologies (DO vs. Mondo vs. OMIM)? Evaluate adopting Mondo as unifying layer.
- How to handle evidence conflicts or multiple evidence codes per annotation? Consider storing as structured list.
- Evaluate need for ontology-specific propagation rules beyond `is_a` (e.g., HPO `part_of`, Reactome `hasEvent`).
- Determine storage solution for large pathway/phenotype networks (flat files vs. graph DB).


## Beyond ontologies: ranking what the method uniquely finds

The breadth work above multiplies *what* domains can be associated with. The
complementary question — which associations are worth reading — is answered by
the **surprise score** (`src/surprise_score.py`,
`scripts/rank_surprising_associations.py`, `SURPRISE_SCORE.md`), which ranks
supra-domain associations by how much the combination beats what its
constituents already predict, discounting redundant InterPro signatures and
curated knowledge. Each new ontology gets this ranking for free; only the
novelty factor needs a per-ontology curated reference (InterPro2GO exists for
GO — an equivalent for EC/Reactome would be a small, self-contained addition).
