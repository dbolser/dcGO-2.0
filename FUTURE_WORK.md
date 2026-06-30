# Future Work: Expanding Ontology Coverage for Human Proteins

## Objectives
- Extend dcGO beyond Gene Ontology (GO) annotations to cover a broader ontology landscape relevant to human protein function, disease, phenotypes, and enzymatic activity.
- Establish a unified integration framework so new ontology layers can be ingested, transformed, and queried through existing dcGO interfaces without bespoke code.

## Target Ontologies and Primary Data Sources
1. **Disease Ontology (DO)**
   - *Ontology*: Disease Ontology (OBO Foundry).
   - *Annotations*: 
     - Monarch Initiative (disease-phenotype/protein associations).
     - UniProtKB "Disease" section with DO and Mondo cross-references.
     - DisGeNET gene-disease associations (map to proteins via UniProt/Ensembl).
2. **Human Phenotype Ontology (HPO)**
   - *Ontology*: `HPO OBO/OWL`.
   - *Annotations*:
     - HPOA gene-phenotype annotation files (official releases) keyed by HGNC, OMIM, Ensembl.
     - Monarch Initiative phenogrid exports (for cross-species phenotype integration).
3. **Mammalian Phenotype (MP)**
   - *Ontology*: MGI Mammalian Phenotype.
   - *Annotations*:
     - Mouse Genome Informatics (MGI) gene-to-MP associations.
     - Monarch cross-species mapping to human orthologs (via HGNC/Ensembl).
4. **Enzyme Commission (EC)**
   - *Ontology*: IUBMB EC hierarchy.
   - *Annotations*:
     - UniProtKB enzyme annotations (EC numbers per protein, canonical source).
     - BRENDA (for supplemental enzyme function data; licensing check required).
5. **Pathway/Process Ontologies**
   - Reactome pathways (hierarchical, stable IDs) with UniProt mappings.
   - Pathway Ontology (PW) via Pathway Commons (covers KEGG, Reactome, WikiPathways).
   - KEGG Orthology (KO) to pathway relationships (check licensing for redistribution).
6. **Chemical/Drug Ontologies (stretch goal)**
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
1. Prototype ingestion for DO and HPO with UniProt + HPOA data (MVP multi-ontology support).
2. Integrate MP via orthology mapping; validate cross-species pipeline.
3. Add EC and Reactome pathway layers to cover enzymatic and pathway context.
4. Implement configuration-driven onboarding for additional ontologies (KEGG, ChEBI, etc.).
5. Deploy monitoring dashboards summarizing ontology coverage per release.

## Open Questions
- Best approach for harmonizing conflicting disease ontologies (DO vs. Mondo vs. OMIM)? Evaluate adopting Mondo as unifying layer.
- How to handle evidence conflicts or multiple evidence codes per annotation? Consider storing as structured list.
- Evaluate need for ontology-specific propagation rules beyond `is_a` (e.g., HPO `part_of`, Reactome `hasEvent`).
- Determine storage solution for large pathway/phenotype networks (flat files vs. graph DB).

