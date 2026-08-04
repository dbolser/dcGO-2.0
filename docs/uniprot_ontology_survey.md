# What ontologies are reachable from the UniProt downloads?

UniProt accessions are this pipeline's protein universe — `protein2ipr` domains,
GOA and Expasy ENZYME are all keyed by them — so any vocabulary UniProt already
carries per accession can be tested for domain enrichment with **no identifier
mapping at all**. This note is the survey behind
[`src/ontology_registry.py`](../src/ontology_registry.py): what the Swiss-Prot
flat file actually contains, which parts of it are usable as ontologies, and
which are not.

Survey basis: `uniprot_sprot.dat.gz` (2026-07 release), 575,503 entries, of which
**20,431 are human** (`OX NCBI_TaxID=9606`). Counts below are human-only.

## The four annotation layers in the flat file

| Layer | Line | Example |
| --- | --- | --- |
| Cross-references | `DR` | `DR   Reactome; R-HSA-71384; Ethanol oxidation.` |
| Keywords | `KW` | `KW   Metal-binding; NAD; Oxidoreductase; Zinc.` |
| Curated comments | `CC   -!- TOPIC:` | `CC   -!- SUBCELLULAR LOCATION: Cytoplasm.` |
| Feature qualifiers | `FT   .../ligand_id=` | `FT  /ligand_id="ChEBI:CHEBI:29105"` |

The first two were already wired up; the last two are new (see
`src/uniprot_annotation_source.py`).

## Sorting the ~150 DR databases

The decisive statistic is **proteins per distinct term**. A vocabulary has many
proteins sharing each term; a database that mirrors the protein 1:1 has a ratio
of ~1 and cannot support an enrichment test at all — every "term" would have a
single protein, so no contingency table has signal.

**(a) Vocabularies — registered as first-class `--ontology` values**

Counts from `docs/dr_survey.tsv` (human, 2026-07 release):

| Database | human proteins | distinct terms | proteins/term | `--ontology` |
| --- | ---: | ---: | ---: | --- |
| `KW` (keywords) | 20,431 | 720 | 28.4 | `keyword` |
| Pharos (development level) | 20,191 | 4 | ~5,000 | `pharos` |
| Reactome | 11,392 | 2,284 | 5.0 | `reactome` |
| MIM, phenotype-typed | 5,044 | 6,920 | 0.7 | `disease` |
| Orphanet | 4,482 | 4,119 | 1.1 | `orphanet` |
| CD-CODE (condensates) | 3,752 | 238 | 15.8 | `condensate` |
| ComplexPortal | 3,583 | 2,419 | 1.5 | `complex` |
| DrugBank | 3,388 | 8,486 | 0.4 | `drugbank` |
| TCDB | 2,347 | 1,952 | 1.2 | `tcdb` |
| UniPathway | 1,354 | 203 | 6.7 | `unipathway` |
| MEROPS | 858 | 891 | 1.0 | `merops` |
| CAZy | 230 | 79 | 2.9 | `cazy` |

Two of these need a twist: `Pharos` and `CD-CODE` key the `DR` line by the
protein and put the vocabulary in the **third** field
(`DR   Pharos; P31946; Tbio.`), so those sources take the term from that field
(`term_from_id_type`, also exposed as `--xref-term-from-type`). Their term counts
above are that vocabulary (Pharos: Tbio 11,968 / Tdark 5,549 / Tchem 1,970 /
Tclin 704; CD-CODE: 238 condensate names), not the raw `DR` ids the survey TSV
counts.

The low-ratio disease and drug resources (MIM, Orphanet, DrugBank) sit at the
edge: many terms have a single protein and contribute nothing, but the shared
ones — a disease caused by several paralogues, a drug hitting a family — are
exactly the informative rows. They are registered, with the caveat that their
effective term space is far smaller than the raw count.

**(b) 1:1 accession mirrors — deliberately not registered**

AlphaFoldDB, STRING, PaxDb, PeptideAtlas, PhosphoSitePlus, iPTMnet, GeneCards,
BioMuta, CTD, DisGeNET, GeneID, KEGG (`hsa:7529` is a *gene*, not a pathway),
MassIVE, jPOST, Bgee, ExpressionAtlas, OpenTargets, SignaLink, PathwayCommons,
SMR, DMDM, Antibodypedia, ChiTaRS, GenomeRNAi, BioGRID-ORCS, DNASU, RNAct, PRO,
PAN-GO, AGR, HPA, ClinPGx, Agora, VEuPathDB, MalaCards, ChEMBL, BindingDB,
CORUM, GuidetoPHARMACOLOGY, DrugCentral … — all ~1.0 proteins per term. They are
still reachable through the escape hatch (`--ontology xref --xref-db NAME`) for
deliberate experiments, but they are not offered as ontologies.

Borderline: the orthology resources (eggNOG 2.0, OrthoDB 1.9, HOGENOM 1.4,
GeneTree 1.3) *are* genuine many-to-one groupings, but an orthologous group is
defined by sequence similarity, so associating it with domain content is close
to circular. Excluded by default for that reason, not for lack of signal.

**(c) Domain databases — excluded as circular**

InterPro, Pfam, PANTHER, SUPFAM, Gene3D, FunFam, CDD, PROSITE, SMART, PRINTS,
PIRSF, NCBIfam, HAMAP, SFLD. Associating InterPro domains with domain-database
signatures measures the signature redundancy of InterPro, not biology. (It is
occasionally useful as a positive control — hence, again, `xref`.)

**(d) Already covered by a dedicated source**

`DR GO` duplicates GOA (which additionally carries evidence codes, so it stays
the GO source), and `DR BRENDA` carries EC numbers already obtained from Expasy
ENZYME.

## Layers curated into the entry body

These are not `DR` lines, so each needs its own extractor. They are the most
interesting additions here, because they are *curated, evidence-tagged*
annotation rather than a pointer to another database.

| `--ontology` | Source in the flat file | Vocabulary | Hierarchy |
| --- | --- | --- | --- |
| `subcellular` | `CC   -!- SUBCELLULAR LOCATION` | `subcell.txt` `SL-` terms | `HI`/`HP` lines |
| `ligand` | `FT   …/ligand_id="ChEBI:…"` | ChEBI | `chebi_lite.obo` |
| `cofactor` | `CC   -!- COFACTOR` `Xref=ChEBI:…` | ChEBI | `chebi_lite.obo` |
| `rhea` | `CC   -!- CATALYTIC ACTIVITY` `Xref=Rhea:…` | Rhea reactions | none |

`subcellular` is prose parsed against a controlled vocabulary: the comment is a
run of `.`-terminated statements, optionally prefixed with `[Isoform n]:`,
carrying `;`-separated topology qualifiers and a free-text `Note=` tail. Each
piece is matched against the `SL` content lines, entry names and synonyms in
`subcell.txt`; unmatched strings are counted and reported rather than silently
dropped.

`ligand` deserves emphasis: binding sites are annotated at the *residue* level,
which is precisely the scale of a domain, so a domain→ligand association is the
most mechanistically direct claim in this whole set.

## Hierarchies, and where they come from

The True Path Rule needs ancestors. Four different mechanisms cover everything
registered:

| Mechanism | Ontologies |
| --- | --- |
| OBO graph | GO (`go-basic.obo`), ChEBI (`chebi_lite.obo`) |
| Implicit in the id | EC (`1.1.1.1`), TCDB (`8.A.98.1.10`), MEROPS (`S01.151` → `S01` → `S`), CAZy (`GT32` → `GT`) |
| Companion hierarchy file | Reactome, UniProt keywords, subcellular locations |
| None | OMIM, Orphanet, Rhea, DrugBank, ComplexPortal, UniPathway, Pharos, CD-CODE, `xref` |

`--enable-true-path` now **fails loudly** for the last group instead of running
without propagation.

## Coverage on the human domain set

Proteins here are those among the 18,908 human proteins carrying InterPro
domains that the layer annotates — the rows that can actually enter a
contingency table. "Ancestor terms added" counts term ids that True Path
propagation introduces but that no protein is annotated with directly; a zero
means the vocabulary is either flat or already used at every level (as with
keywords, where 720 terms in use already cover their own parents), not that
propagation does nothing per protein.

| `--ontology` | proteins (with domains) | terms | ancestor terms added |
| --- | ---: | ---: | ---: |
| keyword | 18,859 | 720 | 0 |
| pharos | 18,744 | 4 | — |
| subcellular | 16,750 | 261 | 8 |
| disease | 5,029 | 6,904 | — |
| ligand | 4,627 | 448 | 849 |
| orphanet | 4,461 | 4,103 | — |
| rhea | 4,135 | 10,728 | — |
| condensate | 3,711 | 316 | — |
| complex | 3,575 | 2,419 | — |
| drugbank | 3,381 | 8,483 | — |
| tcdb | 2,322 | 1,931 | 1,232 |
| cofactor | 1,801 | 46 | 232 |
| unipathway | 1,353 | 203 | — |
| merops | 843 | 880 | 114 |
| cazy | 230 | 79 | 3 |

(`—` = no hierarchy, so no propagation. `go`, `ec` and `reactome` are omitted:
they were already covered before this sweep.)

Two things worth reading off this table. **Term/protein ratio decides
statistical power, not protein count**: `rhea` annotates 4,135 proteins with
10,728 reactions, so most reactions have one or two proteins and only the
well-populated ones will yield associations; `pharos` is the opposite extreme,
18,744 proteins over 4 classes. **The ChEBI layers propagate the hardest** —
`ligand` gains 849 ancestor terms over 448 direct ones, because ChEBI's chemical
hierarchy is deep, so True Path there produces many generic terms
("metal cation", "molecular entity") alongside the specific ones.

## Do these layers actually predict anything?

Coverage is not usefulness. `validation/temporal_breadth.py` trains each layer on
an archived 2021 Swiss-Prot release and scores it against 2026 curation
(`VALIDATION_PLAN.md` §2, breadth subsection). Enrichment is the hit rate over
the term's own acquisition rate. **Corrected 2026-08-04**: the first published
version of this table scored association sets built before `restrict_to_universe`
(#26) and understated every layer — see the banner in `VALIDATION_PLAN.md` §2.

| `--ontology` | enrichment (95% CI) | verdict |
| --- | --- | --- |
| `go` (anchor) | 11.5× [11.1, 12.0] | strongest |
| `reactome` | 11.4× [8.9, 14.6] | indistinguishable from GO at this power (1,430 hits vs GO's 106,224) |
| `subcellular` | 3.7× [3.5, 3.9] | real |
| `cofactor` | 3.5× [2.1, 4.5] | real, term-specific |
| `keyword` | 3.4× [3.2, 3.6] | real; its 0.44% base rate is the highest of any layer, which caps the ratio |
| `complex` | no demonstrated signal | 2 hits / 9,904; CI [0.00, 130.93] includes zero |
| `disease` | undefined (0 hits / 369) | too sparse to test — 6,904 raw OMIM ids over 5,029 proteins. Superseded by `doid`, which re-keys the same curation onto the Disease Ontology (`src/disease_ontology.py`); re-running the breadth test on `doid` is outstanding |
| `ligand` | untestable at this split | `/ligand_id` postdates 2021 (see below) |

The layers not listed were not trained at t0. Note that `ligand` is *entirely
post-2022 annotation*: UniProt used free-text `/note="ATP"` on binding sites
until the structured ChEBI qualifier arrived, so that layer rests on recent, and
still-growing, curation.

## Reproducing this survey

```bash
# distinct terms and proteins per DR database, human only
uv run python scripts/survey_uniprot_ontologies.py --output docs/dr_survey.tsv
```
