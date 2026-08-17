# dcGO Pipeline

> **Scope & safety — safe to review.** This is a purely computational
> bioinformatics project: it statistically associates protein **domains** with
> Gene Ontology **function terms** using public databases (UniProt/GOA, InterPro,
> the Gene Ontology), reimplementing a published academic method (Fang & Gough
> 2013). No dual-use content — no sequence/organism design, no pathogen or toxin
> work, no wet-lab or synthesis protocols; just statistics over public data.

A working implementation of the domain-centric Gene Ontology (dcGO) methodology
for protein function prediction. The pipeline transforms protein-level GO
annotations into statistically validated **domain → GO term** associations using
Fisher's exact tests with FDR correction and hypergeometric association scoring.

> **Status:** the human analysis path (below) runs end-to-end and is the
> supported entry point. It uses **pre-computed InterPro domain annotations**
> (`protein2ipr.dat`) — it does **not** run InterProScan or scan sequences
> itself. A general multi-organism orchestrator is future work (see
> [FUTURE_WORK.md](FUTURE_WORK.md)); there is no `main_pipeline` module yet.

---

## What it does

1. **Parses GO annotations** for a species from GOA (`goa_human.gaf.gz`).
2. **Maps proteins to InterPro domains** from pre-computed `protein2ipr.dat`.
3. **Builds domain × GO contingency tables** and runs **Fisher's exact tests**
   (vectorized via the Cython `fisher` package), with **Benjamini–Hochberg FDR**
   correction and a **hypergeometric association score (1–100)**.
4. **Generates supra-domains** (contiguous domain combinations up to triplets)
   (contiguous combinations up to triplets).
5. *(optional)* **Applies the True Path Rule** to propagate associations up the
   GO DAG, when run with `--enable-true-path`; and *(optional)* applies
   the paper's **relative inference** — a parental-background filter — when run
   with `--enable-relative-inference`. These are separate steps of the published
   method and are selected separately; see [Two hierarchy stages](#two-hierarchy-stages).

Output: a TSV of significant domain–GO associations (FDR < 0.01 by default) with
p-values, FDR q-values, odds ratios, and hypergeometric scores.

---

## Quick Start (human, ~1 hour + download time)

See **[QUICKSTART.md](QUICKSTART.md)** for the full walkthrough.

```bash
# 1. Install dependencies
uv sync

# 2. Download the required datasets
#    (GOA human ~11 MB, GO ontology ~30 MB, InterPro protein2ipr ~20 GB)
uv run python scripts/download_data.py

# 3. Extract the species subset of protein2ipr (one-time; the 20 GB file is
#    streamed once and filtered down to the proteins for this species)
uv run python extract_human_interpro.py

# 4. Run the statistical inference
uv run python run_dcgo_human.py --num-cores 8

# 4b. (optional) also propagate annotations up the GO DAG (True Path Rule)
uv run python run_dcgo_human.py --num-cores 8 \
    --enable-true-path --go-ontology data/raw/go_ontology/go-basic.obo
```

`scripts/download_data.py --list` shows every dataset it knows about;
`--datasets NAME` selects a subset and `--all` grabs the optional sources too.

### Other species

The whole path is species-parameterised — download, extract and run all take a
`--species` flag (default `human`) and share the `goa_<species>` /
`protein2ipr_<species>` file naming. GOA publishes per-species annotations, so
e.g. mouse is:

```bash
uv run python scripts/download_data.py --species mouse   # GOA mouse + shared inputs
uv run python extract_human_interpro.py --species mouse  # mouse subset of protein2ipr
uv run python run_dcgo_human.py --species mouse --num-cores 8
```

The InterPro `protein2ipr.dat.gz` download is shared across species, so it is
only fetched once.

### Other ontologies (Enzyme Commission)

Domains can be associated with ontologies other than GO via `--ontology`. The
first non-GO ontology is **Enzyme Commission (EC)**, sourced from the Expasy
ENZYME database (`enzyme.dat`), which is already keyed by UniProt accession — so
no identifier mapping is needed:

```bash
uv run python scripts/download_data.py --datasets enzyme   # Expasy enzyme.dat
uv run python run_dcgo_human.py --ontology ec              # domain → EC associations
```

EC results are written to `results/domain_ec_associations_*.tsv` (with an
`ec_term` column), leaving the GO outputs untouched. `--enable-true-path` also
works for EC: associations are propagated up the EC hierarchy
(`1.1.1.1 → 1.1.1.- → 1.1.-.- → 1.-.-.-`) into
`results/domain_ec_annotations_propagated.tsv` — no ontology file needed, since
the hierarchy is implicit in the numbering.

### UniProt-native ontologies (Reactome, keywords, …)

Because the protein universe *is* UniProt, the cheapest term annotations are the
ones UniProt already carries per accession — no identifier mapping required. The
UniProt Swiss-Prot flat file (`uniprot_sprot.dat.gz`) cross-references many
resources (`DR` lines: Reactome, KEGG, GO, disease DBs, …) and carries a keyword
vocabulary (`KW` lines):

```bash
uv run python scripts/download_data.py --datasets uniprot_sprot_dat  # ~1 GB
uv run python run_dcgo_human.py --ontology reactome                 # Reactome pathway
uv run python run_dcgo_human.py --ontology keyword                  # UniProt keyword
uv run python run_dcgo_human.py --ontology disease                  # OMIM phenotype (DR MIM)
uv run python run_dcgo_human.py --ontology xref --xref-db KEGG      # any DR database
```

`disease` uses `DR MIM` restricted to `phenotype` entries (dropping the `gene`
links). `xref` opens **any** DR database by name, so KEGG / BRENDA / PANTHER / …
need no code change — add `--xref-type phenotype` to filter a typed database, or
`--xref-term-from-type` for databases that key the `DR` line by accession and put
the vocabulary in the third field. Results land in
`results/domain_<vocab>_associations_*.tsv` (e.g. `domain_disease_*`,
`domain_kegg_*`), leaving GO/EC outputs untouched.

#### The full ontology list

A survey of every `DR` database in the human subset of Swiss-Prot
(`docs/uniprot_ontology_survey.md`, reproducible with
`scripts/survey_uniprot_ontologies.py`) separates real vocabularies from 1:1
accession mirrors (AlphaFoldDB, STRING, GeneCards …, useless for an enrichment
test) and from domain databases (Pfam, PANTHER …, circular). The vocabularies,
plus the layers curated into the entry body rather than into `DR` lines, are
registered in `src/ontology_registry.py`:

| `--ontology` | Terms | Source | Hierarchy |
| --- | --- | --- | --- |
| `go` | Gene Ontology | GOA GAF | `go-basic.obo` |
| `ec` | Enzyme Commission | Expasy `enzyme.dat` | implicit in the number |
| `reactome` | pathways | `DR Reactome` | `ReactomePathwaysRelation.txt` |
| `keyword` | UniProt keywords | `KW` | `keywlist.txt` |
| `subcellular` | `SL-` locations | `CC SUBCELLULAR LOCATION` | `subcell.txt` |
| `ligand` | ChEBI ligands | `FT …/ligand_id` | `chebi_lite.obo` |
| `cofactor` | ChEBI cofactors | `CC COFACTOR` | `chebi_lite.obo` |
| `rhea` | Rhea reactions | `CC CATALYTIC ACTIVITY` | — |
| `tcdb` | transporter classes | `DR TCDB` | implicit (`8.A.98.1.10`) |
| `merops` | peptidase families | `DR MEROPS` | implicit (`S01.151 → S01 → S`) |
| `cazy` | CAZy families | `DR CAZy` | implicit (`GT32 → GT`) |
| `disease` | OMIM phenotypes | `DR MIM` (phenotype) | — |
| `doid` | Disease Ontology | `DR MIM` re-keyed via `doid.obo` | `doid.obo` |
| `orphanet` | rare diseases | `DR Orphanet` | — |
| `orphanet_doid` | Disease Ontology | `DR Orphanet` re-keyed via `doid.obo` | `doid.obo` |
| `hpo` | HPO phenotypes | `genes_to_phenotype.txt`, GeneID re-keyed via `DR GeneID` | `hp.obo` |
| `syngo` | SynGO synaptic terms | SynGO release zip, HGNC re-keyed via `DR HGNC` | same zip (`ontologies.xlsx`) |
| `unipathway` | metabolic pathways | `DR UniPathway` | — |
| `complex` | protein complexes | `DR ComplexPortal` | — |
| `drugbank` | drugs | `DR DrugBank` | — |
| `pharos` | target development level | `DR Pharos` (3rd field) | — |
| `condensate` | biomolecular condensates | `DR CD-CODE` (3rd field) | — |
| `xref` | anything else | `DR <--xref-db>` | — |

```bash
uv run python run_dcgo_human.py --ontology subcellular   # where the domain puts the protein
uv run python run_dcgo_human.py --ontology ligand        # what chemistry the domain binds
uv run python run_dcgo_human.py --ontology tcdb          # transporter classification
```

`--enable-true-path` propagates every ontology in the "Hierarchy" column above
through the shared engine in `src/hierarchy.py` — OBO graphs, hierarchies
implicit in the term id, and companion hierarchy files alike:

```bash
uv run python scripts/download_data.py --datasets reactome_relations uniprot_keywlist \
    uniprot_subcell chebi
uv run python run_dcgo_human.py --ontology reactome    --enable-true-path
uv run python run_dcgo_human.py --ontology subcellular --enable-true-path
uv run python run_dcgo_human.py --ontology ligand      --enable-true-path
```

For an ontology with no hierarchy (`disease`, `rhea`, `xref`, …),
`--enable-true-path` now **fails with an explicit error** rather than running
without propagation, and any missing input is reported before the expensive
stages start.

### Two hierarchy stages

The dcGO paper (Fang & Gough 2013) uses the hierarchy in two distinct places,
and this pipeline exposes them as two flags:

| Flag | Paper step | What it does | Direction |
|------|-----------|--------------|-----------|
| `--enable-relative-inference` | Step 2, "relative inference" | Also tests each association within the background of proteins annotated to its term's **direct parents**, and corrects the larger of the two p-values | **Removes** associations |
| `--enable-true-path` | Step 3, "true-path rule" | Propagates each association to its **ancestor** terms | **Adds** annotations |

Only Step 3 is the True Path Rule. The parental-background test is a separate
statistical inference that happens to consult the hierarchy — the paper is
explicit about this ("*two types of statistical inference followed by FDR
calculation*", then "*following the true-path rule to obtain the complete
domain-centric GO annotations*"), but it is easy to miss in a dense methods
section.

#### Why the *larger* p-value

The paper "*first took the larger one of the overall and relative p-values*",
then applied BH to that. Each inference alone fails in a characteristic way:
the overall test cannot locate which level of the DAG the signal lives at (a
domain associated with a broad term makes every descendant look enriched for
free), while the relative test rests on a background that can be small and
idiosyncratic, and roots have no parents at all. Requiring **both** is the
claim "real globally, *and* specific to this term rather than inherited".

Taking the maximum is also the statistically valid way to say that. It is an
**intersection-union test**: the null is "fails at least one inference",
rejection requires rejecting both, and `max(p₁, p₂)` is itself a valid p-value
for that null with no multiplicity correction (Berger 1982). `min(p₁, p₂)` is
not. So the maximum is exactly the quantity BH is entitled to correct — which
is why the combination happens *before* the correction. Consistently, the
h-score reported is the **smaller** of the overall and relative scores: the
weaker evidence governs in both directions.

#### What each stage does to the numbers

Human GO, `--disable-supra-domains`, FDR < 0.01:

| Configuration | Significant associations |
|---|---:|
| overall inference only | 44,453 |
| `--enable-relative-inference` | **3,876** |

The drop is large because the relative p-value is now held to the same FDR
standard as the overall one (a threshold of 1.2e-07 on this run), rather than
to a loose uncorrected cutoff. 623,092 of the 655,659 candidate pairs are
governed by the relative inference — that is, their relative p-value exceeds
their overall one — which is the inherited-association cascade the test exists
to remove.

Propagation, being additive, works the other way:

```
--enable-true-path                                44,453 -> 328,486 annotations
--enable-relative-inference --enable-true-path     3,876 ->  21,634 annotations
```

Two caveats worth knowing:

- **`--enable-relative-inference` works for every ontology with a hierarchy** —
  all 12 of them (`go`, `ec`, `reactome`, `keyword`, `doid`, `orphanet_doid`,
  `tcdb`, `merops`, `cazy`, `subcellular`, `ligand`, `cofactor`). The 9 flat
  cross-reference layers (`disease`, `orphanet`, `unipathway`, `complex`,
  `drugbank`, `pharos`, `condensate`, `rhea`, `xref`) have no parental
  background to test against and are rejected with an explicit error.
- **The relative inference is not usable yet, and both stages are off by
  default.** Terms with no parents skip the test and pass on the overall
  inference alone, so the GO roots dominate the output and enabling the layer
  makes results *less* specific — the opposite of its purpose. It needs an
  information-content floor. See `VALIDATION_PLAN.md` next-steps item 2.
- **Neither stage is on by default yet.** Both are opt-in while the numbers
  they change are re-measured.

### Disease: re-keying OMIM onto the Disease Ontology

OMIM is a catalogue, not an ontology. Its ids have no DAG, and it splits one
disease across many locus-specific entries, so `--ontology disease` is both
un-propagatable and badly underpowered. The Human Disease Ontology supplies the
missing structure: it is an `is_a` DAG *and* it cross-references OMIM
(`xref: MIM:<id>`) and Orphanet (`xref: ORDO:<id>`).

```bash
uv run python scripts/download_data.py --datasets disease_ontology
uv run python run_dcgo_human.py --ontology doid --enable-true-path
uv run python run_dcgo_human.py --ontology orphanet_doid --enable-true-path
```

The translation happens **at parse time**, in the protein→term map the Fisher
engine consumes, so sparse OMIM phenotypes pool into a better-supported DO class
*before* any test is run; a post-hoc re-labelling of the output could only
rename terms that had already reached significance. `disease` and `orphanet`
keep emitting the raw ids, so the two hypothesis universes stay comparable.

The mapping is not one-to-one, and `src/disease_ontology.py` documents (and
counts, at parse time) what happens to each case: **unmapped** ids are dropped
but logged with their annotation counts, **one-to-many** ids expand to every DO
term, and **obsolete** DO terms are skipped unless `replaced_by` resolves. It
also reports mapping coverage — over distinct ids *and* over protein
annotations — because a layer that reaches significance by shrinking its own
term space has not learned anything.

The `disease_ontology` dataset is pinned to an immutable OBO Foundry release
PURL with a SHA-256 checksum in `config/settings.py`, which
`scripts/download_data.py` verifies on every run.

**What it actually bought** (human, current UniProt, FDR<0.01): 74% annotation
coverage, a hierarchy where there was none — 160 True Path annotations, 16
direct + 144 propagated — and interpretable term labels. It did *not* buy more
significant associations: 16 for `doid` against 17 for `disease`, on a term
space shrunk from 6,904 to 4,917 and a protein universe from 5,029 to 3,928.
Both layers return **0** associations under `--permute-annotations 7`, so
neither count is an FDR artefact. The remaining sparsity is at the protein level
(≈0.8 proteins per DO term), which pooling through a hierarchy cannot fix. See
`TODO.md` for the full before/after table.

> **Note on evidence:** these are UniProt-native *cross-references*, not GO
> annotations, so there is no IEA/evidence code to filter. (GO is the exception —
> UniProt's `DR GO` lines include IEA, which is why GO is pulled from GOA with its
> evidence filter instead.) No source is "preferred"; if the same annotation
> appears in UniProt and a primary DB, they deduplicate to the union.

Adding any ontology means one `OntologyEntry` in `src/ontology_registry.py`,
backed by an `AnnotationSource` subclass — see `src/annotation_source.py`,
`src/ec_annotation_source.py`, and `src/uniprot_annotation_source.py` for the
pattern.

---

## Finding the emergent predictions (surprise score)

The associations that matter most are the ones a *combination* of domains
supports but none of its constituents does — the signal single-domain and
homology methods cannot see. `scripts/rank_surprising_associations.py` ranks
them:

```bash
uv run python scripts/rank_surprising_associations.py --ontology go
```

Each supra-domain association is scored as
`-log10(q_emergence) × distinctness × novelty`: a binomial test of the observed
rate against what the parts already predict (noisy-OR over constituents, floored
by the best sub-combination and by the term's background rate), times a penalty
for constituents that are really one region annotated by redundant InterPro
signatures, times a discount for what InterPro2GO already records.

On the current human GO run this puts textbook multi-domain architectures on top
— SH2 + kinase → non-receptor tyrosine kinase activity, PH + EF-hand → PLC
activity, BTB/POZ + C2H2 → transcriptional repressor — recovered without being
told about them. Output: `results/domain_<ontology>_surprising.tsv`, with every
component in its own column. Full method, results and caveats:
**[SURPRISE_SCORE.md](SURPRISE_SCORE.md)**.

---

## Required inputs

The dcGO methodology needs three inputs for a given set of proteins. All three
are downloaded by `scripts/download_data.py` into `data/raw/<source>/`:

| Input | File | Size | Purpose |
|-------|------|------|---------|
| Domain annotations | `interpro_mappings/protein2ipr.dat.gz` | ~20 GB | Which InterPro domains are in each protein |
| GO annotations | `goa_annotations/goa_<species>.gaf.gz` | ~11 MB (human) | Protein → GO term assignments (GAF 2.2) |
| Ontology structure | `go_ontology/go-basic.obo` | ~30 MB | GO DAG (only needed for `--enable-true-path`) |

`extract_human_interpro.py` filters `protein2ipr.dat.gz` down to the proteins
of the chosen species (`--species`, default `human`) found in the GOA file,
writing `data/interim/protein2ipr_<species>.dat.gz` so subsequent runs are fast.

---

## Key options (`dcgo`)

| Flag | Default | Description |
|------|---------|-------------|
| `--species` | `human` | Species / GOA file to analyze |
| `--ontology` | `go` | Ontology to associate domains with — 19 registered, see the table above or `--help` |
| `--xref-db` | — | UniProt DR database name (required for `--ontology xref`, e.g. `KEGG`) |
| `--xref-type` | — | Optional DR third-field filter for `xref` (e.g. `phenotype`) |
| `--xref-term-from-type` | off | For `xref`, take the term from the DR third field instead of the id |
| `--enzyme-dat` | `data/raw/enzyme/enzyme.dat` | Expasy ENZYME file (used when `--ontology ec`) |
| `--uniprot-dat` | `data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz` | UniProt flat file (every UniProt-native ontology) |
| `--subcell` | `data/raw/uniprot_subcell/subcell.txt` | Subcellular-location vocabulary (`--ontology subcellular`) |
| `--chebi-obo` | `data/raw/chebi/chebi_lite.obo` | ChEBI ontology (True Path for `ligand`/`cofactor`) |
| `--evidence-filter` | `manual` | GO evidence codes: `all`, `manual`, `experimental` |
| `--fdr-threshold` | `0.01` | FDR (q-value) significance cutoff |
| `--num-cores` | `8` | CPU cores for parallel Fisher tests |
| `--batch-size` | `50000` | Fisher test batch size |
| `--enable-supra-domains` / `--disable-supra-domains` | enabled | Test contiguous domain combinations |
| `--enable-true-path` | off | Propagate associations up the term hierarchy, and only that (fails if the ontology has none) |
| `--enable-relative-inference` | off | Parental-background filter: keep an association only if still enriched within its term's direct parents. Any ontology with a hierarchy; it *removes* associations |
| `--go-ontology` | `data/raw/go_ontology/go-basic.obo` | GO OBO file (`--ontology go` only; read by either of the two flags above) |
| `--output-dir` | `results/` | Output directory |

---

## Example output

The ranked TSV has the full header
`rank, domain, go_term, p_value, adj_p_value, odds_ratio, hyper_score, domain_type, constituent_domains, n_observations`
(the unranked export is identical without the leading `rank` column):

```
rank  domain      go_term      p_value        adj_p_value    odds_ratio  hyper_score  domain_type  constituent_domains  n_observations
1     IPR015812   GO:0005178   2.894064e-307  8.791816e-299  inf         100.00       single       -                    42
2     IPR000471   GO:0005125   4.869485e-294  7.396453e-286  inf         100.00       single       -                    35
```

---

## Validation

Two independent checks live in `validation/` (see
[VALIDATION_PLAN.md](VALIDATION_PLAN.md) for the full plan):

- **InterPro2GO coverage (§1)** — treats the curated InterPro2GO map as an
  incomplete *positive* reference and reports recall on shared domains
  (propagated). dcGO recovers **~65%** of curated pairs at FDR<0.01.
  `validation/validate_results.py`.
- **Temporal CAFA-style benchmark (§2)** — trains on the 2021 GOA and scores
  protein-centric predictions against **newly-curated, experimentally-supported**
  2026 annotations (CAFA *no-knowledge* design, per aspect):

  ```bash
  # fetch a dated GOA snapshot (t0), run the pipeline on it, then score
  uv run python scripts/download_data.py --goa-archive 205        # 2021-04
  uv run python validation/temporal_benchmark.py \
      --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \
      --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
      --predictions results_t0_2021/domain_go_associations_significant.tsv \
      --min-ic 0 --min-ic 2 --min-ic 4      # information-content sweep
  ```

  F_max, dcGO (p-score) vs the two baselines, as an **information-content floor**
  removes near-universal low-IC terms (e.g. `protein binding`, 85% of exp. MF):

  | Aspect | IC≥0 dcGO/naive | IC≥2 | IC≥4 | dcGO/random @ IC≥4 |
  |--------|:---------------:|:----:|:----:|:------------------:|
  | BP | **0.248** / 0.115 | **0.170** / 0.071 | **0.115** / 0.031 | 6.1× |
  | MF | 0.360 / **0.464** | **0.365** / 0.053 | **0.337** / 0.045 | 4.7× |
  | CC | **0.380** / 0.343 | **0.239** / 0.153 | **0.134** / 0.099 | 4.3× |

  dcGO beats naive at face value on BP and CC; MF (whose truth is dominated by
  `protein binding`) flips to dcGO the moment uninformative terms are excluded.
  On informative terms dcGO beats both baselines in every aspect, staying
  1.3–25× above the random-domain null. See VALIDATION_PLAN.md §2.

---

## Development

```bash
# Tests (525 tests, ~8 s)
uv run pytest

# Lint + format (the same set CI checks)
uv run ruff check src/ tests/ config/ scripts/ validation/ \
    run_dcgo_human.py extract_human_interpro.py
uv run ruff format --check

# Coverage
uv run pytest --cov=src --cov-report=html
```

CI (`.github/workflows/ci.yml`) runs ruff + pytest on every push and PR, then
builds the wheel and sdist and smoke-tests the installed `dcgo` CLI.

### Reproducible runs

Every analysis run writes `run_manifest_<ontology>.json` in its output
directory. The manifest records the SHA-256 (and embedded release header, where
the format has one) of every input the chosen ontology consumed, the Git
revision and dirty state, the `uv.lock` hash, the full command line, every
effective parameter and threshold, runtime metadata, timestamps, summary counts
and output hashes. A completed run has `"status": "completed"`; an interrupted
one leaves `"status": "running"` for diagnosis.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the publication/release
checklist. Citation metadata is in [CITATION.cff](CITATION.cff).

---

## Repository layout

```
dcGO-2.0/
├── run_dcgo_human.py            # Pipeline implementation + legacy script entry
├── extract_human_interpro.py    # Filter protein2ipr.dat down to human proteins
├── scripts/
│   └── download_data.py         # Download required datasets from source
├── src/
│   ├── runner.py                # Generic request + input-resolution stage
│   ├── goa_parser.py            # Parse GOA GAF files
│   ├── domain_annotation_parser.py  # Parse protein2ipr domain mappings
│   ├── vectorized_fisher.py     # Vectorized Fisher's exact test + BH-FDR
│   ├── sparse_fisher.py         # Sparse contingency-table construction
│   └── ontology_processor.py    # True Path Rule / GO DAG propagation
├── config/settings.py           # Dataset URLs + configuration
├── tests/                       # unit / integration tests
├── validation/                  # validate_results.py (§1) + temporal_benchmark.py (§2)
└── docs/                        # Reference papers
```

> **Downloading data:** use `scripts/download_data.py`, which reads its URLs
> from `config/settings.py`. An older async (aiohttp) downloader and a SQLite
> storage layer used to sit in `src/`; both were unreachable from any supported
> entry point and were removed.

---

## Method references

- Fang & Gough (2013), *dcGO: database of domain-centric ontologies*. The
  statistical framework (per-domain enrichment, FDR, True Path Rule) follows
  this work.
- Gene Ontology Consortium — ontology structure and annotation principles.
- InterPro — domain classification (`protein2ipr`).

---

## Known limitations / not yet done

See [FUTURE_WORK.md](FUTURE_WORK.md) and [TODO.md](TODO.md). In brief:

- Only pre-computed InterPro annotations are consumed; no local domain scanning.
- The True Path Rule is **opt-in**, not part of the default run.
- Validation covers InterPro2GO coverage (§1) and a temporal CAFA-style benchmark
  (§2, above); still open are the ablation study (§4), per-protein score
  calibration (§5), and a comparison to the original dcGO results (§3). See
  [VALIDATION_PLAN.md](VALIDATION_PLAN.md).

---

## License

MIT License — see [LICENSE](LICENSE).
