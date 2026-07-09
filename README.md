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
   with optional empirical-Bayes shrinkage toward their constituent domains.
5. *(optional)* **Applies the True Path Rule** to propagate associations up the
   GO DAG, when run with `--enable-true-path`.

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

# 3. Extract the human subset of protein2ipr (one-time; the 20 GB file is
#    streamed once and filtered down to human proteins)
uv run python extract_human_interpro.py

# 4. Run the statistical inference
uv run python run_dcgo_human.py --num-cores 8

# 4b. (optional) also propagate annotations up the GO DAG (True Path Rule)
uv run python run_dcgo_human.py --num-cores 8 \
    --enable-true-path --go-ontology data/raw/go_ontology/go-basic.obo
```

`scripts/download_data.py --list` shows every dataset it knows about;
`--datasets NAME` selects a subset and `--all` grabs the optional sources too.

---

## Required inputs

The dcGO methodology needs three inputs for a given set of proteins. All three
are downloaded by `scripts/download_data.py` into `data/raw/<source>/`:

| Input | File | Size | Purpose |
|-------|------|------|---------|
| Domain annotations | `interpro_mappings/protein2ipr.dat.gz` | ~20 GB | Which InterPro domains are in each protein |
| GO annotations | `goa_annotations/goa_human.gaf.gz` | ~11 MB | Protein → GO term assignments (GAF 2.2) |
| Ontology structure | `go_ontology/go-basic.obo` | ~30 MB | GO DAG (only needed for `--enable-true-path`) |

`extract_human_interpro.py` filters `protein2ipr.dat.gz` down to the human
proteins found in the GOA file, writing
`data/interim/protein2ipr_human.dat.gz` so subsequent runs are fast.

---

## Key options (`run_dcgo_human.py`)

| Flag | Default | Description |
|------|---------|-------------|
| `--species` | `human` | Species / GOA file to analyze |
| `--evidence-filter` | `manual` | GO evidence codes: `all`, `manual`, `experimental` |
| `--fdr-threshold` | `0.01` | FDR (q-value) significance cutoff |
| `--num-cores` | `8` | CPU cores for parallel Fisher tests |
| `--batch-size` | `50000` | Fisher test batch size |
| `--enable-supra-domains` / `--disable-supra-domains` | enabled | Test contiguous domain combinations |
| `--enable-shrinkage` | off | Empirical-Bayes shrinkage for supra-domains |
| `--enable-true-path` | off | Propagate associations up the GO DAG |
| `--go-ontology` | `data/raw/go_ontology/go-basic.obo` | OBO file (required for True Path Rule) |
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
      --predictions results_t0_2021/domain_go_associations_significant.tsv
  ```

  | Aspect | dcGO F_max | naive | random-domain |
  |--------|-----------:|------:|--------------:|
  | BP | 0.276 | 0.289 | 0.144 |
  | MF | 0.319 | 0.579 | 0.098 |
  | CC | 0.395 | 0.520 | 0.234 |

  dcGO clears the random-domain null by **1.7–3.2×** (the associations are
  genuinely informative) but does not beat the CAFA naive frequency baseline on
  F_max — competitive on BP, below on MF/CC. Honest, precision-capable result;
  see VALIDATION_PLAN.md §2 for the interpretation and next steps (§4/§5).

---

## Development

```bash
# Tests (155 tests, ~5 s)
uv run pytest

# Lint + format (CI uses ruff)
uv run ruff check src/ tests/
uv run ruff format --check

# Coverage
uv run pytest --cov=src --cov-report=html
```

CI (`.github/workflows/ci.yml`) runs ruff + pytest on every push and PR.

---

## Repository layout

```
dcGO-2.0/
├── run_dcgo_human.py            # Main entry point: human dcGO analysis
├── extract_human_interpro.py    # Filter protein2ipr.dat down to human proteins
├── scripts/
│   └── download_data.py         # Download required datasets from source
├── src/
│   ├── goa_parser.py            # Parse GOA GAF files
│   ├── domain_annotation_parser.py  # Parse protein2ipr domain mappings
│   ├── vectorized_fisher.py     # Vectorized Fisher's exact test + BH-FDR
│   ├── sparse_fisher.py         # Sparse contingency-table construction
│   ├── hierarchical_inference.py    # Supra-domains + shrinkage
│   ├── ontology_processor.py    # True Path Rule / GO DAG propagation
│   ├── data_acquisition.py      # (async downloader library — see note below)
│   └── database_manager.py      # SQLite storage/export helpers
├── config/settings.py           # Dataset URLs + configuration
├── tests/                       # unit / integration tests
├── validation/                  # validate_results.py (§1) + temporal_benchmark.py (§2)
└── docs/                        # Reference papers
```

> **Note:** `src/data_acquisition.py` is an older async (aiohttp) download
> library and is **not** on the supported path. Use `scripts/download_data.py`
> instead — it reads the same URLs from `config/settings.py`.

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
