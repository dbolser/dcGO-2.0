# TODO — scratchpad

> **Scope & safety — safe to review.** dcGO-2.0 is a purely computational
> bioinformatics project: it statistically associates protein **domains** with
> Gene Ontology **function terms** using public databases (UniProt/GOA, InterPro,
> the Gene Ontology), reimplementing a published academic method (Fang & Gough
> 2013). No dual-use content — no sequence/organism design, no pathogen or toxin
> work, no wet-lab or synthesis protocols; just statistics (Fisher's exact tests,
> FDR) over existing public data.

The structured roadmap lives in **[VALIDATION_PLAN.md](VALIDATION_PLAN.md)**
(validation/benchmarking, §1–§6, tracked as issues #8–#12) and
**[FUTURE_WORK.md](FUTURE_WORK.md)** (expanding beyond GO to other ontologies).
This file is just loose notes.

## Done
- ~~Switch to pre-computed InterPro download~~ — `scripts/download_data.py`.
- ~~Download the InterPro→GO mapping~~ — `interpro2go` data source.
- ~~Validate results against InterPro2GO~~ — §1 (coverage reframe), ~65% at FDR<0.01.
- ~~Fix pipeline correctness~~ — contingency-table + True Path Rule bugs (#15, #17).
- ~~§2 temporal/CAFA benchmark~~ (#8) — 2021→2026 no-knowledge split.
  `validation/temporal_benchmark.py`. On **informative** terms (IC floor) dcGO
  beats both baselines in every aspect: 4–26× the random-domain null and above
  naive. Naive's raw-F_max lead was base-rate recovery of near-universal terms
  (`protein binding` = 85% of experimental MF). Dated GOA via
  `download_data.py --goa-archive`; IC sweep via `--min-ic`.

## Next (see VALIDATION_PLAN.md) — NOT ablation yet
- **Method-vs-paper audit** — read Fang & Gough 2013 (`docs/`) for how they
  *validated* (likely domain-centric, not protein-centric F_max) and the
  optimal-level / True-Path test (OFF in the §2 run). Align where it matters.
- **Harden the metric** — promote the IC-controlled view; consider a
  domain-centric evaluation (does a domain's predicted GO match known function).
- §4 ablation (#10), §3 original-dcGO comparison (#9), §5–§6 — after the
  yardstick is trusted.

## Loose ideas / nice-to-haves
- Add InterPro names / gene names / GO term descriptions to the output TSVs
  (currently just IDs) — improves readability for downstream users.
- Protein → genome positions from Ensembl (for a genome-browser view).
- Build anomaly/analysis data dir (`BUILD NOMALY DATADIR` — original note).
- Other ontologies (DO/HPO/MP/EC/Reactome) — see FUTURE_WORK.md.
