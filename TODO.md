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

## Next (see VALIDATION_PLAN.md "Next steps")
- ~~Method-vs-paper audit~~ done — their validation is protein-centric CAFA
  PR-RC; restored the two missing pieces (relative inference + p-score) and added
  a domain-centric eval.
- **Temporal domain-centric test** — fetch a dated 2021 `interpro2go`, pass as
  `--reference` to `domain_centric_eval.py`.
- **Fold method into the pipeline** — wire the relative (parental-background)
  test into `run_dcgo_human.py` (combine + FDR<1e-3, per paper) and expose the
  p-score predictor as the standard path.
- §4 ablation (#10), §3 original-dcGO domain re-keying SSF/PF (#9), §5–§6.

## Loose ideas / nice-to-haves
- Add InterPro names / gene names / GO term descriptions to the output TSVs
  (currently just IDs) — improves readability for downstream users.
- Protein → genome positions from Ensembl (for a genome-browser view).
- Build anomaly/analysis data dir (`BUILD NOMALY DATADIR` — original note).
- Other ontologies (DO/HPO/MP/EC/Reactome) — see FUTURE_WORK.md.
