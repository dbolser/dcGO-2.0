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

## Done — 2026-07-28 session (PR #26)
- ~~**Surprise score**~~ — `src/surprise_score.py`,
  `scripts/rank_surprising_associations.py`, `SURPRISE_SCORE.md`. Ranks emergent
  supra-domain associations by `-log10(q_emergence) × distinctness × novelty`:
  a binomial test against a parts-only expectation (noisy-OR over constituents,
  floored by the best sub-combination and by the term's background rate), times a
  redundant-InterPro-signature penalty from matched-region overlap, times an
  InterPro2GO novelty discount. Top of the GO ranking is textbook architectures
  (SH2+kinase, PH+EF-hand, BTB/POZ+C2H2, kinase+SAM).
- ~~**13 new ontologies behind a registry**~~ — `src/ontology_registry.py` is now
  the single `--ontology` dispatch table (19 keys). New: orphanet, tcdb, merops,
  cazy, unipathway, complex, drugbank, pharos, condensate, plus the layers
  curated in the entry body rather than `DR` lines — subcellular (CC → `SL-`),
  ligand (FT `/ligand_id` ChEBI), cofactor, rhea. Survey of all ~150 DR databases
  in `docs/uniprot_ontology_survey.md` + `docs/dr_survey.tsv`.
- ~~**Held-out validation of the surprise score**~~ — `validation/temporal_surprise.py`.
  The *associations* predict future curation (12.5×); the *ranking* is not
  demonstrably better than the dcGO q-value, and at matched prediction budgets
  the comparison is not resolvable at all — 9,923 of 10,136 associations score
  exactly 0.000, so most of a budget slice is an arbitrary tie-break. Verdict,
  the bootstrap diagnosis behind it, and the percentile/basic/BCa intervals
  side by side, in `SURPRISE_SCORE.md`.
- ~~**Predictive power across the breadth**~~ — `validation/temporal_breadth.py`.
  One archived Swiss-Prot release (2021_02) gives t0 for every UniProt-native
  layer. GO 11.3×, reactome 8.0×, cofactor 3.2×, subcellular 2.9×, keyword 1.7×;
  complex degenerate, disease undefined, ligand untestable. See
  `VALIDATION_PLAN.md` §2 breadth subsection.

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

## Queued

- **Surprise score v2: weigh emergence against testability.** The held-out test
  exposed a structural tension — emergence requires that a combination's carriers
  are *already* nearly all annotated, so the most emergent associations leave the
  fewest standing predictions and are the least verifiable (`surprise top-25`
  yields 117 predictions against 0.14 expected hits). A useful score should trade
  emergence off against how much it still predicts.
- **Re-test `ligand` and `cofactor` at a 2023 t0.** UniProt only introduced the
  structured `FT /ligand_id="ChEBI:…"` qualifier after 2021 (April 2021 used
  free-text `/note="ATP"`), so the ligand layer is entirely post-2022 annotation
  and untestable on the 2021→2026 split. Needs `release-2023_01` or similar —
  a shorter, non-comparable window, so report separately.
- **Breadth test with `--supra-only`.** The per-ontology numbers pool single
  domains and supra-domains; isolating supra-domains would say whether the
  *emergent* claim generalises beyond GO, which is the more interesting question.
- **Decide the fate of `agent/reproducible-runs`.** Two unpushed commits, never
  PR'd, worktree in `/tmp/dcgo-repro` (will not survive a reboot). Contains
  `src/run_manifest.py` + tests and `CITATION.cff`, which close P0/P1 review
  items. It conflicts with PR #26 on `run_dcgo_human.py`; cheapest path is to
  cherry-pick those two files onto a fresh branch after #26 merges.
- **`validation/bench_A`–`bench_D` are still untracked** — the review flagged
  their provenance as unclear. Either commit or archive them.

- **Give the disease layers a hierarchy via the Disease Ontology.** `--ontology
  disease` is raw OMIM phenotype ids with no DAG, so `--enable-true-path`
  refuses, and the sparsity shows: the 2021 t0 run yielded only **53**
  significant associations from 6,904 terms over 5,029 proteins. DO supplies the
  missing structure.
  - Download: `http://purl.obolibrary.org/obo/doid.obo` (checked 2026-07-28:
    live, 7.2 MB, 14,735 terms, release 2026-06-30). MONDO —
    `http://purl.obolibrary.org/obo/mondo.obo` — is the alternative, and is what
    dcGO 2023 switched to.
  - The join already exists: DO cross-references OMIM as **`xref: MIM:<id>`**
    (note the prefix is `MIM`, *not* `OMIM`) — 6,467 xrefs over 6,123 distinct
    MIM ids — and Orphanet as `xref: ORDO:<id>` (2,319), which would give
    `--ontology orphanet` a hierarchy too.
  - Two ways to use it, and the choice matters: post-hoc mapping of existing
    output is nearly worthless here (the 53 t0 associations use just 9 distinct
    OMIM terms, 6 mappable), so the value is in **re-keying at parse time** — a
    disease source that emits DOIDs, letting sparse OMIM phenotypes pool into
    better-supported DO classes *before* the Fisher tests, then propagating up
    the DO DAG.
  - Machinery needed: none new. `parse_obo_child_parents` already reads DO's
    `is_a` edges, and the registry takes the ancestors factory. The work is an
    OMIM→DOID translation in the annotation source plus a registry entry.
  - Acceptance: more significant associations than the OMIM-keyed run, and a
    `disease` row in the breadth test that is no longer underpowered.

## Loose ideas / nice-to-haves
- Add InterPro names / gene names / GO term descriptions to the output TSVs
  (currently just IDs) — improves readability for downstream users.
- Protein → genome positions from Ensembl (for a genome-browser view).
- Build anomaly/analysis data dir (`BUILD NOMALY DATADIR` — original note).
- Other ontologies (DO/HPO/MP/EC/Reactome) — see FUTURE_WORK.md.
