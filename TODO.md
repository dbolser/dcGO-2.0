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

## Done — Disease Ontology re-keying

- ~~**Give the disease layers a hierarchy via the Disease Ontology.**~~ —
  `src/disease_ontology.py`, `--ontology doid` and `--ontology orphanet_doid`.
  UniProt's `DR MIM` (phenotype) and `DR Orphanet` ids are re-keyed onto DOID
  terms **at parse time**, using DO's own `xref: MIM:` / `xref: ORDO:` lines,
  and then propagate up DO's `is_a` DAG. `disease`/`orphanet` keep emitting the
  raw ids so the two hypothesis universes stay comparable. The DO release is
  pinned to an immutable OBO Foundry release PURL with a SHA-256 checksum that
  `scripts/download_data.py` now verifies.

  **The result is not the one the acceptance criterion asked for, and that is
  the finding.** Current UniProt, human, FDR<0.01. These are *not* comparable to
  the 53 significant associations the 2021 t0 run produced — different snapshot,
  different term space:

  | | `disease` (OMIM ids) | `doid` (re-keyed) |
  | --- | ---: | ---: |
  | proteins | 5,029 | 3,928 |
  | domain features | 51,311 | 43,019 |
  | terms | 6,904 | 4,917 |
  | Fisher tests | 354,251,144 | 211,524,423 |
  | BH threshold p | 4.72e-10 | 3.96e-10 |
  | significant | 17 | **16** |
  | distinct terms among them | 4 | 2 |
  | permutation control (seed 7) | 0 | 0 |
  | True Path annotations | *impossible* | 160 (16 direct + 144 propagated) |

  So re-keying did **not** buy more significant associations — it bought a
  hierarchy (144 propagated annotations the OMIM layer cannot produce at all)
  and interpretable term labels, at the cost of ~26% of the annotations. Note
  the direction of the arithmetic: a *smaller* term space with a *similar* count
  is not a win, and the honest read is that this layer's sparsity is at the
  protein level (3,928 proteins over 4,917 DO terms ≈ 0.8 proteins/term), which
  pooling through DO does not fix.

  Mapping coverage (the first-class number): 5,087 / 6,920 distinct OMIM ids
  (73.5%) and 5,534 / 7,457 protein–term annotations (74.2%). 1,833 OMIM ids had
  no DO term. One concrete loss: `IPR051503` (complement system regulators) →
  MIM 235400 (atypical HUS) disappears, because DO cross-references
  `DOID:0080301` only as `ORDO:2134`, with no `MIM:` xref — which is exactly why
  `orphanet_doid` exists as a second, complementary route.

- **Still open here:** re-run `validation/temporal_breadth.py` on `doid` at the
  2021 t0 to replace the "undefined (0 hits / 369)" `disease` row; and try
  MONDO (`https://purl.obolibrary.org/obo/mondo.obo`, what dcGO 2023 switched
  to), which has broader OMIM coverage and would use the same machinery — only
  the OBO and the xref prefix change.

## Loose ideas / nice-to-haves
- Add InterPro names / gene names / GO term descriptions to the output TSVs
  (currently just IDs) — improves readability for downstream users.
- Protein → genome positions from Ensembl (for a genome-browser view).
- Build anomaly/analysis data dir (`BUILD NOMALY DATADIR` — original note).
- Other ontologies (DO/HPO/MP/EC/Reactome) — see FUTURE_WORK.md.
