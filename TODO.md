# TODO

Last reviewed against `main` (`d166013`) on 2026-08-31. This is the short,
prioritized project queue. Detailed experiment designs and results live in
[VALIDATION_PLAN.md](VALIDATION_PLAN.md), [MULTISPECIES_BACKGROUND.md](MULTISPECIES_BACKGROUND.md),
and [REPRODUCIBILITY.md](REPRODUCIBILITY.md); ontology expansion in
[FUTURE_WORK.md](FUTURE_WORK.md); the current snapshot in [STATUS.md](STATUS.md).

Since the previous review (`5d3ef91`, 2026-08-11) the paper-parity machinery
landed end to end (#61–#70), the registry grew to 35 ontology keys
(#66/#68/#71), and the full 63-cell production matrix ran clean on 2026-08-18
(`results/production/`, untracked). The remaining blockers are chiefly
publication design, statistical independence, and reproducibility rather than
implementation.

## P0 — publication blockers

- [ ] **Freeze and run a genuinely untouched final evaluation.** Model choices
  (p-score transfer, IC floors, minimum support, True Path policy, and feature
  family) were inspected on the headline 2021→2026 split. Use the nested human
  interval for development, pre-specify the primary endpoint and cohort, then
  evaluate once on a later untouched interval or other untouched dataset.
- [ ] **Remove temporal look-ahead from domain architectures.** The temporal
  benchmarks train on dated GOA but reuse current `protein2ipr`. Obtain an
  archived t0 InterPro/protein2ipr snapshot, or explicitly limit every claim to
  an annotation-temporal benchmark. This also closes the 31.9% universe loss in
  the all-species held-out arm.
- [ ] **Correct phylogenetic non-independence (gene / UniRef50 collapse).**
  Protein-level pooling inflates all-species support by about 2.44× and half of
  associations rest on at most three UniRef50 clusters; the model-organism and
  all-species *surprise rankings* are additionally pseudo-replicated by
  one-gene-to-many-isoform mapping (top fly/worm rows are single-gene stacks —
  human layers are verified clean). Test at an ortholog/cluster level (or
  justify a weighting scheme), report pooled and corrected counts side by side,
  and repeat held-out and permutation analyses after correction. No claims from
  the affected layers until then.
- [ ] **Re-run the stale evaluations post-#67 (regulates-edge fix).** The
  manuscript's §3.2–3.7 blocks are marked PROVISIONAL because they predate the
  fix, and the §4 ablation's "True Path worse in 12/12 cells" predates both the
  background fix (#46) and the relative-inference/True-Path flag split —
  re-measure before citing either.
- [ ] **Choose and justify the primary method configuration.** The (stale)
  ablation finds no protein-centric gain from supra-domains, although they
  remain central to the emergent-association claim. Pre-specify whether the
  primary predictor is single-domain/no-TPR and treat supra-domain discovery as
  a separate analysis if so. Set and justify a minimum-support/effect-size
  policy rather than merely exposing `--min-support`. Includes deciding the
  HPO paper-parity collapse (996 → 38 associations, driven by relative
  inference, not the IC floor): specificity or over-conservatism is decidable
  only by the post-#67 temporal evaluation.
- [ ] **Add a stronger external comparator.** Original dcGO coverage is now
  measured, but the predictive benchmark still needs at least one independent
  domain- or protein-function predictor. Treat InterPro2GO strictly as an
  incomplete positive reference, not a precision benchmark.
- [ ] **Independently verify the evaluator.** Check F_max and AUPRC against a
  standard CAFA implementation; document threshold sampling, interpolation,
  coverage handling, and the changing cohorts induced by IC filtering.

## P1 — reproducibility and release

- [ ] **Pin every reported input to an immutable release.** Current manifests
  record hashes but several production URLs still point at mutable releases.
  Pin the exact t1 GOA and GO ontology as well as GOA, InterPro, UniProt, and all
  ontology inputs; verify expected hashes during download.
- [ ] **Add manifests to downstream analyses.** Cover the temporal benchmarks,
  ablation, breadth, comparison, and surprise-score tools, including random
  seeds and the hashes of upstream association files.
- [ ] **Provide one-command clean-checkout reproduction.** Download and verify
  inputs, run inference, regenerate manuscript tables/figures, and document
  runtime, RAM, disk, and supported platforms. Archive exact snapshots and
  outputs with a DOI.
- [ ] **Decide the fate of the untracked production outputs.**
  `results/production/` (2.7 GB, the 63-cell matrix + surprise rankings) and
  `data/ACQUISITION_MATRIX.md` (acquisition ledger, marked "do not commit")
  need an archival home — release artifact, DOI deposit, or documented
  regeneration path.
- [ ] **Make the calibrated p-score protein-prediction path reproducible**
  without post-hoc scripts (the inference side of paper parity is done; this
  transfer step is the remainder).
- [ ] **Create one generated source of truth for reported metrics.** README,
  `VALIDATION_PLAN.md`, the manuscript, and committed TSVs currently drift.
  Generate prose tables from committed metrics and retain provenance for
  `bench_A`–`bench_D` and their logs.
- [ ] **Finish release essentials.** Add a changelog, semantic-versioning and
  support policy, archived release/DOI, and resource estimates. `CITATION.cff`,
  wheel/sdist builds, installed-CLI smoke tests, and CI are already present.

## P2 — focused follow-up analyses

- [ ] Run the temporal domain-centric evaluation against a dated t0
  `interpro2go`, not the current mapping.
- [ ] Extend the published-dcGO comparison to `--domain-key pfam` and the
  non-GO Domain2EC/Domain2KW/Domain2DO tables.
- [ ] Re-test ligand and cofactor at a 2023 t0, when structured ChEBI ligand
  identifiers are available; report the shorter window separately.
- [ ] Repeat the breadth benchmark with `--supra-only` to test whether emergent
  combinations generalize beyond GO.
- [ ] Redesign the surprise score to balance emergence with testability. The
  current score predicts future curation in aggregate (12.5× enrichment) but
  has no demonstrated ranking advantage over q-value and leaves 9,923/10,136
  associations tied at zero; a graded emergence score is the precondition for
  redoing the head-to-head.
- [ ] Evaluate Mondo as the unifying disease layer. The adapter landed
  (`--ontology mondo` / `orphanet_mondo`, #71); what is open is the analysis —
  whether Mondo re-keying beats DOID for sparse disease annotations (the 2021
  DOID breadth run produced 0/160 hits).
- [ ] Quantify sensitivity to evidence policy, FDR, minimum support, effect-size
  floor, supra-domain length, GO release, and transfer rule; include qualitative
  error analysis of representative successes and failures.
- [ ] Investigate dependence-aware or hierarchical multiple testing across GO
  terms, domains, and nested supra-domain hypotheses — including `elim`-style
  decorrelation (Alexa 2006): even at IC≥1, ~50% of surviving single-domain
  associations still sit on an ancestor chain.

## Blocked on external inputs

- [ ] MAxO: the ontology is on disk but there is no annotation source — needs
  Monarch's `maxo-annotations`.
- [ ] SNOMED CT and MedDRA are licence-gated; OMIM `genemap2.txt` is
  registration-gated. All three need Dan's credentials.

## Documentation and maintenance

- [ ] Merge or close **PR #72** (`agent/manuscript-ontology-update`, 10
  commits) — held for Dan's read since 2026-08-18 — then prune the five stale
  `agent/*` branches and their two worktrees.
- [ ] Reconcile stale status text. README's counts predate the ontology sprint
  (the suite is now 963 tests); `RESULTS.md` contains two claims that
  `VALIDATION_PLAN.md` §4 flags as wrong (the AUPRC "across the board" claim
  and two "dcGO ÷ random" cells); the pre-August design docs
  (`IMPLEMENTATION_GUIDE.md`, `SUMMARY.md`, `THINKING.md`) are stale relative
  to the code.
- [ ] Add InterPro names, gene names, and ontology labels to result TSVs (or a
  deterministic annotation companion table) for usable downstream output. The
  production surprise run already emits per-ontology term-name tables
  (`results/production/surprise/names/`) — fold that into the supported path.
- [ ] Remove or formalize `scratch_allspecies/` scripts after the multi-species
  workflow is captured by supported, tested commands; likewise the ad-hoc
  `collect_stats.py` / `gene_collapse_check.py` under
  `results/production/surprise/`.

## Recently completed

- [x] **Paper-parity inference end to end (#61–#70):** relative inference split
  from the True Path flag and generalized to every hierarchical ontology,
  combined with the overall test *before* BH as the paper does, union-of-
  direct-parents background, True-Path input propagation, memory bounded by
  domain block, and the `--min-ic` IC floor with the exported `ic` column.
- [x] **Ontology expansion to 35 registry keys:** HPO + SynGO gene-keyed layers
  (#66), model-organism phenotypes mp/wbphenotype/zfa/fbcv/fbbt (#68), and
  wave 3 mondo/efo/celltype/wbbt/ncit/oncotree (#71) — each with counted
  id-remap policies.
- [x] **Regulates-edge propagation fix (#67).**
- [x] **Production matrix run (2026-08-18):** 63/63 cells, 0 failures, via
  `scripts/run_production_matrix.py` (restartable, per-cell manifests);
  surprise rankings for 59/63 cells.
- [x] Correct contingency tables, Fisher testing, BH families, and True Path
  direction/background; remove the invalid shrinkage heuristic.
- [x] Package and smoke-test the installed `dcgo` CLI; validate inputs and fail
  explicitly when requested hierarchy support is unavailable.
- [x] Add run manifests, `CITATION.cff`, CI coverage, end-to-end tests, and
  committed benchmark artifacts. Current verification: **963 tests pass**.
- [x] Complete InterPro2GO coverage, the 2021→2026 temporal benchmark,
  uncertainty/permutation analysis, component ablation, mouse validation,
  original-dcGO SSF comparison, and 19-ontology breadth evaluation.
- [x] Run the all-species background experiment and make supra-domain inference
  tractable by enumerating only co-occurring combinations while preserving the
  full BH hypothesis count.
- [x] Add and evaluate the surprise score and DOID/Orphanet-to-DOID re-keying;
  retain their negative or inconclusive results rather than treating them as
  successes.
