# TODO

Last reviewed against `main` (`5d3ef91`) on 2026-08-11. This is the short,
prioritized project queue. Detailed experiment designs and results live in
[VALIDATION_PLAN.md](VALIDATION_PLAN.md), [MULTISPECIES_BACKGROUND.md](MULTISPECIES_BACKGROUND.md),
and [REPRODUCIBILITY.md](REPRODUCIBILITY.md); longer-term ontology expansion
lives in [FUTURE_WORK.md](FUTURE_WORK.md).

The core human pipeline is operational, packaged, tested, and validated. The
remaining blockers are chiefly publication design, statistical independence,
and reproducibility rather than basic implementation.

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
- [ ] **Correct phylogenetic non-independence in the multi-species analysis.**
  Protein-level pooling inflates support by about 2.44× and half of associations
  rest on at most three UniRef50 clusters. Test at an ortholog/cluster level (or
  justify a weighting scheme), report pooled and corrected counts side by side,
  and repeat held-out and permutation analyses after correction.
- [ ] **Complete paper-parity inference end to end.** Integrate the relative
  parental-background test into `run_dcgo_human.py`, combine overall and
  relative evidence with the paper's multiple-testing policy, and make the
  calibrated p-score protein-prediction path reproducible without post-hoc
  scripts.
- [ ] **Choose and justify the primary method configuration.** Current ablation
  finds no protein-centric gain from supra-domains and a loss from True Path
  propagation, although supra-domains remain central to the emergent-association
  claim. Pre-specify whether the primary predictor is single-domain/no-TPR and
  treat supra-domain discovery as a separate analysis if so. Set and justify a
  minimum-support/effect-size policy rather than merely exposing
  `--min-support` and confidence intervals.
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
  current score predicts future curation in aggregate but has no demonstrated
  ranking advantage over q-value and leaves 9,923/10,136 associations tied at
  zero.
- [ ] Test MONDO re-keying as a broader alternative to DOID for sparse disease
  annotations. The 2021 DOID breadth run is complete and produced 0/160 hits.
- [ ] Quantify sensitivity to evidence policy, FDR, minimum support, effect-size
  floor, supra-domain length, GO release, and transfer rule; include qualitative
  error analysis of representative successes and failures.
- [ ] Investigate dependence-aware or hierarchical multiple testing across GO
  terms, domains, and nested supra-domain hypotheses.

## Documentation and maintenance

- [ ] Reconcile stale status text. In particular, README still says §3/§4 and
  score calibration are open and reports 525 tests, while the current suite has
  632 passing tests; `VALIDATION_PLAN.md` still has a 2026-07-09 “next steps”
  section and marks completed work as pending; and
  `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` retains resolved launch blockers.
- [ ] Add InterPro names, gene names, and ontology labels to result TSVs (or a
  deterministic annotation companion table) for usable downstream output.
- [ ] Remove or formalize `scratch_allspecies/` scripts after the multi-species
  workflow is captured by supported, tested commands.

## Recently completed

- [x] Correct contingency tables, Fisher testing, BH families, and True Path
  direction/background; remove the invalid shrinkage heuristic.
- [x] Package and smoke-test the installed `dcgo` CLI; validate inputs and fail
  explicitly when requested hierarchy support is unavailable.
- [x] Add run manifests, `CITATION.cff`, CI coverage, end-to-end tests, and
  committed benchmark artifacts. Current verification: **632 tests pass**.
- [x] Complete InterPro2GO coverage, the 2021→2026 temporal benchmark,
  uncertainty/permutation analysis, component ablation, mouse validation,
  original-dcGO SSF comparison, and 19-ontology breadth evaluation.
- [x] Run the all-species background experiment and make supra-domain inference
  tractable by enumerating only co-occurring combinations while preserving the
  full BH hypothesis count.
- [x] Add and evaluate the surprise score and DOID/Orphanet-to-DOID re-keying;
  retain their negative or inconclusive results rather than treating them as
  successes.
