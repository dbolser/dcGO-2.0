> **Historical snapshot.** Active work is tracked in `TODO.md`; detailed code-debt rationale is in `docs/CODE_QUALITY_ROADMAP.md`.

# Engineering and Scientific Review TODOs

Review date: 2026-07-15

This is a read-only review of the implementation and the validation evidence as
they stood on the review date. The overall conclusion is that the core pipeline
is promising and substantially tested, but it is not yet ready to launch as a
packaged tool or to support publication-grade claims of general superiority.

The reported benchmark values appear to be genuine outputs of the implemented
methods. At present, the evidence supports the narrower conclusion that dcGO
contains predictive signal in this retrospective human benchmark. It does not
yet establish robust general performance, calibration, or superiority.

## P0: Software launch blockers

- [x] Fix `run_analysis.sh`. It passes the nonexistent options
  `--disable-true-path` and `--disable-shrinkage`; with `set -e`, the advertised
  eight-configuration analysis exits immediately.
- [ ] Fix and test the installed CLI. `pyproject.toml` declares
  `dcgo = "run_dcgo_human:main"`, but the wheel packages only `src`, not the
  top-level `run_dcgo_human.py`. Move the entry point into a packaged module and
  add a wheel-install smoke test.
- [x] Expand CI to cover the actual product. CI now lints `src/ tests/ config/
  scripts/ validation/ run_dcgo_human.py extract_human_interpro.py`, and already
  builds both distributions and smoke-tests the installed `dcgo` CLI. The
  reported Ruff errors are fixed (the last three, `E402` in
  `validation/sprint1_validation.py`, by moving the logger setup below the
  imports rather than suppressing the rule).
- [x] Add a small end-to-end test that invokes the installed CLI on fixtures and
  verifies schema-valid output. The unit and integration tests do not currently
  demonstrate that the installed command works.
- [ ] Fail rather than silently degrade when requested functionality is
  unavailable. In particular, a run using `--enable-true-path` currently logs an
  error and succeeds without propagation if the ontology is absent.
- [x] Validate CLI parameters: FDR range, positive batch size and core count,
  shrinkage range, supported input/species selection, non-empty protein
  intersection, and non-empty feature spaces.
- [x] Pin every production input using a release identifier, source URL, and
  checksum. Each run now records the SHA-256, byte size, source URL and embedded
  release header of every input the selected ontology consumed (GAF
  `!gaf-version`/`!date-generated`, OBO `data-version`, UniProt vocabulary
  `Release:`, Expasy `Release of`), so a reported number can always be traced to
  the exact bytes behind it. See `REPRODUCIBILITY.md`.
- [ ] Still mutable upstream: recording the hash identifies what was consumed
  but does not make it re-fetchable. Switch the production runs to dated or
  archived releases (`scripts/download_data.py --goa-archive` already fetches
  dated GOA snapshots) and verify the recorded hashes on download.

## P0: Publication blockers

- [ ] Separate model selection from final evaluation. The p-score transfer, IC
  thresholds, shrinkage settings, and other choices appear to have been compared
  using the same 2021-to-2026 benchmark used for the headline result. Freeze
  choices on a development interval or nested temporal split, then evaluate once
  on an untouched interval.
- [x] Add protein-level bootstrap confidence intervals for F-max, AUPRC, and
  paired differences versus baselines. **Done 2026-08-04** —
  `validation/resampling.py` (`paired_bootstrap`, one protein resample per
  replicate shared by every method), driven by `validation/ablation.py`;
  1,000 replicates; results in `validation/ablation_paired_bootstrap.tsv` and
  VALIDATION_PLAN §4. The paired test changes a headline: dcGO's AUPRC advantage
  over naive does **not** survive it at IC≥0 in any aspect.
- [x] Replace the single random-domain shuffle with many seeded permutations.
  Report the null distribution, confidence interval, and empirical p-value.
  **Done 2026-08-04** — `validation/temporal_benchmark.py --n-permutations`
  (default 100; permutation 0 is the old single shuffle) and the same in
  `validation/ablation.py` (200). Null mean/sd/percentile interval and empirical
  p in `validation/temporal_benchmark_permutation_null.tsv` and
  `validation/ablation_permutation_null.tsv`. dcGO clears the null in all 24
  aspect × IC × metric cells at the attainable floor 1/(n+1). The exchangeability
  assumption, what the permutation preserves (term base rates, architectures,
  prediction coverage) and what it cannot test (the Fisher+BH stage) are stated
  in VALIDATION_PLAN §2.
- [x] Complete the planned ablation: single domains; plus supra-domains; plus
  shrinkage; plus True Path Rule; and the full method. **Done 2026-08-04** —
  `validation/ablation.py`; VALIDATION_PLAN §4. **The result is negative:**
  supra-domains 0/12 cells improved, shrinkage 0/12, True Path 12/12 *worse*.
  The best configuration on this benchmark is single domains only.
- [x] Establish the statistical validity of the claimed empirical-Bayes
  shrinkage, or rename it as a heuristic. **Resolved by removal (2026-08-05,
  #44)** — neither of the two options offered, because the step did not survive
  examination well enough to be worth renaming.

  What was measured: the implementation geometrically interpolated each
  supra-domain p-value toward the geometric mean of its constituents' p-values
  using a hand-set decay. Because a well-supported constituent has a far smaller
  p-value than the thin combination containing it, that pulled combinations
  *toward* significance — a 3-protein combination at p=0.01 became p=1e-24 at the
  default strength — and enabling it took the FDR<0.01 count from 163,277 to
  463,924 (+184%), with 56.4% of supra-domain p-values decreased. A shrinkage
  toward a null prior cannot increase rejections.

  On validity: the transformed quantities were never shown to be valid p-values,
  so BH over them was not shown to control FDR. The stronger claim that they
  cannot be valid under *any* null is not something this work established, and is
  not asserted. The §4 ablation separately found the step moved 0 of 12 cells, so
  nothing was lost by removing it. A statistically defensible version would shrink
  the observed *rate* toward a pooled rate and recompute Fisher on the adjusted
  counts; that is a different method and was not attempted.

- [ ] Pre-specify the primary endpoint. The unfiltered result loses to the naive
  method for MF, while the headline superiority depends on IC filtering. The
  IC-filtered analysis is scientifically reasonable, but should be a
  pre-specified secondary analysis or be confirmed on untouched data.
- [ ] Document that IC thresholds can alter the evaluation cohort. Proteins whose
  truth becomes empty are dropped, so IC >= 0, 2, and 4 results do not necessarily
  measure the same proteins. Report cohort sizes and paired analyses explicitly.
- [ ] Pin the exact 2026 t1 GOA snapshot and GO ontology used for propagation.
  The current example pins t0 but uses mutable current files for t1 and GO.
- [ ] Address temporal look-ahead from current InterPro architectures. Either
  reconstruct historical InterPro inputs or state clearly that this is an
  annotation-temporal benchmark rather than a fully prospective simulation.
- [ ] Add stronger external comparators. The naive frequency baseline and
  shuffled-domain null establish signal, but not state-of-the-art utility.
  Include original dcGO output where feasible and at least one independent
  protein-function or domain-based predictor.
- [x] Add an external validation axis: another species, a later untouched time
  interval, or both. *(Mouse, matched 2021→2026 window: performance improves in
  all 9 aspect × IC cells and clears the permutation null in all 24. Nested
  human split 205→215→current in progress.)*

## P1: Engineering and reproducibility

- [x] Generate a machine-readable run manifest containing the Git commit,
  command line, dependency-lock hash, input releases and checksums, timestamps,
  evidence filter, ontology relations, and every threshold.
  `run_manifest_<ontology>.json`, written by `src/run_manifest.py`. The recorded
  inputs are whatever the selected registry entry declares, so all 19 ontologies
  are covered, including their True Path hierarchy files. Not yet done for the
  surprise-score driver or the `validation/` benchmarks.
- [ ] Commit or archive the exact benchmark artifacts used by the manuscript.
  The `bench_A` through `bench_D` outputs and logs were untracked at review time,
  leaving their provenance unclear.
- [ ] Consolidate the README tables, validation plan, committed metrics, and
  benchmark directories into one generated source of truth.
- [ ] Add provenance and interpretation fields to output: input releases, GO
  aspect, evidence policy, contingency cells `a/b/c/d`, direct/propagated status,
  and scoring method.
- [x] Replace `calculate_hypergeometric_score`'s fallback value of `50.0` with an
  explicit failure or missing-value state. A numerical error must not become a
  plausible medium-confidence score.
- [x] Avoid overstating parallelism. `--num-cores` is logged but currently unused
  by the Cython Fisher implementation.
- [ ] Remove unsupported legacy modules and dependencies or bring them under
  testing. The old data-acquisition and database code increases the maintenance
  surface without being part of the supported path.
- [ ] Provide clean-checkout reproduction automation, ideally a container or
  workflow that downloads pinned inputs, verifies hashes, runs inference, and
  regenerates every manuscript table.
- [ ] Add release essentials. `CITATION.cff` is done (software metadata plus the
  Fang & Gough method reference). Still missing: semantic-versioning policy,
  changelog, archived DOI, supported-platform statement, and resource estimates.

## P1: Scientific analysis and reporting

- [x] Report counts at every selection stage: t0 proteins, t1 proteins,
  no-knowledge candidates, proteins with domains, aspect-specific cohorts, and
  cohorts retained at each IC threshold. **Done 2026-08-04** —
  `validation/ablation_selection_counts.tsv`. The concern is confirmed: MF keeps
  only 41% of its cohort at IC≥2, CC only 27% at IC≥6, so cross-floor
  comparisons are not paired. Every paired test in §4 is within one (aspect, IC)
  cell.
- [x] Report prediction coverage alongside F-max. CAFA-style precision can omit
  proteins without predictions while recall includes them. **Done 2026-08-04** —
  `coverage_at_fmax` and `coverage_any` columns of
  `validation/ablation_metrics.tsv`. dcGO covers 42–71% of the cohort; the naive
  baseline covers 100%.
- [ ] Verify F-max and AUPRC against an independent CAFA evaluation
  implementation. The current evaluator samples 51 score quantiles and computes
  trapezoidal AUPRC using an upper envelope; document these choices and check
  them against standard tooling.
- [ ] Report sensitivity to evidence filters, FDR threshold, minimum domain
  support, effect-size floor, supra-domain length, GO release, and transfer rule.
- [~] Introduce and justify a minimum-support/effect-size policy. FDR significance
  alone can retain biologically fragile associations from sparse tables. Report
  contingency cells and odds-ratio confidence intervals.
  *(Cells `a/b/c/d` and a Haldane-corrected Woolf odds-ratio interval are now
  exported per association, so fragility is visible and every row is
  recomputable. The **policy** — an actual minimum-support or effect-size
  cutoff — is still an open method decision, tracked in VALIDATION_PLAN §5.)*
- [ ] Treat InterPro2GO recovery strictly as coverage/recall against an incomplete
  positive reference. It does not measure precision or constitute fully
  independent validation.
- [ ] Discuss dependence between GO terms, domains, and supra-domains. BH is
  applied across a highly dependent hierarchical hypothesis family; simulations
  or hierarchical multiple testing would strengthen the error-control claim.
- [ ] Include qualitative error analysis covering representative true positives,
  false positives, incomplete-domain failures, generic GO terms, and
  multifunctional proteins.

## Verified strengths

- [x] The test suite passed: 162 tests in 6.63 seconds.
- [x] Contingency-table overflow and GO directionality have dedicated regression
  tests.
- [x] The no-knowledge gate uses the same non-IEA evidence space as training,
  reducing direct label leakage.
- [x] GO aspects are evaluated separately, roots are excluded, and truth and
  predictions are propagated consistently.
- [x] InterPro2GO coverage is framed as coverage rather than precision.

## Scientific interpretation

The implementation supports the descriptive claim that, on this retrospective
human split, dcGO achieved the recorded F-max values and outperformed the supplied
naive baseline for BP and CC without IC filtering, and for all three aspects after
filtering to higher-information terms.

It does not yet justify the claim that the complete method is generally superior,
calibrated, controls FDR after shrinkage, or will generalize. A defensible current
summary is:

> Domain-derived associations contain predictive signal for later human GO
> annotations and outperform simple baselines in several retrospective settings,
> particularly for higher-information GO terms.

## Blocker summary

Software launch is blocked by the broken comparison workflow, uncertain packaged
CLI, incomplete CI coverage, and unpinned inputs.

A publication-ready paper is blocked by the absence of an untouched evaluation
set, uncertainty estimates, a repeated random-null analysis, component ablation,
exact reproducibility metadata, and a statistically defensible treatment of the
shrinkage procedure.
