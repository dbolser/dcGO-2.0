# Validation & Benchmarking Plan

This document defines what must be done to turn dcGO-2.0 from a working codebase
into a defensible methods contribution. It is the "(b)" companion to the
engineering cleanup: the code runs; this is about showing the results are
*reliable*.

## Status snapshot (2026-07-09)

| Item | State |
|------|-------|
| §0 Pipeline correctness | ✅ done — 4 correctness bugs found & fixed with tests |
| §1 Reframe InterPro2GO comparison | ✅ done (#14, #15) — ~65% coverage at FDR<0.01 |
| **§2 Temporal held-out benchmark (CAFA-style)** | ✅ **done (#8)** — 2021→2026 CAFA split; see results below |
| §3 Compare to original dcGO | ✅ **done (2026-08-04)** — SSF re-keying + published-dcGO join; precision 0.54–0.63, recall uninterpretable (see §3.1) |
| **§4 Current component ablation + uncertainty** | ✅ **done (2026-09-01)** — nine current-code configurations, 1,000 paired bootstraps, 200 permutations |
| **§5 Pre-paper method decisions** | 🟡 partly done (2026-08-05) — shrinkage removed, True Path background fixed, BH families split, `--min-support` added; minimum-support *policy* still open |
| **Untouched evaluation axis** | ✅ **done (2026-08-05)** — mouse, matched 2021→2026 window; performance improves in all 9 cells. Nested human split (205→215→current) in progress |
| §6 Reproducibility | ⬜ open (#12) |

**Where §2 landed:** with the calibrated **p-score** predictor (default) and a
leak-free no-knowledge gate, dcGO **beats the CAFA naive baseline on F_max at face
value in BP (0.248 vs 0.115) and CC (0.380 vs 0.343)**, and beats it in **every
aspect on informative terms** (IC-filtered), while staying **1.3–25× above the
random-domain null**. MF is the one aspect naive leads at IC≥0 (0.360 vs 0.464),
because its truth is dominated by `protein binding` (84.6% of experimental MF
annotation *lines*, and carried by 87.7% of proteins with any experimental MF
term — 34.5% have nothing else; `validation/protein_binding_dominance.py`) — and dcGO overtakes it decisively the moment that
noise is excluded (IC≥2: 0.365 vs 0.053). This confirms the concern that a
temporal CAFA benchmark rewards recovery of the popularity-weighted curation
frontier; **on informative function, and measured by F_max, dcGO is ahead**.

**All of that is F_max.** AUPRC, computed in the same run, points the other way
at IC≥0 in every aspect (naive BP 0.314 vs 0.137, MF 0.325 vs 0.195, CC 0.513
vs 0.240) and keeps pointing the other way in CC at IC≥2 and IC≥4. So the
claim above is about one metric at one operating point, and the metric was
chosen after both were computed. The defensible summary is the external
review's: domain-derived associations contain predictive signal for later human
GO annotations and outperform simple baselines in several retrospective
settings, particularly for higher-information GO terms. See the AUPRC note in
§2.

> **Update (2026-08-04, §4).** Those §2 numbers now carry uncertainty, a real
> null and a component ablation. Three things changed how they should be read:
> (a) the random-domain null is now 100–200 seeded permutations, and dcGO clears
> it in every aspect × IC × metric cell at the attainable floor of the empirical
> p-value; (b) the F_max advantage over naive is confirmed by a **paired**
> bootstrap, but the **AUPRC** advantage is not — naive wins AUPRC at IC≥0 in all
> three aspects and at every floor in CC; (c) the hierarchy stages separate cleanly: **relative inference is the main
> negative component, output propagation is mostly neutral-to-helpful, and input
> propagation rescues relative-inference configurations.** See §4.

### Next steps (as of 2026-07-09, after the §2 benchmark + method audit)

Done this round: temporal CAFA benchmark (`temporal_benchmark.py`), the paper's
two missing pieces — relative inference (`apply_relative_inference.py`) and the
per-target p-score (`--transfer pscore`, now default) — and a domain-centric
evaluation (`domain_centric_eval.py`). Both baselines cleared on informative
terms. Remaining, in rough priority:

1. **Temporal domain-centric test.** The domain-centric eval currently scores
   against the *current* InterPro2GO. Fetch a dated (~2021) `interpro2go` and pass
   it as `--reference` so it becomes a true held-out temporal test (mirrors §2 on
   the domain side). Small.
2. **Relative inference (paper Step 2).** Implemented; usable **with the
   `--min-ic` floor** (see the sweep below), still opt-in.

   Matching the paper: the parental-background p-value is computed for every
   candidate pair *before* the correction, BH corrects
   `max(overall_p, relative_p)` (an intersection-union statistic, so a valid
   p-value without further adjustment) and the h-score is
   `min(overall, relative)`; the background is the **union** of a term's direct
   parents' proteins, per the Figure 1 caption's `N_pa`; and the input
   protein→term map is true-path propagated, which the paper states explicitly.

   **The defect (now floored away): parentless terms skip the test.** They
   score `relative_p = 0` and pass on the overall inference alone, so the DAG
   roots dominate the output — `biological_process` is the single most
   frequent term on an allspecies run (9,623 domains). Nothing in Steps 1–3
   prevents it; `--min-ic` removes it at reporting time.

   The stated purpose is to report each domain at the most specific level it
   supports. It currently does the opposite (human GO single domains, FDR<0.01):

   | configuration | significant | mean #ancestors | on an ancestor chain |
   |---|---:|---:|---:|
   | overall only | 44,453 | 11.8 | 28.6% |
   | + relative + input propagation | 33,185 | 6.8 | 55.2% |

   **Ruled out:** a wider universe (allspecies is *worse* — 82.4% of
   associations on a chain, the three GO roots as the top three results); the
   attainable-p floor (binds only below ~100-protein backgrounds); degenerate
   parents whose propagated protein set equals their child's (6.3% of
   (term, parent) pairs). The per-parent-maximum background was also wrong and
   is fixed, which improved specificity but not the cascade.

   **Done (2026-08-17): the information-content floor** — dcGO's own
   IC/"meta-GO" analysis (their Figure 3) is plausibly where level selection
   actually happens. `--min-ic` floors reported terms on the shared
   annotation-frequency IC (`src/information_content.py`, estimated from the
   propagated map, exported as the `ic` column), applied after BH exactly like
   `--min-support` so the hypothesis family is untouched. Sweep on human GO
   single domains, paper-parity config (`--propagate-annotations
   --enable-relative-inference --enable-true-path`, post-regulates-fix;
   computed with `validation/specificity_metrics.py`, which makes the ad hoc
   chain numbers above reproducible; re-derived unchanged from a
   full-precision `ic` export after the column moved from 4 decimals to
   `%.10g`, so no row sits on a rounding boundary):

   | --min-ic | significant | mean #ancestors | on an ancestor chain | GO roots |
   |---:|---:|---:|---:|---|
   | 0 (off) | 30,655 | 6.0 | 52.7% | all three |
   | 1 | 30,302 | 6.1 | 50.2% | none |
   | 2 | 28,348 | 6.3 | 46.4% | none |
   | 3 | 26,401 | 6.4 | 42.8% | none |
   | 5 | 18,888 | 7.3 | 33.3% | none |

   The GO roots sit at IC 0.09–0.17 (each aspect covers most but not all of
   the universe), so any floor ≥ 0.2 removes them; **IC ≥ 1** — keep terms
   carried by under half the universe — also removes the near-universal band
   just above them ("biological regulation" ×84, "cytoplasm" ×70, "organelle"
   ×57, "cellular process", "binding") at a cost of 353 associations (1.2%).
   That is the recommended setting: the floor's job is killing vacuous terms,
   not making the method conservative. What it does *not* fix is the residual
   cascade — half the surviving associations still sit on an ancestor chain —
   so the remaining idea if that matters is `elim`-style decorrelation (Alexa
   et al. 2006): test specific terms first, remove their proteins, retest
   ancestors on the residue.

   *(The p-score protein-prediction path is deliberately **not** pursued: the
   deliverable is domain→term annotation, not protein-centric prediction.)*

3. **§4 ablation** — now that the yardstick is trusted, run the temporal benchmark
   per pipeline config (single-domain / +supra / +shrinkage / +TPR).
4. ~~**§3 original-dcGO comparison**~~ — **done** (2026-08-04), see §3.1.
   `--domain-key ssf` re-keys on the SUPERFAMILY/SCOP signature; the remaining
   §3 work is the `pfam` key and the non-GO reference tables.
5. **§5 method decisions** — TPR default on/off; species scope (human-only vs
   multi-species); minimum-evidence / effect-size floor; Haldane-corrected odds
   ratio. **§6 reproducibility** — pin GOA/InterPro/GO versions, one-command
   repro, archive the exact input snapshots.

Report BP/MF as the headline (as the original did; CC is least domain-relevant),
and keep **both** protein-centric and domain-centric evaluations.

---

## 0. Pipeline correctness (done — the results can now be trusted)

Before any validation number means anything, the underlying associations have to
be computed correctly. Four correctness bugs were found and fixed (each with a
regression test), so this section is now **closed**:

1. **Fisher-tail vs. hypergeometric disagreement** was the symptom that exposed
   three contingency-table bugs (#17):
   - **int8 overflow** — overlap counts computed as an int8 sparse matmul wrapped
     at 127 (true `a=300` → `44`).
   - **Non-binary matrices** — a `(protein, InterPro)` pair is listed once per
     member signature in `protein2ipr`, so duplicates were summed (matrix cells
     reached 125), inflating overlaps and driving the `d` cell negative.
   - **occurrence-vs-protein counts** — `observation_count` counted domain
     occurrences, not proteins (reported `n_observations` exceeded the proteome).
2. **True Path Rule DAG direction** (#15) — `OntologyProcessor` traversed the GO
   graph the wrong way (obonet emits child→parent), so `get_ancestors` returned
   descendants and the optimal-level filter tested against children.

These were caught by adversarial code review (Gemini + Codex) on the first
results PR, *after* an initial run reported an inflated ~79% coverage. Lesson
retained: sanity-check the extreme rows (`odds_ratio`/`hyper_score`/`p` triples)
before quoting any number.

**Why the naive InterPro2GO comparison was also misleading** (the original §0,
now addressed by §1): it treated the deliberately-incomplete, positive-only
InterPro2GO map as complete ground truth, without GO-DAG propagation, across
domains InterPro2GO does not even cover — so "absence" was miscounted as "wrong"
and precision was meaningless. §1 fixed the *comparison*; §0 fixed the *pipeline*.

---

## 1. Fix and reframe the InterPro2GO comparison  *(quick win)*

Treat InterPro2GO as an **incomplete positive reference**, and report only what
it can legitimately support.

- [x] De-duplicate the threshold sweep; verify each threshold changes the
      prediction set (the current file has identical rows — a bug).
- [x] **Propagate both** predictions and the reference to their GO ancestor
      closure before intersecting (use `ontology_processor` + `go-basic.obo`).
- [x] Report **recall / coverage of InterPro2GO** as the headline (what fraction
      of curated pairs we recover), and explicitly label "novel" pairs as
      *not-in-reference* rather than *false*.
- [x] Restrict the comparison to the **domains actually present in both** sets
      (a domain absent from InterPro2GO contributes only noise).

**Status: done** (PRs #14, #15). Implemented in `validation/validate_results.py`
with unit tests in `tests/unit/test_validation_metrics.py`.

### Results — first human run (2026-07-08, after the contingency-table fixes #17)

> An earlier version of this run reported ~79% coverage over 457,939 significant
> associations. That was **corrupted** by three contingency-table bugs (int8
> overflow, non-binary matrices, occurrence-vs-protein counts) caught in review
> and fixed in #17. The numbers below are the corrected run.

Full human pipeline: 165,823 significant associations (FDR<0.01). InterPro2GO
reference: 30,190 curated pairs → 47,393 after propagation, on the 2,747 domains
shared with our predictions.

| Threshold | Reference coverage (recall) | Recovered / 47,393 |
|-----------|-----------------------------|--------------------|
| p ≤ 1e-10 | 29.3% | 13,881 |
| p ≤ 1e-8  | 42.3% | 20,051 |
| p ≤ 1e-6  | 64.7% | 30,673 |
| **FDR < 0.01 (significance)** | **64.7%** | 30,673 |

**Headline:** at the FDR<0.01 significance cutoff dcGO recovers ~65% of curated
InterPro2GO associations (on shared domains, propagated) — versus the old,
misleading "3–4% precision" that treated InterPro2GO as complete truth without
propagation. This is the reframing working, not a new result about accuracy.

**Caveats (do not oversell):**
- `precision_lower_bound` ≈ 22–26% is a *floor*, not precision. The ~107k
  "candidate" pairs are curation-gap candidates (InterPro2GO is incomplete), not
  errors. Real precision needs the §2 temporal/CAFA benchmark.
- The loose end of the sweep plateaus because the input is already FDR<0.01
  filtered; the informative signal is the *tightening* end (1e-4 → 1e-10).

---

## 2. Temporal held-out benchmark (CAFA-style)  *(core evidence — DONE)*

This is the single most important addition. It measures whether dcGO predicts
annotations that were **later** confirmed — a real, precision-capable test.

**Status: done.** Implemented in `validation/temporal_benchmark.py` (pure metric
functions, unit-tested in `tests/unit/test_temporal_benchmark.py`); dated GOA
snapshots fetched via `scripts/download_data.py --goa-archive <version>`.

- [x] Obtain two dated GOA snapshots — **t0 = release 205 (2021-04-21)**,
      **t1 = current (2026-06)**, a ~5-year gap. EBI dated GOA archive.
- [x] Train domain→GO associations using only annotations available at `t0`
      (ran the standard pipeline on the 2021 GOA → 164,549 significant
      associations at FDR<0.01).
- [x] Define the benchmark as **no-knowledge proteins per aspect**: proteins with
      *no annotation known to training in that aspect at t0* that gained
      experimental annotation by t1. Score against their **full** propagated t1
      experimental terms. The gate uses the **same evidence filter the pipeline
      trains on** (`manual`/non-IEA), not experimental-only — gating on anything
      narrower leaks already-seen labels into the held-out set (a bug caught in
      review; fixed). (An earlier delta-only truth — `t1 minus t0` — was also
      wrong: it scored correct predictions of already-known terms as
      misinformation. Both fixed before any number was quoted — the §0 lesson.)
- [x] Score with the **CAFA protein-centric metric**: transfer predicted GO terms
      to each protein via its domains, then **F_max** over a threshold sweep,
      plus **S_min** (marginal-IC weighted, IC from t0) and **AUPRC**. Rank by
      −log10(p) — `hyper_score` saturates (37% at exactly 100) and collapses the
      sweep. Default transfer is the **p-score** (Fang & Gough: sum of scores,
      min-max normalised per protein); `--transfer max` for the simpler variant.
- [x] Report separately for BP / MF / CC.
- [x] Sweep an **information-content floor** (`--min-ic`) that excludes
      near-universal, low-IC terms from truth and all methods alike — the fair,
      principled way to stop rewarding base-rate recovery of terms like
      GO:0005515 `protein binding`, near-zero IC. Quantified by
      `validation/protein_binding_dominance.py` (metrics:
      `validation/protein_binding_dominance.tsv`): it is **84.6%** of
      experimental MF annotation *lines*, but the measure that explains the
      benchmark is protein coverage — **87.7%** of the 15,260 proteins with any
      experimental MF term carry it, and **34.5%** carry nothing else, so a
      baseline that always predicts it is right about the entire truth for a
      third of the cohort.

### Results — 2021→2026 temporal split (2026-07-09, p-score transfer)

No-knowledge benchmark sizes (IC≥0): **BP 324 / MF 418 / CC 572** proteins (a
leak-free gate on training evidence — much smaller and cleaner than the earlier
experimental-only gate, which wrongly admitted proteins whose t0 computational
labels the model had already seen).

**dcGO (p-score) beats the naive baseline at face value on BP and CC, and beats
it in every aspect once uninformative terms are excluded**, staying well above the
random-domain null throughout:

| Aspect | IC floor | dcGO F_max | naive F_max | random F_max | dcGO / random | dcGO AUPRC | naive AUPRC |
|--------|:--------:|-----------:|------------:|-------------:|--------------:|-----------:|------------:|
| BP | ≥0 | **0.248** | 0.115 | 0.158 | 1.6× | 0.137 | **0.314** |
| BP | ≥2 | **0.170** | 0.071 | 0.053 | 3.2× | **0.069** | 0.049 |
| BP | ≥4 | **0.115** | 0.031 | 0.019 | 6.1× | **0.032** | 0.012 |
| BP | ≥6 | **0.077** | 0.010 | 0.003 | 24× | **0.017** | 0.003 |
| MF | ≥0 | 0.360 | **0.464** | 0.262 | 1.4× | 0.195 | **0.325** |
| MF | ≥2 | **0.365** | 0.053 | 0.088 | 4.2× | **0.227** | 0.025 |
| MF | ≥4 | **0.337** | 0.045 | 0.072 | 4.7× | **0.210** | 0.011 |
| MF | ≥6 | **0.217** | 0.018 | 0.009 | 25× | **0.121** | 0.006 |
| CC | ≥0 | **0.380** | 0.343 | 0.291 | 1.3× | 0.240 | **0.513** |
| CC | ≥2 | **0.239** | 0.153 | 0.072 | 3.3× | 0.073 | **0.117** |
| CC | ≥4 | **0.134** | 0.099 | 0.031 | 4.3× | 0.029 | **0.038** |
| CC | ≥6 | **0.124** | 0.044 | 0.015 | 8.1× | **0.025** | 0.012 |

(IC in bits; IC≥2 ⇒ term in ≤25% of proteins, IC≥6 ⇒ ≤1.6%. Full table with
S_min in `validation/temporal_benchmark_metrics.tsv`.)

> **The AUPRC columns disagree with the F_max columns, and that has to be said
> up front.** At IC≥0 the naive baseline has the higher AUPRC in *all three*
> aspects — including BP and CC, the two where dcGO wins F_max — and in CC it
> keeps the higher AUPRC at IC≥2 and IC≥4 as well, where dcGO's F_max lead is
> otherwise clear. Only in BP and MF does dcGO take AUPRC once the IC floor
> rises. So "dcGO beats naive" is a statement about F_max at a chosen operating
> point, not a statement about the whole precision-recall curve, and CC is the
> aspect where the two metrics most persistently disagree. The primary endpoint
> was never pre-specified (a P0 item in
> `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md`), so choosing F_max after seeing both
> is exactly the freedom that review objects to. These AUPRC numbers were
> present in `temporal_benchmark_metrics.tsv` from the start but were not
> surfaced in any summary until 2026-08-04.

**Read-out:**
- **dcGO beats naive on F_max at IC≥0 in BP (0.248 vs 0.115) and CC (0.380 vs
  0.343)**; MF is the only aspect where naive leads at IC≥0 (0.360 vs 0.464) —
  because MF truth is dominated by `protein binding`. At IC≥2 dcGO wins MF
  decisively (0.365 vs 0.053) and every aspect thereafter.
- **naive is a base-rate mirage.** Its F_max collapses toward the random null as
  informative terms are required (BP 0.115→0.010, MF 0.464→0.018, CC 0.343→0.044
  over IC 0→6): all it ever had were high-frequency generic terms.
- **dcGO stays 1.3–1.6× the random null at IC≥0, rising to 8–25× at IC≥6** — it
  degrades gracefully because its predictions are specific. On MF it even holds
  flat 0→2 (0.360→0.365) as protein-binding noise is removed.
- This confirms Julian's point: a temporal CAFA benchmark rewards recovery of the
  **popularity-weighted** curation frontier, so raw F_max flatters a frequency
  baseline. On informative terms — what a domain→function method is for — dcGO is
  clearly ahead. The calibrated p-score (default) is what lifts it above naive
  even at IC≥0.

**Acceptance: met on F_max, qualified on AUPRC.** dcGO clears both mandatory
baselines (random null and naive) on F_max across all three aspects on
informative terms, and beats naive at face value on 2 of 3 aspects. It does
*not* clear naive on AUPRC at IC≥0 in any aspect, nor in CC at IC≥2 or ≥4. The
defensible summary is the review's: domain-derived associations contain
predictive signal for later human GO annotations and outperform simple
baselines in several retrospective settings, particularly for higher-information
GO terms — not that the method is generally superior.

### Method audit vs the original dcGO Predictor (2026-07-09)

Read of Fang & Gough 2013, *A domain-centric solution to functional genomics via
dcGO Predictor* (BMC Bioinformatics 14(S3):S9 — the CAFA paper) and the database
paper (`gks1080.pdf`). **Their validation is protein-centric CAFA PR-RC**
(per-target precision/recall, averaged over targets), so our protein-centric
choice matches theirs. But three concrete divergences likely explain our raw
gap to naive — and two of them are *method*, not metric:

| Dimension | Original (Fang & Gough) | Ours (§2 run) | Likely impact |
|-----------|-------------------------|---------------|---------------|
| Aspects scored | **BP + MF only** (found MF > BP; CC not scored — least domain-relevant) | BP + MF + CC | report BP/MF as headline, CC secondary |
| Gold-standard evidence | newly-added **EXP + TAS + IC** | strict experimental set | minor |
| **Statistical inference** | **two tests**: *overall* (background = all UniProt) **and** *relative* (background = only proteins with all direct-parent terms); keep the **larger p / smaller h-score** | **overall only** | **big** — the relative/parental test is what enforces specificity; omitting it lets generic low-IC associations through |
| **Prediction transfer / score** | **sum** of h-scores over supporting domains, then **min-max normalised _per target_** to 0–1 | **max** over domains, ranked by global −log10(p) | **big** — per-target normalisation is the per-protein calibration we lacked |
| Propagation / true-path | always on | opt-in, **off** in the run | medium |
| Domain universe | SCOP/SUPERFAMILY (+Pfam) | InterPro entries | see §3 |
| FDR threshold | < 1e-3 | < 0.01 | minor |
| Score saturation | **explicitly noted**: FDR-based p-scores "collapse to FDR=0" for top hits → switched to h-score | we hit the same with `hyper_score` (37% at 100) → switched to −log10(p) | corroborates our fix |

**Method-piece experiment (2026-07-09, leak-free gate).** We restored both
missing pieces and measured each (`apply_relative_inference.py`;
`temporal_benchmark.py --transfer`). dcGO F_max / AUPRC at IC≥0 (naive: BP
0.115/0.314, MF 0.464/0.325, CC 0.343/0.513):

| Config | BP F_max | MF F_max | CC F_max | BP AUPRC | MF AUPRC | CC AUPRC |
|--------|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| A base (max, overall) | 0.218 | 0.234 | 0.349 | 0.083 | 0.106 | 0.201 |
| B **+ p-score transfer** | **0.248** | 0.360 | **0.380** | **0.137** | 0.195 | **0.240** |
| C + relative inference | 0.219 | 0.273 | 0.356 | 0.065 | 0.164 | 0.180 |
| D both | 0.224 | **0.385** | 0.356 | 0.100 | **0.217** | 0.196 |

**The per-target p-score (2) is the main lever; the relative inference (1) helps
MF.** Per-protein calibration (B) lifts F_max/AUPRC across the board (BP
0.218→0.248, MF 0.234→0.360, CC 0.349→0.380). The relative inference (C) helps
MF (F_max 0.234→0.273, AUPRC 0.106→0.164) and is ~neutral on BP/CC F_max; the two
together (D) give the best MF (0.385). So (2) is the broad protein-centric win and
(1) adds specificity that shows up most in MF — complementary, as their
domain-centric result below confirms.

**Correction to an earlier claim:** the relative inference is *not* equivalent to
the evaluation-time IC filter. The IC filter changes the *scoring universe* (only
informative terms count, so naive collapses); the relative inference only prunes
dcGO's *predictions* while the evaluation still scores every term — a different
operation with a much smaller effect. Its clearest payoff is association-set
quality — best seen in the **domain-centric** evaluation below, not protein-centric
F_max. Recommended default: **p-score transfer** (now the default).

### Domain-centric evaluation — where the relative inference earns its keep

`validation/domain_centric_eval.py` scores the **(domain, GO) associations
directly** against InterPro2GO (propagated, shared single-domain space — the §1
frame). This is the right instrument for the relative inference, whose job is
association-set specificity. On the 2021 t0 associations:

| Config | pairs | precision (lb) | recall | F1 |
|--------|------:|---------------:|-------:|---:|
| base (overall only) | 134,610 | 0.218 | 0.631 | 0.324 |
| **+ relative inference** | 69,206 | **0.253** | 0.430 | 0.319 |

The relative inference **raises domain-level precision 0.218 → 0.253** (a real
gain — these are curated-pair rates) while halving the set (recall 0.63 → 0.43).
So the two restored pieces earn their keep on **different** metrics, and are
complementary: **(1) relative inference** cleans the *associations* (domain-centric
precision ↑), **(2) p-score transfer** calibrates the *protein* predictions
(protein-centric F_max/AUPRC ↑). Keep both evaluations; keep both method pieces.

(Caveat: precision here is a *lower bound* vs the incomplete, *current*
InterPro2GO — good for comparing configs, not an absolute precision. A dated
2021 InterPro2GO would make this a fully temporal domain-centric test — a small
follow-up: fetch an archived `interpro2go` and pass it as `--reference`.)

### Baselines
- [x] **Naive baseline**: predict each GO term at its (propagated) t0 frequency
      (CAFA's standard `Naive`).
- [x] **Random-domain null — now a distribution, not one draw** (2026-08-04).
      See "Permutation null" immediately below.
- [ ] **BLAST/annotation-transfer baseline** *(optional, still open)*: transfer
      GO terms from the most similar annotated protein.
- [ ] **Full-pipeline shuffle** *(stronger null, still open)*: shuffle labels and
      re-run Fisher end-to-end to confirm the FDR itself is calibrated (the
      current shuffle is at the association level, not the whole pipeline).

### Permutation null (2026-08-04) — replaces the single shuffle

The published "dcGO ÷ random" ratios came from **one** seeded shuffle. That is an
anecdote: it gives no spread, no interval and no p-value, and a single draw can
land anywhere in the null. `validation/temporal_benchmark.py` now takes
`--n-permutations` (default 100; permutation 0 uses `--seed`, so the old
`random_domain` row is the first draw of the reported distribution) and writes
`validation/temporal_benchmark_permutation_null.tsv`.

**What is permuted, and the exchangeability assumption.** The null hypothesis is
*"a domain's identity carries no information about the functions of the proteins
that carry it"*. The randomisation that hypothesis licenses is a permutation of
the labels of the domain→GO map, and nothing else. It **preserves** every
protein's real architecture, the multiset of GO-sets (so each GO term's marginal
frequency across the table is exactly preserved — the null still predicts
`protein binding` as often as the real method does), the number of domains making
any prediction, and therefore the set of proteins that receive any prediction at
all. It **destroys** the pairing between a domain and its terms, and with it the
correlation between a domain's prevalence and the size/specificity of its term
set. It is a *base-rate-preserving* null, which is why it scores in the same
range as the naive baseline at IC≥0 rather than at zero.

**What it does not test — stated plainly.** It permutes the *surviving,
FDR-significant* table, so it inherits the real pipeline's decisions about which
domains are predictive at all and how many terms each gets. It is a null for the
**transfer step**, not for the whole method, and it cannot say whether Fisher+BH
is calibrated. The end-to-end label shuffle above is still open and would be a
harder null, because permuted training labels would yield far fewer significant
associations rather than the same number of scrambled ones.

F_max, 100 permutations, seeds 0–99, `--transfer pscore`:

| Aspect | IC | observed | null mean ± sd | null 95% pct | obs ÷ null mean | *z* | empirical *p* |
|---|--:|--:|--:|--:|--:|--:|--:|
| BP | ≥0 | 0.2484 | 0.1552 ± 0.0086 | [0.140, 0.173] | 1.6× | 10.8 | 0.0099 |
| BP | ≥2 | 0.1701 | 0.0536 ± 0.0070 | [0.042, 0.068] | 3.2× | 16.7 | 0.0099 |
| BP | ≥4 | 0.1153 | 0.0167 ± 0.0050 | [0.009, 0.027] | 6.9× | 19.9 | 0.0099 |
| BP | ≥6 | 0.0774 | 0.0057 ± 0.0031 | [0.002, 0.013] | 13.7× | 23.4 | 0.0099 |
| MF | ≥0 | 0.3597 | 0.2010 ± 0.0198 | [0.166, 0.240] | 1.8× | 8.0 | 0.0099 |
| MF | ≥2 | 0.3654 | 0.0488 ± 0.0120 | [0.032, 0.077] | **7.5×** | 26.3 | 0.0099 |
| MF | ≥4 | 0.3367 | 0.0373 ± 0.0102 | [0.023, 0.060] | **9.0×** | 29.2 | 0.0099 |
| MF | ≥6 | 0.2168 | 0.0097 ± 0.0053 | [0.003, 0.022] | 22.3× | 39.0 | 0.0099 |
| CC | ≥0 | 0.3796 | 0.2758 ± 0.0149 | [0.242, 0.299] | 1.4× | 7.0 | 0.0099 |
| CC | ≥2 | 0.2385 | 0.0708 ± 0.0102 | [0.058, 0.096] | 3.4× | 16.5 | 0.0099 |
| CC | ≥4 | 0.1345 | 0.0302 ± 0.0078 | [0.017, 0.047] | 4.5× | 13.3 | 0.0099 |
| CC | ≥6 | 0.1245 | 0.0176 ± 0.0070 | [0.006, 0.032] | 7.1× | 15.3 | 0.0099 |

(AUPRC in the same file; its ratios are larger still — 1.5× to 500×.)

**Read-out.**
- **The result holds, and holds strongly.** In all 24 aspect × IC × metric cells,
  *no* permutation out of 100 ever reached the observed value, so every empirical
  *p* is at the attainable floor of 1/(n+1) = 0.0099. With 100 permutations that
  is the smallest p obtainable; it is a limit of the design, not a measurement of
  0.0099, and is reported as such.
- **Two published ratios were understated by a single unlucky draw.** The one
  shuffle behind `RESULTS.md` gave MF IC≥2 random = 0.088 and MF IC≥4 = 0.072 —
  both in the *upper tail* of the null (95th percentiles 0.077 and 0.060; maxima
  0.098 and 0.088). Against the null **mean** the correct ratios are **7.5× and
  9.0×**, not the reported 4.2× and 4.7×. BP ≥0 (1.6×) and CC ≥0 (1.3× → 1.38×)
  were essentially unaffected. **`RESULTS.md`'s "1.3–25×" range should be read as
  "1.4–22× against the null mean, with the low end at IC≥0 and the high end at
  IC≥6"; the individual per-cell ratios in that table are single draws and the
  table above supersedes them.**
- The null is **not** trivially weak. At IC≥0 it scores 0.155–0.276 F_max,
  because it preserves each term's marginal frequency even though it destroys
  every domain→term pairing. In **BP** it is *above* the naive baseline (0.155
  vs 0.115); in MF and CC naive is still ahead (0.201 vs 0.464, 0.276 vs
  0.343). dcGO's margin over the null grows with the IC floor exactly as the
  specificity argument predicts.

  That BP row is worth stating plainly, because it is the sharpest available
  measurement of what a CAFA-style F_max actually rewards: a shuffle carrying
  **no domain information whatsoever** outscores the standard baseline, purely
  by reproducing how often each term is annotated. Any method evaluated this
  way is being graded substantially on base-rate recovery.

### Implementation sketch (concrete next steps)

**Data.** EBI keeps dated GOA archives at
`https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/` (files like
`goa_human.gaf.<version>.gz`, plus a `submitted/` history). Pick `t0` and `t1`
≥ ~2 years apart (e.g. a 2021 release vs the current one) so enough newly-curated
annotations accumulate. Extend `scripts/download_data.py` with a `--goa-archive
<version>` option and record exact versions (feeds §6).

**Split.** Domain architectures come from `protein2ipr` and are *not* time-
varying in our data, so hold them fixed; only the GOA annotations move in time.
- Train set: `(protein, GO)` present at `t0`, EXP/curated evidence.
- Eval set: `(protein, GO)` present at `t1` but **not** at `t0` (newly curated),
  EXP evidence, propagated to ancestor closure. Exclude proteins with zero
  domains (unpredictable by any domain method).

**Predict.** Run the (now-correct) pipeline on `t0` to get domain→GO
associations, then transfer to each eval protein via its domains: a protein's
predicted GO set = union over its domains of that domain's associations, scored
by the max `hyper_score` (or min FDR) across contributing domains.

**Score (protein-centric CAFA).** For a sweep of score thresholds τ, compute
per-protein precision/recall against the eval set (propagated), average over
proteins, then take **F_max** = max over τ of the harmonic mean. Also compute
**S_min** (information-content weighted, IC from `t0` term frequencies) and
**AUPRC**. Report per aspect (BP/MF/CC) separately — mixing them is misleading.

**Where it lives.** New `validation/temporal_benchmark.py` reusing
`ontology_processor.get_ancestors` (propagation) and the `validate_results.py`
helpers; new `tests/unit/test_temporal_benchmark.py` for the metric maths on a
synthetic split (mirror the §1 test approach — pure functions, tiny fixtures).

**Reuse note.** The propagation + shared-space logic from §1 is already on
`main`; the CAFA metric (F_max/S_min) is the genuinely new code here.

---

### Breadth: does the predictive signal hold beyond GO? (2026-08-04, corrected) — DONE

> **These numbers replace the 2026-07-28 table, which was wrong.** The earlier
> run scored association tables built before `restrict_to_universe` (#26), when
> a `--species human` run of a UniProt-native ontology built its contingency
> tables over the whole of Swiss-Prot rather than the human proteome —
> `protein2ipr` is extracted per species, `uniprot_sprot.dat` is not. The
> t0 universes were inflated 3–70×, and every non-GO association set was
> substantially wrong:
>
> | Ontology | t0 universe before | after | t0 associations before | after |
> | --- | ---: | ---: | ---: | ---: |
> | reactome | 35,849 | 10,736 | 93,172 | 59,426 |
> | keyword | 555,986 | 18,797 | 166,924 | 83,739 |
> | subcellular | ~362k | 16,750 | 22,106 | 10,867 |
> | complex | — | — | 4,482 | 1,385 |
> | cofactor | ~124k | 1,801 | 2,284 | 1,725 |
> | disease | — | — | 53 | 53 |
> | **go** (anchor) | 18,735 | 18,382 | 164,549 | 163,277 |
>
> **The contamination was suppressing the signal, not creating it.** Every
> enrichment below is *higher* than the superseded figure. GO is the control
> that shows the correction is real rather than a change of protocol: GOA is
> already species-specific, so its universe moves 1.9% and its enrichment does
> not move at all (11.3× → 11.5×).

`validation/temporal_breadth.py` applies the §2 split to the ontologies added in
`src/ontology_registry.py`. One archived Swiss-Prot release (**2021_02**,
07-Apr-2021 — the same month as GOA release 205) supplies t0 for every
UniProt-native layer at once, so seven ontologies are trained and scored under
one protocol. Metrics: `validation/temporal_breadth_metrics.tsv`,
`validation/temporal_breadth_go.tsv`.

Statistic as in §2: predictions are proteins carrying a feature that lacked the
term at t0; hits are those carrying it at t1; enrichment is the hit rate over the
term's own acquisition rate. GO is included **under the identical protocol** as
the anchor — without it, "weaker than GO" could not be said, because the 12.5×
figure came from a different subset (supra-domains only) and a stricter truth
(experimental evidence only).

| Ontology | assoc. | predictions | hits | hit rate | expected | enrichment (95% CI) | superseded |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| **go** (anchor) | 91,830 | 3,106,234 | 106,224 | 3.42% | 0.30% | **11.5× [11.1, 12.0]** | 11.3× |
| reactome | 34,397 | 760,073 | 1,430 | 0.19% | 0.02% | **11.4× [8.9, 14.6]** | 8.0× |
| cofactor | 291 | 1,404 | 235 | 16.74% | 4.72% | **3.5× [2.1, 4.5]** | 3.2× |
| subcellular | 6,340 | 172,679 | 3,798 | 2.20% | 0.59% | **3.7× [3.5, 3.9]** | 2.9× |
| keyword | 35,981 | 889,117 | 13,140 | 1.48% | 0.44% | **3.4× [3.2, 3.6]** | 1.7× |
| complex | 958 | 9,904 | 2 | 0.02% | ~0% | *artefact, see below* | 265× |
| disease | 27 | 90 | 0 | 0% | ~0% | *undefined* | undefined |
| ligand | — | — | — | — | — | *not testable, see below* | — |

(In-scope universe where it differs from all-domains; keyword and GO annotate
essentially every protein, so the two coincide.)

**The signal generalises, and Reactome is no longer behind GO.** Every ontology
with enough data enriches above 1 with the interval excluding it — dcGO's
association-finding is not a GO artefact. The correction changes one conclusion
outright: Reactome at **11.4× [8.9, 14.6]** now overlaps GO at **11.5×
[11.1, 12.0]**, so the earlier "GO remains strongest and Reactome follows" is
not supported. Reactome's interval is much wider (1,430 hits against GO's
106,224), so the honest statement is that the two are *indistinguishable at this
power*, not that Reactome matches GO.

The remaining ordering — subcellular 3.7×, keyword 3.4×, cofactor 3.5× — is a
band well clear of 1 but well below the two pathway/function vocabularies.
Keywords are no longer the outlier they appeared to be (1.7× → 3.4×): that gap
was contamination, not a property of the vocabulary. Their base rate is still
the highest of any layer (0.44%), which is what caps the achievable ratio.

**Three results that must not be quoted at face value:**

- **`complex` was an artefact, and the correction exposed it.** The superseded
  run reported 265× on 24 hits and an interval excluding 1, which was already
  flagged as a degenerate ratio (ComplexPortal averages ~1.5 proteins per
  complex, so the base rate rounds to zero and any hit divides by nearly
  nothing). With the contaminating proteins removed it falls to **2 hits** on
  9,904 predictions and the interval becomes **[0.00, 130.93]** — it now
  includes zero. There is no demonstrated signal in this layer at all. The
  earlier "detectable, magnitude meaningless" was too generous: it was not
  detectable either.
- **`disease` scored 0 hits on 369 predictions** — no signal, and enrichment is
  undefined because expected hits ≈ 0 as well. This is the 53-association t0 run
  playing out as predicted: 6,904 OMIM phenotypes over 5,029 proteins is too thin
  for contingency tables. The queued Disease-Ontology re-keying (see `TODO.md`)
  is the concrete fix — pool sparse OMIM terms into DO classes *before* testing.
- **`ligand` cannot be tested at this split at all.** UniProt introduced the
  structured `FT /ligand_id="ChEBI:…"` qualifier *after* 2021; the April 2021
  release used free-text `/note="ATP"`. The layer is therefore entirely post-2022
  annotation. Testing it needs a later t0 (2022_05 or 2023_01), which would be a
  shorter, non-comparable split.

**Controls and caveats.**

- **In-scope control.** A sparse ontology could enrich merely because its
  proteins are the ones that get curated at all, so every ontology is also
  scored against only the proteins it reaches by t1. `cofactor` is the sharp
  case: the hit rate rises 1.24% → 12.67% while enrichment holds at ~3× with the
  interval still excluding 1, so that signal is term-specific rather than
  scope-driven.
- **No evidence filter exists off GO.** GO trains on non-IEA and scores against
  experimental only; the UniProt layers have no such split, so an automated
  annotation added between snapshots counts as a hit. The anchor row bounds how
  much that matters: GO scores **11.5×** under the loose protocol against
  **12.5×** under the strict one, i.e. the looser truth did not inflate GO — but
  that is one ontology's worth of reassurance, not a general guarantee.
- Same look-ahead caveat as §2: architectures come from the current
  `protein2ipr`, so this is annotation-temporal, not prospective.

- [ ] **Open:** re-test `ligand`/`cofactor` at a 2023 t0, once the structured
      binding annotation exists in both snapshots.
- [ ] **Open:** repeat with `--supra-only` to isolate the emergent claim per
      ontology (the numbers above pool single domains and supra-domains).

---

## 3. Comparison to the original dcGO  *(reproducibility)*

A method claiming to implement dcGO must relate its output to the published one.

**What the original did** (from `docs/gks1080.pdf`, `docs/1471-2105-14-S3-S9.pdf`):
all completely-sequenced genomes (2,414 genomes + UniProt, >80M sequences);
**SCOP superfamily/family** domains via SUPERFAMILY HMMs, plus Pfam — *not*
InterPro directly; hypergeometric ≈ Fisher with BH-FDR at **FDR < 1e-3**
(stricter than our default 0.01); true-path propagation. Our run is human-only,
InterPro-entry-keyed, FDR<0.01.

**The identifier mapping is easier than it looks.** InterPro *integrates* both
SUPERFAMILY (SCOP-based) and Pfam as member databases, and `protein2ipr` column 4
carries the member signature (`SSFxxxxx` = SUPERFAMILY/SCOP, `PFxxxxx` = Pfam)
for every match. So this is not a disjoint domain universe.

**Resource check (2026-07-28).** Both dcGO websites are live, so this task is
unblocked:

| Resource | Status |
| --- | --- |
| `http://www.protdomainonto.pro/dcGO` (2023 site) | up — faceted search, hierarchy browser, enrichment tool. HTTP only, no HTTPS |
| `https://supfam.org/SUPERFAMILY/dcGO/` (original) | up — bulk downloads still served by `cgi-bin/dcdownload.cgi` (incl. `?outer=PFAM`) |
| `github.com/hfang-bristol/dcGO` | exists, last pushed 2023-05-01 |

There is also a **third dcGO paper** we had not been comparing against:
`docs/EMS185259.pdf` — Bao *et al.*, *The dcGO Domain-Centric Ontology Database
in 2023* (JMB 435:168093). It supersedes the 2013 database paper and changes what
"the original" means for this comparison:

* **Domains.** 2023 adds Pfam and InterPro alongside SCOP superfamily/family —
  so an InterPro-keyed comparison is now possible without re-keying. But the
  paper reports only *~1,000 Pfam and ~800 InterPro domains* against our 19,449
  InterPro entries in human; **verify that figure against the actual download
  before citing it**, since it may be an increment rather than a total.
* **Ontologies.** The 2023 set both grew and shrank. Added: MONDO, EFO, KEGG,
  Reactome, PANTHER, WikiPathways, MitoCarta, DGIdb, Open Targets tractability
  buckets, ENRICHR/TRRUST transcription factors, MSigDB hallmarks. **Dropped**
  from the 2013 set: EC, UniPathway, UniProt keywords, DrugBank ATC codes.
* **Consequence for §3.** Compare against 2023 for GO/Reactome, but the 2013
  release is the only comparator for EC / UniPathway / keywords. Our
  chemistry-level layers (`ligand`, `cofactor`, `rhea`) and `subcellular`,
  `tcdb`, `merops`, `cazy`, `complex`, `condensate` have **no counterpart in any
  dcGO release**, so for those the InterPro2GO-style and temporal tests are the
  only available validation.

- [x] Download the original dcGO / SUPFAM domain–GO associations.
      (`scripts/download_data.py --group dcgo-reference`; SHA-256 pinned in
      `config/settings.py`.)
- [x] Re-key our domain parser on the `SSF` signature
      (`run_dcgo_human.py --domain-key ssf`).
- [x] Report agreement on the shared domain space and characterize the deltas.

**Acceptance:** a table quantifying overlap with the original dcGO on mappable
domains, with a written explanation of the deltas.

---

### 3.1 Result (2026-08-04) — done, with the confounds stated first

Driver: `validation/compare_original_dcgo.py`. Metrics:
`validation/dcgo_comparison_{metrics,variants,by_aspect,by_domain,supra}_{manual,iea}.tsv`.

#### Read the confounds before the numbers

They are large, they are quantified, and they push in one direction — **against
recall**. Anyone reading "recall 0.07" as a statement about method quality has
been misled by the experiment, not by the method.

| # | Confound | Quantified |
|---|---|---|
| 1 | **Vintage.** Their tables are stamped 3 Apr 2016; our GOA is 2026. | **1,795** of the 18,666 GO terms they scored are now obsolete in GO. **4,518** of our 14,517 terms appear nowhere in their table, so 792 of our 6,829 significant pairs (11.6%) are not comparable in either direction. |
| 2 | **Species scope.** They used 2,414 genomes + all UniProt (>80M sequences); we use human only. | **69.4%** of the pairs they call significant and we do not (30,485 / 43,946) have **zero co-occurring human proteins** — no threshold, method or statistic could recover them. This is the single largest driver of the recall figure. |
| 3 | **Evidence policy.** Their IEA/evidence-code policy is **not stated** in any paper available to us. We do not know it; we did not assume it. | Bracketed by running both ways. `--evidence-filter manual` → precision 0.537, recall 0.069; `--evidence-filter all` (IEA included) → precision 0.552, recall 0.082. The unknown is worth ~1.5 precision points, so it does not change any conclusion. |
| 4 | **Reachable domain space.** InterPro integrates SUPERFAMILY at superfamily level only; there is no `fa`. | Only **1,561 of their 4,355** domains (**35.8%**) are reachable at all. The 2,794 SCOP-family entries are structurally out of scope. |
| 5 | **Their table is not a full test matrix.** `GO_mapping` covers 1.9% of the shared sunid × GO space. | Absence from their table is reported as its own category (`extra_not_in_their_table`), never folded into "they disagree". |

Two further confounds surfaced while building the comparison. Neither is in the
plan, and a naive implementation would have got both wrong.

**(a) Their published set is not `all_score < 1e-3`.** At the `sf` level they
score 404,288 pairs, of which 134,665 sit below their 1e-3 threshold — but they
*ship* 108,612 direct + 85,683 inherited. Every direct row bar one is below 1e-3
(direct ⊂ significant); **26,054 significant pairs (19.3%) are not shipped at
all**, removed by the further filtering their method applies on top of the FDR
(the relative/parental-background test and most-specific-level selection); and
**none** of the inherited rows is itself below 1e-3 — they are pure true-path
propagation onto ancestors. The comparison therefore uses **what they
published**, split direct vs. inherited, as the primary comparator, and reports
the naive re-threshold as one variant among three.

**(b) 80,826 of their 404,288 `sf` rows (20%) carry `all_score` exactly 1.0**,
which is the column's SQL `DEFAULT`. Those are indistinguishable from "never
actually tested". The clearest illustration is that P-loop NTPase
(`SSF52540`) → ATP binding (`GO:0005524`) — about as canonical as
domain–function associations get, and one we call at p = 2e-31 — carries
`all_score = 1.0` on their side and is not shipped, while its parent
`GO:0016887` (ATPase activity) *is* shipped as inherited. This caps how much any
rank correlation against their FDR column can mean, and it is reported rather
than worked around.

#### The two keyings, run end to end

`uv run python run_dcgo_human.py [--domain-key ssf] --num-cores 8`, human, GO,
`--evidence-filter manual`, FDR < 0.01:

| | `--domain-key interpro` (default) | `--domain-key ssf` |
|---|---:|---:|
| Proteins in the universe | 18,908 | **12,316** (−34.9%) |
| Distinct single domains | 19,449 | **911** (−95%) |
| Domain features (single + supra) | 103,167 | 3,553 |
| GO terms | 16,389 | 14,517 |
| Fisher tests | 1,690,803,963 | 51,578,901 |
| Significant (FDR < 0.01) | 165,618 | 13,771 |
| BH raw-p cutoff | 9.79e-07 | 2.66e-06 |
| Wall clock (Fisher / BH / total) | 15.2 / 50.5 / **68.8** min | 0.5 / 1.3 / **1.8** min |

(Side observation, not part of §3: at 1.69 B tests the run is dominated by
`benjamini_hochberg_correction`, which walks a Python loop over every p-value —
50 of the 69 minutes. It is 3× the cost of the 1.7 M tests/s Fisher stage and is
the obvious next optimisation. It affects the default keying only.)

The trade is exactly the one the original made: **breadth collapses ~20×, depth
per domain rises** (median 4 proteins per superfamily vs 1 per InterPro entry).
SSF keying is the right lens for this comparison and the wrong default for the
pipeline. SSF↔InterPro is a 1:1 bijection over human data (911 ↔ 911, no
many-to-one either way), so every SSF-keyed row carries its InterPro entry in a
trailing `interpro_id` column for cross-referencing.

#### Agreement, single domains

Shared space: **904 / 911 (99.2%)** of our superfamilies are in their table
(the 7 absentees are small superfamilies they never associated with any GO term
— CENP-B dimerisation, HBS1-like, GTP cyclohydrolase I feedback regulatory
protein, …; listed in full in `dcgo_comparison_metrics_manual.tsv`);
9,999 shared GO terms.

Primary comparison — ours at FDR < 0.01 vs **their published direct** set:

| | value |
|---|---:|
| Our significant pairs (in shared space) | 6,037 |
| Their significant pairs (in shared space) | 47,189 |
| Shared | 3,243 |
| **Precision (ours also called by them)** | **0.537** |
| Recall (theirs also called by us) | 0.069 |
| Jaccard | 0.065 |

All six threshold × definition variants, none selected after the fact:

| our threshold | their definition | ours | theirs | shared | precision | recall |
|---|---|---:|---:|---:|---:|---:|
| FDR<0.01 | published, direct | 6,037 | 47,189 | 3,243 | 0.537 | 0.069 |
| FDR<0.001 | published, direct | 4,478 | 47,189 | 2,530 | 0.565 | 0.054 |
| FDR<0.01 | published, direct+inherited | 6,037 | 77,652 | 3,591 | 0.595 | 0.046 |
| FDR<0.001 | published, direct+inherited | 4,478 | 77,652 | 2,799 | 0.625 | 0.036 |
| FDR<0.01 | `all_score` < 1e-3 | 6,037 | 58,529 | 3,561 | 0.590 | 0.061 |
| FDR<0.001 | `all_score` < 1e-3 | 4,478 | 58,529 | 2,724 | 0.608 | 0.047 |

Precision sits in a narrow **0.54–0.63** band across every variant, and tightening
our threshold to theirs moves it by ~3 points. Per aspect (primary variant):
CC 0.623 > MF 0.544 > BP 0.496. Per domain, over the 190 superfamilies where we
make ≥10 calls, median precision is 0.549 (IQR 0.455–0.634); only 2 have zero
agreement and 7 exceed 0.8 — disagreement is spread, not concentrated in a few
pathological families.

#### Where the disagreement actually is

This is the part a set intersection cannot tell you, and it is why the
per-pair FDR in `Domain2GO.sql.gz` was worth parsing.

**Pairs they call and we do not (43,946)** — our Fisher p recomputed for each,
on the same universe the run used (recomputation verified against the pipeline's
own p-values to 4.9e-7 relative):

| bucket | n | share |
|---|---:|---:|
| Zero co-occurring human proteins — unreachable | 30,485 | 69.4% |
| Supported, p < 0.05 but above our BH cutoff (2.66e-6) | 9,096 | 20.7% |
| Supported, p ≥ 0.05 — we genuinely see nothing | 4,365 | 9.9% |

So **90% of our "misses" are explained by the species-scope confound plus
sub-threshold signal in the right direction** (median p for supported misses:
0.015). Only ~10% are pairs where human data is available and shows nothing.

**Pairs we call and they do not (2,794):**

| bucket | n |
|---|---:|
| Never scored by them (absent from `GO_mapping`) | 735 |
| They published it, but as a true-path-**inherited** annotation | 348 |
| They scored it below 1e-3 but did not ship it | 318 |
| They scored it in 1e-3 … 0.05 ("nearly called it too") | 586 |
| They scored it > 0.5 | 569 |

#### Calibration (threshold-independent)

Spearman of our −log10(p) against their −log10(FDR), our p recomputed for every
pair they scored in the shared space:

| population | n | ρ |
|---|---:|---:|
| All pairs they scored | 171,006 | 0.127 |
| … with ≥1 co-occurring human protein (drops the a=0 mass confound 2 makes invisible to us) | 48,034 | **0.467** |
| … and excluding their `all_score = 1.0` column default (confound b) | 41,065 | 0.398 |

So on pairs where both sides genuinely have a number, our ranking and theirs
correlate at ρ ≈ 0.4–0.47: clearly related, far from interchangeable. Against
their `all_hscore_max`: ρ = 0.019 all / 0.224 supported — but that column is a
max over the GO subtree on published rows only, not the same quantity as a
p-value, so little should be read into it.

#### Supra-domains — the novel part

Nothing else in this plan gives our supra-domain machinery an external
reference. `SP2GO.txt` does.

| | value |
|---|---:|
| Supra-domain architectures we observe in human | 2,642 |
| … with ≥1 significant association | 1,107 |
| Their supra-domain architectures | 7,163 |
| Shared, exact N→C order | 1,335 (precision 0.505) |
| Shared, order-insensitive | 1,267 (precision 0.539) |
| Architectures compared (both sides, shared GO terms) | 709 |
| Our associations on those / theirs (direct) | 4,320 / 41,574 |
| Shared associations | 2,226 |
| **Supra precision (ours in theirs, direct)** | **0.515** |
| Supra recall | 0.054 |
| Supra precision vs their direct + inherited | 0.553 |

Order convention on their side is undocumented, so both the exact-order and
order-insensitive variants are reported; the difference is small (0.505 vs
0.539), which suggests our N→C convention broadly matches theirs. **Our
supra-domain associations agree with the published ones at essentially the same
rate as our single-domain ones (0.515 vs 0.537)** — the supra-domain machinery
is not producing a different quality of output from the single-domain path.

#### Verdict: §3 is closed for GO, on the `sf` half of their domain space

What is now established: the domain universes join cleanly (99.2% of our
superfamilies), our associations agree with the published ones at ~54–63%
precision on every reasonable definition of "their significant set", the
disagreement is diagnosed rather than merely counted (90% of misses are
species-scope or sub-threshold), and the supra-domain machinery is externally
corroborated for the first time.

What is **not** established, and should not be claimed:

* **Nothing about recall.** The comparison cannot measure it — a human-only
  universe cannot reach 69% of their calls at any threshold.
* **Nothing about the SCOP-family (`fa`) half** of their release, or about the
  Pfam-keyed release. Reachable only by adding a `pfam` domain key (the parser
  seam now supports it: add one entry to `SIGNATURE_PREFIXES`).
* **Nothing about the 2023 dcGO release.** It publishes no bulk download; only
  the 2013 tables are comparable at scale.
* **Nothing about the non-GO ontologies.** `Domain2EC.txt`, `Domain2KW.txt`,
  `Domain2UP.txt` exist and are the only comparator for our EC / keyword /
  UniPathway layers; they were not used here.
* The surprise score's **novelty discount is InterPro2GO-based** and reports
  `no-reference` for every SSF-keyed candidate. `scripts/rank_surprising_associations.py
  --domain-key ssf` runs, but its novelty factor is inert.

- [ ] **Open:** extend to `--domain-key pfam` and to `Domain2EC` / `Domain2KW` /
      `Domain2UP` for the non-GO layers.

---

## 4. Current component ablation *(2026-09-01; current main)*

This section reports only the ablation run made from commit `3f37558` after the
hierarchy and relative-inference fixes in #61–#67. Earlier ablation measurements
have been removed because they describe an implementation that no longer
exists.

### Design

Training uses archived human GOA release 205; evaluation uses the current human
GOA with experimental evidence only and the CAFA no-knowledge cohort. All
inference configurations were rerun from scratch. The experiment is a factorial
over the three independent hierarchy stages on the supra-domain model:

| rung | input annotations propagated | relative inference | inferred associations propagated |
|---|---:|---:|---:|
| `supra` | no | no | no |
| `supra_input` | yes | no | no |
| `supra_relative` | no | yes | no |
| `supra_output` | no | no | yes |
| `supra_input_relative` | yes | yes | no |
| `supra_input_output` | yes | no | yes |
| `supra_relative_output` | no | yes | yes |
| `full` | yes | yes | yes |

A separate `single` rung measures the effect of adding supra-domains. Every
hierarchy combination is a complete pipeline run; no post-hoc legacy rung is
constructed. Scores are `-log10(q)` for every primary comparison. Results cover
BP, MF and CC at IC floors 0, 2 and 4, with 1,000 paired protein bootstraps.
The supra rung also uses 200 seeded domain-label permutations.

The scored cohorts contain 324 BP, 418 MF and 572 CC proteins before IC
filtering. Run manifests are in `validation/ablation_manifests/`; machine-readable
results are:

- `validation/ablation_metrics.tsv`
- `validation/ablation_paired_bootstrap.tsv`
- `validation/ablation_permutation_null.tsv`
- `validation/ablation_selection_counts.tsv`
- `validation/ablation_provenance.tsv`

### Results

The domain signal is real: the supra model beats all 200 shuffled mappings for
both F_max and AUPRC in all nine aspect × IC cells (18/18 empirical
`p = 1/201`).

Supra-domains themselves have little protein-centric effect. Relative to single
domains, their F_max change is not significant in any of nine cells; AUPRC
improves significantly in one cell.

The three hierarchy stages do not behave as one component:

- **Input propagation is context-dependent but important.** Added directly to
  supra it changes F_max significantly in five of nine cells (six positive,
  three negative). Added to relative inference plus output propagation, it
  improves both F_max and AUPRC in all nine cells, significantly in eight.
- **Relative inference is the main negative component on this benchmark.**
  Added directly to supra, it lowers F_max and AUPRC in all nine cells,
  significantly in eight. Added to input plus output propagation, it lowers
  F_max in seven of nine cells (six significant) and AUPRC in eight of nine
  (seven significant).
- **Output True Path propagation is not uniformly harmful.** Added directly to
  supra, F_max rises in seven of nine cells and is significant in one; AUPRC is
  mixed. Added after input propagation and relative inference, it raises F_max
  and AUPRC in six of nine cells, significantly in three and four respectively.

There is no single winning configuration across all cells. The full
paper-parity pipeline is best only for MF at IC 0. Input propagation without
relative inference is best in several informative-term cells; input plus output
propagation is best for BP at IC 0/2 and CC at IC 0. Single domains remain best
for CC F_max at IC 2/4.

### Decision

Do not describe “True Path” as beneficial or harmful as a unit. Input
propagation, relative inference and output propagation are separate operations
with different effects. On this temporal benchmark, relative inference requires
the strongest justification and should not be enabled by default merely for
paper parity. The primary configuration must be pre-specified against an
untouched final evaluation rather than selected from these cells.

## 5. Decisions to settle before writing the paper

### Held-out validation of the surprise score (2026-07-28) — DONE

`validation/temporal_surprise.py` re-uses the §2 split (t0 = GOA 205, 2021-04;
t1 = 2026-06) to ask whether the emergent-combination ranking predicts *future
curation*. For each t0 association, the proteins carrying the combination that
lacked the term at t0 are its predictions; hits are those annotated by t1; the
control is the term's own acquisition rate over all domain-carrying proteins
that lacked it. Metrics in `validation/temporal_surprise_metrics.tsv`, outcomes
per association in `validation/temporal_surprise_associations.tsv`.

- [x] **The supra-domain associations predict future curation: 12.5×
      enrichment** over the terms' own acquisition rates (2,181 hits on 170,416
      predictions vs ~175 expected), 95% percentile CI [10.87, 14.39]. These
      strata are fixed sets, so this is an ordinary bootstrap and the interval is
      sound. This is the out-of-sample evidence the domain-combination claim was
      missing.
- [x] **The surprise ranking is *not* demonstrably better than ranking by the
      dcGO q-value — and at these budgets the comparison is not resolvable at
      all.** The score assigns exactly 0.000 to 9,923 of the 10,136 evaluated
      associations (`q_emergence` saturates at 1), so a slice sized by prediction
      budget is mostly an arbitrary tie-break: 63% of it at 10,000 predictions,
      85% at 40,000. Re-breaking the tie at random moves the observed +7.96 at
      the 10,000 budget to +0.12 [−5.53, +7.57] — the observed value is the 98th
      percentile of its own tie-break spread. At 2,000 predictions only 17 dcGO
      associations fit, and the percentile and basic intervals disagree about the
      *sign*. Report it as "no demonstrated ranking advantage"; the score's
      contribution is interpretability (redundant-signature and curated-novelty
      filtering) and a bias toward rarer, higher-IC terms, not ranking power.
      Re-running the head-to-head needs a graded emergence score first.
- [ ] **Open: the emergence/testability tension.** The most emergent
      associations leave the fewest standing predictions (emergence requires
      that carriers are already nearly all annotated), so the sharp end of the
      ranking cannot be validated this way — `surprise top-25` yields 117
      predictions against an expected 0.14 hits. A revised score should weigh
      emergence against how many predictions it leaves outstanding.
- [ ] **Open: same look-ahead caveat as §2.** Domain architectures come from the
      current `protein2ipr`, so this is annotation-temporal, not prospective.


- [ ] **True Path Rule default.** The original dcGO makes TPR central; here it is
      opt-in. Decide: make it default (recommended, with `--disable-true-path`
      escape) or justify keeping it optional. Whatever is chosen, benchmark both.
- [ ] **Species scope.** Human-only reduces statistical power. State this as a
      deliberate scope choice, or extend to a multi-species run for the paper.
- [ ] **Significance vs. effect size.** The human run does ~1.7 B domain×GO tests
      and keeps 165,823 at FDR<0.01 — many on thin evidence (low `n_observations`).
      Decide and document a minimum-evidence / effect-size floor
      (`MIN_PROTEINS_PER_ASSOCIATION`, odds-ratio bounds). Note `odds_ratio`
      prints `0.0000` when `d=0` and `inf` when `b*c=0` — cosmetic artifacts of
      `a*d/(b*c)`, not the ranking signal (`hyper_score`/FDR are); consider a
      Haldane correction so the reported odds ratio is interpretable.

---

## 6. Reproducibility (must-haves for submission)

- [x] **Record exact run provenance**: every run writes
      `run_manifest_<ontology>.json` with input/output SHA-256 hashes, embedded
      release headers, source URLs, Git state, the `uv.lock` hash, the command
      line, every effective parameter and threshold, and summary counts. See
      [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
- [ ] **Pin dated dataset versions**: replace the mutable `current_release`
      inputs with archived GOA / InterPro / GO releases for every reported
      number.
- [ ] **Extend provenance to the downstream tools**: the surprise-score driver
      and the `validation/` benchmarks still record nothing.
- [ ] **One-command reproduction** of each table/figure from raw downloads.
- [ ] **Archive** the exact input snapshots (Zenodo/figshare) since
      `current_release` URLs move.

---

## Suggested order of work

1. ~~§0 pipeline correctness~~ — ✅ done (#15, #17).
2. ~~§1 (fix the existing comparison)~~ — ✅ done (#14, #15).
3. ~~§2 (temporal benchmark + baselines) — the core result~~ — ✅ done (#8).
4. **Method-vs-paper audit + metric hardening — ← next** (not ablation yet):
   verify against Fang & Gough 2013 (validation approach, optimal-level test),
   and settle the informative-term / domain-centric evaluation.
5. §4 (ablation), §3 (original-dcGO comparison), §5–§6 after the yardstick is
   trusted.

§2 gives a real precision-capable result: on informative terms dcGO beats the
random null by 4–26× and beats the naive baseline in every aspect. Naive's raw
F_max lead was base-rate recovery of uninformative terms (e.g. `protein binding`).
