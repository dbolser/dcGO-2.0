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
| §3 Compare to original dcGO | 🟡 method audit done (see §2→"Method audit"); domain re-keying open (#9) |
| §4 Supra-domain ablation | ⬜ open (#10) — **now unblocked** (reuses the §2 harness) |
| §5 Pre-paper method decisions | ⬜ open (#11) |
| §6 Reproducibility | ⬜ open (#12) |

**Where §2 landed:** on **informative** GO terms dcGO clearly beats both required
baselines — the random-domain null (by **4–26×** on F_max once low-IC terms are
excluded) and the CAFA naive frequency baseline (in every aspect at every
informative IC floor). At face value (no IC filter) naive's raw F_max looks
higher, but that lead is pure base rate: it evaporates the moment near-universal
terms like `protein binding` (84.6% of experimental MF annotations) are removed,
while dcGO holds up. This is the real *precision* result §2 was for, and it
confirms the concern that a temporal CAFA benchmark rewards recovery of the
popularity-weighted curation frontier.

**Next action:** *not* ablation. The method audit vs Fang & Gough 2013 is **done**
(see §2 → "Method audit" and §3) and points to two concrete, paper-grounded
method changes: (a) add the **relative (parental-background) inference** — the
specificity test we omitted, which demotes generic low-IC associations at
inference time; and (b) adopt their **per-target p-score** (sum of h-scores,
min-max normalised per protein). Metric side: keep **both** protein-centric and
domain-centric evaluations, and report BP/MF as the headline (as the original did;
CC is least domain-relevant).

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
      (ran the standard pipeline on the 2021 GOA → 69,493 significant
      associations at FDR<0.01).
- [x] Define the benchmark as **no-knowledge proteins per aspect**: proteins with
      *no experimental annotation in that aspect at t0* that gained experimental
      annotation by t1. Score against their **full** propagated t1 experimental
      terms. (An earlier delta-only truth — `t1 minus t0` — was wrong: it scored
      correct predictions of already-known terms as misinformation. Fixed before
      any number was quoted — the §0 lesson.)
- [x] Score with the **CAFA protein-centric metric**: transfer predicted GO terms
      to each protein via its domains (union, max score, propagated), then
      **F_max** over a threshold sweep, plus **S_min** (marginal-IC weighted, IC
      from t0) and **AUPRC**. Rank by −log10(p) — `hyper_score` saturates (37% at
      exactly 100) and collapses the sweep.
- [x] Report separately for BP / MF / CC.
- [x] Sweep an **information-content floor** (`--min-ic`) that excludes
      near-universal, low-IC terms from truth and all methods alike — the fair,
      principled way to stop rewarding base-rate recovery of terms like
      GO:0005515 `protein binding` (**84.6%** of human experimental MF
      annotations, near-zero IC).

### Results — 2021→2026 temporal split (2026-07-09)

No-knowledge benchmark sizes (IC≥0): **BP 1,537 / MF 1,124 / CC 2,305** proteins.

**The headline is the information-content sweep.** At face value (IC≥0) dcGO
trails the naive frequency baseline on F_max. But that lead is *entirely* base
rate: the moment low-information terms are excluded, naive collapses toward the
random null while dcGO holds up — and dcGO **beats naive in every aspect at every
informative IC floor**.

| Aspect | IC floor | dcGO F_max | naive F_max | random F_max | dcGO / random |
|--------|:--------:|-----------:|------------:|-------------:|--------------:|
| BP | ≥0 | 0.276 | **0.289** | 0.144 | 1.9× |
| BP | ≥2 | **0.215** | 0.131 | 0.051 | 4.2× |
| BP | ≥4 | **0.175** | 0.042 | 0.023 | 7.6× |
| BP | ≥6 | **0.142** | 0.012 | 0.009 | 15× |
| MF | ≥0 | 0.319 | **0.579** | 0.098 | 3.2× |
| MF | ≥2 | **0.385** | 0.082 | 0.055 | 7.0× |
| MF | ≥4 | **0.373** | 0.081 | 0.040 | 9.4× |
| MF | ≥6 | **0.337** | 0.029 | 0.013 | 26× |
| CC | ≥0 | 0.395 | **0.520** | 0.234 | 1.7× |
| CC | ≥2 | **0.278** | 0.208 | 0.077 | 3.6× |
| CC | ≥4 | **0.203** | 0.103 | 0.038 | 5.4× |
| CC | ≥6 | **0.192** | 0.047 | 0.029 | 6.6× |

(IC in bits; IC≥2 ⇒ term in ≤25% of proteins, IC≥6 ⇒ ≤1.6%. Full table with
S_min/AUPRC in `validation/temporal_benchmark_metrics.tsv`.)

**Read-out:**
- **naive is a base-rate mirage.** Its F_max advantage vanishes with one filter:
  BP 0.289→0.012, MF 0.579→0.029, CC 0.520→0.047 as IC rises 0→6 — it converges to
  the random null because all it ever had were high-frequency generic terms.
  In MF, ~46% of benchmark proteins (1,124→608 at IC≥2) had *only* low-IC
  newly-curated terms — i.e. `protein binding` was their entire MF "truth".
- **dcGO degrades gracefully and dominates on informative terms.** It stays
  **4–26× above the random-domain null** once uninformative terms are removed,
  and beats naive by up to ~12× (BP) / ~11× (MF) / ~4× (CC). On MF it even
  *improves* at IC≥2 (0.319→0.385) — removing the protein-binding noise clarifies
  the signal. AUPRC follows the same flip (dcGO ≥ naive at every IC≥2).
- This is exactly the concern raised by Julian: a temporal CAFA benchmark rewards
  recovery of the **attention-biased, popularity-weighted** curation frontier, so
  raw F_max flatters a frequency baseline. Restricting to *informative* terms —
  what a domain→function method is actually for — shows dcGO clearly ahead.

**Acceptance: met.** dcGO clears both mandatory baselines (random null and naive)
on informative terms across all three aspects. The raw-F_max caveat is now
understood and controlled, not a mystery.

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

**The load-bearing insight:** their *relative (parental-background) inference*
does at **inference time** what our IC filter does at **evaluation time** — demote
associations that aren't stronger than their generic parent. We omitted it, so
our raw predictions are more promiscuous toward low-IC terms, which is exactly
what let naive look competitive at IC≥0. Adding the relative inference + their
per-target p-score is the paper-grounded way to close the gap, and is almost
certainly worth more than any ablation.

### Baselines
- [x] **Naive baseline**: predict each GO term at its (propagated) t0 frequency
      (CAFA's standard `Naive`).
- [x] **Random-domain baseline**: permute the domain→GO association labels
      (seeded) and re-transfer — an empirical null for the transfer step. dcGO
      sits far above it, confirming the associations are non-random.
- [ ] **BLAST/annotation-transfer baseline** *(optional, still open)*: transfer
      GO terms from the most similar annotated protein.
- [ ] **Full-pipeline shuffle** *(stronger null, still open)*: shuffle labels and
      re-run Fisher end-to-end to confirm the FDR itself is calibrated (the
      current shuffle is at the association level, not the whole pipeline).

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

- [ ] Download the original dcGO / SUPFAM domain–GO associations.
- [ ] **Preferred:** re-key our domain parser on the `SSF` (or `PF`) signature
      instead of the integrated `IPR` entry — a near-apples-to-apples
      reproduction of the original's domain definitions (roughly a one-column
      change in `domain_annotation_parser`). Alternatively map `IPR ↔ SSF/PF` via
      `protein2ipr` and document coverage.
- [ ] Report agreement on the shared domain space and characterize where
      dcGO-2.0 differs (newer GOA, SCOP-vs-InterPro granularity, FDR threshold,
      supra-domains).

**Acceptance:** a table quantifying overlap with the original dcGO on mappable
domains, with a written explanation of the deltas.

---

## 4. Ablation study  *(isolates the contribution)*

The supra-domain + shrinkage machinery is the main methodological novelty. It
must be shown to help, not just to exist. Run the temporal benchmark (§2) for
each configuration:

| Config | Flags |
|--------|-------|
| Single domains only | `--disable-supra-domains` |
| + supra-domains | (default) |
| + supra + shrinkage | `--enable-shrinkage` |
| + True Path Rule | `--enable-true-path` |
| Full | `--enable-shrinkage --enable-true-path` |

- [ ] Report F_max / AUPRC per configuration and per GO aspect.
- [ ] Quantify **how many supra-domains produce associations not obtainable
      from their constituents**, and whether low-count supra-domains are
      correctly regularized by shrinkage (they should not dominate the top
      predictions).

**Acceptance:** each enabled stage shows a measurable, explained effect (a
stage that doesn't help is a finding too — report it).

---

## 5. Decisions to settle before writing the paper

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

- [ ] **Pin dataset versions**: record the GOA / InterPro / GO release dates and
      URLs used for every reported number (extend `scripts/download_data.py` to
      log resolved versions).
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
