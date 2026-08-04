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

**Where §2 landed:** with the calibrated **p-score** predictor (default) and a
leak-free no-knowledge gate, dcGO **beats the CAFA naive baseline on F_max at face
value in BP (0.248 vs 0.115) and CC (0.380 vs 0.343)**, and beats it in **every
aspect on informative terms** (IC-filtered), while staying **1.3–25× above the
random-domain null**. MF is the one aspect naive leads at IC≥0 (0.360 vs 0.464),
because its truth is dominated by `protein binding` (84.6% of experimental MF
annotations, near-zero IC) — and dcGO overtakes it decisively the moment that
noise is excluded (IC≥2: 0.365 vs 0.053). This confirms the concern that a
temporal CAFA benchmark rewards recovery of the popularity-weighted curation
frontier; on informative function, dcGO is clearly ahead.

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
2. **Fold the method into the pipeline (paper-parity, end-to-end).** The relative
   inference + p-score are applied post-hoc here. Wire the relative
   (parental-background) test into `run_dcgo_human.py` inference itself (combine
   overall/relative then FDR<1e-3, per the paper — currently a post-hoc
   `alpha<0.05` filter), and expose the p-score predictor as the standard
   protein-prediction path. Makes "dcGO-2.0 == dcGO Predictor" defensible.
3. **§4 ablation** — now that the yardstick is trusted, run the temporal benchmark
   per pipeline config (single-domain / +supra / +shrinkage / +TPR).
4. **§3 original-dcGO comparison** — re-key the domain parser on the `SSF`
   (SUPERFAMILY/SCOP) or `PF` (Pfam) signature from `protein2ipr` col 4 to compare
   on the original's domain universe. See [[original-dcgo-methodology]].
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
      GO:0005515 `protein binding` (**84.6%** of human experimental MF
      annotations, near-zero IC).

### Results — 2021→2026 temporal split (2026-07-09, p-score transfer)

No-knowledge benchmark sizes (IC≥0): **BP 324 / MF 418 / CC 572** proteins (a
leak-free gate on training evidence — much smaller and cleaner than the earlier
experimental-only gate, which wrongly admitted proteins whose t0 computational
labels the model had already seen).

**dcGO (p-score) beats the naive baseline at face value on BP and CC, and beats
it in every aspect once uninformative terms are excluded**, staying well above the
random-domain null throughout:

| Aspect | IC floor | dcGO F_max | naive F_max | random F_max | dcGO / random |
|--------|:--------:|-----------:|------------:|-------------:|--------------:|
| BP | ≥0 | **0.248** | 0.115 | 0.158 | 1.6× |
| BP | ≥2 | **0.170** | 0.071 | 0.053 | 3.2× |
| BP | ≥4 | **0.115** | 0.031 | 0.019 | 6.1× |
| BP | ≥6 | **0.077** | 0.010 | 0.003 | 24× |
| MF | ≥0 | 0.360 | **0.464** | 0.262 | 1.4× |
| MF | ≥2 | **0.365** | 0.053 | 0.088 | 4.2× |
| MF | ≥4 | **0.337** | 0.045 | 0.072 | 4.7× |
| MF | ≥6 | **0.217** | 0.018 | 0.009 | 25× |
| CC | ≥0 | **0.380** | 0.343 | 0.291 | 1.3× |
| CC | ≥2 | **0.239** | 0.153 | 0.072 | 3.3× |
| CC | ≥4 | **0.134** | 0.099 | 0.031 | 4.3× |
| CC | ≥6 | **0.124** | 0.044 | 0.015 | 8.1× |

(IC in bits; IC≥2 ⇒ term in ≤25% of proteins, IC≥6 ⇒ ≤1.6%. Full table with
S_min/AUPRC in `validation/temporal_benchmark_metrics.tsv`.)

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

**Acceptance: met.** dcGO clears both mandatory baselines (random null and naive)
across all three aspects on informative terms, and beats naive at face value on 2
of 3 aspects.

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

### Breadth: does the predictive signal hold beyond GO? (2026-07-28) — DONE

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

| Ontology | assoc. | predictions | hits | hit rate | expected | enrichment (95% CI) |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| **go** (anchor) | 92,377 | 3,139,544 | 106,816 | 3.40% | 0.30% | **11.3× [10.9, 11.8]** |
| reactome | 54,685 | 1,487,359 | 2,056 | 0.14% | 0.02% | **8.0× [6.7, 9.6]** |
| cofactor | 363 | 1,744 | 221 | 12.67% | 3.97% | **3.2× [1.7, 4.0]** |
| subcellular | 13,373 | 507,262 | 7,005 | 1.38% | 0.49% | **2.9× [2.7, 3.0]** |
| keyword | 85,323 | 2,619,504 | 35,350 | 1.35% | 0.79% | **1.7× [1.6, 1.8]** |
| complex | 3,606 | 31,076 | 24 | 0.08% | ~0% | *degenerate, see below* |
| disease | 44 | 369 | 0 | 0% | ~0% | *undefined* |
| ligand | — | — | — | — | — | *not testable, see below* |

(In-scope universe where it differs from all-domains; keyword and GO annotate
essentially every protein, so the two coincide.)

**The signal generalises.** Every ontology with enough data enriches above 1
with the interval excluding it — dcGO's association-finding is not a GO artefact.
GO remains strongest and Reactome follows; the ordering tracks how structured and
curated the vocabulary is. Keywords come last because 720 near-universal terms
sit at a 0.79% base rate, leaving little headroom.

**Three results that must not be quoted at face value:**

- **`complex` is a degenerate ratio, not a 265× triumph.** ComplexPortal averages
  ~1.5 proteins per complex, so the base rate rounds to zero and any hit at all
  divides by nearly nothing. 24 hits. The interval excludes 1, so signal exists;
  the *magnitude* is uninterpretable. Report as "detectable, magnitude
  meaningless".
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
  much that matters: GO scores **11.3×** under the loose protocol against
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
      predictions vs ~175 expected), bootstrap CI [10.9, 14.4]. This is the
      out-of-sample evidence the domain-combination claim was missing.
- [x] **The surprise ranking is *not* demonstrably better than ranking by the
      dcGO q-value.** At matched prediction budgets the point estimate favours
      surprise 3/3 (15.5 vs 5.3, 21.2 vs 13.2, 11.6 vs 10.8) but a **paired**
      bootstrap — re-ranking both ways inside each resample, because the two
      rankings share a candidate pool and their independent intervals are
      therefore correlated — puts zero inside every interval. Report it as "no
      demonstrated ranking advantage"; the score's contribution is
      interpretability (redundant-signature and curated-novelty filtering) and a
      bias toward rarer, higher-IC terms, not ranking power.
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
