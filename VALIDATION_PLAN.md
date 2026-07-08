# Validation & Benchmarking Plan

This document defines what must be done to turn dcGO-2.0 from a working codebase
into a defensible methods contribution. It is the "(b)" companion to the
engineering cleanup: the code runs; this is about showing the results are
*reliable*.

## Status snapshot (2026-07-08)

| Item | State |
|------|-------|
| §0 Pipeline correctness | ✅ done — 4 correctness bugs found & fixed with tests |
| §1 Reframe InterPro2GO comparison | ✅ done (#14, #15) — ~65% coverage at FDR<0.01 |
| **§2 Temporal held-out benchmark (CAFA-style)** | ⏳ **next — the core result** (#8) |
| §3 Compare to original dcGO | ⬜ open (#9) |
| §4 Supra-domain ablation | ⬜ open (#10) |
| §5 Pre-paper method decisions | ⬜ open (#11) |
| §6 Reproducibility | ⬜ open (#12) |

**Next action:** implement §2 (see the implementation sketch there). It reuses
the propagation harness now on `main` and is what converts today's *coverage*
story into a real *precision* story.

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

## 2. Temporal held-out benchmark (CAFA-style)  *(core evidence)*

This is the single most important addition. It measures whether dcGO predicts
annotations that were **later** confirmed — a real, precision-capable test.

- [ ] Obtain two dated GOA snapshots, e.g. a **training** release (`t0`) and a
      **test** release (`t1`, ≥1 year later). EBI keeps dated GOA archives.
- [ ] Train domain→GO associations using only annotations available at `t0`.
- [ ] Define the evaluation set as annotations that appear in `t1` but not `t0`
      (newly curated knowledge), restricted to experimental evidence codes.
- [ ] Score predictions with the **CAFA protein-centric metric**:
      propagate predicted GO terms to each protein via its domains, then compute
      **F_max** over a score threshold sweep, plus **S_min** (information-content
      weighted) and **AUPRC**. Use term information content from `t0`.
- [ ] Report separately for the three GO aspects (BP / MF / CC).

**Acceptance:** F_max clearly above the two mandatory baselines below, with a
plotted precision–recall curve per aspect.

### Baselines (required for any claim of value)
- [ ] **Naive baseline**: predict each GO term with probability = its frequency
      in the training set (CAFA's standard `Naive`).
- [ ] **BLAST/annotation-transfer baseline** *(optional but strong)*: transfer
      GO terms from the most similar annotated protein.
- [ ] **Random-domain baseline**: shuffle domain→GO labels and re-run the whole
      pipeline to get an empirical null for the association scores (confirms the
      FDR is calibrated, not just nominal).

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
3. **§2 (temporal benchmark + baselines) — the core result. ← next.**
4. §4 (ablation) — reuses the §2 harness.
5. §3 (original-dcGO comparison) and §5–§6 in parallel.

The engineering is now correct; §2 is what turns "it runs" into "it works."
