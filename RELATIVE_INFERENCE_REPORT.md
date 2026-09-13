# The relative inference in dcGO-2.0: what we implemented, what we expected, what we measured

*Pre-publication report, 2026-09-12. Prepared for discussion with the original
dcGO authors before any of this is written up for submission.*

---

## Summary

We have reimplemented dcGO's two-step statistical inference (Fang & Gough
2013) from the published descriptions, on current data (UniProtKB/GOA and
InterPro, human-centred). The overall inference reproduces the published
behaviour well. The **relative inference** — the parental-background test that
is Step 2 of the method — does not behave the way the papers led us to
expect, and we would like to compare notes before concluding anything.

Three observations drive this report:

1. **It removes most of the signal.** With everything else held fixed,
   enabling the relative test cuts the significant association set to ~19% of
   its size, in both configurations we ran it in (163,153 → 30,700 and
   486,041 → 91,190 associations).
2. **On a held-out temporal benchmark it makes protein-level prediction
   worse, and not simply by trading sensitivity (recall) for PPV
   (precision).** Training on 2021 GOA and evaluating on annotations gained
   by 2026 (CAFA-style, no-knowledge cohort), adding the relative test to
   the un-propagated model lowers both F_max and AUPRC (area under the
   PPV–sensitivity curve) in 9 of 9 aspect × IC cells (8 of 9 significant,
   paired protein bootstrap); in the configuration closest to the paper's
   (input propagated) the effect is smaller but still mostly negative.
   Decomposing the operating point: sensitivity falls in 26 of 27 contrast
   cells, while PPV rises only in some configurations — and in the
   biological-process aspect the test loses on both axes almost everywhere.
   True-Path propagation of the *input* annotations recovers part of the
   loss but not all of it.
3. **The published description is ambiguous at the exact point our
   implementation had to make a choice.** The Methods text of the 2013
   BMC Bioinformatics paper defines the relative background as proteins
   annotated to **all** direct parents; Figure 1's caption of the same paper
   says **any** direct parents; the worked example's numbers are only
   consistent with **any** (union). We implemented the union. If the original
   implementation did something else — the intersection, a per-parent test, a
   conditional rule — that could explain much of the difference.

Section 6 lists the specific questions. Everything quantitative in this
report is reproducible from the repository; Section 8 gives provenance.

---

## 1. What the papers say

Sources: Fang & Gough, *BMC Bioinformatics* 2013, 14(S3):S9 (the methods
paper, "S9" below); Fang & Gough, *NAR* 2013, 41:D536 (the database paper);
Bao et al., *J Mol Biol* 2023, 435:168093 (the 2023 update).

**The two tests.** S9, Methods ("Two statistical inferences"): an overall
p-value is computed "using all analyzable UniProt proteins (i.e., those
annotated to the root of GO term after applying the true-path rule) as the
background", and a relative p-value "using the background of only those
UniProt proteins annotated to **all** direct parental GO terms."

**The combination.** S9, Methods: "We first took the larger one of the
overall and relative p-values … the Benjamini-Hochberg derived FDR … was used
to determine the statistical significance … A stringent threshold of FDR
(< 10⁻³) was accepted … we also took the smaller of the overall and relative
hypergeometric scores … denoted as h-score."

**The ambiguity.** Within the same S9 paper:

| Place | Wording | Reading |
|---|---|---|
| Methods text, p.9 | "annotated to **all** direct parental GO terms" | intersection |
| Figure 1, flowchart box | "annotated to **all** direct parental GO terms" | intersection |
| Figure 1, caption | "N_pa is the total number of Uniprots that can be annotated by **any** direct parental GO terms" | union |
| NAR paper, D537 | "annotated to direct parents of the term" | either |

The Figure 1 worked example (supra-domain 82199,57667 × GO:0019827 "stem
cell maintenance") prints N_pa = 21,240 against N = 119,685. The three direct
parents include "developmental process", a near-root term; a background of
17.7% of the universe is what the **union** gives. The intersection — proteins
annotated to "developmental process" AND "negative regulation of cell
differentiation" AND "stem cell development" simultaneously — would be far
smaller: it is capped by the smallest of the three parents ("stem cell
development", which cannot plausibly reach 17.7% of the universe) and
bounded below only by the child's own K = 362, so it could be up to ~60-fold
smaller than the printed 21,240. So the example arithmetic supports *any*;
the prose supports *all*.

**The stated purpose.** S9: the dual inferences "aimed to ensure that only
the most relevant GO terms could be retained." The NAR paper is more
operational: "If a GO term and its parent term are both significantly
associated … and if the term is not significantly different from the parent
term using the second background, then it is desirable to only associate the
parent term." Neither paper quantifies the effect — no before/after counts,
no ablation of the relative step appears in any of the three.

**What the papers do not specify** (each of these forced a choice on us):
the behaviour at terms with no parents (the ontology roots); any minimum
background size; whether the relative test runs for every pair or only
conditionally (the NAR prose can be read as a conditional keep-the-parent
rule rather than an always-on second test); the BH hypothesis family
(pooled or split by domain level); the exact evidence-code list behind
"experimental or manual"; and input propagation, stated only via the
background definition ("after applying the true-path rule") rather than as
an explicit step — the shared K = 362 in both Figure 1 tests confirms it.

We also note the 2023 update's method summary describes Fisher's exact test,
FDR and output propagation but does not mention the relative inference at
all — we cannot tell whether it was dropped from later builds or merely
elided in a simplified summary.

## 2. What we implemented

The full recipe, as it runs (`--enable-relative-inference`; code in
`src/relative_inference.py`, wiring in `run_dcgo_human.py`, every parameter
recorded in machine-readable run manifests):

For every (domain D, term T) pair that co-occurs in at least one protein:

1. **Overall test** — one-sided Fisher (enrichment), background = all
   analysable proteins (≥1 domain and ≥1 annotation).
2. **Relative test** — let B = the **union** of the true-path-propagated
   protein sets of T's *direct* parents (GO: `is_a`/`part_of` edges only).
   - no parents → `relative_p = 0` (the pair passes on the overall test alone);
   - |B| < 3 → `relative_p = 1` (counted, conservatively unpassable);
   - otherwise: one-sided Fisher on the 2×2 restricted to B, all counts
     propagated.
3. **Combine** — `p = max(overall_p, relative_p)` *before* correction (a
   valid intersection-union p-value; Berger 1982), then Benjamini-Hochberg at
   FDR < 0.01, single-domain and supra-domain hypotheses corrected as
   separate families against their full dense sizes.
4. **h-score** = min(overall, relative) hypergeometric score, as published.
5. Optional reporting floor `--min-ic` (below).

Independent switches, matching the paper's three steps: `--propagate-annotations`
(true-path propagation of the input protein→term map, before any testing) and
`--enable-true-path` (propagation of the significant output up the DAG,
after selection). The relative test's own backgrounds are internally
propagated regardless — otherwise a parent nobody is directly annotated to
would have an empty background.

**Choices the papers did not force** (candidate divergence points, in order
of expected impact):

- **Union background** ("any"), per Figure 1's caption and example — not the
  intersection ("all") of the Methods text. An earlier version of our code
  ran one test per parent and took the worst p; moving to the union raised
  human single-domain associations (with input propagation) from 18,612 to
  33,185 (56.1% of GO terms have more than one parent), so this choice alone
  is worth tens of percent.
- **Parentless terms pass.** With `relative_p = 0` at the roots, the max()
  combination leaves root terms to the overall test alone — and they dominate
  output (on an all-species run, *biological_process* was the single most
  frequently associated term). We added a reporting-time information-content
  floor (`--min-ic 1`, IC = −log₂ of a term's propagated annotation
  frequency) to remove root and near-universal terms. The papers have no such
  mechanism; we suspect the original's meta-GO/IC stratification plays this
  role downstream.
- **Minimum background of 3 proteins**, rejected conservatively (p = 1) —
  our constant, counted and logged on every run (e.g. 11,073 pairs on the
  human GO run, 9,043 on HPO).
- **FDR < 0.01** rather than the paper's 10⁻³, and **separate BH families**
  for single vs supra domains. (We note S9's own worked example quotes an
  inherited annotation at FDR = 7.15×10⁻³ — above the paper's stated 10⁻³ —
  so we are unsure what the operational threshold was.)
- **Universe**: human-centred (and per-organism), InterPro entries — where
  the original was all-species SCOP superfamilies. All within-run contrasts
  below hold the universe fixed, so this difference cannot explain them.

## 3. What we expected

From the papers' framing we expected the relative test to act as a level
selector: prune child terms whose domain signal is fully explained by their
parents, retaining "only the most relevant" terms — a moderate reduction in
the association set, concentrated on redundant descendants, with protein-level
prediction quality (which is evaluated on propagated truth and so is largely
level-insensitive) roughly unchanged or slightly improved. In the language we
prefer here: a gain in the trustworthiness of what is kept (PPV) at a small
cost in coverage (sensitivity), with the balance neutral or favourable.

## 4. What we measured

All numbers below are from repository artifacts generated on current code
(after an August 2026 fix restricting GO propagation to `is_a`/`part_of`
edges) and shipped alongside this report; Section 8 has exact provenance.

### 4.1 Effect on the association set

Human GO, training data held fixed, FDR < 0.01:

| configuration | significant associations | with relative inference | kept |
|---|---:|---:|---:|
| supra-domains, direct input | 163,153 | 30,700 | 18.8% |
| supra-domains, propagated input | 486,041 | 91,190 | 18.8% |

Domains retaining at least one association: 27,192 → 10,539 and
27,300 → 11,629. The reporting IC floor is *not* the cause: on the
production paper-parity run it removed 963 of 97,382 post-BH associations
(1.0%); on HPO it removed zero.

### 4.2 Held-out temporal benchmark (protein-centric)

Design: train on archived human GOA release 205 (April 2021), transfer
domain associations to proteins (p-score, as published), evaluate against
experimental annotations gained by June 2026 on the CAFA no-knowledge
cohort (324 BP / 418 MF / 572 CC proteins); F_max (the best balanced
PPV/sensitivity score over the threshold sweep) and AUPRC per aspect at
IC floors 0/2/4; paired protein bootstrap (1,000 replicates) for every
contrast; scores are −log₁₀(q) for every configuration.

**Adding the relative test to the supra-domain model** (without input
propagation — the paper-parity configuration includes it): F_max and AUPRC
drop in **9/9 cells each, 8/9 significant**. The largest cell: MF at
IC ≥ 2, F_max 0.352 → 0.242, AUPRC 0.219 → 0.079.

**Adding it with input propagation in place** (the configuration closest to
the paper's): negative in every significant BP cell and in all but one MF
cell (F_max down in 7/9, AUPRC down significantly in 5/9), roughly neutral
for CC, with two isolated significant wins (MF/IC0 AUPRC +0.019; CC/IC2
F_max +0.016). Added to the otherwise-complete pipeline the picture is the
same.

**The rescue.** Adding input propagation to a relative-inference
configuration improves F_max in 9/9 cells (8 significant) and AUPRC in 8/9
(7 significant) — consistent with the view that the original method never
runs the relative test on an unpropagated input map, and that relative
inference without input propagation is a configuration the original never
computed. But even rescued, the paper-parity cells remain behind the
propagation-only configuration in most cells: the full pipeline beats
input+output propagation significantly only at MF/IC0 (AUPRC; its F_max
edge is not significant) and CC/IC2 (F_max).

### 4.3 Sensitivity / PPV decomposition at the operating point

*(Recall = sensitivity. Precision = positive predictive value, PPV — the
stand-in for specificity in this setting, where true negatives are not
countable: the fraction of predictions made that are correct.)*

Our working hypothesis (§3) was that the relative test buys specificity at
the cost of sensitivity — fewer predictions, but more trustworthy. The F_max
summary cannot distinguish that from losing on both axes, so we decomposed
each configuration's F_max operating point into its PPV and sensitivity
components, on the identical evaluation panels.

The decomposition runs on evaluation panels byte-equivalent to the §4.2
benchmark's (recomputed F_max matches all 90 of its cells to < 10⁻¹⁶; full
table in `validation/ablation_pr_decomposition.tsv`). Each configuration is
read at its own F_max operating point; Δ = (with relative test) − (without):

| contrast | cell | ΔPPV | Δsensitivity | ΔF_max |
|---|---|---:|---:|---:|
| added to supra | BP 0 | −0.037 | −0.059 | −0.053 |
| | BP 2 | −0.056 | −0.054 | −0.056 |
| | BP 4 | −0.055 | −0.044 | −0.049 |
| | MF 0 | −0.025 | −0.011 | −0.017 |
| | MF 2 | **−0.279** | −0.017 | −0.110 |
| | MF 4 | **−0.278** | −0.020 | −0.103 |
| | CC 0 | +0.059 | −0.067 | −0.030 |
| | CC 2 | +0.051 | −0.083 | −0.021 |
| | CC 4 | −0.044 | −0.087 | −0.068 |
| added to supra + input propagation | BP 0 | −0.009 | −0.049 | −0.030 |
| | BP 2 | −0.029 | −0.041 | −0.035 |
| | BP 4 | −0.023 | −0.063 | −0.041 |
| | MF 0 | +0.042 | −0.041 | −0.007 |
| | MF 2 | +0.002 | −0.058 | −0.044 |
| | MF 4 | −0.012 | −0.064 | −0.055 |
| | CC 0 | +0.013 | −0.008 | −0.001 |
| | CC 2 | +0.059 | −0.057 | +0.016 |
| | CC 4 | −0.022 | +0.009 | +0.001 |
| added to input + output propagation | BP 0 | +0.006 | −0.053 | −0.036 |
| | BP 2 | −0.035 | −0.039 | −0.037 |
| | BP 4 | −0.025 | −0.059 | −0.040 |
| | MF 0 | +0.095 | −0.030 | +0.017 |
| | MF 2 | −0.017 | −0.038 | −0.033 |
| | MF 4 | −0.016 | −0.043 | −0.038 |
| | CC 0 | +0.030 | −0.023 | −0.003 |
| | CC 2 | +0.061 | −0.034 | +0.016 |
| | CC 4 | +0.042 | −0.038 | −0.002 |

The pattern:

- **Sensitivity always falls** — 26 of 27 cells (the one exception, +0.009,
  is noise-level). The trade's cost side is universal.
- **The specificity payoff depends on the input map.** In the
  un-propagated configuration (no input propagation — one the original
  method arguably never ran), the test lowers **both** axes in 7 of 9 cells;
  in MF at IC ≥ 2 the PPV collapses by −0.28 and the operating threshold
  (the score cutoff τ) degenerates to τ ≈ 0 — the surviving score
  distribution is so thin that no high-PPV operating point exists at all.
- **With input propagation in place, the classic trade appears — sometimes.**
  PPV rises in 9 of 18 cells, most consistently in the cellular-component
  aspect (5/6) and at IC 0 (5/6); the largest net win is MF/IC0 on the
  full pipeline (PPV +0.095 for sensitivity −0.030, F_max +0.017; the two
  CC/IC2 cells also come out ahead, +0.016). But in the biological-process
  aspect the test loses PPV *and* sensitivity in nearly every cell — it
  removes true signal without making the remainder more trustworthy.
- Coverage falls in 20 of 27 cells, by up to ~10 points — fewer proteins get
  any prediction at the operating point. In the degenerate τ ≈ 0 MF cells
  (and one CC cell) it rises instead: the near-zero threshold lets
  everything through.

### 4.4 The extreme case: human phenotype annotations (HPO)

The HPO layer collapses from 996 single-domain associations to 38 under the
full paper-parity configuration. The production log attributes this
measurably: the IC floor removed nothing (43 kept, 0 dropped, counting both
families) and 13,742,781 of 14,829,673 evaluated pairs (92.7%) were governed
by the relative test (relative_p > overall_p); a separate measurement over
the propagated HPO map (recipe recorded in the project's evidence ledger)
found degenerate parents — child protein set equal to its parental union —
rare (4.1% of 4,000 sampled terms).
HPO's gene→phenotype input is disease-block derived and dense (332,599
rows over 5,274 genes), so most term-level signal is echoed across levels —
exactly where a parental-background test bites hardest. Whether 38 is the
correct specific core or over-pruning is precisely the question this report
is trying to settle.

### 4.5 What still works

For completeness: the underlying associations (without the relative test)
clear a 200-permutation domain-label null in 18/18 cells at the attainable
floor (p = 1/201), and the supra-domain combinations predict future
curation at 12.5× enrichment — the machinery around the relative test is in
good health, which is what makes its isolated negative effect stand out.

## 5. Interpretation (provisional)

Two readings fit the data so far:

- **(a) The test is working as designed and the benchmark punishes it.**
  Protein-centric CAFA evaluation propagates truth up the DAG, so predicting
  a parent instead of a child costs little — but *removing* a correct child
  (and having no term at that branch at all) costs sensitivity directly. A pure
  level-selector should relocate signal, not delete it; our max(p) rule
  deletes. If the original implementation *replaced* pruned children with
  their parents (the NAR paper's "only associate the parent term"), its
  output would keep the branch covered where ours leaves a hole. Our output
  propagation step partially compensates, which may be why the full pipeline
  is closer to neutral.
- **(b) Our background is not the original's.** If the original used the
  intersection ("all") background — smaller, hence more conservative
  per-pair p-values but in a *different pattern* — or a per-parent rule, or
  ran the test conditionally rather than always, the selected set would
  differ in ways we cannot reconstruct from the publications.

The decomposition (§4.3) splits the difference, informatively. Our working
hypothesis going in was the classic trade — specificity bought with
sensitivity. That holds only in the configurations closest to the paper's
(input propagated), and only for some aspects: there the test does behave
like a PPV filter, just one whose price in sensitivity usually exceeds the
gain. In the un-propagated configuration, and in the biological-process
aspect nearly everywhere, it loses on **both** axes — which is not what a
correctly-specified level selector should be able to do, and points at (a)'s
deletion-without-replacement mechanism, at (b), or both. The MF PPV collapse
at τ ≈ 0 in the un-propagated configuration is the sharpest single symptom:
after the max(p) combination, too few strong scores survive for any
high-PPV operating point to exist.

## 6. Questions for the original authors

1. **Union or intersection?** For a term with several direct parents, was
   N_pa the number of proteins annotated to *any* of them (union) or to
   *all* of them (intersection)? (S9 Methods and Figure 1's box say "all";
   Figure 1's caption says "any"; the worked example's N_pa = 21,240 looks
   like a union.)
2. **Root terms.** What did the relative test do at terms with no parents —
   skip (pass on overall alone), exclude, or something else? Did root/near
   -root terms appear in the released Domain2GO tables?
3. **Always-on or conditional?** Was the relative p computed for every
   (domain, term) pair and combined by max(), as S9's formulas suggest — or
   applied as a conditional decision rule ("if the child is not
   significantly different from the parent, associate only the parent"), as
   the NAR wording suggests? If the latter: did the child's signal transfer
   to the parent, or was it dropped?
4. **Input propagation.** Was the protein→term matrix explicitly true-path
   propagated before both tests? (The shared K in Figure 1 implies yes.)
5. **Thresholds in practice.** Was FDR < 10⁻³ the operational cutoff of the
   released tables? (S9's own example quotes an inherited annotation at
   FDR = 7.15×10⁻³.) Was there any minimum background size or minimum
   support?
6. **BH family.** Was the correction pooled across single domains and
   supra-domains (and SF/FA levels), or per family / per aspect?
7. **Effect size then.** Do you have (or remember) the before/after counts —
   what fraction of candidate associations the relative test removed in the
   original builds? Ours removes ~81%; does that match your experience?
8. **Is it still there?** The 2023 update's method summary does not mention
   the relative inference. Is it still part of the current dcGO build?

## 7. What we will do with the answers

If the original was conditional/parent-replacing (Q3) or intersection-backed
(Q1), we will implement that variant and re-run the same nine-cell paired
benchmark — the harness makes each variant a one-flag rerun. Answers to Q2
and Q4–Q6 each fix a parameter the harness already exposes (root handling,
input propagation, the FDR threshold, the BH family), so each is likewise a
one-flag rerun of the same benchmark; Q7 and Q8 need no experiment, only
comparison. If our recipe matches the original, the temporal result stands
as an honest finding about the method on modern data, and the paper will
present the relative test with its measured cost rather than as a default.

## 8. Provenance

- **Ablation** (nine configurations, paired bootstrap, permutation null):
  `validation/ablation_metrics.tsv`, `ablation_paired_bootstrap.tsv`,
  `ablation_permutation_null.tsv`, `ablation_selection_counts.tsv`,
  `ablation_provenance.tsv`, plus per-run manifests in
  `validation/ablation_manifests/` — all committed at `0088362` ("Land the
  current (2026-09-01) nine-configuration component ablation"); versions
  before that commit are an older, superseded experiment. The manifests
  internally record git commit `3f37558` as the code state of the runs;
  training GAF sha256 `69ae7d90…` = archived GOA release 205; evaluation on
  current GOA, experimental evidence only. Design and prose:
  `VALIDATION_PLAN.md` §4.
- **Rung regeneration for §4.3** (this report): the nine pipeline runs were
  regenerated from the run manifests' exact command lines and input
  hashes (commit `e5a73bd`, 2026-09-12) and reproduce the originals
  **byte-for-byte** — all 24 output files match the manifests' recorded
  SHA-256s (24/24), so the §4.3 decomposition is a decomposition of the
  original runs' results, not of a re-run that merely resembles them.
- **Production runs** (association-set counts, HPO case): 63-cell matrix of
  2026-08-18 at commit `d166013`, per-cell `run_manifest_*.json` with input
  and output SHA-256s; logs quoted are `results/production/go_paperparity.log`
  and `hpo_paperparity.log`.
- **Paper quotes**: page-referenced from the PDFs in `docs/`
  (`1471-2105-14-S3-S9.pdf`, `gks1080.pdf`, `EMS185259.pdf`), extracted with
  two independent tools; Figure 1 details read from the rendered image.
- **Sensitivity/PPV decomposition** (§4.3): `validation/ablation_pr_decomposition.py`
  (reuses the ablation's own panel machinery; recomputed F_max agrees with
  all 90 cells of the §4.2 benchmark to < 10⁻¹⁶), output
  `validation/ablation_pr_decomposition.tsv` (90 rows: 9 configurations +
  naive × 3 aspects × 3 IC floors, with PPV/sensitivity/coverage at each
  F_max operating point).
- **Implementation**: `src/relative_inference.py` (vectorised path
  cross-validated pair-for-pair against a scalar reference in the test
  suite), stage wiring in `run_dcgo_human.py`; every threshold named here is
  recorded in each run's manifest (`relative_inference_combination`,
  `parental_background_min_size`, `fdr_families`, `ic_source`).
