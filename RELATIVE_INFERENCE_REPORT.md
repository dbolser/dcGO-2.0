# The relative inference in dcGO-2.0: what we ran, what we expected, what we measured

*Pre-publication report, revised 2026-09-15. Prepared for discussion with the
original dcGO authors before any of this is written up for submission.*

---

## Executive summary

**What was run.** dcGO-2.0 is our reimplementation of the two-step
statistical inference of Fang & Gough (2013). Everything in this report is a
human-only run of that reimplementation:

- **Training annotations:** human UniProtKB-GOA, archived release 205
  (April 2021), manual (non-IEA) evidence codes — 18,735 annotated proteins.
- **Domains:** InterPro entries (InterPro `current_release` snapshot of
  2026-07-22), plus contiguous "supra-domain" combinations of up to three
  entries (as in the original method, which used SCOP superfamilies).
- **Ontology:** GO, `go-basic.obo` release 2026-06-15, `is_a`/`part_of`
  edges. Input annotations are true-path propagated before testing, as the
  original did.
- **Universe:** the 18,908 human proteins with at least one domain and at
  least one annotation.
- **Statistics:** one-sided Fisher tests, Benjamini-Hochberg at FDR < 0.01,
  single-domain and supra-domain hypotheses corrected as separate families.

Four configurations were run, identical except for two switches — whether
the **relative (parental-background) test** is on, and whether the
significant output is propagated up the DAG afterwards (the paper's Step 3):

| Configuration | Relative test | Output propagation |
|---|:---:|:---:|
| **Base** | off | off |
| **Base + relative** | on | off |
| **Base + output** | off | on |
| **Full** (the paper's complete method) | on | on |

Each was evaluated **held-out**: predictions transferred to proteins by
p-score (as published) and scored against experimental annotations those
proteins gained by June 2026 (GOA 2026-06-17), on the CAFA "no-knowledge"
cohort (324 BP / 418 MF / 572 CC proteins), by F_max and AUPRC per GO
aspect at information-content floors 0, 2 and 4, with 1,000 paired protein
bootstraps per contrast. Two additional **production** runs on current
(2026) inputs — not held out — supply the floor check in §4.1 and the HPO example in §4.4; they are
defined where used.

**What we found.**

1. **The relative test removes four fifths of the associations.** Base →
   Base + relative: 486,041 → 91,190 significant associations (18.8% kept);
   domains with at least one association 27,300 → 11,629.
2. **Held-out, it makes protein-level prediction worse in most cells.** Base
   → Base + relative: F_max lower in 7 of 9 aspect × IC cells (5
   significantly), AUPRC lower in 8 of 9 (5 significantly); one significant
   gain on each metric (CC at IC ≥ 2 for F_max; MF at IC ≥ 0 for AUPRC). The
   same pattern holds with output propagation on (Base + output → Full).
3. **It is not simply trading sensitivity for specificity.** Decomposing
   each configuration's operating point: sensitivity (recall) falls in 17 of
   18 contrast cells, but PPV (precision, the stand-in for specificity here)
   rises in only 9 of 18. The trade the test is meant to make does appear —
   in the cellular-component aspect and at MF IC ≥ 0 — but in biological
   process the test loses on *both* axes in 5 of 6 cells.
4. **The published description is ambiguous at the exact point our
   implementation had to choose.** The 2013 Methods text defines the
   relative background as proteins annotated to **all** direct parents;
   Figure 1's caption in the same paper says **any**; the worked example's
   numbers are consistent only with **any** (union). We implemented the
   union. If the original did something else — the intersection, a
   per-parent test, a conditional rule — that could explain much of the
   difference.

Section 6 lists seven specific questions; Section 8 gives provenance for
every number.

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
(pooled or split by domain level); and the exact evidence-code list behind
"experimental or manual".

We also note the 2023 update's method summary describes Fisher's exact test,
FDR and output propagation but does not mention the relative inference at
all — we cannot tell whether it was dropped from later builds or merely
elided in a simplified summary.

## 2. What we implemented

The relative test as it runs in every configuration that has it on (code in
`src/relative_inference.py`, wiring in `run_dcgo_human.py`; every parameter
is recorded in a machine-readable manifest per run).

For every (domain D, term T) pair that co-occurs in at least one protein:

1. **Overall test** — one-sided Fisher (enrichment), background = the
   18,908-protein universe.
2. **Relative test** — let B = the **union** of the (propagated) protein sets
   of T's *direct* parents.
   - no parents → `relative_p = 0` (the pair passes on the overall test alone);
   - |B| < 3 → `relative_p = 1` (counted, conservatively unpassable);
   - otherwise: one-sided Fisher on the 2×2 restricted to B.
3. **Combine** — `p = max(overall_p, relative_p)` *before* correction (a
   valid intersection-union p-value; Berger 1982), then Benjamini-Hochberg at
   FDR < 0.01 within each family.
4. **h-score** = min(overall, relative) hypergeometric score, as published.

**Choices the papers did not force** (candidate divergence points, in order
of expected impact):

- **Union background** ("any"), per Figure 1's caption and example — not the
  intersection ("all") of the Methods text. An earlier version of our code
  ran one test per parent and took the worst p; moving to the union raised
  human single-domain associations (Base settings) from 18,612 to 33,185 (56.1% of GO terms
  have more than one parent), so this choice alone is worth tens of percent.
- **Parentless terms pass.** With `relative_p = 0` at the roots, the max()
  combination leaves root terms to the overall test alone — and they dominate
  output (on an all-species run, *biological_process* was the single most
  frequently associated term). Our production runs therefore add a
  reporting-time information-content floor (`--min-ic 1`, IC = −log₂ of a
  term's propagated annotation frequency) that removes root and
  near-universal terms. The papers have no such mechanism; we suspect the
  original's meta-GO/IC stratification plays this role downstream. The four
  held-out configurations above do **not** apply the floor.
- **Minimum background of 3 proteins**, rejected conservatively (p = 1) —
  our constant, counted and logged on every run (11,073 pairs on the
  production human GO run, 9,043 on HPO).
- **FDR < 0.01** rather than the paper's 10⁻³, and **separate BH families**
  for single vs supra domains. (S9's own worked example quotes an inherited
  annotation at FDR = 7.15×10⁻³ — above the paper's stated 10⁻³ — so we are
  unsure what the operational threshold was.)
- **Universe**: human-only, InterPro entries — where the original was
  all-species SCOP superfamilies. Every contrast below holds the universe
  fixed, so this difference cannot explain them.

## 3. What we expected

From the papers' framing we expected the relative test to act as a level
selector: prune child terms whose domain signal is fully explained by their
parents, retaining "only the most relevant" terms — a moderate reduction in
the association set, concentrated on redundant descendants, with protein-level
prediction quality (which is evaluated on propagated truth and so is largely
level-insensitive) roughly unchanged or slightly improved. In the language we
use here: a gain in the trustworthiness of what is kept (PPV) at a small cost
in coverage (sensitivity), with the balance neutral or favourable.

## 4. What we measured

### 4.1 Effect on the association set

Base → Base + relative (same training data, same universe, FDR < 0.01):

| | Base | Base + relative | kept |
|---|---:|---:|---:|
| significant associations | 486,041 | 91,190 | 18.8% |
| domains with ≥ 1 association | 27,300 | 11,629 | 42.6% |

The information-content floor used in production is *not* the cause of the
reduction: on the production Full run (defined in §4.4) it removed 963 of
97,382 post-BH associations (1.0%); on HPO it removed zero.

### 4.2 Held-out temporal benchmark

Δ = (with relative test) − (without), each cell at its own F_max operating
point; \* marks a paired-bootstrap difference significant at 5%.

**Base → Base + relative**

| cell | F_max: Base → +relative | Δ | AUPRC: Base → +relative | Δ |
|---|---:|---:|---:|---:|
| BP 0 | 0.248 → 0.218 | −0.030\* | 0.123 → 0.093 | −0.031\* |
| BP 2 | 0.201 → 0.166 | −0.035\* | 0.075 → 0.054 | −0.021\* |
| BP 4 | 0.142 → 0.102 | −0.041\* | 0.034 → 0.015 | −0.018\* |
| MF 0 | 0.362 → 0.355 | −0.007 | 0.178 → 0.197 | +0.019\* |
| MF 2 | 0.363 → 0.318 | −0.044\* | 0.226 → 0.166 | −0.060\* |
| MF 4 | 0.346 → 0.291 | −0.055\* | 0.207 → 0.134 | −0.073\* |
| CC 0 | 0.376 → 0.375 | −0.001 | 0.220 → 0.220 | −0.000 |
| CC 2 | 0.206 → 0.221 | +0.016\* | 0.067 → 0.063 | −0.004 |
| CC 4 | 0.132 → 0.133 | +0.001 | 0.036 → 0.031 | −0.004 |

F_max lower in 7 of 9 cells (5 significantly), higher in 2 (1
significantly); AUPRC lower in 8 of 9 (5 significantly), higher in 1
(significantly). Every biological-process cell is significantly worse on
both metrics.

**Base + output → Full**

| cell | F_max: Base+output → Full | Δ | AUPRC: Base+output → Full | Δ |
|---|---:|---:|---:|---:|
| BP 0 | 0.258 → 0.222 | −0.036\* | 0.135 → 0.099 | −0.036\* |
| BP 2 | 0.203 → 0.166 | −0.037\* | 0.087 → 0.054 | −0.033\* |
| BP 4 | 0.139 → 0.100 | −0.040\* | 0.029 → 0.015 | −0.014\* |
| MF 0 | 0.368 → 0.385 | +0.017 | 0.198 → 0.215 | +0.017\* |
| MF 2 | 0.358 → 0.325 | −0.033\* | 0.226 → 0.175 | −0.051\* |
| MF 4 | 0.334 → 0.296 | −0.038\* | 0.204 → 0.142 | −0.062\* |
| CC 0 | 0.387 → 0.383 | −0.003 | 0.257 → 0.241 | −0.016\* |
| CC 2 | 0.218 → 0.233 | +0.016\* | 0.070 → 0.062 | −0.008 |
| CC 4 | 0.133 → 0.131 | −0.002 | 0.030 → 0.023 | −0.007 |

F_max lower in 7 of 9 (5 significantly), higher in 2 (1 significantly);
AUPRC lower in 8 of 9 (6 significantly), higher in 1 (significantly). Output
propagation does not change the picture: the relative test costs the same
cells either way.

### 4.3 Sensitivity / PPV decomposition at the operating point

*(Recall = sensitivity. Precision = positive predictive value, PPV — the
stand-in for specificity in this setting, where true negatives are not
countable: the fraction of predictions made that are correct.)*

Our working hypothesis (§3) was that the relative test buys specificity at
the cost of sensitivity — fewer predictions, but more trustworthy. F_max
cannot distinguish that from losing on both axes, so we split each
configuration's F_max operating point into its PPV and sensitivity
components on the identical evaluation panels (recomputed F_max matches
every §4.2 cell to < 10⁻¹⁶; full table in
`validation/ablation_pr_decomposition.tsv`).

| contrast | cell | ΔPPV | Δsensitivity | ΔF_max |
|---|---|---:|---:|---:|
| Base → Base + relative | BP 0 | −0.009 | −0.049 | −0.030 |
| | BP 2 | −0.029 | −0.041 | −0.035 |
| | BP 4 | −0.023 | −0.063 | −0.041 |
| | MF 0 | +0.042 | −0.041 | −0.007 |
| | MF 2 | +0.002 | −0.058 | −0.044 |
| | MF 4 | −0.012 | −0.064 | −0.055 |
| | CC 0 | +0.013 | −0.008 | −0.001 |
| | CC 2 | +0.059 | −0.057 | +0.016 |
| | CC 4 | −0.022 | +0.009 | +0.001 |
| Base + output → Full | BP 0 | +0.006 | −0.053 | −0.036 |
| | BP 2 | −0.035 | −0.039 | −0.037 |
| | BP 4 | −0.025 | −0.059 | −0.040 |
| | MF 0 | +0.095 | −0.030 | +0.017 |
| | MF 2 | −0.017 | −0.038 | −0.033 |
| | MF 4 | −0.016 | −0.043 | −0.038 |
| | CC 0 | +0.030 | −0.023 | −0.003 |
| | CC 2 | +0.061 | −0.034 | +0.016 |
| | CC 4 | +0.042 | −0.038 | −0.002 |

The pattern:

- **Sensitivity always falls** — 17 of 18 cells (the one exception, +0.009,
  is noise-level). The cost side of the trade is universal.
- **PPV rises in only 9 of 18 cells.** Where it does, the classic trade is
  real: cellular component in 5 of 6 cells, MF at IC 0 in both contrasts —
  the largest single gain is MF 0 in the Full configuration (PPV +0.095 for
  sensitivity −0.030, F_max +0.017).
- **In biological process the test loses on both axes** in 5 of 6 cells (and
  in MF at IC ≥ 2 in 3 of 4). There the relative test removes true signal
  without making the remainder more trustworthy — behaviour a level selector
  should not be able to produce.
- Coverage (the share of cohort proteins with any prediction at the
  operating point) falls in 13 of 18 cells, by up to 7.5 points.

### 4.4 The extreme case: human phenotype annotations (HPO)

Two production runs on current inputs (HPO `genes_to_phenotype.txt`, HP
ontology release 2026-06-23; not held out), 5,180 proteins in the universe:

- **Production Full** — the Full configuration plus the `--min-ic 1`
  reporting floor (this is the "paper-parity" production setting).
- **Production Direct** — overall test only: no relative test, no output
  propagation, no floor, and — unlike every configuration above — no
  propagation of the input annotations either. It is the repository's
  production baseline, so the contrast below spans more than the one switch.

Single-domain associations: Production Direct 996 → Production Full 38
(whole run, both families: 2,653 → 43). The production log attributes the
collapse to the relative test: the floor removed nothing (43 kept, 0
dropped), and 13,742,781 of 14,829,673 evaluated pairs (92.7%) were governed
by the relative test (relative_p > overall_p). A separate measurement over
the propagated HPO map found degenerate parents — child protein set equal to
its parental union — rare (4.1% of 4,000 sampled terms). HPO's
gene→phenotype input is disease-block derived and dense (332,599 rows over
5,274 genes), so most term-level signal is echoed across levels — exactly
where a parental-background test bites hardest. Whether 38 is the correct
specific core or over-pruning is precisely the question this report is
trying to settle; the clean single-switch evidence is the GO benchmark of
§4.2.

## 5. Interpretation (provisional)

Two readings fit the data so far:

- **(a) The test is working as designed and the benchmark punishes it.**
  Protein-centric CAFA evaluation propagates truth up the DAG, so predicting
  a parent instead of a child costs little — but *removing* a correct child
  (and having no term at that branch at all) costs sensitivity directly. A
  pure level-selector should relocate signal, not delete it; our max(p) rule
  deletes. If the original implementation *replaced* pruned children with
  their parents (the NAR paper's "only associate the parent term"), its
  output would keep the branch covered where ours leaves a hole. Output
  propagation partially compensates — but Full is still behind Base + output
  in 7 of 9 F_max cells, so it does not compensate enough.
- **(b) Our background is not the original's.** If the original used the
  intersection ("all") background — smaller, hence more conservative
  per-pair p-values but in a *different pattern* — or a per-parent rule, or
  ran the test conditionally rather than always, the selected set would
  differ in ways we cannot reconstruct from the publications.

The decomposition (§4.3) discriminates partially. Where the test behaves as
a PPV filter (cellular component, MF at IC 0) it is doing roughly what the
papers describe, at a price in sensitivity that usually exceeds the gain.
Where it loses on both axes (biological process throughout, MF at IC ≥ 2),
it is doing something the papers do not describe — which points at (a)'s
deletion-without-replacement mechanism, at (b), or both. Biological process
is where multi-parent terms are most common — 65% of BP terms have more
than one direct parent, against 50% in CC and 19% in MF (`go-basic.obo`
2026-06-15, `is_a`/`part_of`) — which is exactly where "all" and "any"
backgrounds diverge most; that is consistent with (b) but does not establish
it.

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
4. **Thresholds in practice.** Was FDR < 10⁻³ the operational cutoff of the
   released tables? (S9's own example quotes an inherited annotation at
   FDR = 7.15×10⁻³.) Was there any minimum background size or minimum
   support?
5. **BH family.** Was the correction pooled across single domains and
   supra-domains (and SF/FA levels), or per family / per aspect?
6. **Effect size then.** Do you have (or remember) the before/after counts —
   what fraction of candidate associations the relative test removed in the
   original builds? Ours removes ~81%; does that match your experience?
7. **Is it still there?** The 2023 update's method summary does not mention
   the relative inference. Is it still part of the current dcGO build?

## 7. What we will do with the answers

If the original was conditional/parent-replacing (Q3) or intersection-backed
(Q1), we will implement that variant and re-run the same nine-cell paired
benchmark — the harness makes each variant a one-flag rerun. Answers to Q2,
Q4 and Q5 each fix a parameter the harness already exposes (root handling,
the FDR threshold, the BH family), so each is likewise a one-flag rerun of
the same benchmark; Q6 and Q7 need no experiment, only comparison. If our
recipe matches the original, the held-out result stands as an honest finding
about the method on modern data, and the paper will present the relative
test with its measured cost rather than as a default.

## 8. Provenance

- **The four configurations** are the ablation runs `supra_input` (Base),
  `supra_input_relative` (Base + relative), `supra_input_output` (Base +
  output) and `full` (Full) of `validation/ablation.py`; metrics and paired
  bootstraps in `validation/ablation_metrics.tsv` and
  `ablation_paired_bootstrap.tsv`, per-run manifests in
  `validation/ablation_manifests/` — all committed at `0088362`. The
  manifests record git commit `3f37558` as the code state of the runs,
  training GAF sha256 `69ae7d90…` (= archived GOA release 205), domain
  annotations `data/interim/protein2ipr_human.dat.gz` sha256 `a932a515…`
  (InterPro `current_release` snapshot of 2026-07-22), and `go-basic.obo` sha256 `c72fc198…`
  (release 2026-06-15). Design and prose: `VALIDATION_PLAN.md` §4.
- **Regeneration check** (this report): the nine ablation pipeline runs were
  regenerated from the manifests' exact command lines and inputs (commit
  `e5a73bd`, 2026-09-12) and reproduce the originals **byte-for-byte** — all
  24 output files match the manifests' recorded SHA-256s (24/24).
- **Sensitivity/PPV decomposition** (§4.3):
  `validation/ablation_pr_decomposition.py` (reuses the ablation's own panel
  machinery), output `validation/ablation_pr_decomposition.tsv`.
- **Production runs** (§4.1 floor count, §4.4 HPO): 63-cell matrix of
  2026-08-18 at commit `d166013`, per-cell `run_manifest_*.json` with input
  and output SHA-256s; `results/production/go_paperparity.log` and
  `hpo_paperparity.log` are the logs quoted.
- **Paper quotes**: page-referenced from the PDFs in `docs/`
  (`1471-2105-14-S3-S9.pdf`, `gks1080.pdf`, `EMS185259.pdf`), extracted with
  two independent tools; Figure 1 details read from the rendered image.
- **Implementation**: `src/relative_inference.py` (vectorised path
  cross-validated pair-for-pair against a scalar reference in the test
  suite), stage wiring in `run_dcgo_human.py`; every threshold named here is
  recorded in each run's manifest (`relative_inference_combination`,
  `parental_background_min_size`, `fdr_families`, `ic_source`).
