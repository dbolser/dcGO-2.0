# Knowledge projection: making dcGO ~10× better

Goal (Dan, 2026-11-05): make the method about ten times better, by any mix of
a new approach on the same core idea, new data (BioEmbed), techniques from the
literature, and measures that reward *robustly projecting knowledge*.

## 1. Pre-registration (written 2026-11-05, before any number below existed)

### 1.1 What "projecting knowledge" means here

dcGO's job is to carry what is known about studied proteins to proteins
nobody has studied, through the parts they share. The measure that matches
that job is: **how many correct new facts can the method hand over while
keeping its error rate fixed.** A list where half the claims are wrong is
still usable by a curator; a list where 95% are wrong is not, however good
its F_max.

### 1.2 Primary criterion: knowledge yield at 50% precision (KY50)

On the CAFA no-knowledge temporal cohort (the §4 cohort: `manual` training
evidence, experimental truth at t1, True-Path propagated, aspect roots
excluded, cohort = proteins with a domain in the cohort `protein2ipr`):

- restrict both truth and predictions to informative terms, **IC ≥ 2**
  (IC from propagated t0 training frequencies, as `floor_cell` does);
- **pool** every (protein, term) prediction across the aspect's cohort and
  sweep a threshold over the distinct scores;
- `TP(τ)` = predictions with score ≥ τ that are in the protein's truth,
  `PP(τ)` = all predictions with score ≥ τ;
- **KY@π = max TP(τ) over thresholds with TP(τ)/PP(τ) ≥ π.** Primary π = 0.5.

Reported alongside, not tallied: KY@0.3 and KY@0.7; the IC-weighted yield
(bits) and the number of distinct proteins with ≥ 1 correct prediction at the
KY50 threshold; pooled (micro) AUPRC; the existing protein-centric F_max,
S_min and AUPRC under CAFA conventions.

Aspects: **BP and MF** primary; CC reported as sensitivity.
Uncertainty: paired bootstrap over cohort proteins (1,000 resamples), CI on
the ratio new / old per aspect.

### 1.3 Baseline

"Old dcGO" = the **best** of every rung of the 2026-10-15 §4 ladder
(`results/ablation-t0ipr-2026-10-15/ipr_manual/*`) under both transfers
(`pscore`, `max`) and both score columns (q, p), chosen **per aspect, on the
split being evaluated, by KY50**. This selection favours the baseline.
The naive baseline is reported too.

### 1.4 What counts as 10×

KY50(new) / KY50(best old) ≥ 10 in both BP and MF on the final split, point
estimate, with the bootstrap interval reported next to it. A baseline with
fewer than 10 true positives at a level is reported as "≈ 0", not as a ratio.

Two arms are always reported separately, because they project different
kinds of knowledge:

- **domain-only** — sequence-derived units only (InterPro entries and their
  combinations); it applies to any protein in any species;
- **domain + literature** — adds BioEmbed literature evidence published
  before t0; it applies only to proteins with papers.

### 1.5 Splits and protocol

- **Development:** t0 ≈ mid-2018 (GOA human + the InterPro release of the
  same month), t1 = GOA 205 (2021-04). Literature: papers with year ≤ 2017.
  Every design choice (shrinkage strength, combination rule, features,
  stacker settings) is made here or by cross-validation inside t0.
- **Final:** t0 = GOA 205 + InterPro 85.0 (2021-04), t1 = current GOA
  (2026-07), cohort `protein2ipr` current — the split all existing §4 numbers
  are on. Literature: papers with year ≤ 2020. Run **once** with the frozen
  configuration. Anything changed afterwards is reported as a post-hoc
  deviation.
- The project's own pre-registered endpoints (`TODO.md` top) are reported as
  well, unchanged.

### 1.6 Addendum, frozen after development and before the final run (2026-11-06)

Everything here was decided on the 2018 → 2021 development split or by
cross-validation inside t0, before the final split was run. It does not
change §1.2–1.5.

**Frozen configuration** (`validation/projection_benchmark.py` flags):

| Setting | Value | Why (dev evidence in §2.3–2.4) |
|---|---|---|
| Units | InterPro entries only (`--unit-sets ipr`) | signatures hurt BP (IBA circularity); supra/architecture no gain |
| θ | posterior mean, `s = 0.5` | lower bound worse under CV; s 0.5–1 equal |
| Training | human t0 `manual` (as old dcGO) **plus** all-species **experimental** annotations at t0 (`--extra-train-gaf`, `--extra-evidence-filter experimental`) | best summed BP+MF KY50 on dev; avoids IBA/ISS circularity |
| Sequence neighbours | ProtT5 (UniProt per-protein, human + Swiss-Prot), k = 25, T = 0.1 | best single-channel MF/CC on dev |
| Literature neighbours | MedEmbed, papers ≤ t0 year − 1, centred, same-gene vectors dropped, k = 25, T = 0.1 | |
| Combination | equal-weight mean of the channels | beat noisy-OR in BP on dev; equal in MF |

**Arms reported** (all three, plus each channel alone):
domain-only = `proj`; sequence = `mix(proj=1,seq=1)`; full =
`mix(proj=1,seq=1,lit=1)`. The pre-registered domain-only arm is the one
that applies to any protein in any species; the full arm applies to
proteins with papers.

**Analyses added** (reported, not tallied):
- *Deployable cutoff:* each method's KY50 score cutoff on dev is frozen and
  applied on the final split; precision and yield there are reported
  (`--frozen-tau-from`).
- *Non-redundant yield* (`ky50_leaf`): correct predictions with no correct
  descendant, so an ancestor chain counts once.
- *Comparators:* InterPro2GO(t0) projected through each protein's entries,
  and every t0 IEA annotation.
- *High-throughput sensitivity:* truth without papers annotating > 99
  proteins (`--max-proteins-per-paper 99`).
- The project's association-level endpoint, with the new estimator ranked
  by lower credible bound over each unit's most specific confident terms.

**Expected outcome, stated now:** on dev the old method's KY50 was ≈ 0 in
BP (7) and CC (0) and 89 in MF; the frozen full arm gave 656, 241 and 370.
So the pre-registered "≥ 10× in BP and MF" is not expected to hold for MF
(≈ 4×), where 10× would need ~80% recall at 50% precision.

### 1.7 Pre-final adversarial review (2026-11-06, before the final run)

Four reviewers (leakage, metric, fairness, algorithm), each finding checked
by a separate skeptic who tried to refute it, on the dev split only. All
confirmed findings were rated minor; none changes a pre-registered verdict.
What changed as a result, before the final run:

- **t0-existence gate.** Embedding files are current, so the sequence and
  literature channels could score cohort proteins that had no record at t0
  (demerged or new accessions: AMY1A/B/C, UGT2A1/2, H3Y1/2 …), which old
  dcGO cannot. Every channel now scores only cohort proteins present in the
  t0 InterPro subset or the t0 GAF.
- **Same-gene guard.** Human accessions in the all-species background that
  share a cohort protein's gene symbol are dropped, so a protein's own gene
  never trains for it (HLA-A at t0 on allele entries, ARHGEF10L under an
  older accession). Gate and guard together moved the dev full arm's KY50
  BP 656 → 681, MF 370 → 328, CC 232 → 196.
- **Sensitivity rows added:** the cohort without proteins whose *gene* had a
  same-aspect manual annotation at t0 under any accession
  (`--drop-gene-known`; on dev it moves MF's best old rung from 89 to 133),
  and without high-throughput papers (§1.6).
- **Bugs fixed:** bootstrap ratio replicates where both methods yield 0 were
  counted as +∞ (now dropped); the kNN weight normaliser underflowed to NaN
  in float32 when every neighbour was excluded; the driver now refuses an
  old ladder trained on other inputs (as `ablation.py` does); the final run
  emits only the frozen arms and single channels.
- **Kept as frozen, after testing:** the weighted mean counts a missing
  channel as 0. Renormalising over present channels lowered the full arm's
  dev KY50 (BP 656 → 537, MF 370 → 344), so it is a sensitivity option
  (`--mix-skip-missing`), not the rule.
- **Known and left:** the bootstrap holds "best old" fixed at its point
  selection (as pre-registered); the association-level path does not remap
  retired GO ids on either side (as `ablation.py`); KY50 is jumpy (§2.5).

**Dev numbers of the frozen configuration after these changes**
(`validation/projection/dev2018_frozen/`, IC ≥ 2, 1,000 paired resamples):
KY50 best old → full arm BP 7 → 681, MF 89 → 328 (3.7× [2.2, ∞]), CC 0 →
196; pooled AUPRC ratio BP 5.0× [3.8, 6.9], MF 3.9× [2.7, 6.0], CC 6.9×
[4.6, 11.0]. Without gene-known proteins MF's best old falls to 6; without
high-throughput truth it rises to 141 (ratio 2.4×). The final run follows.

## 2. Method and results

*(Filled in after the runs. Sections below this line were written later.)*

### 2.1 What changed, in one paragraph

dcGO asks "is this domain–term pair more common than chance?" and projects
the answer (a p-value) onto proteins. Projection needs a different number:
"if a protein carries this unit, how often does it have this term?" The new
method estimates exactly that, per unit and term, on the proteins annotated
in that GO aspect, shrunk towards the term's base rate when the evidence is
thin. A protein's score for a term is the best estimate among its units, so
the most specific unit it shares with annotated proteins speaks for it.
Scores mean the same thing on every protein, so all predictions can be put
in one list and cut at a chosen precision. Two more channels feed the same
scale: a protein's nearest neighbours by sequence embedding (ProtT5) and by
the papers written about it (BioEmbed's MedEmbed vectors, papers before t0
only).

### 2.2 Components (code)

| Component | Where | Idea | Source |
|---|---|---|---|
| Unit projection `proj` | `src/knowledge_projection.py` | `θ(u,t) = (k + s·π)/(n + s)`, per aspect, max over a protein's units | dcGO's units; empirical-Bayes shrinkage |
| Lower-bound variant | `UnitTermEvidence.posterior(quantile=…)` | Beta lower credible bound instead of the mean | robust to thin support |
| Unit families | `unit_families`, `read_restricted_units` | InterPro entries, bag supra-domains, whole architecture, member signatures incl. PANTHER subfamilies | dcGO; PANTHER/PAINT |
| LR-InterPro `lr` | `lr_project` | one logistic regression per term over unit indicators | GOLabeler (Zhu et al. 2018) |
| Literature neighbours `lit` | `knn_project`, `load_literature` | similarity-weighted vote of the k nearest annotated genes by pre-t0 papers | NetGO 2.0 LR-Text; GORetriever retrieval stage |
| Sequence neighbours `seq` | `load_h5_embeddings` | same vote in ProtT5 per-protein embedding space | goPredSim (Littmann et al. 2021) |
| All-species background | `--extra-train-gaf` | estimate `θ` on every species' manual annotations at t0 | `MULTISPECIES_BACKGROUND.md` |
| Combination | `noisy_or`, `--mix` | noisy-OR or weighted mean of channel scores | GOLabeler consensus / weighted vote |
| Learned stacker | `src/projection_stacker.py` | monotone gradient boosting over evidence features | GOLabeler / NetGO LTR — **dropped on dev, see 2.4** |
| Knowledge yield | `validation/knowledge_yield.py` | pooled precision-controlled yield, bootstrap | this document §1.2 |
| Driver | `validation/projection_benchmark.py` | old ladder, naive and new methods on one cohort | — |

### 2.4 What did not work on dev (2018 → 2021), and was dropped

All numbers: KY50 at IC ≥ 2, dev cohort (BP 337, MF 206, CC 220 proteins).

- **A learned stacker** (GOLabeler-style; monotone gradient boosting over 17
  evidence features, cross-fitted over cohort proteins): BP 19 vs 188 for
  the fixed rule, MF 207 vs 248. Dropping protein-level features or using
  logistic regression did not close the gap. With ~170 training proteins
  per fold the model learns proteins, not evidence. Code kept
  (`src/projection_stacker.py`) but not used.
- **Member-database signatures as units** (incl. PANTHER subfamilies): BP
  108 vs 240 for InterPro entries alone; MF level. Likely cause, from the
  literature sweep: training includes IBA annotations made by PAINT on
  PANTHER families, and the no-knowledge cohort is exactly the proteins
  PAINT did not annotate, so siblings' IBA labels overstate θ for them.
- **Supra-domains and whole architectures as units**: no gain in KY50
  (BP 156 vs 236 at s = 1), as in every earlier protein-level ablation.
  They stay in the code; the frozen configuration does not use them.
- **The Beta lower credible bound** instead of the posterior mean: worse
  under cross-validation inside t0 (BP 155k vs 169k KY50 at q = 0.1), mixed
  on dev. The posterior mean is kept.
- **Existing electronic annotation as a shortcut**: InterPro2GO(t0)
  projected through each protein's entries, and every IEA annotation in the
  t0 GAF, both give KY50 = 0 in BP and MF (one block, precision below 0.5).

### 2.5 KY50 is a jumpy statistic

KY50 reads the largest yield whose pooled precision clears 0.5. When the
top of a list sits near 0.5, a resample can flip it across and the yield
collapses: on dev, domain-only BP has a bootstrap interval of [4, 378]
around 108. Pooled AUPRC moves smoothly and is the steadier companion;
both are reported, and the ratio claims below quote both.

## 3. Final split (2021 → 2026), run once on 2026-11-06

`validation/projection/final/` (primary), `final_geneknown/`,
`final_htfilter/` (sensitivities), `final_assoc/` (association level);
`scripts/run_projection_benchmark.sh final …` with the §1.6 configuration
and `--frozen-tau-from validation/projection/dev2018_frozen/…`. The old-ladder
input guard passed; 11–13 cohort proteins had no t0 record and got no
embedding-based score. Tables: `uv run python validation/summarize_projection.py
validation/projection/final`.

### 3.1 Verdict on the pre-registered criterion

**Not met.** KY50, IC ≥ 2, best old rung (chosen on this split) vs the full
arm, 1,000 paired protein resamples:

| aspect | best old | full arm | ratio [95% CI] | pre-registered ≥ 10×? |
|---|---|---|---|---|
| BP | 54 | 438 | 8.1× [3.4, ∞] | no (point 8.1) |
| MF | 259 | 360 | 1.4× [1.1, 2.0] | no |
| CC (sensitivity) | 43 | 198 | 4.6× [1.7, ∞] | — |

MF could not have passed with any method: old dcGO already yields 259 of the
887 truth pairs, so even perfect recall at 50% precision caps the ratio at
3.4×. The pre-registration should have checked that ceiling on the final
cohort's baseline; it did not, and that is a design error in §1.4.

### 3.2 Everything else measured (reported, not tallied)

| aspect | arm | KY30 | KY50 | KY70 | pooled AUPRC (ratio [CI]) | F_max |
|---|---|---|---|---|---|---|
| BP | best old | 147 | 54 | 0 | 0.028 | 0.184 |
| BP | domain-only | 343 | 1 | 1 | 0.085 (3.0× [2.0, 5.2]) | 0.244 |
| BP | sequence | 696 | 94 | 1 | 0.137 (4.9× [3.2, 8.1]) | 0.280 |
| BP | full | 982 | 438 | 217 | 0.175 (6.2× [4.1, 10.3]) | 0.293 |
| MF | best old | 293 | 259 | 0 | 0.214 | 0.333 |
| MF | domain-only | 456 | 382 | 239 | 0.377 (1.8× [1.4, 2.3]) | 0.429 |
| MF | sequence | 511 | 377 | 242 | 0.426 (2.0× [1.6, 2.6]) | 0.446 |
| MF | full | 487 | 360 | 276 | 0.435 (2.0× [1.6, 2.6]) | 0.431 |
| CC | best old | 54 | 43 | 0 | 0.037 | 0.174 |
| CC | domain-only | 273 | 32 | 24 | 0.117 (3.2× [2.0, 5.9]) | 0.310 |
| CC | sequence | 638 | 114 | 19 | 0.220 (6.0× [3.8, 11.1]) | 0.409 |
| CC | full | 652 | 198 | 17 | 0.247 (6.8× [4.2, 12.2]) | 0.413 |

F_max: the old row is the best old rung *by F_max* under the p-score
transfer (BP and CC: `supra_input` and `supra_input_output` tie; MF:
`supra_input`, p-value column), which reproduces the committed §4 cell
exactly (BP 0.185, MF 0.333, CC 0.174); CAFA conventions, IC ≥ 2. Naive: BP
0.139, MF 0.076, CC 0.230.

- **At 70% precision the selected old rung yields nothing** in any aspect
  (the best of all 28 old variants yields 20 in BP, none in MF or CC); the
  full arm yields 217 / 276 / 17 correct informative predictions.
- **Deployable cutoff** (dev KY50 cutoff applied unchanged): full arm
  precision 0.46 BP (512 correct), 0.54 MF (333), 0.50 CC (213) — the dev
  calibration carried over. Domain-only BP did not (0.25).
- **Comparators:** InterPro2GO(t0) and every t0 IEA annotation have KY50 = 0
  in BP and MF (one block below 50% precision); t0 IEA reaches 275 in CC.
- **Sensitivities** change little: without gene-known proteins the KY50
  ratios are BP 6.6×, MF 1.3×, CC 4.3×; without high-throughput papers BP
  8.1×, MF 1.4×, CC 4.5×. Pooled AUPRC ratios stay ~6× (BP, CC) and ~2× (MF).
- **Domain-only BP KY50 collapsed to 1** while its AUPRC tripled: the
  knife-edge of §2.5. KY50 alone should not be used to rank arms.

### 3.3 Association level (the project's own pre-registered primary)

Precision at the fixed budget (4,024) against InterPro2GO pairs added 85.0 →
current, single InterPro entries:

| table | raw ids | live ids |
|---|---|---|
| new, lower-bound leaves | 0.113 | 0.114 |
| new, posterior mean | 0.540 | 0.042 |
| old `supra_input_relative` (the pre-registered reference) | 0.090 | 0.105 |
| best old by live-id precision | `supra_input_relative` 0.090 | same, 0.105 |

The raw 0.540 is an artefact, and so is most of the raw endpoint: of 5,534
"added" pairs only 2,054 remain once both InterPro2GO snapshots are put on
live GO ids — the rest are mappings to a GO id retired since t0, re-made on
its successor. On live ids the new method ties the old one. **This affects
the committed §4 primary-endpoint results too** (`validation/ablation.py`
scores the raw ids), which should be rescored on live ids.

### 3.4 What this means

The new method is reliably better at handing over correct knowledge on
informative terms — about 6× by pooled AUPRC in BP and CC, 2× in MF, and it
produces usable lists at 70% precision where old dcGO produces at most 20
correct predictions (BP) and none in MF or CC. The paper write-up is
`paper/knowledge_projection.md`. Most of
the gain comes from three things old dcGO did not have: scores that mean the
same thing on every protein, all-species experimental evidence, and two new
evidence channels (ProtT5 sequence neighbours; pre-t0 literature via
BioEmbed). It is not ten times better by the criterion fixed in advance.

## 4. Corrections to this record (2026-11-06, after the final run)

§1 and its addenda are left as written; these notes correct them. Where
this file and `paper/knowledge_projection.md` disagree, the paper (whose
tables are generated from the result files) is right.

- **§1.6, expected outcome.** The dev full arm before the review was
  656 / 370 / **232** (BP / MF / CC; `dev2018_pass5_allsp_exp`), not
  "656, 241 and 370". §1.7's "CC 232 → 196" is the correct value.
- **§1.7, gene-known sensitivity.** "Moves MF's best old rung from 89 to
  133" came from a run before the existence gate and same-gene guard and
  has no surviving result file. The frozen dev run
  (`dev2018_frozen_geneknown`) gives MF's best old KY50 = 6, as the end of
  §1.7 says.
- **§3.2, F_max.** The replayed ladder reproduces the earlier F_max to
  within 0.004, not exactly; the best old BP F_max on the final split is
  0.184.
- **§3.2 table and §3.4, AUPRC ratios.** Those ratios compare every metric
  with the old variant best by KY50. Against the old variant best on each
  metric (as the paper reports), full-method pooled AUPRC is 4.5× in BP,
  2.0× in MF and 6.8× in CC, not 6.2× / "about 6×" in BP.
- **The headline arm.** The pre-registered arms (§1.4) are domain-only and
  domain + literature. The full arm (domain + sequence + literature) was
  added in §1.6; the paper lists this as a departure and reports the
  pre-registered arm next to it (BP 7.1×, MF 1.3×).
- **BP domain-only failure.** §2–3 attribute it to thinly supported units.
  Post hoc, human-only training gives BP KY50 241 on the final split, so
  the all-species background is the cause; units carried only by other
  species' proteins are about half of the top BP predictions and mostly
  wrong (`final/reliability_by_origin.tsv`), but see the follow-up below.
- **§2.4, learned combiner.** "BP 19 vs 188, MF 207 vs 248" compares with
  one fixed rule (all unit families + literature, s = 1, T = 0.2). The best
  fixed rule over the same units and literature in `dev2018_pass1` gives
  227 (BP) and 272 (MF); with InterPro entries alone, 274 and 280.
- **BP domain-only failure, follow-up.** Removing every unit with no human
  carrier (`paper/background_origin_ablation.py`,
  `final/background_origin_ablation.tsv`) raises BP KY50 only from 1 to 15
  (human-only 241; paired difference human-only − frozen 240 [51, 403]), so
  the other-species-only units are part of the cause, not all of it; the
  background also lowers the scores of the predictions human training ranks
  highest.
