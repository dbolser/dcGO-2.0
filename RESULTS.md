# dcGO-2.0 — Results Summary (§2 temporal benchmark)

*One-page summary of the validation findings. Full detail in
[VALIDATION_PLAN.md](VALIDATION_PLAN.md) §2; numbers from
`validation/temporal_benchmark_metrics.tsv` and
`validation/domain_centric_metrics.tsv`.*

## Headline

Trained on a **2021** GOA snapshot and tested against **newly-curated,
experimentally-supported 2026** annotations (a real, precision-capable CAFA-style
test), **dcGO beats the naive frequency baseline on F_max at face value in BP and
CC**, and beats it in **every aspect once uninformative terms are excluded**,
while staying **1.3–25× above the random-domain null**. Molecular Function is the
one aspect where naive leads at face value — because its "truth" is dominated by a
single near-universal term (`protein binding`: 84.6% of experimental MF
annotation lines, carried by 87.7% of annotated proteins, and the *only* MF
term for 34.5% of them — see `validation/protein_binding_dominance.py`) — and dcGO overtakes it decisively the moment that noise is removed.

**The MF pattern does not reproduce in mouse.** Evaluated on mouse — an
organism that influenced no method choice, over a matched 2021-04 → 2026 window
— dcGO *beats* naive on MF F_max at face value (0.471 vs 0.359) and takes MF
AUPRC too (0.304 vs 0.256). The two species also differ in how concentrated
their MF annotation is: `protein binding` is 84.6% of human experimental MF
annotation lines but only 48.7% of mouse's, and across the pair naive's MF
F_max falls (0.464 → 0.359) as dcGO's rises (0.358 → 0.471).

That is *consistent with* the saturation explanation but does not establish it.
Two species is an association, not a demonstration of cause, and mouse differs
from human in curation depth and practice as well as in term concentration. The
supportable statement is narrower: the MF result is not a stable property of the
method — it does not survive a change of organism.

**This headline is about F_max only.** On AUPRC, from the same run, naive leads
in all three aspects at IC≥0 and in CC at IC≥2 and IC≥4 as well
([VALIDATION_PLAN.md](VALIDATION_PLAN.md) §2). The primary endpoint was not
pre-specified, so read the above as a descriptive retrospective result at one
operating point, not as a general claim of superiority.

## Setup

- **Split:** train = GOA release 205 (2021-04); test = current (2026-06), ~5-year gap.
- **Benchmark:** CAFA *no-knowledge* — per aspect, proteins with no annotation
  known to training at t0 that gained experimental annotation by t1, scored
  against their full propagated t1 experimental terms. Gate uses the **training
  evidence** (`manual`/non-IEA) to avoid leaking already-seen labels.
- **Predictor:** domain→GO associations transferred to proteins via the per-target
  **p-score** (Fang & Gough 2013: sum of scores, min-max normalised per protein).
- **Metric:** protein-centric **F_max** over a threshold sweep, plus an
  **information-content (IC) floor** applied identically to truth and every method,
  to separate informative predictions from base-rate recovery of generic terms.

Benchmark sizes (IC≥0): **BP 324 / MF 418 / CC 572** proteins. The IC floors
change the cohort as well as the truth — MF keeps 170 of 418 proteins at IC≥2 —
so rows at different floors are not paired comparisons.

## Protein-centric result — F_max (dcGO / naive / permutation null)

*Regenerated 2026-08-05 on the final code (#44–#47, #50): shrinkage removed,
single and supra-domain BH families corrected separately, True Path background
propagated. Every row here and in the mouse section below comes from that one
pipeline version. The "random" column is now the mean of a **100-permutation**
random-domain null rather than a single shuffle, so the ratios supersede the
earlier single-draw ones — MF IC≥2 is 7.4×, not the 4.2× previously published.*

| Aspect | IC floor | dcGO | naive | permutation null (mean) | dcGO ÷ null |
|--------|:--------:|-----:|------:|------------------------:|------------:|
| **BP** | ≥0 | **0.251** | 0.115 | 0.156 | 1.6× |
| BP | ≥2 | **0.172** | 0.071 | 0.054 | 3.2× |
| BP | ≥4 | **0.115** | 0.031 | 0.017 | 6.8× |
| BP | ≥6 | **0.077** | 0.010 | 0.006 | 13.6× |
| **MF** | ≥0 | 0.358 | **0.464** | 0.193 | 1.9× |
| MF | ≥2 | **0.366** | 0.053 | 0.050 | 7.4× |
| MF | ≥4 | **0.338** | 0.045 | 0.037 | 9.1× |
| MF | ≥6 | **0.218** | 0.018 | 0.010 | 21.3× |
| **CC** | ≥0 | **0.382** | 0.343 | 0.279 | 1.4× |
| CC | ≥2 | **0.242** | 0.153 | 0.069 | 3.5× |
| CC | ≥4 | **0.138** | 0.099 | 0.030 | 4.7× |
| CC | ≥6 | **0.124** | 0.044 | 0.019 | 6.6× |

*IC in bits; IC≥2 ⇒ term in ≤25% of proteins, IC≥6 ⇒ ≤1.6%. Bold = winner
between dcGO and naive.*
**Naive is a base-rate mirage:** its F_max collapses toward the null as
informative terms are required (BP 0.115→0.010, MF 0.464→0.018, CC 0.343→0.044
over IC 0→6). dcGO degrades gracefully because its predictions are specific.
Note that at IC≥0 in BP the null itself (0.156) *exceeds* naive (0.115) — a
shuffle carrying no domain information at all outscores the standard baseline,
purely by reproducing each term's marginal frequency.

> **Four corrections from the §4 uncertainty work (2026-08-04).** Read them
> before citing this table. Detail in `VALIDATION_PLAN.md` §2 (permutation null)
> and §4 (ablation).
>
> 1. **The `random` column is a single shuffle.** It is now a distribution over
>    100 seeded permutations. dcGO clears it everywhere with the smallest
>    attainable empirical p (1/101), but two of the ratios above were drawn from
>    the null's upper tail and understate the margin: **MF ≥2 is 7.5×, not 4.2×;
>    MF ≥4 is 9.0×, not 4.7×** (against the null mean). Use
>    `validation/temporal_benchmark_permutation_null.tsv`.
> 2. **The F_max wins survive a paired protein-level bootstrap; the AUPRC wins
>    do not.** Naive *significantly beats* dcGO on AUPRC at IC≥0 in all three
>    aspects, and at every IC floor in CC — it predicts every term for every
>    protein and so owns the high-recall end of the curve. "Beats naive in every
>    aspect on informative terms" is an **F_max** statement only.
> 3. **Coverage.** dcGO makes no prediction at all for ~45% of the benchmark
>    cohort (BP 0.55, MF 0.63–0.71, CC 0.42–0.52 coverage); naive covers 100%.
>    CAFA recall already charges for this, but the F_max values above should be
>    read against it.
> 4. **These numbers are from a pipeline version that has since changed.**
>    Re-running the t0 training command on current `main` gives 163,277
>    significant associations, not the 164,549 behind this table, because
>    `restrict_to_universe` now narrows the Fisher protein universe to
>    domain-annotated proteins (18,735 → 18,382). The current behaviour is the
>    correct one; the table is an artefact of the earlier version.

## What each pipeline component contributes (ablation, 2026-08-04)

Paired protein-level bootstrap, 1,000 replicates, 12 aspect × IC-floor cells.
Full detail and the mechanisms in `VALIDATION_PLAN.md` §4.

| Component | Cells where it helps | Cells where it hurts | Verdict |
|---|---:|---:|---|
| Supra-domains | **0 / 12** | 1 / 12 | No measurable effect, at 5.3× the feature space |
| Hierarchical shrinkage | **0 / 12** | **0 / 12** | No effect on prediction — but it takes the "significant" association count from 163,277 to **463,924** by *decreasing* 56% of supra p-values. Not a shrinkage; BH on its output does not control FDR. |
| True Path Rule | 0 / 12 | **12 / 12** | Significantly **worse** (−0.04 to −0.24 F_max). Its parental-background filter keeps 14% of associations, and 54,951 parent tests are rejected *untested* because the background uses the unpropagated annotation map. |

**On this benchmark the best configuration is the simplest one: single domains,
no shrinkage, no True Path.** `full − single` is significantly negative in all
twelve cells. The supra-domain machinery's demonstrated value is elsewhere — the
*emergent* combinations validated in `SURPRISE_SCORE.md` — not in
protein-centric F_max.

## Held-out evaluation on an untouched species (2026-08-05)

Every method choice — transfer rule, IC thresholds, the since-removed shrinkage
— was compared on the human 2021→2026 split that also carries the numbers
above. Mouse influenced none of them, so it is genuinely untouched, and it
keeps the same five-year window. GOA release numbers are per-species: mouse
release **191** (2021-04-08) is the match for human **205** (2021-04-21), not
mouse 205, which is dated 2023-09.

| Aspect | IC | human dcGO | mouse dcGO | mouse naive |
|--------|:--:|-----------:|-----------:|------------:|
| BP | ≥0 | 0.251 | **0.275** | 0.114 |
| BP | ≥2 | 0.172 | **0.189** | 0.078 |
| BP | ≥4 | 0.115 | **0.147** | 0.036 |
| MF | ≥0 | 0.358 | **0.471** | 0.359 |
| MF | ≥2 | 0.366 | **0.481** | 0.072 |
| MF | ≥4 | 0.338 | **0.469** | 0.045 |
| CC | ≥0 | 0.382 | **0.403** | 0.403 |
| CC | ≥2 | 0.242 | **0.294** | 0.125 |
| CC | ≥4 | 0.138 | **0.299** | 0.069 |

**Performance does not degrade off the tuned split — it improves in all nine
cells**, and clears the 100-permutation random-domain null in all 24 aspect ×
IC × metric cells at the attainable floor 1/(n+1).

This closes the review's **external validation axis** item ("another species, a
later untouched time interval, or both"). It does **not** close the separate
model-selection P0, which asks for choices to be frozen on a development
interval and evaluated once on an untouched one. Mouse was never tuned on, but
the choices applied to it were still selected using human 2021→2026. The nested
human split (205 → 215 → current) is what addresses that, and is in progress.

**It is not a claim of general superiority.** Mouse is a different organism with
shallower curation (9,136 proteins with experimental MF terms against human's
15,260), so "better on mouse" is not "better in general". The defensible claim
is the narrow one: the method does not depend on having been tuned on its
evaluation data.

Metrics: `validation/mouse/temporal_benchmark_metrics.tsv`,
`validation/mouse/temporal_benchmark_permutation_null.tsv`. Full provenance —
the exact releases, the release-matching trap, the commands and the code version
— in `validation/mouse/PROVENANCE.md`.

## What each restored method piece contributes

Both were part of the original dcGO Predictor and had been dropped; restoring them
was worth more than any parameter ablation. They are **complementary**:

- **Per-target p-score** (now the default transfer) is the broad protein-centric
  lever — it lifts dcGO's own F_max/AUPRC across the board relative to `max`
  transfer, and puts dcGO above naive at IC≥0 on BP and CC **on F_max**. It does
  not put dcGO above naive on AUPRC, which naive still wins at IC≥0 in all three
  aspects; see the AUPRC note in `VALIDATION_PLAN.md` §2.
- **Relative (parental-background) inference** barely moves protein-centric F_max
  but earns its keep on the **domain-centric** metric: scored directly against
  InterPro2GO it raises association precision **0.218 → 0.253** (recall 0.63 →
  0.43, halving the set by pruning generic, parent-driven associations).

## Caveats (not oversold)

- MF's face-value gap to naive is real but an artifact of `protein binding`
  dominating the curated truth; a temporal CAFA benchmark rewards recovery of the
  popularity-weighted curation frontier.
- The domain-centric precision is a *lower bound* vs the incomplete, *current*
  InterPro2GO — good for comparing configs, not an absolute precision. A dated
  2021 InterPro2GO would make it fully temporal (planned).
- Information content is the marginal `−log2 P(t)` approximation, not full
  information-accretion. Human-only, single split.

## Reproduce

Assumes the human InterPro subset is already extracted
(`data/interim/protein2ipr_human.dat.gz`; see [README](README.md) Quick Start).

```bash
# 1. Fetch the 2021 t0 GOA snapshot
uv run python scripts/download_data.py --goa-archive 205

# 2. Stage t0 inputs under the --species name run_dcgo_human.py expects.
#    Domains are fixed in time (only GOA moves), so symlink the current subset.
ln -sf "$PWD/data/raw/goa_archive/goa_human.gaf.205.gz" \
       data/raw/goa_annotations/goa_human_t0_2021.gaf.gz
ln -sf "$PWD/data/interim/protein2ipr_human.dat.gz" \
       data/interim/protein2ipr_human_t0_2021.dat.gz

# 3. Train domain→GO associations on t0
uv run python run_dcgo_human.py --species human_t0_2021 --num-cores 32 \
    --output-dir results_t0_2021

# 4. Protein-centric temporal benchmark (t1 = current GOA)
uv run python validation/temporal_benchmark.py \
    --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \
    --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
    --predictions results_t0_2021/domain_go_associations_significant.tsv \
    --min-ic 0 --min-ic 2 --min-ic 4 --min-ic 6

# 5. Domain-centric eval — build the relative-filtered set, then compare both
uv run python validation/apply_relative_inference.py \
    --predictions results_t0_2021/domain_go_associations_significant.tsv \
    --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \
    --output results_t0_2021/domain_go_associations_relative.tsv
uv run python validation/domain_centric_eval.py \
    --predictions base=results_t0_2021/domain_go_associations_significant.tsv \
    --predictions relative=results_t0_2021/domain_go_associations_relative.tsv
```
