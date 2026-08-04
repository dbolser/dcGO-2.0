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
single near-universal term (`protein binding`, 84.6% of experimental MF
annotations) — and dcGO overtakes it decisively the moment that noise is removed.

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

Benchmark sizes (IC≥0): **BP 324 / MF 418 / CC 572** proteins.

## Protein-centric result — F_max (dcGO / naive / random)

| Aspect | IC floor | dcGO | naive | random | dcGO ÷ random |
|--------|:--------:|-----:|------:|-------:|--------------:|
| **BP** | ≥0 | **0.248** | 0.115 | 0.158 | 1.6× |
| BP | ≥2 | **0.170** | 0.071 | 0.053 | 3.2× |
| BP | ≥4 | **0.115** | 0.031 | 0.019 | 6.1× |
| BP | ≥6 | **0.077** | 0.010 | 0.003 | 24× |
| **MF** | ≥0 | 0.360 | **0.464** | 0.262 | 1.4× |
| MF | ≥2 | **0.365** | 0.053 | 0.088 | 4.2× |
| MF | ≥4 | **0.337** | 0.045 | 0.072 | 4.7× |
| MF | ≥6 | **0.217** | 0.018 | 0.009 | 25× |
| **CC** | ≥0 | **0.380** | 0.343 | 0.291 | 1.3× |
| CC | ≥2 | **0.239** | 0.153 | 0.072 | 3.3× |
| CC | ≥4 | **0.134** | 0.099 | 0.031 | 4.3× |
| CC | ≥6 | **0.124** | 0.044 | 0.015 | 8.1× |

*IC in bits; IC≥2 ⇒ term in ≤25% of proteins, IC≥6 ⇒ ≤1.6%. Bold = winner.*
**Naive is a base-rate mirage:** its F_max collapses toward the random null as
informative terms are required (BP 0.115→0.010, MF 0.464→0.018, CC 0.343→0.044
over IC 0→6). dcGO degrades gracefully because its predictions are specific.

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

## What each restored method piece contributes

Both were part of the original dcGO Predictor and had been dropped; restoring them
was worth more than any parameter ablation. They are **complementary**:

- **Per-target p-score** (now the default transfer) is the broad protein-centric
  lever — it lifts F_max/AUPRC across the board and puts dcGO above naive at IC≥0
  on BP and CC.
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
