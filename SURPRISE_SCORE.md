# The surprise score: ranking emergent domain-combination predictions

## Why

dcGO's distinctive claim is **emergence**: a *combination* of domains predicts a
term that none of its constituents predicts alone. That is the part of the
method homology transfer and single-domain annotation cannot reproduce, so it is
where the interesting predictions live.

Ranking those by raw association significance does not work. Three failure modes
dominate the top of that list:

1. **Redundant signatures.** InterPro annotates one region with several
   signatures (family + domain + homologous superfamily + site). "GPCR
   rhodopsin-like + 7TM" looks like a two-domain architecture in the
   architecture string but is one region described twice, and it reaches
   spectacular p-values because it is effectively a single very common domain.
2. **Restated curated knowledge.** A combination "predicting" what InterPro2GO
   already maps to one of its constituents is not a discovery.
3. **The support/novelty tension.** Genuinely novel combinations sit on two to
   eight proteins; the well-supported ones tend to be cases 1 and 2.

The surprise score addresses all three, and — importantly — is *statistical
about emergence itself* rather than reusing the original association p-value.

## Definition

For a supra-domain `S = d1,d2[,d3]` associated with term `t`:

```
surprise = -log10(q_emergence) × distinctness × novelty
```

### 1. Emergence — does the whole beat the sum of its parts?

Let `rate(f) = P(t | protein carries f)`. Rates are shrunk toward the term's
background rate with one pseudo-observation,
`(a + α·p_bg) / (n + α)`, so a domain seen in three proteins that all carry `t`
does not get to claim `rate = 1.0` and thereby declare every combination
containing it unremarkable.

The **parts-only expectation** is the largest of:

* the **noisy-OR** over the constituent single domains, `1 - Π(1 - rate(d_i))` —
  the rate expected if each domain contributed independently;
* the best **proper sub-combination** rate (a triple whose contained pair already
  predicts `t` is not surprising);
* the term's **background rate** in the protein universe (nothing that merely
  matches random sampling is evidence, and this keeps the null positive).

The observed count is tested against that expectation with a one-sided binomial
tail, `P(X ≥ a | n, p_expected)`, and Benjamini–Hochberg corrected across all
candidates — a second, separate hypothesis family from the original Fisher tests,
so it needs its own correction.

Support is handled honestly by construction rather than by a hand-set threshold:
two supporting proteins reach significance only when the parts-only expectation
is genuinely tiny.

### 2. Distinctness — is it really two domains?

`1 - overlap`, where `overlap` is the median (over supporting proteins) largest
pairwise overlap between the constituents' matched regions, taken from the
`protein2ipr` coordinates. Redundant signatures for one region score ~0 and are
filtered out (`--max-overlap`, default 0.5); genuinely separate domains score ~1.

### 3. Novelty — do curators already say this?

Against a curated domain→term reference (InterPro2GO, GO only):

| Status | Meaning | Weight |
| --- | --- | ---: |
| `curated` | the term is already mapped to a constituent | 0.1 |
| `implied` | the term is *more general* than curated constituent knowledge | 0.3 |
| `refines` | the term is *more specific* than curated knowledge | 0.6 |
| `novel` | outside the curated closure | 1.0 |
| `no-reference` | the reference says nothing about these constituents | 1.0 |

Only GO has such a reference here; other ontologies score `no-reference`, and the
ranking then rests on emergence and distinctness alone.

Every component is written to the output, so results can be re-ranked or
re-weighted without recomputing anything.

## Usage

```bash
uv run python scripts/rank_surprising_associations.py --ontology go
uv run python scripts/rank_surprising_associations.py --ontology ec
```

Output: `results/domain_<ontology>_surprising.tsv`, ranked, with columns
`rank, surprise, domain, domain_names, term, term_name, n_feature, n_both,
observed_rate, expected_rate, expectation_source, lift, p_emergence,
q_emergence, region_overlap, distinctness, novelty, novelty_status,
uninformative_constituents, dcgo_adj_p_value`.

Useful knobs: `--max-overlap` (redundancy filter), `--min-support`,
`--pseudo-count` (rate shrinkage), `--alpha` (emergence FDR).

## Results on the current human run

### Gene Ontology (18,908 proteins, `results/domain_go_associations_significant.tsv`)

| Stage | Count |
| --- | ---: |
| Supra-domain associations scored | 123,203 |
| Dropped as redundant signatures (overlap > 0.5) | 100,960 (82%) |
| Ranked candidates retained | 22,243 |
| Emergent beyond their parts at FDR ≤ 0.05 | 24 |

Novelty breakdown of the retained candidates: 10,353 `novel`, 6,014
`no-reference`, 3,840 `refines`, 1,716 `curated`, 320 `implied`.

Top of the ranking (all with distinctness 1.0):

| Surprise | Architecture | Predicted term | Support | Lift |
| ---: | --- | --- | --- | ---: |
| 6.7 | SH2 domain + protein kinase-like | non-membrane spanning protein tyrosine kinase activity (GO:0004715) | 22/22 | 3× |
| 3.4 | Ig C1-set + Ig-like fold | immunoglobulin receptor binding (GO:0034987) | 9/9 | 7× |
| 3.2 | PH domain + EF-hand pair | PIP2 phospholipase C activity (GO:0004435) | 7/7 | 11× |
| 2.0 | Tyr kinase catalytic + SAM/pointed | ephrin receptor signaling (GO:0048013) | 12/14 | 4× |
| 1.8 | BTB/POZ + C2H2 zinc finger | DNA-binding transcription repressor activity (GO:0001227) | 17/22 | 3× |
| 1.5 | SH3 + P-loop NTPase | receptor localization to synapse (GO:0097120) | 4/7 | 26× |

These are textbook multi-domain architectures recovered *without* being told
about them — Src-family kinases (SH2 + kinase), phospholipase C (PH + EF-hand),
Eph receptors (kinase + SAM), BTB-ZF repressors, and the postsynaptic
SH3 + P-loop case, whose two constituents are individually the most promiscuous
domains in the whole set and therefore useless on their own. That the ranking is
dominated by *known* architectures is the validation: the score puts real
emergent biology on top, so the lower-ranked `novel` rows are worth reading as
hypotheses.

### Enzyme Commission (3,704 enzymes)

| Stage | Count |
| --- | ---: |
| Supra-domain associations scored | 8,637 |
| Dropped as redundant signatures | 7,236 (84%) |
| Ranked candidates retained | 1,401 |
| Emergent beyond their parts at FDR ≤ 0.05 | 1 |

This is a deflationary result and worth stating plainly. Earlier exploratory
work counted 96 EC associations that were significant for a combination but for
none of its constituents. Under the parts-baseline test almost all of them
disappear, for two reasons visible in the output:

* **A catalytic constituent usually explains the enzyme.** The tankyrase example
  (SAM + PARP catalytic → EC 2.4.2.30) has a parts-only expectation of 0.82,
  because the PARP catalytic domain alone predicts the activity. The
  architecture is real; the *prediction* is not emergent.
* **The remainder lack support.** Ankyrin repeat + SAM → EC 2.4.2.30 gets
  `p = 0.019` on two proteins, which does not survive correction across 8,637
  candidates.

"Significant for the pair, not significant for either part" is a much weaker
criterion than "the pair statistically beats what the parts predict", and the EC
run is where the two diverge most, because the EC universe is enzymes only, so
constituent rates are high. The one survivor is the Ser/Thr kinase
active site + catalytic domain pair for EC 2.7.11.25 (8/11, 6.6× lift).

## Held-out temporal validation (2026-08-04)

Everything above is computed on the proteins that produced the associations, so
it measures internal consistency, not predictive power. `validation/temporal_surprise.py`
is the held-out test: **score the associations dcGO found in 2021, then ask how
often curators added the predicted term by 2026.**

For each t0 association (supra-domain `S`, term `t`): *predictions* are proteins
carrying `S` that had no `t` at t0 (propagated, non-IEA — the same evidence space
the pipeline trains on, so the gate cannot leak a seen label); *hits* are those
annotated `t` (propagated, experimental) at t1; and the *base rate* is the term's
own acquisition rate among all domain-carrying proteins lacking it at t0. The
statistic is `enrichment = hit rate / base rate`, which removes the fact that
popular terms accumulate annotations regardless of any prediction.

Pool: 22,376 t0 candidates, of which 10,136 make at least one prediction —
170,416 protein-term predictions in total.

Every number below comes from `validation/temporal_surprise_metrics.tsv`, which
is regenerated by the command at the end of this section. Interval types are
named explicitly, because on one of these statistics the choice of interval
changes the answer.

### Result 1 — the associations are strongly predictive (12.5×)

Percentile bootstrap over associations, 10,000 resamples:

| Stratum | assoc. | predictions | hits | hit rate | expected | enrichment (95% percentile CI) |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| **all candidates** | 10,136 | 170,416 | 2,181 | 1.3% | 0.1% | **12.49× [10.87, 14.39]** |
| surprise top-100 | 100 | 481 | 12 | 2.5% | 0.2% | 16.47× [5.73, 26.04] |
| surprise top-500 | 500 | 4,440 | 56 | 1.3% | 0.1% | 14.78× [10.68, 21.00] |
| dcGO-q top-100 | 100 | 7,293 | 109 | 1.5% | 0.1% | 10.30× [5.79, 21.25] |
| dcGO-q top-500 | 500 | 18,651 | 393 | 2.1% | 0.2% | 13.23× [9.39, 19.53] |

2,181 predictions came true where the terms' own rates predict ~175. This is
out-of-sample evidence that the supra-domain associations anticipate curation,
and it is the strongest claim in this document.

These strata are *fixed* sets: nothing is selected inside a resample, so the
bootstrap here is an ordinary one and the percentile interval is fine. That is
not true of Result 2.

### Result 2 — the ranking comparison is mostly unresolvable, and surprise is still not shown to beat plain significance

Top-K is not a fair comparison: surprise favours tight architectures (≈5
predictions each) while significance favours common domains (≈73), so the same K
exposes 15× more predictions for one than the other. Matching on **prediction
budget** — "given capacity to check N predictions, which ranking finds more that
come true?" — is the right idea. Executing it on this data is not possible, for
a reason that only became visible once the bootstrap was diagnosed properly.

#### The blocking problem: the surprise score does not rank most of the pool

Of the 10,136 evaluated associations, **213 have a surprise score above zero and
9,923 are tied at exactly 0.000** (their emergence test returns `q = 1`, so
`-log10(q) = 0` and the product is zero regardless of distinctness and novelty).
A budget slice therefore consists of the genuinely ranked head plus however much
of that 9,923-way tie the input file happened to list first:

| Budget | surprise slice | of which decided by the tie-break | dcGO-q slice | of which tied |
| --- | ---: | ---: | ---: | ---: |
| 2,000 predictions | 363 assoc. | 150 (41%) | 17 assoc. | 0 |
| 10,000 predictions | 580 assoc. | 367 (63%) | 169 assoc. | 0 |
| 40,000 predictions | 1,391 assoc. | 1,178 (85%) | 1,630 assoc. | 0 |

The 213 genuinely ranked associations make only 940 predictions between them, at
an enrichment of 15.3×. Everything beyond a 940-prediction budget is the tie
block. So at the 10,000 budget, 9,372 of the surprise arm's 10,312 predictions
come from associations the score does not order at all.

Re-breaking that tie at random (2,000 shuffles, same shuffle applied to both
rankings so the comparison stays paired) shows what that is worth:

| Budget | observed difference | difference under random tie-breaks | observed sits at |
| --- | ---: | --- | ---: |
| 2,000 predictions | +10.19 | +9.82 [+0.57, +20.00] | 54th percentile |
| 10,000 predictions | +7.96 | **+0.12 [−5.53, +7.57]** | **98th percentile** |
| 40,000 predictions | +0.81 | +1.77 [−0.94, +5.71] | 27th percentile |

The +7.96 at the 10,000 budget — the number the previous version of this
document treated as a real if non-significant edge — is a property of the order
of the input file. Under a random ordering of associations the score says
nothing about, the same statistic is +0.12.

#### The three intervals, side by side

Two paired resampling designs are reported, both driven by one draw per
replicate so that the two rankings always see identical data. **reselect** is
the original design: resample the pool, then re-fill the budget from scratch, so
the composition of the top of each ranking is resampled too. **fixed** conditions
on the selection: the two budget slices are the ones the observed ranking
produced, and a resample only re-weights associations (shared members move
together). Neither dominates — `fixed` understates selection uncertainty,
`reselect` is barely a bootstrap when few associations fit the budget.

10,000 resamples, seed 0, BCa acceleration from a full 10,136-point jackknife:

| Budget | design | difference | percentile | basic | BCa | favours surprise |
| --- | --- | ---: | --- | --- | --- | ---: |
| 2,000 | reselect | +10.19 | [−83.81, +18.50] | [+1.88, +104.20] | [−42.78, +21.50] | 75% |
| 2,000 | fixed | +10.19 | [−3.29, +18.84] | [+1.55, +23.68] | [−0.71, +20.02] | 95% |
| 10,000 | reselect | +7.96 | [−8.05, +17.57] | [−1.65, +23.96] | [−4.57, +21.17] | 78% |
| 10,000 | fixed | +7.96 | [−2.33, +18.11] | [−2.19, +18.25] | [−2.30, +18.15] | 94% |
| 40,000 | reselect | +0.81 | [−3.00, +7.25] | [−5.62, +4.63] | [−4.63, +5.28] | 72% |
| 40,000 | fixed | +0.81 | [−2.99, +5.93] | [−4.30, +4.61] | [−3.87, +5.03] | 67% |

Read the 2,000/reselect row first. Percentile says the difference could be −84;
basic — the exact same replicates, reflected about the point estimate — says it
is at least +1.9 and therefore *significant in favour of surprise*. Two textbook
intervals from one set of numbers, disagreeing about the sign. That happens
because the replicate distribution is violently left-skewed (sd 24.9 against a
point estimate of +10.2, skew −1.93): only 17 dcGO associations fit a
2,000-prediction budget, the largest of them carries 333 of its 2,031
predictions, and one duplicated copy of that association swings the dcGO arm by
tens of enrichment units. BCa splits the difference and should not be believed
either. **We report no interval at this budget.**

#### Verdict per budget

| Budget | resolvable? | why |
| --- | --- | --- |
| 2,000 | **no** | only 17 dcGO associations fit; percentile and basic disagree about the sign |
| 10,000 | **no** | 63% of the surprise slice is an arbitrary tie-break, and the observed value is the 98th percentile of the tie-break spread |
| 40,000 | no (and no effect anyway) | 85% of the surprise slice is an arbitrary tie-break; every interval, on either design, comfortably spans zero |

So the conclusion of the previous version stands — **the surprise score is not
demonstrated to rank better than the dcGO q-value** — but the reason is stronger
and less flattering than "the intervals span zero". At the budgets a curator
would care about, the surprise score is not ranking: it has 213 opinions and
9,923 abstentions, and the comparison is measuring the abstentions.

What the score demonstrably does is still real, and unaffected by any of this:

* Its **213 scored associations** hit at 15.3× the base rate on 940 predictions
  (16 hits against 1.04 expected) — better than the pool's 12.5×, on a set small
  enough that the interval is wide.
* It applies two filters a q-value cannot: redundant-signature removal and
  novelty against curated knowledge. Those are interpretability work, not
  ranking work, and they are why the top of the list is readable at all.

The actionable consequence is a design fix, not a statistical one: the emergence
FDR needs to produce a *graded* score over the whole candidate set (e.g. keep the
raw `p_emergence`, or rank the `q = 1` tail by lift) before any budget-matched
comparison against the q-value can mean anything.

### Two problems this section used to have

Both were found by cross-checking this document against the committed artifact,
and both are worth recording so they are not reintroduced.

**1. The prose did not match the committed data.** The intervals quoted here were
not the ones in `validation/temporal_surprise_metrics.tsv` (the point estimates
matched; the intervals and the favouring fractions did not). Replaying the
committed `validation/temporal_surprise_associations.tsv` at the documented
defaults reproduces the *prose* exactly, so the prose was the honest half; the
committed metrics file was stale, and no combination of `--seed` (0–11) and
`--bootstrap` (500/1,000/2,000) reproduces its comparison block from any
committed input. It was written by a state of the code or the inputs that no
longer exists.

That was possible because `validation/temporal_surprise_associations.tsv` — the
file the documented `--replay` command re-analyses — was written lossily:
`base_rate` to 5 decimal places, which is one or two significant figures for the
smallest acquisition rates (1/18,908 ≈ 5.3e-5), and `dcgo_score` to 2. Replay
therefore could not reproduce a run even in principle, so the two artifacts could
drift apart unnoticed. It now writes 12 significant figures, and the replay
command below is bit-identical to the run.

**2. The point estimate fell outside its own confidence interval.** The stale
artifact reported a difference of +7.96 at the 10,000 budget with a percentile
interval of [−8.93, +7.06] — the statistic above its own upper bound. That exact
run is not reproducible, but the mechanism that pushes the distribution away from
the point estimate is, and it is measured rather than guessed:

| Budget | tie-break rule inside the resample | bootstrap median | z0 (bias correction) |
| --- | --- | ---: | ---: |
| 2,000 | matched to the point estimate | +8.54 | +0.30 |
| 2,000 | random (the old behaviour) | +7.89 | +0.28 |
| 10,000 | matched to the point estimate | +6.34 | +0.25 |
| 10,000 | random (the old behaviour) | **+0.64** | **+1.48** |
| 40,000 | matched to the point estimate | +1.45 | −0.26 |
| 40,000 | random (the old behaviour) | +1.59 | −0.42 |

`sorted(..., reverse=True)` is stable in the order of the list it is given, so
re-sorting the *drawn* list inside each resample ordered tied associations at
random, while the point estimate used the input order. With 98% of the pool tied
that is not a nuisance term — it is the dominant one, and it means the point
estimate and its replicates were estimating different quantities. At the 10,000
budget it moved the bootstrap median from +6.34 to +0.64 and pushed the bias
correction to z0 = +1.48 — the observed value at the 93rd percentile of its own
resampling distribution, most of the way to falling outside a 95% interval, which
is what the stale artifact recorded. The tie order is now fixed and explicit, and
tie-break variability is measured on purpose instead.

A smaller, genuine bias survives the fix: a with-replacement resample contains
only ~63% of the distinct associations, so filling a *fixed prediction budget*
forces it deeper into the ranking than the observed slice, picking up worse
associations. That is the residual z0 ≈ +0.25 to +0.30 on the `reselect` design,
and it is why the `fixed` design (z0 ≈ 0.02) is also reported.

`validation/temporal_surprise.py` now refuses to quote an interval when the
bootstrap cannot support one — when the observed statistic falls outside its own
percentile interval, when |z0| > 0.5, when fewer than 30 associations fit the
budget, or when over half a slice is an arbitrary tie-break. Every comparison row
in the metrics file carries `trustworthy` (can this bootstrap support an
interval?), `resolvable` (is the comparison answerable at all at this budget?),
a `verdict` and a `note`.

### A structural tension the test exposed

`surprise top-25` scored **0 hits on 117 predictions** — but the base rate
predicts 0.09 hits there, so this cell is uninformative rather than negative.
The reason it is so small is the interesting part: emergence *requires* that a
combination's carriers are nearly all annotated with the term, which by
construction leaves almost nothing unannotated to predict. The most emergent
associations are therefore the least testable, and the score's sharpest end
cannot be validated this way. Any future version should weigh emergence against
the number of standing predictions it leaves.

Reproduce with:

```bash
uv run python scripts/rank_surprising_associations.py --ontology go \
    --gaf data/raw/goa_archive/goa_human.gaf.205.gz \
    --associations results_t0_2021/domain_go_associations_significant.tsv \
    --output results_t0_2021/domain_go_surprising.tsv
# writes both committed artifacts; ~7 min
uv run python validation/temporal_surprise.py \
    --bootstrap 10000 --tie-shuffles 2000 --seed 0
# re-analyse without re-parsing the snapshots — bit-identical to the above
uv run python validation/temporal_surprise.py \
    --replay validation/temporal_surprise_associations.tsv \
    --bootstrap 10000 --tie-shuffles 2000 --seed 0
```

1,000 resamples (the default) is enough for the fixed-stratum intervals of
Result 1 but not for the BCa tails of Result 2; the committed artifact uses
10,000.

## Caveats

* The emergence test conditions on associations that already passed the dcGO FDR
  filter, so it is a *re-ranking* of a selected set, not an unbiased screen of
  all domain combinations.
* Rates are estimated on the same proteins used to call the original
  associations, so the score itself measures internal consistency of the
  evidence. The held-out section above is what tests predictive power — and its
  verdict is that the *associations* predict future curation strongly (12.5×)
  while the *ranking* is not demonstrably better than the q-value.
* The score is **not a total order**. On the GO run 9,923 of the 10,136 evaluated
  associations score exactly 0.000, because `q_emergence` saturates at 1 and the
  product collapses. Any comparison that reaches deeper than the ~940 predictions
  those 213 scored associations make is comparing an arbitrary ordering, and the
  budget-matched head-to-head above says so explicitly rather than reporting a
  number. Fixing this — a graded emergence score over the whole candidate set —
  is the precondition for repeating that comparison.
* `--max-overlap` is a heuristic threshold on a continuous quantity; the overlap
  value is reported for every row, so a stricter or looser cut needs no re-run.
* Weights (0.1/0.3/0.6/1.0) for the novelty statuses are a prioritisation
  convention, not an inference. The statistically meaningful quantity is
  `q_emergence`.
