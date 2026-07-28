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
`surprise, domain, domain_names, term, term_name, n_feature, n_both,
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

## Caveats

* The emergence test conditions on associations that already passed the dcGO FDR
  filter, so it is a *re-ranking* of a selected set, not an unbiased screen of
  all domain combinations.
* Rates are estimated on the same proteins used to call the original
  associations; the score measures internal consistency of the evidence, not
  out-of-sample predictive performance. A temporal or held-out evaluation of the
  top-ranked hypotheses is the natural follow-up
  (see `VALIDATION_PLAN.md`).
* `--max-overlap` is a heuristic threshold on a continuous quantity; the overlap
  value is reported for every row, so a stricter or looser cut needs no re-run.
* Weights (0.1/0.3/0.6/1.0) for the novelty statuses are a prioritisation
  convention, not an inference. The statistically meaningful quantity is
  `q_emergence`.
