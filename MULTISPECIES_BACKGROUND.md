# The multi-species background

*Run 2026-08-06. Implements the "Expand the background dimension" item queued in
`TODO.md` on the same day.*

Every dcGO-2.0 run before this one was human-only: the Fisher universe was the
18,908 human proteins carrying both a GOA annotation and an InterPro domain.
This run widens the background to **every species GOA annotates**, and reports
what that buys and what it costs.

The short version:

* **It predicts human function better.** On the held-out 2021 → 2026 split, with
  the evaluation set held fixed and only the training universe changed, the
  all-species background wins **8 of 9 F_max cells and 9 of 9 AUPRC cells** (§5).
  This is the criterion that matters most, it is scored against future
  annotation rather than another method's opinion, and the all-species arm won
  it while carrying a handicap.
* **It fixes the defect that motivated the work.** Agreement with the published
  (all-species) dcGO improves 2.3–3.0×, because 566 SCOP superfamilies that
  human data cannot test become testable (§4).
* **But it is not free, and it is not clean.** Three quarters of the annotation
  it adds is inferred rather than observed, and over half of *that* was inferred
  from human (§2) — so the `manual` universe is substantially human function
  copied outward. Support is inflated ~2.44× by orthology, with half the
  associations standing on three or fewer ortholog groups. Precision against the
  published reference falls out of the band the acceptance criteria required.
* **The one thing it was most wanted for is untested.** The emergent n = 2–8
  domain combinations cannot be evaluated at this universe size (§1).

Scored strictly against `TODO.md`'s criteria: two met, one met decisively, one
half met, one unreachable.

---

## 1. What "all species as background" means here

| | human-only baseline | all-species |
|---|---|---|
| annotation source | `goa_human.gaf.gz` | `goa_uniprot_all.gaf.gz` (2026-07-29, 11.7 GB) |
| annotated proteins | 19,089 | 1,503,592 |
| Fisher universe (∩ domains) | 18,908 | **1,464,355** (77×) |
| distinct taxa | 1 | 9,074 |
| InterPro domains testable | 19,449 | 40,141 |
| GO terms testable | 16,389 | 28,112 |

The pipeline needed no new code path. The species-parameterised inputs already
support this: the universe is dropped in as `--species allspecies`, reading
`data/raw/goa_annotations/goa_allspecies.gaf.gz` and
`data/interim/protein2ipr_allspecies.dat.gz` like any other organism.

Two build details matter for reproducing it:

* The all-species GAF is **pre-filtered to non-IEA** before it is dropped in.
  This is lossless for the runs below — the pipeline's `manual` preset is a
  strict subset of non-IEA (it also drops `HTP/HDA/HMP/HGI/HEP`) and
  `experimental` is a subset of that — and it takes the file from 11.7 GB to
  258 MB.
* `goa_uniprot_all` also carries ComplexPortal (`CPX-`) and RNAcentral (`URS`)
  identifiers and 4,716 isoform-form accessions. None of them appear in
  `protein2ipr`, so all are dropped at the domain intersection: 1,503,592 →
  1,464,355 proteins (2.6%). This is the same behaviour as the human path, not a
  new approximation.

One code change was needed, in `config/settings.py`: `goa_url_for` composed a
per-species URL for every input, which for this universe would have written a
plausible-looking `.../ALLSPECIES/goa_allspecies.gaf.gz` into every run manifest
as the input's origin. That URL does not exist. It now resolves the all-species
aliases to EBI's actual cross-organism release.

### Supra-domains are out of reach at this scale

The Fisher engine tests the **dense** domain × term product and materialises one
int32 2×2 table per test. Measured by `scratch_allspecies/07_probe_scale.py`:

| configuration | domains | tests | peak memory |
|---|---|---|---|
| all-species, single domains | 40,141 | 1.13 B | ~34 GB |
| all-species, single + supra | 464,490 | **13.1 B** | **~389 GB** |

389 GB is not impossible on this machine (754 GB) but the Benjamini–Hochberg
sort needs a further ~104 GB of index on top, and this is a shared machine.
**Every all-species run below is single-domain only**, and the human baselines
were re-run single-domain to match. That is not a free choice: the emergent
domain combinations in `SURPRISE_SCORE.md` are exactly the n = 2–8 cases that a
wider background was supposed to rescue, so the acceptance criterion about the
emergent-combination support distribution **cannot be evaluated from this run**.

---

## 2. The three traps, measured

`TODO.md` recorded three traps before any of this ran, on the grounds that each
can invalidate the result silently. All three are real, and two are larger than
the entry assumed.

### Trap 1 — annotation-transfer circularity is the big one

Of the 6,758,017 non-IEA annotations in the all-species universe:

| class | annotations | share |
|---|---|---|
| projected (`ISS ISO ISA ISM IGC IBA IBD IKR IRD RCA`) | 5,123,562 | **75.8%** |
| experimental (`EXP IDA IPI IMP IGI IEP`) | 1,260,457 | 18.7% |

For comparison, the human-only universe is ~16% projected. Widening the
background makes the projected share roughly **five times worse in proportion**.
`IBA` alone — phylogenetic inference, which assigns a term across an ortholog
family by construction — is the single largest evidence code.

Worse, the projection largely points back at us. Of the 5,018,650 projected
annotations on **non-human** proteins, **2,771,461 (55.2%)** name a human
protein in their With/From field. Over half of the non-human annotation this
background adds was itself derived from human annotation.

So a `manual` all-species background is substantially human function, copied
outward and counted again. This is why both evidence policies were run, and why
the `experimental` run is the one that supports a non-circular claim.

Restricting to experimental evidence costs an order of magnitude of universe:
**1,464,355 → 142,653 proteins**. That is still 8.8× the human experimental
baseline (16,242), so the widening is real even under the strict policy — it is
just far smaller than the headline 77×.

### Trap 2 — phylogenetic non-independence is real but moderate

Fisher counts proteins and assumes each is an independent observation; orthologs
are not. Measured against UniRef50 clusters as an ortholog-group proxy
(`scratch_allspecies/05_non_independence.py`), on both keyed runs:

| | InterPro-keyed | SCOP-keyed |
|---|---|---|
| associations | 535,133 | 68,291 |
| pooled protein support | 22,107,423 | 4,429,322 |
| UniRef50-collapsed support | 9,051,092 | 1,848,624 |
| **overall inflation** | **2.44×** | **2.40×** |
| median per-association | 2.00× | 2.00× |
| p90 / p99 / max | 4.0× / 8.8× / 92× | 4.0× / 8.5× / 69× |

The two agree closely on the aggregate, which is reassuring — the inflation is a
property of the universe, not of the domain vocabulary. The bulk of the
distribution is only modestly inflated. The tail is where the problem is, and it
is fatter for the finer-grained InterPro vocabulary:

| support concentration | InterPro | SCOP |
|---|---|---|
| entire support is **1** UniRef50 cluster | **98,937 (18.5%)** | 8,812 (12.9%) |
| ≤ 2 clusters | 207,557 (38.8%) | 19,953 (29.2%) |
| ≤ 3 clusters | **271,176 (50.7%)** | 27,582 (40.4%) |
| entire support is 1 taxon | 42,891 (8.0%) | 5,469 (8.0%) |

So the pooled p-value is optimistic by ~2.4× in aggregate — but **half of the
535,133 headline associations rest on three or fewer ortholog groups**, and
18.5% on a single one, where the optimism is whatever that group's size is (up
to 92-fold). Those associations have one observation behind them, not fifty.

This is the number to quote when the significant count is quoted. **Pooled and
collapsed counts are written per association** to
`scratch_allspecies/out/non_independence_{interpro,ssf}_allspecies.tsv`, so the
inflation is auditable case by case rather than assumed away.

### Trap 3 — "all species" is really about two dozen organisms

Nominally 9,074 taxa. In practice:

* **19 taxa carry 50%** of the annotations
* **79 taxa carry 90%**
* human alone is 11.1% (747,710 annotations)

Ranked: human, mouse, rat, *Arabidopsis*, zebrafish, *Drosophila*, wheat,
*Xenopus*, *S. cerevisiae*, cow. This is a well-curated model-organism panel
with a long thin tail, not a phylogenetically balanced sweep. Any claim about
"conservation across species" from this background is a claim about those
organisms.

---

## 3. What the wider background produces

All runs single-domain, InterPro-keyed, FDR < 0.01. The significant count is
reported next to the size of the hypothesis space that produced it, because a
77× universe grows the count arithmetically and the raw number is not evidence
of anything on its own.

| universe | evidence | proteins | domains | terms | tests | significant | per 1e6 tests |
|---|---|---|---|---|---|---|---|
| human | manual | 18,908 | 19,449 | 16,389 | 318,749,661 | 44,453 | 139.5 |
| human | experimental | 16,242 | 18,962 | 12,966 | 245,861,292 | 18,310 | 74.5 |
| all-species | manual | 1,464,355 | 40,141 | 28,112 | 1,128,443,792 | **535,133** | **474.2** |
| all-species | experimental | 142,653 | 36,460 | 27,504 | 1,002,795,840 | **158,385** | **157.9** |
| all-species | manual, **permuted null** | 1,464,355 | 40,141 | 28,112 | 1,128,443,792 | **0** | **0.0** |

The rate — significant calls per million tests — is the comparable quantity, and
it does rise: 139.5 → 474.2 under `manual` (3.4×) and 74.5 → 157.9 under
`experimental` (2.1×). So the wider background is finding more than a bigger
haystack would explain, under both policies.

But the gap between the two policies is the circularity of §2 showing up in the
results: the rate gain is **3.4× when projected annotation is allowed in and
2.1× when it is not**. Read as excess over the matched human baseline, +2.4×
becomes +1.1× — **more than half the apparent gain from widening the background
does not survive excluding annotation that was inferred rather than observed**,
and per §2 over half of that inferred annotation was inferred from human.

That split is indicative, not a clean decomposition: the two policies differ in
universe size (1,464,355 vs 142,653 proteins) and term count as well as in
evidence, so the comparison attributes to circularity everything that changes
between them.

**The `experimental` row is the one that supports a claim about human function.**
The `manual` row is the one comparable to every previous dcGO-2.0 run.

### The wider universe is calibrated

A 77× universe changes what "random" means, so the permutation control was
re-run rather than inherited: term sets shuffled across all 1,464,355 proteins
(seed 7), which preserves every marginal and the whole hypothesis space and
destroys only the domain↔term link.

It returns **0 significant associations out of 1,128,443,792 tests** — against
535,133 on the real data. The false-positive floor for this universe is at or
below one in a billion, so the significant count is not an artefact of the
hypothesis space having grown. This is the cleanest of the criteria: it passes
outright, and it is what licenses reading the rate column above at all.

Note what it does *not* license. Permutation breaks the domain↔term link but
keeps each protein a separate row, so it cannot detect the phylogenetic
non-independence of §2 — fifty orthologs are still fifty rows under the null as
they are in the real data. **The null being clean and the support being inflated
2.44× are both true, and neither answers the other.**

## 4. Does it help? The published-dcGO comparison (§3)

This is the test the widening was for. `VALIDATION_PLAN.md` §3 compares against
the published dcGO, which is all-species and 2016; our human-only recall against
it was uninterpretable precisely because **69.4% of the pairs they call that we
miss have zero co-occurring human proteins** — unreachable at any threshold. An
all-species universe removes that structural blocker.

Both runs below are SCOP-superfamily-keyed (`--domain-key ssf`), single-domain,
`manual` evidence — the published dcGO's own domain universe.

| | human-only | all-species | change |
|---|---|---|---|
| SCOP superfamilies testable | 911 | **1,477** | +566 |
| shared domains / GO terms | 904 / 9,999 | 1,266 / 14,392 | wider |
| **recall** (theirs found by us) | 0.0739 | **0.2610** | **3.5×** |
| **precision** (ours also called by them) | 0.5259 | **0.2984** | −43% |
| **Jaccard** (overall agreement) | 0.0693 | **0.1617** | **2.3×** |

Per aspect:

| aspect | human P / R | all-species P / R |
|---|---|---|
| molecular function | 0.531 / 0.139 | 0.347 / 0.370 |
| biological process | 0.487 / 0.048 | 0.275 / 0.217 |
| cellular component | 0.611 / 0.112 | 0.327 / 0.353 |

**Read this carefully.** Overall agreement with the published method more than
doubles, and it does so by fixing exactly the defect that motivated the work:
566 superfamilies that human data cannot test become testable, and recall
triples in every aspect. That is a genuine win and it is the strongest evidence
so far that the background dimension was the right one to push.

But precision falls from 0.53 to 0.30, well outside the 0.54–0.63 band the
acceptance criterion required it to hold.

### Most of the precision collapse is the circularity, not the widening

Repeating the whole comparison under `--evidence-filter experimental` — the
policy that admits no projected annotation at all — separates the two effects:

| | human | all-species | change |
|---|---|---|---|
| **manual** precision | 0.5259 | 0.2984 | **−43%** |
| **experimental** precision | 0.5644 | **0.4381** | **−22%** |
| **manual** recall | 0.0739 | 0.2610 | 3.5× |
| **experimental** recall | 0.0427 | **0.1482** | **3.5×** |
| **manual** Jaccard | 0.0693 | 0.1617 | 2.3× |
| **experimental** Jaccard | 0.0414 | **0.1245** | **3.0×** |

Under the non-circular policy the recall gain is **undiminished** (still 3.5×),
the precision loss is **halved** (−22% rather than −43%), and overall agreement
improves **more** (3.0× rather than 2.3×). In other words, the projected
annotation that dominates the `manual` universe was buying recall we already had
and paying for it in precision.

**The experimental all-species run is the better experiment on every axis here,
not merely the more defensible one.** Its precision, 0.438, still sits below the
0.54–0.63 band the criterion demanded, so the criterion is still not met — but
it is a near miss under a stricter policy rather than a collapse.

The residual gap deserves the standing caveat: we are 2026 GOA and they are
2016, so associations we call that they do not include genuinely newer
annotation, and "precision" against a decade-old reference penalises that. The
note in `CLAUDE.md` that *recall* against them is not interpretable now applies
to *precision* too, in the opposite direction. What survives without
qualification is the Jaccard: overall agreement with an all-species method
improves 2.3–3.0× when we stop being human-only.

---

## 5. The held-out test: the wider background predicts human function better

This is the criterion that matters most, because it is the only one scored
against *future* annotation rather than against another method's opinion.

The design isolates the background: **only the training universe changes.** One
arm trains on human 2021, the other on all-species 2021 (release 205, the same
release the human t0 file comes from). Both are scored on the identical human
t0 (2021) → t1 (2026) split from `VALIDATION_PLAN.md` §2 — the naive baseline
comes out byte-identical in both runs, which is the check that the evaluation
really is held fixed.

The criterion was that held-out enrichment must not **fall**. It does not fall;
it rises almost everywhere.

| aspect | IC | F_max human | F_max all-species | AUPRC human | AUPRC all-species |
|---|---|---|---|---|---|
| MF | ≥0 | 0.3668 | **0.4486** | 0.1944 | **0.3054** |
| MF | ≥2 | 0.3660 | 0.3522 *(−)* | 0.2208 | **0.2641** |
| MF | ≥4 | 0.3312 | 0.3319 | 0.2084 | **0.2447** |
| BP | ≥0 | 0.2458 | **0.3315** | 0.1379 | **0.2330** |
| BP | ≥2 | 0.1662 | **0.2459** | 0.0638 | **0.1267** |
| BP | ≥4 | 0.1129 | **0.1802** | 0.0303 | **0.0883** |
| CC | ≥0 | 0.3845 | **0.4634** | 0.2399 | **0.3628** |
| CC | ≥2 | 0.2442 | **0.2820** | 0.0722 | **0.1357** |
| CC | ≥4 | 0.1432 | **0.1923** | 0.0309 | **0.0801** |

**All-species wins 8/9 cells on F_max and 9/9 on AUPRC.** The single F_max loss
(MF IC≥2, −0.014) is the only cell that moves the wrong way in either metric.
AUPRC roughly doubles or triples in the high-IC cells — BP IC≥4 goes 0.030 →
0.088, CC IC≥4 goes 0.031 → 0.080 — which is where a function-prediction method
earns its keep.

And this is **conservative**, for the reason in the caveat below: the
all-species arm was handicapped by losing 31.9% of its 2021 universe, against
the human arm's 1.9%. It won anyway.

Both arms clear the random-domain permutation null in all 18 cells at the
attainable p-floor (p = 0.048 with 20 permutations).

### It partly closes a documented limitation

`CLAUDE.md` records that dcGO "**loses to naive on AUPRC at IC≥0 in all
aspects**". Comparing against the naive baseline (identical in both runs):

| cell | naive | dcGO human | dcGO all-species |
|---|---|---|---|
| MF IC≥0 | 0.3251 | 0.1944 ✗ | 0.3054 ✗ |
| BP IC≥0 | 0.3135 | 0.1379 ✗ | 0.2330 ✗ |
| CC IC≥0 | 0.5130 | 0.2399 ✗ | 0.3628 ✗ |
| CC IC≥2 | 0.1169 | 0.0722 ✗ | **0.1357 ✓** |
| CC IC≥4 | 0.0375 | 0.0309 ✗ | **0.0801 ✓** |

The all-species background takes dcGO from beating naive in 4/9 AUPRC cells to
**6/9**, newly winning both informative cellular-component cells. **The IC≥0
limitation itself survives** — naive still wins there in all three aspects,
though the gap narrows in every one. That line in `CLAUDE.md` should be kept and
qualified, not retired.

### A caveat specific to the held-out arm

The t0 runs reuse the **current** `protein2ipr` rather than a 2021 InterPro
release. That is the existing convention — `protein2ipr_human_t0_2021.dat.gz` is
a symlink to the current file — and for human it costs almost nothing: the human
t0 universe loses 353 of 18,735 proteins (1.9%) to accessions with no current
domain data.

For all species it costs a great deal. The 2021 all-species t0 universe loses
**416,635 of 1,306,108 proteins (31.9%)**, because a third of the accessions
annotated in 2021 have since been merged or deleted out of UniProt and so are
absent from today's `protein2ipr`.

The comparison is still meaningful — the surviving t0 universe is 889,473
proteins against the human arm's 18,382, still 48× — but the attrition is
**16× heavier on the all-species arm**, and it is one-directional: it can only
weaken that arm. So a held-out win here would be conservative, while a held-out
loss would be partly confounded with this attrition and should not be read as a
clean negative. Removing the confound needs an archived 2021 `protein2ipr`, not
a re-run.

---

## 6. Against the acceptance criteria

`TODO.md` deliberately set criteria that a bigger universe cannot satisfy by
arithmetic. Scored honestly:

| criterion | verdict |
|---|---|
| held-out enrichment does not fall vs human-only, on the same human evaluation set | **met, decisively** — all-species wins 8/9 F_max and 9/9 AUPRC cells (§5), and does so despite a handicap (caveat above) |
| emergent-combination support shifts up, no more redundant signatures | **cannot be evaluated** — supra-domains do not fit in memory at this universe size (§1) |
| recall vs published dcGO rises materially from 0.069 **while precision holds at 0.54–0.63** | **half met** — recall 0.074 → 0.261; precision 0.526 → 0.298, criterion fails |
| non-independence measured, pooled vs collapsed side by side | **met** — 2.44× overall, 18.5% single-cluster (§2) |
| permutation null re-run on the wider universe | **met** — 0 significant of 1.13 B tests (§3) |

Three criteria met (one decisively), one half met, one unreachable from this
configuration.

**The widening is established as an improvement in predictive power**, on the
strongest evidence available in this project: a held-out split, scored against
annotation that did not exist when the model was trained, with the evaluation
set held byte-identical and only the training universe varied. 8/9 and 9/9 is
not a marginal result, and the losing arm had a 16× lighter handicap.

Three things keep that from being the whole story:

* **The `manual` universe is partly circular** (§2). The held-out test is the
  one result circularity cannot fake — future human annotation is not in the
  2021 training data under any evidence code — but the *size* of the win is
  still measured on a background that is 75.8% projected. The
  experimental-evidence held-out arm was not run and is the obvious next step.
* **Support is inflated ~2.44× by orthology**, with half the associations
  standing on three or fewer ortholog groups (§2). This does not touch the
  held-out result, which never uses the p-values as probabilities, but it does
  mean the 535,133 count and its FDR should not be quoted unqualified.
* **The emergent domain combinations remain untested** (§1) — the case a wider
  background was most wanted for.

---

## 7. Reproducing

Build scripts are in `scratch_allspecies/`, numbered in run order. The heavy
steps are the two streaming scans (the 11.7 GB GAF and the 13 GB
`protein2ipr.dat.gz`).

```bash
bash scratch_allspecies/01_build_gaf.sh        # non-IEA all-species GAF
bash scratch_allspecies/06_build_universe.sh   # protein list + protein2ipr subset
uv run python scratch_allspecies/03_characterise_gaf.py   # §2 traps 1 and 3
uv run python scratch_allspecies/07_probe_scale.py        # §1 memory arithmetic
bash scratch_allspecies/08_run_allspecies.sh   # the three InterPro-keyed runs
bash scratch_allspecies/09_run_ssf.sh          # the two SCOP-keyed runs for §3
uv run python scratch_allspecies/05_non_independence.py \
    --associations results_ssf_allspecies_manual/domain_ssf_go_associations_significant.tsv \
    --domain-key ssf --out scratch_allspecies/out/non_independence_ssf_allspecies.tsv
```

Every run writes its own `run_manifest_go.json` with input SHA-256s, so the
universe a number came from is recoverable from the number.
