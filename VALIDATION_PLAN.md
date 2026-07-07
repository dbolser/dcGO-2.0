# Validation & Benchmarking Plan

This document defines what must be done to turn dcGO-2.0 from a working codebase
into a defensible methods contribution. It is the "(b)" companion to the
engineering cleanup: the code runs; this is about showing the results are
*reliable*.

## 0. Why the current validation is not enough

`validation/validate_results.py` compares predictions to **InterPro2GO** with a
plain set-overlap: any predicted `(domain, GO)` pair not present in InterPro2GO
is counted as a false positive. The committed metrics
(`validation/performance_metrics.tsv`) therefore show ~3–4% "precision" and
~96% "novel" predictions.

This number is **largely an artifact of the comparison, not a measured error
rate**, for three reasons:

1. **InterPro2GO is deliberately incomplete.** It is a conservative, manually
   curated mapping. It is a *positive-only, partial* set — absence of a pair
   does not mean the pair is wrong. So it can bound **recall**, but it cannot be
   used to compute **precision**.
2. **No propagation alignment.** Predictions may be propagated up the GO DAG
   (True Path Rule) while InterPro2GO entries are at specific levels, or vice
   versa. Without propagating both sides to a common closure, true matches are
   miscounted as misses.
3. **Duplicated / plateaued rows** in the metrics file suggest the sweep itself
   is buggy (identical rows for several thresholds), so even the recall trend is
   not trustworthy as-is.

**Conclusion:** the first task is not "improve the score" — it is to build a
validation design where the numbers *mean* something.

---

## 1. Fix and reframe the InterPro2GO comparison  *(quick win)*

Treat InterPro2GO as an **incomplete positive reference**, and report only what
it can legitimately support.

- [ ] De-duplicate the threshold sweep; verify each threshold changes the
      prediction set (the current file has identical rows — a bug).
- [ ] **Propagate both** predictions and the reference to their GO ancestor
      closure before intersecting (use `ontology_processor` + `go-basic.obo`).
- [ ] Report **recall / coverage of InterPro2GO** as the headline (what fraction
      of curated pairs we recover), and explicitly label "novel" pairs as
      *not-in-reference* rather than *false*.
- [ ] Restrict the comparison to the **domains actually present in both** sets
      (a domain absent from InterPro2GO contributes only noise).

**Acceptance:** a corrected `performance_metrics.tsv` where recall of
InterPro2GO rises monotonically as thresholds loosen, and "novel" is framed as
candidate discoveries, not errors.

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

---

## 3. Comparison to the original dcGO  *(reproducibility)*

A method claiming to implement dcGO must relate its output to the published one.

- [ ] Download the original dcGO / SUPFAM domain–GO associations.
- [ ] Map identifier spaces (InterPro ↔ SUPERFAMILY/Pfam) as far as possible;
      document coverage of the mapping.
- [ ] Report agreement on the shared domain space and characterize where
      dcGO-2.0 differs (and why — newer GOA, InterPro vs SUPERFAMILY domains,
      supra-domains, etc.).

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
- [ ] **Significance vs. effect size.** With ~300M tests, many pairs reach tiny
      p-values on thin evidence (note the `odds_ratio = inf`, `hyper_score = 100`
      rows). Decide and document the minimum-evidence / effect-size floor
      (`MIN_PROTEINS_PER_ASSOCIATION`, odds-ratio bounds).

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

1. §1 (fix the existing comparison) — days, unblocks honest reporting.
2. §2 (temporal benchmark + baselines) — the core result.
3. §4 (ablation) — reuses the §2 harness.
4. §3 (original-dcGO comparison) and §5–§6 in parallel.

The engineering is sound; §2 is what turns "it runs" into "it works."
