# Evidence ledger for `paper/dcgo2_manuscript.md`

Every quantitative claim that appears (or was considered for) the manuscript is
listed here with its provenance. **Rule applied:** a number may appear in the
manuscript only if it has a row here with verdict *supported* or *provisional*
(and provisional numbers must be labelled as such in the text). Anything marked
*unsupported* has been kept out of the manuscript body.

Compiled: 2026-08-04. Repository: `/home/danbolser/Build/dcGO-2.0`,
branch `feature/surprise-score-and-ontology-breadth`, HEAD `a334b0b`.

**Verdict key**

| Verdict | Meaning |
| --- | --- |
| **supported** | Traceable to a committed data file, log, or code line in the repository; reproduces on inspection. |
| **provisional** | Traceable only to a prose document (no machine-readable artifact found), **or** known to be affected by an identified defect, **or** contradicted by a second source in the repository. Must be labelled in the text. |
| **unsupported** | Cannot be traced to any file. **Excluded from the manuscript.** |

A cross-cutting caveat applies to every row: none of the runs below carries a
run manifest pinning input releases and checksums
(`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` P0/P1), so "the 2026 GOA snapshot" and
"the current InterPro release" are not reproducible identifiers today.

---

## A. Datasets and pipeline scale

| ID | Claim | Value | Source (file : locator) | What was actually measured | Verdict |
| --- | --- | --- | --- | --- | --- |
| A1 | Proteins in the t0 (2021) training run | 18,382 | `logs/t0_2021_run.log` (2026-07-09 11:16:47, `main:275` "Dataset prepared") | Proteins present in **both** the GOA release-205 non-IEA annotation set and the human `protein2ipr` subset | supported |
| A2 | Domain features tested in the t0 run | 102,206 | `logs/t0_2021_run.log` : same line | Single domains **plus** supra-domains (contiguous combinations up to length 3) that survive the feature build | supported |
| A3 | Distinct single InterPro domains in the t0 run | 19,230 | `logs/t0_2021_run.log` (`src.sparse_fisher:build_domain_metadata:140`) | Single-domain columns of the protein × feature matrix | supported |
| A4 | GO terms in the t0 run | 16,055 | `logs/t0_2021_run.log` : same line as A1 | Distinct GO terms annotated to at least one protein in the intersection | supported |
| A5 | Fisher tests performed, t0 run | 1,640,917,330 | `logs/t0_2021_run.log` (`main:278` "Total tests") | `n_features × n_terms` (102,206 × 16,055) | supported |
| A6 | Single-domain annotations / supra-domain annotations, t0 run | 230,275 / 405,928 | `logs/t0_2021_run.log` (`main:260`, `main:262`) | (protein, feature) incidences | supported |
| A7 | Significant domain–GO associations, t0 run, FDR < 0.01 | 164,549 | `logs/t0_2021_run.log` (`main:512`); `results_t0_2021/domain_go_associations_significant.tsv` = 164,550 lines incl. header | BH-adjusted Fisher p < 0.01 | supported |
| A8 | Significant domain–GO associations, current-release run, FDR < 0.01 | 165,823 | `results/domain_go_associations_significant.tsv` (165,824 lines incl. header); also `VALIDATION_PLAN.md` §1 | As A7, on the current GOA/InterPro release | supported |
| A9 | Human proteins carrying InterPro domains | 18,908 | `SURPRISE_SCORE.md` L103; `docs/uniprot_ontology_survey.md` L136 | Protein universe for the current-release run | provisional — stated identically in two prose documents; no run log located that emits this figure |
| A10 | GO DAG size used for propagation | 38,245 terms / 71,895 relationships | `logs/bench_A.log` … `logs/bench_primary.log` (`_prepare_graph:143`) | `go-basic.obo` after obsolete-term removal | supported |
| A11 | Swiss-Prot flat file scale | 575,503 entries, of which 20,431 human (`OX NCBI_TaxID=9606`) | `docs/uniprot_ontology_survey.md` L11–12 | Survey of the 2026-07 `uniprot_sprot.dat.gz` | supported |
| A12 | Registered `--ontology` values | 19 | `README.md` L243; `src/ontology_registry.py`; `TODO.md` L38 | Keys in the dispatch table | supported |
| A13 | Test suite at external-review date | 162 tests, 6.63 s | `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` L133 | Reviewer's `pytest` run, 2026-07-15 | provisional — `CLAUDE.md`/`README.md` state "155 tests, ~5 s"; the two disagree and neither is a committed artifact |

---

## B. Domain-centric evaluation against InterPro2GO (§1)

Source file for B1–B6: `validation/performance_metrics.tsv`. Semantics from
`validation/validate_results.py` and `validation/domain_centric_eval.py`.

| ID | Claim | Value | Source (file : locator) | What was actually measured | Verdict |
| --- | --- | --- | --- | --- | --- |
| B1 | InterPro2GO reference coverage at the FDR < 0.01 cutoff | 0.6472 (64.7%) | `validation/performance_metrics.tsv`, row `adj_p≤1.00e-02`, col `reference_coverage` | Fraction of curated (domain, GO) pairs — propagated to ancestor closure, restricted to domains present in **both** sets — that the prediction set recovers. **Recall against an incomplete positive-only reference; not precision, not independent validation.** | supported |
| B2 | Pairs recovered / reference size | 30,673 / 47,393 | same row, cols `recovered`, `n_reference` | Propagated curated pairs on shared domains | supported |
| B3 | "Precision lower bound" at the same cutoff | 0.2231 | same row, col `precision_lower_bound` | Recovered ÷ predicted. A **floor**, because non-recovered pairs are curation-gap candidates, not established errors | supported |
| B4 | Coverage at a tighter cutoff, p ≤ 1e-10 | 0.2929 (13,881 pairs) | `validation/performance_metrics.tsv`, row `p_value≤1.00e-10` | As B1 | supported |
| B5 | The threshold sweep plateaus at the loose end | Rows `p≤1e-6`, `p≤1e-4`, `p≤1e-2`, `adj_p≤1e-2`, `adj_p≤5e-2`, `adj_p≤1e-1`, `score≥20…60` all give 0.6472 / 30,673 / 137,490 | `validation/performance_metrics.tsv` | The input file is already FDR < 0.01-filtered, so looser cutoffs cannot add predictions | supported |
| B6 | Reference before propagation; shared-domain count | 30,190 curated pairs; 2,747 shared domains | `VALIDATION_PLAN.md` §1 L121–123 | Pre-propagation InterPro2GO size and shared-domain count for the current-release run | provisional — prose only; the TSV records only the post-propagation figure (47,393) |
| B7 | Domain-centric metrics, t0 associations, base configuration | 164,549 assoc.; 2,693 shared domains; 134,610 predicted pairs; 29,382 recovered; precision-lb 0.2183; recall 0.6306; F1 0.3243 | `validation/domain_centric_metrics.tsv`, row `base` | Association set scored directly against propagated InterPro2GO on the shared single-domain space | supported |
| B8 | Domain-centric metrics, t0 associations, + relative (parental-background) inference | 86,772 assoc.; 2,386 shared domains; 69,206 predicted pairs; 17,525 recovered; precision-lb 0.2532; recall 0.4301; F1 0.3188 | `validation/domain_centric_metrics.tsv`, row `relative` | As B7, after post-hoc relative-inference filtering (`validation/apply_relative_inference.py`) | supported |
| B9 | Relative inference raises domain-level precision-lb from 0.218 to 0.253 while halving the set | +0.035 precision-lb; recall 0.63 → 0.43 | derived from B7, B8 | Difference between two rows of the same file. **No confidence interval and no paired test** were computed for this difference | supported (derived); the *difference* is uncertainty-free and therefore uninterpretable as a significance claim |
| B10 | The reference is the **current** InterPro2GO, not a 2021 snapshot | — | `VALIDATION_PLAN.md` L299–302; `RESULTS.md` L74–76 | Means B7–B9 are not a temporal test on the domain side | supported |

---

## C. Protein-centric temporal (CAFA-style) benchmark (§2)

Source file for all of C: `validation/temporal_benchmark_metrics.tsv`.
Semantics from `validation/temporal_benchmark.py`
(`build_nk_benchmark_by_aspect`, `transfer_predictions_pscore`,
`naive_predictions`, `shuffle_domain_go`, `f_max`, `s_min`, `auprc`,
`filter_by_ic`).

### C.0 Design

| ID | Claim | Value | Source | What was actually measured | Verdict |
| --- | --- | --- | --- | --- | --- |
| C0a | Split | t0 = GOA release 205 (2021-04); t1 = current GOA (2026-06) | `VALIDATION_PLAN.md` §2 L155–156; `RESULTS.md` L21; `validation/temporal_benchmark.py` CLI defaults | ~5-year annotation gap | supported (t0 release id pinned; **t1 is a mutable "current" file and is not pinned**) |
| C0b | Benchmark definition | CAFA *no-knowledge*, per aspect | `validation/temporal_benchmark.py:80–125` | A protein is a target for an aspect iff it had **no annotation known to training** (non-IEA/`manual` evidence) in that aspect at t0 and gained experimental annotation by t1; truth = its **full** propagated t1 experimental term set, aspect roots excluded | supported |
| C0c | Predictor | dcGO p-score (sum of association scores over the protein's domains and their propagated ancestors, min–max normalised per protein) | `validation/temporal_benchmark.py:189–221` | Fang & Gough per-target transfer, restored | supported |
| C0d | Baseline 1 | CAFA `naive` — every protein predicted every term at its propagated t0 frequency | `validation/temporal_benchmark.py:223–235` | — | supported |
| C0e | Baseline 2 | `random_domain` — a **single seeded permutation** (seed 0) reassigning each domain another domain's whole term set | `validation/temporal_benchmark.py:237–260` | A null for the *transfer* step only; Fisher is not re-run | supported — and note the review's demand for many permutations (`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` L49–50) |
| C0f | IC definition | marginal `IC(t) = −log2 P(t)` from t0 frequencies | `validation/temporal_benchmark.py:131–157` | Not information accretion | supported |
| C0g | AUPRC construction | 51 score quantiles + a predict-nothing sentinel; trapezoidal integration over an upper envelope of precision-vs-recall | `validation/temporal_benchmark.py:293–318, 396–418` | Non-standard relative to reference CAFA tooling; unverified against an independent implementation (`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` L112–115) | supported (as a description of what the code does) |
| C0h | Precision/recall averaging | Precision averaged over proteins with ≥1 prediction at τ; recall over all benchmark proteins | `validation/temporal_benchmark.py:262–292` | CAFA convention; **prediction coverage (the count `m`) is computed but not reported in the committed metrics file** | supported |

### C.1 Benchmark cohort sizes (col `n_eval_proteins`)

| ID | IC floor | BP | MF | CC | Verdict |
| --- | --- | --- | --- | --- | --- |
| C1a | ≥ 0 | 324 | 418 | 572 | supported |
| C1b | ≥ 2 | 324 | **170** | 405 | supported |
| C1c | ≥ 4 | 318 | 162 | 252 | supported |
| C1d | ≥ 6 | 289 | 145 | 154 | supported |

C1e — **The IC filter changes the evaluation cohort.** MF drops from 418 to 170
proteins (−59%) between IC ≥ 0 and IC ≥ 2; CC from 572 to 405 (−29%); BP is
unchanged at IC ≥ 2. Derived from C1a–C1d. **supported (derived).** This is the
review item at `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` L63–65: the IC-filtered
rows do not measure the same proteins as the unfiltered rows, so the comparison
across IC floors is unpaired.

### C.2 F_max (col `f_max`)

All rows **supported**; values as stored (manuscript rounds to 3 d.p.).

| ID | Aspect | IC | dcGO | naive | random | dcGO ÷ random (derived) | Winner |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C2a | BP | ≥0 | 0.24841 | 0.11541 | 0.15771 | 1.58 | dcGO |
| C2b | BP | ≥2 | 0.17013 | 0.07062 | 0.05318 | 3.20 | dcGO |
| C2c | BP | ≥4 | 0.11535 | 0.03085 | 0.01893 | 6.09 | dcGO |
| C2d | BP | ≥6 | 0.07745 | 0.00978 | 0.00320 | 24.2 | dcGO |
| C2e | MF | ≥0 | 0.35966 | **0.46394** | 0.26200 | 1.37 | **naive** |
| C2f | MF | ≥2 | 0.36541 | 0.05303 | 0.08780 | 4.16 | dcGO |
| C2g | MF | ≥4 | 0.33672 | 0.04523 | 0.07179 | 4.69 | dcGO |
| C2h | MF | ≥6 | 0.21677 | 0.01821 | 0.00868 | 25.0 | dcGO |
| C2i | CC | ≥0 | 0.37963 | 0.34264 | 0.29146 | 1.30 | dcGO |
| C2j | CC | ≥2 | 0.23854 | 0.15282 | 0.07217 | 3.31 | dcGO |
| C2k | CC | ≥4 | 0.13447 | 0.09913 | 0.03148 | 4.27 | dcGO |
| C2l | CC | ≥6 | 0.12445 | 0.04353 | 0.01534 | 8.11 | dcGO |

C2m — "dcGO beats naive on F_max at IC ≥ 0 in BP and CC but **loses in MF**"
— derived from C2a, C2e, C2i. **supported.**
C2n — "naive's F_max collapses as the IC floor rises: BP 0.115→0.010,
MF 0.464→0.018, CC 0.343→0.044" — derived from the naive column.
**supported.**
C2o — "dcGO ÷ random ranges 1.3× to 25× across aspects and IC floors" — derived.
**supported.** No confidence interval exists for any of these ratios.

### C.3 AUPRC (col `auprc`) — the counter-result

All rows **supported**.

| ID | Aspect | IC | dcGO | naive | random | Winner |
| --- | --- | --- | --- | --- | --- | --- |
| C3a | BP | ≥0 | 0.13728 | **0.31353** | 0.07073 | **naive** |
| C3b | MF | ≥0 | 0.19483 | **0.32514** | 0.09889 | **naive** |
| C3c | CC | ≥0 | 0.23965 | **0.51302** | 0.17031 | **naive** |
| C3d | BP | ≥2 / ≥4 / ≥6 | 0.06893 / 0.03161 / 0.01679 | 0.04891 / 0.01161 / 0.00275 | 0.00789 / 0.00087 / 0.00003 | dcGO |
| C3e | MF | ≥2 / ≥4 / ≥6 | 0.22703 / 0.20986 / 0.12057 | 0.02464 / 0.01142 / 0.00559 | 0.01516 / 0.00708 / 0.00018 | dcGO |
| C3f | CC | ≥2 | 0.07298 | **0.11693** | 0.01008 | **naive** |
| C3g | CC | ≥4 | 0.02891 | **0.03754** | 0.00215 | **naive** |
| C3h | CC | ≥6 | 0.02497 | 0.01187 | 0.00039 | dcGO |

C3i — "naive attains higher AUPRC than dcGO in **all three aspects** at IC ≥ 0,
and in CC also at IC ≥ 2 and IC ≥ 4" — derived from C3a–C3g. **supported.**
This result is not reported in `README.md`, `RESULTS.md` or `VALIDATION_PLAN.md`,
which foreground F_max only. It materially qualifies the headline and is
therefore reported in the manuscript Results.

### C.4 S_min (col `s_min`; lower is better)

| ID | Aspect | IC | dcGO | naive | random | Winner | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C4a | BP | ≥0 | 99.70 | 104.86 | 103.85 | dcGO | supported |
| C4b | MF | ≥0 | 12.64 | 14.25 | 14.51 | dcGO | supported |
| C4c | CC | ≥0 | 20.04 | 21.34 | 21.29 | dcGO | supported |

C4d — "dcGO attains the lowest S_min of the three methods in all three aspects
at IC ≥ 0, including MF where it loses on F_max" — derived from C4a–C4c.
**supported.** Margins are small (≤ 5%) and have no confidence intervals.

### C.5 Method-piece experiment (configurations A–D)

| ID | Claim | Value | Source | Verdict |
| --- | --- | --- | --- | --- |
| C5a | A (max transfer, overall test): BP/MF/CC F_max 0.218 / 0.234 / 0.349; AUPRC 0.083 / 0.106 / 0.201 | as stated | `VALIDATION_PLAN.md` L259 | provisional — prose table only; `logs/bench_A.log` exists but `validation/bench_A/` is **untracked** (`git status`; `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` L83–85) |
| C5b | B (+ p-score transfer): 0.248 / 0.360 / 0.380; AUPRC 0.137 / 0.195 / 0.240 | as stated | `VALIDATION_PLAN.md` L260; F_max/AUPRC agree with the committed `temporal_benchmark_metrics.tsv` dcGO rows at IC ≥ 0 | supported for the B row (it is the primary configuration); the A/C/D rows remain provisional |
| C5c | C (+ relative inference): 0.219 / 0.273 / 0.356 | as stated | `VALIDATION_PLAN.md` L261 | provisional (as C5a) |
| C5d | D (both): 0.224 / 0.385 / 0.356 | as stated | `VALIDATION_PLAN.md` L262 | provisional (as C5a) |
| C5e | "The p-score transfer is the main protein-centric lever; relative inference helps MF" | qualitative | `VALIDATION_PLAN.md` L264–270 | provisional — rests on C5a/C5c/C5d, and no paired test or CI was computed for any A-vs-B difference |
| C5f | Earlier reported figures for the same experiment (BP 0.276→0.365, MF 0.319→0.446, CC 0.395→0.458) | as stated | memory `original-dcgo-method-audit.md` | **unsupported for the manuscript** — superseded by C5a–C5d and by the leak-free gate; excluded |

### C.6 Ancillary §2 claims

| ID | Claim | Value | Source | Verdict |
| --- | --- | --- | --- | --- |
| C6a | GO:0005515 *protein binding* accounts for 84.6% of human experimental MF annotations | 84.6% | `VALIDATION_PLAN.md` L179–180; `RESULTS.md` L17; memory `temporal-benchmark-s2-result.md` | provisional — prose only in three places; no committed artifact computes it |
| C6b | IC ≥ 2 ⇒ term in ≤ 25% of proteins; IC ≥ 6 ⇒ ≤ 1.6% | as stated | `RESULTS.md` L51; `VALIDATION_PLAN.md` L208 | supported (arithmetic identity of `−log2 P`: 2^−2 = 25%, 2^−6 = 1.56%) |
| C6c | `hyper_score` saturates at exactly 100 for ~37% of associations, which is why ranking uses −log10(p) | 37% | `VALIDATION_PLAN.md` L173–174, L250 | provisional — prose only |
| C6d | Earlier, leaky no-knowledge gate gave cohorts of 1537 / 1124 / 2305 | as stated | memory `temporal-benchmark-s2-result.md` | provisional — memory file only; used in the manuscript only as a narrative point about a fixed bug, without quoting the numbers |

---

## D. Surprise score — internal (same-data) results

| ID | Claim | Value | Source | What was actually measured | Verdict |
| --- | --- | --- | --- | --- | --- |
| D1 | Score definition | `−log10(q_emergence) × distinctness × novelty` | `SURPRISE_SCORE.md` L28–33; `src/surprise_score.py` | One-sided binomial tail of the observed co-occurrence count against a parts-only expectation (max of noisy-OR over constituents, best proper sub-combination rate, and the term's background rate; rates shrunk by one pseudo-observation), BH-corrected; × (1 − median region overlap); × a fixed novelty weight vs InterPro2GO | supported |
| D2 | Novelty weights 0.1 / 0.3 / 0.6 / 1.0 / 1.0 | as stated | `SURPRISE_SCORE.md` L71–78 | A prioritisation convention, explicitly **not** an inference | supported (as a description) |
| D3 | Ranked supra-domain candidates retained, current GO run | 22,243 | `results/domain_go_surprising.tsv` = 22,243 data rows (verified by line count) | Candidates surviving the `--max-overlap 0.5` redundancy filter | supported |
| D4 | Emergent at FDR ≤ 0.05, current GO run | 24 | `results/domain_go_surprising.tsv`, count of rows with `q_emergence ≤ 0.05` (verified) | — | supported |
| D5 | Supra-domain associations scored / dropped as redundant signatures, GO | 123,203 scored; 100,960 dropped (82%) | `SURPRISE_SCORE.md` L108–110 | — | provisional — prose only; consistent with D3 by arithmetic (123,203 − 100,960 = 22,243) but the pre-filter count is not in any committed file |
| D6 | Novelty breakdown of retained GO candidates | 10,353 novel / 6,014 no-reference / 3,840 refines / 1,716 curated / 320 implied | `SURPRISE_SCORE.md` L112–113 | Sums to 22,243 (= D3) | provisional — prose only, though the sum checks out |
| D7 | EC: ranked candidates retained | 1,401 | `results/domain_ec_surprising.tsv` = 1,401 data rows (verified) | — | supported |
| D8 | EC: emergent at FDR ≤ 0.05 | 1 | `results/domain_ec_surprising.tsv`, `q_emergence ≤ 0.05` (verified) | — | supported |
| D9 | EC: scored / dropped | 8,637 scored; 7,236 dropped (84%) | `SURPRISE_SCORE.md` L138–140 | Sums to 1,401 (= D7) | provisional — prose only |
| D10 | EC significant associations, current run | 11,745 | `results/domain_ec_associations_significant.tsv` (11,746 lines incl. header) | — | supported |
| D11 | Top-ranked GO architectures (SH2 + kinase → GO:0004715, 22/22, 3× lift; Ig C1-set + Ig-like → GO:0034987, 9/9, 7×; PH + EF-hand → GO:0004435, 7/7, 11×; Tyr-kinase + SAM → GO:0048013, 12/14, 4×; BTB/POZ + C2H2 → GO:0001227, 17/22, 3×; SH3 + P-loop → GO:0097120, 4/7, 26×) | as stated | `SURPRISE_SCORE.md` L117–124; underlying rows in `results/domain_go_surprising.tsv` | Support counts and lift for the six highest-surprise rows | provisional — the table is prose; individual rows exist in the TSV but were not re-verified field-by-field for this ledger |
| D12 | The tankyrase case (SAM + PARP catalytic → EC 2.4.2.30) has a parts-only expectation of 0.82 | 0.82 | `SURPRISE_SCORE.md` L150–153 | Why an architecturally real combination is not an *emergent prediction* | provisional — prose only |
| D13 | Earlier counts of "emergent" associations: 96 in EC, ~6,000 in GO | as stated | memory `emergent-domain-combination-predictions.md` | Superseded by the parts-baseline test (D4, D8) | **unsupported for the manuscript** — explicitly retracted in the same memory file; excluded except as a narrative "superseded" note without the figures |
| D14 | The score re-ranks a set already selected by the dcGO FDR filter, on the same proteins | — | `SURPRISE_SCORE.md` L252–261; `CLAUDE.md` Known Limitations | Internal consistency, not predictive power | supported |

---

## E. Surprise score — held-out temporal test

Source for E2–E12: `validation/temporal_surprise_metrics.tsv`; semantics from
`validation/temporal_surprise.py` (`score_association`, `acquisition_base_rates`,
`pool`, `_bootstrap_enrichment_ci`, `_budget_enrichment`, `compare_rankings`).

| ID | Claim | Value | Source | What was actually measured | Verdict |
| --- | --- | --- | --- | --- | --- |
| E1 | Candidate pool | 22,376 t0 supra-domain candidates, of which 10,136 make ≥1 standing prediction | `results_t0_2021/domain_go_surprising.tsv` = 22,376 rows (verified); `validation/temporal_surprise_associations.tsv` = 10,136 rows (verified); `validation/temporal_surprise_metrics.tsv` rows `surprise all (10136)` | An association contributes only if ≥1 carrier lacked the term at t0 | supported |
| E2 | Whole-pool result | 10,136 assoc.; 170,416 predictions; 2,181 hits; hit rate 1.28%; expected 0.10%; **enrichment 12.49×**, bootstrap CI **[10.61, 14.13]** | `validation/temporal_surprise_metrics.tsv`, row `surprise all (10136)` | Predictions = (carrier, term) pairs lacking the term at t0 under non-IEA propagated evidence; hits = annotated at t1 under experimental propagated evidence; expected = prediction-weighted mean of each term's own acquisition rate among domain-carrying proteins lacking it at t0 | supported |
| E3 | Published CI for E2 | [10.9, 14.4] | `SURPRISE_SCORE.md` L185; `VALIDATION_PLAN.md` L531 | — | provisional — **does not match either committed row** (`surprise all` = [10.61, 14.13]; `dcgo all` = [10.94, 13.75]). Manuscript quotes the committed file and flags the discrepancy |
| E4 | Expected hits, whole pool | ~175 | `SURPRISE_SCORE.md` L192 | 170,416 × 0.10% ≈ 175 | supported (derived, consistent with E2) |
| E5 | Bootstrap unit | associations, 500 or 1000 resamples, percentile CI | `validation/temporal_surprise.py:177–211`; `temporal_breadth.py` default `--bootstrap 500` | Resampling associations, not proteins — the conservative choice given within-association protein correlation | supported |
| E6 | Rank-slice results | surprise top-100: 16.46× [6.61, 24.34]; top-500: 14.78× [11.00, 21.82]; dcGO-q top-100: 10.31× [5.73, 21.33]; top-500: 13.23× [8.93, 19.08] | `validation/temporal_surprise_metrics.tsv` | Top-K by each ranking | supported. Note `SURPRISE_SCORE.md` L186–189 prints CIs of [6.2, 26.4], [10.6, 21.6], [5.9, 21.4], [9.5, 19.2] — **all differ from the committed file** |
| E7 | Top-K is not a fair comparison | surprise ≈ 4.8 predictions per association (481/100), dcGO-q ≈ 73 (7,293/100) | derived from `temporal_surprise_metrics.tsv` rows `surprise top-100`, `dcgo top-100` | Motivates the prediction-budget matching | supported (derived) |
| E8 | Budget-matched enrichments | @2,000: 15.48 vs 5.28; @10,000: 21.17 vs 13.22; @40,000: 11.63 vs 10.81 | `validation/temporal_surprise_metrics.tsv`, rows `surprise @Npreds` / `dcgo @Npreds` | Enrichment of the top associations under each ranking until the prediction budget is exhausted | supported |
| E9 | **Paired bootstrap of the difference (surprise − dcGO-q)** | @2,000: +10.19, 95% CI **[−87.28, +18.95]**, 70.0% of resamples favour surprise; @10,000: +7.96, CI **[−8.93, +7.06]**, 32.0%; @40,000: +0.82, CI **[−1.01, +4.54]**, 80.0% | `validation/temporal_surprise_metrics.tsv`, "# paired ranking comparison" block | Re-ranks both ways inside each resample of the shared candidate pool | supported |
| E10 | **Zero lies inside every paired interval at every budget** | — | derived from E9 | The negative result: no demonstrated ranking advantage | supported (derived) |
| E11 | Published paired figures | @2,000 [−85.26, +20.87] 75%; @10,000 [−9.39, +10.23] 54%; @40,000 [−2.15, +5.25] 82% | `SURPRISE_SCORE.md` L206–208; memory `paired-bootstrap-for-ranking-comparisons.md` | — | provisional — **differs from the committed file at every budget** (E9). Manuscript quotes E9 and flags it |
| E12 | Internal inconsistency in the committed paired result at budget 10,000 | point estimate +7.96 lies **outside** its own CI [−8.93, +7.06], and only 32.0% of resamples favour surprise despite a positive point estimate | derived from E9 | Percentile bootstrap intervals need not contain the point estimate, but combined with the sign disagreement this indicates a heavy-tailed, unstable statistic | supported (derived) — raised as an open question rather than interpreted |
| E13 | Selection bias toward rare terms | @10,000 budget, surprise takes 580 associations at 0.05% base rate and 0.98% hit rate; dcGO-q takes 169 at 0.14% base rate and 1.83% hit rate | `validation/temporal_surprise_metrics.tsv`, rows `surprise @10000 preds`, `dcgo @10000 preds` | Lower raw hit rate, higher enrichment ⇒ rarer, more specific terms | supported. `SURPRISE_SCORE.md` L221–223 states "1.0% vs 1.8%" and base rates "0.0% vs 0.1%", consistent to rounding |
| E14 | The sharp end of the ranking is untestable | surprise top-25: 117 predictions, 0 hits, expected rate 0.0008 ⇒ ≈0.09 expected hits | `validation/temporal_surprise_metrics.tsv`, row `surprise top-25` | An uninformative cell, not a negative result | supported. `SURPRISE_SCORE.md` L232 says "0.14 expected hits"; the committed file's 4-d.p. expected rate implies ≈0.09. Manuscript uses "fewer than 0.2 expected hits" and flags the discrepancy |
| E15 | Structural emergence/testability tension | qualitative | `SURPRISE_SCORE.md` L229–238; `VALIDATION_PLAN.md` L542–547 | Emergence requires that carriers are already nearly all annotated, which leaves few standing predictions | supported (as a stated design tension) |

---

## F. Multi-ontology breadth (§2 breadth subsection) — ALL PROVISIONAL

**Blocking defect (identified 2026-08-04, after the numbers below were produced).**
The UniProt-native annotation sources parse the entire Swiss-Prot flat file with
no taxonomic restriction, and the Fisher protein universe is built as the
**union** of the annotation map's proteins and the domain map's proteins:

- `src/uniprot_annotation_source.py` — no `OX` / `NCBI_TaxID` handling anywhere in
  the module (only a Reactome `species_prefix` filter for hierarchy *edges*, at
  L590–615); the `parse()` implementations return every entry in the file.
- `src/sparse_fisher.py:169` — `all_proteins = sorted(set(protein_domains.keys()) | set(protein_go.keys()))`.
- `run_dcgo_human.py:400–401` intersects the two maps only to choose which
  proteins get *features*, not to set the matrix row space.

Consequence: for every UniProt-native ontology the contingency tables are built
over a universe on the order of the whole of Swiss-Prot (575,503 entries, A11)
rather than the ~20k human proteins, inflating the `d` cell and therefore the
significance of every association. GO is sourced from the species-specific GOA
GAF and Expasy `enzyme.dat`/other non-UniProt sources are unaffected by this
particular defect, but every ontology trained from `uniprot_sprot.dat.gz` is.

**Therefore every row in F is marked provisional and pending recomputation.**
Ledger IDs F1–F9 record the value as currently committed, so that the corrected
run can be diffed against it.

Source file: `validation/temporal_breadth_metrics.tsv` (and
`validation/temporal_breadth_go.tsv` for the GO anchor). Statistic identical to
E2. t0 = archived Swiss-Prot 2021_02 (07-Apr-2021) for the UniProt layers,
GOA release 205 for GO.

| ID | Ontology (universe) | assoc. | predictions | hits | hit rate | expected | enrichment [95% CI] | Verdict |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| F1 | go (all-domains) — anchor | 91,830 | 3,106,234 | 106,224 | 3.420% | 0.298% | **11.49 [11.12, 11.99]** | supported (corrected 2026-08-04; was 11.33) — GO was never affected by the species defect, and its stability across the correction is the control. Shares the remaining §2 caveats |
| F2 | reactome (in-scope) | 34,397 | 760,073 | 1,430 | 0.188% | 0.016% | 11.44 [8.85, 14.60] | supported (corrected; was 7.99). Interval now overlaps F1, so "weaker than GO" is no longer supported — but on 1,430 hits vs GO's 106,224, read as indistinguishable at this power |
| F3 | reactome (all-domains) | 38,187 | 1,194,388 | 1,430 | 0.120% | 0.012% | 9.66 [7.81, 11.96] | supported (corrected; was 7.33) |
| F4 | keyword (all-domains) | 35,981 | 889,117 | 13,140 | 1.478% | 0.440% | 3.36 [3.16, 3.57] | supported (corrected; was 1.70 — the largest single change, and the intervals do not overlap) |
| F5 | subcellular (in-scope) | 6,340 | 172,679 | 3,798 | 2.199% | 0.592% | 3.72 [3.50, 3.93] | supported (corrected; was 2.85) |
| F6 | cofactor (in-scope) | 291 | 1,404 | 235 | 16.738% | 4.722% | 3.54 [2.10, 4.45] | supported (corrected; was 3.19) |
| F7 | cofactor (all-domains) | 747 | 12,838 | 235 | 1.831% | 0.416% | 4.40 [1.81, 6.89] | supported (corrected; was 3.85) |
| F8 | complex (in-scope) | 958 | 9,904 | 2 | 0.020% | ~0% | 61.14 [0.00, 130.93] | **no demonstrated signal** (corrected; was 265.45 [164.55, 455.34]). The interval now includes zero on 2 hits. ~1.5 proteins per ComplexPortal entry makes the denominator round to zero, so the magnitude was never interpretable — but the corrected run shows the *existence* claim was not supported either |
| F9 | disease (in-scope / all-domains) | 27 / 44 | 90 / 369 | 0 / 0 | 0% | 0% | `nan` (undefined) | supported — unchanged by the correction (OMIM phenotype xrefs are essentially human-only, so this layer had no contamination). 0 hits *and* ~0 expected: uninformative, not a negative result |
| F10 | ligand | — | — | — | — | — | not testable at this split | supported as a *design* statement: UniProt's structured `FT /ligand_id="ChEBI:…"` postdates the 2021_02 release (`VALIDATION_PLAN.md` L402–405; `TODO.md` L73–76) |
| F11 | Published rounding of F1–F6 as "GO 11.5×, reactome 11.4×, subcellular 3.7×, cofactor 3.5×, keyword 3.4×" | as stated | `VALIDATION_PLAN.md` §2 breadth; `TODO.md`; `docs/uniprot_ontology_survey.md` | supported — consistent with the corrected F1–F7 to rounding. The pre-correction rounding (GO 11.3, reactome 8.0, cofactor 3.2, subcellular 2.9, keyword 1.7) is superseded |
| F12 | The breadth rows pool single domains **and** supra-domains | — | `validation/temporal_breadth.py:288`; `VALIDATION_PLAN.md` L426–427 | supported — so the breadth test does **not** isolate the emergent (supra-domain) claim |
| F13 | No evidence filter exists off GO; an automated annotation added between snapshots counts as a hit | — | `validation/temporal_breadth.py` docstring L34–38; `VALIDATION_PLAN.md` L416–420 | supported |
| F14 | The anchor bounds how much the looser truth matters: GO scores 11.49× under the loose protocol vs 12.49× under the strict one | derived from F1 and E2 | `VALIDATION_PLAN.md` §2 breadth | provisional — the two figures also differ in candidate set (E2 is supra-domains only, F1 pools all associations), so this is not a clean protocol contrast |
| F15 | Per-ontology coverage on the human domain set (e.g. keyword 18,859 proteins / 720 terms; subcellular 16,750 / 261; disease 5,029 / 6,904; ligand 4,627 / 448; cofactor 1,801 / 46) | as stated | `docs/uniprot_ontology_survey.md` L144–160 | supported — this is a coverage survey, independent of the Fisher universe defect |
| F16 | DR-database survey ratios (keyword 28.4 proteins/term; Reactome 5.0; ComplexPortal 1.5; MEROPS 1.0; DrugBank 0.4) | as stated | `docs/uniprot_ontology_survey.md` L37–50; `docs/dr_survey.tsv` | supported |
| F17 | The `disease` t0 run produced only 53 significant associations from 6,904 OMIM phenotypes over 5,029 proteins | 53 | `TODO.md` L90–92; `results_t0_2021/domain_disease_associations_significant.tsv` = 54 lines incl. header (verified) | supported |

---

## G. Comparison to the original dcGO

| ID | Claim | Value | Source | Verdict |
| --- | --- | --- | --- | --- |
| G1 | The original dcGO covered all completely sequenced genomes: 2,414 genomes + UniProt, > 80M sequences | as stated | memory `original-dcgo-methodology.md`; `VALIDATION_PLAN.md` L437–438, attributed to `docs/gks1080.pdf` and `docs/1471-2105-14-S3-S9.pdf` | provisional — **the PDFs could not be read in this environment** (no `pdftotext`/`pypdf`/`poppler`); the figures come from a prior reading recorded in prose |
| G2 | Original domain universe: SCOP superfamily/family via SUPERFAMILY HMMs, plus Pfam | as stated | as G1 | provisional (same reason) |
| G3 | Original significance threshold: FDR < 1e-3 (vs 0.01 here) | as stated | as G1 | provisional (same reason) |
| G4 | Original combined two hypergeometric tests — *overall* (all-UniProt background) and *relative* (background = proteins annotated to all direct parent terms) — keeping the larger p / smaller h-score | as stated | as G1 | provisional (same reason) |
| G5 | Original transfer: sum of h-scores over supporting domains, min–max normalised per target | as stated | as G1 | provisional (same reason) |
| G6 | Original validation was protein-centric CAFA precision–recall, BP + MF only | as stated | memory `original-dcgo-method-audit.md`; `VALIDATION_PLAN.md` L243 | provisional (same reason) |
| G7 | InterPro integrates SUPERFAMILY and Pfam, so `protein2ipr` col. 4 carries `SSF*`/`PF*` signatures — the domain universes are not disjoint | as stated | `VALIDATION_PLAN.md` L441–445; memory `original-dcgo-methodology.md` | provisional (rests on a reading of `protein2ipr` recorded in prose; the claim about InterPro's member databases is independently well established) |
| G8 | Human subset contains 8,256 distinct Pfam and 911 distinct SUPERFAMILY signatures rolled into 19,534 InterPro entries | as stated | memory `original-dcgo-methodology.md` | provisional — memory only; not recomputed |
| G9 | The 2023 dcGO release reports ~1,000 Pfam and ~800 InterPro domains | as stated | `VALIDATION_PLAN.md` L465–468 | **unsupported** — the source document itself instructs "verify that figure against the actual download before citing it". Excluded from the manuscript |
| G10 | The 2023 dcGO release added MONDO, EFO, KEGG, Reactome, PANTHER, WikiPathways, MitoCarta, DGIdb, Open Targets, ENRICHR/TRRUST, MSigDB, and dropped EC, UniPathway, UniProt keywords, DrugBank ATC | as stated | `VALIDATION_PLAN.md` L469–473, attributed to `docs/EMS185259.pdf` | provisional — PDF unreadable here |
| G11 | No quantitative comparison against original dcGO output has been run | — | `VALIDATION_PLAN.md` §3 checkboxes L477–485 all unchecked; `TODO.md` L62 | supported |

---

## H. Limitations carried from the external review

All rows sourced from `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` (review date
2026-07-15) and **supported** as accurate statements of that document.

| ID | Limitation | Locator |
| --- | --- | --- |
| H1 | Model selection and final evaluation used the same 2021→2026 split; no untouched evaluation interval | L42–46 |
| H2 | No protein-level bootstrap CIs for F_max, AUPRC, or paired differences vs baselines | L47–48 |
| H3 | The random-domain null is a **single** seeded shuffle, not a permutation distribution with an empirical p-value | L49–50 (corroborated by `validation/temporal_benchmark.py:237–260`) |
| H4 | Component ablation (single domains / +supra / +shrinkage / +True Path) not run | L51–52; `VALIDATION_PLAN.md` §4 checkboxes unchecked |
| H5 | The "empirical-Bayes shrinkage" geometrically interpolates observed and constituent p-values with a hand-set decay; it is not a fitted model, the outputs are not shown to be valid p-values, and BH over them does not guarantee nominal FDR | L53–58 |
| H6 | The primary endpoint was not pre-specified; the unfiltered result loses to naive for MF and the headline depends on IC filtering | L59–62 |
| H7 | IC thresholds change the evaluation cohort (see C1e) | L63–65 |
| H8 | The t1 GOA snapshot and the GO release used for propagation are mutable, unpinned files | L66–67, L37–38 |
| H9 | Temporal look-ahead: domain architectures come from the **current** `protein2ipr`, so this is an annotation-temporal, not a prospective, simulation | L68–70; `VALIDATION_PLAN.md` L421–423, L548–549 |
| H10 | Comparators are weak: no original dcGO output, no independent protein-function or domain-based predictor, no BLAST transfer baseline | L71–74; `VALIDATION_PLAN.md` L310–311 |
| H11 | No external validation axis: one species, one time interval | L75–76 |
| H12 | Benchmark artifacts `validation/bench_A`–`bench_D` were untracked at review time and remain untracked | L83–85; `git status` |
| H13 | `calculate_hypergeometric_score` falls back to a value of 50.0 on numerical error — a plausible medium-confidence score produced by a failure | L91–93 |
| H14 | `--num-cores` is logged but unused by the Cython Fisher implementation | L94–95 |
| H15 | BH is applied across a highly dependent hierarchical hypothesis family (GO terms × domains × supra-domains) without simulation or hierarchical multiple-testing correction | L124–126 |
| H16 | InterPro2GO recovery is coverage against an incomplete positive reference — not precision, not independent validation | L121–123, L140 |
| H17 | Prediction coverage is not reported alongside F_max | L110–111 |
| H18 | F_max/AUPRC not verified against an independent CAFA evaluation implementation; 51-quantile sweep and upper-envelope trapezoidal AUPRC are undocumented choices | L112–115 |
| H19 | No minimum-support / effect-size policy; FDR significance alone retains associations from sparse tables; contingency cells and odds-ratio CIs are not reported | L116–120; `VALIDATION_PLAN.md` L557–563 |
| H20 | Reviewer's overall verdict: the evidence supports predictive signal in this retrospective benchmark, **not** general performance, calibration, or superiority | L10–13, L143–155 |
| H21 | `odds_ratio` prints `0.0000` when `d = 0` and `inf` when `b·c = 0`; no Haldane correction | `VALIDATION_PLAN.md` L560–563; visible in `README.md` L271–272 example output | 

---

## I. Claims considered and REJECTED (not in the manuscript)

| ID | Rejected claim | Why |
| --- | --- | --- |
| I1 | "dcGO outperforms the naive baseline" (unqualified) | Contradicted by C2e (MF F_max) and C3a–C3g (AUPRC in all three aspects at IC ≥ 0). |
| I2 | "dcGO beats both baselines in every aspect" | True only for F_max under an IC floor (C2f–C2l); false for AUPRC (C3f, C3g) and for MF F_max unfiltered. |
| I3 | "The method controls FDR" | H5: the shrinkage-transformed quantities have not been shown to be valid p-values, and H15: BH is applied across a dependent family. |
| I4 | "The surprise score ranks better than the dcGO q-value" | E9/E10: paired interval spans zero at every budget. |
| I5 | "The associations generalise across 19 ontologies" | Only 7 were trained at t0; 1 was untestable, 1 undefined, 1 degenerate (F1–F10); and all UniProt-native rows are affected by the species defect. |
| I6 | Any statement of *calibration* | No calibration analysis exists in the repository. |
| I7 | "165,823 / 164,549 associations are validated predictions" | They are hypotheses passing an FDR filter with no minimum-support or effect-size policy (H19). |
| I8 | "96 emergent EC associations" / "~6,000 emergent GO associations" | D13: retracted in the source memory and superseded by D4/D8. |
| I9 | The 2023 dcGO domain counts | G9: the repository's own note forbids citing them unverified. |
| I10 | Any speedup or parallel-scaling claim | H14: `--num-cores` does not drive the Fisher implementation; no timing benchmark with a controlled comparison exists. |

---

## J. References — verification status

| ID | Reference | How verified | Verdict |
| --- | --- | --- | --- |
| J1 | Fang H, Gough J. *A domain-centric solution to functional genomics via dcGO Predictor.* BMC Bioinformatics 2013;14(Suppl 3):S9. doi:10.1186/1471-2105-14-S3-S9 | `docs/1471-2105-14-S3-S9.pdf` XMP metadata (`dc:title`, `prism:doi`, `prism:publicationName`, volume 14, Suppl 3) + web search | verified |
| J2 | Fang H, Gough J. *dcGO: database of domain-centric ontologies on functions, phenotypes, diseases and more.* Nucleic Acids Res 2013;41(D1):D536–D544. | `docs/gks1080.pdf` (NAR pagination "536..544" in metadata; embedded links to `supfam.org/SUPERFAMILY/dcGO`) + web search | verified |
| J3 | *The dcGO Domain-Centric Ontology Database in 2023: New Website and Extended Annotations for Protein Structural Domains.* J Mol Biol 2023;435(14):168093. doi:10.1016/j.jmb.2023.168093 | `docs/EMS185259.pdf` XMP metadata (title, journal, volume 435, issue 14, cover date 2023-07-15, DOI) | title/venue verified; **author list not independently verified** (recorded as "Bao *et al.*" in `VALIDATION_PLAN.md` L459) |
| J4 | The Gene Ontology Consortium. *Expansion of the Gene Ontology knowledgebase and resources.* Nucleic Acids Res 2017;45(D1):D331–D338. doi:10.1093/nar/gkw1108 | `docs/gkw1108.pdf` decompressed content stream contains `academic.oup.com/nar/article/45/D1/D331`; DOI confirmed by web search | verified |
| J5 | Benjamini Y, Hochberg Y. *Controlling the false discovery rate.* J R Stat Soc B 1995;57(1):289–300. | web search (OUP/Wiley records) | verified; not present in `docs/` |
| J6 | Zhou N, Jiang Y, Bergquist T, *et al.* *The CAFA challenge reports improved protein function prediction…* Genome Biol 2019;20:244. doi:10.1186/s13059-019-1835-8 | web search | verified; not present in `docs/` |
| J7 | Paysan-Lafosse T, Blum M, Chuguransky S, *et al.* *InterPro in 2022.* Nucleic Acids Res 2023;51(D1):D418–D427. doi:10.1093/nar/gkac993 | web search | verified; not present in `docs/`. Cited as the InterPro reference; **the exact InterPro release used by the runs is not recorded anywhere in the repository** |
| J8 | The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2025.* Nucleic Acids Res 2025;53(D1):D609–D617. doi:10.1093/nar/gkae1010 | web search | verified; not present in `docs/`. The runs used a 2026-07 release and an archived 2021_02 release, neither of which has a citable release identifier recorded |
| J9 | Any CAFA F_max / S_min formal definition beyond J6 | — | not cited; the manuscript describes the implementation in `validation/temporal_benchmark.py` instead |
