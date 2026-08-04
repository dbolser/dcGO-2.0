# Benchmark artifacts — what each committed metrics file is

`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` (P1) flagged that the `bench_A`–`bench_D`
outputs and their logs were untracked, "leaving their provenance unclear". This
file commits them and records everything that can be established about how they
were produced — including, honestly, the parts that cannot.

**Read the last section before citing any of these numbers.** Four of the five
runs are variant selection carried out on the evaluation split, which is a
known P0 publication blocker, not a set of independent results.

## The runs

All five were produced by `validation/temporal_benchmark.py` on 2026-07-09
against the §2 split — t0 = GOA release 205 (2021-04), t1 = GOA 2026-06 — over
the human `protein2ipr` subset, with `go-basic.obo` for propagation, supra-domains
enabled, `--score-column p_value` ranked as `-log10`, and the CAFA no-knowledge
truth design.

| Artifact | Log | Ran | Transfer | n_eval (BP) | BP F_max | Status |
| --- | --- | --- | --- | ---: | ---: | --- |
| `validation/temporal_benchmark_metrics.tsv` | `logs/bench_primary.log` | 23:13 | `pscore` | 324 | 0.2484 | **The reported §2 result.** |
| `validation/bench_A/…` | `logs/bench_A.log` | 23:18 | `max` | 324 | 0.2183 | Transfer-rule variant |
| `validation/bench_C/…` | `logs/bench_C.log` | 23:22 | `max` | 324 | 0.2187 | Transfer-rule variant |
| `validation/bench_D/…` | `logs/bench_D.log` | 23:25 | `pscore` | 324 | 0.2237 | Transfer-rule variant |
| `validation/bench_B/…` | `logs/bench_B.log` | 22:00 | — | 1537 | 0.3654 | **Superseded — leaky gate** |

A sixth artifact was added on 2026-08-04:

| Artifact | Log | Ran | What |
| --- | --- | --- | --- |
| `validation/temporal_benchmark_permutation_null.tsv` | `logs/bench_permutation_null.log` | 2026-08-04 | The random-domain null as a **distribution** over 100 seeded permutations (seeds 0–99), replacing the single shuffle. Produced by the same command as `bench_primary` plus `--n-permutations 100`; its `temporal_benchmark_metrics.tsv` reproduced the archived one exactly (see below). |

## What is established, and what is inferred

**Established from the logs.** The t0/t1 GAF files, the GO release (38,245 terms,
71,895 relationships), the supra-domain count (413,251), the aspects, and the
transfer rule (the runner logs `Transferring dcGO predictions … (max|pscore)`)
are all recorded directly.

**`bench_primary`'s input is now established by reproduction (2026-08-04).**
An earlier version of this file *inferred*, from file timestamps, that
`bench_primary` had read `results_t0_2021/domain_go_associations_relative.tsv`.
**That inference was wrong.** Re-running the harness against the plain
`results_t0_2021/domain_go_associations_significant.tsv` with the documented
defaults —

```bash
uv run python validation/temporal_benchmark.py \
    --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \
    --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
    --predictions results_t0_2021/domain_go_associations_significant.tsv \
    --min-ic 0 --min-ic 2 --min-ic 4 --min-ic 6
```

— reproduces `validation/temporal_benchmark_metrics.tsv` **exactly**: every
`f_max`, `f_max_tau`, `auprc`, `n_eval_proteins` and every `naive` and
`random_domain` row is bit-identical, with differences only in the last one or
two digits of `s_min` (floating-point summation order). The §2 headline is
therefore the plain significant table under `--transfer pscore`, exactly as
`RESULTS.md`'s reproduce block says, and the relative-inference table plays no
part in it.

**Still inferred, not recorded: `bench_A`, `bench_C`, `bench_D`.** The runner
did not log its command line or its `--predictions` path, so which association
table fed each of the three *variant* runs remains a reconstruction from
timestamps. `bench_A` and `bench_C` differ only in their score scale (F_max τ
7.92 vs 1.31) at near-identical F_max (0.2183 vs 0.2187), which is what two
different input tables under the same transfer rule would look like. Treat them
as archived exploratory output whose exact inputs are not certain, and do not
cite them as independent measurements. Closing that gap for future runs is the
P1 run-manifest item.

**`bench_B` is superseded.** Its 1,537-protein BP cohort is the *leaky*
no-knowledge gate: it gated on experimental-only evidence while the pipeline
trains on all non-IEA evidence, so proteins whose labels the model had already
seen entered the benchmark. Fixing the gate to use the same evidence space as
training cut the cohort to BP 324 / MF 418 / CC 572. `bench_B`'s higher F_max
(0.3654 BP) is that leak, not a better result. It is kept only so the
correction is auditable.

## The methodological caveat that matters most

`bench_primary`, `A`, `C` and `D` are four configurations of the same method
compared on the same evaluation split, and the configuration reported in
`RESULTS.md` is the one that scored best on it. That is model selection on the
test set. The P0 remedy — freeze the transfer rule and the IC threshold on a
development interval or a nested temporal split, then evaluate once on an
untouched interval — has **not** been done. Until it is, the §2 numbers should
be read as a retrospective descriptive result, not as an unbiased estimate of
performance.
