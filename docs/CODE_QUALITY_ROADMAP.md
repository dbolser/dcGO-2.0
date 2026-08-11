# Code quality roadmap

Review date: 2026-08-11

This is the maintained engineering companion to `TODO.md`. The scientific
roadmap belongs in `VALIDATION_PLAN.md`; this document tracks code structure,
testability, consistency, and maintenance debt.

## Current assessment

The numerical core is stronger than the repository surface suggests. Sparse
contingency construction, Fisher testing, BH correction, parsers, hierarchy
helpers, manifests, and surprise-score primitives have extensive regression
coverage. Most debt is concentrated in orchestration, manuscript-producing
validation programs, duplicated representations, and compatibility layers.

Baseline measured on 2026-08-11:

- 632 tests pass.
- Expanded coverage across `src`, the runner, scripts, and validation is 63%.
- `src` is generally well covered; the main runner and several analysis CLIs
  are not exercised in-process.
- The installed CLI works end to end, but most of its implementation lives in
  one large `main()` transaction.

## Work streams

### 1. Canonical propagation semantics

Status: in progress.

GO and non-GO paths currently use separate propagation implementations with
different ordering rules. A directly significant parent can be hidden by an
earlier propagation from a child, which makes `annotation_type`,
`direct_source_term`, and direct/propagated counts order-dependent.

Acceptance criteria:

- One explicit merge policy is used by every ontology.
- Direct evidence takes precedence over propagated evidence.
- Competing propagated sources resolve deterministically using evidence strength.
- Collision and input-order regression tests cover the contract.

### 2. Decompose the runner and adopt a thin Typer CLI

Status: in progress as a prototype.

Typer is appropriate for user-facing parsing, validation, discoverable help,
and future subcommands. It should not merely decorate the existing monolithic
`main()`. Extract stages with typed inputs and outputs first or alongside the
CLI migration:

1. resolve configuration and inputs;
2. load and restrict the protein universe;
3. infer association records;
4. apply reporting policy;
5. propagate annotations;
6. export results and finalize provenance.

Preserve existing command flags and exit behavior during migration. The
installed `dcgo` entry point must continue to pass a wheel-install smoke test.

### 3. Test scientific command-line programs

Status: in progress.

The core unit suite is strong, but expanded coverage exposes weakly tested
manuscript-producing programs:

| Module | Baseline coverage |
| --- | ---: |
| `run_dcgo_human.py` | 24% |
| `scripts/download_data.py` | 19% |
| `scripts/rank_surprising_associations.py` | 25% |
| `validation/compare_original_dcgo.py` | 31% |
| `validation/temporal_benchmark.py` | 52% |
| `validation/validate_results.py` | 33% |
| `validation/apply_relative_inference.py` | 0% |
| `validation/domain_centric_eval.py` | 0% |
| `validation/sprint1_validation.py` | 0% |

Prioritize fixture-based CLI tests that assert output schemas, counts, failure
behavior, and provenance—not coverage percentage alone.

### 4. Retire the `OntologyProcessor` grab bag

Status: queued after propagation consolidation.

`OntologyProcessor` currently owns OBO loading, graph direction, caches,
parental-background inference, propagation, validation, statistics, DataFrame
export, test fixtures, and an executable demonstration. Separate this into GO
graph access, parental inference, and shared propagation. Move synthetic OBO
fixtures and demonstrations into tests/examples.

### 5. Consolidate analysis utilities and schemas

Status: starting with validation tests.

Repeated implementations exist for association TSV loading, compressed text
opening, hashing, GO propagation, Fisher-table reconstruction, and result
serialization. Consolidate only after behavior is pinned by tests. Introduce a
canonical `AssociationRow` so significant, top-100, propagated, validation, and
downstream readers share named fields rather than parallel arrays and repeated
index arithmetic.

### 6. Replace broad exception swallowing

Status: queued.

Catch expected data errors explicitly. Aggregate rejected rows with counts and
examples, and fail when a configurable rejection threshold is exceeded.
Unexpected exceptions must escape rather than silently deleting scientific
results.

### 7. Make typing and dependency policy real

Status: queued.

The strict mypy configuration is not enforced and mypy is absent from the
normal development group. `joblib` and `statsmodels` remain production
dependencies although the current implementation imports neither. Align the
declared toolchain with CI, remove obsolete dependencies and comments, and
either deprecate or remove no-op parallelism parameters.

### 8. Remove misleading human/GO-specific internal names

Status: queued behind the CLI seam.

The runner supports multiple species and ontologies but still uses names such
as `run_dcgo_human.py`, `protein_go_map`, and `go_list`. Introduce generic
internal names (`protein_terms`, `term_ids`) while retaining compatibility
aliases and established filenames at public boundaries until a versioned
deprecation is possible.

### 9. Move historical narrative out of implementation modules

Status: queued.

Long bug postmortems and exact historic run counts are useful evidence but make
current contracts difficult to see. Preserve them in design records or issue
links; keep code comments focused on invariants, failure policy, and surprising
current constraints.

## Change strategy

Use small, independently reviewable PRs. Correctness contracts and tests land
before structural consolidation. Avoid combining numerical changes with CLI,
packaging, or naming migrations, and demonstrate byte- or schema-equivalence
where compatibility is promised.
