# Reproducibility checklist

Use this checklist for any result intended for a report, release, or paper. A
run is reproducible only when its code, exact inputs, parameters, environment,
and outputs can all be identified.

## Before the run

- [ ] Start from a clean Git checkout and record the commit or release tag.
- [ ] Install from the committed `uv.lock`: `uv sync --frozen`.
- [ ] Record the Python and uv versions: `python --version` and `uv --version`.
- [ ] Prefer dated upstream releases over `current_release` URLs.
- [ ] Preserve the raw downloaded files; do not replace them in place.
- [ ] Record licenses and citation requirements for every upstream dataset.
- [ ] Choose a new output directory so an older manifest cannot be overwritten.
- [ ] Record the complete command, including defaults that affect the analysis.

## Run and verify

```bash
uv run python run_dcgo_human.py \
  --species human \
  --ontology go \
  --enable-true-path \
  --output-dir results/reproduction_run
```

The run writes `run_manifest_<ontology>.json` (e.g. `run_manifest_go.json`) in
its output directory — one manifest per ontology, so parallel runs sharing a
directory do not clobber each other's provenance. Check that:

- [ ] `status` is `completed`; `running` means the run did not finalize.
- [ ] `git.commit` is the intended revision and `git.dirty` is `false`
      (`git.dirty` covers tracked modifications only; `git.untracked_files`
      counts everything else in the checkout).
- [ ] `software.dependency_lock.sha256` matches the `uv.lock` you installed from.
- [ ] Every input has a SHA-256 digest and byte size.
- [ ] Release headers are present under `release_metadata` when the format
      supplies one (GAF `!gaf-version`/`!date-generated`, OBO `data-version`,
      UniProt vocabulary `Release:`, Expasy `Release of`).
- [ ] Each mutable upstream URL is paired with the hash of the exact local file.
- [ ] `analysis.ontology` names the ontology actually resolved, whether True Path
      propagation ran, and which files supplied the term hierarchy.
- [ ] `analysis.thresholds` matches the intended method (evidence filter, FDR
      threshold and method, Fisher alternative, supra-domain length and minimum
      domain length, shrinkage settings, parental-background test settings).
- [ ] `parameters` and `command` match the intended invocation.
- [ ] Output hashes match the files submitted for analysis or publication.
- [ ] `summary` counts agree with the log and the manuscript tables.

Hashing is intentionally part of the run. For large inputs it adds sequential
I/O, but prevents a mutable filename from silently referring to different data.
The species-specific InterPro file is a derived input: its own hash identifies
what was analyzed, while `derived_from` records the upstream mapping source.

The inputs recorded are exactly those the selected ontology declares in
`src/ontology_registry.py` (`needs`, plus `hierarchy_needs` when
`--enable-true-path` is set), so every registered ontology is covered and a
newly added one is covered without editing the manifest code.

## Validate the software checkout

Run these commands independently from a clean checkout:

```bash
uv run ruff check src/ tests/ config/ scripts/ validation/ \
    run_dcgo_human.py extract_human_interpro.py
uv run ruff format --check
uv run pytest
```

- [ ] Record pass/fail, test count, operating system, and Python version.
- [ ] Retain CI links for the exact commit.
- [ ] If a release artifact is used, install that artifact in a fresh environment
      and smoke-test `dcgo --help` rather than relying on an editable checkout.

## Archive a result

- [ ] Archive the completed manifest with the result tables and figures.
- [ ] Archive raw-input hashes and release headers even if raw files cannot be
      redistributed.
- [ ] Archive the environment lock file and analysis scripts.
- [ ] Give the archive an immutable version or DOI.
- [ ] Cite dcGO using `CITATION.cff` and cite GO, GOA, InterPro, UniProt, Expasy,
      Reactome, ChEBI, or whichever other data sources the run actually used.
- [ ] Document exclusions, evidence filters, ontology propagation, score
      transfer, FDR threshold, and any post hoc analyses in the manuscript.

## Known limitations

- Most configured download URLs point to mutable `current_release` locations.
  The manifest makes the consumed file identifiable, but it does not make that
  upstream file permanently retrievable. Publication runs should archive or use
  dated source releases (`scripts/download_data.py --goa-archive` fetches dated
  GOA snapshots).
- Not every source embeds a machine-readable release identifier — the UniProt
  flat file and the Reactome relations file do not. Missing `release_metadata`
  means "not present or not recognized", not "current". Do not infer a release
  from file modification time.
- `analysis.thresholds.min_proteins_per_association` is `null` on purpose: the
  pipeline keeps an association on FDR significance alone and applies no
  minimum-support or effect-size floor. Introducing one is an open item in
  `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md`.
- The manifest covers the supported `run_dcgo_human.py` analysis path. The
  downstream tools that consume its output — `scripts/rank_surprising_associations.py`
  and the `validation/` benchmarks — do not yet emit equivalent manifests, so
  their commands and inputs must be archived separately.
