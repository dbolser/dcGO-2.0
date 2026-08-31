# dcGO-2.0 — Status

*2026-08-29. Compiled from the repo, VALIDATION_PLAN.md, and session history.*

## Where we are

The mandate is essentially built. 35 ontologies are registered (34 production
layers), the path is species-parameterised, and the paper-parity machinery
(`--propagate-annotations --enable-relative-inference --enable-true-path
--min-ic 1`) all landed in PRs #61–#70. The full production matrix ran on
2026-08-18: **63/63 cells, 0 failures** (`results/production/`, 2.7 GB,
uncommitted), with surprise rankings for 59 of them. CI is green with 963
tests. A manuscript draft exists in `paper/` with an evidence ledger.

## Headline numbers

- Human GO: 165,687 significant associations (baseline); 96,419 (paper-parity);
  allspecies 2.9M / 5.6M.
- Temporal benchmark (2021→2026): dcGO beats the naive baseline on F_max on
  informative terms (MF IC≥2: 0.365 vs 0.053) and clears the random null in
  every cell — but naive wins AUPRC at IC≥0 in all aspects. Report both.
- All-species background wins 8/9 F_max and 9/9 AUPRC held-out; 9/9 and 9/9
  when trained on experimental evidence only.
- Supra-domain associations predict future curation at **12.5×** enrichment.
  The surprise *ranking* is no better than q-value (98% of scores tie at 0).
- Show-and-tell emergent hits: SH3+PDZ → social behaviour, kinase+SAM →
  Eph–Ephrin, RRM+DAZ, tankyrase/PARP, srGAP.

## Waiting on you

1. **PR #72** (manuscript ontology update, 10 commits) — deliberately unmerged
   until you've read it.
2. **Credentials**: SNOMED CT and MedDRA (licences), OMIM `genemap2.txt`
   (registration key). MAxO also needs Monarch's `maxo-annotations` file —
   ontology is on disk, annotations aren't.
3. **Decide what to keep**: `results/production/` (2.7 GB) and
   `data/ACQUISITION_MATRIX.md` (marked "do not commit") are untracked.

## Next work, in order

1. **Gene / UniRef50 collapse.** Model-organism and allspecies surprise
   rankings are pseudo-replicated (isoform / ortholog stacks) — no claims from
   those layers until collapsed. Human layers are verified clean.
2. **Re-run manuscript evaluations §3.2–3.7 post-#67** (regulates-edge fix) —
   those blocks are marked PROVISIONAL.
3. **Re-run the §4 ablation post-#67.** The "True Path worse 12/12" number
   predates two fixes and conflates filter with propagation — don't cite it.
4. **HPO paper-parity collapse (996 → 38)** is driven by relative inference,
   not the IC floor. Specificity vs over-conservatism is decidable only by the
   post-#67 temporal evaluation.
5. **elim-style decorrelation.** Even at IC≥1, ~50% of associations still sit
   on an ancestor chain.
6. **§6 reproducibility** — the only VALIDATION_PLAN section still fully open:
   pin dated dataset versions, manifests for the surprise/validation tools,
   one-command reproduction.

## Standing caveats

- Allspecies support is inflated ~2.44× by orthology; half its associations
  rest on ≤3 UniRef50 clusters. Don't quote its counts with FDR unqualified.
- The `manual` allspecies universe is 75.8% projected annotation.
- The 24 "Error" strings in paper-parity logs are benign
  (`InsufficientBackgroundError` counts, handled conservatively).

## Housekeeping

- `TODO.md` (reviewed 2026-08-11) and `FUTURE_WORK.md` (says 21 ontologies;
  it's 35) both predate the sprint — refresh or retire.
- Five `agent/*` branches still local, two with live worktrees — prune after
  #72 lands.
