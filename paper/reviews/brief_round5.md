You are an expert peer reviewer for a computational-biology paper. This is review ROUND 5 of a revised manuscript. Read-only: do not edit any files in the repository.

Files (checkout /home/danbolser/.t3/worktrees/dcGO-2.0/t3-a5ab0cd0):
- Manuscript: paper/knowledge_projection.md (Figure 1: paper/figures/projection_precision_yield.png)
- Protocol with pre-registration and dated addenda: KNOWLEDGE_PROJECTION.md
- Result data behind every table: validation/projection/*/
  - test split: final_v2, final_geneknown_v2, final_htfilter_v2, final_strictnk_v2 (these are the runs the tables use; final_rescored etc. are the earlier re-runs with identical metrics), final_humanonly, final_assoc, final_curves_dense
  - final/ also holds training_sizes.tsv, curation_lag.tsv, run_facts.tsv, reliability_by_support.tsv, reliability_by_origin.tsv; final_humanonly/ and dev2018_frozen/ hold reliability tables too; final_rescored/cluster_sizes.tsv; final/background_origin_ablation.tsv, final/background_origin_curves.tsv; reliability_contrasts.tsv in final/, final_humanonly/, dev2018_frozen/
  - dev split: dev2018_frozen, dev2018_humanonly, dev2018_skipmissing, dev2018_pass*, cv2018_quantile, dev2018_assoc
- Table generator: paper/make_projection_tables.py
- Code: src/knowledge_projection.py, validation/projection_benchmark.py, validation/knowledge_yield.py

Rules: Do NOT penalise the paper for analyses it has not done if they are clearly, coherently and accurately described in Section 7 (Future work). DO penalise an obvious follow-up or limitation that is missing. Be specific and fair; do not invent problems. Spot-check as many numbers against the TSVs as you reasonably can.
ROUND 5. Round-4 scores (clarity/coherence/accuracy): Gemini 10/10/10, Codex 9/8/8, Claude 8/8/8; statistics specialist 8/10.

Round-4 criticisms and what was done (verify each; do not assume):
- The stale Discussion sentence ("mainly because ...") is replaced; Discussion, abstract and Section 5.3 now agree.
- Change-of-sign claims now name the cells that actually flipped (BP sequence neighbours and domain + sequence; CC domain + sequence, full method and literature neighbours), and note that the BP domain-unit loss appeared on both splits. MF literature neighbours are described as unaffected. The reason the background was frozen is quoted from protocol §1.6.
- New paired intervals (paper/background_origin_ablation.py; paper/reliability_by_support.py -> reliability_contrasts.tsv; Supplementary Table S4): removal vs human-only is −226 [−360, 100], so the paper now says it cannot tell how much the other-species-only units explain; the BP AUPRC change from the background is −0.013 [−0.033, 0.011], so the collapse is a top-of-list (KY50) effect; origin and support contrasts have paired intervals; the support paragraph is replaced by a calibration-gap summary (observed − mean score by support band, both training sets, both splits).
- The estimator's contribution on old dcGO's own evidence (human-only domain units vs old dcGO) is now the first point of the Discussion; Future Work item 11 is one-directional.
- Main-table differences and AUPRC ratios now all use the reference re-chosen in each resample (new diff_reselect columns; validation/projection_benchmark.py); the fixed-reference KY50 ratio is kept as the pre-registered column. Old dcGO's fixed-cutoff precision now has an interval (MF 0.55 [0.43, 0.69]).
- CAFA comparison corrected (three concrete differences), "headline" wording, the domain-only units departure (no supra-domains) in 4.7, the literature channel's training set, Table 4's "baseline ≈ 0" rule, the 14-protein count, garbled sentences, strict-cohort notes, the Table 8 tie-break by AUPRC.
- Three detail tables moved to a supplement (S1 support reliability, S2 companion measures, S3 no IC floor; S4 new). Main text has Tables 1–14.
- Figure 1 adds the human-only domain-unit curve (dotted).
- Future Work (31 items, four headings) adds: re-testing post-hoc claims on a fresh split, coverage of the ratio-to-max interval, gene study-level stratification and intermediate IC floors, background error analysis and pre-registering the background choice, observed-vs-mean-score by support, other embedding models and neighbour settings.

Judge this version on its own merits, against the same rules.
