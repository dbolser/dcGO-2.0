You are an expert peer reviewer for a computational-biology paper. This is review ROUND 3 of a revised manuscript. Read-only: do not edit any files in the repository.

Files (checkout /home/danbolser/.t3/worktrees/dcGO-2.0/t3-a5ab0cd0):
- Manuscript: paper/knowledge_projection.md (Figure 1: paper/figures/projection_precision_yield.png)
- Protocol with pre-registration and dated addenda: KNOWLEDGE_PROJECTION.md
- Result data behind every table: validation/projection/*/
  - test split: final_rescored, final_geneknown_rescored, final_htfilter_rescored, final_strictnk, final_humanonly, final_assoc, final_curves_dense
  - final/ also holds training_sizes.tsv, curation_lag.tsv, run_facts.tsv, reliability_by_support.tsv, reliability_by_origin.tsv; final_humanonly/ and dev2018_frozen/ hold reliability tables too; final_rescored/cluster_sizes.tsv
  - dev split: dev2018_frozen, dev2018_humanonly, dev2018_skipmissing, dev2018_pass*, cv2018_quantile, dev2018_assoc
- Table generator: paper/make_projection_tables.py
- Code: src/knowledge_projection.py, validation/projection_benchmark.py, validation/knowledge_yield.py

Rules: Do NOT penalise the paper for analyses it has not done if they are clearly, coherently and accurately described in Section 7 (Future work). DO penalise an obvious follow-up or limitation that is missing. Be specific and fair; do not invent problems. Spot-check as many numbers against the TSVs as you reasonably can.
ROUND 3. Round-2 scores (clarity/coherence/accuracy): Gemini 9/9/9, Codex 9/9/8, Claude 8/7/8; statistics specialist 7/10.

Round-2 criticisms and what was done (verify each; do not assume):
- BP domain-unit collapse was blamed on thinly supported units. A new carrier-origin analysis (Table 6, dev and test) shows the main cause on the test split is units carried only by other species' proteins, which take over the top of the BP list. Human-only training gives BP KY50 241. Thin support is now a secondary, "consistent with winner's curse" effect (Table 7: support bands for dev, test human-only, and test with background). Both tables have protein-bootstrap intervals, protein counts, and omit cells under 20 predictions. Abstract, Discussion and Future Work (items 3 and 4) were updated to match.
- Discussion attribution corrected: the background's effects on each arm; the cutoff transfer is a property of the combined score, not of θ alone; channel and background ablations were measured, the estimator's own contribution was not.
- The domain-only arm's weakness is stated (CC AUPRC no better than naive; naive beats it without an IC floor).
- The strict cohort is now "no curated annotation in any aspect", with a note that it is not identical to CAFA no-knowledge.
- Table 11 (sensitivity) now has re-chosen intervals, KY50 differences, and each cohort's own best old variant; the CC gene-known sentence was corrected.
- Every interval column says whether its reference is fixed or re-chosen. The pre-registered domain + literature arm is in Tables 3 and 10, and the change of headline arm is the first item of Section 4.7.
- Section 4.4 now says: the share of resamples with ratio ≥ 10 counts all-zero-baseline resamples; it is not a probability; the re-chosen interval leans conservative; bootstrap coverage for KY50 is not guaranteed.
- Table 9 shows old dcGO's dev-best variant at its own dev cutoff; the text says the dev cutoff is a selected point.
- The "71 CC proteins" was wrong; cluster sizes are now computed (cluster_sizes.tsv): in CC, 51 of 322 proteins share a paper, largest cluster 17.
- Title no longer says "frequency-scale"; the abstract no longer says old predictions "cannot be pooled"; the q-value wording is now "largest −log₁₀ q-value"; the figure caption describes the axis range; p-values were recomputed (3.5 × 10⁻⁴⁴, 1.6 × 10⁻⁷).
- Curation forecasting softened to "may".
- Future work adds: interval-coverage checks, interaction networks (STRING), an audit of top-ranked misses, a test of literature overlap with evidence papers, and the cross-validation evidence against a uniformly stronger prior.
- KNOWLEDGE_PROJECTION.md has a dated corrections section (§4) for its stale numbers.
- Tables were renumbered (there are now 17).

Judge this version on its own merits, against the same rules.
