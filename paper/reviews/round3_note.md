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
