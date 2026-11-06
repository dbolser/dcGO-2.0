ROUND 6. Round-5 scores (clarity/coherence/accuracy): Codex 9/9/9, Claude 8/8/9, Gemini 10/10/10 in round 4 (its round-5 call timed out); statistics specialist 8/10.

Round-5 criticisms and what was done (verify each; do not assume):
- Discussion and Future Work item 13 now use Section 5.3's wording on the removal ablation (too uncertain to say how much those units explain).
- The Discussion opening is restructured into three sources (estimator; neighbour channels incl. coverage and threshold transfer; background). The estimator claim now has paired intervals from a new human-only bootstrap run (final_humanonly_boot): AUPRC 2.5x [2.0, 3.3], 1.5x [1.2, 1.8], 2.8x [1.8, 4.7]; KY50 clear only in BP (+187 [29, 315]); MF +47 [-15, 121]; CC -19 [-69, 135].
- The support claim was re-checked with a new dev human-only reliability run: thin units are less precise in all four BP runs and more over-confident in three; the earlier "mainly with the background" claim was wrong and is withdrawn; whether the background changes the effect is stated as untested. Table S4 has a dev human-only row and applies the 10-protein rule.
- The origin contrast is quoted with its pooled value; AUPRC wording "no detectable effect"; "post hoc" added in the abstract; MF "raised or kept (point estimates)".
- Reference-variant ties and the tie rule are disclosed (Section 4.1); naive-vs-old BP AUPRC now "by point estimate (0.90-2.00)"; CAFA codes corrected (TAS/IC removed from the stricter list); the literature-training sentence moved; combiner sentence rewritten; lower-bound result stated (8.5% lower); reproduction tolerances harmonised; uncertainty paragraph de-duplicated.
- Supplementary tables are in order S1-S4; the header note says S1-S4.
- Figure 1: the human-only curve has short dashes and hollow markers.
- Future Work (32 items) adds a human base-rate prior with all-species counts, term-centric results and NOT-qualified negatives, a logistic reliability model, the support-by-background interaction among claims to re-test, and separating human from non-human background annotations.

Judge this version on its own merits, against the same rules.
