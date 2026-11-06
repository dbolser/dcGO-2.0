ROUND 7. Round-6 scores (clarity/coherence/accuracy): Codex 9/9/9, Gemini 10/10/9, Claude 8/8/9; statistics specialist 8/10.

Round-6 criticisms and what was done (verify each; do not assume):
- Codex: thin-support over-confidence is now "in all four by point estimate", with no intervals for the gap differences stated.
- Gemini: the human-only bootstrap and naive-vs-old rows are now in the data appendix.
- Claude: InterPro2GO's near-50% blocks (49.8% BP, 49.0% MF) and its KY50 of 204 in the gene-known BP cohort are now reported (Section 5.4) and tied to the domain-only weakness (Discussion); the Discussion opening calls the three items "ways the method differs", and the background is discussed once in full; paired intervals for all three domain-unit background rows are cited; the 0-of-50 cell cites the TSV; the Figure 1 caption explains hollow markers and the invisible BP domain-unit KY50; the lower-bound result is stated exactly (BP +0.3%/-8.5%, MF lower at both); "one got worse" is replaced; the tie rule says reverse alphabetical name; a repeated multiplicity sentence is removed; Section 4.6 states how the review changes moved development numbers (656->681, 370->328, 232->196).
- Statistics: the threshold-transfer conclusion is corrected (domain units did transfer in MF and matched in CC); "no detectable difference" wording with pooled dev values; unpaired comparisons are labelled as such; "unstable" is marked as a point-estimate conclusion.
- Future Work (33 items) adds Monte Carlo error and KY50 concentration on few proteins, paired background intervals for all arms, graded comparators and InterPro2GO as a prior for domain units, and sequence embeddings for unreviewed background proteins.

Judge this version on its own merits, against the same rules.
