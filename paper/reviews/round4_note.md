ROUND 4. Round-3 scores (clarity/coherence/accuracy): Gemini 10/10/9, Codex 9/9/9, Claude 8/8/8; statistics specialist 8/10.

Round-3 criticisms and what was done (verify each; do not assume):
- "Main cause" of the BP domain-unit collapse was too strong. A new ablation (paper/background_origin_ablation.py, final/background_origin_ablation.tsv) removes every unit with no human carrier: BP KY50 rises only 1 -> 15 [0, 313], against 241 with human training (paired difference 240 [51, 403]). Section 5.3, the abstract, the Discussion and Future Work now say the background is the cause, other-species-only units are part of it (about half the top predictions, mostly wrong), and the background also lowers the scores of the predictions human training ranks highest (241 at >= 0.9 at 50-73% correct, against 97 with 18 correct).
- The background's change of sign between splits is now stated (Section 5.3, and a new Discussion paragraph "The background"), with the post-hoc human-only sequence CC KY50 of 380 beside the IEA comparison.
- The thin-support claim is now band-specific (holds at 0.7-0.9, not at 0.5-0.7 or >= 0.9), with the theta ceiling for 1-3 carriers stated as a condition on the term's training frequency.
- Reliability tables: cells need >= 10 distinct proteins; all-right or all-wrong cells show counts with no interval; "none" marks empty cells.
- Both pre-registered arms (domain only 0.0x/1.5x, + literature 7.1x/1.3x) are in Section 4.7 and Table 10. Section 4.4 names the pre-registered contrasts. The cohort name in the protocol ("CAFA no-knowledge") is addressed in 4.7.
- KY50 threshold rule: the largest yield, ties to the highest threshold (Sections 4.4 and 5.4). 0/0 resamples are dropped, and positive/0 counts as reaching 10.
- Cluster sizes for all three aspects. Strict-cohort CAFA wording fixed, with notes that more channels did not always help there.
- Learned-combiner comparison numbers corrected (227/272, or 274/280 with InterPro entries alone). Coverage wording corrected. Table 9: old dcGO's dev-best variant, the BP tie-break and the q-value floor explained. The cutoff-transfer claim restated (full and two-channel combinations transfer approximately; domain units alone do not).
- Abstract: the BP ratio interval also given as unbounded with the fixed reference, plus the KY50 difference.
- Future Work is grouped under four headings (30 items), with new items: background per channel and aspect with an error analysis, old dcGO as a channel, where the BP gain lies, PubTator time-matching.
- The protocol corrections section (§4) is extended.

Judge this version on its own merits, against the same rules.
