ROUND 2. This is a revised manuscript. Round-1 scores were: Gemini 9/10/9, Codex 8/8/7, Claude 7/6/6 (clarity/coherence/accuracy), statistics specialist 6/10.

Main round-1 criticisms and what was done (check that each is actually fixed; do not assume it):
- Calibration was overclaimed, including in the title. Now: scores are "comparable across proteins", and calibration is tested only via the winner's-curse table (precision by score band x support band, new Table 6), which shows low-support units are over-confident.
- The claim about 70% precision was wrong: the best old variant reaches KY70 = 20 in BP, not 0. Corrected.
- Ratios were quoted against one baseline for all metrics. Now each metric is compared with the old variant that is best on that metric (BP AUPRC 4.5x, not 6.2x). The bootstrap also re-chooses the best old variant inside each resample (all 28 variants), reports P(ratio >= 10), the KY50 difference with its interval, the ratio to the naive baseline, a cluster bootstrap over proteins linked by shared evidence papers, and precision at the frozen development cutoff with intervals.
- "reproduces exactly" -> "within 0.004".
- Curation-lag analysis recomputed on the scored truth (IC >= 2). The earlier MF 5% figure was an artefact.
- Several exact counts were corrected (supra-domain CV count 16,028; 152/486; IEA CC precision 0.56; InterPro 69.0 dated 3 days after the dev t0).
- The high-throughput filter is now described as uninformative (it removes 0/1/0 truth pairs).
- Added: strict no-knowledge cohort (no annotation in any aspect), human-only background on the test split (Table 7), old dcGO "max" transfer alongside p-score, and a full list of deviations and post-hoc analyses (Section 4.7).
- The association-level endpoint is scored on current GO ids, because most "added" InterPro2GO pairs were GO-id renames.
- Structural embeddings and other follow-ups were added to Future Work (Section 7).
- Tables were renumbered into the order they first appear.

Judge this version on its own merits, against the same rules.
