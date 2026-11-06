**SCORES**
- **Clarity 9/10**: The manuscript is exceptionally well-written, with precise language, clear definitions of novel metrics (like KY50), and a logical flow that makes complex methodological choices easy to understand.
- **Coherence 10/10**: The narrative is highly logical, building a compelling argument from the flaws of the old method through the development of the new one, with the abstract, results, and discussion perfectly aligned.
- **Accuracy 9/10**: The numbers in the text and tables match the data appendix almost perfectly, demonstrating a highly reproducible and rigorous analytical pipeline.

**ISSUES**
1. **Severity: minor**. Location: Section 5.1, "A lower credible bound". What is wrong: The text quotes a KY50 of "155,010", but the data appendix (`cv2018_quantile/projection_metrics.tsv`) reports this value as `1.55e+05` (155,000). Concrete fix: Change "155,010" to "155,000" to match the precision available in the data appendix, or update the appendix to output exact integers.
2. **Severity: minor**. Location: Section 4.6, "An existence rule" and "The same-gene rule". What is wrong: The exact numbers of proteins removed by the existence rule (20 and 13) and the same-gene rule (261 and 135) cannot be verified from the provided data appendix. Concrete fix: Add these exclusion counts to a summary file in the data appendix so they can be independently checked.

**MISSING FUTURE WORK**
- **Structural embeddings:** Given the explosion of predicted protein structures (e.g., AlphaFold DB) and structure-based language models (e.g., Foldseek, ProstT5, ESMFold) between 2021 and 2026, using structural embeddings as an additional or alternative nearest-neighbour channel is a glaringly obvious follow-up. Structure is more conserved than sequence and often provides stronger functional signals for unannotated proteins.

**WHAT WOULD RAISE EACH SCORE TO 9 OR 10**
- **Clarity**: To reach 10/10, ensure that every single number cited in the text (such as the exclusion counts in Section 4.6) has a direct, traceable counterpart in the data appendix.
- **Coherence**: Already at 10/10.
- **Accuracy**: To reach 10/10, fix the minor discrepancy regarding the `1.55e+05` vs `155,010` value and ensure the exclusion counts are verifiable.