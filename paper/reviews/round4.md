**SCORES**
- **Clarity 10/10**: The manuscript is exceptionally well-written, with precise definitions, plain language, and clear explanations of complex metrics like KY50 and the projection estimator.
- **Coherence 10/10**: The narrative flows logically from the limitations of the old method to the design of the new one, followed by a rigorous, transparent evaluation that directly addresses the pre-registered goals and unexpected findings.
- **Accuracy 10/10**: Every checked claim, number, and interval in the text and tables matches the provided data appendix perfectly, including the complex post-hoc ablations and bootstrap intervals.

**ISSUES**
None. The authors have meticulously addressed all previous criticisms, corrected the numbers, clarified the text, and ensured the data appendix perfectly supports the manuscript. The handling of the background ablation and the association-level scoring artefact is particularly commendable.

**MISSING FUTURE WORK**
The Future Work section (Section 7) is extraordinarily comprehensive, covering 30 well-reasoned items. One minor addition could be:
- **Alternative embedding models and hyperparameters.** The current method relies on ProtT5 and MedEmbed-small with fixed kNN parameters (k=25, T=0.1). Future work could explore whether larger or newer models (e.g., ESM-2 for sequences, PubMedBERT for literature) or tuned retrieval hyperparameters improve the neighbour channels' precision.

**WHAT WOULD RAISE EACH SCORE TO 9 OR 10**
All scores are already at 10. The manuscript is flawless, transparent, and ready for submission.