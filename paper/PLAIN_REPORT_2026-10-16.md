# dcGO-2.0 in plain language

*What we built, what we tested it on, and what we found. Written 2026-10-16 so that someone with no background can follow it.*

## 1. The problem in one paragraph

Proteins do jobs in cells. Scientists write those jobs down using a shared dictionary called the Gene Ontology (GO). Each entry in the dictionary is a "term", such as "ATP binding" or "cell division". Most proteins in most species have never been studied directly, so nobody has written down their jobs. We want a computer method that guesses a protein's jobs from its parts.

## 2. The key idea: proteins are built from domains

A protein is a long chain folded into a few compact, reusable units called domains. Think of domains as Lego bricks. The same brick turns up in thousands of different proteins across all of life. If a brick is almost always found in proteins that do a particular job, then a new protein carrying that brick probably does that job too.

The original dcGO method (Fang and Gough, 2013) turned this into a statistical test. For every brick and every job, it asks: "Among proteins whose jobs we know, do the ones carrying this brick do this job more often than chance would predict?" If yes, and strongly, the brick and the job are linked. We call such a link an **association**. The method also tests combinations of two or three bricks, in case a job appears only when bricks are used together.

dcGO-2.0 is our rebuild of that method on today's databases.

## 3. The inputs

- **Which bricks each protein has.** From InterPro, a public database that has already scanned every known protein and listed its domains. We never scan sequences ourselves. We can name a brick in three ways: as an InterPro entry, as a SUPERFAMILY structural family (what the original dcGO used), or as a Pfam family.
- **Which jobs each protein has.** From the GO annotation database for human proteins. We only use jobs that a human curator wrote down, not ones a computer guessed.
- **The dictionary itself.** The GO file, which also says which terms are more general versions of which others ("ATP binding" is a kind of "nucleotide binding").

The numbers below are for the training set used in the test of section 6: one species (human), GO annotations as of April 2021 (GOA release 205, curator-made only, so no IEA or ND codes) and InterPro release 85.0 from the same month. A protein enters the statistics only if it has at least one domain *and* at least one job, so the final universe is smaller than either input file.

| Data type | Count | Source |
|---|---|---|
| Species | 1 (human) | UniProt |
| Proteins with at least one curated GO job | 18,713 | GOA release 205 (April 2021) |
| Proteins with at least one InterPro domain | 19,394 | InterPro 85.0 (April 2021), human subset |
| Proteins in the final universe (both) | 18,422 | intersection |
| Distinct InterPro domain entries | 17,061 | InterPro 85.0 |
| Protein-to-domain annotations | 81,314 pairs (232,836 signature rows) | InterPro 85.0 |
| Distinct GO terms used | 16,059 | GOA release 205 |
| Protein-to-GO annotations | 220,978 protein–term pairs (542,872 evidence rows) | GOA release 205 |
| GO dictionary | 2026-06-15 release of go-basic.obo | Gene Ontology |

The full production run (section 8) uses the current releases instead: GOA and InterPro as downloaded in July 2026.

## 4. The statistics, briefly

For each brick–job pair we build a 2 × 2 count: proteins with the brick and the job, with the brick only, with the job only, and with neither. Fisher's exact test turns that into a p-value. Because we test hundreds of millions of pairs, we correct for multiple testing and keep only pairs at a 1% false discovery rate.

The original method adds a second, stricter test. Instead of comparing against all proteins, it compares against only the proteins that already have the job's more general parent term. This asks: "Is this brick linked to the *specific* job, or just to the broad category it sits in?" We call this the **relative test**. The method keeps whichever of the two tests is weaker.

Each kept association gets a strength score called the **h-score**, which measures how far the observed count is above what chance would give.

To guess the jobs of a new protein, we add up the h-scores of all its bricks for each job, then rescale within that protein so its best guess scores 1 and its worst scores 0. This is the **p-score**.

## 5. What we fixed this week

A review of our code against the 2013 paper found three real mistakes.

1. **The h-score was the wrong formula.** We had used a transformed p-value that maxes out, so 37% of associations shared the top score. Now it is the paper's formula.
2. **The p-score double-counted.** When a brick had several specific associations under one general term, each one was added to the general term. Now each brick contributes once per term.
3. **Brick combinations were built wrongly.** We joined neighbouring database rows, but InterPro lists overlapping entries for the same region, so 40% of our "combinations" were the same brick listed twice. Now each protein is treated as a bag of distinct bricks and we test every pair and triple.

## 6. The test: predict the future

The honest way to test a predictor is to train it on the past and check it against what was discovered later.

- **Training.** Human job annotations as they stood in April 2021, and the InterPro release from the same month.
- **Test.** Annotations added by June 2026.

We score in two ways, because they answer different questions.

**Association-level.** Curators also maintain a hand-made list of brick-to-job links, called InterPro2GO. We take the links curators *added* between 2021 and 2026 (5,534 of them) and ask what fraction of our 2021 associations are among them. This tests the associations directly. It is strict, because many of our links may be right but not yet curated.

**Protein-level (the CAFA way).** We take human proteins that had no known job in 2021 but gained one by 2026 (318 for "biological process", 173 for "molecular function" after removing uninformative terms). For each, we predict jobs with the p-score and compare against the jobs they really gained. Two numbers come out: **sensitivity** (of the true jobs, how many did we predict?) and **PPV** (of the jobs we predicted, how many were true?). F_max is the best balance of the two over all possible cutoffs.

We also compare with a "naive" predictor that guesses the most common jobs for every protein. Any method worth having must beat it.

## 7. What we found

**Both methods beat naive on molecular function** by a wide margin (F_max 0.33 and 0.30 versus 0.08). On biological process the margin is small and only the simpler configuration clears it.

**The relative test is a trade.** With it, more of the associations were later curated (7.4% of them, against 5.6% without it, and 9.0% against 6.9% when both are limited to their best 4,024). But it keeps far fewer associations, so it recovers less than half as many of the newly curated links (7.4% versus 16.6%). At the protein level this costs F_max about 0.03 in both aspects. In short: fewer predictions, more of which curators have since confirmed. Because curation is incomplete, an association not yet curated is not shown to be wrong, so this is evidence of higher reliability, not proof of it.

**Adding ancestor terms beforehand helps.** Giving the test every protein's general terms as well as its specific ones improves both scores in every configuration.

**Brick combinations add nothing to the protein-level score.** They matter for a different question: which combinations predict a job that none of the single bricks does. Those combinations are 8 to 15 times more likely than chance to have their prediction confirmed by 2026.

**About half the test proteins get no prediction at all**, because none of their bricks had any association in 2021. That caps sensitivity before method quality even enters.

## 8. Which configuration is "best" is a judgement call

The two scores disagree about the relative test. The association-level score favours keeping it, because a larger share of its predictions has since been curated. The protein-level score says drop it, because you find more. Neither is wrong. Which matters more depends on what the predictions are for: if each one will be checked by an expensive experiment, reliability wins; if they feed a cheap downstream screen, coverage wins.

We ran the full production set (34 ontologies, three ways of naming bricks) with the relative test on, matching the original method. Switching it off is one flag.

## 9. Sources

- Fang H, Gough J. A domain-centric solution to functional genomics via dcGO Predictor. BMC Bioinformatics 2013.
- Fang H, Gough J. dcGO: database of domain-centric ontologies. Nucleic Acids Research 2013.
- Bao L et al. The dcGO domain-centric ontology database in 2023. Journal of Molecular Biology 2023.
- Zhou N et al. The CAFA challenge reports improved protein function prediction. Genome Biology 2019.
- This project's review (`REVIEW_2026-10-15.md`) and manuscript draft.
