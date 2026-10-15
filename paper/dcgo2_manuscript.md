# dcGO-2.0: domain-centric ontology inference reimplemented on InterPro, with a temporal evaluation

Authors: [to be completed]
Affiliations: [to be completed]
Corresponding author: [to be completed]

*Draft of 2026-10-15. Text in square brackets is a placeholder for a value or
passage that depends on the frozen method and the final runs. Methods describe
the method as implemented at the time of the final runs; where a choice is
still open it is marked.*

## Abstract

**Background.** Protein domains are the units at which molecular function is
conserved and recombined. The dcGO method (Fang and Gough, 2013) associates
ontology terms with domains and with combinations of domains by enrichment
tests over annotated proteins, and transfers those associations to any
protein that carries the domain. The original implementation used SCOP
superfamilies assigned by SUPERFAMILY and predates a decade of growth in
InterPro, UniProt and the Gene Ontology.

**Results.** We present dcGO-2.0, an open reimplementation of the method on
precomputed InterPro assignments. Domains may be keyed as InterPro entries,
SUPERFAMILY signatures or Pfam signatures. Each protein is treated as a set of
distinct domains, and every combination of two or three of them is tested as
a supra-domain. Associations are tested by Fisher's exact test against the
whole protein universe and against the proteins annotated to the term's
direct parents, with the weaker of the two p-values corrected by the
Benjamini–Hochberg procedure. Annotation input is generalised behind a single
interface, so [N] vocabularies beyond the Gene Ontology can be used without
changing the statistics. We evaluate the method on a temporal split of human
annotation (training on GOA release 205 of April 2021 and InterPro 85.0,
testing on annotation added by June 2026) with two instruments: an
association-level test against domain–term pairs curated into InterPro2GO
after the training date, and a protein-level CAFA-style benchmark. [One
sentence per primary result.] We compare the association set with the
published dcGO release on the shared SCOP-superfamily space.

**Conclusions.** [To be written from the frozen run.]

**Keywords.** protein domains; Gene Ontology; functional annotation; InterPro;
supra-domains; temporal benchmark.

## 1. Introduction

Function prediction for uncharacterised proteins remains a standing problem in
computational biology. Community assessments have shown that a predictor which
reproduces the annotation frequency distribution can score well on
protein-centric metrics without making an informative statement [6]. The Gene
Ontology [4] is the dominant target vocabulary, and its term frequencies are
highly uneven.

Domain-centric annotation is one response. Domains are the units at which
structure and function are conserved, so an association learned at the domain
level transfers to every protein carrying that domain, including proteins with
no close annotated homologue. The dcGO method of Fang and Gough [1,2]
formalised this. For each domain and term, it tests whether proteins carrying
the domain are enriched for the term against all annotated proteins, and
against the proteins annotated to the term's direct parents. It corrects for
multiple testing, propagates associations up the ontology, and extends the
feature space from single domains to combinations of domains, so that a
function carried only by an architecture can be captured. Predictions are
transferred to a target protein by summing association scores over the
domains it carries.

Three things motivate a reimplementation. InterPro [7] now integrates
SUPERFAMILY, Pfam and other member databases and distributes precomputed
protein-to-domain assignments, so no sequence scanning is required. The
evaluation standard has moved: a claim about a function predictor now needs a
temporal held-out design, a frequency baseline, and explicit treatment of
low-information terms. And the statistical machinery is indifferent to the
vocabulary that supplies the terms, which has not been exercised
systematically.

This paper describes the reimplementation, states where it departs from the
original, and reports its performance on a temporal benchmark with the
evaluation design fixed before the final runs.

## 2. Methods

### 2.1 Inputs

The pipeline consumes three public inputs and performs no sequence analysis.

- Domain assignments: precomputed InterPro matches (`protein2ipr`) keyed by
  UniProt accession, with match coordinates. For the temporal evaluation the
  training architectures are taken from InterPro 85.0 (April 2021), the last
  release before the training annotation snapshot.
- Annotations: for the Gene Ontology, the species GOA file in GAF 2.2 format.
  By default annotations with evidence code IEA or ND are excluded. Other
  vocabularies are described in §2.6.
- Ontology structure: `go-basic.obo` for GO; other vocabularies supply their
  own hierarchy.

Every run writes a manifest recording the SHA-256 digest and release header
of each input, the Git state, the dependency lock hash, the command line and
every threshold.

### 2.2 Domain features and supra-domains

For each protein the parser records the set of distinct domains it carries,
under one of three keys: the integrated InterPro entry, the SUPERFAMILY
signature (whose accession is the SCOP superfamily identifier used by the
original dcGO), or the Pfam signature. A supra-domain is any combination of
two or three distinct domains present on the same protein, keyed in a
canonical order. Positional adjacency is not required. This differs from the
original method, which required successive domains in the SUPERFAMILY
architecture; the reasons are that InterPro entries nest and overlap, so
adjacency of matches does not correspond to adjacency of domains, and that a
set-based definition is the same under every domain key.

Each feature is recorded once per protein regardless of how many matches
support it.

### 2.3 Statistical inference

For each feature f and term t, a 2 × 2 contingency table is formed over the
protein universe, defined as the proteins with at least one domain and at
least one annotation. A one-sided Fisher exact test gives the overall
p-value. Only feature–term pairs that co-occur in at least one protein are
evaluated; every other pair has p = 1 under the one-sided test and
contributes to the hypothesis count only.

The relative test uses the same table restricted to the background of
proteins annotated to at least one direct parent of t (the union of the
parents' protein sets, the quantity the original paper's Figure 1 calls
N_pa). A term with no parent, or whose parental background has fewer than
[3] proteins, is not tested relatively and passes on the overall test alone.

The two p-values are combined as their maximum, and the Benjamini–Hochberg
procedure [5] is applied to the combined p-values, with single domains and
supra-domains corrected as separate hypothesis families. Pairs with adjusted
p below [0.01 | 0.001; decision pending] are retained.

Each retained association carries an h-score, defined as in the original
method: the hypergeometric standard score (observed count minus its
expectation, divided by its standard deviation) of the overall table and of
the relative table, and the smaller of the two.

Input annotations are propagated to all ancestor terms before the tables are
built. GO identifiers that have been merged or replaced are remapped to their
current identifier; identifiers with no live replacement are dropped and
counted.

### 2.4 Propagation and reporting floors

Significant associations are propagated to ancestor terms, an ancestor
carrying the minimum adjusted p and the maximum h-score of its significant
descendants. Propagation follows `is_a` and `part_of` edges only.

Two reporting floors may be applied after correction, so that neither alters
the hypothesis family. A minimum-support floor drops associations supported by
fewer than N proteins. An information-content floor drops associations to
terms whose annotation-frequency information content, IC(t) = −log₂ P(t)
estimated from the propagated annotation map, is below X bits. The
information floor exists because a parentless term passes the relative test
by default, so without it the ontology roots dominate the output. [State the
floor used in the final runs.]

### 2.5 Prediction for a target protein

A target protein's predictions are formed from the associations of the
features it carries. For each feature and term, the feature contributes one
score, the h-score of its propagated association to that term. Scores are
summed across features, and the sums are scaled to [0, 1] within the protein
by min–max normalisation. This is the p-score of the original method.

### 2.6 Annotation sources beyond the Gene Ontology

Annotations enter through a single interface that yields a map from protein
to term set, so the inference engine never sees anything vocabulary-specific.
A registry holds one entry per vocabulary with its annotation source, its
hierarchy (or none) and its required input files. [N] vocabularies are
registered: the Gene Ontology; Enzyme Commission from Expasy; UniProt-native
layers from the Swiss-Prot flat file (cross-reference databases, keywords,
subcellular location, ChEBI ligands and cofactors, Rhea reactions); disease
layers re-keyed from OMIM and Orphanet onto the Disease Ontology or Mondo;
gene-keyed human layers (HPO, SynGO, GWAS Catalog traits, HPA cell-type
expression); and model-organism phenotype and anatomy layers (mouse, worm,
zebrafish, fly). Gene-keyed sources are re-keyed to UniProt accessions at
parse time; unmapped identifiers are dropped and counted, and one-to-many
identifiers credit all of their accessions. Model-organism associations are
learned on the organism's own proteins; because domains are species-agnostic,
the association applies to any protein carrying the domain.

[Supplementary Table S1: registry, source files, hierarchy mechanism and
identifier-mapping coverage per vocabulary.]

### 2.7 Evaluation design

All evaluation uses one temporal split of human data. Training uses GOA
release 205 (April 2021) and InterPro 85.0 (April 2021). Testing uses
annotation added by GOA of June 2026 and InterPro2GO of [release]. One GO
release ([date]) is used for training, prediction and scoring.

The primary endpoint, fixed before the final runs, is association-level:
the precision of a configuration's single-domain associations against
domain–term pairs added to the curated InterPro2GO mapping between the two
InterPro releases, at a fixed number of associations, for biological process
and molecular function. Both sides are propagated to ancestor closure, and
the comparison is restricted to domains present in both mappings. Uncertainty
is a paired bootstrap resampling domains. Supra-domain associations are
assessed by their enrichment for annotations gained by 2026 on carrier
proteins that lacked the term in 2021, relative to the term's own acquisition
rate.

The secondary endpoint is protein-level. For each GO aspect, a protein enters
the benchmark if it had no non-IEA annotation in that aspect in 2021 and
gained experimental annotation by 2026. Its truth set is the propagated 2026
experimental term set, roots removed. F_max is computed following
CAFA-evaluator conventions: precision averaged over proteins with at least one
prediction at the threshold, recall over all benchmark proteins, thresholds
at 0.001 and above. Area under the precision–recall curve is reported without
an origin anchor. An information-content floor of 2 bits, estimated from the
2021 training annotations, is applied identically to truth and predictions;
proteins left without truth leave the cohort, and cohort sizes are reported.
Baselines are the CAFA naive predictor, scored both at full coverage and
restricted to the proteins the method covers, and a random-domain null in
which each domain is reassigned another domain's term set under [200] seeded
permutations. Uncertainty is a paired bootstrap resampling proteins.

Configuration choices (domain key, transfer rule, information floor,
threshold) are made on the primary endpoint over a nested development split
[GOA 205 → 215] and evaluated once on the full interval.

### 2.8 Comparison with the published dcGO

The published dcGO tables (Domain2GO, SP2GO, April 2016) are keyed by SCOP
superfamily identifier. With the SUPERFAMILY domain key, our associations lie
in the same identifier space, and the fraction of our associations present
in the published set is reported. Recall against the published set is not
interpretable because the published set was learned on all sequenced genomes
and ours on human proteins; the fraction of published pairs with no
co-occurring human protein is reported instead.

## 3. Results

### 3.1 Association sets

[Table 1: for each domain key (InterPro, SUPERFAMILY, Pfam) and for the
baseline and full configurations: proteins in the universe, single domains,
supra-domains, hypotheses, associations retained, associations after floors.]

### 3.2 Component ablation on the primary endpoint

[Table 2: nine configurations (single domains; supra-domains; with and
without input propagation, relative test and output propagation) scored on
the association-level temporal endpoint, SUPERFAMILY key, with paired
bootstrap intervals. Text: which components change precision at fixed
budget, and in which direction.]

### 3.3 Protein-level benchmark

[Table 3: F_max, AUPRC and cohort size per aspect at IC ≥ 0 and IC ≥ 2 for
the chosen configuration, the naive baseline at full and matched coverage,
and the permutation null, with paired bootstrap intervals. Text: the result
and the sensitivity and positive predictive value at the F_max threshold.]

### 3.4 Supra-domain associations

[Enrichment of supra-domain associations for later annotation, with
interval; the number of associations emergent beyond their constituent
domains; two or three examples.]

### 3.5 Comparison with the published dcGO

[Precision against the published set under the SUPERFAMILY key; fraction of
published pairs unreachable from a human universe.]

### 3.6 Other vocabularies

[Association counts per vocabulary for the chosen configuration, labelled as
hypothesis-set sizes, with the temporal enrichment for the vocabularies that
can be split in time.]

## 4. Discussion

[To be written from the results. Points to cover: what the association-level
and protein-level endpoints agree and disagree on; the effect of the relative
test, expressed as sensitivity against positive predictive value; the
set-based supra-domain definition and its consequences; the relation to the
original method, including the h-score, threshold and domain universe; the
dependence of the Benjamini–Hochberg family on nested terms and nested
features; limitations: human training universe, single time interval,
look-ahead through the GO release, no independent predictor as comparator.]

## 5. Conclusions

[To be written.]

## 6. Data and code availability

The implementation, evaluation code and manuscript are in the dcGO-2.0
repository [URL, licence, archived DOI]. All inputs are public: InterPro
`protein2ipr` and InterPro2GO (releases [85.0 and current]); UniProt-GOA
(release 205 and [current]); the Gene Ontology (`go-basic.obo`, [date]); the
published dcGO tables (supfam.org, April 2016, archived with the release);
and the annotation sources listed in Supplementary Table S1. Run manifests
with input and output digests are distributed with the results. [Resource
requirements and one-command reproduction instructions.]

## References

1. Fang H, Gough J. A domain-centric solution to functional genomics via dcGO
   Predictor. BMC Bioinformatics. 2013;14(Suppl 3):S9.
2. Fang H, Gough J. dcGO: database of domain-centric ontologies on functions,
   phenotypes, diseases and more. Nucleic Acids Res. 2013;41(D1):D536–D544.
3. Bao L, et al. The dcGO domain-centric ontology database in 2023: new
   website and extended annotations for protein structural domains. J Mol
   Biol. 2023;435(14):168093. [author list to be verified]
4. The Gene Ontology Consortium. [current GO reference]. Nucleic Acids Res.
5. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical
   and powerful approach to multiple testing. J R Stat Soc B.
   1995;57(1):289–300.
6. Zhou N, et al. The CAFA challenge reports improved protein function
   prediction and new functional annotations for hundreds of genes through
   experimental screens. Genome Biol. 2019;20:244.
7. Paysan-Lafosse T, et al. InterPro in 2022. Nucleic Acids Res.
   2023;51(D1):D418–D427.
8. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in
   2025. Nucleic Acids Res. 2025;53(D1):D609–D617.
9. [CAFA-evaluator reference]
10. [InterPro2GO reference]
11. [SUPERFAMILY reference]
12. [Pfam reference]
13. [HPO, SynGO, MGI, WormBase, ZFIN, FlyBase, Mondo, Reactome references as
    used]
