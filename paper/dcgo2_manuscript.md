# dcGO-2.0: a reimplementation of domain-centric ontology inference, with a retrospective assessment of what its predictions are worth

**Authors.** *[To be completed.]*

**Affiliations.** *[To be completed.]*

**Corresponding author.** *[To be completed.]*

---

> **Draft conventions.** Bracketed identifiers such as `[C2a]` point to rows in
> `paper/EVIDENCE_LEDGER.md`, which records the file and locator for every number
> below. They are a drafting aid and should be stripped before submission. Every
> number in this manuscript has a ledger row; no number appears here that could
> not be traced to a repository artifact. Passages marked **PROVISIONAL** rest on
> results that are known to be affected by an identified defect and must be
> regenerated before submission.
>
> **Era marking (added 2026-08-17).** The pipeline now restricts GO annotation
> propagation to `is_a` and `part_of` edges; previously it also traversed the
> DAG's 7,799 regulates-family edges `[K3]`. Every number in this draft that
> passed through GO propagation — on the inference side or inside the evaluation
> code — was produced before that fix and is marked **PROVISIONAL —
> pre-regulates-fix era**. Run manifests record the policy as
> `analysis.ontology.propagation_relations`; a manifest without that key (or an
> artifact with no manifest at all, which includes every run in blocks A–H of
> the ledger — the 2026-08-18 production matrix, block M, is manifest-carrying
> and post-fix) identifies a pre-fix artifact `[K12]`. Numbers that never touch
> GO propagation
> (dataset scale counts, download and survey figures, non-GO hierarchies read by
> the `is_a`/`part_of`-only OBO reader) are not marked.

---

## Abstract

**Background.** Protein domains are the units at which much of molecular
function is conserved, and domain-centric ontology inference — associating an
ontology term with a domain or a domain combination rather than with a whole
protein — is an established route to transferable functional annotation. The
original dcGO method (Fang & Gough, 2013) established the statistical framework
but predates a decade of growth in InterPro, UniProt and the Gene Ontology, and
its accompanying evaluation was not designed to separate genuine predictive
signal from recovery of the annotation base rate.

**Results.** We present dcGO-2.0, an open reimplementation that associates
InterPro entries and contiguous InterPro combinations ("supra-domains") with
ontology terms by Fisher's exact tests with Benjamini–Hochberg control, and that
generalises the annotation input behind a single interface: thirty-three
vocabularies beyond the Gene Ontology are registered — UniProt-native layers,
Enzyme Commission, disease, expression, and human and model-organism phenotype
ontologies whose gene-keyed annotations are re-keyed to protein accessions
under an audited mapping policy — without touching the statistics `[M8]`. The method
now reproduces the original pipeline's structure: the relative
(parental-background) test is combined with the overall test before FDR
correction, input annotations are propagated by the True Path Rule, GO
propagation is restricted to `is_a`/`part_of`, and an information-content
floor removes vacuous near-root terms `[K3, K11]`. On current human data the
baseline configuration performs 1.69 × 10⁹ tests and reports 165,687
associations at FDR < 0.01; the paper-parity configuration reports 96,419, of
which 30,302 are single-domain under the information floor — both from
manifest-carrying post-fix runs `[M1, M2]`. We assess the associations in
three retrospective settings; all three were evaluated against the 2021
training run, predate the propagation fix, and are marked provisional in the
text pending regeneration. (i) Against the curated InterPro2GO map, treated
strictly as an incomplete positive reference, the association set recovers 64.7%
of curated pairs on shared domains after propagation `[B1]`. (ii) In a CAFA-style
no-knowledge benchmark trained on a 2021 GO annotation snapshot and scored
against experimentally supported 2026 annotations, dcGO exceeds a naive
frequency baseline on F_max for biological process (0.248 vs 0.115) and cellular
component (0.380 vs 0.343) but **loses** for molecular function (0.360 vs 0.464)
`[C2a, C2e, C2i]`; under an information-content floor that removes near-universal
terms, dcGO leads on F_max in all three aspects `[C2b–C2l]`. The naive baseline
nevertheless attains a higher area under the precision–recall curve than dcGO in
all three aspects without the filter, and in cellular component up to a floor of
four bits `[C3a–C3g]`. (iii) Supra-domain associations found in 2021 anticipate
curation added by 2026 at 12.5-fold the terms' own acquisition rates (2,181 hits
on 170,416 predictions; bootstrap CI 10.6–14.1) `[E2]`; however, a proposed
"surprise" ranking of emergent combinations shows **no demonstrated advantage**
over ranking by the dcGO q-value, with a paired bootstrap interval spanning zero
at every prediction budget `[E9, E10]`. Applying the unchanged statistics to the expansion vocabularies —
HPO, SynGO, disease and trait re-keys, expression layers, and mouse, worm,
zebrafish and fly phenotype layers learned on each organism's own proteins,
with domains transferring species-agnostically — yields per-layer association
sets now measured by a 63-cell, manifest-carrying production matrix
`[M0, M7]`; none of these layers is yet validated.

**Conclusions.** Domain-derived associations contain predictive signal for later
human GO annotations and outperform simple baselines in several retrospective
settings, particularly for higher-information GO terms. The evidence assembled
here does not establish general superiority or score calibration; the
information-content filter that carries the headline was not pre-specified; and
the relative-inference layer, though now structurally faithful to the original
method, does not yet achieve its stated specificity purpose `[K5]`. We report
the analysis, its negative results, and the open items required before any
stronger claim would be defensible.

**Keywords.** protein domains; Gene Ontology; functional annotation; InterPro;
supra-domains; phenotype ontologies; retrospective benchmark.

---

## 1. Introduction

Function prediction for uncharacterised proteins remains one of the standing
problems in computational biology, and the community's own assessments have
repeatedly shown how easily apparent performance can be an artefact of the
annotation base rate rather than of biological insight [6]. The Gene Ontology [4]
is the dominant target vocabulary, and its highly uneven term frequencies mean
that a predictor which simply reproduces the annotation frequency distribution
can score respectably on protein-centric metrics without making a single
informative statement.

Domain-centric annotation is one response. Because domains are the units at
which structure and function are conserved and recombined, an association learnt
at the domain level transfers to every protein carrying that domain, including
proteins with no close homologue among annotated sequences. The dcGO method of
Fang and Gough [1,2] formalised this: for each (domain, term) pair, test whether
proteins carrying the domain are enriched for the term against a background of
all analysable proteins; correct for multiple testing; propagate up the ontology
by the True Path Rule; and, distinctively, extend the feature space from single
domains to *supra-domains* — contiguous combinations — so that a function
emerging only from an architecture can be captured. The original work also
included a second, *relative* test against a background restricted to proteins
already annotated to all direct parent terms, which forces an association to beat
its own generic parent, and transferred predictions to proteins through a
per-target normalised score `[G4, G5]`.

Three things motivate a reimplementation now. First, the underlying resources
have changed substantially: InterPro [7] integrates SUPERFAMILY and Pfam among
many member databases and offers precomputed protein-to-domain assignments at a
scale that removes the need to run sequence scanning at all, while UniProt [8]
and the GO annotation corpus have grown and been re-curated. Second, and more
importantly, the evaluation standard has moved. A modern claim about a
function-prediction method requires a held-out design that is temporal rather
than random, baselines that include the frequency-only predictor, and explicit
treatment of the fact that low-information terms dominate the curated record.
Third, the domain-centric statistical machinery is indifferent to which
vocabulary supplies the terms, and that generality has never been systematically
exercised.

This paper describes dcGO-2.0 and — this is the substance of the contribution —
reports what its output does and does not support. We deliberately foreground
three negative or qualified results: the unfiltered molecular-function comparison
that dcGO loses; the area-under-curve comparison that the naive baseline wins;
and the emergent-combination ranking that fails to beat plain significance. We do
this because the alternative — reporting only the filtered F_max table — would
misrepresent the evidence.

---

## 2. Methods

### 2.1 Inputs

The pipeline consumes three public inputs and performs no sequence analysis of
its own.

- **Domain assignments.** Precomputed InterPro matches from `protein2ipr.dat.gz`
  [7], keyed by UniProt accession, including per-match sequence coordinates. A
  one-off extraction step streams the full file and retains the subset of
  proteins for the species under study.
- **Annotations.** For the Gene Ontology, the species GO annotation file in GAF
  2.2 format, with evidence-code filtering (default: all non-IEA, i.e. `manual`).
  For other vocabularies, see §2.6–2.7.
- **Ontology structure.** `go-basic.obo` for GO propagation; other vocabularies
  supply their own hierarchy in one of five forms (§2.6).

All three are downloaded by a single script from URLs held in a central
configuration module. Non-GO vocabularies declare their own annotation and
hierarchy inputs in the registry (§2.6), so a run fails before the statistics
if an input is missing. Retrieved releases are still referenced through
mutable "current" URLs `[H8]`, but every run now writes a machine-readable
manifest (`run_manifest_<ontology>.json`) recording a SHA-256 digest and any
release header for each input, the Git state, the dependency lock hash, the
full command line and every threshold `[K12]`. The 2026-08-18 production
matrix (§3.9–3.10) is fully manifest-carrying; **every number from earlier
runs — including everything the evaluations of §3.2–3.7 rest on — predates
the manifest machinery** and is not covered by it. See §5.

### 2.2 Domain features and supra-domains

For each protein the parser produces a domain architecture: the set of InterPro
entries matched, and the ordered sequence of matched regions. A **supra-domain**
is a contiguous run of 2 or 3 entries in that order; supra-domains are treated as
additional columns of the feature matrix and tested identically to single
domains. In the human 2021 run this expands 19,230 single-domain features
`[A3]` into 102,206 total features `[A2]`, contributing 230,275 single-domain and
405,928 supra-domain incidences `[A6]`.

A domain–protein pair is recorded once per supporting member signature in
`protein2ipr`, so the raw incidence matrix contains duplicate entries; the
implementation collapses the matrix to binary presence/absence before any
counting. (Failing to do so was one of four correctness defects found and fixed
during development; see §5.)

### 2.3 Statistical inference

For each (feature *f*, term *t*) pair a 2 × 2 contingency table is formed over the
protein universe — proteins carrying *f* and annotated *t*; carrying *f* only;
annotated *t* only; neither — using sparse matrix products, and a one-sided
Fisher exact test is evaluated with a vectorised Cython implementation.
*p*-values are adjusted by the Benjamini–Hochberg procedure [5], with single
domains and supra-domains corrected as **separate hypothesis families**, each
against its own dense pair count — a supra-domain is not an exchangeable
sibling of its constituent domains — and pairs with adjusted *p* < 0.01 are
retained; run manifests record the two families and their BH thresholds
`[M9]`. (Runs made before 2026-08-05 used a single pooled family; the 2021
training run of §2.8 is one of them.) A hypergeometric association score
rescaled to 1–100 is also emitted.

Two properties of this design are load-bearing for the interpretation of the
results and are stated here rather than in the discussion. (i) The hypothesis
family is strongly dependent — GO terms are nested, supra-domains contain their
constituents, and domains co-occur — so BH control over it is a convention rather
than a demonstrated guarantee `[H15]`. (ii) An earlier version of the method
included an optional "empirical-Bayes shrinkage" step that geometrically
interpolated a sparse supra-domain's observed *p*-value toward its
constituents' with a hand-set decay constant. The resulting quantities were
never shown to be valid *p*-values under any null — enabling the step took
FDR < 0.01 rejections from 163,277 to 463,924 (+184%) by pulling thin evidence
toward its well-supported parts — and the component ablation found no effect on
prediction quality in any of 12 aspect × IC cells. **The step has been removed
from the codebase** (2026-08-05) `[K7]`; it was already off in every run
reported here `[H5]`, so no number in this draft depends on it.

The current human baseline run — manifest-carrying, at the production commit —
reports 165,687 associations at FDR < 0.01 (44,453 single-domain / 121,234
supra-domain) from 1,690,803,963 tests over 18,908 proteins and 16,389 GO
terms `[M1]`; the paper-parity configuration of §2.5 reports 96,419, of which
30,302 are single-domain under the `--min-ic 1` floor `[M2]`. The 2021
training run used for all temporal evaluation reports 164,549 `[A7]` from
1,640,917,330 tests `[A5]` over 18,382 proteins and 16,055 GO terms
`[A1, A4]`; that run predates the manifest machinery and the per-family BH
split, and its era is discussed in the draft conventions.

By default no minimum-support or effect-size floor is applied beyond FDR
significance. Two opt-in **reporting** filters exist, both applied *after* the
BH correction so the hypothesis family is never altered: `--min-support N`
drops associations supported by fewer than *N* proteins, and `--min-ic X`
drops associations to terms whose annotation-frequency information content —
IC(*t*) = −log₂ P(*t*), estimated from the propagated annotation map and
exported as an `ic` column — falls below *X* bits `[K11]`. The IC floor exists
for a specific reason connected to relative inference and is motivated in
§2.5. Odds ratios are emitted uncorrected, so they print as `0` when the *d*
cell is empty and as `inf` when *b·c* = 0 `[H21]`; the ranking signal is the
*p*-value, not the odds ratio.

### 2.4 True Path Rule

The True Path Rule appears in three distinct places, and the method now keeps
them separate.

**Input-side propagation** (`--propagate-annotations`) closes the input
protein → term map over the hierarchy before any counting, which the original
paper states as part of its design. Terms the hierarchy does not contain are
handled explicitly rather than silently: GO `alt_id`s are remapped to their
primary identifiers, and terms still unknown after the remap are dropped from
the tested universe, with every case counted and recorded in the run manifest
`[K11]`.

**Output-side propagation** (`--enable-true-path`) adds ancestor associations
to the significant set. It applies to every registered vocabulary that has a
hierarchy, through a shared engine that computes a transitive ancestor closure
from a child→parent map, and it now **fails explicitly** for a vocabulary with
no hierarchy instead of degrading silently. It is propagation and nothing
else: the parental-background filter that used to run under the same flag for
GO is the paper's separate *relative inference* step and has its own flag
(§2.5).

**Evaluation-side propagation** is applied within the evaluation code, to both
predictions and truth, wherever the metric requires it (§2.8).

For GO, all three paths propagate over **`is_a` and `part_of` edges only** —
the relations GO's own annotation-propagation rules license. Until 2026-08-17
the propagation graph also traversed the ontology's 7,799 regulates-family
edges, which put proteins annotated to "negative regulation of X" into X's
propagated background `[K3]`. **This changes every previously quoted number
that passed through GO propagation**, which is why the affected results below
carry era markers. Run manifests record the edge policy as
`analysis.ontology.propagation_relations` `[K12]`. The runs reported in this
draft were made before the fix; output-side propagation was disabled in them,
but evaluation-side propagation was not.

### 2.5 The relative (parental-background) test

The original method combined the overall enrichment test with a *relative* test
whose background is restricted to proteins annotated to the direct parent terms
of *t*, keeping the weaker of the two `[G4]`. This is the mechanism that is
supposed to force an association to be more specific than its parent. In
dcGO-2.0 it is now **folded into the inference** (`--enable-relative-inference`)
rather than applied as a post-hoc filter, matching the published pipeline's
structure `[K11]`:

- the relative *p*-value is computed for every candidate pair before any
  correction, against a background formed as the **union of the direct
  parents'** protein sets (the paper's Figure-1 *N_pa*);
- BH corrects `max(overall_p, relative_p)` — an intersection–union statistic,
  and therefore a valid *p*-value without further adjustment — and the
  reported h-score is `min(overall, relative)`;
- the input annotation map is True-Path propagated (`--propagate-annotations`),
  as the paper states;
- the test is available for every registered vocabulary with a hierarchy
  (19 of the 28 registry entries `[K2]`), not only GO.

One structural defect remains and is handled by a reporting floor rather than
by the test itself: a **parentless term skips the relative test entirely** and
passes on the overall inference alone, so the DAG roots and the near-universal
band just below them dominate the raw output. The `--min-ic` floor of §2.3
exists for exactly this case: the GO aspect roots sit at 0.09–0.17 bits of
annotation-frequency information `[K6]`, so any floor clear of that band
removes them. The floor is applied to every deliverable, including the
propagated annotations file, so output-side propagation cannot re-derive a
floored-away ancestor. What the floor does and does not fix is measured in
§3.8.

The **paper-parity configuration** referred to throughout this draft is
`--propagate-annotations --enable-relative-inference --enable-true-path`, with
`--min-ic 1` as the recommended reporting floor; on human GO single domains it
reports 30,655 significant associations at FDR < 0.01 before the floor `[K4]`,
and the manifest-carrying production run reproduces the floored count exactly
— 30,302 single-domain associations under `--min-ic 1` `[M2]`, matching the
sweep's floor-1 row `[K5]`. The measurement in §3.4 predates this design and
evaluated the earlier post-hoc variant.

### 2.6 The annotation seam and the vocabulary registry

Annotations enter through a single interface that yields a `{protein → {term}}`
map, so the inference engine never sees anything ontology-specific. A registry
holds one entry per vocabulary: a source factory, an ancestors factory and a
direct-parents factory (or `None` when the vocabulary has no hierarchy), and
the input files required, so a run fails early on a missing input.
**Thirty-five `--ontology` keys are registered** `[M8]` (nineteen at the time
this draft was first compiled `[A12]`, twenty-eight at the 2026-08-17 update
`[K1]`), in seven kinds:

- the **Gene Ontology**, from the species GAF file;
- **Enzyme Commission**, from Expasy `enzyme.dat` (UniProt-keyed already);
- **UniProt-native layers** from three parts of the Swiss-Prot flat file —
  twelve `DR` cross-reference vocabularies (Reactome, OMIM phenotypes,
  Orphanet, TCDB, MEROPS, CAZy, UniPathway, ComplexPortal, DrugBank, Pharos,
  CD-CODE, plus `KW` keywords) and four layers curated into the entry body
  (`CC SUBCELLULAR LOCATION` mapped to `SL-` terms, `FT /ligand_id` ChEBI
  ligands, `CC COFACTOR` ChEBI cofactors, `CC CATALYTIC ACTIVITY` Rhea
  reactions);
- **term re-keys**: the OMIM and Orphanet disease layers re-keyed at parse
  time onto Disease Ontology terms (`doid`, `orphanet_doid`) or Mondo terms
  (`mondo`, `orphanet_mondo`), which is what gives the disease curation a
  hierarchy at all;
- **gene-keyed human layers** (HPO, SynGO, GWAS-Catalog EFO traits, HPA
  cell-type expression) and **model-organism layers** (MP, WBPhenotype, ZFA,
  FBcv, FBbt phenotypes; WBbt anatomy expression), whose identifier mapping
  is the subject of §2.7;
- **cross-reference-chained cancer layers** (`ncit`, `oncotree`): CIViC's
  CC0 gene → disease evidence carried through DOID's NCI Thesaurus
  cross-references, and onward to OncoTree codes via shared NCI/UMLS ids —
  two mapping hops, every stage counted `[L10]`;
- a generic `xref` escape hatch reaching any other `DR` database by name.

UniProt-derived vocabularies were selected by a survey of all ~150 `DR`
databases in the human Swiss-Prot subset, using proteins-per-term as the
decisive statistic: a resource that mirrors the protein 1:1 (ratio ≈ 1) cannot
support an enrichment test at all. Registered vocabularies range from 28.4
proteins/term (keywords) through 5.0 (Reactome) to 1.5 (ComplexPortal) and 0.4
(DrugBank) `[F16]`; ~1:1 accession mirrors and domain databases (which would be
circular) are deliberately excluded `[F15, F16]`.

Hierarchies come from five mechanisms: OBO graphs (GO via its own processor;
ChEBI, DOID, Mondo, HP, MP, WBPhenotype, WBbt, ZFA, FBcv, FBbt, EFO and NCIt
via a light `is_a`/`part_of` reader), structure implicit in the identifier
(EC, TCDB, MEROPS, CAZy), companion hierarchy files (Reactome, UniProt
keywords, subcellular locations, OncoTree's JSON parent tree), a hierarchy
sheet shipped inside the annotation release itself (SynGO), and none at all —
for which `--enable-true-path` and `--enable-relative-inference` fail
explicitly rather than degrading silently. Twenty-five of the thirty-five
entries have a hierarchy `[M8]` (nineteen of twenty-eight at the previous
update `[K2]`).

> **Resolved defect (fixed 2026-08-04).** An earlier version of the
> UniProt-native sources parsed the entire Swiss-Prot flat file with no
> taxonomic restriction, and the contingency-table row space was built as the
> union of the annotation and domain maps, inflating the protein universe for
> those vocabularies toward all of Swiss-Prot (575,503 entries `[A11]`) rather
> than the ~20k human proteins. The annotation map is now restricted to the
> domain-annotated intersection, and every §3.7 row was regenerated after the
> fix; the superseded values are retained there for comparison. GO (sourced
> from the species-specific GAF) and EC were never affected. See `[F]` in the
> evidence ledger for the audit trail.

### 2.7 Gene-keyed annotation sources and identifier mapping

The statistics join domains and annotations on UniProt accessions
(`protein2ipr`'s key space), but the phenotype resources added in this round
are keyed by *gene*: HPO's `genes_to_phenotype.txt` by NCBI GeneID, SynGO's
bulk release by HGNC id, and the model-organism databases by their own gene
identifiers (MGI, WBGene, ZDB-GENE, FBal/FBgn). Each layer re-keys
gene → accession **at parse time**, before the statistics see anything, under
one audited policy: an unmapped gene id is dropped, counted and logged, never
silently discarded; a one-to-many id (real, for readthrough loci and
unresolved paralogs) credits *all* of its accessions, since choosing one
arbitrarily is not reproducible. Every count is exposed by the source and will
be recorded in run manifests `[L0]`.

The two human layers map through the Swiss-Prot flat file's own `DR GeneID`
and `DR HGNC` lines, so no separate id-mapping download is needed: HPO maps
5,196 of 5,274 gene ids (**98.5%**) and SynGO 1,787 of 1,789 (**99.9%**)
`[L1, L2]`. The model-organism layers instead translate through each
database's own mapping files — essential because model-organism proteomes are
mostly unreviewed, so the map must include TrEMBL — with coverage of **85.0%**
of MGI genes (14,033/16,509), **82.5%** of WBGene ids (8,704/10,550), and
**68.6%** of ZDB-GENE ids (5,139/7,487) `[L3, L4, L5]`. The fly layers map
**35.1%** (FBcv; 6,121/17,438) and **42.0%** (FBbt; 3,358/7,987) of FBgn ids
`[L6]` — low by design rather than by defect: FlyBase issues gene identifiers
to transgenic constructs and drivers that have no *Drosophila* protein
product, and those dominate the unmapped tail `[L6a]`.

Genotype-to-gene attribution is handled per database, and restrictively. MGI
rows are kept only for **single-gene genotypes** (a phenotype observed on an
*Atm*/*Rad50* double mutant cannot be attributed to either gene alone; the
file is already overwhelmingly single-gene, with 2 of 283,003 rows dropped)
`[L3]`. FlyBase rows are kept only for **single-allele genotypes**, which
retains 199,921 of 399,972 data rows (50.0%) `[L6]`. WormBase publishes
gene-level associations directly, with `NOT`-qualified rows dropped and
counted as negative evidence `[L4]`.

One layer deserves an explicit honesty note. ZFIN curates phenotypes as EQ
post-compositions rather than pre-composed phenotype-ontology identifiers, and
composing them to ZP terms is **not derivable from the released files**: the
annotation file carries no ZP ids, and the OBO edition of `zp.obo` carries
logical definitions for 0 of 43,521 non-obsolete terms (only ~2% have a
machine-readable EQ signature in a comment) `[L5]`. The zebrafish layer is
therefore an **affected-anatomy** layer — "proteins whose mutation produces an
abnormality of anatomical structure S are enriched for domain D", propagated
up the ZFA DAG — which loses only the quality dimension: all 169,887 rows in
the current file are tagged `abnormal` `[L5]`.

The model-organism layers are the original dcGO's central trick restated on
this infrastructure: the domain → phenotype association is learned on the
model organism's *own* protein universe (`--species mouse/worm/zebrafish/fly`,
with the domain universe built by the same species-parameterised extraction),
and because domains are species-agnostic the association annotates any protein
carrying the domain — including human ones.

### 2.8 Evaluation designs

Three evaluations are reported, measuring different things.

**(a) Coverage of a curated domain-level reference (domain-centric).** The
association set is compared with InterPro2GO, treated as an **incomplete,
positive-only** reference. Both sides are propagated to ancestor closure and the
comparison is restricted to domains present in both. The reported quantity is
*coverage* (recall of curated pairs). A "precision lower bound" is also emitted,
but non-recovered pairs are curation-gap candidates rather than demonstrated
errors, so it is a floor and not an estimate of precision `[B3, H16]`.

**(b) Temporal, protein-centric, CAFA-style benchmark.** Training uses GO
annotation release 205 (April 2021); testing uses the current release (June 2026)
`[C0a]`. Per GO aspect, a protein enters the benchmark only if it had **no
annotation known to training** in that aspect at t0 — the gate uses the same
non-IEA evidence space the pipeline trains on, so a computational label the model
had already seen cannot re-enter the held-out set — and gained experimental
annotation by t1. Its truth set is the **full** propagated t1 experimental term
set for that aspect, roots excluded `[C0b]`. Predictions are transferred by the
per-target p-score: the sum of association scores over the protein's features and
their propagated ancestors, min–max normalised within the protein `[C0c]`.
Baselines are CAFA `naive` (every protein predicted every term at its propagated
t0 frequency) `[C0d]` and a random-domain null in which each domain is reassigned
another domain's entire term set under a **single seeded permutation** — a null
for the transfer step only, since Fisher is not re-run `[C0e]`.

Metrics are F_max, S_min (with marginal information content from t0), and AUPRC,
computed per aspect. Precision is averaged over proteins with at least one
prediction at the threshold and recall over all benchmark proteins, following the
CAFA convention; prediction coverage is computed internally but is not reported
in the committed metrics `[C0h, H17]`. The threshold sweep uses 51 score
quantiles plus a predict-nothing sentinel, and AUPRC is a trapezoidal integral
over an upper envelope of precision as a function of recall — choices that have
**not** been verified against an independent CAFA evaluation implementation
`[C0g, H18]`.

An **information-content floor** is applied identically to truth and to every
method: terms with marginal IC = −log₂P(t) below a threshold are dropped, and
proteins left with no truth terms are dropped from the cohort `[C0f]`. This
filter was introduced during analysis and was **not pre-specified** `[H6]`; its
consequences for cohort composition are reported in §3.3.

**(c) Held-out test of emergent combinations.** For each supra-domain
association found at t0, its *predictions* are the proteins carrying the
combination that lacked the term at t0 (propagated, non-IEA); its *hits* are
those annotated with the term at t1 (propagated, experimental); and the control
is the term's own acquisition rate among all domain-carrying proteins that lacked
it at t0. Enrichment is the hit rate divided by that base rate, which removes the
tendency of popular terms to accumulate annotations regardless of any prediction.
Confidence intervals are percentile bootstraps resampling **associations** (the
conservative unit, since proteins within one association share an architecture)
`[E5]`. When two rankings are compared, the difference is estimated by a
**paired** bootstrap that re-ranks both ways inside each resample, because both
rankings act on one shared candidate pool and their independent intervals are
correlated `[E9]`.

### 2.9 The surprise score for emergent combinations

Ranking supra-domain associations by raw significance is dominated by three
artefacts: redundant InterPro signatures describing one region as if it were two
domains; restatements of what InterPro2GO already curates; and a tension whereby
genuinely novel combinations rest on very few proteins. We therefore score each
supra-domain association *S* → *t* as

```
surprise = −log10(q_emergence) × distinctness × novelty
```

where **emergence** is a one-sided binomial tail of the observed co-occurrence
count against a parts-only expectation — the maximum of a noisy-OR over the
constituents' individual rates, the best proper sub-combination's rate, and the
term's background rate, with rates shrunk toward the background by one
pseudo-observation — BH-corrected across all candidates as its own hypothesis
family; **distinctness** is 1 − the median largest pairwise overlap between the
constituents' matched regions, from the `protein2ipr` coordinates, with
candidates above 0.5 overlap discarded; and **novelty** is a fixed weight
(0.1 curated / 0.3 implied / 0.6 refines / 1.0 novel or no-reference) against
InterPro2GO `[D1, D2]`. The novelty weights are a prioritisation convention, not
an inference; the statistically meaningful quantity is `q_emergence`. Every
component is written to the output so results can be re-weighted without
recomputation.

---

## 3. Results

### 3.1 Scale of the association set

On the 2021 human training snapshot the pipeline tested 1.64 × 10⁹ feature–term
pairs and retained 164,549 at FDR < 0.01 `[A5, A7]` — about 0.01% of tests.
The current, manifest-carrying baseline run retains 165,687 (44,453
single-domain / 121,234 supra-domain) `[M1]`, and the paper-parity
configuration retains 96,419, of which 30,302 are single-domain under the
`--min-ic 1` floor `[M2]`. Feature space is dominated by combinations: 19,230
single domains expand to 102,206 features once contiguous pairs and triples
are included `[A2, A3]`.

> **Era status (updated 2026-08-18).** The current-release headline is now
> **resolved**: the production matrix regenerated it post-fix under a
> manifest (165,687, superseding the earlier manifest-less 165,823 `[A8]`).
> The 164,549 figure remains the **pre-fix, manifest-less 2021 training run**
> on which every temporal evaluation below rests; it stays quoted because it
> is what those evaluations used, and it carries their era marking.

We stress that "significant" here means only that: no minimum-support or
effect-size policy is applied, so the set includes associations resting on very
small contingency tables `[H19]`. The counts above should be read as the size of
a hypothesis set, not as a count of validated findings.

### 3.2 Coverage of the curated InterPro2GO reference

> **PROVISIONAL — pre-regulates-fix era.** Both sides of this comparison are
> propagated to ancestor closure through the GO processor, which at the time
> of these runs traversed the regulates-family edges `[K3]`. Every figure in
> this subsection must be regenerated from post-#67 runs.

At the FDR < 0.01 operating point the association set recovers **64.7%** of
curated InterPro2GO pairs — 30,673 of 47,393 propagated pairs on the domains
present in both sets `[B1, B2]`. Tightening the cutoff trades coverage for
selectivity in the expected direction: 29.3% (13,881 pairs) at *p* ≤ 10⁻¹⁰
`[B4]`. The loose end of the sweep is uninformative because the input is already
FDR-filtered, so every threshold from *p* ≤ 10⁻⁶ downward in stringency returns
the identical set `[B5]`.

This is a **coverage** result and nothing more. InterPro2GO is deliberately
incomplete and positive-only, so the 22.3% "precision lower bound" at the same
cutoff `[B3]` is a floor: the ~107,000 non-recovered pairs are candidates in a
curation gap, not demonstrated errors. This evaluation cannot establish precision
and is not independent validation `[H16]`.

### 3.3 Temporal benchmark: what the headline depends on

> **PROVISIONAL — pre-regulates-fix era.** The benchmark's no-knowledge gate,
> its truth sets, the naive baseline's term frequencies and the p-score
> transfer all propagate through the GO processor, which at the time of these
> runs traversed the regulates-family edges `[K3]`. Cohort sizes and every
> value in Tables 1 and 2 must be regenerated from post-#67 runs before
> submission; the qualitative reading below is retained because the defect
> applies identically to dcGO and to both baselines, but no number in this
> subsection should be quoted as current.

**Cohorts.** The no-knowledge benchmark comprises 324 biological-process, 418
molecular-function and 572 cellular-component proteins at IC ≥ 0 `[C1a]`. These
are small cohorts; no protein-level confidence intervals were computed for any
metric below `[H2]`.

**The information-content filter changes the cohort, not just the scoring.**
This must be stated before the results, because it conditions their
interpretation. Applying an IC floor drops truth terms and then drops proteins
whose truth becomes empty. Molecular function falls from 418 to **170** proteins
(−59%) at IC ≥ 2, cellular component from 572 to 405 (−29%), while biological
process is unchanged `[C1a, C1b, C1e]`. At IC ≥ 6 the cohorts are 289 / 145 / 154
`[C1d]`. The IC ≥ 0 and IC ≥ 2 rows therefore **do not measure the same
proteins**, and the comparison across floors is unpaired `[H7]`.

**F_max.** Table 1 gives F_max for dcGO and both baselines.

**Table 1.** Protein-centric F_max, 2021 → 2026 no-knowledge benchmark. Bold
marks the better of dcGO and naive. Values from
`validation/temporal_benchmark_metrics.tsv` `[C2a–C2l]`; the final column is
derived.

| Aspect | IC floor | *n* proteins | dcGO | naive | random | dcGO ÷ random |
|---|:--:|--:|--:|--:|--:|--:|
| BP | ≥ 0 | 324 | **0.248** | 0.115 | 0.158 | 1.6× |
| BP | ≥ 2 | 324 | **0.170** | 0.071 | 0.053 | 3.2× |
| BP | ≥ 4 | 318 | **0.115** | 0.031 | 0.019 | 6.1× |
| BP | ≥ 6 | 289 | **0.077** | 0.010 | 0.003 | 24× |
| MF | ≥ 0 | 418 | 0.360 | **0.464** | 0.262 | 1.4× |
| MF | ≥ 2 | 170 | **0.365** | 0.053 | 0.088 | 4.2× |
| MF | ≥ 4 | 162 | **0.337** | 0.045 | 0.072 | 4.7× |
| MF | ≥ 6 | 145 | **0.217** | 0.018 | 0.009 | 25× |
| CC | ≥ 0 | 572 | **0.380** | 0.343 | 0.291 | 1.3× |
| CC | ≥ 2 | 405 | **0.239** | 0.153 | 0.072 | 3.3× |
| CC | ≥ 4 | 252 | **0.134** | 0.099 | 0.031 | 4.3× |
| CC | ≥ 6 | 154 | **0.124** | 0.044 | 0.015 | 8.1× |

Read plainly, this table says three things.

1. **Without any filter, dcGO loses to the naive frequency baseline in molecular
   function** (0.360 vs 0.464) `[C2e]`, and wins in biological process and
   cellular component `[C2a, C2i]`. The margin in cellular component is modest
   (0.380 vs 0.343) and has no confidence interval.
2. **The claim that dcGO leads in all three aspects holds only under an IC floor
   that was not pre-specified** `[H6]`, and that floor simultaneously changes the
   molecular-function cohort by 59% `[C1e]`. We report the IC-filtered analysis
   because we regard it as scientifically motivated — a domain-to-function method
   should be judged on informative terms — but it is a post hoc secondary
   analysis, and we do not present it as the primary endpoint.
3. **The naive baseline's apparent strength is base-rate recovery.** Its F_max
   collapses toward the random null as informative terms are required: BP
   0.115 → 0.010, MF 0.464 → 0.018, CC 0.343 → 0.044 across IC 0 → 6 `[C2n]`.
   dcGO degrades far more gently, sitting 1.3–25× above the random-domain null
   throughout `[C2o]`. The molecular-function case is the sharpest: dcGO's F_max
   is essentially unchanged from IC ≥ 0 to IC ≥ 2 (0.360 → 0.365) while naive
   falls by an order of magnitude. The mechanism is that molecular-function truth
   is dominated by a single near-universal term, GO:0005515 *protein binding*,
   reported elsewhere in this project as 84.6% of human experimental
   molecular-function annotations — a figure we could trace only to project
   documentation and therefore mark **provisional** `[C6a]`.

**AUPRC tells a different story, and we report it.** Table 2 gives the area under
the precision–recall curve for the same runs.

**Table 2.** AUPRC, same benchmark. Bold marks the better of dcGO and naive
`[C3a–C3h]`.

| Aspect | IC ≥ 0 | IC ≥ 2 | IC ≥ 4 | IC ≥ 6 |
|---|---|---|---|---|
| BP dcGO / naive | 0.137 / **0.314** | **0.069** / 0.049 | **0.032** / 0.012 | **0.017** / 0.003 |
| MF dcGO / naive | 0.195 / **0.325** | **0.227** / 0.025 | **0.210** / 0.011 | **0.121** / 0.006 |
| CC dcGO / naive | 0.240 / **0.513** | 0.073 / **0.117** | 0.029 / **0.038** | **0.025** / 0.012 |

The naive baseline attains a **higher AUPRC than dcGO in all three aspects
without the IC filter**, and in cellular component it continues to do so at
IC ≥ 2 and IC ≥ 4, losing only at IC ≥ 6 `[C3i]`. This is a genuine counter-result
to the F_max table and is not reported in the project's existing summaries. Its
most likely explanation is structural rather than substantive — a predictor that
emits every term for every protein sweeps out a full recall range and so
integrates favourably, whereas dcGO's predictions are sparse and its curve
truncates — but that explanation is a hypothesis we have not tested, and the
AUPRC implementation itself uses a quantile sweep and an upper-envelope
trapezoidal rule that have not been checked against reference CAFA tooling
`[C0g, H18]`. Either way, "dcGO outperforms the naive baseline" is not a
statement the evidence supports without qualification.

**S_min.** On the information-content-weighted semantic distance, lower being
better, dcGO attains the best value of the three methods in all three aspects at
IC ≥ 0: BP 99.70 vs naive 104.86, MF 12.64 vs 14.25, CC 20.04 vs 21.34 `[C4a–C4c]`.
Notably this includes molecular function, the aspect dcGO loses on F_max. The
margins are ≤ 5% and carry no confidence intervals `[C4d]`, so we present this as
consistent with, rather than as independent confirmation of, the F_max picture.

**Restored method pieces.** Two components of the original method were absent
from an earlier version of this implementation and were restored: the per-target
p-score transfer and the relative (parental-background) test — the latter in its
earlier post-hoc form, since superseded by the in-inference design of §2.5. The
p-score is the configuration used throughout Tables 1 and 2 `[C5b]`. A
four-configuration comparison of the two pieces exists `[C5a, C5c, C5d]`, and
its benchmark outputs have since been committed (`validation/bench_A`–`bench_D`
`[H12]`), but no paired test or confidence interval was computed for any
pairwise difference and the runs are pre-regulates-fix era, so we mark that
comparison **provisional** and draw no conclusion about which component matters
more from protein-centric metrics alone.

### 3.4 Domain-centric effect of the relative test

> **PROVISIONAL — superseded machinery, pre-regulates-fix era.** This
> measurement evaluated the earlier *post-hoc* relative filter, which has
> since been replaced by the in-inference combination of §2.5, and its
> propagation ran over the pre-#67 edge set `[K3]`. It is retained as the
> only committed measurement of the relative test's domain-centric direction;
> it must be re-run under the current design (`--enable-relative-inference`)
> before any number here is quoted.

Scored directly against propagated InterPro2GO on the shared single-domain
space, the base association set gives coverage 0.631 at a precision floor of
0.218 (134,610 predicted pairs, 29,382 recovered); adding the relative
parental-background filter gives 0.430 coverage at a precision floor of 0.253
(69,206 pairs, 17,525 recovered) `[B7, B8]`. The filter thus raises the precision
floor by 0.035 while halving the association set `[B9]`, consistent with its
intended role of pruning generic, parent-driven associations. F1 is essentially
unchanged (0.324 vs 0.319).

Two caveats bound this. The reference is the **current** InterPro2GO rather than a
2021 snapshot, so this is not a temporal test on the domain side `[B10]`. And no
uncertainty was computed for the 0.218 → 0.253 difference `[B9]`, so it should be
read as a direction, not as a significant improvement.

### 3.5 Emergent domain combinations

> **Era note.** The candidate set derives from the current-release
> association run, a manifest-less pre-parity artifact (§3.1). The emergence
> statistic itself does not propagate, so the edge fix does not touch the
> arithmetic below, but the counts inherit the association set's era and will
> be regenerated with it.

Of the supra-domain associations in the current human GO run, 22,243 survive the
redundant-signature filter `[D3]` and **24** are emergent beyond their
constituents at FDR ≤ 0.05 under the binomial parts-baseline test `[D4]`. For
Enzyme Commission the corresponding figures are 1,401 retained and **1** emergent
`[D7, D8]`.

That the emergent count is two orders of magnitude smaller than the retained
candidate count is the substantive finding, and it is deflationary. The weaker
and more common criterion — "significant for the combination but for none of its
constituents" — is not equivalent to "the combination beats what its parts
predict", and the two diverge most sharply for Enzyme Commission, where the
protein universe is enzymes and a catalytic constituent usually already explains
the activity `[D12]`. Earlier counts of emergent associations produced under the
weaker criterion have been retracted within this project and are not repeated
here `[D13]`.

The 24 surviving GO cases are dominated by textbook multi-domain architectures
recovered without supervision — SH2 + protein-kinase-like → non-membrane-spanning
tyrosine kinase activity; PH + EF-hand pair → phosphatidylinositol
phospholipase C activity; BTB/POZ + C2H2 zinc finger → DNA-binding transcription
repressor activity; tyrosine-kinase catalytic + SAM → ephrin receptor signalling
`[D11]`. We treat this as a sanity check on the ranking rather than as a result:
recovering known architectures shows the score is not measuring noise, but it is
not evidence about the unknown cases. All figures in this paragraph derive from a
prose table and are marked **provisional** in the ledger pending field-level
re-verification `[D11]`.

Critically, every quantity in this subsection is computed on the same proteins
that produced the associations, so it measures internal consistency of the
evidence, not predictive power `[D14]`. That is what the next subsection tests.

### 3.6 Held-out test: the associations predict, the ranking does not

> **PROVISIONAL — pre-regulates-fix era.** The prediction and hit sets in
> this subsection are defined over propagated non-IEA and experimental
> annotation closures computed through the pre-#67 GO processor `[K3]`, and
> the underlying 2021 association set is a manifest-less pre-parity run. The
> verdicts below (the associations predict; the ranking shows no demonstrated
> advantage) are qualitative and paired, but every number must be regenerated
> from post-fix runs.

Applying the score to the 2021 association set yields 22,376 candidates, of which
10,136 leave at least one standing prediction `[E1]`.

**The associations anticipate curation.** Across the whole pool, 170,416
(protein, term) predictions produced 2,181 hits by 2026 — a 1.28% hit rate
against a base-rate expectation of 0.10%, i.e. **12.5-fold enrichment**
(bootstrap CI over associations: 10.6–14.1) `[E2]`. Roughly 175 hits would be
expected from the terms' own acquisition rates alone `[E4]`. This is the
strongest evidence in this paper: it is genuinely out-of-sample, the control
removes term popularity, and the effect is large.

**The ranking does not beat plain significance.** Comparing the surprise ranking
with ranking by the dcGO q-value requires matching on prediction budget rather
than on rank, because the two rankings expose very different numbers of
predictions per association (≈ 4.8 vs ≈ 73 at K = 100) `[E7]`. At matched
budgets the point estimates favour surprise at all three budgets — 15.5 vs 5.3 at
2,000 predictions, 21.2 vs 13.2 at 10,000, 11.6 vs 10.8 at 40,000 `[E8]` — but a
paired bootstrap of the *difference* puts zero inside every interval:
+10.19 [−87.28, +18.95]; +7.96 [−8.93, +7.06]; +0.82 [−1.01, +4.54] `[E9, E10]`.

**The honest conclusion is that the surprise score has no demonstrated ranking
advantage over the dcGO q-value.** We record two further observations. First, the
paired result is itself unstable: at the 10,000-prediction budget the point
estimate (+7.96) lies outside its own percentile interval and only 32% of
resamples favour surprise despite the positive point estimate `[E12]` — a
signature of a heavy-tailed statistic that we flag rather than interpret. Second,
comparing the two rankings' *independent* intervals would have looked
substantially more favourable; because both act on one shared candidate pool
their estimates are correlated, and the paired test is the correct one. We note
this explicitly because it is the kind of comparison that is easy to get wrong.

What the score demonstrably does instead is select **rarer, more specific**
predictions. At a matched 10,000-prediction budget it takes 580 associations at a
0.05% base rate and a 0.98% hit rate, against significance's 169 associations at a
0.14% base rate and a 1.83% hit rate `[E13]`: a lower raw hit rate with higher
enrichment, i.e. harder and more informative predictions. This echoes the §3.3
finding that dcGO's advantage lives in high-information terms. Its other
contributions — redundant-signature removal and novelty discounting against
curated knowledge — are interpretability work, not ranking work.

**A structural limit.** The sharp end of the ranking cannot be validated this
way. The top 25 associations by surprise yield 117 predictions and 0 hits, but the
base rate predicts fewer than 0.2 hits there `[E14]`, so the cell is
uninformative rather than negative. The reason is intrinsic: emergence requires
that a combination's carriers are already nearly all annotated with the term,
which by construction leaves almost nothing outstanding to predict `[E15]`. Any
revised score should trade emergence off against how many standing predictions it
leaves.

### 3.7 Breadth across vocabularies

> **Recomputed 2026-08-04.** The values first committed for this subsection were
> produced before the species defect described in §2 was fixed: the
> UniProt-native annotation sources parsed all of Swiss-Prot rather than the
> selected species, and the contingency-table row space was the union of the
> annotation and domain protein sets, so the Fisher protein universe for every
> UniProt-sourced vocabulary was inflated from ~20k human proteins toward the
> 575,503 entries of the flat file. Restricting the annotation map to the
> domain-annotated intersection corrects it. Every row below has been
> regenerated; the superseded values are retained in the last column, and the
> Gene Ontology anchor — drawn from the species-specific annotation file and
> therefore never affected — serves as the control showing that the change is
> the correction rather than a change of protocol.

Applying the held-out enrichment statistic of §3.6 to seven vocabularies trained
on an archived April 2021 Swiss-Prot release gives the values in Table 3.

> **Era note.** The non-GO rows propagate over their own hierarchies, which
> are read by the light `is_a`/`part_of`-only OBO reader and its companions
> and were therefore never exposed to the regulates-edge defect. The **GO
> anchor row is pre-regulates-fix era** `[K3]` and is PROVISIONAL with the
> rest of the GO-propagated results; since it is the control the correction
> narrative leans on, the whole table should be re-anchored after the post-fix
> regeneration.

**Table 3.** Held-out enrichment by vocabulary, 2021 → 2026, after the
species correction. Values from `validation/temporal_breadth_metrics.tsv`.

| Vocabulary | assoc. | predictions | hits | hit rate | expected | enrichment [95% CI] | superseded |
|---|--:|--:|--:|--:|--:|---|--:|
| GO (anchor) | 91,830 | 3,106,234 | 106,224 | 3.42% | 0.30% | 11.5 [11.1, 12.0] | 11.3 |
| Reactome | 34,397 | 760,073 | 1,430 | 0.19% | 0.02% | 11.4 [8.9, 14.6] | 8.0 |
| Subcellular location | 6,340 | 172,679 | 3,798 | 2.20% | 0.59% | 3.7 [3.5, 3.9] | 2.9 |
| Cofactor (ChEBI) | 291 | 1,404 | 235 | 16.74% | 4.72% | 3.5 [2.1, 4.5] | 3.2 |
| UniProt keyword | 35,981 | 889,117 | 13,140 | 1.48% | 0.44% | 3.4 [3.2, 3.6] | 1.7 |
| ComplexPortal | 958 | 9,904 | 2 | 0.02% | ≈ 0% | 61 [0.0, 131] | 265 |
| OMIM phenotype | 27 | 90 | 0 | 0% | ≈ 0% | *undefined* | undefined |
| ChEBI ligand | — | — | — | — | — | *not testable* | — |

The correction moved the results in a direction worth stating explicitly: the
contamination was **suppressing** the cross-vocabulary signal, not manufacturing
it. Every affected layer enriches more strongly once the non-human proteins
leave the universe, and the UniProt keyword layer doubles (1.7 → 3.4) with
non-overlapping intervals. The GO anchor does not move (11.3 → 11.5), which is
what a correct correction should look like.

Two conclusions in the superseded analysis do not survive. First, Reactome at
11.4 [8.9, 14.6] now overlaps GO at 11.5 [11.1, 12.0], so the ordering "GO
strongest, Reactome second" is not supported; because Reactome rests on 1,430
hits against GO's 106,224 its interval is far wider, and the defensible
statement is that the two are indistinguishable at this power rather than that
Reactome matches GO. Second, ComplexPortal falls from 24 hits to 2 and its
interval now includes zero, so the layer shows no demonstrated signal at all.

Three rows must not be read at face value. **Complex** never had an
interpretable magnitude: ComplexPortal averages ~1.5 proteins per complex, so
the base rate rounds to zero and any hit at all divides by nearly nothing. The
superseded run's 24 hits gave a nominal enrichment in the hundreds with an
interval excluding one; after correction 2 hits give an interval including zero,
so the earlier reading of "signal exists, magnitude meaningless" was itself too
generous. **Disease** returned 0 hits on 369 predictions with an expected count
also near zero — an uninformative cell, not a negative result, and a direct
consequence of 6,904 OMIM phenotypes spread over 5,029 proteins, which yielded
only 53 significant associations at t0. This layer is unchanged by the species
correction, because OMIM phenotype cross-references are essentially human-only
already. Re-keying the same curation onto the Disease Ontology has since been
implemented and does not repair it: the residual sparsity is at the protein
level, roughly 0.8 proteins per Disease Ontology term, which no amount of term
pooling can fix.
**Ligand** cannot be tested at this split at all: UniProt's structured
`FT /ligand_id="ChEBI:…"` qualifier postdates the April 2021 release, so the layer
consists entirely of post-2022 annotation `[F10]`.

Two protocol caveats apply to the whole table regardless of the defect. First, no
evidence filter exists outside GO, so an automated annotation added between
snapshots counts as a hit `[F13]`; the GO anchor bounds how much this matters but
does so imperfectly, since the anchor and the strict §3.6 figure differ in
candidate set as well as in protocol `[F14]`. Second, the rows pool single
domains and supra-domains, so this table does **not** test whether the emergent
(combination-specific) claim generalises beyond GO `[F12]` — which is the more
interesting question and remains open.

### 3.8 Relative inference needs an information floor, and gets one

The in-inference relative test of §2.5 was built to report each domain at the
most specific term it supports. Run without further measures it does the
opposite: parentless terms skip the test and pass on the overall inference
alone, so the three GO aspect roots — carried by most of the universe, at
0.09–0.17 bits of annotation-frequency information `[K6]` — and the
near-universal band just above them ("biological regulation", "cytoplasm",
"organelle") dominate the output.

The `--min-ic` reporting floor removes exactly that failure mode. On the human
GO single-domain paper-parity run, sweeping the floor gives:

**Table 4.** The `--min-ic` sweep, human GO single domains, paper-parity
configuration, FDR < 0.01. "On a chain" is the share of significant
associations whose term is an ancestor of another significant term for the
same domain — the specificity failure the relative test exists to prevent.
Values from the sweep in `VALIDATION_PLAN.md` (item 2), computed with
`validation/specificity_metrics.py`; the sweep itself has no committed
machine-readable artifact, so the ledger carries it as provisional
`[K4, K5]` — but its floor-1 significant count is now independently
confirmed by the manifest-carrying production run (30,302 `[M2]`).

| `--min-ic` (bits) | significant | mean #ancestors | on an ancestor chain | GO roots present |
|---:|---:|---:|---:|---|
| 0 (off) | 30,655 | 6.0 | 52.7% | all three |
| 1 | 30,302 | 6.1 | 50.2% | none |
| 2 | 28,348 | 6.3 | 46.4% | none |
| 3 | 26,401 | 6.4 | 42.8% | none |
| 5 | 18,888 | 7.3 | 33.3% | none |

`--min-ic 1` — keep terms carried by under half the universe — removes the
roots and the near-universal band for 353 associations, 1.2% of the set
`[K6]`. That is the recommended setting: the floor's job is killing vacuous
terms, not making the method conservative.

What the floor does **not** fix should be stated as plainly. The residual
ancestor-chain cascade barely moves — 52.7% of associations sit on a chain
without the floor, 50.2% with it — so half the surviving output still restates
a more specific association at a more general level. A wider universe does not
help (an all-species run is *worse*, at 82.4% on a chain), and degenerate
parents and the attainable-*p* floor have been ruled out as causes `[K5]`. The
open idea is `elim`-style decorrelation in the manner of Alexa et al.: test
specific terms first, remove their proteins, and retest ancestors on the
residue. Until something of that kind is implemented and measured, relative
inference should be described as *structurally faithful to the original
method but not yet achieving its stated purpose* — which is how we describe
it here.

### 3.9 The multi-vocabulary expansion

The registry now holds thirty-five `--ontology` keys `[M8]`, and the round of
work reported here added fourteen annotation layers in two waves. The first
wave is the gene-keyed phenotype infrastructure of §2.7: the Human Phenotype
Ontology, SynGO, and five model-organism phenotype vocabularies learned on
their own organisms' proteomes. The second (merged after the first draft of
this section) extends the same machinery `[L9, L10]`: Mondo and
Orphanet–Mondo term re-keys of the UniProt disease layer; GWAS-Catalog EFO
traits (genetic-association evidence — the loosest evidence type in the
registry, and labelled as such); HPA single-cell **expression** layers
(`celltype`, flat, on HPA's own cell-type names after Cell Ontology name
matching fell below the pre-set usefulness floor); WormBase anatomy
**expression** (`wbbt` — "expressed in", unlike the phenotype layers'
"abnormal when mutated"); and the CIViC-chained cancer layers (`ncit`,
`oncotree`). A systematic acquisition pass preceded the adapters: **39 open
annotation and ontology sources were fetched, integrity-checked and
versioned** against UniProt release 2026_02 `[L8]`. The boundary of the open
data is itself a finding: SNOMED CT and MedDRA are licence-gated and were not
attempted; OMIM's `genemap2.txt` is registration-gated, with UniProt's `MIM`
cross-references covering the open part; and MAxO publishes an ontology but
no released annotation source, so there is nothing to test `[L8]`.

**The production matrix.** Every runnable registry ontology has now been run
under manifests: 63 cells, zero failures — one *baseline* cell per ontology,
one *paper-parity* cell wherever a hierarchy exists, plus all-species and
experimental-evidence GO variants (§3.10) — each writing
`results/production/<cell>/run_manifest_<ontology>.json` with the
`is_a`/`part_of` edge policy recorded, i.e. post-fix era throughout `[M0]`.
Table 5 reads the expansion layers from those artifacts.

**Table 5.** Significant single-domain associations per expansion layer, from
the production matrix `[M7]`. *Baseline* is the default configuration (no
propagation, no relative inference, no floors); *paper-parity* is
`--propagate-annotations --enable-relative-inference --enable-true-path
--min-ic 1` (§2.5). Both columns count rows with `domain_type = single` in
the cell's `domain_<ontology>_associations_significant.tsv`; single domains
and supra-domains are corrected as separate BH families (§2.3), so each count
is FDR < 0.01 within the single-domain family. `celltype` is a flat
vocabulary and has no paper-parity cell.

| Layer | Organism | Vocabulary (semantics) | Baseline | Paper-parity |
|---|---|---|---:|---:|
| `hpo` | human | Human Phenotype Ontology | 996 | 38 |
| `syngo` | human | SynGO | 484 | 241 |
| `mondo` | human | Mondo (re-keyed OMIM) | 9 | 112 |
| `orphanet_mondo` | human | Mondo (re-keyed Orphanet) | 205 | 116 |
| `efo` | human | EFO (GWAS traits) | 1,504 | 814 |
| `celltype` | human | HPA cell-type expression | 1,236 | — |
| `ncit` | human | NCI Thesaurus (via CIViC) | 9 | 10 |
| `oncotree` | human | OncoTree (via CIViC) | 5 | 0 |
| `mp` | mouse | Mammalian Phenotype | 873 | 261 |
| `wbphenotype` | worm | WBPhenotype | 26,624 | 13,809 |
| `wbbt` | worm | WormBase anatomy expression | 67,379 | 50,259 |
| `zfa` | zebrafish | Zebrafish anatomy (affected) | 37,776 | 31,994 |
| `fbcv` | fly | FlyBase phenotype class | 5,033 | 2,704 |
| `fbbt` | fly | Drosophila anatomy (manifesting) | 10,791 | 8,316 |

Three readings, and their limits. First, the spread follows curation shape,
not biology: the expression and phenotype screens that annotate tens of
thousands of genes (worm, zebrafish, fly) support the largest sets, while
HPO's disease-derived universe (5,180 analysable proteins `[M7]`), SynGO's
deliberately narrow expert curation (1,781), and the two-hop cancer layers
(435 and 300) support correspondingly few. Second, the paper-parity column
moves in **both directions**: sparse re-keyed layers grow (input propagation
pools annotations onto testable terms before the statistics — `mondo` 9 →
112, `cofactor` and `tcdb` behave likewise `[M7]`), while richly annotated
layers shrink as the intersection–union test prunes parent-driven
associations. The extreme is `hpo`, which collapses 996 → 38; we flag that
collapse as **unexplained** and requiring investigation before the layer's
paper-parity output is used for anything. Third, and most important: **none
of these counts is validated**. No held-out temporal test, no permutation
control and no reference comparison has yet been run for any of these
layers; the counts are the size of hypothesis sets, exactly as §3.1 cautions
for GO. The model-organism layers carry the original dcGO's transfer
argument — the association is learned on the model organism's proteins and
the domain is species-agnostic — but that argument has not yet been tested
here either (for instance, against human disease annotation via the HPO
layer). We record the expansion as infrastructure with verified identifier
mapping (§2.7), honest per-layer semantics, and manifest-carrying first
runs — not as validated predictive breadth.

### 3.10 Training-universe and evidence-mode variants

The production matrix also regenerated, under manifests and post-fix, the two
GO variants this project has studied outside the primary human/`manual`
configuration.

**Training universe.** The all-species background (`--species allspecies`;
1,464,355 proteins, 28,112 testable GO terms) yields 2,911,662 significant
associations at baseline (535,133 single-domain) `[M3]` and 5,597,840 under
paper-parity (1,051,061 single-domain) `[M4]`. Two things are worth stating
plainly. First, the paper-parity configuration moves in the **opposite
direction** from human GO — it roughly doubles the all-species set where it
cut the human set by 42% — which is consistent with input propagation adding
testable annotation mass faster than the relative test prunes it at this
scale, though we have not isolated that mechanism. Second, the *evaluative*
claims for the wider universe — 8/9 F_max and 9/9 AUPRC cells won on the
held-out split, 9/9 and 9/9 under an experimental-evidence filter — still
come from the pre-fix evaluation `[K10]` and remain provisional; what the
matrix supplies is the manifest-carrying training runs a post-fix
re-evaluation can now be built on. The known caveats travel with the counts:
the all-species `manual` universe is majority-projected annotation and
support is inflated ~2.44× by orthology `[K10]`, so the raw significant
counts above should not be read as 17× the human evidence.

**Evidence mode.** Restricting training to experimentally supported
annotations (`--evidence-filter experimental`; 16,242 proteins, 12,966 terms)
yields 62,426 significant associations at baseline (18,310 single-domain)
`[M5]` and 28,250 under paper-parity (10,127 single-domain) `[M6]` — about
38% of the `manual` baseline's count, the price of dropping curated-inference
evidence from the training signal. No evaluation of the experimental-mode
associations has been run at these settings; the cells exist so that the
evidence-policy sensitivity of any future result can be quoted from
manifest-carrying runs rather than re-derived.

---

## 4. Discussion

The defensible summary of this work is narrow and we state it as such:
**domain-derived associations contain predictive signal for later human GO
annotations and outperform simple baselines in several retrospective settings,
particularly for higher-information GO terms.**

Three strands support it. Associations recover about two-thirds of the curated
InterPro2GO reference on shared domains `[B1]`, which shows the inference is
finding the same kind of thing curators do. Predictions transferred to
no-knowledge proteins beat a frequency baseline on F_max in two of three aspects
unfiltered and in all three once uninformative terms are excluded, while staying
1.3–25× above a random-domain null `[C2]`. And the supra-domain associations
found in 2021 anticipate 2026 curation at 12.5× the terms' own acquisition rates
`[E2]` — the cleanest out-of-sample evidence in the paper, and the one that
speaks directly to the domain-combination claim that distinguishes this method
from homology transfer.

Equally, three strands limit it, and we regard reporting them as part of the
contribution.

**The information-content filter is doing visible work.** Unfiltered, dcGO loses
molecular function to a frequency baseline `[C2e]`. The filter is defensible on
its face — near-universal terms carry no information and rewarding their recovery
measures curation demographics rather than prediction — but it was chosen after
seeing the data, and it changes the evaluated cohort as well as the scoring
universe `[C1e]`. A pre-specified endpoint, or confirmation on an untouched
interval, would settle this; neither exists yet.

**The metric choice matters more than it should.** F_max and AUPRC disagree,
with the naive baseline winning AUPRC in all three aspects unfiltered `[C3i]`. A
method whose apparent superiority depends on which of two standard metrics is
reported has not established superiority. Our reading is that the disagreement is
partly structural — an all-terms predictor integrates favourably over recall —
but this is untested, and our AUPRC implementation has not been checked against
reference tooling `[H18]`.

**Ranking emergent predictions remains unsolved.** The surprise score was built
to make the combination-specific predictions readable, and it does: it removes
redundant InterPro signatures, discounts what curators already record, and
biases toward rare, specific terms `[E13]`. It does not rank better than a plain
q-value `[E10]`. There is a structural reason to expect this to be hard:
emergence, as defined, requires that the combination's carriers are already
almost all annotated, which leaves nothing to predict `[E15]`. A score that
explicitly trades emergence against outstanding predictions is the obvious next
step.

**Relation to the original dcGO.** Our reading of the original papers —
recorded in project documentation rather than re-verified here, and therefore
marked provisional — is that it operated over all sequenced genomes with
SCOP/SUPERFAMILY and Pfam domains at FDR < 10⁻³, combined an overall with a
relative parental-background test, and transferred predictions by a per-target
normalised sum `[G1–G5]`. Three of those pieces are now present here (the
p-score transfer; the relative test, folded into the inference before FDR
correction as published; input-map propagation), one differs (InterPro entries
rather than SCOP/Pfam signatures by default, a coarser but not disjoint
universe, since InterPro integrates both `[G7]`), and one is looser
(FDR < 0.01 against their 10⁻³, bracketed in the comparison below).

The quantitative comparison the first draft called the single most valuable
missing analysis has since been run: re-keying the parser on the `SSF`
member-signature (`--domain-key ssf`) puts our associations in the published
release's SCOP-superfamily space, where **precision — the fraction of our
calls the published dcGO also made — sits in a narrow 0.537–0.625 band across
six pre-declared threshold × definition variants** `[K9]`. Recall against them
is not interpretable: they are all-species and 2016 where we are human-only
and 2026, and 69.4% of their pairs have zero co-occurring human proteins, so
a human universe cannot reach them at any threshold `[K9]`. Their SCOP-family
half and their Pfam-keyed release remain unreachable pending a `pfam` domain
key, and their non-GO tables (EC, keywords, UniPathway) have not been used.

---

## 5. Limitations

We list these as blocking or qualifying items rather than as a formality. All are
carried from an internal engineering and scientific review of this repository
`[H]`, whose own conclusion was that the work "does not yet establish robust
general performance, calibration, or superiority" `[H20]`.

**Evaluation design.**

1. **No untouched evaluation set.** Model-selection choices — the transfer rule,
   the IC thresholds, shrinkage settings — were compared on the same 2021 → 2026
   split that produced the headline `[H1]`. Choices should be frozen on a
   development interval or a nested temporal split and evaluated once on an
   untouched interval.
2. **The primary endpoint was not pre-specified**, and the unfiltered result
   loses to naive for molecular function `[H6]`. The IC-filtered analysis is a
   post hoc secondary analysis.
3. **IC thresholds alter the evaluation cohort** (MF −59% at IC ≥ 2), so results
   across floors are not paired and cohort sizes must always be reported
   alongside `[H7, C1e]`.
4. **No uncertainty estimates.** There are no protein-level bootstrap confidence
   intervals for F_max, AUPRC, or the paired differences against baselines
   `[H2]`; every entry in Tables 1 and 2 is a point estimate on cohorts of
   145–572 proteins.
5. **The random-domain null is a single seeded shuffle**, not a permutation
   distribution `[H3]`. No empirical p-value or null interval exists, so
   "1.3–25× above random" is a ratio of two point estimates.
6. **The component ablation has now been run, and it is negative for two of
   three components.** Over 12 aspect × IC cells with a protein-level paired
   bootstrap: supra-domains improve 0/12 cells; the shrinkage rung moved
   0/12 cells (and the step has since been removed `[K7]`); and the True Path
   rung is significantly worse in 12/12 — but that last figure is a statement
   about the then-combined filter-plus-propagation, measured with a
   background defect since fixed and over the pre-fix edge set, so it cannot
   be attributed to propagation and must be re-measured against the split
   flags `[K8]`. On this benchmark the best configuration is single domains
   only; the supra-domain machinery's demonstrated value is the emergent
   combinations of §3.5–3.6, not protein-centric F_max.
7. **Weak comparators.** A frequency baseline and a shuffled null establish that
   signal exists; they do not establish utility. No comparison exists against
   original dcGO output, against a homology-transfer baseline, or against any
   independent function predictor `[H10]`.
8. **One species, one interval — partially addressed.** The reported results
   remain human-trained over a single 2021 → 2026 interval `[H11]`. An
   all-species training universe (1,464,355 proteins across 9,074 taxa) has
   since been run and wins 8/9 F_max and 9/9 AUPRC cells on the same held-out
   split — 9/9 and 9/9 under an experimental-evidence filter — but its
   evaluation is pre-regulates-fix era and its support counts are inflated
   ~2.44× by orthology, so we cite it as direction, not as a result of this
   paper `[K10]`. Manifest-carrying post-fix training runs for both the
   all-species and experimental-evidence variants now exist (§3.10
   `[M3–M6]`); the missing piece is the post-fix re-evaluation on top of
   them.
9. **Prediction coverage is not reported** alongside F_max, although CAFA
   precision omits proteins with no predictions while recall includes them
   `[H17]`; and the F_max/AUPRC implementation has not been verified against
   independent CAFA tooling `[H18]`.

**Statistical validity.**

10. **The shrinkage step was a heuristic, not a fitted empirical-Bayes model,
    and has been removed.** It geometrically interpolated observed and
    constituent *p*-values with a hand-set decay; the transformed quantities
    were not valid *p*-values, so BH over them did not control FDR `[H5]`, and
    enabling it inflated rejections by 184% with no measurable effect on
    prediction quality `[K7]`. It was disabled in every run reported here and
    no longer exists in the codebase.
11. **BH is applied across a strongly dependent hypothesis family** — nested GO
    terms, supra-domains containing their constituents, co-occurring domains —
    without simulation or hierarchical multiple-testing correction `[H15]`.
12. **No minimum-support or effect-size policy.** FDR significance alone retains
    associations from very sparse tables; contingency cells and odds-ratio
    confidence intervals are not reported, and odds ratios are emitted without a
    Haldane correction so they degenerate to 0 or ∞ at empty cells `[H19, H21]`.
13. **InterPro2GO recovery is coverage, not precision, and is not independent
    validation** `[H16]`.

**Temporal validity.**

14. **Look-ahead through domain architectures.** Domain assignments come from the
    **current** `protein2ipr`, not a 2021 reconstruction `[H9]`. This is an
    annotation-temporal benchmark, not a fully prospective simulation, and we
    describe it as such throughout.
15. **The t1 annotation snapshot and the GO release used for propagation are
    mutable, unpinned files** `[H8]`; only the t0 release (GOA 205) is pinned.
16. **The domain-centric evaluation scores against the current InterPro2GO**, so
    it is not temporal on the domain side `[B10]`.

**Scope and implementation.**

17. **Every GO-propagated result in this draft is pre-regulates-fix era.** The
    species-filtering defect that previously invalidated §3.7 has been fixed
    and those rows regenerated `[F]`; the outstanding regeneration debt is now
    the propagation-edge fix: §3.2–3.4, §3.6 and the §3.7 GO anchor must be
    reproduced from post-#67, manifest-carrying runs `[K3, K12]`.
18. **The A–D configuration-comparison artifacts are now committed**
    (`validation/bench_A`–`bench_D`), closing the provenance gap the review
    raised `[H12]`; the comparison itself still lacks a paired test or
    confidence interval and stays provisional (§3.3).
19. **Run manifests cover the association counts but none of the
    evaluations.** Every run writes `run_manifest_<ontology>.json` with
    input/output SHA-256s, release headers, Git state, dependency-lock hash,
    command line and thresholds `[K12]`, and the production matrix regenerated
    every association count under one `[M0]` — but the evaluation results
    (§3.2–3.7) all rest on pre-manifest runs and are not exactly reproducible
    today. Regenerating the evaluations from the manifest-carrying runs is a
    prerequisite for submission.
20. **Implementation caveats that touch reported behaviour:** a numerical failure
    in the hypergeometric score falls back to a value of 50.0, i.e. a plausible
    medium-confidence score `[H13]`; and `--num-cores` is logged but does not
    drive the Fisher implementation, so no parallel-scaling claim is made `[H14]`.

Four correctness defects were found and fixed during development — an int8
overflow in contingency counts, non-binary incidence matrices, occurrence-versus-
protein counting, and an inverted GO DAG traversal — each with a regression test.
An earlier reported coverage figure was corrupted by the first three. We mention
this because it is the reason the project's conventions now require sanity-checking
extreme rows before quoting any number, and because it bears on how much
confidence to place in numbers that have not yet been independently regenerated.

---

## 6. Data and code availability

**Code.** The implementation, evaluation harnesses and this manuscript live in
the `dcGO-2.0` repository, MIT-licensed. The supported entry point is
`run_dcgo_human.py`; evaluation code is under `validation/`. Committed metrics
files are `validation/performance_metrics.tsv`,
`validation/domain_centric_metrics.tsv`,
`validation/temporal_benchmark_metrics.tsv`,
`validation/temporal_surprise_metrics.tsv`,
`validation/temporal_surprise_associations.tsv`,
`validation/temporal_breadth_metrics.tsv` and
`validation/temporal_breadth_go.tsv`. *[A public URL and archived DOI are
required before submission and do not yet exist.]*

**Data.** All inputs are public: InterPro `protein2ipr` [7]; UniProt-GOA
annotation files, including the dated archive from which release 205 (April 2021)
was taken; the Gene Ontology `go-basic.obo` [4]; the UniProt Swiss-Prot flat file
and its companion vocabularies [8]; Expasy ENZYME; InterPro2GO; and, for the
layers of §2.7, the HPO `genes_to_phenotype` release with `hp.obo`, the SynGO
bulk release, MGI's `MGI_GenePheno`/`MRK_SwissProt_TrEMBL` reports with
`mp.obo`, WormBase's phenotype-association file with `wbphenotype.obo`, ZFIN's
`phenoGeneCleanData_fish`/`uniprot.txt` with `zfa.obo`, FlyBase's
genotype–phenotype and identifier-mapping tables with `fbcv.obo`/`fbbt.obo`,
and UniProt's per-organism idmapping files. The acquisition record — source
URLs, releases, integrity checks and licence status for all 39 fetched
sources — is `data/ACQUISITION_MATRIX.md` `[L8]`.

**Reproducibility status — association counts regenerated, evaluations not
yet.** Every run emits a machine-readable manifest
(`run_manifest_<ontology>.json`: input/output SHA-256 digests, release headers
where the format supplies one, Git state, dependency-lock hash, command line,
thresholds, and the propagation-edge policy as an era marker), with a
completion checklist in `REPRODUCIBILITY.md` `[K12]`, and the 2026-08-18
production matrix regenerated every association count in this manuscript
under one (63 cells, driver `scripts/run_production_matrix.py`) `[M0]`.
However, only the t0 annotation release is pinned by identifier; the t1
snapshot, the GO release, the InterPro release and the Swiss-Prot releases
are still referenced through mutable "current" URLs `[H8]`; and **every
evaluation result (§3.2–3.7) predates the manifest machinery**, so those
numbers are not exactly reproducible from a clean checkout. Closing this —
the post-fix, manifest-based regeneration of every evaluation table, pinned
releases, an archived input snapshot, and one-command regeneration — is a
prerequisite for submission.

---

## References

Verification status for each reference is recorded in
`paper/EVIDENCE_LEDGER.md` §J.

1. Fang H, Gough J. A domain-centric solution to functional genomics via dcGO
   Predictor. *BMC Bioinformatics* 2013;**14**(Suppl 3):S9.
   doi:10.1186/1471-2105-14-S3-S9. *(Verified from the PDF in `docs/`.)*
2. Fang H, Gough J. dcGO: database of domain-centric ontologies on functions,
   phenotypes, diseases and more. *Nucleic Acids Res* 2013;**41**(D1):D536–D544.
   *(Verified from the PDF in `docs/` and by search.)*
3. The dcGO Domain-Centric Ontology Database in 2023: New Website and Extended
   Annotations for Protein Structural Domains. *J Mol Biol*
   2023;**435**(14):168093. doi:10.1016/j.jmb.2023.168093. *(Title, venue and DOI
   verified from the PDF in `docs/`; **author list not independently verified** —
   recorded in project notes as Bao et al. and to be confirmed before
   submission.)*
4. The Gene Ontology Consortium. Expansion of the Gene Ontology knowledgebase and
   resources. *Nucleic Acids Res* 2017;**45**(D1):D331–D338.
   doi:10.1093/nar/gkw1108. *(Verified from the PDF in `docs/` and by search. A
   more recent GO Consortium reference should be substituted before submission.)*
5. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
   powerful approach to multiple testing. *J R Stat Soc Series B* 1995;**57**(1):
   289–300. *(Verified by search.)*
6. Zhou N, Jiang Y, Bergquist T, et al. The CAFA challenge reports improved
   protein function prediction and new functional annotations for hundreds of
   genes through experimental screens. *Genome Biol* 2019;**20**:244.
   doi:10.1186/s13059-019-1835-8. *(Verified by search.)*
7. Paysan-Lafosse T, Blum M, Chuguransky S, et al. InterPro in 2022. *Nucleic
   Acids Res* 2023;**51**(D1):D418–D427. doi:10.1093/nar/gkac993. *(Verified by
   search. **The InterPro release actually used by these runs is not recorded in
   the repository**, so this citation identifies the resource, not the version.)*
8. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025.
   *Nucleic Acids Res* 2025;**53**(D1):D609–D617. doi:10.1093/nar/gkae1010.
   *(Verified by search. The runs used a 2026-07 release and an archived 2021_02
   release; neither has a citable release identifier recorded in the
   repository.)*

*No reference is cited that could not be verified from a PDF in `docs/` or by
search. In particular, no citation is made to the reported domain counts of the
2023 dcGO release, because the project's own notes require that figure to be
checked against the actual download first `[G9]`.*

---

## Open questions for the authors

These are the points where the repository evidence was ambiguous, contradictory,
or absent, and where an author decision is needed before this draft can be
completed.

**Discrepancies between committed data and prose (must be resolved).**

1. **The surprise-score paired bootstrap disagrees with its own write-up.**
   `validation/temporal_surprise_metrics.tsv` records intervals
   [−87.28, +18.95] / [−8.93, +7.06] / [−1.01, +4.54] and
   favouring-fractions 0.700 / 0.320 / 0.800, while `SURPRISE_SCORE.md` states
   [−85.26, +20.87] / [−9.39, +10.23] / [−2.15, +5.25] and 75% / 54% / 82%
   `[E9, E11]`. Which run is authoritative? The verdict (no demonstrated
   advantage) is unchanged either way, but the numbers must agree.
2. **Related, and more concerning:** in the committed file, at the
   10,000-prediction budget the point estimate (+7.96) lies **outside** its own
   percentile interval, and only 32% of resamples favour surprise despite a
   positive point difference `[E12]`. Is this heavy-tailed instability, a seeding
   inconsistency between the point estimate and the resamples, or a bug in
   `_budget_enrichment` when a resample contains duplicate associations?
3. **Whole-pool and rank-slice confidence intervals also disagree** between the
   committed metrics file and `SURPRISE_SCORE.md` (e.g. 12.5× CI [10.61, 14.13]
   vs [10.9, 14.4]) `[E3, E6]`, and the "expected hits" for the top-25 slice is
   ≈0.09 by the file and 0.14 in the prose `[E14]`. Same question.
4. **Two different test-suite sizes are documented** (162 tests in the review,
   155 in `CLAUDE.md`/`README.md`) `[A13]`. Trivial, but it should be a generated
   number.

**Missing artifacts behind reported numbers.**

5. **The A–D configuration comparison** (base / +p-score / +relative / both) is
   quoted only in prose; `validation/bench_A`–`bench_D` are untracked `[C5a–C5e,
   H12]`. Can these be regenerated and committed? Without them, the claim that
   the p-score transfer is "the main lever" cannot be reported.
   *Partially resolved (2026-08-17): the bench_A–D metrics files are now
   committed; the paired test and confidence intervals remain missing, and the
   runs are pre-regulates-fix era.*
6. **`protein binding` = 84.6% of human experimental MF annotations** appears in
   three documents but is computed nowhere `[C6a]`. Likewise the `hyper_score`
   saturation figure (~37% at exactly 100) `[C6c]`, the pre-filter supra-domain
   counts (123,203 GO / 8,637 EC) `[D5, D9]`, and the pre-propagation
   InterPro2GO size (30,190 pairs) and shared-domain count (2,747) `[B6]`. Should
   these be emitted by the code, or dropped from the paper?
7. **The count of human proteins carrying InterPro domains (18,908)** `[A9]` is
   documented but not present in any located run log. Which run produced it?

**Design and interpretation decisions.**

8. **What is the primary endpoint?** The review requires one to be pre-specified
   `[H6]`. Our recommendation is: F_max at IC ≥ 2, per aspect, on an untouched
   interval, with IC ≥ 0 reported as a secondary analysis — but this is the
   authors' call and it determines how §3.3 is framed.
9. **Why does the naive baseline win AUPRC?** `[C3i]` We hypothesise that an
   all-terms predictor integrates favourably over the recall range while dcGO's
   curve truncates, but we have not tested this, and the AUPRC implementation's
   quantile sweep and upper-envelope trapezoidal rule are unverified `[C0g]`.
   This needs either an explanation backed by prediction-coverage numbers or a
   verified reimplementation.
10. ~~**Is the relative test going to be folded into inference, or stay post
    hoc?**~~ *Resolved (2026-08-17): folded into the inference before BH as
    `max(overall_p, relative_p)`, per the original method (§2.5) `[K11]`. The
    threshold divergence (0.01 vs their 10⁻³) remains and is bracketed in the
    published-dcGO comparison `[K9]`.*
11. **True Path Rule default.** The original made propagation central; here it is
    opt-in and was off in every reported run. Should the paper report a
    propagated configuration as primary? *Updated context: the paper-parity
    configuration (§2.5) enables it, but the ablation's True Path rung
    (confounded, pre-fix) is the current evidence against making it the
    protein-centric primary `[K8]`; the re-measurement must come first.*
12. **Species scope.** Human-only is a real limitation `[H11]`. Is a second
    species in scope for this paper, or is it explicitly deferred?
    *Updated context: an all-species background run `[K10]` and the
    model-organism phenotype layers (§3.9) both exist now; what remains open
    is which of them this paper reports as results rather than direction.*

**The breadth section.**

13. ~~**The species-filtering defect must be fixed and every UniProt-native
    vocabulary retrained** before §3.7 can say anything `[F]`.~~ *Resolved
    (2026-08-04): the sources are species-restricted, the row space is the
    intersection, and every §3.7 row was regenerated — the corrected values
    and the superseded ones are both in the table and the ledger `[F]`.*
14. **Should the breadth test be rerun with `--supra-only`?** As it stands it
    pools single domains and supra-domains, so it cannot address whether the
    *emergent* claim generalises beyond GO `[F12]` — which is the more
    interesting question and the one a reviewer will ask.
15. **How should the degenerate rows be presented?** ComplexPortal's
    hundreds-fold ratio and OMIM's undefined ratio are both artefacts of
    near-zero base rates `[F8, F9]`. Our draft reports them as "detectable,
    magnitude meaningless" and "uninformative" respectively; an alternative is to
    omit them from the table and describe them in text.

**Comparison to prior work.**

16. **The dcGO papers in `docs/` could not be read in this environment** (no PDF
    text extraction available), so all statements about the original method's
    scope, thresholds and validation are second-hand from project notes `[G1–G6,
    G10]`. Every such statement needs a direct check against the papers before
    submission.
17. ~~**The §3 comparison against published dcGO output has not been run**
    `[G11]`.~~ *Resolved (2026-08-04): run via `--domain-key ssf`; precision
    0.537–0.625 across six pre-declared variants, recall not interpretable —
    see Discussion and `[K9]`. Still open within it: the `pfam` domain key and
    the published non-GO tables.*
18. **Author list for reference [3]** was not independently verified.

**Framing.**

19. **What is this paper for?** As it stands the strongest, cleanest result is
    the 12.5× held-out enrichment of supra-domain associations `[E2]` — an
    out-of-sample test of the domain-combination claim that homology transfer
    cannot make. The F_max comparisons are weaker, metric-dependent, and rest on
    a post hoc filter. One option is to restructure around the emergence result
    and demote the CAFA-style table to a supporting analysis. That would be a
    smaller paper, but a more defensible one.
