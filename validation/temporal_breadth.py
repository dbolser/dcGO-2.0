#!/usr/bin/env python3
"""Does dcGO predict future curation *beyond GO*? — VALIDATION_PLAN §2/§7.

`validation/temporal_surprise.py` showed that the GO associations dcGO found in
2021 anticipate 2026 curation at 12.5x the terms' own acquisition rates. That is
one ontology. This asks the same question of the whole breadth added in
`src/ontology_registry.py` — pathways, localisation, chemistry, disease,
complexes — using the identical statistic, so the numbers are comparable across
ontologies.

**Why this is cheap to do.** Every UniProt-native ontology is harvested from one
file, so a single archived Swiss-Prot release supplies t0 for all of them at
once. UniProt 2021_02 (07-Apr-2021) lines up with GOA release 205 (2021-04), so
the split matches the GO benchmark's and the results are directly comparable.

For one ontology:

* Train at t0: ``run_dcgo_human.py --ontology X --uniprot-dat <archived>``.
* **Predictions** — proteins carrying feature ``S`` that lacked term ``t`` at t0.
* **Hits** — those carrying ``t`` at t1.
* **Base rate** — the term's acquisition rate over the whole universe, which is
  what makes enrichment comparable between a 46-term vocabulary (cofactors) and
  a 6,904-term one (OMIM).

Two design points that matter for correctness:

* **The universe is restricted to proteins that existed at t0.** A protein with
  no Swiss-Prot entry in 2021 cannot have been "predicted" to gain a term; it
  simply had no entry. Counting new entries as successful predictions would
  inflate every ontology, and unevenly — the layers that grew fastest would look
  best. :func:`src.uniprot_annotation_source.parse_uniprot_accessions` supplies
  the t0 membership that an ``AnnotationSource`` cannot.
* **No evidence filter is available off GO.** GO has evidence codes, so the GO
  benchmark trains on non-IEA and scores against experimental only. The
  UniProt-native layers are curated cross-references with no such split, so t0
  and t1 are simply the two snapshots. That is a weaker guarantee against
  automated annotation being counted as truth, and it is stated rather than
  hidden.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Set

project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from validation.temporal_surprise import (  # noqa: E402
    AssociationOutcome,
    StratumResult,
    acquisition_base_rates,
    pool,
    propagate,
    score_association,
)

#: Ontologies worth testing: broad enough coverage to yield associations, and
#: spanning distinct kinds of claim (function, pathway, localisation,
#: chemistry, disease, assembly).
DEFAULT_ONTOLOGIES = [
    "reactome",
    "keyword",
    "subcellular",
    "ligand",
    "cofactor",
    "disease",
    "complex",
]


def load_associations(
    path: Path, term_column: str, supra_only: bool = False
) -> List[tuple]:
    """Read ``(feature, term, q_value)`` from a significant-associations TSV."""
    import csv

    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if supra_only and row.get("constituent_domains") == "-":
                continue
            rows.append((row["domain"], row[term_column], float(row["adj_p_value"])))
    return rows


def evaluate_ontology(
    ontology: str,
    associations: List[tuple],
    carriers_of: Dict[str, Set[str]],
    t0_map: Dict[str, Set[str]],
    t1_map: Dict[str, Set[str]],
    universe: Set[str],
    bootstrap: int = 500,
    seed: int = 0,
) -> StratumResult:
    """Pool one ontology's associations into a single enrichment figure."""
    base_rates = acquisition_base_rates(
        {term for _f, term, _q in associations}, t0_map, t1_map, universe
    )
    outcomes: List[AssociationOutcome] = []
    for feature, term, q_value in associations:
        outcome = score_association(
            feature=feature,
            term=term,
            carriers=carriers_of.get(feature, set()),
            t0_map=t0_map,
            t1_map=t1_map,
            base_rate=base_rates.get(term, 0.0),
            rank_scores={"dcgo": -q_value},
        )
        if outcome.n_predicted:
            outcomes.append(outcome)
    return pool(ontology, outcomes, bootstrap, seed)


def main() -> int:  # pragma: no cover - I/O wiring
    import argparse

    from loguru import logger

    from src.domain_annotation_parser import DomainAnnotationParser
    from src.ontology_registry import get_ontology
    from src.uniprot_annotation_source import parse_uniprot_accessions

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ontologies", nargs="+", default=DEFAULT_ONTOLOGIES)
    parser.add_argument(
        "--t0-uniprot",
        type=Path,
        default=Path("data/raw/uniprot_archive/uniprot_sprot.dat.gz"),
        help="Archived Swiss-Prot flat file (t0)",
    )
    parser.add_argument(
        "--t1-uniprot",
        type=Path,
        default=Path("data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz"),
    )
    parser.add_argument(
        "--t0-gaf",
        type=Path,
        default=Path("data/raw/goa_archive/goa_human.gaf.205.gz"),
        help="Archived GOA GAF, so 'go' can be included as the anchor row",
    )
    parser.add_argument(
        "--t1-gaf",
        type=Path,
        default=Path("data/raw/goa_annotations/goa_human.gaf.gz"),
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--t0-results",
        type=Path,
        default=Path("results_t0_2021"),
        help="Directory of t0 association TSVs (one per ontology)",
    )
    parser.add_argument(
        "--interpro", type=Path, default=Path("data/interim/protein2ipr_human.dat.gz")
    )
    parser.add_argument(
        "--subcell", type=Path, default=Path("data/raw/uniprot_subcell/subcell.txt")
    )
    parser.add_argument(
        "--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_lite.obo")
    )
    parser.add_argument(
        "--reactome-relations",
        type=Path,
        default=Path("data/raw/reactome_relations/ReactomePathwaysRelation.txt"),
    )
    parser.add_argument(
        "--keyword-list",
        type=Path,
        default=Path("data/raw/uniprot_keywlist/keywlist.txt"),
    )
    parser.add_argument(
        "--doid-obo",
        type=Path,
        default=Path("data/raw/disease_ontology/doid.obo"),
        help="Disease Ontology, for --ontologies doid|orphanet_doid",
    )
    parser.add_argument(
        "--enzyme-dat",
        type=Path,
        default=Path("data/raw/enzyme/enzyme.dat"),
        help="Expasy ENZYME, for --ontologies ec",
    )
    parser.add_argument(
        "--supra-only",
        action="store_true",
        help="Score only supra-domain associations (the emergent claim), "
        "instead of every significant association",
    )
    parser.add_argument(
        "--no-propagate",
        dest="propagate",
        action="store_false",
        default=True,
        help="Skip True Path propagation of the two snapshots",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("validation/temporal_breadth_metrics.tsv")
    )
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    required = [args.t0_uniprot, args.t1_uniprot, args.interpro]
    if "go" in args.ontologies:
        # GO reads GAFs rather than the UniProt flat file, so without this the
        # missing input only surfaced as a parse failure after the expensive
        # architecture pass.
        required += [args.t0_gaf, args.t1_gaf, args.go_ontology]
    for path in required:
        if not path.exists():
            logger.error(f"Missing required input: {path}")
            return 1

    logger.info("Parsing domain architectures...")
    architectures = DomainAnnotationParser(
        max_supra_domain_length=3, min_domain_length=10
    ).parse_protein2ipr_file(args.interpro)

    logger.info("Reading the t0 accession set (proteins that existed in 2021)...")
    t0_present = parse_uniprot_accessions(args.t0_uniprot)
    universe = set(architectures) & t0_present
    logger.info(
        f"  Universe: {len(universe):,} proteins with domains and a 2021 entry "
        f"(of {len(architectures):,} with domains)"
    )

    paths_t0 = {
        "gaf": args.t0_gaf,
        "go_obo": args.go_ontology,
        "uniprot_dat": args.t0_uniprot,
        "subcell": args.subcell,
        "chebi_obo": args.chebi_obo,
        "reactome_relations": args.reactome_relations,
        "keywlist": args.keyword_list,
        "doid_obo": args.doid_obo,
        "enzyme_dat": args.enzyme_dat,
    }
    paths_t1 = dict(paths_t0, uniprot_dat=args.t1_uniprot, gaf=args.t1_gaf)

    # This dict is a second copy of the one in run_dcgo_human.build_ontology_paths,
    # and a registry entry added after it was written will not be in it. That is
    # not hypothetical: adding `doid` to --ontologies crashed a two-hour run with
    # KeyError('doid_obo') at the first parse, after the expensive architecture
    # pass. Check every requested ontology up front instead.
    unresolvable = {
        ontology: sorted(
            (
                set(get_ontology(ontology).needs)
                | set(get_ontology(ontology).hierarchy_needs)
            )
            - set(paths_t0)
        )
        for ontology in args.ontologies
    }
    unresolvable = {k: v for k, v in unresolvable.items() if v}
    if unresolvable:
        for ontology, missing in unresolvable.items():
            logger.error(
                f"--ontologies {ontology} needs input(s) {missing}, which this "
                "script does not resolve. Add them to paths_t0 above."
            )
        return 1

    results: List[StratumResult] = []
    for ontology in args.ontologies:
        entry = get_ontology(ontology)
        assoc_path = args.t0_results / f"domain_{ontology}_associations_significant.tsv"
        if not assoc_path.exists():
            logger.warning(f"[{ontology}] no t0 associations at {assoc_path}; skipping")
            continue

        logger.info(f"[{ontology}] loading t0 and t1 annotations...")
        t0_raw = entry.build_source(paths_t0, {}).parse()
        t1_raw = entry.build_source(paths_t1, {}).parse()

        ancestors_factory = entry.build_ancestors
        if ontology == "go" and args.go_ontology.exists():
            # GO's registry entry defers to OntologyProcessor (which also does
            # optimal-level filtering); for propagation alone the light OBO
            # reader is equivalent and much cheaper.
            from src.hierarchy import closure_ancestors, parse_obo_child_parents

            ancestors_factory = lambda paths: closure_ancestors(  # noqa: E731
                parse_obo_child_parents(args.go_ontology)
            )

        if args.propagate and ancestors_factory is not None:
            ancestors = ancestors_factory(paths_t0)
            t0_map = {p: propagate(t, ancestors) for p, t in t0_raw.items()}
            t1_map = {p: propagate(t, ancestors) for p, t in t1_raw.items()}
        else:
            t0_map = {p: set(t) for p, t in t0_raw.items()}
            t1_map = {p: set(t) for p, t in t1_raw.items()}

        associations = load_associations(
            assoc_path, f"{ontology}_term", supra_only=args.supra_only
        )
        logger.info(f"[{ontology}] {len(associations):,} t0 associations")

        # A sparse ontology (cofactors, complexes) annotates only a slice of the
        # proteome, so scoring against every domain-carrying protein conflates
        # two things: whether the right *term* was predicted, and whether the
        # protein is the kind that gets annotated in this ontology at all. The
        # in-scope universe removes the second by keeping only proteins the
        # ontology reaches by t1 — a stricter, term-specific test.
        in_scope = {p for p in universe if t1_map.get(p)}
        logger.info(
            f"[{ontology}] in-scope universe: {len(in_scope):,} of {len(universe):,} "
            f"proteins are annotated in this ontology by t1"
        )
        universes = [("all-domains", universe)]
        if in_scope and len(in_scope) < len(universe):
            universes.append(("in-scope", in_scope))

        wanted = {feature for feature, _t, _q in associations}
        for label, eval_universe in universes:
            carriers_of: Dict[str, Set[str]] = {}
            for protein in eval_universe:
                arch = architectures[protein]
                for feature in list(arch.single_domains) + list(arch.supra_domains):
                    if feature in wanted:
                        carriers_of.setdefault(feature, set()).add(protein)

            result = evaluate_ontology(
                f"{ontology} ({label})",
                associations,
                carriers_of,
                t0_map,
                t1_map,
                eval_universe,
                args.bootstrap,
                args.seed,
            )
            results.append(result)
            logger.info(
                f"[{ontology}/{label}] {result.n_hit:,}/{result.n_predicted:,} hit "
                f"({100 * result.hit_rate:.2f}%), expected "
                f"{100 * result.expected_rate:.2f}% → enrichment "
                f"{result.enrichment:.2f} [{result.ci_low:.2f}, {result.ci_high:.2f}]"
            )

    if not results:
        logger.error("No ontology could be evaluated — are the t0 runs missing?")
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        f.write(
            "ontology\tn_associations\tn_predictions\tn_hits\thit_rate\t"
            "expected_rate\tenrichment\tci_low\tci_high\n"
        )
        for r in results:
            f.write(
                f"{r.name}\t{r.n_associations}\t{r.n_predicted}\t{r.n_hit}\t"
                f"{r.hit_rate:.5f}\t{r.expected_rate:.5f}\t{r.enrichment:.2f}\t"
                f"{r.ci_low:.2f}\t{r.ci_high:.2f}\n"
            )
    logger.info(f"✓ Wrote {len(results)} ontologies to {args.output}")

    logger.info("")
    logger.info(
        f"{'ontology':<14} {'assoc':>7} {'preds':>9} {'hits':>7} {'hit%':>7} "
        f"{'exp%':>7} {'enrich':>7} {'95% CI':>16}"
    )
    for r in results:
        logger.info(
            f"{r.name:<14} {r.n_associations:>7,} {r.n_predicted:>9,} {r.n_hit:>7,} "
            f"{100 * r.hit_rate:>6.2f}% {100 * r.expected_rate:>6.2f}% "
            f"{r.enrichment:>7.2f} {f'[{r.ci_low:.2f}, {r.ci_high:.2f}]':>16}"
        )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
