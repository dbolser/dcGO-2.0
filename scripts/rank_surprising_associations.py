#!/usr/bin/env python
"""Rank dcGO's emergent domain-combination predictions by surprise score.

Reads a ``domain_<ontology>_associations_significant.tsv`` produced by
``run_dcgo_human.py``, keeps the **supra-domain** rows, and scores each one for
how much the *combination* predicts that its constituent domains do not — see
:mod:`src.surprise_score` for the definition of the score and its three
components (emergence, distinctness, novelty).

The inputs of the original run are re-read (domain architectures and the
annotation source) because the emergence test needs the contingency counts for
the combination *and* for each of its parts, which the association table does
not carry.

Usage
-----
    # Gene Ontology, with the curated InterPro2GO novelty reference
    uv run python scripts/rank_surprising_associations.py --ontology go

    # Enzyme Commission
    uv run python scripts/rank_surprising_associations.py --ontology ec \
        --output results/domain_ec_surprising.tsv

Output
------
One TSV row per candidate, ranked by surprise score, with every component
exposed (observed and parts-only expected rate, lift, emergence p/q, region
overlap, novelty status, support counts) so results can be re-ranked without
recomputation.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from loguru import logger

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.domain_annotation_parser import DomainAnnotationParser  # noqa: E402
from src.hierarchy import closure_ancestors, parse_obo_child_parents  # noqa: E402
from src.ontology_registry import get_ontology, missing_inputs  # noqa: E402
from src.surprise_score import (  # noqa: E402
    EmergenceEvidence,
    conditional_rate,
    SurpriseResult,
    apply_fdr,
    locate_feature_regions,
    max_pairwise_overlap,
    median,
    novelty_factor,
    parse_interpro2go,
    proper_subfeatures,
    score_candidate,
)

#: Proteins sampled per feature when measuring constituent region overlap.
#: The overlap of a signature pair is a property of the signatures, so a sample
#: is enough; the median over it is stable well below this many proteins.
OVERLAP_SAMPLE = 25


def parse_associations(
    path: Path, term_column: str
) -> List[Tuple[str, str, float, Tuple[str, ...]]]:
    """Read the supra-domain rows of a significant-associations TSV.

    Returns ``(feature, term, q_value, constituent_domains)`` per row. Rows for
    single domains (``constituent_domains == "-"``) are skipped: a single domain
    has no parts, so it cannot be emergent.
    """
    rows: List[Tuple[str, str, float, Tuple[str, ...]]] = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_domain = header.index("domain")
            i_term = header.index(term_column)
            i_q = header.index("adj_p_value")
            i_parts = header.index("constituent_domains")
        except ValueError as exc:
            raise SystemExit(
                f"{path}: unexpected header {header} ({exc}). Expected the "
                f"columns written by run_dcgo_human.py (term column "
                f"{term_column!r})."
            ) from None

        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(i_domain, i_term, i_q, i_parts):
                continue
            constituents = fields[i_parts]
            if constituents == "-":
                continue
            rows.append(
                (
                    fields[i_domain],
                    fields[i_term],
                    float(fields[i_q]),
                    tuple(constituents.split(",")),
                )
            )
    return rows


def build_feature_index(
    architectures: Dict[str, object],
    universe: Set[str],
    wanted: Set[str],
) -> Dict[str, Set[str]]:
    """Map each wanted domain feature to the proteins carrying it.

    Only ``wanted`` features are materialised — the full supra-domain space is
    far larger than the candidate set and does not fit comfortably in memory.
    """
    index: Dict[str, Set[str]] = defaultdict(set)
    for protein in universe:
        arch = architectures[protein]
        for feature in arch.single_domains:
            if feature in wanted:
                index[feature].add(protein)
        for feature in arch.supra_domains:
            if feature in wanted:
                index[feature].add(protein)
    return index


def measure_region_overlap(
    feature: str,
    parts: Sequence[str],
    proteins: Iterable[str],
    architectures: Dict[str, object],
) -> float:
    """Median largest pairwise overlap between the combination's matched regions.

    High values flag the redundant-signature artefact: several InterPro
    signatures annotating one region, which looks like a domain combination in
    the architecture string but is not one biologically.
    """
    overlaps: List[float] = []
    for n_seen, protein in enumerate(proteins):
        if n_seen >= OVERLAP_SAMPLE:
            break
        arch = architectures[protein]
        annotations = arch.domain_annotations
        domain_ids = [a.interpro_id for a in annotations]
        intervals = [(a.start, a.end) for a in annotations]
        regions = locate_feature_regions(domain_ids, intervals, parts)
        if regions:
            overlaps.append(max_pairwise_overlap(regions))
    return median(overlaps)


def load_two_column_names(path: Optional[Path]) -> Dict[str, str]:
    """Read an ``id<TAB>name`` table (extra columns ignored, header tolerated)."""
    if path is None or not path.exists():
        return {}
    names: Dict[str, str] = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0] and not fields[0].startswith("#"):
                names.setdefault(fields[0], fields[1])
    return names


def load_obo_names(path: Optional[Path]) -> Dict[str, str]:
    """Read ``id: … / name: …`` pairs out of an OBO file."""
    if path is None or not path.exists():
        return {}
    names: Dict[str, str] = {}
    term_id: Optional[str] = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("["):
                term_id = None
            elif line.startswith("id:"):
                term_id = line[3:].strip()
            elif line.startswith("name:") and term_id:
                names.setdefault(term_id, line[5:].strip())
    return names


def load_interpro_names(path: Optional[Path]) -> Dict[str, str]:
    """Read InterPro ``entry.list`` (``ENTRY_AC<TAB>TYPE<TAB>NAME``)."""
    if path is None or not path.exists():
        return {}
    names: Dict[str, str] = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0].startswith("IPR"):
                names[fields[0]] = fields[2]
    return names


def build_novelty_reference(
    interpro2go: Path, ontology: str
) -> Tuple[Dict[str, Set[str]], str]:
    """Curated domain→term reference for the novelty factor, if one exists.

    Only GO has one here (InterPro2GO). Other ontologies score every prediction
    as ``"no-reference"``, which leaves the novelty factor at 1.0 — the score
    then ranks purely on emergence and distinctness.
    """
    if ontology != "go":
        return {}, f"no curated domain→term reference for --ontology {ontology}"
    if not interpro2go.exists():
        return {}, f"InterPro2GO not found at {interpro2go}; novelty factor disabled"
    with open(interpro2go) as f:
        mapping = parse_interpro2go(f)
    return mapping, f"InterPro2GO: {len(mapping):,} domains with curated terms"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rank emergent domain-combination predictions by surprise score",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ontology", default="go", help="Ontology the run used")
    parser.add_argument("--species", default="human")
    parser.add_argument(
        "--associations",
        type=Path,
        default=None,
        help="Significant-associations TSV (default: "
        "results/domain_<ontology>_associations_significant.tsv)",
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--term-column",
        default=None,
        help="Term column in the associations TSV (default: <ontology>_term). "
        "Needed for --ontology xref runs, whose outputs are named after the DR "
        "database, e.g. 'kegg_term'",
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=2,
        help="Minimum proteins carrying the combination AND annotated with the term",
    )
    parser.add_argument(
        "--max-overlap",
        type=float,
        default=0.5,
        help="Drop combinations whose constituent regions overlap more than this "
        "(redundant InterPro signatures for one region); 1.0 keeps everything",
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05, help="FDR level for the emergence test"
    )
    parser.add_argument(
        "--pseudo-count",
        type=float,
        default=1.0,
        help="Pseudo-observations shrinking each part's rate toward the term's "
        "background rate, so a domain seen in three proteins cannot claim a "
        "predictive rate of 1.0",
    )
    parser.add_argument(
        "--top", type=int, default=25, help="How many rows to print to the log"
    )
    parser.add_argument(
        "--enzyme-dat", type=Path, default=Path("data/raw/enzyme/enzyme.dat")
    )
    parser.add_argument(
        "--uniprot-dat",
        type=Path,
        default=Path("data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz"),
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--interpro2go", type=Path, default=Path("data/raw/interpro2go/interpro2go")
    )
    parser.add_argument(
        "--interpro-entries",
        type=Path,
        default=Path("data/raw/interpro_entry.list"),
        help="InterPro entry.list, for domain names in the output",
    )
    parser.add_argument(
        "--term-names",
        type=Path,
        default=None,
        help="Optional id<TAB>name table for term names (e.g. ReactomePathways.txt)",
    )
    parser.add_argument("--evidence-filter", default="manual")
    parser.add_argument(
        "--gaf",
        type=Path,
        default=None,
        help="GOA GAF to score against (default: the current release for "
        "--species). Point this at an archived release to score a historical "
        "run, e.g. a t0 snapshot for the temporal validation",
    )
    parser.add_argument(
        "--subcell", type=Path, default=Path("data/raw/uniprot_subcell/subcell.txt")
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
        "--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_lite.obo")
    )
    args = parser.parse_args()

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    ontology_label = args.ontology
    associations_path = args.associations or Path(
        f"results/domain_{ontology_label}_associations_significant.tsv"
    )
    output_path = args.output or Path(f"results/domain_{ontology_label}_surprising.tsv")
    if not associations_path.exists():
        logger.error(f"Associations file not found: {associations_path}")
        logger.error("Run run_dcgo_human.py for this ontology first.")
        return 1

    entry = get_ontology(args.ontology)
    paths = {
        "gaf": args.gaf or Path(f"data/raw/goa_annotations/goa_{args.species}.gaf.gz"),
        "go_obo": args.go_ontology,
        "enzyme_dat": args.enzyme_dat,
        "uniprot_dat": args.uniprot_dat,
        "reactome_relations": args.reactome_relations,
        "keywlist": args.keyword_list,
        "subcell": args.subcell,
        "chebi_obo": args.chebi_obo,
    }
    missing = missing_inputs(entry, paths)
    if missing:
        logger.error(
            f"Missing input(s) for --ontology {args.ontology}: " + "; ".join(missing)
        )
        return 1

    logger.info("=" * 70)
    logger.info(f"SURPRISE SCORE — emergent {ontology_label.upper()} predictions")
    logger.info("=" * 70)

    term_column = args.term_column or f"{ontology_label}_term"
    candidates = parse_associations(associations_path, term_column)
    logger.info(
        f"Supra-domain associations read from {associations_path}: {len(candidates):,}"
    )
    if not candidates:
        logger.error(
            "No supra-domain rows in the input — was the run made with "
            "--disable-supra-domains?"
        )
        return 1

    # Rebuild the run's universe: proteins that have both domains and terms.
    interpro_file = Path(f"data/interim/protein2ipr_{args.species}.dat.gz")
    if not interpro_file.exists():
        logger.error(f"Domain annotations not found: {interpro_file}")
        return 1
    annotations = entry.build_source(
        paths, {"evidence_filter": args.evidence_filter}
    ).parse()
    domain_parser = DomainAnnotationParser(
        max_supra_domain_length=3, min_domain_length=10
    )
    architectures = domain_parser.parse_protein2ipr_file(interpro_file)
    universe = set(annotations) & set(architectures)
    n_universe = len(universe)
    logger.info(f"Protein universe (domains ∩ annotations): {n_universe:,}")

    # Features whose protein sets the emergence test needs: each candidate
    # combination, its constituent single domains and its contained sub-combinations.
    wanted: Set[str] = set()
    subfeatures_of: Dict[str, List[str]] = {}
    for feature, _term, _q, parts in candidates:
        if feature in subfeatures_of:
            continue
        subs = proper_subfeatures(parts)
        subfeatures_of[feature] = subs
        wanted.add(feature)
        wanted.update(parts)
        wanted.update(subs)
    logger.info(f"Domain features to index: {len(wanted):,}")
    feature_proteins = build_feature_index(architectures, universe, wanted)

    # Proteins per candidate term, over the same universe.
    term_proteins: Dict[str, Set[str]] = defaultdict(set)
    wanted_terms = {term for _f, term, _q, _p in candidates}
    for protein in universe:
        for term in annotations[protein]:
            if term in wanted_terms:
                term_proteins[term].add(protein)

    curated, reference_note = build_novelty_reference(args.interpro2go, args.ontology)
    logger.info(f"Novelty reference: {reference_note}")
    ancestors_fn: Optional[Callable[[str], Iterable[str]]] = None
    if curated and args.go_ontology.exists():
        ancestors_fn = closure_ancestors(parse_obo_child_parents(args.go_ontology))

    logger.info("Scoring candidates...")
    scored: List[SurpriseResult] = []
    dropped_support = 0
    for feature, term, q_value, parts in candidates:
        carriers = feature_proteins.get(feature, set())
        annotated = term_proteins.get(term, set())
        n_both = len(carriers & annotated)
        if n_both < args.min_support:
            dropped_support += 1
            continue

        background = len(annotated) / n_universe if n_universe else 0.0

        def rate_of(sub_feature: str) -> Optional[float]:
            """Shrunken ``P(term | sub_feature)``, or None if nothing carries it."""
            sub_carriers = feature_proteins.get(sub_feature, set())
            if not sub_carriers:
                return None
            return conditional_rate(
                len(sub_carriers & annotated),
                len(sub_carriers),
                background,
                args.pseudo_count,
            )

        single_rates = [
            rate
            for part in dict.fromkeys(parts)  # de-duplicate repeated domains
            if (rate := rate_of(part)) is not None
        ]
        part_rates = [
            rate
            for sub in subfeatures_of[feature]
            if (rate := rate_of(sub)) is not None
        ]

        evidence = EmergenceEvidence(
            feature=feature,
            term=term,
            n_feature=len(carriers),
            n_both=n_both,
            single_rates=tuple(single_rates),
            part_rates=tuple(part_rates),
            background_rate=background,
            q_value=q_value,
        )
        overlap = measure_region_overlap(
            feature, parts, carriers & annotated, architectures
        )
        novelty, status = novelty_factor(
            term,
            set().union(*(curated.get(part, set()) for part in parts))
            if curated
            else set(),
            ancestors_fn,
        )
        scored.append(score_candidate(evidence, overlap, novelty, status))

    logger.info(f"  Scored: {len(scored):,} (dropped for support: {dropped_support:,})")
    scored = apply_fdr(scored, alpha=args.alpha)

    kept = [r for r in scored if r.region_overlap <= args.max_overlap]
    logger.info(
        f"  Distinct-region filter (overlap ≤ {args.max_overlap}): "
        f"{len(kept):,} kept, {len(scored) - len(kept):,} dropped as redundant signatures"
    )
    significant = [r for r in kept if r.q_emergence <= args.alpha]
    logger.info(
        f"  Emergent beyond their parts (q_emergence ≤ {args.alpha}): {len(significant):,}"
    )
    kept.sort(key=lambda r: -r.surprise)

    domain_names = load_interpro_names(args.interpro_entries)
    term_names = load_two_column_names(args.term_names)
    if args.ontology == "go":
        term_names = {**load_obo_names(args.go_ontology), **term_names}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(
            "rank\tsurprise\tdomain\tdomain_names\tterm\tterm_name\t"
            "n_feature\tn_both\tobserved_rate\texpected_rate\texpectation_source\t"
            "lift\tp_emergence\tq_emergence\tregion_overlap\tdistinctness\t"
            "novelty\tnovelty_status\tuninformative_constituents\tdcgo_adj_p_value\n"
        )
        for rank, r in enumerate(kept, 1):
            parts = r.feature.split(",")
            names = " + ".join(domain_names.get(p, p) for p in parts)
            f.write(
                f"{rank}\t{r.surprise:.3f}\t{r.feature}\t{names}\t{r.term}\t"
                f"{term_names.get(r.term, '')}\t{r.n_feature}\t{r.n_both}\t"
                f"{r.observed_rate:.4f}\t{r.expected_rate:.3e}\t{r.expectation_source}\t"
                f"{r.lift:.1f}\t{r.p_emergence:.3e}\t{r.q_emergence:.3e}\t"
                f"{r.region_overlap:.2f}\t{r.distinctness:.2f}\t{r.novelty:.2f}\t"
                f"{r.novelty_status}\t{r.uninformative_constituents}\t{r.q_value:.3e}\n"
            )
    logger.info(f"✓ Wrote {len(kept):,} ranked candidates to {output_path}")

    logger.info("")
    logger.info(f"Top {min(args.top, len(kept))} by surprise score:")
    for rank, r in enumerate(kept[: args.top], 1):
        names = " + ".join(domain_names.get(p, p) for p in r.feature.split(","))
        term_name = term_names.get(r.term, r.term)
        logger.info(
            f"  {rank:>3}. {r.surprise:6.1f}  {names} → {term_name} "
            f"[{r.term}] n={r.n_both}/{r.n_feature} lift={r.lift:.0f}× "
            f"q={r.q_emergence:.1e} {r.novelty_status}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
