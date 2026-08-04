#!/usr/bin/env python
"""Compare dcGO-2.0 against the published dcGO (Fang & Gough 2013) — §3.

The published dcGO is keyed by SCOP **sunid**; running our pipeline with
``--domain-key ssf`` keys it by the SUPERFAMILY signature ``SSF<sunid>``, so the
two sides live in the same domain namespace and the comparison is a join rather
than a re-derivation.

READ THE CONFOUNDS FIRST — they dominate any disagreement
--------------------------------------------------------
1. **Vintage.** Their tables are stamped 3 April 2016; our GOA is 2026. A decade
   of curation separates the two, in both directions (terms they could not have
   seen; annotations since retracted).
2. **Species scope.** They used all completely-sequenced genomes plus all of
   UniProt (>80M sequences); we use human only. Their evidence per domain is
   orders of magnitude larger, so **recall against them is not a measure of our
   method** — many of their pairs are simply unreachable in a human-only
   universe. This script quantifies that directly, by counting how many of their
   significant pairs have *zero* co-occurring human proteins. **Precision (how
   much of what we call significant they also call significant) is the
   interpretable headline; recall is reported but should not be read as
   quality.**
3. **Evidence policy.** Their IEA/evidence-code policy is not stated in the
   papers available to us. It is **unknown**, not "the same as ours". Run this
   script against pipeline runs made with both ``--evidence-filter manual`` and
   ``--evidence-filter all`` to bracket it.
4. **Reachable domain space.** InterPro integrates SUPERFAMILY at superfamily
   level only, so their ``fa`` (family) half is structurally out of reach: only
   the ``sf`` rows can be compared at all.
5. **Not a full test matrix.** Their ``GO_mapping`` is not exhaustive over
   sf × GO. A pair absent from their table was not necessarily tested and found
   non-significant, so "absent" is reported as its own category and never folded
   into "they disagree".

What it reports
---------------
* Agreement over the shared (sunid, GO) space, at their FDR threshold (1e-3) and
  ours (0.01), overall and per GO aspect.
* For every pair they call significant, **our recomputed Fisher p-value** — so a
  pair we miss narrowly is distinguished from one we assign p≈1 and from one no
  human protein could support. And for every pair we call significant, **their**
  FDR — likewise split into "they tested it and disagreed" vs "not in their
  table".
* Spearman correlation of our −log10(p) against their FDR and h-score over the
  shared pairs (a calibration result, which also feeds §5).
* The same agreement for **supra-domains** against ``SP2GO.txt`` — the only
  external reference for our supra-domain machinery anywhere in the plan.

Usage
-----
    uv run python run_dcgo_human.py --domain-key ssf --output-dir results_ssf
    uv run python scripts/download_data.py --group dcgo-reference
    uv run python validation/compare_original_dcgo.py \
        --associations results_ssf/domain_ssf_go_associations_significant.tsv \
        --output-dir validation
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
from loguru import logger
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.annotation_source import GAFAnnotationSource  # noqa: E402
from src.domain_annotation_parser import (  # noqa: E402
    DomainAnnotationParser,
    superfamily_sunid,
)
from src.vectorized_fisher import fisher_exact_parallel  # noqa: E402

logger.remove()
logger.add(sys.stderr, level="INFO")

#: One complete ``GO_mapping`` row inside an INSERT statement:
#: ``(id,'level',go,single_score,all_score,'inherited_from','inherited_from_all',
#: all_fdr_min,all_hscore_max)``. Parsing the dump with a regex avoids needing a
#: MySQL server; this matches all 841,270 rows, so nothing is silently skipped.
_GO_MAPPING_ROW = re.compile(
    r"\((\d+),'(sf|fa)',(\d{7}),([^,]+),([^,]+),'([^']*)','([^']*)',([^,]+),([^)]+)\)"
)

ASPECTS = ("molecular_function", "biological_process", "cellular_component")


# --------------------------------------------------------------------------
# Reference-side parsing
# --------------------------------------------------------------------------


def parse_go_mapping(
    sql_gz: Path,
) -> Tuple[
    Dict[Tuple[int, str], float], Dict[Tuple[int, str], float], Set[int], Set[int]
]:
    """Read ``GO_mapping`` out of the published MySQL dump.

    Returns ``(fdr, hscore, sf_sunids, fa_sunids)`` for the ``sf`` level only —
    ``fa`` (SCOP family) is unreachable through InterPro, which integrates
    SUPERFAMILY at superfamily level, so only its domain *count* is kept, to
    quantify that confound.

    ``fdr`` is the ``all_score`` column: the FDR computed over all proteins,
    including multi-domain ones. That is the right comparator for us — we impose
    no single-domain restriction, unlike their ``single_score``. ``hscore`` holds
    ``all_hscore_max`` where the dump provides it.
    """
    fdr: Dict[Tuple[int, str], float] = {}
    hscore: Dict[Tuple[int, str], float] = {}
    sf_sunids: Set[int] = set()
    fa_sunids: Set[int] = set()
    with gzip.open(sql_gz, "rt", encoding="latin-1") as fh:
        for line in fh:
            if not line.startswith("INSERT INTO `GO_mapping`"):
                continue
            for match in _GO_MAPPING_ROW.finditer(line):
                sunid = int(match.group(1))
                if match.group(2) != "sf":
                    fa_sunids.add(sunid)
                    continue
                sf_sunids.add(sunid)
                key = (sunid, f"GO:{match.group(3)}")
                fdr[key] = float(match.group(5))
                raw = match.group(9).strip()
                if raw != "NULL":
                    hscore[key] = float(raw)
    logger.info(
        f"Published GO_mapping (sf): {len(fdr):,} scored pairs over "
        f"{len(sf_sunids):,} superfamilies, {len(hscore):,} with an h-score"
    )
    return fdr, hscore, sf_sunids, fa_sunids


def parse_domain2go_flat(
    path: Path, level: str = "sf"
) -> Tuple[Set[Tuple[int, str]], Set[Tuple[int, str]]]:
    """Read the published association *set* (``Domain2GO_supported_only_by_all``).

    This — not ``all_score < 1e-3`` — is what dcGO actually publishes, and the
    two are **not** the same set. At the ``sf`` level, of 404,288 scored pairs,
    134,665 sit below their 1e-3 threshold, but they ship 108,612 direct +
    85,683 inherited:

    * every direct row bar one is below 1e-3, so *direct ⊂ significant*;
    * **26,054 significant pairs (19.3%) are not shipped at all** — removed by
      the further filtering their method applies on top of the FDR (the
      relative/parental-background test and most-specific-level selection);
    * **none** of the inherited rows is itself below 1e-3 — they are pure
      true-path propagation onto ancestor terms.

    Comparing only against a re-thresholded FDR would therefore be comparing
    against something they never shipped, which is why the caller reports both.

    Returns ``(direct, inherited)`` sets of ``(sunid, GO)``. ``direct`` is the
    like-for-like comparator for a run made without ``--enable-true-path``.
    """
    direct: Set[Tuple[int, str]] = set()
    inherited: Set[Tuple[int, str]] = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Domain_type"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or fields[0] != level:
                continue
            try:
                key = (int(fields[1]), fields[2])
            except ValueError:
                continue
            (direct if fields[6] == "1" else inherited).add(key)
    logger.info(
        f"Published association set ({level}): {len(direct):,} direct + "
        f"{len(inherited):,} true-path-inherited"
    )
    return direct, inherited


def parse_sp2go(path: Path, direct_only: bool) -> Dict[Tuple[str, ...], Set[str]]:
    """Read the published supra-domain associations (``SP2GO.txt``).

    Returns ``{sunid_tuple: {GO term}}``. Rows with a single component are kept
    out — those are single domains wearing a supra-domain label, and comparing
    them here would double-count the single-domain table.

    ``direct_only`` keeps ``Annotation_origin == 1``. A run made without
    ``--enable-true-path`` produces direct associations only, so mixing their
    inherited rows in would inflate their side of every count.
    """
    supra: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Type\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or not fields[0].startswith("supra"):
                continue
            if direct_only and fields[6] != "1":
                continue
            parts = tuple(fields[1].split(","))
            if len(parts) < 2:
                continue
            supra[parts].add(fields[2])
    logger.info(
        f"Published SP2GO ({'direct only' if direct_only else 'direct + inherited'}): "
        f"{len(supra):,} supra-domain architectures, "
        f"{sum(len(v) for v in supra.values()):,} associations"
    )
    return supra


def parse_scop_des(path: Path) -> Dict[int, Tuple[str, str, str]]:
    """``dir.des.scop.1.75.txt`` → ``{sunid: (type, sccs, description)}``."""
    des: Dict[int, Tuple[str, str, str]] = {}
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                des[int(fields[0])] = (fields[1], fields[2], fields[4])
            except ValueError:
                continue
    return des


def parse_go_namespaces(obo: Path) -> Tuple[Dict[str, str], Set[str]]:
    """``({GO id (incl. alt_ids): namespace}, {obsolete ids})`` from a go-basic.obo.

    Obsolete terms are still carried in the OBO with a namespace, so "present in
    the file" is not the same as "still a usable term" — the vintage confound
    needs the obsolete set called out separately.
    """
    namespaces: Dict[str, str] = {}
    obsolete: Set[str] = set()
    ids: List[str] = []
    namespace: Optional[str] = None
    is_obsolete = False

    def flush() -> None:
        for term_id in ids:
            if namespace:
                namespaces[term_id] = namespace
            if is_obsolete:
                obsolete.add(term_id)

    with open(obo) as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("["):
                flush()
                ids, namespace, is_obsolete = [], None, False
            elif line.startswith("id: GO:"):
                ids.append(line[4:].strip())
            elif line.startswith("alt_id: GO:"):
                ids.append(line[8:].strip())
            elif line.startswith("namespace:"):
                namespace = line[10:].strip()
            elif line.startswith("is_obsolete: true"):
                is_obsolete = True
    flush()
    return namespaces, obsolete


# --------------------------------------------------------------------------
# Our side
# --------------------------------------------------------------------------


def parse_our_associations(
    path: Path, term_column: str = "go_term"
) -> Tuple[Dict[Tuple[str, str], Tuple[float, float]], Dict[Tuple[str, ...], Set[str]]]:
    """Split a significant-associations TSV into single-domain and supra rows.

    Returns ``({(domain, term): (p, q)}, {domain_tuple: {term}})``.
    """
    singles: Dict[Tuple[str, str], Tuple[float, float]] = {}
    supra: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for name in ("domain", term_column, "p_value", "adj_p_value", "domain_type"):
            if name not in idx:
                raise SystemExit(f"{path}: missing column {name!r} (header {header})")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            domain = fields[idx["domain"]]
            term = fields[idx[term_column]]
            p_value = float(fields[idx["p_value"]])
            q_value = float(fields[idx["adj_p_value"]])
            if fields[idx["domain_type"]] == "single":
                singles[(domain, term)] = (p_value, q_value)
            else:
                supra[tuple(domain.split(","))].add(term)
    logger.info(
        f"Our run: {len(singles):,} significant single-domain associations, "
        f"{sum(len(v) for v in supra.values()):,} supra-domain associations "
        f"over {len(supra):,} architectures"
    )
    return singles, supra


def build_universe(
    interpro_file: Path, gaf: Path, evidence_filter: str
) -> Tuple[
    Dict[str, Set[str]], Dict[str, Set[str]], int, Dict[Tuple[str, ...], Set[str]]
]:
    """Rebuild the SSF-keyed protein universe the run used.

    Returns ``(domain → proteins, term → proteins, |universe|,
    architecture → proteins)``. The last one covers supra-domains, needed for the
    ``SP2GO`` comparison.
    """
    parser = DomainAnnotationParser(
        max_supra_domain_length=3, min_domain_length=10, domain_key="ssf"
    )
    architectures = parser.parse_protein2ipr_file(interpro_file)
    annotations = GAFAnnotationSource(gaf, evidence_filter=evidence_filter).parse()

    universe = set(architectures) & set(annotations)
    logger.info(f"Protein universe (SSF domains ∩ GO annotations): {len(universe):,}")

    domain_proteins: Dict[str, Set[str]] = defaultdict(set)
    supra_proteins: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)
    term_proteins: Dict[str, Set[str]] = defaultdict(set)
    for protein in universe:
        arch = architectures[protein]
        for domain in arch.single_domains:
            domain_proteins[domain].add(protein)
        for feature in arch.supra_domains:
            supra_proteins[tuple(feature.split(","))].add(protein)
        for term in annotations[protein]:
            term_proteins[term].add(protein)
    return domain_proteins, term_proteins, len(universe), supra_proteins


def fisher_p_for_pairs(
    pairs: Sequence[Tuple[str, str]],
    domain_proteins: Dict[str, Set[str]],
    term_proteins: Dict[str, Set[str]],
    n_universe: int,
    n_jobs: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Right-tailed Fisher p-value and overlap count ``a`` for arbitrary pairs.

    Identical construction to ``src/sparse_fisher.compute_contingency_tables_sparse``
    — the pipeline reports p-values only for pairs that passed FDR, but the
    interesting question here is what we say about the pairs we *did not* report.
    ``verify_recomputation`` checks this against the pipeline's own output.
    """
    tables = np.empty((len(pairs), 2, 2), dtype=np.int64)
    overlaps = np.empty(len(pairs), dtype=np.int64)
    for i, (domain, term) in enumerate(pairs):
        carriers = domain_proteins.get(domain, set())
        annotated = term_proteins.get(term, set())
        a = len(carriers & annotated)
        tables[i, 0, 0] = a
        tables[i, 0, 1] = len(carriers) - a
        tables[i, 1, 0] = len(annotated) - a
        tables[i, 1, 1] = n_universe - len(carriers) - len(annotated) + a
        overlaps[i] = a
    _odds, pvalues = fisher_exact_parallel(
        tables, alternative="greater", n_jobs=n_jobs, batch_size=200_000
    )
    return pvalues, overlaps


def verify_recomputation(
    ours: Dict[Tuple[str, str], Tuple[float, float]],
    domain_proteins: Dict[str, Set[str]],
    term_proteins: Dict[str, Set[str]],
    n_universe: int,
    n_jobs: int,
    sample: int = 2000,
) -> float:
    """Max relative deviation between recomputed and pipeline-reported p-values.

    A non-trivial deviation means the universe rebuilt here is not the one the
    run used, which would invalidate every "our p-value for a pair we did not
    report" number below. Reported, not silently assumed.
    """
    pairs = sorted(ours)[:sample]
    if not pairs:
        return 0.0
    recomputed, _ = fisher_p_for_pairs(
        pairs, domain_proteins, term_proteins, n_universe, n_jobs
    )
    reported = np.array([ours[pair][0] for pair in pairs])
    with np.errstate(divide="ignore", invalid="ignore"):
        deviation = np.abs(recomputed - reported) / np.maximum(reported, 1e-300)
    # The TSV stores p-values at 6 significant figures, so agreement is expected
    # to ~1e-6 relative, not exactly.
    return float(np.nanmax(np.where(reported > 0, deviation, 0.0)))


# --------------------------------------------------------------------------
# Reporting helpers
# --------------------------------------------------------------------------


def write_tsv(
    path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")
    logger.info(f"✓ wrote {path}")


def fraction(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.4f}" if denominator else ""


def agreement(ours: Set[object], theirs: Set[object]) -> Tuple[int, str, str, str]:
    shared = ours & theirs
    union = ours | theirs
    return (
        len(shared),
        fraction(len(shared), len(ours)),
        fraction(len(shared), len(theirs)),
        fraction(len(shared), len(union)),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare dcGO-2.0 (SSF-keyed) against the published dcGO",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--associations",
        type=Path,
        default=Path("results_ssf/domain_ssf_go_associations_significant.tsv"),
        help="Significant associations from a --domain-key ssf run",
    )
    parser.add_argument(
        "--dcgo-sql",
        type=Path,
        default=Path("data/raw/dcgo_reference/Domain2GO.sql.gz"),
    )
    parser.add_argument(
        "--dcgo-flat",
        type=Path,
        default=Path("data/raw/dcgo_reference/Domain2GO_supported_only_by_all.txt"),
        help="Their published association set (the 'all proteins' flavour, which "
        "matches our lack of a single-domain restriction)",
    )
    parser.add_argument(
        "--dcgo-supra",
        type=Path,
        default=Path("data/raw/dcgo_reference/SP2GO.txt"),
    )
    parser.add_argument(
        "--scop-des", type=Path, default=Path("data/raw/scop/dir.des.scop.1.75.txt")
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--gaf", type=Path, default=Path("data/raw/goa_annotations/goa_human.gaf.gz")
    )
    parser.add_argument(
        "--interpro", type=Path, default=Path("data/interim/protein2ipr_human.dat.gz")
    )
    parser.add_argument("--evidence-filter", default="manual")
    parser.add_argument(
        "--their-fdr",
        type=float,
        default=1e-3,
        help="Threshold the published dcGO used (Fang & Gough 2013)",
    )
    parser.add_argument(
        "--our-fdr", type=float, default=0.01, help="Threshold our run used"
    )
    parser.add_argument("--num-cores", type=int, default=8)
    parser.add_argument("--output-dir", type=Path, default=Path("validation"))
    parser.add_argument(
        "--label",
        default="manual",
        help="Suffix for the output filenames, so evidence-filter variants do "
        "not overwrite each other",
    )
    args = parser.parse_args()

    metrics: List[Tuple[str, object, str]] = []

    def record(key: str, value: object, note: str = "") -> None:
        metrics.append((key, value, note))
        logger.info(f"  {key}: {value}  {note}")

    logger.info("=" * 72)
    logger.info("VALIDATION §3 — dcGO-2.0 vs the published dcGO (Fang & Gough 2013)")
    logger.info("=" * 72)

    # ---- load both sides -------------------------------------------------
    their_fdr, their_hscore, their_sf_domains, their_fa_domains = parse_go_mapping(
        args.dcgo_sql
    )
    their_direct, their_inherited = parse_domain2go_flat(args.dcgo_flat)

    ours_single, ours_supra = parse_our_associations(args.associations)
    namespaces, obsolete_terms = parse_go_namespaces(args.go_ontology)
    scop = parse_scop_des(args.scop_des)

    domain_proteins, term_proteins, n_universe, supra_proteins = build_universe(
        args.interpro, args.gaf, args.evidence_filter
    )

    # ---- confound 4: what fraction of their domain space is reachable ----
    our_domains = {
        sunid
        for sunid in (superfamily_sunid(d) for d in domain_proteins)
        if sunid is not None
    }
    shared_domains = our_domains & their_sf_domains
    record("published_domains_sf", len(their_sf_domains))
    record("published_domains_fa", len(their_fa_domains), "unreachable via InterPro")
    record(
        "published_domains_reachable_fraction",
        fraction(len(their_sf_domains), len(their_sf_domains | their_fa_domains)),
        "InterPro integrates SUPERFAMILY at superfamily level only",
    )
    record("our_domains_ssf", len(our_domains))
    record("shared_domains", len(shared_domains))
    record(
        "our_domains_absent_from_theirs",
        len(our_domains - their_sf_domains),
        "; ".join(
            f"{sunid} ({scop.get(sunid, ('', '', '?'))[2]})"
            for sunid in sorted(our_domains - their_sf_domains)
        ),
    )

    # ---- confound 1: GO vintage -----------------------------------------
    their_terms = {term for _sunid, term in their_fdr}
    our_terms = set(term_proteins)
    record("published_go_terms", len(their_terms))
    record("our_go_terms", len(our_terms))
    record("shared_go_terms", len(their_terms & our_terms))
    record(
        "published_go_terms_absent_from_current_obo",
        len(their_terms - set(namespaces)),
        "gone entirely since their 2016 snapshot",
    )
    record(
        "published_go_terms_now_obsolete",
        len(their_terms & obsolete_terms),
        "confound 1 (vintage): terms they scored that GO has since retired",
    )
    record(
        "our_go_terms_absent_from_their_table",
        len(our_terms - their_terms),
        "terms they never associated with any superfamily — no comparison possible",
    )

    shared_terms = their_terms & our_terms

    def in_shared_space(sunid: Optional[int], term: str) -> bool:
        return sunid in shared_domains and term in shared_terms

    # ---- single-domain agreement ----------------------------------------
    ours_pairs = {(superfamily_sunid(domain), term) for domain, term in ours_single}
    ours_pairs = {
        (sunid, term) for sunid, term in ours_pairs if in_shared_space(sunid, term)
    }
    # Three defensible readings of "what the original called significant". They
    # are reported together rather than the flattering one being picked:
    #
    #   published_direct    what they actually shipped, direct associations only
    #                       — the like-for-like comparator for a run made without
    #                       --enable-true-path, and therefore the PRIMARY one;
    #   published_all       + their true-path-inherited rows;
    #   fdr_lt_1e-3         a naive re-threshold of GO_mapping.all_score at the
    #                       paper's cutoff. This is NOT their published set: it
    #                       both adds 63,551 pairs they did not ship and drops
    #                       174,375 they did, because their release applies
    #                       further filtering on top of the FDR.
    their_sets = {
        "published_direct": {p for p in their_direct if in_shared_space(*p)},
        "published_all": {
            p for p in (their_direct | their_inherited) if in_shared_space(*p)
        },
        "fdr_lt_%g" % args.their_fdr: {
            pair
            for pair, value in their_fdr.items()
            if value < args.their_fdr and in_shared_space(*pair)
        },
    }
    primary = "published_direct"
    theirs_sig = their_sets[primary]

    # How their shipped set relates to a naive re-threshold of the FDR column.
    # Recorded because getting this wrong silently changes every number below.
    their_below = {pair for pair, value in their_fdr.items() if value < args.their_fdr}
    record("their_scored_pairs_sf", len(their_fdr))
    record(
        "their_pairs_below_their_fdr",
        len(their_below),
        f"all_score < {args.their_fdr:g}",
    )
    record("their_published_direct_sf", len(their_direct))
    record("their_published_inherited_sf", len(their_inherited))
    record(
        "their_direct_rows_below_their_fdr",
        len(their_direct & their_below),
        "their direct set is a strict SUBSET of the FDR-significant set",
    )
    record(
        "their_significant_but_unpublished",
        len(their_below - their_direct - their_inherited),
        "significant by FDR yet not shipped — their release applies further "
        "filtering (most-specific-level selection / the relative test)",
    )
    record(
        "their_inherited_rows_below_their_fdr",
        len(their_inherited & their_below),
        "0 expected: inherited rows are pure true-path propagation to ancestors, "
        "not separately significant",
    )
    record(
        "their_scored_rows_at_column_default_1.0",
        sum(1 for value in their_fdr.values() if value == 1.0),
        "the DEFAULT of all_score is 1 — these are indistinguishable from "
        "'never actually tested', which caps how much the rank correlation below "
        "can be trusted",
    )

    # Their table is not an exhaustive sf × GO matrix: absence is "untested",
    # not "tested and rejected". Kept as its own category throughout.
    their_tested = {pair for pair in their_fdr if in_shared_space(*pair)}

    record("testable_space_pairs", len(shared_domains) * len(shared_terms))
    record("their_pairs_tested_in_shared_space", len(their_tested))
    record(
        "their_table_coverage_of_testable_space",
        fraction(len(their_tested), len(shared_domains) * len(shared_terms)),
        "their GO_mapping is sparse — absence ≠ tested-and-rejected",
    )

    # How much of our own output the shared space even covers — the reachability
    # cost of the comparison, before any agreement number is quoted.
    ours_all = {(superfamily_sunid(domain), term) for domain, term in ours_single}
    record("our_significant_single_domain_pairs_total", len(ours_all))
    record(
        "our_pairs_outside_the_shared_space",
        len(ours_all) - len({p for p in ours_all if in_shared_space(*p)}),
        "domain or term absent from their table; not comparable either way",
    )

    # Our side at our FDR (0.01) and at theirs (1e-3), against each of their
    # three definitions — nine numbers, all reported, none selected after the fact.
    ours_strict = {
        pair
        for pair in (
            (superfamily_sunid(domain), term)
            for (domain, term), (_p, q) in ours_single.items()
            if q < args.their_fdr
        )
        if in_shared_space(*pair)
    }
    variant_rows = []
    for their_name, their_set in their_sets.items():
        for our_name, our_set in (
            (f"fdr<{args.our_fdr:g}", ours_pairs),
            (f"fdr<{args.their_fdr:g}", ours_strict),
        ):
            n_v, prec_v, rec_v, jac_v = agreement(our_set, their_set)
            variant_rows.append(
                [
                    our_name,
                    their_name,
                    len(our_set),
                    len(their_set),
                    n_v,
                    prec_v,
                    rec_v,
                    jac_v,
                ]
            )

    n_shared, precision, recall, jaccard = agreement(ours_pairs, theirs_sig)
    record("ours_significant_pairs", len(ours_pairs))
    record("their_significant_pairs", len(theirs_sig), f"definition: {primary}")
    record("shared_significant_pairs", n_shared)
    record("precision_ours_in_theirs", precision, "HEADLINE — see confound 2")
    record("recall_theirs_in_ours", recall, "depressed by species scope; not quality")
    record("jaccard", jaccard)
    for row in variant_rows:
        record(
            f"variant__ours_{row[0]}__theirs_{row[1]}",
            f"precision={row[5]} recall={row[6]} jaccard={row[7]}",
            f"ours={row[2]:,} theirs={row[3]:,} shared={row[4]:,}",
        )

    # ---- per-aspect ------------------------------------------------------
    aspect_rows = []
    for aspect in ASPECTS:
        ours_a = {p for p in ours_pairs if namespaces.get(p[1]) == aspect}
        theirs_a = {p for p in theirs_sig if namespaces.get(p[1]) == aspect}
        n_a, prec_a, rec_a, jac_a = agreement(ours_a, theirs_a)
        aspect_rows.append(
            [aspect, len(ours_a), len(theirs_a), n_a, prec_a, rec_a, jac_a]
        )

    # ---- diagnosis: pairs they call significant that we do not -----------
    missed = sorted(theirs_sig - ours_pairs)
    logger.info(f"Recomputing our p-value for {len(missed):,} pairs we did not report")
    # Map sunid back to the SSF accession our tables are keyed by.
    ssf_of = {superfamily_sunid(d): d for d in domain_proteins}
    missed_pairs = [(ssf_of[sunid], term) for sunid, term in missed]
    if missed_pairs:
        missed_p, missed_a = fisher_p_for_pairs(
            missed_pairs, domain_proteins, term_proteins, n_universe, args.num_cores
        )
    else:
        missed_p, missed_a = np.array([]), np.array([])

    # The effective raw-p cutoff of our run: the largest raw p that passed BH.
    our_p_cutoff = max((p for p, _q in ours_single.values()), default=0.0)
    record("our_bh_raw_p_cutoff", f"{our_p_cutoff:.3e}")

    # Three mutually exclusive buckets, so the reader can see at a glance which
    # kind of "miss" dominates. A pair we miss narrowly is a very different
    # finding from one we assign p≈1, and both differ from one that no human
    # protein could ever support.
    supported_missed = missed_a > 0
    missed_zero_overlap = int((~supported_missed).sum())
    missed_near = int((supported_missed & (missed_p < 0.05)).sum())
    missed_flat = int((supported_missed & (missed_p >= 0.05)).sum())
    record("missed_pairs", len(missed))
    record(
        "missed_zero_human_cooccurrence",
        missed_zero_overlap,
        "no human protein carries both — unreachable, confound 2, not a method miss",
    )
    record(
        "missed_zero_human_cooccurrence_fraction",
        fraction(missed_zero_overlap, len(missed)),
    )
    record(
        "missed_supported_p_below_0.05",
        missed_near,
        "signal present in human, just under our FDR bar",
    )
    record(
        "missed_supported_p_above_0.05",
        missed_flat,
        "some human co-occurrence, but no enrichment we can see",
    )
    if len(missed_p) and supported_missed.any():
        record(
            "missed_supported_median_our_p",
            f"{float(np.median(missed_p[supported_missed])):.3e}",
        )

    # ---- diagnosis: pairs we call significant that they do not -----------
    extra = sorted(ours_pairs - theirs_sig)
    extra_untested = [pair for pair in extra if pair not in their_fdr]
    extra_tested = [pair for pair in extra if pair in their_fdr]
    extra_inherited = [pair for pair in extra if pair in their_inherited]
    extra_scores = np.array([their_fdr[pair] for pair in extra_tested])
    record("extra_pairs", len(extra))
    record(
        "extra_not_in_their_table",
        len(extra_untested),
        "they never scored this pair — untested, not rejected",
    )
    record("extra_they_scored_but_did_not_publish_as_direct", len(extra_tested))
    record(
        "extra_they_published_as_inherited",
        len(extra_inherited),
        "same association, assigned to a different level of the GO DAG — "
        "agreement in substance, disagreement in placement",
    )
    if len(extra_scores):
        record(
            f"extra_their_fdr_between_{args.their_fdr:g}_and_0.05",
            int(((extra_scores >= args.their_fdr) & (extra_scores < 0.05)).sum()),
            "they nearly called it too",
        )
        record(
            "extra_their_fdr_below_their_threshold",
            int((extra_scores < args.their_fdr).sum()),
            "they scored it significant but did not ship it (post-FDR filtering)",
        )
        record("extra_their_fdr_above_0.5", int((extra_scores > 0.5).sum()))
        record("extra_median_their_fdr", f"{float(np.median(extra_scores)):.3e}")

    # ---- rank correlation over everything they tested --------------------
    joint = sorted(their_tested)
    joint_pairs = [(ssf_of[sunid], term) for sunid, term in joint]
    logger.info(f"Recomputing our p-value for all {len(joint):,} pairs they tested")
    joint_p, joint_a = fisher_p_for_pairs(
        joint_pairs, domain_proteins, term_proteins, n_universe, args.num_cores
    )
    our_score = -np.log10(np.maximum(joint_p, 1e-300))
    their_score = -np.log10(np.maximum(np.array([their_fdr[p] for p in joint]), 1e-300))
    rho_all = spearmanr(our_score, their_score)
    record("spearman_our_logp_vs_their_logfdr_all_tested", f"{rho_all.statistic:.4f}")
    record("spearman_p", f"{rho_all.pvalue:.3e}")

    # Restricted to pairs with any human support at all — the correlation above
    # is partly driven by the mass of a=0 pairs that we structurally cannot see.
    supported = joint_a > 0
    rho_supported = spearmanr(our_score[supported], their_score[supported])
    record("pairs_with_human_cooccurrence", int(supported.sum()))
    record(
        "spearman_restricted_to_supported_pairs",
        f"{rho_supported.statistic:.4f}",
        "removes the a=0 mass that is a species-scope artefact",
    )

    # ... and additionally dropping their column-default 1.0 rows, which cannot
    # be distinguished from "never tested". This is the cleanest reading of the
    # calibration: both sides have a real number for every pair in it.
    their_raw = np.array([their_fdr[p] for p in joint])
    scored_and_supported = supported & (their_raw < 1.0)
    if scored_and_supported.sum() > 2:
        rho_clean = spearmanr(
            our_score[scored_and_supported], their_score[scored_and_supported]
        )
        record(
            "pairs_supported_and_not_at_their_default", int(scored_and_supported.sum())
        )
        record(
            "spearman_supported_and_genuinely_scored",
            f"{rho_clean.statistic:.4f}",
            "cleanest calibration comparison",
        )

    hs_idx = [i for i, pair in enumerate(joint) if pair in their_hscore]
    if len(hs_idx) > 2:
        hs_values = np.array([their_hscore[joint[i]] for i in hs_idx])
        rho_hs = spearmanr(our_score[hs_idx], hs_values)
        record("pairs_with_their_hscore", len(hs_idx))
        record(
            "spearman_our_logp_vs_their_hscore",
            f"{rho_hs.statistic:.4f}",
            "all_hscore_max is populated only on their published rows and is a "
            "max over the GO subtree, so it is not the same quantity as our p",
        )
        hs_supported = [i for i in hs_idx if joint_a[i] > 0]
        if len(hs_supported) > 2:
            rho_hs_sup = spearmanr(
                our_score[hs_supported],
                np.array([their_hscore[joint[i]] for i in hs_supported]),
            )
            record(
                "spearman_our_logp_vs_their_hscore_supported",
                f"{rho_hs_sup.statistic:.4f}",
            )

    # ---- honesty check: is our rebuilt universe the run's universe? -------
    deviation = verify_recomputation(
        ours_single, domain_proteins, term_proteins, n_universe, args.num_cores
    )
    record(
        "recomputation_max_relative_deviation",
        f"{deviation:.2e}",
        "vs the pipeline's own p-values; >1e-5 would invalidate the diagnostics",
    )

    # ---- supra-domains ---------------------------------------------------
    # Direct rows only, to match a run without --enable-true-path (their
    # inherited rows would otherwise inflate their side ~5×).
    their_supra = parse_sp2go(args.dcgo_supra, direct_only=True)
    their_supra_all = parse_sp2go(args.dcgo_supra, direct_only=False)

    # Their tuples are bare sunids; ours are SSF accessions in N→C order.
    def to_sunids(feature: Tuple[str, ...]) -> Optional[Tuple[str, ...]]:
        out = []
        for part in feature:
            sunid = superfamily_sunid(part)
            if sunid is None:
                return None
            out.append(str(sunid))
        return tuple(out)

    our_arch = {
        key: value
        for key, value in (
            (to_sunids(feature), terms) for feature, terms in ours_supra.items()
        )
        if key is not None
    }
    our_all_arch = {
        key
        for key in (to_sunids(feature) for feature in supra_proteins)
        if key is not None
    }

    record("our_supra_architectures_observed", len(our_all_arch))
    record("our_supra_architectures_with_associations", len(our_arch))
    record("their_supra_architectures", len(their_supra))

    n_arch, prec_arch, rec_arch, jac_arch = agreement(our_all_arch, set(their_supra))
    record("shared_supra_architectures_exact_order", n_arch)
    record("supra_architecture_precision_ours_in_theirs", prec_arch)

    # Order convention is not documented on their side, so report the
    # order-insensitive variant too rather than assuming ours matches.
    def canon(t: Tuple[str, ...]) -> Tuple[str, ...]:
        return tuple(sorted(t))

    our_arch_canon = {canon(t) for t in our_all_arch}
    their_arch_canon = {canon(t) for t in their_supra}
    n_arch_c, prec_arch_c, rec_arch_c, jac_arch_c = agreement(
        our_arch_canon, their_arch_canon
    )
    record("shared_supra_architectures_order_insensitive", n_arch_c)
    record("supra_architecture_precision_order_insensitive", prec_arch_c)

    ours_supra_pairs = {
        (arch, term) for arch, terms in our_arch.items() for term in terms
    }
    theirs_supra_pairs = {
        (arch, term) for arch, terms in their_supra.items() for term in terms
    }
    # Restrict to architectures both sides observed, and to shared GO terms —
    # otherwise the comparison is dominated by architectures that simply do not
    # occur in human.
    common_arch = {a for a in our_arch if a in their_supra}
    ours_supra_r = {
        (a, t) for a, t in ours_supra_pairs if a in common_arch and t in shared_terms
    }
    theirs_supra_r = {
        (a, t) for a, t in theirs_supra_pairs if a in common_arch and t in shared_terms
    }
    n_sp, prec_sp, rec_sp, jac_sp = agreement(ours_supra_r, theirs_supra_r)
    record("supra_architectures_compared", len(common_arch))
    record("our_supra_associations_compared", len(ours_supra_r))
    record("their_supra_associations_compared", len(theirs_supra_r))
    record("shared_supra_associations", n_sp)
    record("supra_precision_ours_in_theirs", prec_sp)
    record("supra_recall_theirs_in_ours", rec_sp)
    record("supra_jaccard", jac_sp)

    # Same, against their direct + inherited rows — the variant a reader would
    # otherwise ask for, reported rather than left out.
    theirs_supra_all_r = {
        (a, t)
        for a, terms in their_supra_all.items()
        if a in common_arch
        for t in terms
        if t in shared_terms
    }
    n_sp_all, prec_sp_all, rec_sp_all, jac_sp_all = agreement(
        ours_supra_r, theirs_supra_all_r
    )
    record("supra_precision_vs_their_direct_plus_inherited", prec_sp_all)
    record("supra_recall_vs_their_direct_plus_inherited", rec_sp_all)

    # ---- write it out ----------------------------------------------------
    out = args.output_dir
    write_tsv(
        out / f"dcgo_comparison_metrics_{args.label}.tsv",
        ["metric", "value", "note"],
        metrics,
    )
    write_tsv(
        out / f"dcgo_comparison_variants_{args.label}.tsv",
        [
            "our_threshold",
            "their_definition",
            "ours",
            "theirs",
            "shared",
            "precision_ours_in_theirs",
            "recall_theirs_in_ours",
            "jaccard",
        ],
        variant_rows,
    )
    write_tsv(
        out / f"dcgo_comparison_by_aspect_{args.label}.tsv",
        [
            "aspect",
            "ours_significant",
            "theirs_significant",
            "shared",
            "precision_ours_in_theirs",
            "recall_theirs_in_ours",
            "jaccard",
        ],
        aspect_rows,
    )

    # Per-domain precision, so a reader can see whether disagreement is spread
    # evenly or concentrated in a few superfamilies.
    per_domain: Dict[int, Counter] = defaultdict(Counter)
    for pair in ours_pairs:
        per_domain[pair[0]]["ours"] += 1
        if pair in theirs_sig:
            per_domain[pair[0]]["shared"] += 1
    for pair in theirs_sig:
        per_domain[pair[0]]["theirs"] += 1
    write_tsv(
        out / f"dcgo_comparison_by_domain_{args.label}.tsv",
        [
            "sunid",
            "ssf",
            "sccs",
            "description",
            "ours",
            "theirs",
            "shared",
            "precision",
        ],
        sorted(
            (
                [
                    sunid,
                    f"SSF{sunid}",
                    scop.get(sunid, ("", "", ""))[1],
                    scop.get(sunid, ("", "", ""))[2],
                    counts["ours"],
                    counts["theirs"],
                    counts["shared"],
                    fraction(counts["shared"], counts["ours"]),
                ]
                for sunid, counts in per_domain.items()
            ),
            key=lambda row: -row[4],
        ),
    )

    write_tsv(
        out / f"dcgo_comparison_supra_{args.label}.tsv",
        ["metric", "value"],
        [
            ["our_architectures_observed", len(our_all_arch)],
            ["their_architectures", len(their_supra)],
            ["shared_architectures_exact_order", n_arch],
            ["shared_architectures_order_insensitive", n_arch_c],
            ["architecture_precision_exact_order", prec_arch],
            ["architecture_recall_exact_order", rec_arch],
            ["architecture_jaccard_exact_order", jac_arch],
            ["architecture_precision_order_insensitive", prec_arch_c],
            ["architecture_recall_order_insensitive", rec_arch_c],
            ["architecture_jaccard_order_insensitive", jac_arch_c],
            ["architectures_compared", len(common_arch)],
            ["our_associations_compared", len(ours_supra_r)],
            ["their_associations_compared", len(theirs_supra_r)],
            ["shared_associations", n_sp],
            ["precision_ours_in_theirs", prec_sp],
            ["recall_theirs_in_ours", rec_sp],
            ["jaccard", jac_sp],
            ["their_associations_compared_incl_inherited", len(theirs_supra_all_r)],
            ["shared_associations_incl_inherited", n_sp_all],
            ["precision_vs_direct_plus_inherited", prec_sp_all],
            ["recall_vs_direct_plus_inherited", rec_sp_all],
            ["jaccard_vs_direct_plus_inherited", jac_sp_all],
        ],
    )

    logger.info("")
    logger.info("Per-aspect agreement (precision = ours also called by them):")
    for row in aspect_rows:
        logger.info(
            f"  {row[0]:<20} ours={row[1]:>6,}  theirs={row[2]:>6,}  "
            f"shared={row[3]:>6,}  precision={row[4]}  recall={row[5]}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
