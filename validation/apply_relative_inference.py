#!/usr/bin/env python3
"""Apply the *relative (parental-background) inference* to dcGO associations.

VALIDATION_PLAN §2/§3 method-audit follow-up. The original dcGO (Fang & Gough
2013) runs **two** hypergeometric tests per (domain, GO term): an *overall* test
against the whole proteome (what our pipeline already does) **and** a *relative*
test whose background is only the proteins annotated to the term's direct parents
— then keeps the more conservative result. The relative test is what enforces
specificity: it asks whether the domain is enriched for the *child* term beyond
what the *parent* term already explains, so associations that merely ride a
popular parent are demoted. Our pipeline shipped the overall test only; this
script adds the missing relative test.

It runs *only on the already-significant associations* (a subset restriction can
only make a p-value larger, never smaller, so nothing that failed the overall
test could pass), so it is cheap — no re-running of the full Fisher sweep.

For each significant (domain, GO):
  * direct parents = the term's immediate parents in the GO DAG;
  * background N = proteins annotated (true-path) to *all* direct parents;
  * within N: K = domain carriers, n = child-term carriers, k = both;
  * relative p = hypergeom.sf(k - 1, N, K, n)  (enrichment of the child).
Keep the association iff the relative test is still significant (p < ``--alpha``).
Root terms, terms absent from the DAG, and under-powered backgrounds
(N < ``--min-background``) are kept (conservative — untestable, not refuted).

Output: a filtered associations TSV (same columns, plus ``relative_p``) that is a
drop-in ``--predictions`` for ``temporal_benchmark.py``.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from loguru import logger
from scipy.stats import hypergeom

logger.remove()
logger.add(sys.stderr, level="INFO")

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

from src.domain_annotation_parser import DomainAnnotationParser  # noqa: E402
from src.goa_parser import parse_goa_human  # noqa: E402
from src.ontology_processor import OntologyProcessor  # noqa: E402
from validation.association_io import load_associations  # noqa: E402


def build_protein_domain_map(
    interpro_file: Path, enable_supra: bool
) -> dict[str, list]:
    parser = DomainAnnotationParser(max_supra_domain_length=3, min_domain_length=10)
    architectures = parser.parse_protein2ipr_file(interpro_file)
    out: dict[str, list] = {}
    for protein, arch in architectures.items():
        domains = list(arch.single_domains)
        if enable_supra:
            domains.extend(arch.supra_domains)
        if domains:
            out[protein] = domains
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions", type=Path, required=True)
    ap.add_argument("--t0-gaf", type=Path, required=True)
    ap.add_argument(
        "--interpro", type=Path, default=Path("data/interim/protein2ipr_human.dat.gz")
    )
    ap.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-background", type=int, default=10)
    ap.add_argument("--enable-supra-domains", action="store_true", default=True)
    ap.add_argument(
        "--disable-supra-domains", dest="enable_supra_domains", action="store_false"
    )
    args = ap.parse_args()

    logger.info("Loading GO ontology...")
    processor = OntologyProcessor(args.go_ontology)

    logger.info("Parsing t0 GOA (manual evidence) and propagating (true-path)...")
    protein_go = parse_goa_human(args.t0_gaf, evidence_filter="manual")
    # Propagate each protein's annotations to ancestors so parent backgrounds are
    # correct (a protein annotated to a child is annotated to its parents too).
    go2prot: dict[str, set] = defaultdict(set)
    for protein, terms in protein_go.items():
        closed = set(terms)
        for t in terms:
            closed |= processor.get_ancestors(t)
        for t in closed:
            go2prot[t].add(protein)

    logger.info("Parsing t0 domain architectures...")
    protein_domain = build_protein_domain_map(args.interpro, args.enable_supra_domains)
    dom2prot: dict[str, set] = defaultdict(set)
    for protein, domains in protein_domain.items():
        for d in domains:
            dom2prot[d].add(protein)

    logger.info(f"Loading significant associations: {args.predictions}")
    try:
        df = load_associations(args.predictions, required_columns={"p_value"})
    except (OSError, pd.errors.ParserError, ValueError) as exc:
        logger.error(str(exc))
        return 1
    logger.info(f"  {len(df):,} associations to test")

    relative_p = []
    kept_mask = []
    n_root_or_unknown = 0
    n_underpowered = 0
    for domain, go in zip(df["domain"], df["go_term"]):
        if go not in processor.go_graph:
            relative_p.append(float("nan"))
            kept_mask.append(True)
            n_root_or_unknown += 1
            continue
        parents = list(processor.go_graph.predecessors(go))
        if not parents:
            relative_p.append(float("nan"))
            kept_mask.append(True)
            n_root_or_unknown += 1
            continue
        # Background: proteins annotated to ALL direct parents.
        background = set(go2prot.get(parents[0], set()))
        for p in parents[1:]:
            background &= go2prot.get(p, set())
        big_n = len(background)
        if big_n < args.min_background:
            relative_p.append(float("nan"))
            kept_mask.append(True)
            n_underpowered += 1
            continue
        dom_p = dom2prot.get(domain, set())
        child_p = go2prot.get(go, set())
        big_k = len(background & dom_p)  # domain carriers in background
        small_n = len(background & child_p)  # child-annotated in background
        small_k = len(background & dom_p & child_p)  # both
        if big_k == 0 or small_n == 0:
            relative_p.append(1.0)
            kept_mask.append(False)
            continue
        # P(overlap >= small_k) under hypergeometric(N, K, n).
        p = float(hypergeom.sf(small_k - 1, big_n, big_k, small_n))
        relative_p.append(p)
        kept_mask.append(p < args.alpha)

    # Match the paper: the reported p-value is the *more conservative* of the
    # overall and relative tests, so downstream scoring (which ranks on p_value)
    # reflects both. Untested rows (root/absent/under-powered) keep the overall p.
    df["overall_p"] = df["p_value"]
    df["relative_p"] = relative_p
    tested = df["relative_p"].notna()
    df.loc[tested, "p_value"] = df.loc[tested, ["overall_p", "relative_p"]].max(axis=1)
    kept = df[kept_mask].copy()
    logger.info(
        f"Relative inference: kept {len(kept):,} / {len(df):,} "
        f"({100 * len(kept) / len(df):.1f}%) at alpha={args.alpha}"
    )
    logger.info(
        f"  ({n_root_or_unknown:,} root/absent kept untested; "
        f"{n_underpowered:,} under-powered backgrounds kept untested)"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    kept.to_csv(args.output, sep="\t", index=False)
    logger.info(f"✓ Saved filtered associations: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
