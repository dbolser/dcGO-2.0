#!/usr/bin/env python3
"""
Validation Script - Compare dcGO predictions against the InterPro2GO reference.

InterPro2GO is a *manually curated, deliberately incomplete, positive-only*
mapping. It is therefore a reference for **recall / coverage**, not a source of
truth for precision: a predicted (domain, GO) pair that is absent from
InterPro2GO is a *candidate*, not a demonstrated false positive.

This script reflects that. It:
  * propagates BOTH predictions and the reference up the GO DAG to their
    ancestor closure before comparing (True-Path-aware matching), so a prediction
    of a specific term is credited against a more general curated term and vice
    versa;
  * restricts the comparison to the domains present in BOTH sets (a domain absent
    from InterPro2GO can only add noise);
  * reports **reference coverage (recall)** as the headline metric and labels the
    remainder as *candidate* predictions (not "false positives");
  * de-duplicates the threshold sweep.

See VALIDATION_PLAN.md §1.
"""

import sys
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd
from loguru import logger

logger.remove()
logger.add(sys.stderr, level="INFO")

Pair = tuple[str, str]  # (domain_id, go_term)


# --------------------------------------------------------------------------- #
# Pure, unit-testable core                                                     #
# --------------------------------------------------------------------------- #
def propagate_pairs(
    pairs: Iterable[Pair], get_ancestors: Callable[[str], set[str]]
) -> set[Pair]:
    """Expand each (domain, go) to include (domain, ancestor) for every ancestor.

    The original term is always retained. ``get_ancestors`` returns the set of
    more-general ancestor terms for a GO id (excluding the term itself); a
    no-propagation run can pass ``lambda _go: set()``.
    """
    out: set[Pair] = set()
    for domain, go in pairs:
        out.add((domain, go))
        for ancestor in get_ancestors(go):
            out.add((domain, ancestor))
    return out


def restrict_to_shared_domains(
    pred_set: set[Pair], ref_set: set[Pair]
) -> tuple[set[Pair], set[Pair], set[str]]:
    """Keep only pairs whose domain appears in both the prediction and reference.

    A domain that InterPro2GO does not cover at all cannot inform precision, so
    including it only depresses the metric artificially.
    """
    shared = {d for d, _ in pred_set} & {d for d, _ in ref_set}
    pred = {(d, g) for d, g in pred_set if d in shared}
    ref = {(d, g) for d, g in ref_set if d in shared}
    return pred, ref, shared


def compute_metrics(pred_set: set[Pair], ref_set: set[Pair], name: str) -> dict:
    """Coverage-first metrics for one prediction set against the reference.

    ``pred_set`` and ``ref_set`` are expected to be already propagated and
    restricted to the shared domain space.
    """
    recovered = pred_set & ref_set
    candidates = pred_set - ref_set  # NOT false positives — curation gaps
    reference_coverage = len(recovered) / len(ref_set) if ref_set else 0.0
    # A precision *lower bound*: fraction of predictions confirmed by curation,
    # on the shared domain space. Real precision is higher (curation is partial).
    precision_lower_bound = len(recovered) / len(pred_set) if pred_set else 0.0
    return {
        "threshold": name,
        "n_predictions": len(pred_set),
        "n_reference": len(ref_set),
        "recovered": len(recovered),
        "reference_coverage": reference_coverage,
        "candidate_predictions": len(candidates),
        "precision_lower_bound": precision_lower_bound,
    }


# --------------------------------------------------------------------------- #
# I/O                                                                          #
# --------------------------------------------------------------------------- #
def parse_interpro2go(interpro2go_file: Path) -> set[Pair]:
    """Parse the InterPro2GO reference into a set of (domain_id, go_id) pairs.

    Line format:
        InterPro:IPR000003 Retinoid X receptor > GO:DNA binding ; GO:0003677
    """
    logger.info(f"Parsing InterPro2GO reference: {interpro2go_file}")
    reference: set[Pair] = set()

    with open(interpro2go_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split(">")
            if len(parts) != 2:
                continue
            try:
                interpro_id = parts[0].split()[0].replace("InterPro:", "")
                _go_name, go_id = parts[1].rsplit(";", 1)
                reference.add((interpro_id, go_id.strip()))
            except (ValueError, IndexError):
                logger.debug(f"Skipping malformed line: {line.strip()}")
                continue

    domains = {d for d, _ in reference}
    go_terms = {g for _, g in reference}
    logger.info(f"✓ Loaded {len(reference):,} reference pairs")
    logger.info(
        f"  Unique domains: {len(domains):,}; unique GO terms: {len(go_terms):,}"
    )
    return reference


def load_dcgo_predictions(predictions_file: Path) -> pd.DataFrame:
    """Load dcGO predicted associations."""
    logger.info(f"Loading dcGO predictions: {predictions_file}")
    df = pd.read_csv(predictions_file, sep="\t")
    required = {"domain", "go_term", "p_value", "adj_p_value", "hyper_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Predictions file is missing required columns: {missing}")
    logger.info(f"✓ Loaded {len(df):,} predicted associations")
    logger.info(
        f"  Unique domains: {df['domain'].nunique():,}; "
        f"unique GO terms: {df['go_term'].nunique():,}"
    )
    return df


def build_get_ancestors(obo_file: Path | None) -> Callable[[str], set[str]]:
    """Return a get_ancestors(go) function, propagating via the GO DAG if available."""
    if obo_file is None or not obo_file.exists():
        logger.warning(
            "No GO ontology available — comparing without propagation. "
            "Pass a go-basic.obo to enable True-Path-aware matching."
        )
        return lambda _go: set()

    # Imported lazily so the metric helpers don't require obonet/networkx.
    project_root = str(Path(__file__).resolve().parents[1])
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    import networkx as nx

    from src.ontology_processor import OntologyProcessor

    logger.info(f"Loading GO ontology for propagation: {obo_file}")
    processor = OntologyProcessor(obo_file)
    graph = processor.go_graph

    # obonet builds edges child -> parent (verified empirically), so the more
    # GENERAL ancestor terms of a node are reachable *forward* from it:
    # nx.descendants(graph, term). NOTE: do NOT use OntologyProcessor.get_ancestors
    # here — with this edge direction it returns descendants (see issue tracker,
    # the True-Path direction bug). Memoize since we propagate repeatedly.
    cache: dict[str, set[str]] = {}

    def general_ancestors(go: str) -> set[str]:
        if go not in cache:
            cache[go] = nx.descendants(graph, go) if go in graph else set()
        return cache[go]

    return general_ancestors


# --------------------------------------------------------------------------- #
# Sweep                                                                        #
# --------------------------------------------------------------------------- #
def _pairs_from_df(df: pd.DataFrame) -> set[Pair]:
    return set(zip(df["domain"], df["go_term"]))


def calculate_overlap_at_thresholds(
    predictions: pd.DataFrame,
    reference: set[Pair],
    thresholds: dict,
    get_ancestors: Callable[[str], set[str]],
) -> pd.DataFrame:
    """Coverage/precision-lower-bound across a de-duplicated threshold sweep.

    The reference is propagated once. The shared domain space is fixed by the
    *full* prediction set so the coverage denominator is stable across thresholds
    (a domain we never predict is outside the testable universe).
    """
    logger.info("Calculating coverage at various thresholds (propagated)...")

    # Propagation only adds ancestor GO terms; it never changes the domain set.
    # So compute the shared domain space from the RAW domains and filter to it
    # *before* propagating — propagating non-shared domains would be wasted work.
    shared = set(predictions["domain"]) & {d for d, _ in reference}
    ref_fixed = propagate_pairs(
        {pair for pair in reference if pair[0] in shared}, get_ancestors
    )
    logger.info(
        f"  Shared domains: {len(shared):,}; reference pairs on shared domains "
        f"(propagated): {len(ref_fixed):,}"
    )

    results = []

    def sweep(column: str, values: Iterable[float], op: str, label: str):
        for v in sorted(set(values)):  # de-dup + stable order
            if op == "le":
                filtered = predictions[predictions[column] <= v]
                name = f"{label}≤{v:.2e}"
            else:
                filtered = predictions[predictions[column] >= v]
                name = f"{label}≥{v:g}"
            pred_pairs = {p for p in _pairs_from_df(filtered) if p[0] in shared}
            pred_shared = propagate_pairs(pred_pairs, get_ancestors)
            results.append(compute_metrics(pred_shared, ref_fixed, name))

    sweep("p_value", thresholds["p_value"], "le", "p_value")
    sweep("adj_p_value", thresholds["adj_p_value"], "le", "adj_p")
    sweep("hyper_score", thresholds["hyper_score"], "ge", "score")

    return pd.DataFrame(results)


def analyze_candidate_predictions(
    predictions: pd.DataFrame,
    reference: set[Pair],
    get_ancestors: Callable[[str], set[str]],
    top_n: int = 100,
) -> pd.DataFrame:
    """Top candidate predictions: high-confidence pairs not in the (propagated) reference.

    These are curation-gap candidates, not errors. Only pairs on the shared
    domain space are considered (a domain absent from the reference can't be
    called 'novel' relative to it).
    """
    logger.info(f"Selecting top {top_n} candidate predictions (not in reference)...")

    # Filter the reference to domains we actually predict before propagating.
    pred_domains = set(predictions["domain"])
    ref_lookup = propagate_pairs(
        {pair for pair in reference if pair[0] in pred_domains}, get_ancestors
    )
    ref_domains = {d for d, _ in ref_lookup}

    def is_candidate(d: str, g: str) -> bool:
        if d not in ref_domains:
            return False  # outside shared domain space
        if (d, g) in ref_lookup:
            return False
        # Also not covered by propagating this prediction up to a curated ancestor
        return not any((d, anc) in ref_lookup for anc in get_ancestors(g))

    # Vectorized over columns — pandas apply(axis=1) is far slower here.
    mask = [
        is_candidate(d, g)
        for d, g in zip(predictions["domain"], predictions["go_term"])
    ]
    candidates = (
        predictions[mask].sort_values("hyper_score", ascending=False).head(top_n).copy()
    )
    logger.info(f"✓ {sum(mask):,} candidate pairs on shared domains")
    return candidates


def create_visualizations(metrics_df: pd.DataFrame, output_dir: Path) -> None:
    """Coverage-focused validation plots."""
    logger.info("Creating visualizations...")
    # Lazy import: metric computation must not require a plotting backend.
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    ax.scatter(
        metrics_df["n_predictions"],
        metrics_df["reference_coverage"],
        c=metrics_df["precision_lower_bound"],
        cmap="viridis",
        s=90,
        alpha=0.8,
    )
    ax.set_xscale("log")
    ax.set_xlabel("Predictions (propagated, shared domains)")
    ax.set_ylabel("Reference coverage (recall)")
    ax.set_title("Coverage of InterPro2GO vs prediction count", fontweight="bold")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.scatter(
        metrics_df["reference_coverage"],
        metrics_df["precision_lower_bound"],
        s=90,
        alpha=0.8,
    )
    ax.set_xlabel("Reference coverage (recall)")
    ax.set_ylabel("Precision lower bound")
    ax.set_title(
        "Coverage vs precision lower bound\n(true precision is higher — curation is partial)",
        fontweight="bold",
    )
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = output_dir / "validation_metrics.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"✓ Saved visualization: {out}")


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def main() -> None:
    logger.info("=" * 70)
    logger.info("dcGO VALIDATION — InterPro2GO coverage (propagated, reframed)")
    logger.info("=" * 70)

    reference_file = Path("data/raw/interpro2go/interpro2go")
    # Back-compat: older layout kept it under interpro_mappings/
    if not reference_file.exists():
        alt = Path("data/raw/interpro_mappings/interpro2go")
        if alt.exists():
            reference_file = alt
    predictions_file = Path("results/domain_go_associations_significant.tsv")
    obo_file = Path("data/raw/go_ontology/go-basic.obo")
    output_dir = Path("validation")
    output_dir.mkdir(exist_ok=True)

    # The GO ontology is REQUIRED here: the whole point of §1 is a propagated
    # comparison. Without it build_get_ancestors would silently fall back to
    # no-propagation and reproduce the unpropagated comparison this fixes.
    missing_inputs = [
        str(p) for p in (reference_file, predictions_file, obo_file) if not p.exists()
    ]
    if missing_inputs:
        logger.error(
            "Missing required inputs: " + ", ".join(missing_inputs) + ". Download "
            "the reference and ontology with: uv run python scripts/download_data.py "
            "--datasets interpro2go --datasets go_ontology"
        )
        sys.exit(1)

    reference = parse_interpro2go(reference_file)
    predictions = load_dcgo_predictions(predictions_file)
    get_ancestors = build_get_ancestors(obo_file)

    thresholds = {
        "p_value": [1e-10, 1e-8, 1e-6, 1e-4, 1e-2],
        "adj_p_value": [1e-6, 1e-4, 1e-2, 0.05, 0.1],
        "hyper_score": [90, 80, 70, 60, 50, 40, 30, 20],
    }
    metrics_df = calculate_overlap_at_thresholds(
        predictions, reference, thresholds, get_ancestors
    )
    metrics_df.to_csv(output_dir / "performance_metrics.tsv", sep="\t", index=False)
    logger.info(f"✓ Saved metrics: {output_dir / 'performance_metrics.tsv'}")

    best = metrics_df.loc[metrics_df["reference_coverage"].idxmax()]
    logger.info(
        f"✓ Best reference coverage: {best['reference_coverage']:.3f} "
        f"at {best['threshold']} (recovered {best['recovered']:,.0f} pairs)"
    )

    candidates = analyze_candidate_predictions(predictions, reference, get_ancestors)
    candidates.to_csv(
        output_dir / "candidate_predictions_top100.tsv", sep="\t", index=False
    )
    logger.info(
        f"✓ Saved candidates: {output_dir / 'candidate_predictions_top100.tsv'}"
    )

    create_visualizations(metrics_df, output_dir)

    logger.info("=" * 70)
    logger.info("VALIDATION COMPLETE")
    logger.info(
        "Reminder: 'candidate_predictions' are curation-gap candidates, not "
        "false positives. Coverage is the defensible metric against InterPro2GO."
    )


if __name__ == "__main__":
    main()
