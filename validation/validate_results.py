#!/usr/bin/env python3
"""
Validation Script - Compare dcGO predictions against InterPro gold standard

This script validates dcGO domain-GO associations against manually curated
InterPro2GO mappings, analyzing:
- Precision/Recall at various thresholds
- Overlap statistics
- Novel predictions
- Single vs multi-domain protein contributions
"""

from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from loguru import logger

logger.remove()
logger.add(sys.stderr, level="INFO")


def parse_interpro2go(interpro2go_file: Path) -> dict:
    """
    Parse InterPro2GO gold standard file.

    Returns:
        dict: {(domain_id, go_term): go_name}
    """
    logger.info(f"Parsing InterPro gold standard: {interpro2go_file}")

    gold_standard = {}

    with open(interpro2go_file, "r") as f:
        for line in f:
            # Skip comments
            if line.startswith("!"):
                continue

            # Format: InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677
            try:
                parts = line.strip().split(">")
                if len(parts) != 2:
                    continue

                # Extract InterPro ID
                interpro_part = parts[0].strip()
                interpro_id = interpro_part.split()[0].replace("InterPro:", "")

                # Extract GO term
                go_part = parts[1].strip()
                go_name, go_id = go_part.split(";")
                go_id = go_id.strip()

                gold_standard[(interpro_id, go_id)] = go_name.strip()

            except Exception:
                logger.debug(f"Skipping malformed line: {line.strip()}")
                continue

    logger.info(f"✓ Loaded {len(gold_standard):,} gold standard domain-GO associations")

    # Count unique domains and GO terms
    domains = set(d for d, g in gold_standard.keys())
    go_terms = set(g for d, g in gold_standard.keys())
    logger.info(f"  Unique domains: {len(domains):,}")
    logger.info(f"  Unique GO terms: {len(go_terms):,}")

    return gold_standard


def load_dcgo_predictions(predictions_file: Path) -> pd.DataFrame:
    """Load dcGO predicted associations."""
    logger.info(f"Loading dcGO predictions: {predictions_file}")

    df = pd.read_csv(predictions_file, sep="\t")

    logger.info(f"✓ Loaded {len(df):,} predicted associations")
    logger.info(f"  Unique domains: {df['domain'].nunique():,}")
    logger.info(f"  Unique GO terms: {df['go_term'].nunique():,}")
    logger.info(
        f"  p-value range: {df['p_value'].min():.2e} to {df['p_value'].max():.2e}"
    )
    logger.info(
        f"  Hyper score range: {df['hyper_score'].min():.2f} to {df['hyper_score'].max():.2f}"
    )

    return df


def load_domain_architectures(protein2ipr_file: Path) -> tuple:
    """
    Load protein domain architectures to identify single vs multi-domain proteins.

    Returns:
        (protein_domains, domain_protein_counts)
    """
    logger.info(f"Loading protein domain architectures: {protein2ipr_file}")

    from src.domain_annotation_parser import DomainAnnotationParser

    parser = DomainAnnotationParser(max_supra_domain_length=1, min_domain_length=10)
    architectures = parser.parse_protein2ipr_file(protein2ipr_file)

    # Build domain -> protein count mapping (only single domains, not supra-domains)
    domain_proteins = defaultdict(set)

    for protein_id, arch in architectures.items():
        for domain in arch.single_domains:
            domain_proteins[domain].add(protein_id)

    # Count single-domain proteins per domain
    domain_stats = {}
    for domain, proteins in domain_proteins.items():
        single_domain_proteins = []
        for protein in proteins:
            if len(architectures[protein].single_domains) == 1:
                single_domain_proteins.append(protein)

        domain_stats[domain] = {
            "total_proteins": len(proteins),
            "single_domain_proteins": len(single_domain_proteins),
            "multi_domain_proteins": len(proteins) - len(single_domain_proteins),
            "fraction_single": len(single_domain_proteins) / len(proteins)
            if proteins
            else 0,
        }

    logger.info(f"✓ Analyzed {len(domain_stats):,} domain architectures")

    return domain_proteins, domain_stats


def calculate_overlap_at_thresholds(
    predictions: pd.DataFrame, gold_standard: dict, thresholds: dict
) -> pd.DataFrame:
    """
    Calculate precision, recall, F1 at various thresholds.

    Args:
        predictions: dcGO predictions DataFrame
        gold_standard: Gold standard dict {(domain, go): name}
        thresholds: dict of threshold lists for different metrics

    Returns:
        DataFrame with performance metrics
    """
    logger.info("Calculating performance metrics at various thresholds...")

    results = []

    # Test p-value thresholds
    for pval in thresholds["p_value"]:
        filtered = predictions[predictions["p_value"] <= pval]
        metrics = calculate_metrics(filtered, gold_standard, f"p_value≤{pval:.2e}")
        results.append(metrics)

    # Test adjusted p-value thresholds
    for adj_pval in thresholds["adj_p_value"]:
        filtered = predictions[predictions["adj_p_value"] <= adj_pval]
        metrics = calculate_metrics(filtered, gold_standard, f"adj_p≤{adj_pval:.2e}")
        results.append(metrics)

    # Test hyper score thresholds
    for score in thresholds["hyper_score"]:
        filtered = predictions[predictions["hyper_score"] >= score]
        metrics = calculate_metrics(filtered, gold_standard, f"score≥{score}")
        results.append(metrics)

    df = pd.DataFrame(results)
    return df


def calculate_metrics(
    predictions: pd.DataFrame, gold_standard: dict, threshold_name: str
) -> dict:
    """Calculate precision, recall, F1 for a set of predictions."""

    # Create prediction set
    pred_set = set(zip(predictions["domain"], predictions["go_term"]))
    gold_set = set(gold_standard.keys())

    # Calculate overlap
    true_positives = pred_set & gold_set
    false_positives = pred_set - gold_set
    false_negatives = gold_set - pred_set

    # Calculate metrics
    precision = len(true_positives) / len(pred_set) if pred_set else 0
    recall = len(true_positives) / len(gold_set) if gold_set else 0
    f1 = (
        2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    )

    return {
        "threshold": threshold_name,
        "n_predictions": len(pred_set),
        "n_gold_standard": len(gold_set),
        "true_positives": len(true_positives),
        "false_positives": len(false_positives),
        "false_negatives": len(false_negatives),
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
    }


def analyze_novel_predictions(
    predictions: pd.DataFrame, gold_standard: dict, top_n: int = 100
) -> pd.DataFrame:
    """Identify novel high-confidence predictions not in gold standard."""

    logger.info(f"Analyzing top {top_n} novel predictions...")

    # Filter to predictions not in gold standard
    pred_set = set(zip(predictions["domain"], predictions["go_term"]))
    gold_set = set(gold_standard.keys())
    novel_pairs = pred_set - gold_set

    # Get top novel predictions by hyper score
    novel_df = predictions[
        predictions.apply(
            lambda row: (row["domain"], row["go_term"]) in novel_pairs, axis=1
        )
    ].copy()

    novel_df = novel_df.sort_values("hyper_score", ascending=False).head(top_n)

    logger.info(f"✓ Found {len(novel_pairs):,} novel predictions")
    logger.info(f"  Top novel prediction score: {novel_df['hyper_score'].max():.2f}")

    return novel_df


def analyze_domain_architecture_contribution(
    predictions: pd.DataFrame, gold_standard: dict, domain_stats: dict
) -> pd.DataFrame:
    """Analyze how single vs multi-domain proteins contribute to predictions."""

    logger.info("Analyzing domain architecture contributions...")

    results = []

    # Get overlapping predictions
    pred_set = set(zip(predictions["domain"], predictions["go_term"]))
    gold_set = set(gold_standard.keys())
    overlapping = pred_set & gold_set

    for domain, go_term in overlapping:
        if domain in domain_stats:
            stats = domain_stats[domain]
            results.append(
                {
                    "domain": domain,
                    "go_term": go_term,
                    "total_proteins": stats["total_proteins"],
                    "single_domain_proteins": stats["single_domain_proteins"],
                    "multi_domain_proteins": stats["multi_domain_proteins"],
                    "fraction_single": stats["fraction_single"],
                }
            )

    df = pd.DataFrame(results)

    if not df.empty:
        logger.info(f"✓ Analyzed {len(df)} overlapping associations")
        logger.info(f"  Avg fraction single-domain: {df['fraction_single'].mean():.2%}")
        logger.info(
            f"  Median proteins per association: {df['total_proteins'].median():.0f}"
        )

    return df


def create_visualizations(
    metrics_df: pd.DataFrame, architecture_df: pd.DataFrame, output_dir: Path
):
    """Create validation visualizations."""

    logger.info("Creating visualizations...")

    sns.set_style("whitegrid")

    # 1. Precision-Recall curve
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Precision vs Recall
    ax = axes[0, 0]
    ax.scatter(metrics_df["recall"], metrics_df["precision"], alpha=0.6, s=100)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("Recall", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_title(
        "Precision vs Recall at Different Thresholds", fontsize=14, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)

    # F1 Score by number of predictions
    ax = axes[0, 1]
    ax.scatter(
        metrics_df["n_predictions"],
        metrics_df["f1_score"],
        alpha=0.6,
        s=100,
        c=metrics_df["precision"],
        cmap="viridis",
    )
    ax.set_xlabel("Number of Predictions", fontsize=12)
    ax.set_ylabel("F1 Score", fontsize=12)
    ax.set_title("F1 Score vs Prediction Count", fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    ax.grid(True, alpha=0.3)

    # True Positives vs False Positives
    ax = axes[1, 0]
    ax.scatter(
        metrics_df["false_positives"], metrics_df["true_positives"], alpha=0.6, s=100
    )
    ax.set_xlabel("False Positives", fontsize=12)
    ax.set_ylabel("True Positives", fontsize=12)
    ax.set_title("True Positives vs False Positives", fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    # Architecture contribution
    if not architecture_df.empty:
        ax = axes[1, 1]
        ax.hist(
            architecture_df["fraction_single"], bins=20, alpha=0.7, edgecolor="black"
        )
        ax.axvline(
            architecture_df["fraction_single"].median(),
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Median: {architecture_df['fraction_single'].median():.2%}",
        )
        ax.set_xlabel("Fraction Single-Domain Proteins", fontsize=12)
        ax.set_ylabel("Number of Associations", fontsize=12)
        ax.set_title(
            "Single-Domain Protein Contribution", fontsize=14, fontweight="bold"
        )
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(output_dir / "validation_metrics.png", dpi=300, bbox_inches="tight")
    logger.info(f"✓ Saved visualization: {output_dir / 'validation_metrics.png'}")
    plt.close()


def main():
    """Run validation analysis."""

    logger.info("=" * 70)
    logger.info("dcGO VALIDATION - InterPro Gold Standard Comparison")
    logger.info("=" * 70)

    # Setup paths
    gold_standard_file = Path("data/raw/interpro_mappings/interpro2go")
    predictions_file = Path(
        "results_no_true_path/domain_go_associations_significant.tsv"
    )
    protein2ipr_file = Path("data/interim/protein2ipr_human.dat.gz")
    output_dir = Path("validation")
    output_dir.mkdir(exist_ok=True)

    # 1. Load data
    logger.info("")
    logger.info("STEP 1: Loading Data")
    logger.info("─" * 70)

    gold_standard = parse_interpro2go(gold_standard_file)
    predictions = load_dcgo_predictions(predictions_file)

    # 2. Calculate basic overlap
    logger.info("")
    logger.info("STEP 2: Basic Overlap Analysis")
    logger.info("─" * 70)

    pred_set = set(zip(predictions["domain"], predictions["go_term"]))
    gold_set = set(gold_standard.keys())

    overlap = pred_set & gold_set
    novel = pred_set - gold_set
    missed = gold_set - pred_set

    logger.info(f"Predictions: {len(pred_set):,}")
    logger.info(f"Gold standard: {len(gold_set):,}")
    logger.info(
        f"✓ Overlap: {len(overlap):,} ({len(overlap) / len(pred_set) * 100:.2f}% of predictions)"
    )
    logger.info(
        f"  Novel predictions: {len(novel):,} ({len(novel) / len(pred_set) * 100:.2f}%)"
    )
    logger.info(
        f"  Missed gold standard: {len(missed):,} ({len(missed) / len(gold_set) * 100:.2f}%)"
    )

    # 3. Performance at thresholds
    logger.info("")
    logger.info("STEP 3: Performance at Various Thresholds")
    logger.info("─" * 70)

    thresholds = {
        "p_value": [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.05],
        "adj_p_value": [1e-6, 1e-4, 1e-2, 0.01, 0.05, 0.1],
        "hyper_score": [90, 80, 70, 60, 50, 40, 30, 20],
    }

    metrics_df = calculate_overlap_at_thresholds(predictions, gold_standard, thresholds)

    # Find best F1 score
    best_row = metrics_df.loc[metrics_df["f1_score"].idxmax()]
    logger.info(
        f"✓ Best F1 score: {best_row['f1_score']:.4f} at {best_row['threshold']}"
    )
    logger.info(f"  Precision: {best_row['precision']:.4f}")
    logger.info(f"  Recall: {best_row['recall']:.4f}")
    logger.info(f"  True positives: {best_row['true_positives']:,.0f}")

    # Save metrics
    metrics_df.to_csv(output_dir / "performance_metrics.tsv", sep="\t", index=False)
    logger.info(f"✓ Saved metrics: {output_dir / 'performance_metrics.tsv'}")

    # 4. Novel predictions
    logger.info("")
    logger.info("STEP 4: Novel High-Confidence Predictions")
    logger.info("─" * 70)

    novel_df = analyze_novel_predictions(predictions, gold_standard, top_n=100)
    novel_df.to_csv(output_dir / "novel_predictions_top100.tsv", sep="\t", index=False)
    logger.info(
        f"✓ Saved novel predictions: {output_dir / 'novel_predictions_top100.tsv'}"
    )

    # 5. Domain architecture analysis
    if protein2ipr_file.exists():
        logger.info("")
        logger.info("STEP 5: Domain Architecture Analysis")
        logger.info("─" * 70)

        domain_proteins, domain_stats = load_domain_architectures(protein2ipr_file)
        architecture_df = analyze_domain_architecture_contribution(
            predictions, gold_standard, domain_stats
        )

        if not architecture_df.empty:
            architecture_df.to_csv(
                output_dir / "architecture_contribution.tsv", sep="\t", index=False
            )
            logger.info(
                f"✓ Saved architecture analysis: {output_dir / 'architecture_contribution.tsv'}"
            )
    else:
        logger.warning(f"Protein2IPR file not found: {protein2ipr_file}")
        logger.warning("Skipping domain architecture analysis")
        architecture_df = pd.DataFrame()

    # 6. Visualizations
    logger.info("")
    logger.info("STEP 6: Creating Visualizations")
    logger.info("─" * 70)

    create_visualizations(metrics_df, architecture_df, output_dir)

    # 7. Summary report
    logger.info("")
    logger.info("=" * 70)
    logger.info("VALIDATION COMPLETE!")
    logger.info("=" * 70)
    logger.info("Summary:")
    logger.info(
        f"  Overlap with gold standard: {len(overlap):,} / {len(pred_set):,} ({len(overlap) / len(pred_set) * 100:.2f}%)"
    )
    logger.info(f"  Best F1 score: {best_row['f1_score']:.4f}")
    logger.info(f"  Novel high-confidence predictions: {len(novel):,}")
    if not architecture_df.empty:
        logger.info(
            f"  Avg single-domain contribution: {architecture_df['fraction_single'].mean():.2%}"
        )
    logger.info("")
    logger.info("Output files:")
    logger.info(f"  {output_dir / 'performance_metrics.tsv'}")
    logger.info(f"  {output_dir / 'novel_predictions_top100.tsv'}")
    if not architecture_df.empty:
        logger.info(f"  {output_dir / 'architecture_contribution.tsv'}")
    logger.info(f"  {output_dir / 'validation_metrics.png'}")


if __name__ == "__main__":
    main()
