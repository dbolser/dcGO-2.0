#!/usr/bin/env python3
"""
dcGO Pipeline - Protein Domain-GO Association Analysis

This script runs the complete dcGO statistical inference pipeline for any species.
It performs domain-GO association analysis using sparse matrix operations and parallel
Fisher's exact tests.

Usage:
    uv run python run_dcgo_human.py [OPTIONS]

Options:
    --species STR            Species to analyze: 'human', 'mouse', etc. (default: human)
    --evidence-filter STR    Evidence code filter: 'all', 'manual', 'experimental' (default: manual)
    --fdr-threshold FLOAT    FDR significance threshold (default: 0.01)
    --num-cores INT          Number of CPU cores for parallel processing (default: 8)
    --output-dir PATH        Output directory for results (default: results/)
    --batch-size INT         Batch size for Fisher tests (default: 50000)
    --enable-true-path       Enable True Path Rule for GO annotation propagation
    --go-ontology PATH       Path to GO ontology file (default: data/raw/go_ontology/go.obo)

Examples:
    # Run for human proteins
    uv run python run_dcgo_human.py --species human --num-cores 16

    # Run for mouse proteins with experimental evidence only
    uv run python run_dcgo_human.py --species mouse --evidence-filter experimental

    # Run with True Path Rule propagation
    uv run python run_dcgo_human.py --enable-true-path --go-ontology data/raw/go_ontology/go.obo
"""

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from loguru import logger
from scipy.stats import hypergeom

logger.remove()
logger.add(sys.stderr, level="INFO")

from src.domain_annotation_parser import DomainAnnotationParser
from src.goa_parser import parse_goa_human
from src.hierarchical_inference import HierarchicalInferenceEngine
from src.ontology_processor import OntologyProcessor
from src.sparse_fisher import (
    DomainMetadata,
    build_sparse_matrices,
    compute_contingency_tables_sparse,
    parse_domain_id,
)
from src.vectorized_fisher import benjamini_hochberg_correction, fisher_exact_parallel


@dataclass
class AssociationResult:
    """Simple dataclass to hold association results for True Path Rule."""

    domain: str
    go_term: str
    p_value: float
    q_value: float
    hyper_score: float
    a: int  # proteins with both
    b: int  # proteins with domain only
    c: int  # proteins with GO only
    d: int  # proteins with neither


def calculate_hypergeometric_score(a: int, b: int, c: int, d: int) -> float:
    """
    Calculate hypergeometric-based association score on 1-100 scale.

    Args:
        a: Proteins with both domain and GO term
        b: Proteins with GO term only (not domain)
        c: Proteins with domain only (not GO)
        d: Proteins with neither

    Returns:
        float: Association score between 1.0 and 100.0
    """
    n = a + b + c + d  # total proteins
    k = a + c  # proteins with domain
    m = a + b  # proteins with GO term
    x = a  # proteins with both

    if k == 0 or m == 0 or x == 0:
        return 0.0

    try:
        # Calculate hypergeometric survival function (1 - CDF)
        # P(X ≥ x) where X ~ Hypergeometric(n, k, m)
        p_hyper = hypergeom.sf(x - 1, n, k, m)

        if p_hyper > 0 and not np.isnan(p_hyper):
            # Convert to -log10 scale
            score = -np.log10(p_hyper)
            # Scale to 1-100 range (typical values 1e-50 to 1e-1 give scores 1-500)
            scaled_score = min(100.0, max(1.0, score * 10))
        else:
            scaled_score = 100.0  # Maximum score for p ≈ 0

        return scaled_score

    except (ValueError, OverflowError, ZeroDivisionError):
        return 50.0  # Neutral score for edge cases


def main():
    parser = argparse.ArgumentParser(
        description="dcGO Pipeline - Human Protein Analysis"
    )
    parser.add_argument(
        "--evidence-filter",
        default="manual",
        choices=["all", "manual", "experimental"],
        help="GO annotation evidence filter (default: manual)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.01,
        help="FDR significance threshold (default: 0.01)",
    )
    parser.add_argument(
        "--num-cores", type=int, default=8, help="Number of CPU cores (default: 8)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Output directory (default: results/)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50000,
        help="Batch size for Fisher tests (default: 50000)",
    )
    parser.add_argument(
        "--species",
        default="human",
        help="Species to analyze: 'human', 'mouse', or specific GOA file name (default: human)",
    )
    parser.add_argument(
        "--enable-true-path",
        action="store_true",
        help="Enable True Path Rule for GO annotation propagation",
    )
    parser.add_argument(
        "--go-ontology",
        type=Path,
        default=Path("data/raw/go_ontology/go.obo"),
        help="Path to GO ontology file (default: data/raw/go_ontology/go.obo)",
    )
    parser.add_argument(
        "--enable-supra-domains",
        action="store_true",
        default=True,
        help="Include supra-domain (multi-domain) combinations in analysis (default: True)",
    )
    parser.add_argument(
        "--disable-supra-domains",
        dest="enable_supra_domains",
        action="store_false",
        help="Disable supra-domain analysis (single domains only)",
    )
    parser.add_argument(
        "--enable-shrinkage",
        action="store_true",
        help="Enable hierarchical shrinkage for supra-domains (empirical Bayes regularization)",
    )
    parser.add_argument(
        "--shrinkage-strength",
        type=float,
        default=0.5,
        help="Shrinkage strength factor 0-1 (default: 0.5). Higher = more regularization",
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info(f"dcGO PIPELINE - {args.species.upper()} PROTEIN ANALYSIS")
    logger.info("=" * 70)
    logger.info("Configuration:")
    logger.info(f"  Species: {args.species}")
    logger.info(f"  Evidence filter: {args.evidence_filter}")
    logger.info(f"  FDR threshold: {args.fdr_threshold}")
    logger.info(f"  CPU cores: {args.num_cores}")
    logger.info(
        f"  Supra-domains: {'ENABLED' if args.enable_supra_domains else 'DISABLED'}"
    )
    if args.enable_supra_domains and args.enable_shrinkage:
        logger.info(
            f"  Hierarchical shrinkage: ENABLED (strength={args.shrinkage_strength})"
        )
    logger.info(f"  Output directory: {args.output_dir}")

    # File paths - support different species
    goa_file = Path(f"data/raw/goa_annotations/goa_{args.species}.gaf.gz")
    interpro_file = Path(f"data/interim/protein2ipr_{args.species}.dat.gz")

    # Check files exist
    if not interpro_file.exists():
        logger.error(f"{args.species.title()} InterPro file not found: {interpro_file}")
        logger.error(f"Please extract {args.species} data from protein2ipr.dat.gz")
        return 1

    if not goa_file.exists():
        logger.error(f"{args.species.title()} GOA file not found: {goa_file}")
        logger.error(f"Please download GOA file for {args.species}")
        return 1

    # Load data
    logger.info("")
    logger.info("STAGE 1: Loading Data")
    logger.info("─" * 70)

    logger.info("Parsing GO annotations...")
    protein_go_map = parse_goa_human(
        goa_file, evidence_filter=args.evidence_filter, aspects={"P", "F", "C"}
    )

    logger.info("Parsing domain annotations...")
    parser_obj = DomainAnnotationParser(max_supra_domain_length=3, min_domain_length=10)
    domain_architectures = parser_obj.parse_protein2ipr_file(interpro_file)

    # Get intersection
    proteins_with_both = set(protein_go_map.keys()) & set(domain_architectures.keys())

    # Build protein-domain map (using lists for compatibility with ontology processor)
    # CRITICAL: Include both single domains AND supra-domains as per dcGO methodology
    protein_domain_map = {}
    all_domains = set()
    single_domain_count = 0
    supra_domain_count = 0

    for protein_id in proteins_with_both:
        arch = domain_architectures[protein_id]

        # Always include single domains
        domains = list(arch.single_domains)
        single_domain_count += len(domains)

        # Include supra-domains if enabled (default: True)
        if args.enable_supra_domains:
            domains.extend(arch.supra_domains)
            supra_domain_count += len(arch.supra_domains)

        if domains:
            protein_domain_map[protein_id] = domains
            all_domains.update(domains)

    logger.info(f"  Single domains: {single_domain_count:,} annotations")
    if args.enable_supra_domains:
        logger.info(f"  Supra-domains: {supra_domain_count:,} annotations")
        logger.info(
            f"  Total domain features: {single_domain_count + supra_domain_count:,}"
        )

    # Get all GO terms
    all_go_terms = set()
    for protein_id in proteins_with_both:
        all_go_terms.update(protein_go_map[protein_id])

    domain_list = sorted(all_domains)
    go_list = sorted(all_go_terms)

    logger.info(
        f"✓ Dataset prepared: {len(proteins_with_both):,} proteins, {len(domain_list):,} domains, {len(go_list):,} GO terms"
    )
    logger.info(f"  Total tests: {len(domain_list) * len(go_list):,}")

    # Build sparse matrices
    logger.info("")
    logger.info("STAGE 2: Building Sparse Matrices")
    logger.info("─" * 70)
    start_time = time.time()

    protein_domain_matrix, protein_go_matrix, domain_metadata = build_sparse_matrices(
        protein_domain_map, protein_go_map, domain_list, go_list
    )

    matrix_time = time.time() - start_time
    logger.info(f"✓ Sparse matrices built in {matrix_time:.2f}s")

    # Compute contingency tables
    logger.info("")
    logger.info("STAGE 3: Computing Contingency Tables")
    logger.info("─" * 70)
    start_time = time.time()

    tables = compute_contingency_tables_sparse(protein_domain_matrix, protein_go_matrix)

    table_time = time.time() - start_time
    logger.info(
        f"✓ Contingency tables computed in {table_time:.2f}s ({table_time / 60:.1f} min)"
    )

    # Run Fisher's exact tests
    logger.info("")
    logger.info("STAGE 4: Running Fisher's Exact Tests")
    logger.info("─" * 70)
    logger.info(f"Processing {len(tables):,} tests with {args.num_cores} cores...")
    start_time = time.time()

    # Progress callback
    def progress_callback(completed, total):
        progress_pct = (completed / total) * 100
        elapsed = time.time() - start_time
        rate = completed / elapsed if elapsed > 0 else 0
        eta = (total - completed) / rate if rate > 0 else 0
        logger.info(
            f"  Progress: {completed:,} / {total:,} ({progress_pct:.1f}%) | {rate:,.0f} tests/s | ETA: {eta / 60:.1f} min"
        )

    odds_ratios, pvalues = fisher_exact_parallel(
        tables,
        alternative="greater",
        n_jobs=args.num_cores,
        batch_size=args.batch_size,
        progress_callback=progress_callback,
    )

    test_time = time.time() - start_time
    logger.info(
        f"✓ Fisher tests completed in {test_time:.2f}s ({test_time / 60:.1f} min)"
    )
    logger.info(f"  Rate: {len(pvalues) / test_time:,.0f} tests/second")

    # STAGE 4.5: Hierarchical Shrinkage (Optional)
    if args.enable_supra_domains and args.enable_shrinkage:
        logger.info("")
        logger.info("STAGE 4.5: Hierarchical Shrinkage")
        logger.info("─" * 70)
        start_time = time.time()

        # Initialize shrinkage engine
        shrinkage_engine = HierarchicalInferenceEngine(
            shrinkage_strength=args.shrinkage_strength,
            min_observations=3,  # From config
        )

        # Apply shrinkage to p-values
        original_pvalues = pvalues.copy()
        pvalues = shrinkage_engine.shrink_pvalues(
            pvalues, domain_list, go_list, domain_metadata
        )

        # Report shrinkage statistics
        stats = shrinkage_engine.get_shrinkage_statistics(
            original_pvalues,
            pvalues,
            domain_list,
            domain_metadata,
            significance_threshold=args.fdr_threshold,
        )

        shrinkage_time = time.time() - start_time
        logger.info(f"✓ Hierarchical shrinkage completed in {shrinkage_time:.2f}s")
        logger.info(f"  Supra-domain tests affected: {stats['n_supra_tests']:,}")
        logger.info(
            f"  P-values increased (regularized): {stats['n_pvalues_increased']:,} ({stats['pct_pvalues_increased']:.1f}%)"
        )
        logger.info(f"  Median p-value ratio: {stats['median_pvalue_ratio']:.3f}")

    # Apply FDR correction
    logger.info("")
    logger.info("STAGE 5: FDR Correction")
    logger.info("─" * 70)
    start_time = time.time()

    adjusted_pvalues, threshold = benjamini_hochberg_correction(
        pvalues, alpha=args.fdr_threshold
    )

    fdr_time = time.time() - start_time
    logger.info(f"✓ FDR correction completed in {fdr_time:.2f}s")
    logger.info(f"  Threshold p-value: {threshold:.2e}")

    # Count significant associations
    significant = adjusted_pvalues <= args.fdr_threshold
    n_significant = int(significant.sum())

    # STAGE 5.5: True Path Rule (Optional)
    propagated_annotations = []
    if args.enable_true_path and n_significant > 0:
        logger.info("")
        logger.info("STAGE 5.5: True Path Rule Propagation")
        logger.info("─" * 70)

        # Check GO ontology file exists
        if not args.go_ontology.exists():
            logger.error(f"GO ontology file not found: {args.go_ontology}")
            logger.error("Skipping True Path Rule propagation")
        else:
            start_time = time.time()

            # Load GO ontology
            logger.info(f"Loading GO ontology from: {args.go_ontology}")
            ontology_processor = OntologyProcessor(args.go_ontology)

            # Create AssociationResult objects for significant associations
            logger.info("Preparing significant associations for True Path filtering...")
            significant_associations = []
            significant_indices = np.where(significant)[0]

            for idx in significant_indices:
                domain_idx = idx // len(go_list)
                go_idx = idx % len(go_list)
                table = tables[idx]
                a, b = int(table[0, 0]), int(table[0, 1])
                c, d = int(table[1, 0]), int(table[1, 1])

                assoc = AssociationResult(
                    domain=domain_list[domain_idx],
                    go_term=go_list[go_idx],
                    p_value=float(pvalues[idx]),
                    q_value=float(adjusted_pvalues[idx]),
                    hyper_score=calculate_hypergeometric_score(a, b, c, d),
                    a=a,
                    b=b,
                    c=c,
                    d=d,
                )
                significant_associations.append(assoc)

            logger.info(
                f"Applying True Path Rule to {len(significant_associations):,} significant associations..."
            )

            # Apply optimal level filtering
            # Note: alpha_threshold is for raw p-values from Fisher tests, not FDR-corrected
            # Using 0.05 as recommended threshold for parent-child comparison tests
            filtered_associations = ontology_processor.apply_optimal_level_filter(
                significant_associations,
                protein_domain_map,
                protein_go_map,
                min_background_size=3,
                alpha_threshold=0.05,
            )

            logger.info(
                f"✓ Optimal level filtering: {len(filtered_associations):,} associations retained"
            )

            # Propagate annotations up the GO hierarchy
            propagated_annotations = ontology_processor.propagate_annotations(
                filtered_associations
            )

            direct_count = sum(
                1 for ann in propagated_annotations if ann.annotation_type == "direct"
            )
            propagated_count = len(propagated_annotations) - direct_count

            logger.info(
                f"✓ Generated {len(propagated_annotations):,} total annotations:"
            )
            logger.info(f"  - Direct: {direct_count:,}")
            logger.info(f"  - Propagated: {propagated_count:,}")

            true_path_time = time.time() - start_time
            logger.info(f"✓ True Path Rule completed in {true_path_time:.2f}s")

    # Export results
    logger.info("")
    logger.info("STAGE 6: Exporting Results")
    logger.info("─" * 70)

    # Calculate hypergeometric scores for significant associations
    logger.info("Calculating hypergeometric scores for significant associations...")
    significant_indices = np.where(significant)[0]

    # Export significant associations with hypergeometric scores and domain types
    output_file = args.output_dir / "domain_go_associations_significant.tsv"
    with open(output_file, "w") as f:
        f.write(
            "domain\tgo_term\tp_value\tadj_p_value\todds_ratio\thyper_score\t"
            "domain_type\tconstituent_domains\tn_observations\n"
        )
        for idx in significant_indices:
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)
            domain_id = domain_list[domain_idx]

            # Get contingency table values for hypergeometric score
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            hyper_score = calculate_hypergeometric_score(a, b, c, d)

            # Get domain metadata
            meta = domain_metadata[domain_id]
            constituents = (
                ",".join(meta.constituent_domains) if meta.constituent_domains else "-"
            )

            f.write(
                f"{domain_id}\t{go_list[go_idx]}\t"
                f"{pvalues[idx]:.6e}\t{adjusted_pvalues[idx]:.6e}\t{odds_ratios[idx]:.4f}\t{hyper_score:.2f}\t"
                f"{meta.domain_type.value}\t{constituents}\t{meta.observation_count}\n"
            )

    logger.info(f"✓ Exported significant associations to: {output_file}")
    logger.info(f"  {n_significant:,} associations (FDR < {args.fdr_threshold})")

    # Export top associations with hypergeometric scores and domain types
    top_file = args.output_dir / "domain_go_associations_top100.tsv"
    top_indices = np.argsort(pvalues)[:100]
    with open(top_file, "w") as f:
        f.write(
            "rank\tdomain\tgo_term\tp_value\tadj_p_value\todds_ratio\thyper_score\t"
            "domain_type\tconstituent_domains\tn_observations\n"
        )
        for rank, idx in enumerate(top_indices, 1):
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)
            domain_id = domain_list[domain_idx]

            # Get contingency table values for hypergeometric score
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            hyper_score = calculate_hypergeometric_score(a, b, c, d)

            # Get domain metadata
            meta = domain_metadata[domain_id]
            constituents = (
                ",".join(meta.constituent_domains) if meta.constituent_domains else "-"
            )

            f.write(
                f"{rank}\t{domain_id}\t{go_list[go_idx]}\t"
                f"{pvalues[idx]:.6e}\t{adjusted_pvalues[idx]:.6e}\t{odds_ratios[idx]:.4f}\t{hyper_score:.2f}\t"
                f"{meta.domain_type.value}\t{constituents}\t{meta.observation_count}\n"
            )

    logger.info(f"✓ Exported top 100 associations to: {top_file}")

    # Export propagated annotations if True Path Rule was applied
    if propagated_annotations:
        annotations_file = args.output_dir / "domain_go_annotations_propagated.tsv"
        with open(annotations_file, "w") as f:
            f.write(
                "domain\tgo_term\tq_value\tassociation_score\tannotation_type\tdirect_source_term\n"
            )
            for ann in propagated_annotations:
                f.write(
                    f"{ann.domain}\t{ann.go_term}\t{ann.q_value:.6e}\t{ann.association_score:.2f}\t"
                    f"{ann.annotation_type}\t{ann.direct_source_term}\n"
                )

        logger.info(
            f"✓ Exported {len(propagated_annotations):,} propagated annotations to: {annotations_file}"
        )

    # Performance summary
    total_time = matrix_time + table_time + test_time + fdr_time

    logger.info("")
    logger.info("=" * 70)
    logger.info("PIPELINE COMPLETE!")
    logger.info("=" * 70)
    logger.info("Results Summary:")
    logger.info(f"  Total domain-GO tests: {len(pvalues):,}")
    logger.info(
        f"  Significant associations (FDR < {args.fdr_threshold}): {n_significant:,} ({n_significant / len(pvalues) * 100:.2f}%)"
    )
    if propagated_annotations:
        direct_count = sum(
            1 for ann in propagated_annotations if ann.annotation_type == "direct"
        )
        propagated_count = len(propagated_annotations) - direct_count
        logger.info(
            f"  True Path Rule annotations: {len(propagated_annotations):,} total ({direct_count:,} direct + {propagated_count:,} propagated)"
        )
    logger.info(f"  Total runtime: {total_time:.1f}s ({total_time / 60:.1f} minutes)")
    logger.info("")
    logger.info("Output files:")
    logger.info(f"  {output_file}")
    logger.info(f"  {top_file}")
    if propagated_annotations:
        logger.info(f"  {annotations_file}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
