#!/usr/bin/env python3
"""
dcGO Pipeline - Human Protein Analysis

This script runs the complete dcGO statistical inference pipeline for human proteins.
It performs domain-GO association analysis using sparse matrix operations and parallel
Fisher's exact tests.

Usage:
    uv run python run_dcgo_human.py [OPTIONS]

Options:
    --evidence-filter STR    Evidence code filter: 'all', 'manual', 'experimental' (default: manual)
    --fdr-threshold FLOAT    FDR significance threshold (default: 0.01)
    --num-cores INT          Number of CPU cores for parallel processing (default: 8)
    --output-dir PATH        Output directory for results (default: results/)
    --batch-size INT         Batch size for Fisher tests (default: 50000)
"""

import time
import argparse
from pathlib import Path
from loguru import logger
import sys
import numpy as np
from scipy.stats import hypergeom

logger.remove()
logger.add(sys.stderr, level="INFO")

from src.sparse_fisher import build_sparse_matrices, compute_contingency_tables_sparse
from src.vectorized_fisher import fisher_exact_parallel, benjamini_hochberg_correction
from src.domain_annotation_parser import DomainAnnotationParser
from src.goa_parser import parse_goa_human


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
    k = a + c          # proteins with domain
    m = a + b          # proteins with GO term
    x = a              # proteins with both

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
    parser = argparse.ArgumentParser(description='dcGO Pipeline - Human Protein Analysis')
    parser.add_argument('--evidence-filter', default='manual',
                       choices=['all', 'manual', 'experimental'],
                       help='GO annotation evidence filter (default: manual)')
    parser.add_argument('--fdr-threshold', type=float, default=0.01,
                       help='FDR significance threshold (default: 0.01)')
    parser.add_argument('--num-cores', type=int, default=8,
                       help='Number of CPU cores (default: 8)')
    parser.add_argument('--output-dir', type=Path, default=Path('results'),
                       help='Output directory (default: results/)')
    parser.add_argument('--batch-size', type=int, default=50000,
                       help='Batch size for Fisher tests (default: 50000)')

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("dcGO PIPELINE - HUMAN PROTEIN ANALYSIS")
    logger.info("=" * 70)
    logger.info("Configuration:")
    logger.info(f"  Evidence filter: {args.evidence_filter}")
    logger.info(f"  FDR threshold: {args.fdr_threshold}")
    logger.info(f"  CPU cores: {args.num_cores}")
    logger.info(f"  Output directory: {args.output_dir}")

    # File paths
    goa_file = Path("data/raw/goa_annotations/goa_human.gaf.gz")
    interpro_file = Path("data/interim/protein2ipr_human.dat.gz")

    # Check files exist
    if not interpro_file.exists():
        logger.error(f"Human InterPro file not found: {interpro_file}")
        logger.error("Please run: uv run python extract_human_interpro.py")
        return 1

    # Load data
    logger.info("")
    logger.info("STAGE 1: Loading Data")
    logger.info("─" * 70)

    logger.info("Parsing GO annotations...")
    protein_go_map = parse_goa_human(goa_file, evidence_filter=args.evidence_filter, aspects={'P', 'F', 'C'})

    logger.info("Parsing domain annotations...")
    parser_obj = DomainAnnotationParser(max_supra_domain_length=3, min_domain_length=10)
    domain_architectures = parser_obj.parse_protein2ipr_file(interpro_file)

    # Get intersection
    proteins_with_both = set(protein_go_map.keys()) & set(domain_architectures.keys())

    # Build protein-domain map
    protein_domain_map = {}
    all_domains = set()
    for protein_id in proteins_with_both:
        arch = domain_architectures[protein_id]
        domains = set(arch.single_domains)
        if domains:
            protein_domain_map[protein_id] = domains
            all_domains.update(domains)

    # Get all GO terms
    all_go_terms = set()
    for protein_id in proteins_with_both:
        all_go_terms.update(protein_go_map[protein_id])

    domain_list = sorted(all_domains)
    go_list = sorted(all_go_terms)

    logger.info(f"✓ Dataset prepared: {len(proteins_with_both):,} proteins, {len(domain_list):,} domains, {len(go_list):,} GO terms")
    logger.info(f"  Total tests: {len(domain_list) * len(go_list):,}")

    # Build sparse matrices
    logger.info("")
    logger.info("STAGE 2: Building Sparse Matrices")
    logger.info("─" * 70)
    start_time = time.time()

    protein_domain_matrix, protein_go_matrix = build_sparse_matrices(
        protein_domain_map,
        protein_go_map,
        domain_list,
        go_list
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
    logger.info(f"✓ Contingency tables computed in {table_time:.2f}s ({table_time / 60:.1f} min)")

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
        logger.info(f"  Progress: {completed:,} / {total:,} ({progress_pct:.1f}%) | {rate:,.0f} tests/s | ETA: {eta/60:.1f} min")

    odds_ratios, pvalues = fisher_exact_parallel(
        tables,
        alternative='greater',
        n_jobs=args.num_cores,
        batch_size=args.batch_size,
        progress_callback=progress_callback
    )

    test_time = time.time() - start_time
    logger.info(f"✓ Fisher tests completed in {test_time:.2f}s ({test_time / 60:.1f} min)")
    logger.info(f"  Rate: {len(pvalues) / test_time:,.0f} tests/second")

    # Apply FDR correction
    logger.info("")
    logger.info("STAGE 5: FDR Correction")
    logger.info("─" * 70)
    start_time = time.time()

    adjusted_pvalues, threshold = benjamini_hochberg_correction(pvalues, alpha=args.fdr_threshold)

    fdr_time = time.time() - start_time
    logger.info(f"✓ FDR correction completed in {fdr_time:.2f}s")
    logger.info(f"  Threshold p-value: {threshold:.2e}")

    # Count significant associations
    significant = adjusted_pvalues <= args.fdr_threshold
    n_significant = int(significant.sum())

    # Export results
    logger.info("")
    logger.info("STAGE 6: Exporting Results")
    logger.info("─" * 70)

    # Calculate hypergeometric scores for significant associations
    logger.info("Calculating hypergeometric scores for significant associations...")
    significant_indices = np.where(significant)[0]

    # Export significant associations with hypergeometric scores
    output_file = args.output_dir / "domain_go_associations_significant.tsv"
    with open(output_file, 'w') as f:
        f.write("domain\tgo_term\tp_value\tadj_p_value\todds_ratio\thyper_score\n")
        for idx in significant_indices:
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)

            # Get contingency table values for hypergeometric score
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            hyper_score = calculate_hypergeometric_score(a, b, c, d)

            f.write(f"{domain_list[domain_idx]}\t{go_list[go_idx]}\t"
                   f"{pvalues[idx]:.6e}\t{adjusted_pvalues[idx]:.6e}\t{odds_ratios[idx]:.4f}\t{hyper_score:.2f}\n")

    logger.info(f"✓ Exported significant associations to: {output_file}")
    logger.info(f"  {n_significant:,} associations (FDR < {args.fdr_threshold})")

    # Export top associations with hypergeometric scores
    top_file = args.output_dir / "domain_go_associations_top100.tsv"
    top_indices = np.argsort(pvalues)[:100]
    with open(top_file, 'w') as f:
        f.write("rank\tdomain\tgo_term\tp_value\tadj_p_value\todds_ratio\thyper_score\n")
        for rank, idx in enumerate(top_indices, 1):
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)

            # Get contingency table values for hypergeometric score
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            hyper_score = calculate_hypergeometric_score(a, b, c, d)

            f.write(f"{rank}\t{domain_list[domain_idx]}\t{go_list[go_idx]}\t"
                   f"{pvalues[idx]:.6e}\t{adjusted_pvalues[idx]:.6e}\t{odds_ratios[idx]:.4f}\t{hyper_score:.2f}\n")

    logger.info(f"✓ Exported top 100 associations to: {top_file}")

    # Performance summary
    total_time = matrix_time + table_time + test_time + fdr_time

    logger.info("")
    logger.info("=" * 70)
    logger.info("PIPELINE COMPLETE!")
    logger.info("=" * 70)
    logger.info("Results Summary:")
    logger.info(f"  Total domain-GO tests: {len(pvalues):,}")
    logger.info(f"  Significant associations (FDR < {args.fdr_threshold}): {n_significant:,} ({n_significant / len(pvalues) * 100:.2f}%)")
    logger.info(f"  Total runtime: {total_time:.1f}s ({total_time / 60:.1f} minutes)")
    logger.info("")
    logger.info("Output files:")
    logger.info(f"  {output_file}")
    logger.info(f"  {top_file}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
