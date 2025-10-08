"""
Optimized sparse matrix implementation for dcGO statistical inference.

This approach uses sparse matrices to efficiently compute contingency tables
for all domain-GO combinations simultaneously.
"""

import numpy as np
from scipy import sparse
from typing import Dict, Set, Tuple, List
from loguru import logger


def build_sparse_matrices(
    protein_domains: Dict[str, Set[str]],
    protein_go: Dict[str, Set[str]],
    domain_list: List[str],
    go_list: List[str]
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """
    Build sparse binary matrices for protein-domain and protein-GO relationships.

    Args:
        protein_domains: Dict mapping protein_id -> set of domain IDs
        protein_go: Dict mapping protein_id -> set of GO term IDs
        domain_list: Ordered list of all domain IDs
        go_list: Ordered list of all GO term IDs

    Returns:
        Tuple of (protein_domain_matrix, protein_go_matrix)
        Both matrices have shape (n_proteins, n_features) and are binary (0/1)
    """
    # Get all proteins
    all_proteins = sorted(set(protein_domains.keys()) | set(protein_go.keys()))
    protein_to_idx = {p: i for i, p in enumerate(all_proteins)}
    domain_to_idx = {d: i for i, d in enumerate(domain_list)}
    go_to_idx = {g: i for i, g in enumerate(go_list)}

    n_proteins = len(all_proteins)
    n_domains = len(domain_list)
    n_go_terms = len(go_list)

    logger.info(f"Building sparse matrices: {n_proteins:,} proteins × {n_domains:,} domains")
    logger.info(f"                         {n_proteins:,} proteins × {n_go_terms:,} GO terms")

    # Build protein-domain matrix
    rows_d, cols_d = [], []
    for protein_id, domains in protein_domains.items():
        if protein_id in protein_to_idx:
            p_idx = protein_to_idx[protein_id]
            for domain in domains:
                if domain in domain_to_idx:
                    d_idx = domain_to_idx[domain]
                    rows_d.append(p_idx)
                    cols_d.append(d_idx)

    protein_domain_matrix = sparse.csr_matrix(
        (np.ones(len(rows_d), dtype=np.int8), (rows_d, cols_d)),
        shape=(n_proteins, n_domains),
        dtype=np.int8
    )

    # Build protein-GO matrix
    rows_g, cols_g = [], []
    for protein_id, go_terms in protein_go.items():
        if protein_id in protein_to_idx:
            p_idx = protein_to_idx[protein_id]
            for go_term in go_terms:
                if go_term in go_to_idx:
                    g_idx = go_to_idx[go_term]
                    rows_g.append(p_idx)
                    cols_g.append(g_idx)

    protein_go_matrix = sparse.csr_matrix(
        (np.ones(len(rows_g), dtype=np.int8), (rows_g, cols_g)),
        shape=(n_proteins, n_go_terms),
        dtype=np.int8
    )

    logger.info(f"Protein-domain matrix: {protein_domain_matrix.nnz:,} non-zero entries")
    logger.info(f"Protein-GO matrix: {protein_go_matrix.nnz:,} non-zero entries")

    return protein_domain_matrix, protein_go_matrix


def compute_contingency_tables_sparse(
    protein_domain_matrix: sparse.csr_matrix,
    protein_go_matrix: sparse.csr_matrix
) -> np.ndarray:
    """
    Compute all 2x2 contingency tables using sparse matrix operations.

    For each domain-GO pair, we compute:
                    Has GO | No GO
        Has domain:   a   |   b
        No domain:    c   |   d

    Where:
    - a = proteins with both domain and GO term
    - b = proteins with domain but not GO term
    - c = proteins with GO term but not domain
    - d = proteins with neither

    Args:
        protein_domain_matrix: Binary matrix (n_proteins, n_domains)
        protein_go_matrix: Binary matrix (n_proteins, n_go_terms)

    Returns:
        Array of shape (n_domains * n_go_terms, 2, 2) containing contingency tables
    """
    n_proteins = protein_domain_matrix.shape[0]
    n_domains = protein_domain_matrix.shape[1]
    n_go_terms = protein_go_matrix.shape[1]

    logger.info(f"Computing contingency tables for {n_domains:,} × {n_go_terms:,} = {n_domains * n_go_terms:,} pairs")

    # Compute a: proteins with both (domain AND GO)
    # This is the dot product of transposed matrices
    logger.info("  Computing overlap counts (a)...")
    a_matrix = protein_domain_matrix.T @ protein_go_matrix  # Shape: (n_domains, n_go_terms)
    a_matrix = a_matrix.toarray().astype(np.int32)  # Convert to dense for faster access

    # Compute marginal counts
    logger.info("  Computing marginal counts...")
    domain_counts = np.array(protein_domain_matrix.sum(axis=0)).flatten().astype(np.int32)  # Proteins per domain
    go_counts = np.array(protein_go_matrix.sum(axis=0)).flatten().astype(np.int32)  # Proteins per GO term

    # Build all contingency tables
    logger.info("  Building contingency table array...")
    n_tests = n_domains * n_go_terms
    tables = np.zeros((n_tests, 2, 2), dtype=np.int32)

    idx = 0
    for d_idx in range(n_domains):
        n_with_domain = domain_counts[d_idx]
        n_without_domain = n_proteins - n_with_domain

        for g_idx in range(n_go_terms):
            n_with_go = go_counts[g_idx]

            # Fill contingency table
            a = a_matrix[d_idx, g_idx]  # Both domain and GO
            b = n_with_domain - a  # Domain but not GO
            c = n_with_go - a  # GO but not domain
            d = n_proteins - a - b - c  # Neither

            tables[idx] = [[a, b], [c, d]]
            idx += 1

        # Progress logging
        if (d_idx + 1) % 1000 == 0:
            logger.info(f"    Processed {d_idx + 1:,} / {n_domains:,} domains")

    logger.info(f"✓ Built {n_tests:,} contingency tables")

    return tables
