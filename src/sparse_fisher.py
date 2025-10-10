"""
Optimized sparse matrix implementation for dcGO statistical inference.

This approach uses sparse matrices to efficiently compute contingency tables
for all domain-GO combinations simultaneously.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Set, Tuple

import numpy as np
from loguru import logger
from scipy import sparse


class DomainType(Enum):
    """Classification of domain features."""

    SINGLE = "single"  # Individual domain (e.g., "IPR000001")
    SUPRA_PAIR = "supra_pair"  # 2-domain combination (e.g., "IPR000001,IPR000002")
    SUPRA_TRIPLE = "supra_triple"  # 3-domain combination


@dataclass
class DomainMetadata:
    """
    Metadata about domains in the analysis.

    This tracks information needed for hierarchical statistical inference,
    allowing supra-domains to "borrow strength" from their constituent domains.
    """

    domain_id: str
    domain_type: DomainType
    constituent_domains: List[str]  # Empty for single domains
    observation_count: int  # Number of proteins with this domain
    index: int  # Position in domain_list array

    @property
    def is_single_domain(self) -> bool:
        """Check if this is a single domain (not a supra-domain)."""
        return self.domain_type == DomainType.SINGLE

    @property
    def is_supra_domain(self) -> bool:
        """Check if this is a supra-domain (multi-domain combination)."""
        return self.domain_type in (DomainType.SUPRA_PAIR, DomainType.SUPRA_TRIPLE)


def parse_domain_id(domain_id: str) -> Tuple[DomainType, List[str]]:
    """
    Parse a domain ID to determine its type and constituent domains.

    Args:
        domain_id: Domain identifier (e.g., "IPR000001" or "IPR000001,IPR000002")

    Returns:
        Tuple of (domain_type, constituent_domains)

    Examples:
        >>> parse_domain_id("IPR000001")
        (DomainType.SINGLE, ["IPR000001"])

        >>> parse_domain_id("IPR000001,IPR000002")
        (DomainType.SUPRA_PAIR, ["IPR000001", "IPR000002"])
    """
    if "," not in domain_id:
        return DomainType.SINGLE, [domain_id]

    constituents = domain_id.split(",")
    num_constituents = len(constituents)

    if num_constituents == 2:
        return DomainType.SUPRA_PAIR, constituents
    elif num_constituents == 3:
        return DomainType.SUPRA_TRIPLE, constituents
    else:
        # Fallback for edge cases (shouldn't happen with max_supra_domain_length=3)
        logger.warning(
            f"Unexpected supra-domain length {num_constituents} for {domain_id}"
        )
        return DomainType.SUPRA_TRIPLE, constituents


def build_domain_metadata(
    domain_list: List[str], protein_domains: Dict[str, Set[str]]
) -> Dict[str, DomainMetadata]:
    """
    Build comprehensive metadata for all domains in the analysis.

    Args:
        domain_list: Ordered list of all domain IDs (single + supra)
        protein_domains: Dict mapping protein_id -> set of domain IDs

    Returns:
        Dictionary mapping domain_id -> DomainMetadata
    """
    logger.info("Building domain metadata...")

    metadata = {}

    # Count observations for each domain
    domain_counts: Dict[str, int] = {}
    for domains in protein_domains.values():
        for domain in domains:
            domain_counts[domain] = domain_counts.get(domain, 0) + 1

    # Build metadata for each domain
    single_count = 0
    supra_pair_count = 0
    supra_triple_count = 0

    for idx, domain_id in enumerate(domain_list):
        domain_type, constituents = parse_domain_id(domain_id)

        # Track constituent domains (empty list for single domains)
        constituent_list = constituents if domain_type != DomainType.SINGLE else []

        metadata[domain_id] = DomainMetadata(
            domain_id=domain_id,
            domain_type=domain_type,
            constituent_domains=constituent_list,
            observation_count=domain_counts.get(domain_id, 0),
            index=idx,
        )

        # Count by type
        if domain_type == DomainType.SINGLE:
            single_count += 1
        elif domain_type == DomainType.SUPRA_PAIR:
            supra_pair_count += 1
        elif domain_type == DomainType.SUPRA_TRIPLE:
            supra_triple_count += 1

    logger.info(f"  Total domains: {len(metadata):,}")
    logger.info(f"    Single domains: {single_count:,}")
    logger.info(f"    Supra-domain pairs: {supra_pair_count:,}")
    logger.info(f"    Supra-domain triples: {supra_triple_count:,}")

    return metadata


def build_sparse_matrices(
    protein_domains: Dict[str, Set[str]],
    protein_go: Dict[str, Set[str]],
    domain_list: List[str],
    go_list: List[str],
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix, Dict[str, DomainMetadata]]:
    """
    Build sparse binary matrices for protein-domain and protein-GO relationships.

    Args:
        protein_domains: Dict mapping protein_id -> set of domain IDs
        protein_go: Dict mapping protein_id -> set of GO term IDs
        domain_list: Ordered list of all domain IDs (includes both single and supra-domains)
        go_list: Ordered list of all GO term IDs

    Returns:
        Tuple of (protein_domain_matrix, protein_go_matrix, domain_metadata)
        - protein_domain_matrix: Binary matrix (n_proteins, n_domains)
        - protein_go_matrix: Binary matrix (n_proteins, n_go_terms)
        - domain_metadata: Metadata for hierarchical inference
    """
    # Get all proteins
    all_proteins = sorted(set(protein_domains.keys()) | set(protein_go.keys()))
    protein_to_idx = {p: i for i, p in enumerate(all_proteins)}
    domain_to_idx = {d: i for i, d in enumerate(domain_list)}
    go_to_idx = {g: i for i, g in enumerate(go_list)}

    n_proteins = len(all_proteins)
    n_domains = len(domain_list)
    n_go_terms = len(go_list)

    logger.info(
        f"Building sparse matrices: {n_proteins:,} proteins × {n_domains:,} domains"
    )
    logger.info(
        f"                         {n_proteins:,} proteins × {n_go_terms:,} GO terms"
    )

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
        dtype=np.int8,
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
        dtype=np.int8,
    )

    logger.info(
        f"Protein-domain matrix: {protein_domain_matrix.nnz:,} non-zero entries"
    )
    logger.info(f"Protein-GO matrix: {protein_go_matrix.nnz:,} non-zero entries")

    # Build domain metadata for hierarchical inference
    domain_metadata = build_domain_metadata(domain_list, protein_domains)

    return protein_domain_matrix, protein_go_matrix, domain_metadata


def compute_contingency_tables_sparse(
    protein_domain_matrix: sparse.csr_matrix, protein_go_matrix: sparse.csr_matrix
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

    logger.info(
        f"Computing contingency tables for {n_domains:,} × {n_go_terms:,} = {n_domains * n_go_terms:,} pairs"
    )

    # Compute a: proteins with both (domain AND GO)
    # This is the dot product of transposed matrices
    logger.info("  Computing overlap counts (a)...")
    a_matrix = (
        protein_domain_matrix.T @ protein_go_matrix
    )  # Shape: (n_domains, n_go_terms)
    a_matrix = a_matrix.toarray().astype(np.int32)  # Convert to dense for faster access

    # Compute marginal counts
    logger.info("  Computing marginal counts...")
    domain_counts = (
        np.array(protein_domain_matrix.sum(axis=0)).flatten().astype(np.int32)
    )  # Proteins per domain
    go_counts = (
        np.array(protein_go_matrix.sum(axis=0)).flatten().astype(np.int32)
    )  # Proteins per GO term

    # Build all contingency tables
    logger.info("  Building contingency table array...")
    n_tests = n_domains * n_go_terms
    tables = np.zeros((n_tests, 2, 2), dtype=np.int32)

    idx = 0
    for d_idx in range(n_domains):
        n_with_domain = domain_counts[d_idx]

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
