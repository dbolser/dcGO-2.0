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

    # Count the number of distinct PROTEINS carrying each domain. The per-protein
    # domain list contains repeats (one entry per member signature), so dedupe
    # within each protein — otherwise observation_count is an occurrence count
    # that can exceed the number of proteins.
    domain_counts: Dict[str, int] = {}
    for domains in protein_domains.values():
        for domain in set(domains):
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
    # A (protein, InterPro entry) pair is listed once per supporting member
    # signature in protein2ipr, so the same pair appears many times and the CSR
    # constructor SUMS the duplicates (cells reach values like 125). The
    # contingency analysis is presence/absence, so collapse to 1 — otherwise
    # overlap counts are inflated and a domain's "count" can exceed the number
    # of proteins, driving the `d` cell negative.
    protein_domain_matrix.data[:] = 1

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
    # Same presence/absence guard as the domain matrix (defensive; GO maps are
    # usually already de-duplicated).
    protein_go_matrix.data[:] = 1

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

    # The indicator matrices are stored as int8 to save memory, but int8
    # *accumulation* overflows at 127. A common domain/GO pair co-occurs in far
    # more than 127 proteins, and marginal counts run to the whole proteome, so
    # both the overlap matmul and the marginal sums must accumulate in a wider
    # dtype. Upcast to int32 BEFORE any accumulation. (Casting the int8 product
    # afterwards is too late — the overflow has already happened.)
    protein_domain_matrix = protein_domain_matrix.astype(np.int32)
    protein_go_matrix = protein_go_matrix.astype(np.int32)

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

    # Build all contingency tables via NumPy broadcasting (no Python loop).
    # b, c, d follow from a and the marginals; ravel() in C-order keeps the
    # domain-major / GO-minor ordering the rest of the pipeline expects.
    logger.info("  Building contingency table array...")
    n_tests = n_domains * n_go_terms

    b_matrix = domain_counts[:, np.newaxis] - a_matrix  # Domain but not GO
    c_matrix = go_counts[np.newaxis, :] - a_matrix  # GO but not domain
    d_matrix = (  # Neither
        n_proteins - domain_counts[:, np.newaxis] - go_counts[np.newaxis, :] + a_matrix
    )

    tables = np.empty((n_tests, 2, 2), dtype=np.int32)
    tables[:, 0, 0] = a_matrix.ravel()
    tables[:, 0, 1] = b_matrix.ravel()
    tables[:, 1, 0] = c_matrix.ravel()
    tables[:, 1, 1] = d_matrix.ravel()

    logger.info(f"✓ Built {n_tests:,} contingency tables")

    return tables


def compute_cooccurring_contingency_tables(
    protein_domain_matrix: sparse.csr_matrix,
    protein_go_matrix: sparse.csr_matrix,
    domain_block: int = 4096,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """Contingency tables for the domain-term pairs that actually co-occur.

    ``compute_contingency_tables_sparse`` materialises the **dense** domain x
    term product: one int32 2x2 table per pair whether or not the pair is ever
    seen together. That is what puts supra-domains out of reach on a large
    universe — an all-species run is 464,490 x 28,112 = 13.1e9 tables, ~389 GB
    of tables alone before the p-values and the BH sort.

    Nearly all of those tables have ``a = 0``, and dcGO's Fisher test is
    one-sided in the ``greater`` direction (``FISHER_ALTERNATIVE``). For a=0 the
    one-sided p-value is ``P(X >= 0) = 1`` **exactly**, at any marginals. So the
    a=0 tables are not approximated or discarded here: their p-value is known in
    closed form, and enumerating only the co-occurring pairs is exact rather
    than a heuristic.

    Two things are needed to keep it exact, and both are the caller's job:

    * **BH must still divide by the full hypothesis count.** The returned
      ``n_hypotheses`` is ``n_domains * n_terms``; pass it to
      ``benjamini_hochberg_correction`` as ``n_hypotheses``. Correcting against
      the enumerated subset instead would inflate every rejection.
    * **This is only valid for ``alternative='greater'``.** A two-sided or
      ``less`` test gives a=0 pairs p < 1, and dropping them would change the
      answer. Callers using another alternative must use the dense path.

    The product is computed in blocks of ``domain_block`` columns rather than in
    one matmul: SciPy's sparse-sparse product allocates its full result up
    front, and on the all-species supra design that single allocation is the
    thing most likely to kill the run.

    Args:
        protein_domain_matrix: Binary (n_proteins, n_domains) CSR.
        protein_go_matrix: Binary (n_proteins, n_terms) CSR.
        domain_block: Domain columns per block. Trades peak memory for overhead.

    Returns:
        ``(tables, pair_index, n_hypotheses)``. ``tables`` has shape (k, 2, 2)
        for the k co-occurring pairs. ``pair_index`` is int64 and holds each
        pair's index in the **dense** domain-major layout
        (``domain_idx * n_terms + term_idx``), so callers decompose it exactly
        as they did the dense index. ``n_hypotheses`` is the dense test count.
    """
    n_proteins, n_domains = protein_domain_matrix.shape
    n_terms = protein_go_matrix.shape[1]
    n_hypotheses = n_domains * n_terms

    logger.info(
        f"Enumerating co-occurring pairs from {n_domains:,} x {n_terms:,} = "
        f"{n_hypotheses:,} hypotheses (a=0 pairs have p=1 exactly under a "
        f"one-sided 'greater' test and are not enumerated)"
    )

    # int8 accumulates to 127 and overflows; a common domain/term pair
    # co-occurs in far more proteins than that. Upcast BEFORE any product.
    domain_i32 = protein_domain_matrix.astype(np.int32)
    go_i32 = protein_go_matrix.astype(np.int32)

    domain_counts = np.asarray(domain_i32.sum(axis=0)).ravel().astype(np.int64)
    go_counts = np.asarray(go_i32.sum(axis=0)).ravel().astype(np.int64)

    domain_csc = domain_i32.tocsc()

    table_blocks: List[np.ndarray] = []
    index_blocks: List[np.ndarray] = []
    n_pairs = 0

    for start in range(0, n_domains, domain_block):
        stop = min(start + domain_block, n_domains)
        # (block_domains, n_terms), sparse: only co-occurring pairs are stored.
        overlap = (domain_csc[:, start:stop].T @ go_i32).tocoo()
        if overlap.nnz == 0:
            continue

        a = overlap.data.astype(np.int64, copy=False)
        rows = overlap.row.astype(np.int64, copy=False) + start
        cols = overlap.col.astype(np.int64, copy=False)

        block = np.empty((overlap.nnz, 2, 2), dtype=np.int32)
        block[:, 0, 0] = a
        block[:, 0, 1] = domain_counts[rows] - a
        block[:, 1, 0] = go_counts[cols] - a
        block[:, 1, 1] = n_proteins - domain_counts[rows] - go_counts[cols] + a

        table_blocks.append(block)
        index_blocks.append(rows * n_terms + cols)
        n_pairs += overlap.nnz

        if (start // domain_block) % 20 == 0:
            logger.info(
                f"  domains {stop:,}/{n_domains:,} — {n_pairs:,} co-occurring "
                f"pairs so far"
            )

    if not table_blocks:
        return (
            np.empty((0, 2, 2), dtype=np.int32),
            np.empty(0, dtype=np.int64),
            n_hypotheses,
        )

    tables = np.concatenate(table_blocks)
    pair_index = np.concatenate(index_blocks)
    del table_blocks, index_blocks

    # COO from a CSC-major product is already domain-major, but the block loop
    # only guarantees ordering *between* blocks. Sorting makes pair_index
    # ascending overall, so downstream decomposition and any grouping by domain
    # sees the same order the dense path produced.
    order = np.argsort(pair_index, kind="stable")
    tables = tables[order]
    pair_index = pair_index[order]

    logger.info(
        f"✓ Built {n_pairs:,} contingency tables for co-occurring pairs "
        f"({n_pairs / n_hypotheses:.4%} of the hypothesis space); the remaining "
        f"{n_hypotheses - n_pairs:,} have p=1 by construction"
    )
    return tables, pair_index, n_hypotheses
