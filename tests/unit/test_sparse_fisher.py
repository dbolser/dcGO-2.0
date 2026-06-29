"""
Unit tests for sparse matrix Fisher's exact test implementation.

Tests the core sparse matrix operations used in the production pipeline.
"""

import pytest
import numpy as np
from scipy import sparse

from src.sparse_fisher import build_sparse_matrices, compute_contingency_tables_sparse


class TestBuildSparseMatrices:
    """Test suite for building sparse binary matrices."""

    @pytest.fixture
    def sample_data(self):
        """Create sample protein-domain and protein-GO mappings."""
        protein_domains = {
            "P001": {"IPR001", "IPR002"},
            "P002": {"IPR001"},
            "P003": {"IPR002", "IPR003"},
            "P004": {"IPR003"},
            "P005": {"IPR001", "IPR003"},
        }

        protein_go = {
            "P001": {"GO:0001", "GO:0002"},
            "P002": {"GO:0001"},
            "P003": {"GO:0002", "GO:0003"},
            "P004": {"GO:0003"},
            "P005": {"GO:0001", "GO:0003"},
        }

        domain_list = ["IPR001", "IPR002", "IPR003"]
        go_list = ["GO:0001", "GO:0002", "GO:0003"]

        return protein_domains, protein_go, domain_list, go_list

    def test_matrix_shapes(self, sample_data):
        """Test that matrices have correct dimensions."""
        protein_domains, protein_go, domain_list, go_list = sample_data

        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, domain_list, go_list
        )

        n_proteins = 5
        n_domains = 3
        n_go_terms = 3

        assert domain_matrix.shape == (n_proteins, n_domains)
        assert go_matrix.shape == (n_proteins, n_go_terms)

    def test_matrix_is_binary(self, sample_data):
        """Test that matrices contain only 0s and 1s."""
        protein_domains, protein_go, domain_list, go_list = sample_data

        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, domain_list, go_list
        )

        # Convert to dense for checking
        domain_dense = domain_matrix.toarray()
        go_dense = go_matrix.toarray()

        assert np.all((domain_dense == 0) | (domain_dense == 1))
        assert np.all((go_dense == 0) | (go_dense == 1))

    def test_matrix_sparsity(self, sample_data):
        """Test that matrices are actually sparse."""
        protein_domains, protein_go, domain_list, go_list = sample_data

        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, domain_list, go_list
        )

        # Check matrix types
        assert isinstance(domain_matrix, sparse.csr_matrix)
        assert isinstance(go_matrix, sparse.csr_matrix)

        # Check sparsity (should have fewer non-zeros than total elements)
        total_elements_d = domain_matrix.shape[0] * domain_matrix.shape[1]
        total_elements_g = go_matrix.shape[0] * go_matrix.shape[1]

        assert domain_matrix.nnz < total_elements_d
        assert go_matrix.nnz < total_elements_g

    def test_correct_annotations(self, sample_data):
        """Test that matrices encode the correct protein annotations."""
        protein_domains, protein_go, domain_list, go_list = sample_data

        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, domain_list, go_list
        )

        domain_dense = domain_matrix.toarray()
        go_dense = go_matrix.toarray()

        # P001 should have IPR001 and IPR002
        assert domain_dense[0, 0] == 1  # IPR001
        assert domain_dense[0, 1] == 1  # IPR002
        assert domain_dense[0, 2] == 0  # Not IPR003

        # P001 should have GO:0001 and GO:0002
        assert go_dense[0, 0] == 1  # GO:0001
        assert go_dense[0, 1] == 1  # GO:0002
        assert go_dense[0, 2] == 0  # Not GO:0003

    def test_missing_proteins(self):
        """Test handling of proteins with only domains or only GO terms."""
        protein_domains = {
            "P001": {"IPR001"},
            "P002": {"IPR002"},
        }

        protein_go = {
            "P002": {"GO:0001"},
            "P003": {"GO:0002"},
        }

        domain_list = ["IPR001", "IPR002"]
        go_list = ["GO:0001", "GO:0002"]

        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, domain_list, go_list
        )

        # Should have all 3 proteins
        assert domain_matrix.shape[0] == 3
        assert go_matrix.shape[0] == 3

    def test_empty_input(self):
        """Test handling of empty input."""
        protein_domains = {}
        protein_go = {}
        domain_list = []
        go_list = []

        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, domain_list, go_list
        )

        assert domain_matrix.shape == (0, 0)
        assert go_matrix.shape == (0, 0)


class TestComputeContingencyTables:
    """Test suite for computing contingency tables from sparse matrices."""

    @pytest.fixture
    def sample_matrices(self):
        """Create sample sparse matrices for testing."""
        # Create simple test case with known contingency tables
        # 4 proteins, 2 domains, 2 GO terms

        # Protein-domain matrix:
        # P1: D1
        # P2: D1, D2
        # P3: D2
        # P4: neither
        protein_domain = sparse.csr_matrix(
            np.array(
                [
                    [1, 0],  # P1
                    [1, 1],  # P2
                    [0, 1],  # P3
                    [0, 0],  # P4
                ],
                dtype=np.int8,
            )
        )

        # Protein-GO matrix:
        # P1: GO1
        # P2: GO1, GO2
        # P3: GO2
        # P4: neither
        protein_go = sparse.csr_matrix(
            np.array(
                [
                    [1, 0],  # P1
                    [1, 1],  # P2
                    [0, 1],  # P3
                    [0, 0],  # P4
                ],
                dtype=np.int8,
            )
        )

        return protein_domain, protein_go

    def test_table_shape(self, sample_matrices):
        """Test that contingency tables have correct shape."""
        protein_domain, protein_go = sample_matrices

        tables = compute_contingency_tables_sparse(protein_domain, protein_go)

        n_domains = 2
        n_go_terms = 2
        expected_n_tests = n_domains * n_go_terms

        assert tables.shape == (expected_n_tests, 2, 2)

    def test_table_values_sum_to_total(self, sample_matrices):
        """Test that each contingency table sums to total number of proteins."""
        protein_domain, protein_go = sample_matrices

        tables = compute_contingency_tables_sparse(protein_domain, protein_go)
        n_proteins = protein_domain.shape[0]

        for table in tables:
            assert table.sum() == n_proteins

    def test_table_values_nonnegative(self, sample_matrices):
        """Test that all contingency table values are non-negative."""
        protein_domain, protein_go = sample_matrices

        tables = compute_contingency_tables_sparse(protein_domain, protein_go)

        assert np.all(tables >= 0)

    def test_specific_contingency_values(self, sample_matrices):
        """Test specific contingency table values against manual calculation."""
        protein_domain, protein_go = sample_matrices

        tables = compute_contingency_tables_sparse(protein_domain, protein_go)

        # Test D1-GO1 (first table, index 0)
        # D1: P1, P2 (2 proteins)
        # GO1: P1, P2 (2 proteins)
        # Both: P1, P2 (2 proteins)
        # a=2, b=0, c=0, d=2
        table_d1_go1 = tables[0]
        assert table_d1_go1[0, 0] == 2  # a: both
        assert table_d1_go1[0, 1] == 0  # b: domain only
        assert table_d1_go1[1, 0] == 0  # c: GO only
        assert table_d1_go1[1, 1] == 2  # d: neither

    def test_marginal_consistency(self, sample_matrices):
        """Test that marginals in contingency tables match matrix row/column sums."""
        protein_domain, protein_go = sample_matrices

        tables = compute_contingency_tables_sparse(protein_domain, protein_go)

        domain_counts = np.array(protein_domain.sum(axis=0)).flatten()
        go_counts = np.array(protein_go.sum(axis=0)).flatten()

        n_domains = protein_domain.shape[1]
        n_go_terms = protein_go.shape[1]

        for d_idx in range(n_domains):
            for g_idx in range(n_go_terms):
                table_idx = d_idx * n_go_terms + g_idx
                table = tables[table_idx]

                # Domain marginal: a + b
                assert table[0, 0] + table[0, 1] == domain_counts[d_idx]

                # GO marginal: a + c
                assert table[0, 0] + table[1, 0] == go_counts[g_idx]

    def test_large_sparse_matrix(self):
        """Test with larger sparse matrices to verify scalability."""
        n_proteins = 1000
        n_domains = 100
        n_go_terms = 100

        # Create random sparse matrices with ~1% density
        np.random.seed(42)
        density = 0.01

        protein_domain = sparse.random(
            n_proteins, n_domains, density=density, format="csr", dtype=np.int8
        )
        protein_go = sparse.random(
            n_proteins, n_go_terms, density=density, format="csr", dtype=np.int8
        )

        # Make binary
        protein_domain.data = np.ones_like(protein_domain.data)
        protein_go.data = np.ones_like(protein_go.data)

        tables = compute_contingency_tables_sparse(protein_domain, protein_go)

        # Verify shape
        assert tables.shape == (n_domains * n_go_terms, 2, 2)

        # Verify all tables sum to n_proteins
        assert np.all(tables.sum(axis=(1, 2)) == n_proteins)

        # Verify all values non-negative
        assert np.all(tables >= 0)


class TestIntegration:
    """Integration tests for the complete sparse matrix workflow."""

    def test_full_pipeline(self):
        """Test complete workflow from protein annotations to contingency tables."""
        # Create realistic test data
        protein_domains = {
            f"P{i:04d}": {f"IPR{j:03d}" for j in range(i % 5, i % 5 + 3)}
            for i in range(100)
        }

        protein_go = {
            f"P{i:04d}": {f"GO:{j:07d}" for j in range(i % 7, i % 7 + 2)}
            for i in range(100)
        }

        # Get all unique domains and GO terms
        all_domains = sorted(
            set(d for domains in protein_domains.values() for d in domains)
        )
        all_go_terms = sorted(
            set(g for go_terms in protein_go.values() for g in go_terms)
        )

        # Build matrices
        domain_matrix, go_matrix, _ = build_sparse_matrices(
            protein_domains, protein_go, all_domains, all_go_terms
        )

        # Compute contingency tables
        tables = compute_contingency_tables_sparse(domain_matrix, go_matrix)

        # Verify results
        assert domain_matrix.shape[0] == 100  # 100 proteins
        assert go_matrix.shape[0] == 100
        assert tables.shape == (len(all_domains) * len(all_go_terms), 2, 2)
        assert np.all(tables >= 0)
        assert np.all(tables.sum(axis=(1, 2)) == 100)
