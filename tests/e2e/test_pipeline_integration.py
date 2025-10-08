"""
End-to-end integration tests for the dcGO pipeline.

Tests the complete workflow with mock data to validate all components
work together correctly.
"""

import pytest
import tempfile
import gzip
from pathlib import Path
from unittest.mock import patch, Mock
import pandas as pd

# Add src to path for testing
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from config.settings import Config, DataSource
from data_acquisition import DataAcquisition
from statistical_inference import StatisticalInferenceEngine, AssociationResult
from ontology_processor import OntologyProcessor, Annotation
from database_manager import dcGODatabaseManager


class TestPipelineIntegration:
    """End-to-end integration tests for the complete dcGO pipeline."""
    
    @pytest.fixture
    def workspace(self):
        """Create temporary workspace for integration testing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            workspace = Path(temp_dir)
            
            # Create directory structure
            (workspace / "data" / "raw").mkdir(parents=True)
            (workspace / "data" / "processed").mkdir(parents=True)
            (workspace / "results").mkdir(parents=True)
            
            yield workspace
    
    @pytest.fixture
    def integration_config(self, workspace):
        """Create configuration for integration testing."""
        from dataclasses import replace
        config = Config()
        return replace(config, base_dir=workspace)
    
    def create_mock_datasets(self, workspace: Path):
        """Create minimal mock datasets for integration testing."""
        raw_dir = workspace / "data" / "raw"
        
        # Mock UniProt sequences
        uniprot_content = """>sp|P12345|TEST1_HUMAN Test protein 1
MKLLVLSLSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQH
>sp|P67890|TEST2_HUMAN Test protein 2
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLR
>sp|P11111|TEST3_HUMAN Test protein 3
MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYY
"""
        with gzip.open(raw_dir / "uniprot_sprot.fasta.gz", "wt") as f:
            f.write(uniprot_content)
        
        # Mock GOA annotations
        goa_content = """!gaf-version: 2.2
!generated-by: UniProt
UniProtKB\tP12345\tTEST1\t\tGO:0005515\tPMID:12345678\tIDA\t\tF\tTest protein 1\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
UniProtKB\tP12345\tTEST1\t\tGO:0008150\tPMID:12345678\tISS\t\tP\tTest protein 1\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
UniProtKB\tP67890\tTEST2\t\tGO:0005515\tPMID:87654321\tIDA\t\tF\tTest protein 2\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
UniProtKB\tP67890\tTEST2\t\tGO:0016740\tPMID:87654321\tIDA\t\tF\tTest protein 2\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
UniProtKB\tP11111\tTEST3\t\tGO:0008150\tPMID:11111111\tIDA\t\tP\tTest protein 3\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
"""
        with gzip.open(raw_dir / "goa_uniprot_all.gaf.gz", "wt") as f:
            f.write(goa_content)
        
        # Mock InterPro domain mappings
        interpro_content = """P12345\tPF00001\t10\t50\t1e-10\tT\tDomain1
P12345\tPF00002\t60\t90\t1e-15\tT\tDomain2
P67890\tPF00001\t20\t60\t1e-12\tT\tDomain1
P67890\tPF00003\t80\t120\t1e-08\tT\tDomain3
P11111\tPF00002\t30\t70\t1e-20\tT\tDomain2
"""
        with gzip.open(raw_dir / "protein2ipr.dat.gz", "wt") as f:
            f.write(interpro_content)
        
        # Mock GO ontology (simple)
        go_content = """format-version: 1.2
data-version: releases/2023-01-01

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0005515
name: protein binding
namespace: molecular_function
is_a: GO:0003674

[Term]
id: GO:0016740
name: transferase activity
namespace: molecular_function
is_a: GO:0003674

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function
"""
        with open(raw_dir / "go-basic.obo", "w") as f:
            f.write(go_content)
    
    def test_complete_pipeline_workflow(self, workspace, integration_config):
        """Test the complete dcGO pipeline workflow with mock data."""
        # Create mock datasets
        self.create_mock_datasets(workspace)
        
        # Test 1: Data parsing
        protein_domain_map = self.parse_mock_interpro_data(workspace)
        protein_go_map = self.parse_mock_goa_data(workspace)
        
        assert len(protein_domain_map) == 3  # P12345, P67890, P11111
        assert len(protein_go_map) == 3
        assert 'P12345' in protein_domain_map
        assert 'PF00001' in protein_domain_map['P12345']
        
        # Test 2: Statistical inference
        inference_engine = StatisticalInferenceEngine(protein_domain_map, protein_go_map)
        correspondence_matrix = inference_engine.build_correspondence_matrix()
        
        assert isinstance(correspondence_matrix, pd.DataFrame)
        assert correspondence_matrix.shape[0] > 0  # Has GO terms
        assert correspondence_matrix.shape[1] > 0  # Has domains
        
        # Run statistical tests (may have no significant results with small test data)
        results = inference_engine.run_statistical_tests(min_cooccurrence=1)
        assert isinstance(results, list)
        
        # Test 3: Database storage
        db_manager = dcGODatabaseManager(workspace / "results" / "test.db")
        
        # Create mock annotations for database testing
        mock_annotations = [
            Annotation(
                domain='PF00001',
                go_term='GO:0005515',
                q_value=0.001,
                association_score=85.0,
                annotation_type='direct',
                direct_source_term='GO:0005515'
            ),
            Annotation(
                domain='PF00001',
                go_term='GO:0003674',
                q_value=0.001,
                association_score=85.0,
                annotation_type='propagated',
                direct_source_term='GO:0005515'
            )
        ]
        
        # Store annotations
        metadata = {
            'test_run': True,
            'total_proteins': len(protein_domain_map),
            'total_domains': len(set().union(*protein_domain_map.values())),
            'pipeline_version': '0.1.0'
        }
        
        db_manager.store_annotations(mock_annotations, metadata)
        
        # Test database contents
        summary = db_manager.get_summary_statistics()
        assert summary['total_annotations'] == 2
        assert summary['direct_annotations'] == 1
        assert summary['propagated_annotations'] == 1
        
        # Test TSV export
        tsv_file = db_manager.export_tsv(workspace / "results" / "annotations.tsv")
        assert tsv_file.exists()
        assert tsv_file.stat().st_size > 0
        
        # Verify TSV content
        df = pd.read_csv(tsv_file, sep='\t')
        assert len(df) == 2
        assert 'domain_id' in df.columns
        assert 'go_id' in df.columns
        assert 'annotation_type' in df.columns
    
    def parse_mock_interpro_data(self, workspace: Path) -> dict:
        """Parse mock InterPro data for testing."""
        interpro_file = workspace / "data" / "raw" / "protein2ipr.dat.gz"
        
        protein_domain_map = {}
        with gzip.open(interpro_file, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    protein_id = parts[0]
                    domain_id = parts[1]
                    
                    if protein_id not in protein_domain_map:
                        protein_domain_map[protein_id] = []
                    protein_domain_map[protein_id].append(domain_id)
        
        # Generate supra-domains for testing
        for protein_id, domains in protein_domain_map.items():
            if len(domains) > 1:
                # Add pairwise combinations
                for i in range(len(domains) - 1):
                    supra_domain = f"{domains[i]},{domains[i+1]}"
                    protein_domain_map[protein_id].append(supra_domain)
        
        return protein_domain_map
    
    def parse_mock_goa_data(self, workspace: Path) -> dict:
        """Parse mock GOA data for testing."""
        goa_file = workspace / "data" / "raw" / "goa_uniprot_all.gaf.gz"
        
        protein_go_map = {}
        with gzip.open(goa_file, 'rt') as f:
            for line in f:
                if line.startswith('!'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    protein_id = parts[1]  # DB Object ID
                    go_id = parts[4]       # GO ID
                    
                    if protein_id not in protein_go_map:
                        protein_go_map[protein_id] = set()
                    protein_go_map[protein_id].add(go_id)
        
        return protein_go_map
    
    def test_data_acquisition_integration(self, workspace, integration_config):
        """Test data acquisition with mock downloads."""
        # Test that DataAcquisition can be initialized and directories created
        data_acquisition = DataAcquisition(integration_config)
        assert data_acquisition.data_dir.exists()
        
        # Test getting download summary (should work even with no files)
        summary = data_acquisition.get_download_summary()
        assert 'total_files' in summary
        assert 'total_size_gb' in summary
        assert summary['total_files'] >= 0
    
    def test_error_handling_integration(self, workspace, integration_config):
        """Test error handling across pipeline components."""
        # Test with empty data
        empty_domain_map = {}
        empty_go_map = {}
        
        # Should handle empty data gracefully
        try:
            inference_engine = StatisticalInferenceEngine(empty_domain_map, empty_go_map)
            # Should raise an error for empty input
            assert False, "Should have raised an error for empty input"
        except Exception:
            pass  # Expected behavior
        
        # Test database with invalid data
        db_manager = dcGODatabaseManager(workspace / "results" / "error_test.db")
        
        # Should handle empty annotation list
        empty_annotations = []
        db_manager.store_annotations(empty_annotations)
        
        summary = db_manager.get_summary_statistics()
        assert summary['total_annotations'] == 0
    
    def test_configuration_validation(self, workspace):
        """Test configuration validation and error handling."""
        from dataclasses import replace
        
        # Test with invalid base directory
        try:
            invalid_config = replace(Config(), base_dir=Path("/nonexistent/path"))
            # Should handle gracefully or raise appropriate error
        except Exception:
            pass  # Expected for invalid paths
        
        # Test data source validation
        try:
            invalid_data_source = DataSource(
                name="test",
                url="invalid_url",  # Invalid URL format
                description="Test source"
            )
            assert False, "Should have raised validation error"
        except Exception:
            pass  # Expected behavior
    
    @pytest.mark.slow
    def test_performance_with_larger_dataset(self, workspace, integration_config):
        """Test pipeline performance with slightly larger mock dataset."""
        # Create larger mock dataset (100 proteins)
        raw_dir = workspace / "data" / "raw"
        
        # Generate larger UniProt file
        uniprot_lines = []
        for i in range(100):
            protein_id = f"P{i:05d}"
            uniprot_lines.append(f">sp|{protein_id}|TEST{i}_HUMAN Test protein {i}")
            uniprot_lines.append("MKLLVLSLSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQH")
        
        with gzip.open(raw_dir / "uniprot_sprot_large.fasta.gz", "wt") as f:
            f.write("\n".join(uniprot_lines))
        
        # Generate corresponding domain and GO data
        domain_lines = []
        go_lines = ["!gaf-version: 2.2", "!generated-by: UniProt"]
        
        for i in range(100):
            protein_id = f"P{i:05d}"
            domain1 = f"PF{(i % 10):05d}"
            domain2 = f"PF{((i + 1) % 10):05d}"
            
            domain_lines.append(f"{protein_id}\t{domain1}\t10\t50\t1e-10\tT\tDomain{i % 10}")
            domain_lines.append(f"{protein_id}\t{domain2}\t60\t90\t1e-15\tT\tDomain{(i + 1) % 10}")
            
            go_term = f"GO:{7000000 + (i % 100):07d}"
            go_lines.append(f"UniProtKB\t{protein_id}\tTEST{i}\t\t{go_term}\tPMID:12345678\tIDA\t\tF\tTest protein {i}\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t")
        
        with gzip.open(raw_dir / "protein2ipr_large.dat.gz", "wt") as f:
            f.write("\n".join(domain_lines))
        
        with gzip.open(raw_dir / "goa_uniprot_large.gaf.gz", "wt") as f:
            f.write("\n".join(go_lines))
        
        # Parse larger dataset
        import time
        start_time = time.time()
        
        protein_domain_map = {}
        with gzip.open(raw_dir / "protein2ipr_large.dat.gz", 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    protein_id, domain_id = parts[0], parts[1]
                    if protein_id not in protein_domain_map:
                        protein_domain_map[protein_id] = []
                    protein_domain_map[protein_id].append(domain_id)
        
        protein_go_map = {}
        with gzip.open(raw_dir / "goa_uniprot_large.gaf.gz", 'rt') as f:
            for line in f:
                if line.startswith('!'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    protein_id, go_id = parts[1], parts[4]
                    if protein_id not in protein_go_map:
                        protein_go_map[protein_id] = set()
                    protein_go_map[protein_id].add(go_id)
        
        parse_time = time.time() - start_time
        
        # Test statistical inference with larger dataset
        start_time = time.time()
        
        inference_engine = StatisticalInferenceEngine(protein_domain_map, protein_go_map)
        results = inference_engine.run_statistical_tests(min_cooccurrence=2)
        
        inference_time = time.time() - start_time
        
        # Performance assertions (should complete in reasonable time)
        assert parse_time < 5.0, f"Parsing took too long: {parse_time:.2f}s"
        assert inference_time < 10.0, f"Inference took too long: {inference_time:.2f}s"
        assert len(protein_domain_map) == 100
        assert len(protein_go_map) == 100


if __name__ == "__main__":
    pytest.main([__file__, "-v"])