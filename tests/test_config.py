"""
Test configuration module

Basic tests to ensure configuration is loaded correctly.
"""

from pathlib import Path
from config.settings import Config


class TestConfig:
    """Test configuration settings"""
    
    def test_config_initialization(self):
        """Test that Config class initializes correctly"""
        config = Config()
        
        # Test that basic attributes exist
        assert hasattr(config, 'FDR_THRESHOLD')
        assert hasattr(config, 'MIN_PROTEINS_PER_ASSOCIATION')
        assert hasattr(config, 'DATASOURCES')
        assert hasattr(config, 'BASE_DIR')
        assert hasattr(config, 'DATA_DIR')
        
    def test_fdr_threshold(self):
        """Test FDR threshold is reasonable"""
        config = Config()
        assert 0 < config.FDR_THRESHOLD < 1
        
    def test_min_proteins_per_association(self):
        """Test minimum proteins per association is positive"""
        config = Config()
        assert config.MIN_PROTEINS_PER_ASSOCIATION > 0
        
    def test_data_sources_exist(self):
        """Test that data sources are defined"""
        config = Config()
        expected_sources = [
            'uniprot_sprot',
            'uniprot_trembl',
            'goa_annotations',
            'go_ontology',
            'pfam_hmms',
            'interpro_scan'
        ]
        
        for source in expected_sources:
            assert source in config.DATASOURCES
            assert isinstance(config.DATASOURCES[source], str)
            assert len(config.DATASOURCES[source]) > 0
            
    def test_paths_are_path_objects(self):
        """Test that path attributes are Path objects"""
        config = Config()
        
        assert isinstance(config.BASE_DIR, Path)
        assert isinstance(config.DATA_DIR, Path)
        assert isinstance(config.RESULTS_DIR, Path)
        assert isinstance(config.LOGS_DIR, Path)
        assert isinstance(config.DATABASE_PATH, Path)
        
    def test_num_cores_positive(self):
        """Test that NUM_CORES is positive"""
        config = Config()
        assert config.NUM_CORES > 0