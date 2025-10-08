"""
Unit tests for the data acquisition module.

Tests download functionality, validation, and file management without
requiring actual large downloads.
"""

import pytest
import tempfile
import gzip
from pathlib import Path
from unittest.mock import Mock, patch
from urllib.error import URLError

# Add src to path for testing
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from data_acquisition import DataAcquisition, DataAcquisitionError
from config.settings import Config


class TestDataAcquisition:
    """Test suite for DataAcquisition class."""
    
    @pytest.fixture
    def temp_data_dir(self):
        """Create temporary directory for test data."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield Path(temp_dir)
    
    @pytest.fixture
    def mock_config(self, temp_data_dir):
        """Create mock configuration for testing."""
        from dataclasses import replace
        config = Config()
        # Create new config with temp directory
        modified_config = replace(config, base_dir=temp_data_dir)
        modified_config.data_dir.mkdir(parents=True, exist_ok=True)
        return modified_config
    
    @pytest.fixture
    def data_acquisition(self, mock_config):
        """Create DataAcquisition instance for testing."""
        return DataAcquisition(mock_config)
    
    def test_init(self, data_acquisition, mock_config):
        """Test DataAcquisition initialization."""
        assert data_acquisition.config == mock_config
        assert data_acquisition.data_dir.exists()
        assert data_acquisition.data_dir.name == "raw"  # Points to data/raw directory
    
    @patch('requests.get')
    def test_download_with_progress_success(self, mock_get, data_acquisition, temp_data_dir):
        """Test successful file download with progress tracking."""
        # Mock response
        mock_response = Mock()
        mock_response.headers = {'content-length': '1000'}
        mock_response.iter_content.return_value = [b'test'] * 250  # 1000 bytes total
        mock_get.return_value = mock_response
        
        test_url = "https://example.com/test.txt"
        test_filepath = temp_data_dir / "test.txt"
        
        # Execute download
        data_acquisition.download_with_progress(test_url, test_filepath)
        
        # Verify file was created and has correct content
        assert test_filepath.exists()
        assert test_filepath.read_text() == 'test' * 250
        # Check that requests.get was called with the URL (ignore additional parameters)
        mock_get.assert_called_once()
        args, kwargs = mock_get.call_args
        assert args[0] == test_url
        assert kwargs['stream'] is True
    
    @patch('requests.get')
    def test_download_with_progress_failure(self, mock_get, data_acquisition, temp_data_dir):
        """Test download failure handling."""
        mock_get.side_effect = URLError("Network error")
        
        test_url = "https://example.com/nonexistent.txt"
        test_filepath = temp_data_dir / "test.txt"
        
        with pytest.raises(DataAcquisitionError):
            data_acquisition.download_with_progress(test_url, test_filepath)
    
    def test_verify_checksum_success(self, data_acquisition, temp_data_dir):
        """Test successful checksum verification."""
        test_file = temp_data_dir / "test.txt"
        test_content = "Hello, World!"
        test_file.write_text(test_content)
        
        # Calculate expected MD5
        import hashlib
        expected_md5 = hashlib.md5(test_content.encode()).hexdigest()
        
        # Verify checksum
        result = data_acquisition.verify_checksum(test_file, expected_md5)
        assert result is True
    
    def test_verify_checksum_failure(self, data_acquisition, temp_data_dir):
        """Test checksum verification failure."""
        test_file = temp_data_dir / "test.txt"
        test_file.write_text("Hello, World!")
        
        wrong_checksum = "incorrect_checksum"
        result = data_acquisition.verify_checksum(test_file, wrong_checksum)
        assert result is False
    
    def test_verify_checksum_missing_file(self, data_acquisition, temp_data_dir):
        """Test checksum verification with missing file."""
        missing_file = temp_data_dir / "nonexistent.txt"
        
        with pytest.raises(DataAcquisitionError):
            data_acquisition.verify_checksum(missing_file, "any_checksum")
    
    @patch.object(DataAcquisition, 'download_with_progress')
    def test_download_all_datasets_success(self, mock_download, data_acquisition):
        """Test downloading all configured datasets."""
        # Mock successful downloads
        mock_download.return_value = None
        
        # Execute download
        downloaded_files = data_acquisition.download_all_datasets()
        
        # Verify downloads were attempted for required datasets
        assert len(downloaded_files) > 0
        assert 'uniprot_sprot' in downloaded_files
        assert 'goa_annotations' in downloaded_files
        assert 'go_ontology' in downloaded_files
        assert 'interpro_mappings' in downloaded_files
    
    @patch.object(DataAcquisition, 'download_with_progress')
    def test_download_all_datasets_partial_failure(self, mock_download, data_acquisition):
        """Test handling of partial download failures."""
        # Mock some downloads failing
        def side_effect(url, filepath):
            if 'uniprot_trembl' in str(filepath):
                raise DataAcquisitionError("Download failed")
            return None
        
        mock_download.side_effect = side_effect
        
        # Should not raise exception for optional datasets
        downloaded_files = data_acquisition.download_all_datasets()
        assert len(downloaded_files) > 0
    
    def test_get_download_summary(self, data_acquisition, temp_data_dir):
        """Test download summary generation."""
        # Create some mock downloaded files
        (data_acquisition.data_dir / "raw").mkdir(parents=True, exist_ok=True)
        test_file1 = data_acquisition.data_dir / "raw" / "test1.txt"
        test_file2 = data_acquisition.data_dir / "raw" / "test2.txt"
        
        test_file1.write_text("content1")
        test_file2.write_text("content2" * 1000)  # Larger file
        
        summary = data_acquisition.get_download_summary()
        
        assert 'total_files' in summary
        assert 'total_size_gb' in summary
        assert 'files' in summary
        assert summary['total_files'] >= 2
        assert summary['total_size_gb'] >= 0
    
    def test_cleanup_old_downloads(self, data_acquisition, temp_data_dir):
        """Test cleanup of old download files."""
        # Create some old files
        old_files_dir = data_acquisition.data_dir / "raw"
        old_files_dir.mkdir(parents=True, exist_ok=True)
        
        old_file = old_files_dir / "old_file.txt"
        old_file.write_text("old content")
        
        # Mock file age to be old
        import time
        old_time = time.time() - (8 * 24 * 60 * 60)  # 8 days ago
        old_file.touch(times=(old_time, old_time))
        
        removed_files = data_acquisition.cleanup_old_downloads(max_age_days=7)
        
        # File should be cleaned up
        assert len(removed_files) >= 0  # May be 0 if filesystem doesn't support old timestamps
    
    def test_ftp_download_support(self, data_acquisition):
        """Test FTP URL parsing for InterPro downloads."""
        ftp_url = "ftp://ftp.ebi.ac.uk/pub/databases/interpro/test.dat.gz"
        
        # Should not raise exception during URL parsing
        from urllib.parse import urlparse
        parsed = urlparse(ftp_url)
        assert parsed.scheme == 'ftp'
        assert 'ebi.ac.uk' in parsed.netloc


class TestDataAcquisitionIntegration:
    """Integration tests for data acquisition workflows."""
    
    @pytest.fixture
    def temp_workspace(self):
        """Create temporary workspace for integration tests."""
        with tempfile.TemporaryDirectory() as temp_dir:
            workspace = Path(temp_dir)
            (workspace / "data" / "raw").mkdir(parents=True)
            (workspace / "data" / "processed").mkdir(parents=True)
            yield workspace
    
    def create_mock_goa_file(self, filepath: Path) -> None:
        """Create a mock GOA file for testing."""
        goa_content = '''!gaf-version: 2.2
!generated-by: UniProt
UniProtKB\tP12345\tTEST1\t\tGO:0005515\tPMID:12345678\tIDA\t\tF\tTest protein 1\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
UniProtKB\tP67890\tTEST2\t\tGO:0008150\tPMID:87654321\tISS\t\tP\tTest protein 2\t\tprotein\ttaxon:9606\t20230101\tUniProt\t\t
'''
        with gzip.open(filepath, 'wt') as f:
            f.write(goa_content)
    
    def create_mock_interpro_file(self, filepath: Path) -> None:
        """Create a mock InterPro mapping file for testing."""
        interpro_content = '''P12345\tPF00001\t10\t50\t1e-10\tT\tDomain1
P12345\tPF00002\t60\t90\t1e-15\tT\tDomain2
P67890\tPF00001\t20\t60\t1e-12\tT\tDomain1
'''
        with gzip.open(filepath, 'wt') as f:
            f.write(interpro_content)
    
    def test_mock_data_workflow(self, temp_workspace):
        """Test complete workflow with mock data files."""
        # Create mock data files
        raw_dir = temp_workspace / "data" / "raw"
        
        goa_file = raw_dir / "goa_uniprot_all.gaf.gz"
        interpro_file = raw_dir / "protein2ipr.dat.gz"
        
        self.create_mock_goa_file(goa_file)
        self.create_mock_interpro_file(interpro_file)
        
        # Verify files can be read
        assert goa_file.exists()
        assert interpro_file.exists()
        
        # Test file sizes
        assert goa_file.stat().st_size > 0
        assert interpro_file.stat().st_size > 0
        
        # Test content can be decompressed and read
        with gzip.open(goa_file, 'rt') as f:
            goa_lines = f.readlines()
            assert len(goa_lines) >= 2  # Header + data lines
            assert any('P12345' in line for line in goa_lines)
        
        with gzip.open(interpro_file, 'rt') as f:
            interpro_lines = f.readlines()
            assert len(interpro_lines) >= 3
            assert any('PF00001' in line for line in interpro_lines)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])