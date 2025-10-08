"""
Configuration settings for dcGO Pipeline

This module contains all configuration settings and parameters for the dcGO pipeline,
including data sources, processing parameters, and file paths using Python 3.12 features.
"""

import os
import psutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Self, Union
from urllib.parse import urlparse
import logging

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Exception raised for configuration-related errors."""
    pass


@dataclass(frozen=True)
class DataSource:
    """Configuration for a data source with validation."""
    name: str
    url: str
    description: str
    required: bool = True
    checksum: Optional[str] = None
    
    def __post_init__(self) -> None:
        """Validate data source configuration."""
        if not self.url:
            raise ConfigurationError(f"URL cannot be empty for data source '{self.name}'")
        
        # Basic URL validation
        parsed = urlparse(self.url)
        if not parsed.scheme or not parsed.netloc:
            raise ConfigurationError(f"Invalid URL format for data source '{self.name}': {self.url}")
        
        if parsed.scheme not in {'http', 'https', 'ftp'}:
            raise ConfigurationError(f"Unsupported URL scheme for data source '{self.name}': {parsed.scheme}")


@dataclass(frozen=True)
class ProcessingParameters:
    """Statistical and processing parameters with validation."""
    fdr_threshold: float = 0.01
    min_proteins_per_association: int = 3
    max_supra_domain_length: int = 3
    alpha_threshold: float = 0.05
    min_cooccurrence_threshold: int = 5
    max_sequence_length: int = 50000
    min_sequence_length: int = 20
    chunk_size: int = 65536
    timeout: int = 30
    
    def __post_init__(self) -> None:
        """Validate processing parameters."""
        if not (0 < self.fdr_threshold < 1):
            raise ConfigurationError(f"FDR threshold must be between 0 and 1, got {self.fdr_threshold}")
        
        if self.min_proteins_per_association <= 0:
            raise ConfigurationError(f"Minimum proteins per association must be positive, got {self.min_proteins_per_association}")
        
        if self.max_supra_domain_length <= 0:
            raise ConfigurationError(f"Maximum supra-domain length must be positive, got {self.max_supra_domain_length}")
        
        if not (0 < self.alpha_threshold < 1):
            raise ConfigurationError(f"Alpha threshold must be between 0 and 1, got {self.alpha_threshold}")
        
        if self.min_cooccurrence_threshold <= 0:
            raise ConfigurationError(f"Minimum cooccurrence threshold must be positive, got {self.min_cooccurrence_threshold}")
        
        if self.max_sequence_length <= self.min_sequence_length:
            raise ConfigurationError(f"Maximum sequence length ({self.max_sequence_length}) must be greater than minimum ({self.min_sequence_length})")
        
        if self.chunk_size <= 0:
            raise ConfigurationError(f"Chunk size must be positive, got {self.chunk_size}")
        
        if self.timeout <= 0:
            raise ConfigurationError(f"Timeout must be positive, got {self.timeout}")


@dataclass(frozen=True)
class ComputeResources:
    """Compute resource configuration with auto-detection and validation."""
    num_cores: int = field(default_factory=lambda: psutil.cpu_count())
    memory_limit_gb: Optional[float] = field(default_factory=lambda: psutil.virtual_memory().total / (1024**3))
    java_heap_size: str = "4G"
    temp_dir: Optional[Path] = None
    
    def __post_init__(self) -> None:
        """Validate compute resource settings."""
        available_cores = psutil.cpu_count()
        if self.num_cores <= 0 or self.num_cores > available_cores:
            raise ConfigurationError(f"Number of cores must be between 1 and {available_cores}, got {self.num_cores}")
        
        if self.memory_limit_gb is not None:
            available_memory = psutil.virtual_memory().total / (1024**3)
            if self.memory_limit_gb <= 0 or self.memory_limit_gb > available_memory:
                raise ConfigurationError(f"Memory limit must be between 0 and {available_memory:.1f}GB, got {self.memory_limit_gb}")
        
        # Validate Java heap size format
        if not self.java_heap_size.endswith(('M', 'G', 'm', 'g')):
            raise ConfigurationError(f"Java heap size must end with 'M' or 'G', got {self.java_heap_size}")
        
        try:
            int(self.java_heap_size[:-1])
        except ValueError:
            raise ConfigurationError(f"Invalid Java heap size format: {self.java_heap_size}")


@dataclass
class Config:
    """
    Main configuration class for dcGO Pipeline using Python 3.12 dataclass features.
    
    This class provides comprehensive configuration management with:
    - Type-safe configuration parameters
    - Automatic validation
    - Environment variable override support
    - Path management with automatic directory creation
    - Resource auto-detection
    """
    
    # Core paths (calculated from current file location)
    base_dir: Path = field(default_factory=lambda: Path(__file__).parent.parent.resolve())
    
    # Data sources configuration
    data_sources: Dict[str, DataSource] = field(default_factory=lambda: {
        'uniprot_sprot': DataSource(
            name='uniprot_sprot',
            url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
            description='UniProt Swiss-Prot protein sequences',
            required=True
        ),
        'uniprot_trembl': DataSource(
            name='uniprot_trembl', 
            url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz',
            description='UniProt TrEMBL protein sequences',
            required=False
        ),
        'goa_annotations': DataSource(
            name='goa_annotations',
            url='ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz',
            description='Gene Ontology Annotation (GOA) database',
            required=True
        ),
        'go_ontology': DataSource(
            name='go_ontology',
            url='http://current.geneontology.org/ontology/go-basic.obo',
            description='Gene Ontology basic ontology file',
            required=True
        ),
        'pfam_hmms': DataSource(
            name='pfam_hmms',
            url='https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz',
            description='Pfam HMM profiles for domain detection',
            required=False
        ),
        'interpro_scan': DataSource(
            name='interpro_scan',
            url='https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz',
            description='InterProScan domain analysis software',
            required=True
        )
    })
    
    # Processing parameters
    processing: ProcessingParameters = field(default_factory=ProcessingParameters)
    
    # Compute resources
    compute: ComputeResources = field(default_factory=ComputeResources)
    
    # Environment overrides
    use_env_overrides: bool = True
    
    def __post_init__(self) -> None:
        """Initialize configuration with validation and environment overrides."""
        # Apply environment variable overrides if enabled
        if self.use_env_overrides:
            self._apply_env_overrides()
        
        # Create directory structure
        self._create_directories()
        
        # Log configuration summary
        self._log_configuration()
    
    @property
    def data_dir(self) -> Path:
        """Directory for storing raw and processed data files."""
        return self.base_dir / "data"
    
    @property 
    def results_dir(self) -> Path:
        """Directory for storing analysis results."""
        return self.base_dir / "results"
    
    @property
    def logs_dir(self) -> Path:
        """Directory for storing log files."""
        return self.base_dir / "logs"
    
    @property
    def cache_dir(self) -> Path:
        """Directory for storing cached intermediate results."""
        return self.data_dir / "cache"
    
    @property
    def database_path(self) -> Path:
        """Path to the main SQLite database file."""
        return self.results_dir / "dcgo_database.db"
    
    @property
    def temp_dir(self) -> Path:
        """Temporary directory for processing operations."""
        if self.compute.temp_dir:
            return self.compute.temp_dir
        return self.base_dir / "temp"
    
    def _apply_env_overrides(self) -> None:
        """Apply environment variable overrides to configuration."""
        # Processing parameter overrides
        if fdr_env := os.getenv('DCGO_FDR_THRESHOLD'):
            try:
                fdr_value = float(fdr_env)
                object.__setattr__(self.processing, 'fdr_threshold', fdr_value)
            except ValueError:
                logger.warning(f"Invalid FDR threshold in environment: {fdr_env}")
        
        if cores_env := os.getenv('DCGO_NUM_CORES'):
            try:
                cores_value = int(cores_env)
                object.__setattr__(self.compute, 'num_cores', cores_value)
            except ValueError:
                logger.warning(f"Invalid number of cores in environment: {cores_env}")
        
        if min_proteins_env := os.getenv('DCGO_MIN_PROTEINS'):
            try:
                min_proteins_value = int(min_proteins_env)
                object.__setattr__(self.processing, 'min_proteins_per_association', min_proteins_value)
            except ValueError:
                logger.warning(f"Invalid minimum proteins in environment: {min_proteins_env}")
        
        if java_heap_env := os.getenv('DCGO_JAVA_HEAP'):
            object.__setattr__(self.compute, 'java_heap_size', java_heap_env)
    
    def _create_directories(self) -> None:
        """Create necessary directories if they don't exist."""
        directories = [
            self.data_dir / "raw",
            self.data_dir / "processed", 
            self.data_dir / "interim",
            self.cache_dir,
            self.results_dir,
            self.logs_dir,
            self.temp_dir
        ]
        
        for directory in directories:
            try:
                directory.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Ensured directory exists: {directory}")
            except PermissionError as e:
                raise ConfigurationError(f"Cannot create directory {directory}: {e}")
    
    def _log_configuration(self) -> None:
        """Log configuration summary for debugging."""
        logger.info("dcGO Pipeline Configuration Summary:")
        logger.info(f"  Base directory: {self.base_dir}")
        logger.info(f"  Data directory: {self.data_dir}")
        logger.info(f"  Results directory: {self.results_dir}")
        logger.info(f"  CPU cores: {self.compute.num_cores}")
        logger.info(f"  Memory limit: {self.compute.memory_limit_gb:.1f}GB")
        logger.info(f"  FDR threshold: {self.processing.fdr_threshold}")
        logger.info(f"  Required data sources: {len([ds for ds in self.data_sources.values() if ds.required])}")
    
    @classmethod
    def from_dict(cls, config_dict: Dict) -> Self:
        """Create configuration from dictionary."""
        return cls(**config_dict)
    
    @classmethod
    def from_env(cls) -> Self:
        """Create configuration with environment variable overrides enabled."""
        return cls(use_env_overrides=True)
    
    def get_data_source_url(self, source_name: str) -> str:
        """Get URL for a specific data source."""
        if source_name not in self.data_sources:
            raise ConfigurationError(f"Unknown data source: {source_name}")
        return self.data_sources[source_name].url
    
    def get_required_data_sources(self) -> List[str]:
        """Get list of required data source names."""
        return [name for name, source in self.data_sources.items() if source.required]
    
    def validate_paths(self) -> bool:
        """Validate all configured paths are accessible."""
        paths_to_check = [
            self.base_dir,
            self.data_dir,
            self.results_dir,
            self.logs_dir
        ]
        
        for path in paths_to_check:
            if not path.exists():
                logger.error(f"Required path does not exist: {path}")
                return False
            
            if not os.access(path, os.R_OK | os.W_OK):
                logger.error(f"Insufficient permissions for path: {path}")
                return False
        
        return True
    
    def get_system_info(self) -> Dict[str, Union[str, int, float]]:
        """Get system information for diagnostics."""
        return {
            'python_version': os.sys.version,
            'cpu_cores_available': psutil.cpu_count(),
            'cpu_cores_configured': self.compute.num_cores,
            'memory_total_gb': psutil.virtual_memory().total / (1024**3),
            'memory_configured_gb': self.compute.memory_limit_gb,
            'disk_free_gb': psutil.disk_usage(self.base_dir).free / (1024**3),
            'base_directory': str(self.base_dir)
        }
    
    # Legacy compatibility properties for existing code
    @property
    def DATASOURCES(self) -> Dict[str, str]:
        """Legacy compatibility: return data sources as URL dictionary."""
        return {name: source.url for name, source in self.data_sources.items()}
    
    @property
    def FDR_THRESHOLD(self) -> float:
        """Legacy compatibility: FDR threshold."""
        return self.processing.fdr_threshold
    
    @property
    def MIN_PROTEINS_PER_ASSOCIATION(self) -> int:
        """Legacy compatibility: minimum proteins per association."""
        return self.processing.min_proteins_per_association
    
    @property
    def MAX_SUPRA_DOMAIN_LENGTH(self) -> int:
        """Legacy compatibility: maximum supra-domain length."""
        return self.processing.max_supra_domain_length
    
    @property
    def NUM_CORES(self) -> int:
        """Legacy compatibility: number of CPU cores."""
        return self.compute.num_cores
    
    @property
    def BASE_DIR(self) -> Path:
        """Legacy compatibility: base directory."""
        return self.base_dir
    
    @property
    def DATA_DIR(self) -> Path:
        """Legacy compatibility: data directory."""
        return self.data_dir
    
    @property
    def RESULTS_DIR(self) -> Path:
        """Legacy compatibility: results directory."""
        return self.results_dir
    
    @property
    def LOGS_DIR(self) -> Path:
        """Legacy compatibility: logs directory."""
        return self.logs_dir
    
    @property
    def DATABASE_PATH(self) -> Path:
        """Legacy compatibility: database path."""
        return self.database_path


# Global configuration instance for backward compatibility
config = Config()