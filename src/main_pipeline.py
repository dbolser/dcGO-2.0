"""
Main pipeline orchestrator for the dcGO (Domain-centric Gene Ontology) pipeline.

This module provides the dcGOPipeline class that coordinates all pipeline modules
and a comprehensive CLI interface for running the complete analysis workflow.
Includes checkpoint/resume functionality, progress tracking, and robust error handling.

Author: dcGO Pipeline Team
License: MIT
Python: 3.12+
"""

import asyncio
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum, auto
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Union

import click
from loguru import logger

from config.settings import config
from src.data_acquisition import DataAcquisition, DataAcquisitionError
from src.database_manager import DatabaseError, dcGODatabaseManager
from src.domain_scanning import DomainArchitectureScanner
from src.ontology_processor import OntologyProcessor
from src.statistical_inference import StatisticalInferenceEngine


class PipelineStage(Enum):
    """Enumeration of pipeline stages for checkpoint/resume functionality."""
    
    INITIALIZATION = auto()
    DATA_ACQUISITION = auto()
    DOMAIN_SCANNING = auto()
    ONTOLOGY_PROCESSING = auto()
    STATISTICAL_INFERENCE = auto()
    RESULTS_EXPORT = auto()
    CLEANUP = auto()
    COMPLETED = auto()


@dataclass
class PipelineCheckpoint:
    """Represents a pipeline checkpoint for resume functionality."""
    
    stage: PipelineStage
    timestamp: datetime
    parameters: Dict[str, Any]
    completed_substeps: Set[str] = field(default_factory=set)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert checkpoint to dictionary for JSON serialization."""
        return {
            'stage': self.stage.name,
            'timestamp': self.timestamp.isoformat(),
            'parameters': self.parameters,
            'completed_substeps': list(self.completed_substeps),
            'metadata': self.metadata
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PipelineCheckpoint':
        """Create checkpoint from dictionary loaded from JSON."""
        return cls(
            stage=PipelineStage[data['stage']],
            timestamp=datetime.fromisoformat(data['timestamp']),
            parameters=data['parameters'],
            completed_substeps=set(data.get('completed_substeps', [])),
            metadata=data.get('metadata', {})
        )


@dataclass
class PipelineProgress:
    """Tracks overall pipeline progress and performance metrics."""
    
    current_stage: PipelineStage
    start_time: datetime
    stage_start_time: datetime
    total_stages: int = 7
    completed_stages: int = 0
    estimated_completion: Optional[datetime] = None
    stage_durations: Dict[PipelineStage, timedelta] = field(default_factory=dict)
    
    @property
    def overall_progress(self) -> float:
        """Calculate overall progress percentage (0-100)."""
        return (self.completed_stages / self.total_stages) * 100
    
    @property
    def elapsed_time(self) -> timedelta:
        """Total elapsed time since pipeline start."""
        return datetime.now() - self.start_time
    
    def update_stage(self, stage: PipelineStage) -> None:
        """Update progress to new stage."""
        if hasattr(self, '_prev_stage'):
            # Record duration of previous stage
            duration = datetime.now() - self.stage_start_time
            self.stage_durations[self._prev_stage] = duration
            self.completed_stages += 1
        
        self._prev_stage = self.current_stage
        self.current_stage = stage
        self.stage_start_time = datetime.now()
        
        # Estimate completion time based on average stage duration
        if self.stage_durations:
            avg_duration = sum(self.stage_durations.values(), timedelta()) / len(self.stage_durations)
            remaining_stages = self.total_stages - self.completed_stages
            self.estimated_completion = datetime.now() + (avg_duration * remaining_stages)


class PipelineError(Exception):
    """Base exception for pipeline-specific errors."""
    pass


class dcGOPipeline:
    """
    Main orchestrator for the dcGO (Domain-centric Gene Ontology) analysis pipeline.
    
    This class coordinates all pipeline modules including data acquisition, domain scanning,
    ontology processing, statistical inference, and results export. Features include:
    
    - Comprehensive logging and error handling
    - Progress tracking with time estimation
    - Checkpoint/resume functionality for long-running analyses
    - Async operations where beneficial
    - Modern Python 3.12 features including pattern matching and improved typing
    
    Attributes:
        config: Configuration object with all pipeline parameters
        progress: Progress tracker for the current run
        checkpoint_file: Path to checkpoint file for resume functionality
        max_workers: Maximum number of concurrent workers for parallel operations
    """
    
    def __init__(
        self,
        max_workers: int = 4,
        checkpoint_file: Optional[Path] = None,
        enable_async: bool = True
    ):
        """
        Initialize the dcGO pipeline orchestrator.
        
        Args:
            max_workers: Maximum concurrent workers for parallel operations
            checkpoint_file: Path for checkpoint file (default: logs/pipeline_checkpoint.json)
            enable_async: Whether to use async operations where beneficial
        """
        self.config = config
        self.max_workers = max_workers
        self.enable_async = enable_async
        self.checkpoint_file = checkpoint_file or (self.config.LOGS_DIR / "pipeline_checkpoint.json")
        
        # Initialize components
        self.data_acquisition: Optional[DataAcquisition] = None
        self.domain_scanner: Optional[DomainArchitectureScanner] = None
        self.ontology_processor: Optional[OntologyProcessor] = None
        self.statistical_engine: Optional[StatisticalInferenceEngine] = None
        self.database_manager: Optional[dcGODatabaseManager] = None
        
        # Progress tracking
        self.progress: Optional[PipelineProgress] = None
        self.current_checkpoint: Optional[PipelineCheckpoint] = None
        
        # Configure logging
        self._setup_logging()
        
    def _setup_logging(self) -> None:
        """Configure comprehensive logging for the pipeline."""
        # Remove default handler
        logger.remove()
        
        # Add console handler with color
        logger.add(
            sys.stderr,
            level="INFO",
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
                   "<level>{level: <8}</level> | "
                   "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> | "
                   "<level>{message}</level>",
            colorize=True
        )
        
        # Add file handler for detailed logs
        logger.add(
            self.config.LOGS_DIR / f"dcgo_pipeline_{datetime.now():%Y%m%d_%H%M%S}.log",
            level="DEBUG",
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} | {message}",
            rotation="100 MB",
            retention="30 days",
            compression="gzip"
        )
        
        logger.info("Pipeline logging configured successfully")
    
    def save_checkpoint(self, stage: PipelineStage, parameters: Dict[str, Any], **metadata: Any) -> None:
        """
        Save pipeline checkpoint for resume functionality.
        
        Args:
            stage: Current pipeline stage
            parameters: Pipeline parameters used
            **metadata: Additional metadata to store in checkpoint
        """
        try:
            checkpoint = PipelineCheckpoint(
                stage=stage,
                timestamp=datetime.now(),
                parameters=parameters,
                completed_substeps=set(),  # Can be populated by specific stages
                metadata=metadata
            )
            
            self.current_checkpoint = checkpoint
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint.to_dict(), f, indent=2)
                
            logger.info(f"Checkpoint saved at stage {stage.name}")
            
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")
    
    def load_checkpoint(self) -> Optional[PipelineCheckpoint]:
        """
        Load pipeline checkpoint for resume functionality.
        
        Returns:
            PipelineCheckpoint if found and valid, None otherwise
        """
        if not self.checkpoint_file.exists():
            logger.info("No checkpoint file found")
            return None
            
        try:
            with open(self.checkpoint_file, 'r') as f:
                data = json.load(f)
                
            checkpoint = PipelineCheckpoint.from_dict(data)
            logger.info(f"Loaded checkpoint from stage {checkpoint.stage.name} "
                       f"at {checkpoint.timestamp}")
            return checkpoint
            
        except Exception as e:
            logger.error(f"Failed to load checkpoint: {e}")
            return None
    
    def clear_checkpoint(self) -> None:
        """Remove checkpoint file after successful completion."""
        try:
            if self.checkpoint_file.exists():
                self.checkpoint_file.unlink()
                logger.info("Checkpoint file cleared")
        except Exception as e:
            logger.warning(f"Failed to clear checkpoint: {e}")
    
    async def initialize_components(self, parameters: Dict[str, Any]) -> None:
        """
        Initialize all pipeline components asynchronously.
        
        Args:
            parameters: Pipeline parameters from CLI or checkpoint
        """
        logger.info("Initializing pipeline components...")
        
        # Initialize components concurrently where possible
        async def init_data_acquisition():
            self.data_acquisition = DataAcquisition(
                data_dir=self.config.DATA_DIR / "raw",
                chunk_size=parameters.get('download_chunk_size', 65536),
                timeout=parameters.get('download_timeout', 30)
            )
            
        async def init_database():
            self.database_manager = dcGODatabaseManager(
                db_path=self.config.DATABASE_PATH
            )
            # Initialize database schema
            await asyncio.to_thread(self.database_manager.initialize_database)
            
        # Run initializations concurrently
        await asyncio.gather(
            init_data_acquisition(),
            init_database(),
            return_exceptions=True
        )
        
        # Initialize remaining components that depend on data
        self.domain_scanner = DomainArchitectureScanner(
            data_dir=self.config.DATA_DIR,
            num_cores=parameters.get('num_cores', self.config.NUM_CORES)
        )
        
        self.ontology_processor = OntologyProcessor(
            fdr_threshold=parameters.get('fdr_threshold', self.config.FDR_THRESHOLD)
        )
        
        self.statistical_engine = StatisticalInferenceEngine()
        
        logger.info("All components initialized successfully")
    
    async def run_data_acquisition(self, parameters: Dict[str, Any]) -> Dict[str, Path]:
        """
        Execute data acquisition stage asynchronously.
        
        Args:
            parameters: Pipeline parameters
            
        Returns:
            Dictionary mapping data source names to downloaded file paths
        """
        logger.info("Starting data acquisition stage...")
        
        # Check what datasets to download
        datasets_to_download = parameters.get('datasets', list(self.config.DATASOURCES.keys()))
        
        if parameters.get('skip_existing', True):
            logger.info("Skipping existing files (use --force-download to override)")
        
        try:
            # Download datasets with progress tracking
            if self.enable_async:
                # Use thread pool for I/O bound download operations
                with ThreadPoolExecutor(max_workers=min(len(datasets_to_download), 3)) as executor:
                    download_futures = []
                    
                    for dataset in datasets_to_download:
                        future = executor.submit(
                            self.data_acquisition.download_specific_dataset, 
                            dataset,
                            force=parameters.get('force_download', False)
                        )
                        download_futures.append((dataset, future))
                    
                    downloaded_files = {}
                    for dataset, future in download_futures:
                        try:
                            file_path = future.result(timeout=3600)  # 1 hour timeout per file
                            downloaded_files[dataset] = file_path
                            logger.info(f"Successfully downloaded {dataset}: {file_path}")
                        except Exception as e:
                            logger.error(f"Failed to download {dataset}: {e}")
                            raise PipelineError(f"Data acquisition failed for {dataset}") from e
            else:
                # Sequential download
                downloaded_files = {}
                for dataset in datasets_to_download:
                    try:
                        file_path = self.data_acquisition.download_specific_dataset(
                            dataset,
                            force=parameters.get('force_download', False)
                        )
                        downloaded_files[dataset] = file_path
                    except Exception as e:
                        logger.error(f"Failed to download {dataset}: {e}")
                        raise PipelineError(f"Data acquisition failed for {dataset}") from e
            
            logger.info(f"Data acquisition completed successfully. Downloaded {len(downloaded_files)} datasets.")
            return downloaded_files
            
        except DataAcquisitionError as e:
            logger.error(f"Data acquisition failed: {e}")
            raise PipelineError("Data acquisition stage failed") from e
    
    async def run_domain_scanning(self, parameters: Dict[str, Any], data_files: Dict[str, Path]) -> Dict[str, Any]:
        """
        Execute domain scanning stage.
        
        Args:
            parameters: Pipeline parameters
            data_files: Downloaded data files from previous stage
            
        Returns:
            Dictionary with domain scanning results and statistics
        """
        logger.info("Starting domain scanning stage...")
        
        try:
            # Get protein sequences file (UniProt Swiss-Prot preferably)
            sequence_file = data_files.get('uniprot_sprot')
            if not sequence_file or not sequence_file.exists():
                sequence_file = data_files.get('uniprot_trembl')
                if not sequence_file or not sequence_file.exists():
                    raise PipelineError("No UniProt sequence file available for domain scanning")
            
            logger.info(f"Using sequence file: {sequence_file}")
            
            # Configure scanning parameters
            scan_params = {
                'sequence_file': sequence_file,
                'output_dir': self.config.DATA_DIR / "processed" / "domains",
                'num_cores': parameters.get('num_cores', self.config.NUM_CORES),
                'batch_size': parameters.get('scan_batch_size', 1000),
                'evalue_threshold': parameters.get('evalue_threshold', 1e-3),
                'use_pfam': parameters.get('use_pfam', True),
                'use_interpro': parameters.get('use_interpro', True)
            }
            
            # Create output directory
            scan_params['output_dir'].mkdir(parents=True, exist_ok=True)
            
            # Run domain scanning with progress tracking
            scan_results = {}
            
            if scan_params['use_interpro']:
                logger.info("Running InterProScan analysis...")
                interpro_results = await asyncio.to_thread(
                    self.domain_scanner.run_interproscan_batch,
                    sequence_file,
                    scan_params['output_dir'] / "interpro_results.tsv",
                    batch_size=scan_params['batch_size'],
                    num_cores=scan_params['num_cores']
                )
                scan_results['interpro'] = interpro_results
            
            # Store results in database
            if self.database_manager:
                logger.info("Storing domain scan results in database...")
                await asyncio.to_thread(
                    self.database_manager.store_domain_results,
                    scan_results
                )
            
            logger.info("Domain scanning completed successfully")
            return scan_results
            
        except Exception as e:
            logger.error(f"Domain scanning failed: {e}")
            raise PipelineError("Domain scanning stage failed") from e
    
    async def run_ontology_processing(self, parameters: Dict[str, Any], data_files: Dict[str, Path]) -> Dict[str, Any]:
        """
        Execute ontology processing stage.
        
        Args:
            parameters: Pipeline parameters
            data_files: Downloaded data files from previous stage
            
        Returns:
            Dictionary with processed ontology data and annotations
        """
        logger.info("Starting ontology processing stage...")
        
        try:
            # Get GO ontology and annotation files
            go_ontology_file = data_files.get('go_ontology')
            goa_file = data_files.get('goa_annotations')
            
            if not go_ontology_file or not go_ontology_file.exists():
                raise PipelineError("GO ontology file not available")
            if not goa_file or not goa_file.exists():
                raise PipelineError("GOA annotation file not available")
            
            logger.info(f"Using GO ontology: {go_ontology_file}")
            logger.info(f"Using GOA annotations: {goa_file}")
            
            # Load and process ontology
            await asyncio.to_thread(
                self.ontology_processor.load_ontology,
                go_ontology_file
            )
            
            # Process annotations with filtering
            processing_params = {
                'min_proteins': parameters.get('min_proteins_per_association', 
                                              self.config.MIN_PROTEINS_PER_ASSOCIATION),
                'fdr_threshold': parameters.get('fdr_threshold', self.config.FDR_THRESHOLD),
                'evidence_codes': parameters.get('evidence_codes', None),  # None = use all
                'aspects': parameters.get('go_aspects', ['P', 'F', 'C'])  # Process, Function, Component
            }
            
            # Process annotations
            processed_annotations = await asyncio.to_thread(
                self.ontology_processor.process_annotations,
                goa_file,
                **processing_params
            )
            
            ontology_results = {
                'ontology_graph': self.ontology_processor.go_graph,
                'processed_annotations': processed_annotations,
                'statistics': {
                    'total_terms': len(self.ontology_processor.go_graph.nodes()),
                    'total_annotations': len(processed_annotations),
                    'aspects_processed': processing_params['aspects']
                }
            }
            
            logger.info(f"Ontology processing completed. "
                       f"Processed {ontology_results['statistics']['total_annotations']} annotations "
                       f"across {ontology_results['statistics']['total_terms']} GO terms")
            
            return ontology_results
            
        except Exception as e:
            logger.error(f"Ontology processing failed: {e}")
            raise PipelineError("Ontology processing stage failed") from e
    
    async def run_statistical_inference(
        self, 
        parameters: Dict[str, Any],
        domain_results: Dict[str, Any],
        ontology_results: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Execute statistical inference stage.
        
        Args:
            parameters: Pipeline parameters
            domain_results: Results from domain scanning stage
            ontology_results: Results from ontology processing stage
            
        Returns:
            Dictionary with statistical inference results and associations
        """
        logger.info("Starting statistical inference stage...")
        
        try:
            # Configure statistical parameters
            stats_params = {
                'fdr_threshold': parameters.get('fdr_threshold', self.config.FDR_THRESHOLD),
                'min_proteins': parameters.get('min_proteins_per_association', 
                                              self.config.MIN_PROTEINS_PER_ASSOCIATION),
                'max_supra_domain_length': parameters.get('max_supra_domain_length',
                                                         self.config.MAX_SUPRA_DOMAIN_LENGTH),
                'correction_method': parameters.get('multiple_testing_correction', 'benjamini_hochberg')
            }
            
            logger.info(f"Running statistical tests with FDR threshold: {stats_params['fdr_threshold']}")
            
            # Prepare data for statistical analysis
            protein_domain_map = domain_results.get('protein_domains', {})
            protein_go_map = {}
            
            # Build protein-GO mapping from processed annotations
            for annotation in ontology_results['processed_annotations']:
                protein_id = annotation.protein_id
                go_term = annotation.go_term
                
                if protein_id not in protein_go_map:
                    protein_go_map[protein_id] = set()
                protein_go_map[protein_id].add(go_term)
            
            # Load data into statistical engine
            self.statistical_engine.load_protein_domain_map(protein_domain_map)
            self.statistical_engine.load_protein_go_map(protein_go_map)
            
            # Run statistical tests
            statistical_results = await asyncio.to_thread(
                self.statistical_engine.run_statistical_tests,
                fdr_threshold=stats_params['fdr_threshold'],
                min_proteins=stats_params['min_proteins']
            )
            
            # Process results
            significant_associations = [
                result for result in statistical_results 
                if result.fdr_corrected_p <= stats_params['fdr_threshold']
            ]
            
            inference_results = {
                'all_results': statistical_results,
                'significant_associations': significant_associations,
                'statistics': {
                    'total_tests': len(statistical_results),
                    'significant_associations': len(significant_associations),
                    'fdr_threshold': stats_params['fdr_threshold'],
                    'correction_method': stats_params['correction_method']
                },
                'parameters': stats_params
            }
            
            logger.info(f"Statistical inference completed. "
                       f"Found {len(significant_associations)} significant associations "
                       f"out of {len(statistical_results)} total tests")
            
            return inference_results
            
        except Exception as e:
            logger.error(f"Statistical inference failed: {e}")
            raise PipelineError("Statistical inference stage failed") from e
    
    async def export_results(
        self, 
        parameters: Dict[str, Any],
        inference_results: Dict[str, Any],
        ontology_results: Dict[str, Any]
    ) -> Dict[str, Path]:
        """
        Export pipeline results to various formats.
        
        Args:
            parameters: Pipeline parameters
            inference_results: Results from statistical inference
            ontology_results: Results from ontology processing
            
        Returns:
            Dictionary mapping result types to output file paths
        """
        logger.info("Starting results export stage...")
        
        try:
            # Create results directory
            results_dir = self.config.RESULTS_DIR
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = results_dir / f"dcgo_run_{timestamp}"
            run_dir.mkdir(parents=True, exist_ok=True)
            
            output_files = {}
            
            # Export significant associations
            significant_file = run_dir / "significant_associations.tsv"
            await asyncio.to_thread(
                self._export_associations_tsv,
                inference_results['significant_associations'],
                significant_file
            )
            output_files['significant_associations'] = significant_file
            
            # Export all results if requested
            if parameters.get('export_all_results', False):
                all_results_file = run_dir / "all_results.tsv"
                await asyncio.to_thread(
                    self._export_associations_tsv,
                    inference_results['all_results'],
                    all_results_file
                )
                output_files['all_results'] = all_results_file
            
            # Export summary statistics
            summary_file = run_dir / "pipeline_summary.json"
            summary_data = {
                'run_timestamp': timestamp,
                'parameters': parameters,
                'statistics': {
                    'ontology': ontology_results.get('statistics', {}),
                    'inference': inference_results.get('statistics', {})
                },
                'output_files': {k: str(v) for k, v in output_files.items()}
            }
            
            with open(summary_file, 'w') as f:
                json.dump(summary_data, f, indent=2)
            output_files['summary'] = summary_file
            
            # Store results in database if available
            if self.database_manager:
                logger.info("Storing results in database...")
                await asyncio.to_thread(
                    self.database_manager.store_associations,
                    inference_results['significant_associations']
                )
            
            logger.info(f"Results exported successfully to {run_dir}")
            logger.info(f"Generated {len(output_files)} output files")
            
            return output_files
            
        except Exception as e:
            logger.error(f"Results export failed: {e}")
            raise PipelineError("Results export stage failed") from e
    
    def _export_associations_tsv(self, associations: List[Any], output_file: Path) -> None:
        """Export associations to TSV format."""
        import csv
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            writer.writerow([
                'domain_id', 'go_term', 'p_value', 'fdr_corrected_p',
                'association_score', 'num_proteins', 'description'
            ])
            
            # Write associations
            for assoc in associations:
                writer.writerow([
                    assoc.domain_id,
                    assoc.go_term,
                    assoc.p_value,
                    assoc.fdr_corrected_p,
                    assoc.association_score,
                    assoc.num_proteins,
                    getattr(assoc, 'description', '')
                ])
    
    @asynccontextmanager
    async def _progress_context(self, stage: PipelineStage):
        """Context manager for tracking stage progress."""
        if self.progress:
            self.progress.update_stage(stage)
            logger.info(f"Entering stage: {stage.name} "
                       f"({self.progress.overall_progress:.1f}% complete)")
        
        try:
            yield
        finally:
            if self.progress:
                logger.info(f"Completed stage: {stage.name}")
    
    async def run_full_pipeline(
        self,
        resume: bool = False,
        force_download: bool = False,
        datasets: Optional[List[str]] = None,
        num_cores: Optional[int] = None,
        fdr_threshold: Optional[float] = None,
        export_all_results: bool = False,
        **additional_params: Any
    ) -> Dict[str, Any]:
        """
        Execute the complete dcGO analysis pipeline.
        
        Args:
            resume: Whether to resume from checkpoint if available
            force_download: Force re-download of existing files
            datasets: Specific datasets to download (default: all)
            num_cores: Number of cores for parallel processing
            fdr_threshold: FDR threshold for significance
            export_all_results: Export all results, not just significant ones
            **additional_params: Additional pipeline parameters
            
        Returns:
            Dictionary containing all pipeline results and output paths
            
        Raises:
            PipelineError: If any stage fails
        """
        start_time = datetime.now()
        logger.info(f"Starting dcGO pipeline run at {start_time}")
        
        # Prepare parameters
        parameters = {
            'force_download': force_download,
            'datasets': datasets or list(self.config.DATASOURCES.keys()),
            'num_cores': num_cores or self.config.NUM_CORES,
            'fdr_threshold': fdr_threshold or self.config.FDR_THRESHOLD,
            'export_all_results': export_all_results,
            'skip_existing': not force_download,
            **additional_params
        }
        
        # Initialize progress tracking
        self.progress = PipelineProgress(
            current_stage=PipelineStage.INITIALIZATION,
            start_time=start_time,
            stage_start_time=start_time
        )
        
        # Check for resume
        start_stage = PipelineStage.INITIALIZATION
        loaded_checkpoint = None
        
        if resume:
            loaded_checkpoint = self.load_checkpoint()
            if loaded_checkpoint:
                start_stage = loaded_checkpoint.stage
                # Merge checkpoint parameters with current parameters
                parameters.update(loaded_checkpoint.parameters)
                logger.info(f"Resuming from stage: {start_stage.name}")
        
        try:
            results = {}
            
            # Stage 1: Initialization
            if start_stage.value <= PipelineStage.INITIALIZATION.value:
                async with self._progress_context(PipelineStage.INITIALIZATION):
                    await self.initialize_components(parameters)
                    self.save_checkpoint(PipelineStage.DATA_ACQUISITION, parameters)
            
            # Stage 2: Data Acquisition
            if start_stage.value <= PipelineStage.DATA_ACQUISITION.value:
                async with self._progress_context(PipelineStage.DATA_ACQUISITION):
                    data_files = await self.run_data_acquisition(parameters)
                    results['data_files'] = data_files
                    self.save_checkpoint(PipelineStage.DOMAIN_SCANNING, parameters, 
                                       data_files={str(k): str(v) for k, v in data_files.items()})
            
            # Stage 3: Domain Scanning
            if start_stage.value <= PipelineStage.DOMAIN_SCANNING.value:
                async with self._progress_context(PipelineStage.DOMAIN_SCANNING):
                    # Restore data files if resuming
                    if 'data_files' not in results and loaded_checkpoint:
                        data_files_meta = loaded_checkpoint.metadata.get('data_files', {})
                        results['data_files'] = {k: Path(v) for k, v in data_files_meta.items()}
                    
                    domain_results = await self.run_domain_scanning(parameters, results['data_files'])
                    results['domain_results'] = domain_results
                    self.save_checkpoint(PipelineStage.ONTOLOGY_PROCESSING, parameters)
            
            # Stage 4: Ontology Processing
            if start_stage.value <= PipelineStage.ONTOLOGY_PROCESSING.value:
                async with self._progress_context(PipelineStage.ONTOLOGY_PROCESSING):
                    ontology_results = await self.run_ontology_processing(parameters, results['data_files'])
                    results['ontology_results'] = ontology_results
                    self.save_checkpoint(PipelineStage.STATISTICAL_INFERENCE, parameters)
            
            # Stage 5: Statistical Inference
            if start_stage.value <= PipelineStage.STATISTICAL_INFERENCE.value:
                async with self._progress_context(PipelineStage.STATISTICAL_INFERENCE):
                    inference_results = await self.run_statistical_inference(
                        parameters, 
                        results['domain_results'],
                        results['ontology_results']
                    )
                    results['inference_results'] = inference_results
                    self.save_checkpoint(PipelineStage.RESULTS_EXPORT, parameters)
            
            # Stage 6: Results Export
            if start_stage.value <= PipelineStage.RESULTS_EXPORT.value:
                async with self._progress_context(PipelineStage.RESULTS_EXPORT):
                    output_files = await self.export_results(
                        parameters,
                        results['inference_results'],
                        results['ontology_results']
                    )
                    results['output_files'] = output_files
                    self.save_checkpoint(PipelineStage.CLEANUP, parameters)
            
            # Stage 7: Cleanup
            async with self._progress_context(PipelineStage.CLEANUP):
                # Final progress update
                if self.progress:
                    self.progress.update_stage(PipelineStage.COMPLETED)
                
                # Clear checkpoint on successful completion
                self.clear_checkpoint()
                
                # Final statistics
                total_time = datetime.now() - start_time
                logger.info(f"Pipeline completed successfully in {total_time}")
                
                if results.get('inference_results'):
                    stats = results['inference_results']['statistics']
                    logger.info(f"Found {stats['significant_associations']} significant "
                               f"domain-GO associations (FDR < {stats['fdr_threshold']})")
                
                results['run_metadata'] = {
                    'start_time': start_time.isoformat(),
                    'end_time': datetime.now().isoformat(),
                    'total_duration': str(total_time),
                    'parameters': parameters,
                    'resumed_from': loaded_checkpoint.stage.name if loaded_checkpoint else None
                }
            
            return results
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            # Save checkpoint on failure for potential resume
            if hasattr(self, 'current_checkpoint') and self.current_checkpoint:
                logger.info("Checkpoint saved for potential resume")
            raise PipelineError(f"Pipeline execution failed: {e}") from e


# CLI Interface
@click.command()
@click.option('--resume', is_flag=True, help='Resume from checkpoint if available')
@click.option('--force-download', is_flag=True, 
              help='Force re-download of all datasets')
@click.option('--datasets', multiple=True, 
              help='Specific datasets to download (can specify multiple)')
@click.option('--num-cores', type=int, 
              help=f'Number of cores for parallel processing (default: {config.NUM_CORES})')
@click.option('--fdr-threshold', type=float,
              help=f'FDR threshold for significance (default: {config.FDR_THRESHOLD})')
@click.option('--export-all-results', is_flag=True,
              help='Export all results, not just significant associations')
@click.option('--max-workers', type=int, default=4,
              help='Maximum concurrent workers for pipeline operations')
@click.option('--disable-async', is_flag=True,
              help='Disable asynchronous operations')
@click.option('--checkpoint-file', type=click.Path(),
              help='Custom checkpoint file path')
@click.option('--log-level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Logging level')
@click.option('--dry-run', is_flag=True,
              help='Show what would be done without executing')
def main(
    resume: bool,
    force_download: bool,
    datasets: tuple,
    num_cores: Optional[int],
    fdr_threshold: Optional[float],
    export_all_results: bool,
    max_workers: int,
    disable_async: bool,
    checkpoint_file: Optional[str],
    log_level: str,
    dry_run: bool
) -> None:
    """
    dcGO Pipeline - Domain-centric Gene Ontology functional annotation.
    
    This pipeline performs comprehensive functional annotation of protein domains
    using Gene Ontology terms through statistical inference and ontology processing.
    
    Example usage:
        dcgo-pipeline --num-cores 8 --fdr-threshold 0.01
        dcgo-pipeline --resume --export-all-results
        dcgo-pipeline --datasets uniprot_sprot --datasets go_ontology
    """
    # Configure root logger level
    logger.remove()
    logger.add(sys.stderr, level=log_level)
    
    if dry_run:
        click.echo("DRY RUN MODE - No actual processing will be performed")
        click.echo(f"Parameters that would be used:")
        click.echo(f"  Resume: {resume}")
        click.echo(f"  Force download: {force_download}")
        click.echo(f"  Datasets: {list(datasets) if datasets else 'all'}")
        click.echo(f"  Cores: {num_cores or config.NUM_CORES}")
        click.echo(f"  FDR threshold: {fdr_threshold or config.FDR_THRESHOLD}")
        click.echo(f"  Export all results: {export_all_results}")
        click.echo(f"  Max workers: {max_workers}")
        click.echo(f"  Async enabled: {not disable_async}")
        return
    
    # Convert datasets tuple to list
    datasets_list = list(datasets) if datasets else None
    
    # Validate datasets
    if datasets_list:
        invalid_datasets = set(datasets_list) - set(config.DATASOURCES.keys())
        if invalid_datasets:
            click.echo(f"Error: Invalid datasets specified: {invalid_datasets}", err=True)
            click.echo(f"Available datasets: {list(config.DATASOURCES.keys())}")
            sys.exit(1)
    
    # Initialize and run pipeline
    try:
        pipeline = dcGOPipeline(
            max_workers=max_workers,
            checkpoint_file=Path(checkpoint_file) if checkpoint_file else None,
            enable_async=not disable_async
        )
        
        # Run the pipeline
        results = asyncio.run(pipeline.run_full_pipeline(
            resume=resume,
            force_download=force_download,
            datasets=datasets_list,
            num_cores=num_cores,
            fdr_threshold=fdr_threshold,
            export_all_results=export_all_results
        ))
        
        # Display results summary
        click.echo("\n" + "="*60)
        click.echo("PIPELINE COMPLETED SUCCESSFULLY")
        click.echo("="*60)
        
        if 'inference_results' in results:
            stats = results['inference_results']['statistics']
            click.echo(f"Significant associations found: {stats['significant_associations']}")
            click.echo(f"Total statistical tests: {stats['total_tests']}")
            click.echo(f"FDR threshold used: {stats['fdr_threshold']}")
        
        if 'output_files' in results:
            click.echo(f"\nOutput files generated:")
            for file_type, file_path in results['output_files'].items():
                click.echo(f"  {file_type}: {file_path}")
        
        metadata = results.get('run_metadata', {})
        if metadata:
            click.echo(f"\nRun duration: {metadata.get('total_duration', 'Unknown')}")
            if metadata.get('resumed_from'):
                click.echo(f"Resumed from: {metadata['resumed_from']}")
        
    except PipelineError as e:
        click.echo(f"Pipeline error: {e}", err=True)
        sys.exit(1)
    except KeyboardInterrupt:
        click.echo("\nPipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error occurred")
        click.echo(f"Unexpected error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()