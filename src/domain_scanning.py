"""
Domain Architecture Scanner for dcGO Pipeline

This module implements comprehensive domain architecture scanning using InterProScan,
including single domain detection, supra-domain generation, and batch processing
capabilities for large-scale proteomics datasets.

Author: dcGO Pipeline
"""

import subprocess
import pandas as pd
import psutil
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Union
from loguru import logger
from Bio import SeqIO
from dataclasses import dataclass
import tempfile
import shutil
import time
import json


@dataclass
class DomainHit:
    """Represents a single domain hit from InterProScan"""
    protein_id: str
    domain_id: str
    start: int
    stop: int
    score: float
    e_value: Optional[float] = None
    interpro_id: Optional[str] = None
    interpro_desc: Optional[str] = None
    go_annotations: Optional[str] = None


@dataclass
class DomainArchitecture:
    """Represents the complete domain architecture of a protein"""
    protein_id: str
    single_domains: List[str]
    supra_domains: List[str]
    domain_positions: Dict[str, Tuple[int, int]]
    total_length: int


class DomainArchitectureScanner:
    """
    Comprehensive domain architecture scanner using InterProScan.
    
    Features:
    - InterProScan integration with proper error handling
    - Single domain and supra-domain detection
    - Batch processing with progress tracking
    - Sequence preparation and validation
    - Results parsing and validation
    - Comprehensive logging and monitoring
    """
    
    def __init__(
        self,
        interproscan_path: Path,
        num_cores: Optional[int] = None,
        temp_dir: Optional[Path] = None,
        max_sequence_length: int = 50000,
        min_sequence_length: int = 20
    ):
        """
        Initialize the domain scanner.
        
        Args:
            interproscan_path: Path to InterProScan executable
            num_cores: Number of CPU cores to use (defaults to system cores)
            temp_dir: Directory for temporary files (defaults to system temp)
            max_sequence_length: Maximum allowed sequence length
            min_sequence_length: Minimum allowed sequence length
        """
        self.interproscan_path = Path(interproscan_path)
        self.cores = num_cores or psutil.cpu_count()
        self.temp_dir = Path(temp_dir) if temp_dir else Path(tempfile.gettempdir())
        self.max_sequence_length = max_sequence_length
        self.min_sequence_length = min_sequence_length
        
        # Validate InterProScan installation
        self._validate_interproscan()
        
        # Create working directory
        self.work_dir = self.temp_dir / f"dcgo_domain_scan_{int(time.time())}"
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("DomainArchitectureScanner initialized:")
        logger.info(f"  InterProScan: {self.interproscan_path}")
        logger.info(f"  CPU cores: {self.cores}")
        logger.info(f"  Working directory: {self.work_dir}")

    def _validate_interproscan(self) -> None:
        """Validate InterProScan installation and permissions"""
        if not self.interproscan_path.exists():
            raise FileNotFoundError(f"InterProScan not found at {self.interproscan_path}")
        
        if not self.interproscan_path.is_file():
            raise ValueError(f"InterProScan path is not a file: {self.interproscan_path}")
        
        if not self.interproscan_path.stat().st_mode & 0o111:
            raise PermissionError(f"InterProScan is not executable: {self.interproscan_path}")
        
        logger.info("InterProScan installation validated")

    def prepare_sequences_for_scanning(
        self,
        sequence_files: List[Path],
        output_file: Path,
        max_sequences: Optional[int] = None,
        validate_sequences: bool = True
    ) -> Path:
        """
        Combine and prepare sequence files for InterProScan.
        
        Args:
            sequence_files: List of FASTA files to combine
            output_file: Output path for combined sequences
            max_sequences: Maximum number of sequences to process (for testing)
            validate_sequences: Whether to validate sequence format and content
            
        Returns:
            Path to the prepared sequence file
        """
        if output_file.exists() and output_file.stat().st_size > 0:
            logger.info(f"Combined sequence file already exists at {output_file}")
            return output_file
        
        logger.info("Combining and preparing sequence files for InterProScan")
        
        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        sequence_count = 0
        invalid_sequences = 0
        
        with open(output_file, 'w') as outfile:
            for seq_file in sequence_files:
                if not seq_file.exists():
                    logger.warning(f"Sequence file not found: {seq_file}")
                    continue
                    
                logger.info(f"Processing {seq_file}")
                
                # Determine file opener based on extension
                if seq_file.suffix.lower() == '.gz':
                    open_func = lambda f: gzip.open(f, 'rt')
                else:
                    open_func = open
                
                try:
                    with open_func(seq_file) as infile:
                        for record in SeqIO.parse(infile, 'fasta'):
                            # Clean up sequence ID for InterProScan compatibility
                            clean_id = self._clean_sequence_id(record.id)
                            
                            # Validate sequence if requested
                            if validate_sequences:
                                if not self._validate_sequence(record, clean_id):
                                    invalid_sequences += 1
                                    continue
                            
                            # Create clean record
                            record.id = clean_id
                            record.description = ""  # Remove description to avoid parsing issues
                            
                            # Write to output file
                            SeqIO.write(record, outfile, 'fasta')
                            sequence_count += 1
                            
                            # Check sequence limit
                            if max_sequences and sequence_count >= max_sequences:
                                logger.info(f"Reached sequence limit of {max_sequences}")
                                break
                    
                    if max_sequences and sequence_count >= max_sequences:
                        break
                        
                except Exception as e:
                    logger.error(f"Error processing {seq_file}: {e}")
                    continue
        
        logger.info("Sequence preparation complete:")
        logger.info(f"  Total sequences: {sequence_count}")
        logger.info(f"  Invalid sequences skipped: {invalid_sequences}")
        logger.info(f"  Output file: {output_file}")
        
        if sequence_count == 0:
            raise ValueError("No valid sequences found in input files")
        
        return output_file

    def _clean_sequence_id(self, seq_id: str) -> str:
        """Clean sequence ID for InterProScan compatibility"""
        # Remove pipe-separated UniProt format parts if present
        if '|' in seq_id:
            parts = seq_id.split('|')
            # For UniProt format: sp|P12345|PROTEIN_HUMAN, use the accession (P12345)
            if len(parts) >= 2 and parts[0] in ['sp', 'tr']:
                clean_id = parts[1]
            else:
                clean_id = parts[0]
        else:
            clean_id = seq_id
        
        # Remove any remaining problematic characters
        clean_id = clean_id.replace(' ', '_').replace('/', '_').replace('\\', '_')
        
        # Ensure ID starts with alphanumeric character
        if not clean_id[0].isalnum():
            clean_id = 'seq_' + clean_id
        
        return clean_id

    def _validate_sequence(self, record, clean_id: str) -> bool:
        """Validate sequence for InterProScan processing"""
        seq_str = str(record.seq)
        
        # Check sequence length
        if len(seq_str) < self.min_sequence_length:
            logger.debug(f"Sequence {clean_id} too short ({len(seq_str)} aa)")
            return False
        
        if len(seq_str) > self.max_sequence_length:
            logger.debug(f"Sequence {clean_id} too long ({len(seq_str)} aa)")
            return False
        
        # Check for valid amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWYXBZJ*-')
        invalid_chars = set(seq_str.upper()) - valid_aa
        if invalid_chars:
            logger.debug(f"Sequence {clean_id} contains invalid characters: {invalid_chars}")
            return False
        
        # Check for excessive ambiguous residues
        ambiguous_count = sum(1 for c in seq_str.upper() if c in 'XBZ')
        if ambiguous_count / len(seq_str) > 0.5:
            logger.debug(f"Sequence {clean_id} has too many ambiguous residues ({ambiguous_count}/{len(seq_str)})")
            return False
        
        return True

    def run_interproscan_batch(
        self,
        input_fasta: Path,
        output_dir: Path,
        applications: Optional[List[str]] = None,
        timeout_hours: int = 24,
        retry_attempts: int = 3
    ) -> Path:
        """
        Run InterProScan on protein sequences with robust error handling.
        
        Args:
            input_fasta: Path to input FASTA file
            output_dir: Directory for output files
            applications: List of InterProScan applications to run (defaults to Pfam)
            timeout_hours: Timeout in hours for the InterProScan run
            retry_attempts: Number of retry attempts on failure
            
        Returns:
            Path to the TSV output file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_dir / "interpro_results.tsv"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            logger.info(f"InterProScan results already exist at {output_file}")
            return output_file
        
        # Default to Pfam for domain architecture analysis
        if applications is None:
            applications = ['Pfam']
        
        # Validate input file
        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
        
        # Count sequences for progress estimation
        seq_count = self._count_sequences(input_fasta)
        logger.info(f"Running InterProScan on {seq_count} sequences")
        
        for attempt in range(retry_attempts):
            try:
                success = self._run_interproscan_attempt(
                    input_fasta, output_file, applications, timeout_hours, attempt + 1
                )
                
                if success:
                    logger.info("InterProScan completed successfully")
                    return output_file
                    
            except Exception as e:
                logger.error(f"InterProScan attempt {attempt + 1} failed: {e}")
                if attempt < retry_attempts - 1:
                    logger.info("Retrying in 30 seconds...")
                    time.sleep(30)
                else:
                    raise
        
        raise RuntimeError(f"InterProScan failed after {retry_attempts} attempts")

    def _count_sequences(self, fasta_file: Path) -> int:
        """Count sequences in FASTA file"""
        count = 0
        try:
            if fasta_file.suffix.lower() == '.gz':
                with gzip.open(fasta_file, 'rt') as f:
                    for line in f:
                        if line.startswith('>'):
                            count += 1
            else:
                with open(fasta_file, 'r') as f:
                    for line in f:
                        if line.startswith('>'):
                            count += 1
        except Exception as e:
            logger.warning(f"Could not count sequences in {fasta_file}: {e}")
            count = 0
        
        return count

    def _run_interproscan_attempt(
        self,
        input_fasta: Path,
        output_file: Path,
        applications: List[str],
        timeout_hours: int,
        attempt: int
    ) -> bool:
        """Single attempt to run InterProScan"""
        cmd = [
            str(self.interproscan_path),
            '-i', str(input_fasta),
            '-f', 'TSV',
            '-o', str(output_file),
            '--cpu', str(self.cores),
            '--disable-precalc',
            '--goterms',
            '--pathways',
            '--verbose'
        ]
        
        # Add specific applications
        if applications:
            cmd.extend(['-appl', ','.join(applications)])
        
        logger.info(f"InterProScan attempt {attempt}: {' '.join(cmd)}")
        
        try:
            # Create log file for this attempt
            log_file = output_file.parent / f"interproscan_attempt_{attempt}.log"
            
            with open(log_file, 'w') as log_f:
                result = subprocess.run(
                    cmd,
                    check=True,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    timeout=timeout_hours * 3600
                )
            
            # Verify output file was created and has content
            if not output_file.exists():
                logger.error("InterProScan completed but no output file was created")
                return False
            
            if output_file.stat().st_size == 0:
                logger.error("InterProScan output file is empty")
                return False
            
            # Log basic statistics
            line_count = sum(1 for _ in open(output_file))
            logger.info(f"InterProScan produced {line_count} result lines")
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"InterProScan failed with exit code {e.returncode}")
            self._log_interproscan_error(log_file if 'log_file' in locals() else None)
            return False
            
        except subprocess.TimeoutExpired:
            logger.error(f"InterProScan timed out after {timeout_hours} hours")
            return False
        
        except Exception as e:
            logger.error(f"Unexpected error running InterProScan: {e}")
            return False

    def _log_interproscan_error(self, log_file: Optional[Path]) -> None:
        """Log InterProScan error details"""
        if log_file and log_file.exists():
            try:
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    if log_content:
                        logger.error("InterProScan log output (last 20 lines):")
                        lines = log_content.strip().split('\n')
                        for line in lines[-20:]:
                            logger.error(f"  {line}")
            except Exception as e:
                logger.error(f"Could not read InterProScan log file: {e}")

    def parse_interpro_output(
        self,
        tsv_file: Path,
        applications_filter: Optional[Set[str]] = None
    ) -> Dict[str, List[str]]:
        """
        Parse InterProScan TSV output into protein -> domain mapping.
        
        Args:
            tsv_file: Path to InterProScan TSV output file
            applications_filter: Set of applications to include (defaults to Pfam)
            
        Returns:
            Dictionary mapping protein IDs to lists of domain IDs (including supra-domains)
        """
        logger.info(f"Parsing InterProScan output from {tsv_file}")
        
        if not tsv_file.exists():
            raise FileNotFoundError(f"InterProScan output file not found: {tsv_file}")
        
        # Default to Pfam only for domain architecture
        if applications_filter is None:
            applications_filter = {'Pfam'}
        
        # Define InterProScan TSV column names
        columns = [
            'protein_id', 'seq_md5', 'seq_length', 'analysis', 'signature_id',
            'signature_desc', 'start', 'stop', 'score', 'status', 'date',
            'interpro_id', 'interpro_desc', 'go_annotations', 'pathways'
        ]
        
        try:
            # Read TSV file
            df = pd.read_csv(
                tsv_file,
                sep='\t',
                names=columns,
                comment='#',
                na_values=['-', ''],
                dtype={
                    'protein_id': str,
                    'seq_length': 'Int64',
                    'analysis': str,
                    'signature_id': str,
                    'start': 'Int64',
                    'stop': 'Int64',
                    'score': float
                },
                low_memory=False
            )
            
            logger.info(f"Loaded {len(df)} domain hits from InterProScan output")
            
        except Exception as e:
            logger.error(f"Error reading InterProScan output: {e}")
            raise
        
        if df.empty:
            logger.warning("InterProScan output is empty")
            return {}
        
        # Filter for specified applications
        filtered_hits = df[df['analysis'].isin(applications_filter)].copy()
        
        if filtered_hits.empty:
            logger.warning(f"No hits found for applications: {applications_filter}")
            return {}
        
        logger.info(f"Found {len(filtered_hits)} hits from {applications_filter} across {filtered_hits['protein_id'].nunique()} proteins")
        
        # Parse domain hits
        domain_hits = self._parse_domain_hits(filtered_hits)
        
        # Create domain architecture map
        return self._create_domain_architecture_map(domain_hits)

    def _parse_domain_hits(self, df: pd.DataFrame) -> Dict[str, List[DomainHit]]:
        """Parse DataFrame into DomainHit objects grouped by protein"""
        domain_hits = {}
        
        for _, row in df.iterrows():
            protein_id = row['protein_id']
            
            hit = DomainHit(
                protein_id=protein_id,
                domain_id=row['signature_id'],
                start=int(row['start']) if pd.notna(row['start']) else 0,
                stop=int(row['stop']) if pd.notna(row['stop']) else 0,
                score=float(row['score']) if pd.notna(row['score']) else 0.0,
                e_value=None,  # Not directly available in standard TSV format
                interpro_id=row['interpro_id'] if pd.notna(row['interpro_id']) else None,
                interpro_desc=row['interpro_desc'] if pd.notna(row['interpro_desc']) else None,
                go_annotations=row['go_annotations'] if pd.notna(row['go_annotations']) else None
            )
            
            if protein_id not in domain_hits:
                domain_hits[protein_id] = []
            domain_hits[protein_id].append(hit)
        
        return domain_hits

    def _create_domain_architecture_map(
        self,
        domain_hits: Dict[str, List[DomainHit]]
    ) -> Dict[str, List[str]]:
        """
        Create protein -> domain architecture mapping with supra-domains.
        
        Args:
            domain_hits: Dictionary of protein ID -> list of DomainHit objects
            
        Returns:
            Dictionary mapping protein IDs to lists of domain IDs (single + supra-domains)
        """
        protein_domains = {}
        
        for protein_id, hits in domain_hits.items():
            # Sort domains by start position to maintain order
            hits_sorted = sorted(hits, key=lambda x: x.start)
            
            # Extract single domain IDs
            single_domains = [hit.domain_id for hit in hits_sorted]
            
            # Generate supra-domains (contiguous combinations)
            supra_domains = self._generate_supra_domains(single_domains)
            
            # Combine single and supra-domains
            all_domains = single_domains + supra_domains
            
            # Store domain architecture
            protein_domains[protein_id] = all_domains
        
        logger.info(f"Created domain architectures for {len(protein_domains)} proteins:")
        
        # Log statistics
        total_single = sum(len([d for d in domains if ',' not in d]) 
                          for domains in protein_domains.values())
        total_supra = sum(len([d for d in domains if ',' in d]) 
                         for domains in protein_domains.values())
        
        logger.info(f"  Single domains: {total_single}")
        logger.info(f"  Supra-domains: {total_supra}")
        logger.info(f"  Total domain features: {total_single + total_supra}")
        
        return protein_domains

    def _generate_supra_domains(
        self,
        domain_list: List[str],
        max_length: int = 3
    ) -> List[str]:
        """
        Generate supra-domain identifiers for contiguous domain combinations.
        
        Args:
            domain_list: List of domain IDs in positional order
            max_length: Maximum length of supra-domain combinations
            
        Returns:
            List of supra-domain identifiers (comma-separated domain combinations)
        """
        supra_domains = []
        
        if len(domain_list) < 2:
            return supra_domains
        
        # Generate all contiguous combinations up to max_length
        for length in range(2, min(len(domain_list) + 1, max_length + 1)):
            for i in range(len(domain_list) - length + 1):
                combo = ','.join(domain_list[i:i+length])
                supra_domains.append(combo)
        
        return supra_domains

    def create_domain_architecture_objects(
        self,
        domain_hits: Dict[str, List[DomainHit]]
    ) -> Dict[str, DomainArchitecture]:
        """
        Create detailed DomainArchitecture objects for each protein.
        
        Args:
            domain_hits: Dictionary of protein ID -> list of DomainHit objects
            
        Returns:
            Dictionary mapping protein IDs to DomainArchitecture objects
        """
        architectures = {}
        
        for protein_id, hits in domain_hits.items():
            # Sort by position
            hits_sorted = sorted(hits, key=lambda x: x.start)
            
            # Extract components
            single_domains = [hit.domain_id for hit in hits_sorted]
            supra_domains = self._generate_supra_domains(single_domains)
            
            # Create position mapping
            domain_positions = {
                hit.domain_id: (hit.start, hit.stop) for hit in hits_sorted
            }
            
            # Calculate total protein length (approximate from domain spans)
            if hits_sorted:
                total_length = max(hit.stop for hit in hits_sorted)
            else:
                total_length = 0
            
            architecture = DomainArchitecture(
                protein_id=protein_id,
                single_domains=single_domains,
                supra_domains=supra_domains,
                domain_positions=domain_positions,
                total_length=total_length
            )
            
            architectures[protein_id] = architecture
        
        return architectures

    def export_domain_architectures(
        self,
        architectures: Dict[str, DomainArchitecture],
        output_file: Path
    ) -> None:
        """Export domain architectures to JSON format"""
        logger.info(f"Exporting domain architectures to {output_file}")
        
        # Convert to JSON-serializable format
        export_data = {}
        for protein_id, arch in architectures.items():
            export_data[protein_id] = {
                'single_domains': arch.single_domains,
                'supra_domains': arch.supra_domains,
                'domain_positions': {k: list(v) for k, v in arch.domain_positions.items()},
                'total_length': arch.total_length
            }
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2, sort_keys=True)
        
        logger.info(f"Exported {len(architectures)} domain architectures")

    def get_statistics(self) -> Dict[str, Union[int, float, str]]:
        """Get scanner statistics and configuration"""
        return {
            'interproscan_path': str(self.interproscan_path),
            'cpu_cores': self.cores,
            'work_directory': str(self.work_dir),
            'max_sequence_length': self.max_sequence_length,
            'min_sequence_length': self.min_sequence_length,
            'temp_directory': str(self.temp_dir)
        }

    def cleanup(self) -> None:
        """Clean up temporary files and directories"""
        try:
            if self.work_dir.exists():
                shutil.rmtree(self.work_dir)
                logger.info(f"Cleaned up working directory: {self.work_dir}")
        except Exception as e:
            logger.warning(f"Could not clean up working directory: {e}")

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup"""
        self.cleanup()


def main():
    """Example usage of DomainArchitectureScanner"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Run domain architecture scanning")
    parser.add_argument("--interproscan", required=True, help="Path to InterProScan executable")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--cores", type=int, help="Number of CPU cores")
    
    args = parser.parse_args()
    
    # Setup logging
    logger.add("domain_scanning.log", level="INFO")
    
    with DomainArchitectureScanner(
        interproscan_path=Path(args.interproscan),
        num_cores=args.cores
    ) as scanner:
        
        # Run InterProScan
        results = scanner.run_interproscan_batch(
            input_fasta=Path(args.input),
            output_dir=Path(args.output_dir)
        )
        
        # Parse results
        domain_map = scanner.parse_interpro_output(results)
        
        # Export summary
        summary_file = Path(args.output_dir) / "domain_summary.json"
        with open(summary_file, 'w') as f:
            json.dump({
                'total_proteins': len(domain_map),
                'proteins_with_domains': len([p for p, d in domain_map.items() if d]),
                'total_domain_features': sum(len(d) for d in domain_map.values()),
                'scanner_stats': scanner.get_statistics()
            }, f, indent=2)
        
        logger.info(f"Domain scanning complete. Results in {args.output_dir}")


if __name__ == "__main__":
    main()