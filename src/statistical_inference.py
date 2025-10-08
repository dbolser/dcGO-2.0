# src/statistical_inference.py
"""
Statistical Inference Engine for dcGO Pipeline

This module implements the complete statistical inference system for determining 
significant domain-GO term associations using Fisher's exact test, hypergeometric
scoring, and FDR correction according to the dcGO methodology.

Author: dcGO Pipeline
Version: 1.0.0
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import fdrcorrection
from typing import Dict, List, Set, Tuple, Optional, Union
from loguru import logger
from dataclasses import dataclass, field
from tqdm import tqdm
import warnings
from pathlib import Path


@dataclass
class AssociationResult:
    """
    Data class for storing statistical association results between domains and GO terms.
    
    This class encapsulates all the statistical metrics computed for each domain-GO term
    pair, including contingency table values, Fisher's exact test results, hypergeometric
    scores, and multiple testing correction values.
    
    Attributes:
        domain (str): Domain identifier (e.g., 'PF00001')
        go_term (str): GO term identifier (e.g., 'GO:0003674')
        a (int): Count of proteins with both domain and GO term
        b (int): Count of proteins with GO term but not domain  
        c (int): Count of proteins with domain but not GO term
        d (int): Count of proteins with neither domain nor GO term
        p_value (float): Fisher's exact test p-value (one-tailed, greater)
        odds_ratio (float): Odds ratio from Fisher's exact test
        hyper_score (float): Hypergeometric-based association score (1-100 scale)
        q_value (Optional[float]): FDR-corrected q-value (Benjamini-Hochberg)
        n_total (int): Total number of proteins in the analysis
        expected_overlap (float): Expected overlap under null hypothesis
        enrichment_fold (float): Fold enrichment over expected
    """
    domain: str
    go_term: str
    a: int  # both domain and GO term
    b: int  # GO term only  
    c: int  # domain only
    d: int  # neither
    p_value: float
    odds_ratio: float
    hyper_score: float
    q_value: Optional[float] = None
    n_total: int = field(default=0)
    expected_overlap: float = field(default=0.0)
    enrichment_fold: float = field(default=0.0)
    
    def __post_init__(self):
        """Calculate derived metrics after initialization."""
        self.n_total = self.a + self.b + self.c + self.d
        
        # Calculate expected overlap under null hypothesis
        domain_freq = (self.a + self.c) / self.n_total if self.n_total > 0 else 0
        go_freq = (self.a + self.b) / self.n_total if self.n_total > 0 else 0
        self.expected_overlap = domain_freq * go_freq * self.n_total
        
        # Calculate fold enrichment
        if self.expected_overlap > 0:
            self.enrichment_fold = self.a / self.expected_overlap
        else:
            self.enrichment_fold = float('inf') if self.a > 0 else 0.0
    
    @property
    def is_significant(self, alpha: float = 0.01) -> bool:
        """Check if association is statistically significant at given alpha level."""
        return self.q_value is not None and self.q_value < alpha
    
    def to_dict(self) -> Dict[str, Union[str, int, float]]:
        """Convert to dictionary for easy serialization."""
        return {
            'domain': self.domain,
            'go_term': self.go_term,
            'proteins_both': self.a,
            'proteins_go_only': self.b, 
            'proteins_domain_only': self.c,
            'proteins_neither': self.d,
            'total_proteins': self.n_total,
            'p_value': self.p_value,
            'q_value': self.q_value,
            'odds_ratio': self.odds_ratio,
            'hyper_score': self.hyper_score,
            'expected_overlap': self.expected_overlap,
            'enrichment_fold': self.enrichment_fold
        }


class StatisticalInferenceEngine:
    """
    Complete statistical inference engine for dcGO pipeline.
    
    This class implements the core statistical methods for identifying significant
    associations between protein domains and Gene Ontology terms. It uses Fisher's
    exact test for statistical significance, hypergeometric distributions for 
    association scoring, and Benjamini-Hochberg FDR correction for multiple testing.
    
    The engine follows the dcGO methodology:
    1. Build correspondence matrix of domain-GO co-occurrences
    2. Calculate 2x2 contingency tables for each domain-GO pair
    3. Apply Fisher's exact test (one-tailed, greater alternative)
    4. Compute hypergeometric association scores
    5. Apply FDR correction with Benjamini-Hochberg method
    6. Filter results based on significance thresholds
    
    Attributes:
        protein_domain_map (Dict[str, List[str]]): Protein ID -> list of domains
        protein_go_map (Dict[str, Set[str]]): Protein ID -> set of GO terms
        total_proteins (int): Total number of proteins in analysis
        domain_counts (Dict[str, int]): Domain frequency counts
        go_counts (Dict[str, int]): GO term frequency counts
        correspondence_matrix (Optional[pd.DataFrame]): Domain-GO co-occurrence matrix
    """
    
    def __init__(self, 
                 protein_domain_map: Dict[str, List[str]], 
                 protein_go_map: Dict[str, Set[str]],
                 validate_input: bool = True):
        """
        Initialize the statistical inference engine.
        
        Args:
            protein_domain_map: Dictionary mapping protein IDs to lists of domain IDs
            protein_go_map: Dictionary mapping protein IDs to sets of GO term IDs
            validate_input: Whether to validate input data consistency
            
        Raises:
            ValueError: If input validation fails
            TypeError: If inputs are not of expected types
        """
        if validate_input:
            self._validate_inputs(protein_domain_map, protein_go_map)
            
        self.protein_domain_map = protein_domain_map
        self.protein_go_map = protein_go_map
        self.total_proteins = len(protein_domain_map)
        
        # Initialize derived attributes
        self.domain_counts: Dict[str, int] = {}
        self.go_counts: Dict[str, int] = {}
        self.correspondence_matrix: Optional[pd.DataFrame] = None
        
        # Pre-compute marginal frequencies for efficiency
        self._compute_marginals()
        
        logger.info(f"Initialized StatisticalInferenceEngine with {self.total_proteins} proteins")
        logger.info(f"Found {len(self.domain_counts)} unique domains and {len(self.go_counts)} unique GO terms")
    
    def _validate_inputs(self, 
                        protein_domain_map: Dict[str, List[str]], 
                        protein_go_map: Dict[str, Set[str]]) -> None:
        """
        Validate input data for consistency and completeness.
        
        Args:
            protein_domain_map: Protein-domain mapping to validate
            protein_go_map: Protein-GO mapping to validate
            
        Raises:
            TypeError: If inputs are not of expected types
            ValueError: If data is inconsistent or invalid
        """
        # Type checking
        if not isinstance(protein_domain_map, dict):
            raise TypeError("protein_domain_map must be a dictionary")
        if not isinstance(protein_go_map, dict):
            raise TypeError("protein_go_map must be a dictionary")
            
        # Content validation
        if not protein_domain_map:
            raise ValueError("protein_domain_map cannot be empty")
        if not protein_go_map:
            raise ValueError("protein_go_map cannot be empty")
            
        # Check for common proteins
        domain_proteins = set(protein_domain_map.keys())
        go_proteins = set(protein_go_map.keys())
        common_proteins = domain_proteins & go_proteins
        
        if not common_proteins:
            raise ValueError("No proteins found with both domain and GO annotations")
        
        if len(common_proteins) < 10:
            warnings.warn(f"Only {len(common_proteins)} proteins with both annotations. "
                         "Results may be unreliable with small sample sizes.")
        
        # Validate data structure integrity
        for protein_id, domains in protein_domain_map.items():
            if not isinstance(domains, list):
                raise TypeError(f"Domain list for protein {protein_id} must be a list")
            if not all(isinstance(d, str) for d in domains):
                raise TypeError(f"All domains for protein {protein_id} must be strings")
                
        for protein_id, go_terms in protein_go_map.items():
            if not isinstance(go_terms, (set, list)):
                raise TypeError(f"GO terms for protein {protein_id} must be a set or list")
            if not all(isinstance(g, str) for g in go_terms):
                raise TypeError(f"All GO terms for protein {protein_id} must be strings")
    
    def _compute_marginals(self) -> None:
        """
        Pre-compute domain and GO term frequency counts for efficient lookup.
        
        This method calculates how many proteins contain each domain and each GO term.
        These marginal frequencies are used repeatedly in statistical calculations
        and are cached for performance optimization.
        """
        logger.info("Computing marginal frequencies for domains and GO terms")
        
        # Initialize counters
        self.domain_counts = {}
        self.go_counts = {}
        
        # Count domain occurrences
        for protein_id, domains in self.protein_domain_map.items():
            for domain in domains:
                self.domain_counts[domain] = self.domain_counts.get(domain, 0) + 1
        
        # Count GO term occurrences  
        for protein_id, go_terms in self.protein_go_map.items():
            for go_term in go_terms:
                self.go_counts[go_term] = self.go_counts.get(go_term, 0) + 1
        
        # Log statistics
        total_domain_annotations = sum(self.domain_counts.values())
        total_go_annotations = sum(self.go_counts.values())
        
        logger.info(f"Domain statistics: {len(self.domain_counts)} unique domains, "
                   f"{total_domain_annotations} total annotations")
        logger.info(f"GO statistics: {len(self.go_counts)} unique GO terms, "
                   f"{total_go_annotations} total annotations")
        
        # Log frequency distributions
        if self.domain_counts:
            domain_freq_stats = pd.Series(list(self.domain_counts.values()))
            logger.info(f"Domain frequency stats - Mean: {domain_freq_stats.mean():.1f}, "
                       f"Median: {domain_freq_stats.median():.1f}, "
                       f"Max: {domain_freq_stats.max()}")
        
        if self.go_counts:
            go_freq_stats = pd.Series(list(self.go_counts.values()))
            logger.info(f"GO term frequency stats - Mean: {go_freq_stats.mean():.1f}, "
                       f"Median: {go_freq_stats.median():.1f}, "
                       f"Max: {go_freq_stats.max()}")
    
    def build_correspondence_matrix(self, cache_result: bool = True) -> pd.DataFrame:
        """
        Build the correspondence matrix showing domain-GO term co-occurrences.
        
        This method creates a matrix where rows are GO terms, columns are domains,
        and values are the number of proteins that have both the GO term and domain.
        This matrix is fundamental for calculating statistical associations.
        
        Args:
            cache_result: Whether to cache the matrix for reuse
            
        Returns:
            pd.DataFrame: Correspondence matrix with GO terms as rows and domains as columns
            
        Raises:
            RuntimeError: If matrix construction fails
        """
        if self.correspondence_matrix is not None and cache_result:
            logger.info("Using cached correspondence matrix")
            return self.correspondence_matrix
        
        logger.info("Building correspondence matrix for domain-GO term co-occurrences")
        
        try:
            # Create protein-domain pairs
            domain_pairs = []
            for protein_id, domains in self.protein_domain_map.items():
                for domain in domains:
                    domain_pairs.append({'protein': protein_id, 'domain': domain})
            
            # Create protein-GO pairs
            go_pairs = []
            for protein_id, go_terms in self.protein_go_map.items():
                for go_term in go_terms:
                    go_pairs.append({'protein': protein_id, 'go_term': go_term})
            
            # Convert to DataFrames
            domain_df = pd.DataFrame(domain_pairs)
            go_df = pd.DataFrame(go_pairs)
            
            if domain_df.empty or go_df.empty:
                raise RuntimeError("No valid domain-GO pairs found for matrix construction")
            
            # Merge on protein to get co-occurrences
            logger.info(f"Merging {len(domain_df)} domain pairs with {len(go_df)} GO pairs")
            merged = domain_df.merge(go_df, on='protein', how='inner')
            
            if merged.empty:
                raise RuntimeError("No co-occurrences found between domains and GO terms")
            
            # Create correspondence matrix using crosstab
            correspondence_matrix = pd.crosstab(
                merged['go_term'], 
                merged['domain'], 
                margins=False
            ).fillna(0).astype(int)
            
            logger.info(f"Built correspondence matrix: {correspondence_matrix.shape[0]} GO terms × "
                       f"{correspondence_matrix.shape[1]} domains")
            logger.info(f"Matrix density: {(correspondence_matrix > 0).sum().sum()} / "
                       f"{correspondence_matrix.size} = "
                       f"{100 * (correspondence_matrix > 0).sum().sum() / correspondence_matrix.size:.2f}%")
            
            # Cache result if requested
            if cache_result:
                self.correspondence_matrix = correspondence_matrix
                
            return correspondence_matrix
            
        except Exception as e:
            logger.error(f"Failed to build correspondence matrix: {str(e)}")
            raise RuntimeError(f"Correspondence matrix construction failed: {str(e)}") from e
    
    def calculate_contingency_values(self, 
                                   domain: str, 
                                   go_term: str, 
                                   correspondence_matrix: Optional[pd.DataFrame] = None) -> Tuple[int, int, int, int]:
        """
        Calculate 2x2 contingency table values for a specific domain-GO term pair.
        
        The contingency table has the structure:
                    | GO+  | GO-  |
            --------+------+------+
            Domain+ |  a   |  c   |
            Domain- |  b   |  d   |
            
        Where:
        - a: proteins with both domain and GO term
        - b: proteins with GO term but not domain
        - c: proteins with domain but not GO term  
        - d: proteins with neither domain nor GO term
        
        Args:
            domain: Domain identifier to analyze
            go_term: GO term identifier to analyze
            correspondence_matrix: Pre-computed correspondence matrix (optional)
            
        Returns:
            Tuple[int, int, int, int]: Values (a, b, c, d) for contingency table
            
        Raises:
            ValueError: If domain or GO term is invalid
        """
        if correspondence_matrix is None:
            correspondence_matrix = self.correspondence_matrix or self.build_correspondence_matrix()
        
        # Get co-occurrence count (a)
        if (go_term in correspondence_matrix.index and 
            domain in correspondence_matrix.columns):
            a = correspondence_matrix.loc[go_term, domain]
        else:
            a = 0
        
        # Get marginal totals from pre-computed counts
        domain_total = self.domain_counts.get(domain, 0)
        go_term_total = self.go_counts.get(go_term, 0)
        
        # Calculate other contingency table cells
        b = go_term_total - a      # GO term but not domain
        c = domain_total - a       # domain but not GO term
        d = self.total_proteins - (a + b + c)  # neither
        
        # Validate contingency table
        if any(x < 0 for x in [a, b, c, d]):
            raise ValueError(f"Invalid contingency table values for {domain}-{go_term}: "
                           f"a={a}, b={b}, c={c}, d={d}")
        
        return a, b, c, d
    
    def run_statistical_tests(self, 
                             min_cooccurrence: int = 3,
                             fdr_alpha: float = 0.01,
                             enable_progress_bar: bool = True) -> List[AssociationResult]:
        """
        Run complete statistical testing pipeline for all domain-GO term pairs.
        
        This method performs the following steps:
        1. Build correspondence matrix of co-occurrences
        2. Filter pairs with minimum co-occurrence threshold
        3. Calculate Fisher's exact test for each pair
        4. Compute hypergeometric association scores
        5. Apply Benjamini-Hochberg FDR correction
        6. Filter to significant associations
        
        Args:
            min_cooccurrence: Minimum number of proteins with both domain and GO term
            fdr_alpha: False discovery rate threshold for significance
            enable_progress_bar: Whether to show progress bar during computation
            
        Returns:
            List[AssociationResult]: Significant associations after FDR correction
            
        Raises:
            RuntimeError: If statistical testing fails
        """
        logger.info(f"Starting statistical testing with minimum co-occurrence: {min_cooccurrence}")
        logger.info(f"FDR threshold: {fdr_alpha}")
        
        try:
            # Build correspondence matrix
            correspondence_matrix = self.build_correspondence_matrix()
            
            # Identify testable pairs (meeting minimum co-occurrence)
            testable_pairs = []
            for go_term in correspondence_matrix.index:
                for domain in correspondence_matrix.columns:
                    if correspondence_matrix.loc[go_term, domain] >= min_cooccurrence:
                        testable_pairs.append((go_term, domain))
            
            logger.info(f"Found {len(testable_pairs)} testable domain-GO pairs "
                       f"(≥{min_cooccurrence} co-occurrences)")
            
            if not testable_pairs:
                logger.warning("No testable pairs found - consider lowering min_cooccurrence threshold")
                return []
            
            # Run statistical tests
            results = []
            
            # Setup progress tracking
            iterator = tqdm(testable_pairs, desc="Statistical testing") if enable_progress_bar else testable_pairs
            
            for go_term, domain in iterator:
                try:
                    result = self._test_single_association(
                        domain, go_term, correspondence_matrix
                    )
                    if result is not None:
                        results.append(result)
                        
                except Exception as e:
                    logger.warning(f"Failed to test {domain}-{go_term}: {str(e)}")
                    continue
            
            logger.info(f"Completed statistical testing for {len(results)} associations")
            
            # Apply FDR correction
            if results:
                significant_results = self._apply_fdr_correction(results, fdr_alpha)
                logger.info(f"Found {len(significant_results)} significant associations "
                           f"after FDR correction (α={fdr_alpha})")
                return significant_results
            else:
                logger.warning("No valid statistical test results obtained")
                return []
                
        except Exception as e:
            logger.error(f"Statistical testing failed: {str(e)}")
            raise RuntimeError(f"Statistical testing pipeline failed: {str(e)}") from e
    
    def _test_single_association(self, 
                               domain: str, 
                               go_term: str, 
                               correspondence_matrix: pd.DataFrame) -> Optional[AssociationResult]:
        """
        Perform statistical testing for a single domain-GO term association.
        
        Args:
            domain: Domain identifier
            go_term: GO term identifier  
            correspondence_matrix: Pre-computed correspondence matrix
            
        Returns:
            AssociationResult or None if testing fails
        """
        try:
            # Calculate contingency table
            a, b, c, d = self.calculate_contingency_values(domain, go_term, correspondence_matrix)
            
            # Skip if insufficient data
            if a + b == 0 or a + c == 0:
                return None
            
            # Fisher's exact test (one-tailed for enrichment)
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
            
            # Handle edge cases
            if np.isnan(p_value) or np.isnan(odds_ratio):
                return None
            
            # Calculate hypergeometric score
            hyper_score = self._calculate_hypergeometric_score(a, b, c, d)
            
            # Create result object
            result = AssociationResult(
                domain=domain,
                go_term=go_term,
                a=a, b=b, c=c, d=d,
                p_value=p_value,
                odds_ratio=odds_ratio,
                hyper_score=hyper_score
            )
            
            return result
            
        except Exception as e:
            logger.debug(f"Failed to test {domain}-{go_term}: {str(e)}")
            return None
    
    def _calculate_hypergeometric_score(self, a: int, b: int, c: int, d: int) -> float:
        """
        Calculate hypergeometric-based association score on 1-100 scale.
        
        This score represents the strength of association between a domain and GO term,
        scaled to an intuitive 1-100 range where higher values indicate stronger associations.
        
        Args:
            a: Proteins with both domain and GO term
            b: Proteins with GO term only
            c: Proteins with domain only  
            d: Proteins with neither
            
        Returns:
            float: Association score between 1.0 and 100.0
        """
        n = a + b + c + d  # total proteins
        k = a + c          # proteins with domain
        m = a + b          # proteins with GO term
        x = a              # proteins with both
        
        if k == 0 or m == 0 or x == 0:
            return 0.0
        
        try:
            # Calculate hypergeometric survival function (1 - CDF)
            # This gives P(X ≥ x) where X ~ Hypergeometric(n, k, m)
            p_hyper = hypergeom.sf(x - 1, n, k, m)
            
            if p_hyper > 0 and not np.isnan(p_hyper):
                # Convert to -log10 scale
                score = -np.log10(p_hyper)
                
                # Scale to 1-100 range
                # Typical values range from 1e-50 to 1e-1, giving scores 1-500
                # We scale and cap at 100
                scaled_score = min(100.0, max(1.0, score * 10))
            else:
                scaled_score = 100.0  # Maximum score for p ≈ 0
                
            return scaled_score
            
        except (ValueError, OverflowError, ZeroDivisionError):
            # Handle numerical edge cases
            return 50.0  # Neutral score
    
    def _apply_fdr_correction(self, 
                            results: List[AssociationResult], 
                            alpha: float = 0.01) -> List[AssociationResult]:
        """
        Apply Benjamini-Hochberg FDR correction to statistical results.
        
        Args:
            results: List of association results to correct
            alpha: False discovery rate threshold
            
        Returns:
            List[AssociationResult]: Significant results after FDR correction
        """
        logger.info(f"Applying Benjamini-Hochberg FDR correction with α={alpha}")
        
        if not results:
            return results
        
        # Extract p-values
        p_values = [r.p_value for r in results]
        
        # Apply FDR correction
        rejected, q_values = fdrcorrection(p_values, alpha=alpha, method='indep')
        
        # Update results with q-values
        for result, q_value in zip(results, q_values):
            result.q_value = q_value
        
        # Filter to significant results
        significant_results = [r for r in results if r.q_value < alpha]
        
        logger.info("FDR correction results:")
        logger.info(f"  Input associations: {len(results)}")
        logger.info(f"  Significant after correction: {len(significant_results)}")
        logger.info(f"  Rejection rate: {100 * len(significant_results) / len(results):.2f}%")
        
        # Sort by significance (q-value, then by association score)
        significant_results.sort(key=lambda x: (x.q_value, -x.hyper_score))
        
        return significant_results
    
    def export_results(self, 
                      results: List[AssociationResult], 
                      output_path: Union[str, Path],
                      include_all_stats: bool = True) -> Path:
        """
        Export statistical results to TSV file.
        
        Args:
            results: List of association results to export
            output_path: Path where to save the results file
            include_all_stats: Whether to include all statistical details
            
        Returns:
            Path: Path to the exported file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Exporting {len(results)} results to {output_path}")
        
        # Convert results to DataFrame
        data = []
        for result in results:
            if include_all_stats:
                data.append(result.to_dict())
            else:
                data.append({
                    'domain': result.domain,
                    'go_term': result.go_term,
                    'proteins_both': result.a,
                    'p_value': result.p_value,
                    'q_value': result.q_value,
                    'hyper_score': result.hyper_score,
                    'enrichment_fold': result.enrichment_fold
                })
        
        # Create DataFrame and export
        df = pd.DataFrame(data)
        df.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"Results exported successfully to {output_path}")
        return output_path
    
    def get_summary_statistics(self, results: List[AssociationResult]) -> Dict[str, Union[int, float]]:
        """
        Calculate summary statistics for a set of results.
        
        Args:
            results: List of association results to summarize
            
        Returns:
            Dict: Summary statistics
        """
        if not results:
            return {
                'total_associations': 0,
                'significant_associations': 0,
                'unique_domains': 0,
                'unique_go_terms': 0
            }
        
        significant_results = [r for r in results if r.is_significant()]
        p_values = [r.p_value for r in results]
        q_values = [r.q_value for r in results if r.q_value is not None]
        scores = [r.hyper_score for r in results]
        
        stats = {
            'total_associations': len(results),
            'significant_associations': len(significant_results),
            'unique_domains': len(set(r.domain for r in results)),
            'unique_go_terms': len(set(r.go_term for r in results)),
            'min_p_value': min(p_values) if p_values else None,
            'median_p_value': np.median(p_values) if p_values else None,
            'min_q_value': min(q_values) if q_values else None,
            'median_q_value': np.median(q_values) if q_values else None,
            'mean_hyper_score': np.mean(scores) if scores else None,
            'median_hyper_score': np.median(scores) if scores else None,
            'max_hyper_score': max(scores) if scores else None
        }
        
        return stats


def main():
    """
    Main function for standalone testing of the statistical inference engine.
    
    This function creates sample data and demonstrates the usage of the 
    StatisticalInferenceEngine class.
    """
    # Configure logging
    logger.info("Starting StatisticalInferenceEngine demonstration")
    
    # Create sample data for testing
    sample_protein_domain_map = {
        'P001': ['PF00001', 'PF00002'],
        'P002': ['PF00001', 'PF00003'], 
        'P003': ['PF00002', 'PF00003'],
        'P004': ['PF00001'],
        'P005': ['PF00002'],
        'P006': ['PF00003'],
        'P007': ['PF00001', 'PF00002', 'PF00003'],
        'P008': ['PF00004'],
        'P009': ['PF00004', 'PF00001'],
        'P010': ['PF00005']
    }
    
    sample_protein_go_map = {
        'P001': {'GO:0001', 'GO:0002'},
        'P002': {'GO:0001', 'GO:0003'},
        'P003': {'GO:0002', 'GO:0003'},
        'P004': {'GO:0001'},
        'P005': {'GO:0002'},
        'P006': {'GO:0003'},
        'P007': {'GO:0001', 'GO:0002', 'GO:0003'},
        'P008': {'GO:0004'},
        'P009': {'GO:0004', 'GO:0001'},
        'P010': {'GO:0005'}
    }
    
    try:
        # Initialize engine
        engine = StatisticalInferenceEngine(
            sample_protein_domain_map, 
            sample_protein_go_map
        )
        
        # Run statistical tests
        results = engine.run_statistical_tests(min_cooccurrence=1)
        
        # Display results
        logger.info(f"Found {len(results)} significant associations")
        
        for result in results[:5]:  # Show top 5
            logger.info(f"  {result.domain} - {result.go_term}: "
                       f"p={result.p_value:.2e}, q={result.q_value:.2e}, "
                       f"score={result.hyper_score:.1f}")
        
        # Show summary statistics
        stats = engine.get_summary_statistics(results)
        logger.info("Summary statistics:")
        for key, value in stats.items():
            logger.info(f"  {key}: {value}")
        
        logger.info("StatisticalInferenceEngine demonstration completed successfully")
        
    except Exception as e:
        logger.error(f"Demonstration failed: {str(e)}")
        raise


if __name__ == "__main__":
    main()