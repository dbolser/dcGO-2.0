# src/ontology_processor.py
"""
Complete ontology processing module implementing the True Path Rule and annotation propagation.

This module provides the OntologyProcessor class for handling Gene Ontology (GO) processing
including optimal level filtering, True Path Rule implementation, and hierarchical annotation
propagation as specified in the dcGO methodology.
"""

import networkx as nx
import obonet
import pandas as pd
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Union
from loguru import logger
from dataclasses import dataclass
from scipy.stats import fisher_exact
import gzip


@dataclass
class Annotation:
    """
    Represents a domain-GO term annotation with metadata.
    
    Stores the complete information for a domain-GO term association including
    statistical significance, annotation source, and propagation information.
    
    Attributes:
        domain: Domain identifier (e.g., Pfam ID or supra-domain combination)
        go_term: GO term identifier (e.g., GO:0008150)
        q_value: FDR-corrected p-value from statistical testing
        association_score: Hypergeometric-based association score (1-100 scale)
        annotation_type: Either 'direct' (from statistical inference) or 'propagated'
        direct_source_term: Original GO term that generated this annotation (for propagated terms)
    """
    domain: str
    go_term: str
    q_value: float
    association_score: float
    annotation_type: str  # 'direct' or 'propagated'
    direct_source_term: str

    def __post_init__(self) -> None:
        """Validate annotation data after initialization."""
        if self.annotation_type not in {'direct', 'propagated'}:
            raise ValueError(f"annotation_type must be 'direct' or 'propagated', got: {self.annotation_type}")
        if not (0.0 <= self.q_value <= 1.0):
            raise ValueError(f"q_value must be between 0 and 1, got: {self.q_value}")
        if not (1.0 <= self.association_score <= 100.0):
            raise ValueError(f"association_score must be between 1 and 100, got: {self.association_score}")


class OntologyProcessor:
    """
    Processes Gene Ontology structure for domain-centric functional annotation.
    
    This class implements the dcGO methodology's ontology processing components:
    1. True Path Rule implementation with optimal level filtering
    2. Hierarchical annotation propagation
    3. GO graph management and validation
    
    The processor loads GO ontology files, applies statistical filtering based on
    the True Path Rule, and propagates direct associations up the ontology hierarchy
    to generate comprehensive domain annotations.
    """

    def __init__(self, obo_file: Union[str, Path]) -> None:
        """
        Initialize ontology processor with GO ontology file.
        
        Args:
            obo_file: Path to GO ontology file (.obo format, can be gzipped)
            
        Raises:
            FileNotFoundError: If the ontology file doesn't exist
            ValueError: If the ontology file cannot be parsed
        """
        self.obo_file = Path(obo_file)
        if not self.obo_file.exists():
            raise FileNotFoundError(f"Ontology file not found: {obo_file}")
        
        logger.info(f"Loading GO ontology from {self.obo_file}")
        
        try:
            # Handle gzipped files
            if self.obo_file.suffix == '.gz':
                with gzip.open(self.obo_file, 'rt') as f:
                    self.go_graph = obonet.read_obo(f)
            else:
                self.go_graph = obonet.read_obo(str(self.obo_file))
        except Exception as e:
            raise ValueError(f"Failed to parse ontology file {obo_file}: {e}")
        
        self._prepare_graph()
        self._compute_graph_metrics()

    def _prepare_graph(self) -> None:
        """
        Prepare GO graph for efficient processing.
        
        Removes obsolete terms, validates graph structure, and prepares
        cached data structures for optimal performance.
        """
        logger.info("Preparing GO graph for processing")
        
        # Remove obsolete terms
        obsolete_terms = [
            term for term, data in self.go_graph.nodes(data=True)
            if data.get('is_obsolete', False)
        ]
        
        if obsolete_terms:
            self.go_graph.remove_nodes_from(obsolete_terms)
            logger.info(f"Removed {len(obsolete_terms)} obsolete GO terms")
        
        # Validate graph structure
        if not nx.is_directed_acyclic_graph(self.go_graph):
            logger.warning("GO graph contains cycles - this may indicate parsing issues")
        
        # Cache frequently used structures
        self._ancestors_cache = {}
        self._descendants_cache = {}
        
        logger.info(f"GO graph ready: {len(self.go_graph.nodes)} terms, {len(self.go_graph.edges)} relationships")

    def _compute_graph_metrics(self) -> None:
        """Compute and log basic graph metrics."""
        try:
            # Count terms by aspect
            aspects = {}
            for term, data in self.go_graph.nodes(data=True):
                namespace = data.get('namespace', 'unknown')
                aspects[namespace] = aspects.get(namespace, 0) + 1
            
            logger.info("GO term distribution by namespace:")
            for namespace, count in aspects.items():
                logger.info(f"  {namespace}: {count} terms")
                
            # Compute depth statistics
            root_nodes = [n for n in self.go_graph.nodes() if self.go_graph.in_degree(n) == 0]
            if root_nodes:
                depths = []
                for root in root_nodes:
                    try:
                        paths = nx.single_source_shortest_path_length(self.go_graph, root)
                        depths.extend(paths.values())
                    except Exception:
                        continue
                
                if depths:
                    logger.info(f"Graph depth statistics: max={max(depths)}, avg={sum(depths)/len(depths):.1f}")
                    
        except Exception as e:
            logger.warning(f"Failed to compute graph metrics: {e}")

    def get_ancestors(self, go_term: str) -> Set[str]:
        """
        Get all ancestor terms for a given GO term.
        
        Uses caching for performance optimization.
        
        Args:
            go_term: GO term identifier
            
        Returns:
            Set of ancestor GO term identifiers
        """
        if go_term not in self._ancestors_cache:
            if go_term in self.go_graph:
                self._ancestors_cache[go_term] = set(nx.ancestors(self.go_graph, go_term))
            else:
                self._ancestors_cache[go_term] = set()
        
        return self._ancestors_cache[go_term]

    def get_descendants(self, go_term: str) -> Set[str]:
        """
        Get all descendant terms for a given GO term.
        
        Uses caching for performance optimization.
        
        Args:
            go_term: GO term identifier
            
        Returns:
            Set of descendant GO term identifiers
        """
        if go_term not in self._descendants_cache:
            if go_term in self.go_graph:
                self._descendants_cache[go_term] = set(nx.descendants(self.go_graph, go_term))
            else:
                self._descendants_cache[go_term] = set()
        
        return self._descendants_cache[go_term]

    def apply_optimal_level_filter(self, 
                                 significant_associations: List,
                                 protein_domain_map: Dict[str, List[str]], 
                                 protein_go_map: Dict[str, Set[str]],
                                 min_background_size: int = 10,
                                 alpha_threshold: float = 0.05) -> List:
        """
        Apply optimal level determination filter implementing the True Path Rule.
        
        This method implements Phase 1 of the ontology processing pipeline.
        For each significant domain-GO association, it tests whether the association
        is significantly stronger than associations with parent terms. This ensures
        annotations are made at the most specific (optimal) level in the ontology.
        
        The True Path Rule states that if a protein is annotated to a GO term,
        it is implicitly annotated to all ancestor terms. This method identifies
        associations that are not merely inherited from more general parent terms.
        
        Args:
            significant_associations: List of statistically significant associations
            protein_domain_map: Mapping of protein IDs to domain lists
            protein_go_map: Mapping of protein IDs to GO term sets
            min_background_size: Minimum number of proteins required for background testing
            alpha_threshold: P-value threshold for significance testing
            
        Returns:
            List of associations that pass the optimal level filter
            
        Raises:
            ValueError: If input data is malformed
        """
        logger.info("Applying optimal level filter (True Path Rule implementation)")
        
        if not significant_associations:
            logger.warning("No significant associations provided for filtering")
            return []
        
        # Validate input data
        if not protein_domain_map or not protein_go_map:
            raise ValueError("Protein mapping data cannot be empty")
        
        filtered_associations = []
        total_associations = len(significant_associations)
        
        for i, assoc in enumerate(significant_associations):
            if i % 1000 == 0:
                logger.info(f"Processing association {i+1}/{total_associations}")
            
            try:
                if self._passes_optimal_level_test(
                    assoc.domain, 
                    assoc.go_term, 
                    protein_domain_map, 
                    protein_go_map,
                    min_background_size,
                    alpha_threshold
                ):
                    filtered_associations.append(assoc)
            except Exception as e:
                logger.warning(f"Error testing association {assoc.domain}-{assoc.go_term}: {e}")
                continue
        
        logger.info(f"Optimal level filter: {len(filtered_associations)}/{total_associations} associations retained")
        return filtered_associations

    def _passes_optimal_level_test(self, 
                                 domain: str, 
                                 child_term: str, 
                                 protein_domain_map: Dict[str, List[str]], 
                                 protein_go_map: Dict[str, Set[str]],
                                 min_background_size: int,
                                 alpha_threshold: float) -> bool:
        """
        Test if a domain-GO association is at optimal specificity level.
        
        Implements the core True Path Rule logic by testing whether the
        domain-child term association is significantly stronger than
        domain-parent term associations.
        
        Args:
            domain: Domain identifier
            child_term: Child GO term being tested
            protein_domain_map: Protein to domain mapping
            protein_go_map: Protein to GO term mapping
            min_background_size: Minimum background proteins required
            alpha_threshold: Significance threshold
            
        Returns:
            True if association passes optimal level test
        """
        if child_term not in self.go_graph:
            logger.debug(f"GO term {child_term} not found in ontology, keeping association")
            return True  # Keep if term not in graph (conservative approach)
        
        # Get direct parents (predecessors in the graph)
        parents = list(self.go_graph.predecessors(child_term))
        
        if not parents:
            logger.debug(f"GO term {child_term} has no parents, keeping association")
            return True  # Root terms pass by default
        
        # Test against each parent
        for parent_term in parents:
            try:
                p_value = self._test_against_parent_background(
                    domain, 
                    child_term, 
                    parent_term, 
                    protein_domain_map, 
                    protein_go_map,
                    min_background_size
                )
                
                # If not significantly stronger than any parent, reject
                if p_value >= alpha_threshold:
                    logger.debug(f"Association {domain}-{child_term} not significantly stronger than parent {parent_term} (p={p_value:.4f})")
                    return False
                    
            except Exception as e:
                logger.warning(f"Error testing against parent {parent_term}: {e}")
                # Conservative: if we can't test, reject the association
                return False
        
        return True

    def _test_against_parent_background(self, 
                                      domain: str, 
                                      child_term: str, 
                                      parent_term: str,
                                      protein_domain_map: Dict[str, List[str]], 
                                      protein_go_map: Dict[str, Set[str]],
                                      min_background_size: int) -> float:
        """
        Test domain-child association strength within parent term background.
        
        This implements the statistical test at the core of the True Path Rule.
        We test whether the domain is significantly enriched in proteins annotated
        to the child term, when considering only proteins annotated to the parent term.
        
        Args:
            domain: Domain identifier
            child_term: Child GO term
            parent_term: Parent GO term defining the background
            protein_domain_map: Protein to domain mapping
            protein_go_map: Protein to GO term mapping
            min_background_size: Minimum background size for valid test
            
        Returns:
            P-value from Fisher's exact test
            
        Raises:
            ValueError: If background is too small or data is invalid
        """
        # Get all proteins annotated with parent term (background set)
        parent_proteins = {
            protein for protein, terms in protein_go_map.items()
            if parent_term in terms
        }
        
        if len(parent_proteins) < min_background_size:
            raise ValueError(f"Insufficient background size: {len(parent_proteins)} < {min_background_size}")
        
        # Build 2x2 contingency table within parent background
        # a: proteins with domain AND child term (within parent background)
        a = len([
            p for p in parent_proteins
            if domain in protein_domain_map.get(p, []) and 
               child_term in protein_go_map.get(p, set())
        ])
        
        # b: proteins with child term but NOT domain (within parent background)
        b = len([
            p for p in parent_proteins
            if child_term in protein_go_map.get(p, set()) and 
               domain not in protein_domain_map.get(p, [])
        ])
        
        # c: proteins with domain but NOT child term (within parent background)
        c = len([
            p for p in parent_proteins
            if domain in protein_domain_map.get(p, []) and 
               child_term not in protein_go_map.get(p, set())
        ])
        
        # d: proteins with neither domain nor child term (within parent background)
        d = len(parent_proteins) - (a + b + c)
        
        # Validate contingency table
        if a == 0 or (a + b) == 0 or (a + c) == 0:
            return 1.0  # No association possible
        
        if d < 0:
            raise ValueError(f"Invalid contingency table: a={a}, b={b}, c={c}, d={d}")
        
        try:
            # One-tailed Fisher's exact test for enrichment
            _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
            return p_value
        except (ValueError, ZeroDivisionError) as e:
            logger.warning(f"Fisher's exact test failed: {e}")
            return 1.0

    def propagate_annotations(self, direct_associations: List) -> List[Annotation]:
        """
        Propagate annotations up the ontology hierarchy.
        
        This implements Phase 2 of the ontology processing pipeline.
        For each direct association that passed the optimal level filter,
        propagate the annotation to all ancestor terms in the GO hierarchy.
        
        The propagated annotations inherit the statistical significance and
        association score from the direct source term, but are marked as
        'propagated' to distinguish them from direct statistical associations.
        
        Args:
            direct_associations: List of associations that passed optimal level filtering
            
        Returns:
            List containing both direct and propagated annotations
            
        Raises:
            ValueError: If input associations are malformed
        """
        logger.info("Propagating annotations up ontology hierarchy")
        
        if not direct_associations:
            logger.warning("No direct associations provided for propagation")
            return []
        
        propagated_annotations = []
        processed_pairs = set()  # Track domain-GO pairs to avoid duplicates
        
        for i, assoc in enumerate(direct_associations):
            if i % 1000 == 0:
                logger.info(f"Propagating annotation {i+1}/{len(direct_associations)}")
            
            try:
                # Add direct annotation
                direct_key = (assoc.domain, assoc.go_term)
                if direct_key not in processed_pairs:
                    propagated_annotations.append(Annotation(
                        domain=assoc.domain,
                        go_term=assoc.go_term,
                        q_value=assoc.q_value,
                        association_score=assoc.hyper_score,
                        annotation_type='direct',
                        direct_source_term=assoc.go_term
                    ))
                    processed_pairs.add(direct_key)
                
                # Propagate to ancestors
                if assoc.go_term in self.go_graph:
                    ancestors = self.get_ancestors(assoc.go_term)
                    
                    for ancestor in ancestors:
                        ancestor_key = (assoc.domain, ancestor)
                        if ancestor_key not in processed_pairs:
                            propagated_annotations.append(Annotation(
                                domain=assoc.domain,
                                go_term=ancestor,
                                q_value=assoc.q_value,
                                association_score=assoc.hyper_score,
                                annotation_type='propagated',
                                direct_source_term=assoc.go_term
                            ))
                            processed_pairs.add(ancestor_key)
                else:
                    logger.debug(f"GO term {assoc.go_term} not found in ontology for propagation")
                    
            except Exception as e:
                logger.warning(f"Error propagating annotation for {assoc.domain}-{assoc.go_term}: {e}")
                continue
        
        # Count direct vs propagated
        direct_count = sum(1 for ann in propagated_annotations if ann.annotation_type == 'direct')
        propagated_count = len(propagated_annotations) - direct_count
        
        logger.info(f"Generated {len(propagated_annotations)} total annotations "
                   f"({direct_count} direct, {propagated_count} propagated)")
        
        return propagated_annotations

    def validate_annotations(self, annotations: List[Annotation]) -> Dict[str, int]:
        """
        Validate annotation consistency and compute quality metrics.
        
        Performs comprehensive validation of the annotation set including:
        - Consistency checks between direct and propagated annotations
        - Ontology structure validation
        - Statistical measure validation
        - Duplicate detection
        
        Args:
            annotations: List of annotations to validate
            
        Returns:
            Dictionary containing validation statistics and error counts
        """
        logger.info(f"Validating {len(annotations)} annotations")
        
        validation_stats = {
            'total_annotations': len(annotations),
            'valid_annotations': 0,
            'invalid_go_terms': 0,
            'invalid_scores': 0,
            'invalid_q_values': 0,
            'consistency_errors': 0,
            'duplicate_pairs': 0
        }
        
        seen_pairs = set()
        direct_terms = set()
        propagated_sources = set()
        
        for ann in annotations:
            try:
                # Basic validation
                ann.__post_init__()  # Triggers validation
                
                # Check GO term validity
                if ann.go_term not in self.go_graph:
                    validation_stats['invalid_go_terms'] += 1
                    continue
                
                # Check for duplicates
                pair_key = (ann.domain, ann.go_term)
                if pair_key in seen_pairs:
                    validation_stats['duplicate_pairs'] += 1
                    continue
                seen_pairs.add(pair_key)
                
                # Track direct terms and propagated sources
                if ann.annotation_type == 'direct':
                    direct_terms.add(ann.go_term)
                else:
                    propagated_sources.add(ann.direct_source_term)
                
                validation_stats['valid_annotations'] += 1
                
            except ValueError as e:
                if 'association_score' in str(e):
                    validation_stats['invalid_scores'] += 1
                elif 'q_value' in str(e):
                    validation_stats['invalid_q_values'] += 1
                logger.debug(f"Validation error for {ann.domain}-{ann.go_term}: {e}")
        
        # Check consistency between direct and propagated annotations
        orphaned_propagations = propagated_sources - direct_terms
        validation_stats['consistency_errors'] = len(orphaned_propagations)
        
        if orphaned_propagations:
            logger.warning(f"Found {len(orphaned_propagations)} propagated annotations without direct sources")
        
        logger.info("Validation completed:")
        for metric, value in validation_stats.items():
            logger.info(f"  {metric}: {value}")
        
        return validation_stats

    def get_ontology_statistics(self) -> Dict[str, Union[int, float, List[str]]]:
        """
        Compute comprehensive ontology statistics.
        
        Returns:
            Dictionary containing ontology metrics and statistics
        """
        stats = {
            'total_terms': len(self.go_graph.nodes),
            'total_relationships': len(self.go_graph.edges),
            'namespaces': {},
            'root_terms': [],
            'leaf_terms': [],
            'max_depth': 0,
            'average_depth': 0.0
        }
        
        # Count by namespace
        depths = []
        for term, data in self.go_graph.nodes(data=True):
            namespace = data.get('namespace', 'unknown')
            stats['namespaces'][namespace] = stats['namespaces'].get(namespace, 0) + 1
            
            # Identify root and leaf terms
            if self.go_graph.in_degree(term) == 0:
                stats['root_terms'].append(term)
            if self.go_graph.out_degree(term) == 0:
                stats['leaf_terms'].append(term)
        
        # Compute depth statistics
        try:
            for root in stats['root_terms']:
                paths = nx.single_source_shortest_path_length(self.go_graph, root)
                depths.extend(paths.values())
            
            if depths:
                stats['max_depth'] = max(depths)
                stats['average_depth'] = sum(depths) / len(depths)
        except Exception as e:
            logger.warning(f"Failed to compute depth statistics: {e}")
        
        stats['leaf_term_count'] = len(stats['leaf_terms'])
        stats['root_term_count'] = len(stats['root_terms'])
        
        return stats

    def export_graph_summary(self, output_path: Path) -> Path:
        """
        Export ontology graph summary to TSV file.
        
        Args:
            output_path: Path for output TSV file
            
        Returns:
            Path to the created summary file
        """
        logger.info(f"Exporting ontology summary to {output_path}")
        
        summary_data = []
        for term, data in self.go_graph.nodes(data=True):
            summary_data.append({
                'go_id': term,
                'name': data.get('name', ''),
                'namespace': data.get('namespace', ''),
                'definition': data.get('def', ''),
                'is_obsolete': data.get('is_obsolete', False),
                'in_degree': self.go_graph.in_degree(term),
                'out_degree': self.go_graph.out_degree(term),
                'ancestor_count': len(self.get_ancestors(term)),
                'descendant_count': len(self.get_descendants(term))
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"Exported summary for {len(summary_data)} GO terms")
        return output_path

    def clear_cache(self) -> None:
        """Clear internal caches to free memory."""
        logger.info("Clearing ontology processor caches")
        self._ancestors_cache.clear()
        self._descendants_cache.clear()


# Example usage and testing functions
def create_test_processor(test_obo_content: str) -> OntologyProcessor:
    """
    Create a test ontology processor with minimal GO data.
    
    Args:
        test_obo_content: String content for test .obo file
        
    Returns:
        OntologyProcessor instance for testing
    """
    import tempfile
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.obo', delete=False) as f:
        f.write(test_obo_content)
        temp_path = Path(f.name)
    
    try:
        processor = OntologyProcessor(temp_path)
        return processor
    finally:
        temp_path.unlink()  # Clean up temp file


# Minimal test OBO content for testing
TEST_OBO_CONTENT = """
format-version: 1.2
data-version: test

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0006810
name: transport
namespace: biological_process
is_a: GO:0009987 ! cellular process

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005215
name: transporter activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function
"""


if __name__ == "__main__":
    """Example usage of the OntologyProcessor class."""
    
    # Configure logging
    logger.add("ontology_processor.log", level="INFO", rotation="10 MB")
    
    # Create test processor
    processor = create_test_processor(TEST_OBO_CONTENT)
    
    # Print ontology statistics
    stats = processor.get_ontology_statistics()
    print("Ontology Statistics:")
    for key, value in stats.items():
        if not isinstance(value, (list, dict)):
            print(f"  {key}: {value}")
    
    # Test ancestor/descendant relationships
    test_term = "GO:0006810"  # transport
    ancestors = processor.get_ancestors(test_term)
    print(f"\nAncestors of {test_term}: {ancestors}")
    
    # Clear cache when done
    processor.clear_cache()