"""
Basic usage example for dcGO Pipeline

This script demonstrates how to use the dcGO pipeline for domain-centric
Gene Ontology functional annotation analysis.
"""

from pathlib import Path
import sys

# Add src to path to import modules
sys.path.append(str(Path(__file__).parent.parent / "src"))

from config.settings import Config

def main():
    """Basic usage example"""
    print("dcGO Pipeline - Basic Usage Example")
    print("===================================")
    
    # Initialize configuration
    config = Config()
    print(f"Base directory: {config.BASE_DIR}")
    print(f"Data directory: {config.DATA_DIR}")
    print(f"Results directory: {config.RESULTS_DIR}")
    print()
    
    # Display configuration settings
    print("Configuration Settings:")
    print(f"  FDR Threshold: {config.FDR_THRESHOLD}")
    print(f"  Min proteins per association: {config.MIN_PROTEINS_PER_ASSOCIATION}")
    print(f"  Max supra-domain length: {config.MAX_SUPRA_DOMAIN_LENGTH}")
    print(f"  CPU cores: {config.NUM_CORES}")
    print()
    
    # Display data sources
    print("Data Sources:")
    for source, url in config.DATASOURCES.items():
        print(f"  {source}: {url[:60]}...")
    print()
    
    print("To run the complete pipeline:")
    print("  uv run python -m src.main_pipeline")
    print()
    print("Or with custom parameters:")
    print("  uv run python -m src.main_pipeline --num-cores 16 --output-dir /path/to/results")


if __name__ == "__main__":
    main()