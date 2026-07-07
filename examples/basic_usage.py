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

    print("To run the human dcGO analysis (see README / QUICKSTART.md):")
    print("  1. uv run python scripts/download_data.py")
    print("  2. uv run python extract_human_interpro.py")
    print("  3. uv run python run_dcgo_human.py --num-cores 8")
    print()
    print("Enable GO annotation propagation (True Path Rule):")
    print(
        "  uv run python run_dcgo_human.py --enable-true-path "
        "--go-ontology data/raw/go_ontology/go-basic.obo"
    )


if __name__ == "__main__":
    main()
