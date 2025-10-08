# dcGO Pipeline

A comprehensive pipeline for protein function prediction using dcGO (domain-centric Gene Ontology).

## Overview

The dcGO Pipeline is a bioinformatics tool designed to predict protein functions by leveraging domain-centric approaches to Gene Ontology annotations. This pipeline integrates multiple data sources and machine learning techniques to provide accurate and reliable protein function predictions.

## Features

- **Domain-centric approach**: Utilizes protein domain information for enhanced prediction accuracy
- **Multi-source integration**: Combines data from various biological databases
- **Machine learning models**: Implements state-of-the-art ML algorithms for function prediction
- **Comprehensive validation**: Includes extensive testing and validation frameworks
- **Visualization tools**: Provides rich visualizations for results analysis
- **High performance**: Optimized for speed and scalability
- **Modular design**: Easily extensible and customizable

## Requirements

- Python 3.12+
- CUDA-compatible GPU (optional, for accelerated training)

## Installation

### Using uv (recommended)

```bash
# Install uv if you haven't already
pip install uv

# Clone the repository
git clone https://github.com/yourusername/dcgo-pipeline.git
cd dcgo-pipeline

# Install dependencies
uv sync

# Install in development mode
uv pip install -e .
```

### Using pip

```bash
pip install dcgo-pipeline
```

## Quick Start

```python
from dcgo_pipeline import DCGOPipeline

# Initialize the pipeline
pipeline = DCGOPipeline()

# Load protein sequences
sequences = pipeline.load_sequences("path/to/sequences.fasta")

# Run prediction
results = pipeline.predict(sequences)

# Visualize results
pipeline.visualize_results(results, output_dir="results/")
```

## Project Structure

```
dcGOspeed/
├── src/
│   └── dcgo_pipeline/
│       ├── core/           # Core pipeline functionality
│       ├── utils/          # Utility functions
│       ├── validation/     # Validation and testing utilities
│       ├── integration/    # External tool integrations
│       └── visualization/  # Plotting and visualization
├── config/                 # Configuration files
├── data/
│   ├── raw/               # Raw input data
│   ├── processed/         # Processed data
│   └── interim/           # Intermediate processing results
├── results/               # Output results
├── logs/                  # Log files
├── scripts/               # Utility scripts
├── tests/
│   ├── unit/              # Unit tests
│   ├── integration/       # Integration tests
│   └── e2e/               # End-to-end tests
└── examples/
    ├── basic/             # Basic usage examples
    └── advanced/          # Advanced usage examples
```

## Usage

### Command Line Interface

```bash
# Basic prediction
dcgo predict --input sequences.fasta --output results.json

# With custom configuration
dcgo predict --input sequences.fasta --config config/custom.yaml --output results.json

# Batch processing
dcgo batch --input-dir data/sequences/ --output-dir results/

# Training custom models
dcgo train --data training_data.csv --model-type neural_network --output models/
```

### Python API

```python
import dcgo_pipeline as dcgo

# Load configuration
config = dcgo.load_config("config/default.yaml")

# Initialize pipeline with custom settings
pipeline = dcgo.DCGOPipeline(config=config)

# Process multiple sequences
sequences = ["MKLLVLSLSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQH"]
predictions = pipeline.predict_batch(sequences)

# Access detailed results
for seq_id, prediction in predictions.items():
    print(f"Sequence: {seq_id}")
    print(f"GO terms: {prediction.go_terms}")
    print(f"Confidence: {prediction.confidence}")
```

## Configuration

The pipeline uses YAML configuration files for customization. Key configuration sections include:

- **Data sources**: Database connections and file paths
- **Model parameters**: ML model hyperparameters
- **Processing options**: Parallel processing settings
- **Output formats**: Result formatting preferences

Example configuration:

```yaml
# config/default.yaml
data:
  sequence_db: "data/sequences.fasta"
  domain_db: "data/domains.hmm"
  go_annotations: "data/go_annotations.gaf"

model:
  type: "neural_network"
  hidden_layers: [512, 256, 128]
  dropout: 0.2
  learning_rate: 0.001

processing:
  batch_size: 32
  num_workers: 4
  use_gpu: true

output:
  format: "json"
  include_confidence: true
  threshold: 0.5
```

## Development

### Setting up development environment

```bash
# Clone repository
git clone https://github.com/yourusername/dcgo-pipeline.git
cd dcgo-pipeline

# Install development dependencies
uv sync --group dev

# Install pre-commit hooks
pre-commit install

# Run tests
pytest

# Run with coverage
pytest --cov=src/dcgo_pipeline
```

### Code Quality

The project uses several tools to maintain code quality:

- **Black**: Code formatting
- **isort**: Import sorting
- **flake8**: Linting
- **mypy**: Type checking
- **pytest**: Testing

Run all quality checks:

```bash
# Format code
black src/ tests/
isort src/ tests/

# Lint code
flake8 src/ tests/

# Type checking
mypy src/

# Run tests
pytest
```

## Testing

The project includes comprehensive tests at multiple levels:

```bash
# Run all tests
pytest

# Run specific test categories
pytest -m unit          # Unit tests only
pytest -m integration   # Integration tests only
pytest -m e2e           # End-to-end tests only

# Run with coverage
pytest --cov=src/dcgo_pipeline --cov-report=html
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{dcgo_pipeline,
  title={dcGO Pipeline: A comprehensive pipeline for protein function prediction},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/dcgo-pipeline}
}
```

## Support

- **Documentation**: [https://dcgo-pipeline.readthedocs.io](https://dcgo-pipeline.readthedocs.io)
- **Issues**: [GitHub Issues](https://github.com/yourusername/dcgo-pipeline/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/dcgo-pipeline/discussions)

## Acknowledgments

- dcGO methodology and original implementation
- Gene Ontology Consortium
- UniProt database
- InterPro domain database
- All contributors and supporters of this project