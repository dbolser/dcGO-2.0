#!/bin/bash

# dcGO Pipeline Installation Script
# This script sets up the dcGO pipeline environment and dependencies

set -e

echo "dcGO Pipeline Installation"
echo "=========================="

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo "Error: Please run this script from the dcGO pipeline root directory"
    exit 1
fi

# Check Python version
echo "Checking Python version..."
python3 --version | grep -E "Python 3\.[89]|Python 3\.1[0-9]" || {
    echo "Error: Python 3.8+ required"
    echo "Please install Python 3.8 or higher"
    exit 1
}

# Install uv if not present
if ! command -v uv &> /dev/null; then
    echo "Installing uv package manager..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    source $HOME/.cargo/env
    echo "uv installed successfully"
else
    echo "uv already installed: $(uv --version)"
fi

# Create and sync virtual environment
echo "Creating virtual environment and installing dependencies..."
uv sync

# Install development dependencies
echo "Installing development dependencies..."
uv sync --extra dev

# Verify installation
echo "Verifying installation..."
if uv run python -c "from config.settings import Config; print('✓ Configuration module loaded successfully')"; then
    echo "✓ Basic imports working"
else
    echo "✗ Import test failed"
    exit 1
fi

# Run basic tests
echo "Running basic tests..."
if uv run pytest tests/test_config.py -v; then
    echo "✓ Basic tests passed"
else
    echo "✗ Tests failed"
    exit 1
fi

echo ""
echo "Installation completed successfully! 🎉"
echo ""
echo "Next steps:"
echo "1. Review configuration in config/settings.py"
echo "2. Run the example: uv run python examples/basic_usage.py"
echo "3. Start the pipeline: uv run python -m src.main_pipeline --help"
echo ""
echo "For HPC usage, see scripts/run_dcgo_hpc.sh"