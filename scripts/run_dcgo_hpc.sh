#!/bin/bash

# SLURM job script for dcGO Pipeline on HPC cluster
# Adjust parameters according to your cluster specifications

#SBATCH --job-name=dcgo_pipeline
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/dcgo_%j.out
#SBATCH --error=logs/dcgo_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your-email@domain.com

# Load required modules (adjust according to your cluster)
module load python/3.9
module load java/11
# module load interproscan/5.67-99.0  # If available as module

# Print job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Number of tasks: $SLURM_NTASKS"
echo "Working directory: $SLURM_SUBMIT_DIR"

# Change to submission directory
cd $SLURM_SUBMIT_DIR

# Activate environment (assuming uv is available)
export PATH="$HOME/.cargo/bin:$PATH"
export UV_PYTHON=$(which python3)

# Create necessary directories
mkdir -p logs results/checkpoints data/{raw,processed,interim}

# Set environment variables for optimal performance
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=1

# Run pipeline with appropriate number of cores
echo "Starting dcGO pipeline with $SLURM_NTASKS cores..."

uv run python -m src.main_pipeline \
    --num-cores $SLURM_NTASKS \
    --output-dir results \
    2>&1 | tee logs/pipeline_${SLURM_JOB_ID}.log

exit_code=${PIPESTATUS[0]}

echo "Pipeline finished at: $(date)"
echo "Exit code: $exit_code"

# Optional: Send notification with results summary
if [ $exit_code -eq 0 ]; then
    echo "dcGO Pipeline completed successfully!"
    # Count output files
    if [ -f "results/dcgo_annotations.tsv" ]; then
        num_associations=$(wc -l < results/dcgo_annotations.tsv)
        echo "Generated $num_associations domain-GO associations"
    fi
else
    echo "dcGO Pipeline failed with exit code: $exit_code"
fi

exit $exit_code