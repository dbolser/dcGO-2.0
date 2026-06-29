#!/usr/bin/env bash
#
# Comprehensive dcGO Pipeline Comparison
# Runs all meaningful configuration combinations for validation
#

set -e  # Exit on error

CORES=8
SPECIES="human"

echo "=========================================="
echo "  dcGO Pipeline Configuration Comparison"
echo "=========================================="
echo ""
echo "Running 8 configuration combinations..."
echo "Cores: $CORES"
echo "Species: $SPECIES"
echo ""

# Create timestamp for this analysis run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
ANALYSIS_DIR="analysis_${TIMESTAMP}"
mkdir -p "$ANALYSIS_DIR"

# Configuration 1: Baseline (single domains only, no true path)
echo "[1/8] Running BASELINE (single domains only, no true path)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --disable-supra-domains \
  --disable-true-path \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/01_baseline" \
  > "$ANALYSIS_DIR/01_baseline.log" 2>&1 &
PID1=$!

# Configuration 2: True path only (single domains)
echo "[2/8] Running TRUE PATH ONLY (single domains)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --disable-supra-domains \
  --enable-true-path \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/02_true_path_only" \
  > "$ANALYSIS_DIR/02_true_path_only.log" 2>&1 &
PID2=$!

# Configuration 3: Supra-domains only (no true path, no shrinkage)
echo "[3/8] Running SUPRA-DOMAINS ONLY (no true path, no shrinkage)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --enable-supra-domains \
  --disable-true-path \
  --disable-shrinkage \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/03_supra_only" \
  > "$ANALYSIS_DIR/03_supra_only.log" 2>&1 &
PID3=$!

# Configuration 4: Supra-domains + true path (no shrinkage)
echo "[4/8] Running SUPRA-DOMAINS + TRUE PATH (no shrinkage)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --enable-supra-domains \
  --enable-true-path \
  --disable-shrinkage \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/04_supra_true_path" \
  > "$ANALYSIS_DIR/04_supra_true_path.log" 2>&1 &
PID4=$!

# Configuration 5: Supra-domains + shrinkage (no true path)
echo "[5/8] Running SUPRA-DOMAINS + SHRINKAGE (no true path)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --enable-supra-domains \
  --disable-true-path \
  --enable-shrinkage \
  --shrinkage-strength 0.5 \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/05_supra_shrinkage" \
  > "$ANALYSIS_DIR/05_supra_shrinkage.log" 2>&1 &
PID5=$!

# Configuration 6: Full dcGO (supra + true path + shrinkage, strength=0.5)
echo "[6/8] Running FULL DCGO (supra + true path + shrinkage 0.5)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --enable-supra-domains \
  --enable-true-path \
  --enable-shrinkage \
  --shrinkage-strength 0.5 \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/06_full_dcgo_05" \
  > "$ANALYSIS_DIR/06_full_dcgo_05.log" 2>&1 &
PID6=$!

# Configuration 7: Full dcGO with high shrinkage (strength=0.7)
echo "[7/8] Running FULL DCGO HIGH SHRINKAGE (strength=0.7)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --enable-supra-domains \
  --enable-true-path \
  --enable-shrinkage \
  --shrinkage-strength 0.7 \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/07_full_dcgo_07" \
  > "$ANALYSIS_DIR/07_full_dcgo_07.log" 2>&1 &
PID7=$!

# Configuration 8: Full dcGO with low shrinkage (strength=0.3)
echo "[8/8] Running FULL DCGO LOW SHRINKAGE (strength=0.3)..."
uv run python run_dcgo_human.py \
  --species "$SPECIES" \
  --enable-supra-domains \
  --enable-true-path \
  --enable-shrinkage \
  --shrinkage-strength 0.3 \
  --num-cores "$CORES" \
  --output-dir "$ANALYSIS_DIR/08_full_dcgo_03" \
  > "$ANALYSIS_DIR/08_full_dcgo_03.log" 2>&1 &
PID8=$!

echo ""
echo "All 8 configurations launched in background."
echo ""
echo "Monitor progress:"
echo "  tail -f $ANALYSIS_DIR/*.log"
echo ""
echo "Check running processes:"
echo "  ps aux | grep run_dcgo_human"
echo ""
echo "Results will be saved to: $ANALYSIS_DIR/"
echo ""

# Wait for all processes to complete
echo "Waiting for all jobs to complete..."
wait $PID1 && echo "✓ [1/8] Baseline complete"
wait $PID2 && echo "✓ [2/8] True path only complete"
wait $PID3 && echo "✓ [3/8] Supra-domains only complete"
wait $PID4 && echo "✓ [4/8] Supra-domains + true path complete"
wait $PID5 && echo "✓ [5/8] Supra-domains + shrinkage complete"
wait $PID6 && echo "✓ [6/8] Full dcGO (0.5) complete"
wait $PID7 && echo "✓ [7/8] Full dcGO (0.7) complete"
wait $PID8 && echo "✓ [8/8] Full dcGO (0.3) complete"

echo ""
echo "=========================================="
echo "  All configurations complete!"
echo "=========================================="
echo ""
echo "Results directory: $ANALYSIS_DIR/"
echo ""
echo "Compare results:"
echo "  01_baseline                - Single domains only"
echo "  02_true_path_only          - Single domains + ontology propagation"
echo "  03_supra_only              - Multi-domain combinations"
echo "  04_supra_true_path         - Multi-domain + ontology"
echo "  05_supra_shrinkage         - Multi-domain + hierarchical inference"
echo "  06_full_dcgo_05            - Complete methodology (recommended)"
echo "  07_full_dcgo_07            - Complete with high shrinkage"
echo "  08_full_dcgo_03            - Complete with low shrinkage"
echo ""
echo "Key files to compare:"
echo "  - domain_go_associations_top100.tsv   (top predictions)"
echo "  - validation/performance_metrics.tsv   (AUPR, precision)"
echo "  - validation/architecture_contribution.tsv (supra-domain impact)"
echo ""
