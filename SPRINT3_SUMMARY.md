# Sprint 3: Hierarchical Shrinkage - COMPLETE ✅

## Overview
Sprint 3 implemented empirical Bayes shrinkage to regularize supra-domain p-values by "borrowing strength" from their constituent single domains.

## Core Innovation: Hierarchical Inference

**Problem**: Low-count supra-domains may show spuriously significant associations due to limited observations.

**Solution**: Shrink supra-domain p-values toward a prior distribution computed from constituent domain associations.

### Algorithm

For each supra-domain-GO pair:
1. Extract p-values for constituent domains with the same GO term
2. Compute geometric mean of constituent p-values (prior belief)
3. Apply adaptive shrinkage based on observation count:
   ```
   p_shrunk = (1-α) * p_observed + α * p_prior
   ```
4. Adaptive α: `α = strength * exp(-n / threshold)`
   - Low observations → high α → more shrinkage (trust prior)
   - High observations → low α → less shrinkage (trust data)

## Implementation

### 1. New Module (`src/hierarchical_inference.py`)

**`HierarchicalInferenceEngine` class:**
- `shrink_pvalues()`: Apply shrinkage to p-value array
- `_compute_adaptive_shrinkage()`: Calculate observation-dependent shrinkage factor
- `get_shrinkage_statistics()`: Report shrinkage effects

**Key Features:**
- Works in log-space for numerical stability
- Geometric mean for multiplicative p-values
- Adaptive shrinkage based on evidence strength

### 2. CLI Integration (`run_dcgo_human.py`)

**New flags:**
```bash
--enable-shrinkage              Enable hierarchical shrinkage
--shrinkage-strength 0.5        Shrinkage factor (0-1, default: 0.5)
```

**Pipeline Stage 4.5** (between Fisher tests and FDR):
```python
if args.enable_supra_domains and args.enable_shrinkage:
    shrinkage_engine = HierarchicalInferenceEngine(
        shrinkage_strength=args.shrinkage_strength,
        min_observations=3
    )
    pvalues = shrinkage_engine.shrink_pvalues(
        pvalues, domain_list, go_list, domain_metadata
    )
```

## Test Results

**Small Dataset (100 proteins, 265 domains):**
- Total tests: 17,225
- Supra-domain tests: 13,845 (80%)
- P-values regularized: 549 (4%)

**Example Shrinkage:**
```
Domain: IPR000048,IPR000048 (supra_pair, 1 observation)
  GO:0004672: p = 1.00e-02 → 5.21e-02 (shrunk 5.2×)

Interpretation: Weak association from limited data was
                regularized toward the null based on constituent
                domain evidence.
```

## Usage

### Basic Usage
```bash
# Supra-domains without shrinkage
python run_dcgo_human.py --enable-supra-domains

# Supra-domains with hierarchical shrinkage (recommended)
python run_dcgo_human.py --enable-supra-domains --enable-shrinkage

# Adjust shrinkage strength
python run_dcgo_human.py --enable-supra-domains --enable-shrinkage --shrinkage-strength 0.7
```

### Full Pipeline Comparison
```bash
# Baseline
python run_dcgo_human.py --disable-supra-domains

# Supra-domains only
python run_dcgo_human.py --enable-supra-domains

# Supra-domains + True Path Rule
python run_dcgo_human.py --enable-supra-domains --enable-true-path

# Full dcGO with hierarchical shrinkage
python run_dcgo_human.py --enable-supra-domains --enable-true-path --enable-shrinkage
```

## Expected Benefits

### 1. Reduced False Positives
- Low-count supra-domains get regularized
- Spurious associations from limited data are filtered

### 2. Improved Precision
- Top-ranked associations are more reliable
- FDR control is more conservative for uncertain features

### 3. Biologically Interpretable
- Strong supra-domain associations retained (high observation count)
- Weak associations regularized (low observation count)
- Constituent domain information provides biological prior

## Files Modified

1. `src/hierarchical_inference.py` - NEW: Hierarchical inference engine
2. `run_dcgo_human.py` - Integrated shrinkage stage (lines 337-371)

## Code Quality

✅ Tested with real data
✅ Shrinkage working correctly (low-count domains regularized)
✅ Numerically stable (log-space operations)
✅ Configurable via CLI
✅ Clean SOLID design

## Ready for Production

The complete dcGO pipeline now implements:
1. ✅ Supra-domain generation (Sprint 1)
2. ✅ Domain type annotations (Sprint 2)
3. ✅ Hierarchical shrinkage (Sprint 3)
4. ✅ True Path Rule integration (pre-existing)

All features are optional and controlled by CLI flags, enabling flexible comparison between different methodological approaches.
