#!/usr/bin/env bash
# The held-out test again, with projected annotation excluded from TRAINING.
#
# The manual-evidence held-out arm showed all-species winning 8/9 F_max and 9/9
# AUPRC cells. But that background is 75.8% projected annotation, over half of
# which cites a human protein — so while the held-out result is the one thing
# circularity cannot fake (future human annotation is not in 2021 training data
# under any evidence code), the SIZE of the win was measured on a background
# largely made of copied human function. This asks how big it is without that.
#
# The evaluation does not move: temporal_benchmark.py fixes it internally (t0
# parsed 'manual' for IC and the naive baseline, t1 'experimental' as the gold
# standard) independent of what the predictions were trained on. So the naive
# rows must come out identical to the manual arms — 12_compare_heldout.py
# checks exactly that before reporting any delta.
set -euo pipefail

run () {
    local species=$1 out=$2
    echo "=== training run: $out ==="
    uv run python run_dcgo_human.py \
        --species "$species" --disable-supra-domains \
        --evidence-filter experimental --output-dir "$out" \
        > "scratch_allspecies/run_${out}.log" 2>&1
    echo "  done: $out"
}

bench () {
    local preds=$1 out=$2
    echo "=== benchmark: $out ==="
    uv run python validation/temporal_benchmark.py \
        --t0-gaf data/raw/goa_annotations/goa_human_t0_2021.gaf.gz \
        --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
        --predictions "$preds" \
        --interpro data/interim/protein2ipr_human.dat.gz \
        --disable-supra-domains --n-permutations 20 \
        --output-dir "$out" \
        > "scratch_allspecies/bench_$(basename "$out").log" 2>&1
    echo "  done: $out"
}

run human_t0_2021      results_human_t0_exp
run allspecies_t0_2021 results_allspecies_t0_exp

bench results_human_t0_exp/domain_go_associations_significant.tsv \
      validation/heldout_human_single_exp
bench results_allspecies_t0_exp/domain_go_associations_significant.tsv \
      validation/heldout_allspecies_single_exp

echo "EXPERIMENTAL HELD-OUT COMPLETE"
uv run python scratch_allspecies/12_compare_heldout.py \
    --base validation/heldout_human_single_exp/temporal_benchmark_metrics.tsv \
    --wide validation/heldout_allspecies_single_exp/temporal_benchmark_metrics.tsv
