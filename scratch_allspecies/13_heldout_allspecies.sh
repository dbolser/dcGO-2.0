#!/usr/bin/env bash
# The held-out acceptance criterion, all-species arm.
#
# Only the TRAINING universe differs from validation/heldout_human_single: the
# evaluation is human t0 (2021) -> t1 (2026) in both, so the two runs are
# directly comparable and the criterion ("held-out enrichment does not fall
# against the human-only baseline on the same human evaluation set") is
# answerable cell by cell.
#
# --interpro stays the human domain file: it defines which proteins the
# predictions are transferred onto for scoring, and the evaluation set is human.
set -euo pipefail

until grep -q "T0 ALL-SPECIES RUN COMPLETE" scratch_allspecies/11_build_t0.log 2>/dev/null; do
    sleep 60
done

uv run python validation/temporal_benchmark.py \
    --t0-gaf data/raw/goa_annotations/goa_human_t0_2021.gaf.gz \
    --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
    --predictions results_allspecies_t0_2021/domain_go_associations_significant.tsv \
    --interpro data/interim/protein2ipr_human.dat.gz \
    --disable-supra-domains --n-permutations 20 \
    --output-dir validation/heldout_allspecies_single \
    > scratch_allspecies/bench_allspecies.log 2>&1

echo "HELD-OUT BENCHMARK COMPLETE"
uv run python scratch_allspecies/12_compare_heldout.py
