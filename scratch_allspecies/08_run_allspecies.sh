#!/usr/bin/env bash
# Run the all-species background once the domain subset is extracted.
#
# Single domains only, deliberately: the Fisher engine tests the DENSE
# domain x term product, and the supra-domain space over a 1.4M-protein
# universe does not fit in memory at any sane setting. 07_probe_scale.py
# prints the arithmetic; the runs below are what is actually reachable.
#
# Sequential, not concurrent: each run peaks around 60-70 GB and this is a
# shared machine.
set -euo pipefail

IPR=data/interim/protein2ipr_allspecies.dat.gz
until [ -s "$IPR" ]; do sleep 60; done
echo "domain subset ready: $(stat -c%s "$IPR") bytes"

echo "=== scale probe ==="
uv run python scratch_allspecies/07_probe_scale.py 2>&1 | tee scratch_allspecies/07_probe.log

run () {
    local label=$1; shift
    echo "=== run: $label ==="
    uv run python run_dcgo_human.py \
        --species allspecies --disable-supra-domains "$@" \
        > "scratch_allspecies/run_allspecies_${label}.log" 2>&1
    echo "  done: $label"
}

# The two evidence policies: 'experimental' is the non-circular claim,
# 'manual' is what the existing human baselines used.
run manual       --evidence-filter manual       --output-dir results_allspecies_manual
run experimental --evidence-filter experimental --output-dir results_allspecies_experimental

# Calibration control. A wider universe changes what "random" means, so the
# permutation null has to be re-run rather than inherited from the human runs.
run permuted     --evidence-filter manual       --output-dir results_allspecies_permuted \
                 --permute-annotations 7

echo "ALL-SPECIES RUNS COMPLETE"
