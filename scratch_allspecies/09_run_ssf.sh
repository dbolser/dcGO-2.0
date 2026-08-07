#!/usr/bin/env bash
# SCOP-superfamily-keyed runs, human and all-species.
#
# This is the acceptance criterion that actually motivated widening the
# background: against the published dcGO we score precision 0.54-0.63 but
# recall 0.069, and 69.4% of the pairs they call that we miss have ZERO
# co-occurring human proteins - unreachable at any threshold purely because we
# looked at one species. If a multi-species background is worth anything, that
# recall has to move.
#
# --domain-key ssf keys domains by SUPERFAMILY signature (the published dcGO's
# own domain universe), which is ~2,000 domains rather than ~45,000 InterPro
# entries, so these runs are cheap even all-species.
set -euo pipefail

until [ -s data/interim/protein2ipr_allspecies.dat.gz ]; do sleep 60; done
# Do not contend with the InterPro-keyed runs for memory.
while pgrep -f "run_dcgo_human.py --species allspecies --disable-supra" > /dev/null; do
    sleep 60
done

run () {
    local label=$1; shift
    echo "=== ssf run: $label ==="
    uv run python run_dcgo_human.py --domain-key ssf --disable-supra-domains "$@" \
        > "scratch_allspecies/run_ssf_${label}.log" 2>&1
    echo "  done: $label"
}

run human_manual      --species human      --evidence-filter manual \
                      --output-dir results_ssf_human_manual
run allspecies_manual --species allspecies --evidence-filter manual \
                      --output-dir results_ssf_allspecies_manual

echo "SSF RUNS COMPLETE"
