#!/usr/bin/env bash
#
# Run the pipeline across its configuration space, for side-by-side comparison.
#
# Each configuration writes a full result set (plus its provenance manifest)
# into its own directory under a timestamped analysis directory, so the
# configurations can be diffed against each other afterwards.
#
# NOTE: this produces association tables, not evaluation metrics. To score the
# configurations against held-out annotation, feed the resulting
# `domain_go_associations_significant.tsv` files to
# `validation/temporal_benchmark.py --predictions ...`.
#
set -euo pipefail

CORES="${CORES:-8}"
SPECIES="${SPECIES:-human}"
# Full runs are memory-hungry (contingency-table construction dominates), so
# they are throttled rather than all launched at once as this script used to.
JOBS="${JOBS:-2}"
GO_OBO="${GO_OBO:-data/raw/go_ontology/go-basic.obo}"

# Configuration name, then the flags that define it. The pipeline's defaults are
# supra-domains ON and True Path OFF, so "disabling" True Path means omitting
# its flag — there is deliberately no --disable-true-path option, and an earlier
# version of this script passed one, which made every configuration fail
# immediately.
CONFIGS=(
  "01_baseline|--disable-supra-domains"
  "02_true_path_only|--disable-supra-domains --enable-true-path --go-ontology ${GO_OBO}"
  "03_supra_only|"
  "04_supra_true_path|--enable-true-path --go-ontology ${GO_OBO}"
)

INTERPRO="data/interim/protein2ipr_${SPECIES}.dat.gz"
GAF="data/raw/goa_annotations/goa_${SPECIES}.gaf.gz"

echo "=========================================="
echo "  dcGO Pipeline Configuration Comparison"
echo "=========================================="
echo "  Species:     ${SPECIES}"
echo "  Cores/run:   ${CORES}"
echo "  Concurrency: ${JOBS}"
echo "  Configs:     ${#CONFIGS[@]}"
echo ""

# Check the inputs up front. Discovering a missing 20 GB extract after an hour
# of running seven other configurations is a waste of everybody's afternoon.
missing=0
for path in "${INTERPRO}" "${GAF}"; do
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: missing required input: ${path}" >&2
    missing=1
  fi
done
if [[ ! -f "${GO_OBO}" ]]; then
  echo "ERROR: missing GO ontology: ${GO_OBO}" >&2
  echo "       (needed by the --enable-true-path configurations)" >&2
  missing=1
fi
if (( missing )); then
  echo "" >&2
  echo "Fetch them with: uv run python scripts/download_data.py" >&2
  echo "then:            uv run python extract_human_interpro.py --species ${SPECIES}" >&2
  exit 1
fi

ANALYSIS_DIR="analysis_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${ANALYSIS_DIR}"
echo "Output directory: ${ANALYSIS_DIR}/"
echo ""

# Failures are recorded as marker files, not as a shell variable. Each
# configuration runs in a backgrounded subshell, so an increment inside one is
# lost when it exits — an earlier version of this counted that way and cheerfully
# reported "All 8 configurations complete" after all 8 had failed.
FAILURE_DIR="${ANALYSIS_DIR}/.failures"
mkdir -p "${FAILURE_DIR}"

run_config() {
  local name="$1" flags="$2"
  local log="${ANALYSIS_DIR}/${name}.log"
  local associations="${ANALYSIS_DIR}/${name}/domain_go_associations_significant.tsv"
  echo "  [start] ${name}"
  # shellcheck disable=SC2086  # flags are intentionally word-split
  if uv run python run_dcgo_human.py \
       --species "${SPECIES}" \
       --num-cores "${CORES}" \
       --output-dir "${ANALYSIS_DIR}/${name}" \
       ${flags} > "${log}" 2>&1 && [[ -f "${associations}" ]]; then
    echo "  [ok]    ${name} — $(( $(wc -l < "${associations}") - 1 )) significant associations"
  else
    echo "  [FAIL]  ${name} — see ${log}" >&2
    touch "${FAILURE_DIR}/${name}"
  fi
}

running=0
for entry in "${CONFIGS[@]}"; do
  name="${entry%%|*}"
  flags="${entry#*|}"
  run_config "${name}" "${flags}" &
  running=$(( running + 1 ))
  if (( running >= JOBS )); then
    wait -n
    running=$(( running - 1 ))
  fi
done
wait

failures=$(find "${FAILURE_DIR}" -type f | wc -l)
rmdir "${FAILURE_DIR}" 2>/dev/null || true

echo ""
echo "=========================================="
if (( failures )); then
  echo "  Completed with ${failures} of ${#CONFIGS[@]} configuration(s) FAILED"
else
  echo "  All ${#CONFIGS[@]} configurations complete"
fi
echo "=========================================="
echo ""
echo "Results: ${ANALYSIS_DIR}/"
echo ""
echo "  01_baseline          Single domains only"
echo "  02_true_path_only    Single domains + GO DAG propagation"
echo "  03_supra_only        + contiguous domain combinations (the default)"
echo "  04_supra_true_path   + GO DAG propagation"
echo ""
echo "Per configuration:"
echo "  domain_go_associations_significant.tsv   all associations at FDR < 0.01"
echo "  domain_go_associations_top100.tsv        strongest 100"
echo "  run_manifest_go.json                     inputs, hashes, thresholds, git commit"
echo ""
echo "To score them against held-out annotation, pass each"
echo "significant-associations file to validation/temporal_benchmark.py."
echo ""

exit $(( failures > 0 ))
