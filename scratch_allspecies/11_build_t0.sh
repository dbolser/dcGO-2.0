#!/usr/bin/env bash
# Build the 2021 all-species universe, for the held-out acceptance criterion.
#
# The criterion is: does held-out enrichment FALL against the human-only
# baseline, on the same human evaluation set? So only the training side moves
# to all-species — the evaluation stays human t0 (2021) -> t1 (2026), exactly
# the split VALIDATION_PLAN §2 already uses. Release 205 is the same release the
# human t0 file comes from, so the two training sets are contemporaneous.
#
# protein2ipr is deliberately NOT re-extracted for 2021. The existing human t0
# run does the same thing (protein2ipr_human_t0_2021.dat.gz is a symlink to the
# current file): domain assignments are far more stable than annotations, and
# the pipeline intersects anyway, so the universe is t0-annotated ∩
# current-domain-covered. Re-extracting would change the comparison, not
# improve it.
set -euo pipefail

SRC=data/raw/goa_uniprot_all/goa_uniprot_all.gaf.205.gz
EXPECTED=16592085788
OUT=data/raw/goa_annotations/goa_allspecies_t0_2021.gaf.gz

until [ -f "$SRC" ] && [ "$(stat -c%s "$SRC")" -ge "$EXPECTED" ]; do sleep 60; done
echo "t0 download complete: $(stat -c%s "$SRC") bytes"

pigz -dc -p 8 "$SRC" \
  | awk -F'\t' '/^!/ { if (h < 30) { print; h++ } ; next } $7 != "IEA" { print }' \
  | pigz -p 16 -c > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
echo "wrote $OUT: $(stat -c%s "$OUT") bytes"

ln -sfn "$(pwd)/data/interim/protein2ipr_allspecies.dat.gz" \
        data/interim/protein2ipr_allspecies_t0_2021.dat.gz

# Wait out the InterPro-keyed runs rather than competing with them for memory.
while pgrep -f "run_dcgo_human.py --species allspecies --disable-supra" > /dev/null; do
    sleep 60
done

uv run python run_dcgo_human.py \
    --species allspecies_t0_2021 --disable-supra-domains \
    --evidence-filter manual --output-dir results_allspecies_t0_2021 \
    > scratch_allspecies/run_allspecies_t0.log 2>&1

echo "T0 ALL-SPECIES RUN COMPLETE"
