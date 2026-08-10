#!/usr/bin/env bash
# Build the all-species non-IEA GAF from goa_uniprot_all.gaf.gz.
#
# Filtering to non-IEA here is lossless for both evidence policies we run:
# the pipeline's `manual` preset is a strict subset of non-IEA (it also drops
# HTP/HDA/HMP/HGI/HEP), and `experimental` is a subset of that. So the pipeline
# still applies its own filter; this step only removes the ~97% IEA bulk so the
# file is tractable.
set -euo pipefail

SRC=data/raw/goa_uniprot_all/goa_uniprot_all.gaf.gz
EXPECTED=11714420656
OUT=data/raw/goa_annotations/goa_allspecies.gaf.gz

# 1. Wait for the download to reach the advertised Content-Length.
until [ -f "$SRC" ] && [ "$(stat -c%s "$SRC")" -ge "$EXPECTED" ]; do
    sleep 20
done
echo "download complete: $(stat -c%s "$SRC") bytes"

# 2. Stream once: keep '!' headers and every non-IEA data line.
#    GAF col 7 (1-based) is the evidence code.
mkdir -p "$(dirname "$OUT")"
pigz -dc -p 8 "$SRC" \
  | awk -F'\t' '/^!/ { if (h < 30) { print; h++ } ; next } $7 != "IEA" { print }' \
  | pigz -p 16 -c > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"

echo "wrote $OUT: $(stat -c%s "$OUT") bytes"
