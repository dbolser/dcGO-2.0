#!/usr/bin/env bash
# Reduce idmapping_selected.tab.gz to the four columns the non-independence
# measurement needs, so the 7 GB file is scanned once and never again.
#
# idmapping_selected columns (1-based): 1 UniProtKB-AC, 9 UniRef90,
# 10 UniRef50, 13 NCBI-taxon.
set -euo pipefail

SRC=data/raw/uniprot_idmapping/idmapping_selected.tab.gz
OUT=data/interim/uniref_taxon.tsv.gz
EXPECTED=7066467385

until [ -f "$SRC" ] && [ "$(stat -c%s "$SRC")" -ge "$EXPECTED" ]; do
    sleep 60
done

pigz -dc -p 8 "$SRC" \
  | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $9, $10, $13}' \
  | pigz -p 16 -c > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
echo "wrote $OUT: $(stat -c%s "$OUT") bytes"
