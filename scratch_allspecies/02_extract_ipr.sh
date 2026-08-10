#!/usr/bin/env bash
# Step 2 of the all-species build: pull the protein2ipr rows for our universe.
#
# Same semantics as extract_human_interpro.py's step 2 (keep a row iff its
# column-1 accession is in the protein list), but pigz/awk instead of Python's
# single-threaded gzip, because the input is 13 GB and the protein list is
# ~100x larger than the human one.
set -euo pipefail

LIST=data/interim/allspecies_proteins.txt
SRC=data/raw/interpro_mappings/protein2ipr.dat.gz
OUT=data/interim/protein2ipr_allspecies.dat.gz

test -s "$LIST"
echo "protein list: $(wc -l < "$LIST") accessions"

pigz -dc -p 12 "$SRC" \
  | awk -F'\t' 'NR==FNR { ids[$1]; next } ($1 in ids)' "$LIST" - \
  | pigz -p 24 -c > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"

echo "wrote $OUT: $(stat -c%s "$OUT") bytes"
