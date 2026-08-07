#!/usr/bin/env bash
# Fetch idmapping_selected.tab.gz once the GOA download is out of the way.
#
# We need it for the non-independence measurement the TODO's acceptance criteria
# demand: columns 8/9/10 are UniRef100/90/50 cluster ids, which give an
# ortholog-group proxy so pooled protein counts can be reported next to
# cluster-collapsed ones. Column 13 is the NCBI taxon.
set -euo pipefail

GOA=data/raw/goa_uniprot_all/goa_uniprot_all.gaf.gz
GOA_EXPECTED=11714420656
OUT=data/raw/uniprot_idmapping/idmapping_selected.tab.gz

until [ -f "$GOA" ] && [ "$(stat -c%s "$GOA")" -ge "$GOA_EXPECTED" ]; do
    sleep 30
done

mkdir -p "$(dirname "$OUT")"
curl -sL -C - -o "$OUT" \
  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
echo "idmapping done: $(stat -c%s "$OUT") bytes"
