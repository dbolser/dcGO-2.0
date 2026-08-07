#!/usr/bin/env bash
# Wait for the filtered all-species GAF, then build the rest of the universe:
# the protein list (via the project's own parser, so the set is identical to
# what the run will see) and the protein2ipr subset for it.
#
# The protein list is built with --evidence-filter manual, which is a superset
# of experimental — one domain extraction then serves both planned runs.
set -euo pipefail

GAF=data/raw/goa_annotations/goa_allspecies.gaf.gz

until [ -s "$GAF" ]; do sleep 20; done
echo "GAF ready: $(stat -c%s "$GAF") bytes"

echo "== step 1: protein list (project parser) =="
uv run python - <<'PY'
from pathlib import Path
from src.goa_parser import parse_goa

gaf = Path("data/raw/goa_annotations/goa_allspecies.gaf.gz")
protein_terms = parse_goa(gaf, evidence_filter="manual", aspects={"P", "F", "C"})
out = Path("data/interim/allspecies_proteins.txt")
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as fh:
    for accession in sorted(protein_terms):
        fh.write(f"{accession}\n")
print(f"wrote {out}: {len(protein_terms):,} accessions")
PY

echo "== step 2: protein2ipr subset =="
bash scratch_allspecies/02_extract_ipr.sh

echo "ALL-SPECIES UNIVERSE BUILD COMPLETE"
