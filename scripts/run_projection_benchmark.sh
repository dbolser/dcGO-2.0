#!/usr/bin/env bash
# Knowledge-projection benchmark (KNOWLEDGE_PROJECTION.md): the frozen
# configuration on one split. Usage: run_projection_benchmark.sh dev|final OUTDIR [extra flags]
set -euo pipefail
split=$1; out=$2; shift 2
H5="data/raw/uniprot_embeddings/UP000005640_9606/per-protein.h5 data/raw/uniprot_embeddings/uniprot_sprot/per-protein.h5"
FROZEN="--strength 0.5 --unit-sets ipr --seq-embeddings $H5 --seq-knn-k 25 --seq-knn-temperature 0.1 \
  --knn-k 25 --knn-temperature 0.1 --mix proj=1,lit=1 proj=1,seq=1 proj=1,seq=1,lit=1 \
  --extra-evidence-filter experimental --iea-comparator --no-noisy-or --min-ic 0 2 --protein-metrics-min-ic 2"
case $split in
  dev)
    SPLIT="--t0-gaf data/raw/goa_archive/goa_human.gaf.181.gz --t1-gaf data/raw/goa_archive/goa_human.gaf.205.gz \
      --interpro data/interim/protein2ipr_human_t0_2018_ipr69.dat.gz --cohort-interpro data/interim/protein2ipr_human_t0_ipr85.dat.gz \
      --old-run-dir results/projection-dev2018/ipr_manual --literature data/processed/literature/medembed_protein_y2017.npz \
      --extra-train-gaf data/raw/goa_annotations/goa_allspecies_manual_2018.gaf.gz \
      --extra-train-interpro data/interim/protein2ipr_allspecies_t0_2018_ipr69.dat.gz \
      --interpro2go-t0 data/raw/interpro_archive/69.0/interpro2go --interpro2go-t1 data/raw/interpro_archive/85.0/interpro2go" ;;
  final)
    SPLIT="--t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \
      --interpro data/interim/protein2ipr_human_t0_ipr85.dat.gz --cohort-interpro data/interim/protein2ipr_human.dat.gz \
      --old-run-dir results/ablation-t0ipr-2026-10-15/ipr_manual --literature data/processed/literature/medembed_protein_y2020.npz \
      --extra-train-gaf data/raw/goa_annotations/goa_allspecies_manual_2021.gaf.gz \
      --extra-train-interpro data/interim/protein2ipr_allspecies_t0_2021_ipr85.dat.gz \
      --interpro2go-t0 data/raw/interpro_archive/85.0/interpro2go --interpro2go-t1 data/raw/interpro2go/interpro2go" ;;
  *) echo "split must be dev or final" >&2; exit 2 ;;
esac
exec uv run --group projection python validation/projection_benchmark.py $SPLIT $FROZEN --output-dir "$out" "$@"
