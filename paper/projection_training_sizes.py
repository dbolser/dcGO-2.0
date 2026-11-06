"""Training-set sizes behind the knowledge-projection runs (paper Table 1).

Rebuilds, with the driver's own functions and the frozen configuration's
inputs, the per-aspect training sets (human t0 ``manual`` plus the
all-species ``experimental`` background, queries and same-gene accessions
excluded) and writes their sizes to ``<out>/training_sizes.tsv``.

    uv run --group projection python paper/projection_training_sizes.py dev
    uv run --group projection python paper/projection_training_sizes.py final
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location(
    "projection_benchmark", ROOT / "validation" / "projection_benchmark.py"
)
pb = importlib.util.module_from_spec(_spec)
sys.modules["projection_benchmark"] = pb
_spec.loader.exec_module(pb)

from src import knowledge_projection as kp  # noqa: E402

SPLITS = {
    "dev": [
        "--t0-gaf", "data/raw/goa_archive/goa_human.gaf.181.gz",
        "--t1-gaf", "data/raw/goa_archive/goa_human.gaf.205.gz",
        "--interpro", "data/interim/protein2ipr_human_t0_2018_ipr69.dat.gz",
        "--cohort-interpro", "data/interim/protein2ipr_human_t0_ipr85.dat.gz",
        "--extra-train-gaf", "data/raw/goa_annotations/goa_allspecies_manual_2018.gaf.gz",
        "--extra-train-interpro",
        "data/interim/protein2ipr_allspecies_t0_2018_ipr69.dat.gz",
        "--output-dir", "validation/projection/dev2018_frozen",
    ],
    "final": [
        "--t0-gaf", "data/raw/goa_archive/goa_human.gaf.205.gz",
        "--t1-gaf", "data/raw/goa_annotations/goa_human.gaf.gz",
        "--interpro", "data/interim/protein2ipr_human_t0_ipr85.dat.gz",
        "--cohort-interpro", "data/interim/protein2ipr_human.dat.gz",
        "--extra-train-gaf", "data/raw/goa_annotations/goa_allspecies_manual_2021.gaf.gz",
        "--extra-train-interpro",
        "data/interim/protein2ipr_allspecies_t0_2021_ipr85.dat.gz",
        "--output-dir", "validation/projection/final",
    ],
}  # fmt: skip


def main() -> int:
    split = sys.argv[1]
    args = pb.parse_args(SPLITS[split] + ["--extra-evidence-filter", "experimental"])
    inputs = pb.ab.load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        "interpro",
        args.evidence_filter,
        cohort_interpro=args.cohort_interpro,
    )
    families = pb.unit_families(
        inputs.protein_domains, args.interpro, set(inputs.protein_domains)
    )
    human_trainings = kp.build_aspect_trainings(
        inputs.t0_map,
        inputs.get_ancestors,
        inputs.term_aspect,
        pb.ASPECTS,
        exclude={a: list(inputs.benchmark[a]) for a in pb.ASPECTS},
    )
    merged = pb.add_background(args, inputs, families)
    trainings = kp.build_aspect_trainings(
        merged,
        inputs.get_ancestors,
        inputs.term_aspect,
        pb.ASPECTS,
        exclude={a: list(inputs.benchmark[a]) for a in pb.ASPECTS},
    )
    rows = []
    for a in pb.ASPECTS:
        rows.append(
            {
                "split": split,
                "aspect": a,
                "human_training_proteins": len(human_trainings[a].proteins),
                "training_proteins_with_background": len(trainings[a].proteins),
                "training_terms_with_background": len(trainings[a].terms),
                "cohort_proteins_ic0": len(inputs.benchmark[a]),
            }
        )
    out = Path(args.output_dir) / "training_sizes.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
