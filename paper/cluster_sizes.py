"""Sizes of the evidence clusters used by the cluster bootstrap (test split).

Cohort proteins (IC >= 2 truth) are linked when experimental t1 annotations
of each, in the aspect, cite the same paper (``evidence_clusters``). Writes
``validation/projection/final_rescored/cluster_sizes.tsv``.

    uv run --group projection python paper/cluster_sizes.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from paper.projection_training_sizes import SPLITS, pb  # noqa: E402


def main() -> int:
    args = pb.parse_args(SPLITS["final"])
    inputs = pb.ab.load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        "interpro",
        args.evidence_filter,
        cohort_interpro=args.cohort_interpro,
    )
    rows = []
    for aspect in pb.ASPECTS:
        truth = pb.tb.filter_by_ic(inputs.benchmark[aspect], inputs.ic, 2.0)
        proteins = sorted(p for p, ts in truth.items() if ts)
        labels = pb.evidence_clusters(args.t1_gaf, aspect, proteins)
        sizes = np.bincount(np.unique(labels, return_inverse=True)[1])
        rows.append(
            {
                "aspect": aspect,
                "proteins": len(proteins),
                "clusters": len(sizes),
                "singletons": int((sizes == 1).sum()),
                "largest": int(sizes.max()),
                "proteins_in_shared_clusters": int(sizes[sizes > 1].sum()),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(
        ROOT / "validation/projection/final_rescored/cluster_sizes.tsv",
        sep="\t",
        index=False,
    )
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
