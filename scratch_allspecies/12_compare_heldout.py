#!/usr/bin/env python3
"""Compare two temporal-benchmark runs cell by cell.

The acceptance criterion is directional and specific: held-out enrichment must
not *fall* against the human-only baseline **on the same human evaluation set**.
So only the training universe differs between the two runs being compared here —
both are scored against human t0 (2021) -> t1 (2026), the split
VALIDATION_PLAN §2 already uses.

Prints every aspect x IC cell rather than an average, because the §4 ablation
established that this method's effects are cell-dependent and an average hides
sign flips.
"""

import csv
import sys
from pathlib import Path

BASE = Path("validation/heldout_human_single/temporal_benchmark_metrics.tsv")
WIDE = Path("validation/heldout_allspecies_single/temporal_benchmark_metrics.tsv")


def load(path: Path) -> dict:
    out = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["aspect"], r["min_ic"], r["method"])] = r
    return out


def main() -> int:
    if not WIDE.exists():
        print(f"missing {WIDE} - run the all-species held-out benchmark first")
        return 1
    base, wide = load(BASE), load(WIDE)

    for metric in ("f_max", "auprc"):
        print(f"\n=== {metric} (dcGO), human-trained vs all-species-trained ===")
        print(
            f"{'aspect':<8}{'IC':>5}{'human':>10}{'allspec':>10}{'delta':>10}  verdict"
        )
        wins = losses = 0
        for aspect in ("MF", "BP", "CC"):
            for ic in ("0.0", "2.0", "4.0"):
                b = base.get((aspect, ic, "dcGO"))
                w = wide.get((aspect, ic, "dcGO"))
                if not b or not w:
                    continue
                bv, wv = float(b[metric]), float(w[metric])
                delta = wv - bv
                verdict = "better" if delta > 0 else "WORSE" if delta < 0 else "equal"
                wins += delta > 0
                losses += delta < 0
                print(
                    f"{aspect:<8}{ic:>5}{bv:>10.4f}{wv:>10.4f}{delta:>+10.4f}  {verdict}"
                )
        print(f"  all-species better in {wins}/{wins + losses} cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
