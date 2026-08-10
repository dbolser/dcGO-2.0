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

import argparse
import csv
import sys
from pathlib import Path

DEFAULT_BASE = Path("validation/heldout_human_single/temporal_benchmark_metrics.tsv")
DEFAULT_WIDE = Path(
    "validation/heldout_allspecies_single/temporal_benchmark_metrics.tsv"
)


def load(path: Path) -> dict:
    out = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["aspect"], r["min_ic"], r["method"])] = r
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--wide", type=Path, default=DEFAULT_WIDE)
    args = ap.parse_args()

    for path in (args.base, args.wide):
        if not path.exists():
            print(f"missing {path} - run that held-out benchmark first")
            return 1
    base, wide = load(args.base), load(args.wide)

    # The evaluation is fixed inside temporal_benchmark.py (t0 parsed 'manual'
    # for IC and the naive baseline, t1 'experimental' as the gold standard) and
    # does not depend on what the predictions were trained on. So the naive rows
    # must reproduce across any two arms; if they do not, the two runs are not
    # scored on the same evaluation and no cell-by-cell delta below is valid.
    drift = [
        (k, c)
        for k in base
        if k[2] == "naive" and k in wide
        for c in ("f_max", "auprc", "n_eval_proteins")
        if base[k][c] != wide[k][c]
    ]
    if drift:
        print(
            f"WARNING: evaluation differs between the two runs in {len(drift)} "
            f"naive cells - the comparison below is not like-for-like"
        )
        for k, c in drift[:5]:
            print(f"  {k[0]} IC>={k[1]} {c}: {base[k][c]} vs {wide[k][c]}")
    else:
        print("evaluation held fixed: naive f_max/auprc/n_eval identical in all cells")

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
