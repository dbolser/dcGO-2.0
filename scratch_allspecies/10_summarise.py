#!/usr/bin/env python3
"""Collect every run in this experiment into one comparison table.

Deliberately reports the significant count *next to* the size of the hypothesis
space that produced it. A wider universe grows the count arithmetically, so the
raw number is not evidence of anything on its own — the TODO's acceptance
criteria say so explicitly. What is comparable across runs is the rate: how many
of the tested pairs survived FDR, and how that compares to the permutation null
run on the same universe.
"""

import json
import sys
from pathlib import Path

RUNS = [
    ("human", "manual", "results_human_single_manual"),
    ("human", "experimental", "results_human_single_experimental"),
    ("allspecies", "manual", "results_allspecies_manual"),
    ("allspecies", "experimental", "results_allspecies_experimental"),
    ("allspecies", "manual (PERMUTED null)", "results_allspecies_permuted"),
    ("human ssf", "manual", "results_ssf_human_manual"),
    ("allspecies ssf", "manual", "results_ssf_allspecies_manual"),
]

HEAD = (
    f"{'universe':<16} {'evidence':<22} {'proteins':>10} {'domains':>9} "
    f"{'terms':>8} {'tests':>16} {'signif':>10} {'per 1e6':>9} {'min':>7}"
)


def main() -> int:
    print(HEAD)
    print("-" * len(HEAD))
    rows = []
    for universe, evidence, directory in RUNS:
        manifest = Path(directory) / "run_manifest_go.json"
        if not manifest.exists():
            print(f"{universe:<16} {evidence:<22} {'(not run)':>10}")
            continue
        s = json.load(manifest.open())["summary"]
        rate = s["significant_associations"] / s["tests"] * 1e6
        print(
            f"{universe:<16} {evidence:<22} {s['proteins']:>10,} {s['domains']:>9,} "
            f"{s['terms']:>8,} {s['tests']:>16,} "
            f"{s['significant_associations']:>10,} {rate:>9.1f} "
            f"{s['runtime_seconds'] / 60:>7.1f}"
        )
        rows.append((universe, evidence, s, rate))

    # The only comparison the acceptance criteria allow on counts alone.
    null = next((r for r in rows if "PERMUTED" in r[1]), None)
    real = next(
        (r for r in rows if r[0] == "allspecies" and r[1] == "manual"),
        None,
    )
    if null and real:
        print()
        print(
            f"Calibration: the all-species universe yields "
            f"{real[2]['significant_associations']:,} associations against "
            f"{null[2]['significant_associations']:,} under permutation "
            f"({real[3]:.1f} vs {null[3]:.1f} per million tests)."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
