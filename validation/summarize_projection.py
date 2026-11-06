"""Headline tables for a knowledge-projection run (``KNOWLEDGE_PROJECTION.md``).

Reads ``projection_metrics.tsv`` and ``projection_bootstrap.tsv`` written by
``validation/projection_benchmark.py`` and prints, per aspect at one IC floor,
the best old dcGO rung (chosen by KY50, as pre-registered), the comparators,
naive, and the three frozen arms, with paired-bootstrap ratio intervals.

    uv run python validation/summarize_projection.py validation/projection/final
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

#: The frozen arms (§1.6), matched by method-name pattern.
ARMS = {
    "domain-only": r"^proj\[ipr\|s=0\.5\]@",
    "sequence": r"^mix\(proj=1,seq=1\)\[",
    "full": r"^mix\(proj=1,seq=1,lit=1\)\[",
}
SINGLE = {"seq alone": r"^seq\[", "lit alone": r"^lit\["}


def fmt(x: float, nd: int = 0) -> str:
    if pd.isna(x):
        return "–"
    if x == float("inf"):
        return "∞"
    return f"{x:,.{nd}f}"


def ratio(a: float, b: float) -> str:
    if b < 10:
        return "old ≈ 0"
    return f"{a / b:.1f}×"


def summarise(run: Path, min_ic: float = 2.0) -> str:
    metrics = pd.read_csv(run / "projection_metrics.tsv", sep="\t")
    boot_path = run / "projection_bootstrap.tsv"
    boot = pd.read_csv(boot_path, sep="\t") if boot_path.exists() else None
    out = []
    for aspect in ("BP", "MF", "CC"):
        m = metrics[(metrics.aspect == aspect) & (metrics.min_ic == min_ic)]
        if m.empty:
            continue
        old = m[m.method.str.startswith("old[")].sort_values(
            ["ky50", "pooled_auprc", "method"], ascending=False
        )
        best = old.iloc[0]
        if boot is not None:
            ref = boot[(boot.aspect == aspect) & (boot.min_ic == min_ic)].reference
            if not ref.empty:  # the bootstrap's own reference, for consistency
                best = m[m.method == ref.iloc[0]].iloc[0]
        rows = [("best old dcGO", best)]
        for label in ("comparator[interpro2go_t0]", "comparator[iea_t0]", "naive"):
            hit = m[m.method == label]
            if not hit.empty:
                rows.append(
                    (label.removeprefix("comparator[").rstrip("]"), hit.iloc[0])
                )
        for label, pattern in {**SINGLE, **ARMS}.items():
            hit = m[m.method.str.contains(pattern, regex=True)]
            if not hit.empty:
                rows.append((label, hit.iloc[0]))
        n, truth = int(best.n_proteins), int(best.n_truth)
        out.append(
            f"\n**{aspect}** — {n} cohort proteins, {truth:,} truth pairs "
            f"(IC ≥ {min_ic:g}); best old = `{best.method}`\n"
        )
        out.append(
            "| method | KY30 | KY50 [95% CI] | KY70 | KY50 leaf | proteins@50 "
            "| KY50 vs old [CI] | pooled AUPRC | AUPRC vs old [CI] | F_max |"
        )
        out.append("|---|---|---|---|---|---|---|---|---|---|")
        for label, r in rows:
            ci = ratio_ci = auprc_ci = ""
            if boot is not None:
                b = boot[
                    (boot.aspect == aspect)
                    & (boot.min_ic == min_ic)
                    & (boot.method == r.method)
                ]
                k = b[b.metric == "ky50"]
                a = b[b.metric == "pooled_auprc"]
                if not k.empty:
                    k = k.iloc[0]
                    ci = f" [{fmt(k.ci_lo)}, {fmt(k.ci_hi)}]"
                    if label != "best old dcGO":
                        ratio_ci = (
                            f" [{fmt(k.ratio_ci_lo, 1)}, {fmt(k.ratio_ci_hi, 1)}]"
                        )
                if not a.empty and label != "best old dcGO":
                    a = a.iloc[0]
                    auprc_ci = f" [{fmt(a.ratio_ci_lo, 1)}, {fmt(a.ratio_ci_hi, 1)}]"
            vs = "–" if label == "best old dcGO" else ratio(r.ky50, best.ky50)
            avs = (
                "–"
                if label == "best old dcGO"
                else f"{r.pooled_auprc / best.pooled_auprc:.1f}×"
            )
            out.append(
                f"| {label} | {fmt(r.ky30)} | {fmt(r.ky50)}{ci} | {fmt(r.ky70)} "
                f"| {fmt(r.get('ky50_leaf'))} | {fmt(r.proteins50)} | {vs}{ratio_ci} "
                f"| {r.pooled_auprc:.3f} | {avs}{auprc_ci} "
                f"| {fmt(r.get('f_max'), 3)} |"
            )
        frozen = m[m.method.str.contains("|".join(ARMS.values()), regex=True)]
        if (
            "precision_at_frozen" in frozen.columns
            and frozen.precision_at_frozen.notna().any()
        ):
            out.append("\nAt the dev-frozen KY50 cutoff:")
            for _, r in frozen.iterrows():
                arm = next(k for k, v in ARMS.items() if re.search(v, r.method))
                out.append(
                    f"- {arm}: precision {r.precision_at_frozen:.2f}, "
                    f"{int(r.tp_at_frozen):,} correct of {int(r.pred_at_frozen):,}"
                )
    return "\n".join(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run", type=Path)
    p.add_argument("--min-ic", type=float, default=2.0)
    args = p.parse_args(argv)
    print(summarise(args.run, args.min_ic))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
