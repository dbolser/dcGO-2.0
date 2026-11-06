"""Figure 1 of the knowledge-projection paper: pooled precision against yield.

Reads ``validation/projection/final_curves_dense/projection_curves.tsv.gz`` (written by
``validation/projection_benchmark.py --save-curves`` with the frozen
configuration) and draws one panel per GO aspect at IC >= 2: how many correct
informative predictions each method hands over (x) at what pooled precision
(y), for old dcGO's best variant and the three new arms, plus (dotted) the
domain-unit arm trained on human annotations only, from
``validation/projection/final/background_origin_curves.tsv``
(``paper/background_origin_ablation.py``).

    uv run python paper/make_projection_figures.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
#: Validated categorical slots 1-3 (dataviz reference palette, light mode);
#: the old method is the grey reference series.
SERIES = [
    ("old dcGO (best variant)", r"^old\[", "#898781", (0, (4, 2))),
    ("domain units", r"^proj\[ipr\|s=0\.5\]@", "#2a78d6", "solid"),
    ("+ sequence", r"^mix\(proj=1,seq=1\)\[", "#eb6834", "solid"),
    ("+ sequence + literature", r"^mix\(proj=1,seq=1,lit=1\)\[", "#1baf7a", "solid"),
    (
        "domain units, human training only",
        r"^domain units, human training only$",
        "#2a78d6",
        (0, (2, 2)),
    ),
]
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
TITLES = {
    "BP": "Biological process",
    "MF": "Molecular function",
    "CC": "Cellular component",
}


def draw(curves: pd.DataFrame, out: Path, min_ic: float = 2.0) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), facecolor=SURFACE)
    for ax, aspect in zip(axes, ("BP", "MF", "CC")):
        sub = curves[(curves.aspect == aspect) & (curves.min_ic == min_ic)]
        ax.set_facecolor(SURFACE)
        xmax = 0.0
        for label, pattern, color, style in SERIES:
            hit = sub[sub.method.str.contains(pattern, regex=True)]
            if hit.empty:
                continue
            name = hit.method.iloc[0]
            c = hit[hit.method == name].sort_values("tp")
            keep = c[c.precision >= 0.2]
            xmax = max(xmax, float(keep.tp.max()) if not keep.empty else 0.0)
            ax.plot(c.tp, c.precision, color=color, lw=2, ls=style, label=label)
            # KY50 and KY70: the largest yield still at or above each level.
            for level in (0.5, 0.7):
                ok = c[c.precision >= level - 1e-12]
                if not ok.empty and ok.tp.max() > 0:
                    pt = ok.loc[ok.tp.idxmax()]
                    hollow = style != "solid" and label.endswith("only")
                    ax.plot(
                        pt.tp, pt.precision, "o", ms=6,
                        color=SURFACE if hollow else color,
                        markeredgecolor=color if hollow else SURFACE,
                        markeredgewidth=1.5, zorder=5,
                    )  # fmt: skip
        for level in (0.5, 0.7):
            ax.axhline(level, color=MUTED, lw=0.8, ls=":", zorder=0)
            ax.text(
                xmax * 1.12,
                level,
                f"{level:.0%}",
                color=MUTED,
                fontsize=8,
                va="center",
                ha="right",
            )
        ax.set_xlim(0, xmax * 1.15 if xmax else 1)
        ax.set_ylim(0, 1.0)
        ax.set_title(TITLES[aspect], color=INK, fontsize=11, loc="left")
        ax.set_xlabel(
            "correct informative predictions (yield)", color=MUTED, fontsize=9
        )
        ax.grid(True, color=GRID, lw=0.6)
        ax.tick_params(colors=MUTED, labelsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
    axes[0].set_ylabel("pooled precision", color=MUTED, fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, -0.02),
        labelcolor=INK,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--curves",
        type=Path,
        default=ROOT
        / "validation/projection/final_curves_dense/projection_curves.tsv.gz",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=ROOT / "paper/figures/projection_precision_yield.png",
    )
    args = p.parse_args(argv)
    curves = pd.read_csv(args.curves, sep="\t")
    human = ROOT / "validation/projection/final/background_origin_curves.tsv"
    if human.exists():
        curves = pd.concat([curves, pd.read_csv(human, sep="\t")], ignore_index=True)
    # One old rung per cell is written (that cell's KY50 reference).
    olds = curves[curves.method.str.match(r"^old\[")]
    per_cell = olds.groupby(["aspect", "min_ic"]).method.nunique()
    if (per_cell != 1).any():
        raise SystemExit(f"expected one old reference per cell, got {per_cell}")
    draw(curves, args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
