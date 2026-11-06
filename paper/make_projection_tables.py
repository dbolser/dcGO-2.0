"""Fill the knowledge-projection manuscript's tables from the result files.

Every table in ``paper/knowledge_projection.md`` sits between
``<!-- BEGIN table:NAME -->`` and ``<!-- END table:NAME -->`` and is rewritten
from the TSVs under ``validation/projection/``, so prose and data cannot
drift. Run after any re-run:

    uv run python paper/make_projection_tables.py
    uv run python paper/make_projection_tables.py --check   # fail if stale

Old dcGO is compared metric by metric against its own best variant for that
metric (KY50's is the pre-registered reference); every ratio is read from
the bootstrap file together with the reference it was computed against.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "validation" / "projection"
PAPER = ROOT / "paper" / "knowledge_projection.md"
#: Test-split runs with per-metric old references in the bootstrap. Same
#: frozen configuration and inputs as the one-shot ``final*`` runs, whose
#: metrics they reproduce (paper Section 4.7).
#: The ``*_v2`` runs add the difference with the reference re-chosen in each
#: resample and the fixed-cutoff precision of old variants.
FINAL = "final_v2"
FINAL_GENEKNOWN = "final_geneknown_v2"
FINAL_HTFILTER = "final_htfilter_v2"
FINAL_STRICT = "final_strictnk_v2"
DEV = "dev2018_frozen"
ASPECTS = ("BP", "MF", "CC")
MIN_PROTEINS = 10  # reliability cells on fewer distinct proteins are not shown
ARMS = {
    "domain units": r"^proj\[ipr\|s=0\.5\]@",
    "+ sequence": r"^mix\(proj=1,seq=1\)\[",
    "+ sequence + literature": r"^mix\(proj=1,seq=1,lit=1\)\[",
}
# The primary tables also show the arm named in the original pre-registration.
PRIMARY_ARMS = {
    "domain units": ARMS["domain units"],
    "+ literature (pre-registered arm)": r"^mix\(proj=1,lit=1\)\[",
    "+ sequence": ARMS["+ sequence"],
    "+ sequence + literature": ARMS["+ sequence + literature"],
}
CHANNELS = {
    "domain units": r"^proj\[ipr\|s=0\.5\]",
    "sequence neighbours": r"^seq\[",
    "literature neighbours": r"^lit\[",
    "domain + literature": r"^mix\(proj=1,lit=1\)\[",
    "domain + sequence": r"^mix\(proj=1,seq=1\)\[",
    "domain + sequence + literature": r"^mix\(proj=1,seq=1,lit=1\)\[",
}


# --------------------------------------------------------------------------- #
# Loading                                                                      #
# --------------------------------------------------------------------------- #
def metrics(run: str) -> pd.DataFrame:
    return pd.read_csv(RUNS / run / "projection_metrics.tsv", sep="\t")


def boot(run: str, resampling: str = "proteins") -> pd.DataFrame:
    b = pd.read_csv(RUNS / run / "projection_bootstrap.tsv", sep="\t")
    if "resampling" in b.columns:
        b = b[b.resampling == resampling]
    return b


def cell(frame: pd.DataFrame, aspect: str, min_ic: float = 2.0) -> pd.DataFrame:
    return frame[(frame.aspect == aspect) & (frame.min_ic == min_ic)]


def pick(frame: pd.DataFrame, pattern: str) -> pd.Series:
    hit = frame[frame.method.str.contains(pattern, regex=True)]
    if len(hit) != 1:
        raise ValueError(f"{pattern!r} matched {len(hit)} methods: {list(hit.method)}")
    return hit.iloc[0]


def best_old(m: pd.DataFrame, metric: str = "ky50") -> pd.Series:
    """The old variant best on ``metric`` (ties: AUPRC, or KY50 for AUPRC)."""
    second = "ky50" if metric == "pooled_auprc" else "pooled_auprc"
    old = m[m.method.str.startswith("old[")]
    return old.sort_values([metric, second, "method"], ascending=False).iloc[0]


def reference(
    run: str, aspect: str, metric: str = "ky50", min_ic: float = 2.0
) -> pd.Series:
    """The old variant the bootstrap compared ``metric`` against.

    Read from the bootstrap file and checked against :func:`best_old`, so a
    table can never pair a ratio with the wrong reference.
    """
    m = cell(metrics(run), aspect, min_ic)
    b = cell(boot(run), aspect, min_ic)
    refs = b[b.metric == metric].reference.unique()
    if len(refs) != 1:
        raise ValueError(f"{run} {aspect} {metric}: references {refs}")
    ref = m[m.method == refs[0]].iloc[0]
    if ref.method != best_old(m, metric).method:
        raise ValueError(f"{run} {aspect} {metric}: bootstrap reference is not best")
    return ref


def boot_row(
    run: str, aspect: str, method: str, metric: str, resampling: str = "proteins"
) -> pd.Series:
    b = cell(boot(run, resampling), aspect)
    return b[(b.method == method) & (b.metric == metric)].iloc[0]


# --------------------------------------------------------------------------- #
# Formatting                                                                   #
# --------------------------------------------------------------------------- #
def diff_cell(point: float, k: pd.Series) -> str:
    """KY50 difference with its re-chosen-reference interval."""
    return (
        f"{point:+,.0f}".replace("-", "−")
        + " "
        + ci(k.diff_reselect_ci_lo, k.diff_reselect_ci_hi)
    )


def n(x: float) -> str:
    return f"{x:,.0f}"


def ci(lo: float, hi: float, nd: int = 0) -> str:
    def f(v):
        return "∞" if v == float("inf") else f"{v:,.{nd}f}".replace("-", "−")

    return f"[{f(lo)}, {f(hi)}]"


def ratio(a: float, b: float) -> str:
    return "old < 10" if b < 10 else f"{a / b:.1f}×"


def rung(method: str) -> str:
    """``old[supra_input|p|pscore]`` -> ``supra_input`` (p, p-score)."""
    m = re.match(r"old\[(.+)\|(q|p)\|(\w+)\]", method)
    if not m:
        return method
    transfer = {"pscore": "p-score", "max": "max"}[m.group(3)]
    return f"`{m.group(1)}` ({m.group(2)}, {transfer})"


def table(header: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# Tables                                                                       #
# --------------------------------------------------------------------------- #
def t_cohorts() -> str:
    rows = []
    for run, sizes_run, label in (
        (DEV, DEV, "development (2018 → 2021)"),
        (FINAL, "final", "test (2021 → 2026)"),
    ):
        sizes = pd.read_csv(RUNS / sizes_run / "training_sizes.tsv", sep="\t")
        m = metrics(run)
        for a in ASPECTS:
            s = sizes[sizes.aspect == a].iloc[0]
            c = cell(m, a).iloc[0]
            rows.append(
                [
                    label,
                    a,
                    n(s.human_training_proteins),
                    n(s.training_proteins_with_background),
                    n(s.cohort_proteins_ic0),
                    n(c.n_proteins),
                    n(c.n_truth),
                ]
            )
    return table(
        [
            "split",
            "aspect",
            "training proteins, human",
            "+ all-species background",
            "cohort proteins",
            "cohort proteins, IC ≥ 2",
            "truth pairs, IC ≥ 2",
        ],
        rows,
    )


def _channel_rows(frames, with_old: pd.DataFrame) -> str:
    rows = [
        ["old dcGO, best variant for each column"]
        + [
            v
            for a in ASPECTS
            for v in (
                n(best_old(cell(with_old, a), "ky50").ky50),
                f"{best_old(cell(with_old, a), 'pooled_auprc').pooled_auprc:.3f}",
            )
        ],
        ["naive (term frequency)"]
        + [
            v
            for a in ASPECTS
            for r in [cell(with_old, a)[cell(with_old, a).method == "naive"].iloc[0]]
            for v in (n(r.ky50), f"{r.pooled_auprc:.3f}")
        ],
    ]
    for tag, frame in frames:
        for label, pattern in CHANNELS.items():
            hits = [
                cell(frame, a)[cell(frame, a).method.str.contains(pattern, regex=True)]
                for a in ASPECTS
            ]
            if any(len(h) != 1 for h in hits):
                continue
            vals = []
            for h in hits:
                r = h.iloc[0]
                vals += [n(r.ky50), f"{r.pooled_auprc:.3f}"]
            rows.append([f"{label}{tag}"] + vals)
    header = ["method"] + [f"{a} {k}" for a in ASPECTS for k in ("KY50", "AUPRC")]
    return table(header, rows)


def t_dev_components() -> str:
    human, full = metrics("dev2018_humanonly"), metrics(DEV)
    return _channel_rows(
        [(" — human training", human), (" — + all-species background", full)], full
    )


def t_final_channels() -> str:
    final = metrics(FINAL)
    return _channel_rows([("", final)], final)


def t_final_primary() -> str:
    rows = []
    m = metrics(FINAL)
    for a in ASPECTS:
        ref_k = reference(FINAL, a, "ky50")
        ref_p = reference(FINAL, a, "pooled_auprc")
        rb = boot_row(FINAL, a, ref_k.method, "ky50")
        rows.append(
            [a, f"old dcGO, best by KY50: {rung(ref_k.method)}",
             f"{n(ref_k.ky50)} {ci(rb.ci_lo, rb.ci_hi)}", "—", "—", "—",
             f"{ref_k.pooled_auprc:.3f}", "—", "—"]
        )  # fmt: skip
        if ref_p.method != ref_k.method:
            rows.append(
                [a, f"old dcGO, best by AUPRC: {rung(ref_p.method)}", n(ref_p.ky50),
                 "—", "—", "—", f"{ref_p.pooled_auprc:.3f}", "—", "—"]
            )  # fmt: skip
        for label, pattern in PRIMARY_ARMS.items():
            r = pick(cell(m, a), pattern)
            k = boot_row(FINAL, a, r.method, "ky50")
            p = boot_row(FINAL, a, r.method, "pooled_auprc")
            rows.append(
                [
                    a,
                    label,
                    f"{n(r.ky50)} {ci(k.ci_lo, k.ci_hi)}",
                    f"{ratio(r.ky50, ref_k.ky50)} {ci(k.ratio_ci_lo, k.ratio_ci_hi, 1)}",
                    ci(k.ratio_reselect_ci_lo, k.ratio_reselect_ci_hi, 1),
                    diff_cell(r.ky50 - ref_k.ky50, k),
                    f"{r.pooled_auprc:.3f}",
                    f"{r.pooled_auprc / ref_p.pooled_auprc:.1f}× "
                    f"{ci(p.ratio_reselect_ci_lo, p.ratio_reselect_ci_hi, 1)}",
                    ci(p.ratio_naive_ci_lo, p.ratio_naive_ci_hi, 1),
                ]
            )
    return table(
        [
            "aspect",
            "method",
            "KY50 [95% CI]",
            "KY50 ratio to best old by KY50 [CI, reference fixed]",
            "same ratio, reference re-chosen in each resample [CI]",
            "KY50 difference to best old [CI, reference re-chosen]",
            "pooled AUPRC",
            "AUPRC ratio to best old by AUPRC [CI, reference re-chosen]",
            "AUPRC ratio to naive [CI]",
        ],
        rows,
    )


def t_uncertainty() -> str:
    """The primary contrast (full method, KY50) under four resampling schemes."""
    rows = []
    m = metrics(FINAL)
    arms = {
        "full (named in the frozen addendum)": ARMS["+ sequence + literature"],
        "pre-registered: domain units": ARMS["domain units"],
        "pre-registered: + literature": PRIMARY_ARMS[
            "+ literature (pre-registered arm)"
        ],
    }
    for a, (arm, pattern) in ((a, x) for a in ASPECTS for x in arms.items()):
        r = pick(cell(m, a), pattern)
        ref = reference(FINAL, a, "ky50")
        for scheme, label in (
            ("proteins", "proteins"),
            ("clusters", "evidence clusters"),
        ):
            k = boot_row(FINAL, a, r.method, "ky50", scheme)
            p = boot_row(FINAL, a, r.method, "pooled_auprc", scheme)
            rows.append(
                [
                    a,
                    arm,
                    label,
                    f"{ratio(r.ky50, ref.ky50)} {ci(k.ratio_ci_lo, k.ratio_ci_hi, 1)}",
                    ci(k.ratio_reselect_ci_lo, k.ratio_reselect_ci_hi, 1),
                    f"{k.p_ratio_reselect_ge10:.2f}",
                    diff_cell(r.ky50 - ref.ky50, k),
                    ci(p.ratio_reselect_ci_lo, p.ratio_reselect_ci_hi, 1),
                ]
            )
    return table(
        [
            "aspect",
            "arm",
            "resampling unit",
            "KY50 ratio [CI, reference fixed]",
            "KY50 ratio [CI, reference re-chosen]",
            "share of resamples with ratio ≥ 10 (re-chosen)",
            "KY50 difference [CI, reference re-chosen]",
            "AUPRC ratio [CI, reference re-chosen]",
        ],
        rows,
    )


def _frozen_ci(method: str, aspect: str) -> tuple[float, float]:
    r = boot_row(FINAL, aspect, method, "precision_at_frozen")
    return r.ci_lo, r.ci_hi


def _old_frozen(aspect: str) -> str:
    """Old dcGO's best variant *on dev* (by KY50), at its dev cutoff, on test."""
    dev = best_old(cell(metrics(DEV), aspect), "ky50")
    if not dev.ky50:
        return "no cutoff (every old variant had dev KY50 = 0)"
    r = cell(metrics(FINAL), aspect)
    r = r[r.method == dev.method].iloc[0]
    if not r.pred_at_frozen or pd.isna(r.pred_at_frozen):
        return f"{rung(dev.method)}: nothing scored above its cutoff"
    b = cell(boot(FINAL), aspect)
    b = b[(b.method == r.method) & (b.metric == "precision_at_frozen")]
    if r.pred_at_frozen < 2:
        interval = ""
    elif len(b):
        interval = f" {ci(b.iloc[0].ci_lo, b.iloc[0].ci_hi, 2)}"
    else:
        interval = " (not bootstrapped)"
    return (
        f"{rung(dev.method)}: {r.precision_at_frozen:.2f}{interval} "
        f"({n(r.tp_at_frozen)} of {n(r.pred_at_frozen)})"
    )


def t_final_secondary() -> str:
    rows = []
    m = metrics(FINAL)
    for a in ASPECTS:
        c = cell(m, a)
        ref = reference(FINAL, a, "ky50")
        r30, r70 = best_old(c, "ky30"), best_old(c, "ky70")
        rows.append(
            [
                a,
                "old dcGO, best variant for each column",
                n(r30.ky30),
                n(r70.ky70),
                n(ref.ky50_leaf),
                n(ref.proteins50),
                _old_frozen(a),
            ]
        )
        for label, pattern in ARMS.items():
            r = pick(c, pattern)
            rows.append(
                [
                    a,
                    label,
                    n(r.ky30),
                    n(r.ky70),
                    n(r.ky50_leaf),
                    n(r.proteins50),
                    f"{r.precision_at_frozen:.2f} "
                    f"{ci(*_frozen_ci(r.method, a), 2)} ({n(r.tp_at_frozen)} of "
                    f"{n(r.pred_at_frozen)})",
                ]
            )
    return table(
        [
            "aspect",
            "method",
            "KY30",
            "KY70",
            "KY50, most specific only",
            "proteins with ≥ 1 correct at KY50",
            "precision at the cutoff fixed on dev [CI] (correct of predicted)",
        ],
        rows,
    )


def t_fmax() -> str:
    rows = []
    m = metrics(FINAL)
    for a in ASPECTS:
        c = cell(m, a)
        c = c[c.f_max.notna()]
        old = (
            c[c.method.str.startswith("old[")].sort_values(["f_max", "method"]).iloc[-1]
        )
        naive = c[c.method == "naive"].iloc[0]
        rows.append(
            [a, f"{naive.f_max:.3f}", f"{old.f_max:.3f} ({rung(old.method)})"]
            + [f"{pick(c, p).f_max:.3f}" for p in ARMS.values()]
        )
    return table(
        ["aspect", "naive", "old dcGO, best p-score variant by F_max"] + list(ARMS),
        rows,
    )


def t_companion() -> str:
    """The pre-registered companion summaries (KNOWLEDGE_PROJECTION.md §1.2).

    The IC-weighted yield is read at the KY50 threshold, so old dcGO's value
    is its KY50 reference's; protein-centric AUPRC and S_min take the best
    old p-score variant on each (the only old variants with those metrics).
    """
    rows = []
    m = metrics(FINAL)
    for a in ASPECTS:
        c = cell(m, a)
        ref = reference(FINAL, a, "ky50")
        scored = c[c.method.str.startswith("old[") & c.auprc.notna()]
        a_best = scored.sort_values(["auprc", "method"]).iloc[-1]
        s_best = scored.sort_values(["s_min", "method"]).iloc[0]
        rows.append(
            [
                a,
                "old dcGO (see caption)",
                f"{ref.bits50:,.0f}",
                f"{a_best.auprc:.3f}",
                f"{s_best.s_min:.1f}",
            ]
        )
        for label, pattern in ARMS.items():
            r = pick(c, pattern)
            rows.append(
                [a, label, f"{r.bits50:,.0f}", f"{r.auprc:.3f}", f"{r.s_min:.1f}"]
            )
    return table(
        [
            "aspect",
            "method",
            "IC-weighted yield at KY50 (bits)",
            "protein-centric AUPRC",
            "S_min (bits; lower is better)",
        ],
        rows,
    )


def t_sensitivity() -> str:
    rows = []
    runs = (
        ("primary cohort", FINAL),
        ("without gene-known proteins", FINAL_GENEKNOWN),
        ("without high-throughput papers", FINAL_HTFILTER),
        ("no curated annotation in any aspect", FINAL_STRICT),
    )
    full = ARMS["+ sequence + literature"]
    for a in ASPECTS:
        for label, run in runs:
            ref_k = reference(run, a, "ky50")
            ref_p = reference(run, a, "pooled_auprc")
            r = pick(cell(metrics(run), a), full)
            k = boot_row(run, a, r.method, "ky50")
            p = boot_row(run, a, r.method, "pooled_auprc")
            rows.append(
                [
                    a,
                    label,
                    n(r.n_proteins),
                    f"{n(ref_k.ky50)} ({rung(ref_k.method)})",
                    n(r.ky50),
                    f"{ratio(r.ky50, ref_k.ky50)} {ci(k.ratio_ci_lo, k.ratio_ci_hi, 1)}",
                    ci(k.ratio_reselect_ci_lo, k.ratio_reselect_ci_hi, 1),
                    diff_cell(r.ky50 - ref_k.ky50, k),
                    f"{r.pooled_auprc / ref_p.pooled_auprc:.1f}× "
                    f"{ci(p.ratio_reselect_ci_lo, p.ratio_reselect_ci_hi, 1)}",
                ]
            )
    return table(
        [
            "aspect",
            "cohort",
            "proteins (IC ≥ 2)",
            "old KY50 (best variant in this cohort)",
            "new KY50",
            "KY50 ratio [CI, reference fixed]",
            "KY50 ratio [CI, reference re-chosen]",
            "KY50 difference [CI, reference re-chosen]",
            "AUPRC ratio to best old by AUPRC [CI, reference re-chosen]",
        ],
        rows,
    )


def t_comparators() -> str:
    rows = []
    m = metrics(FINAL)
    for a in ASPECTS:
        c = cell(m, a)
        for label, method in (
            ("InterPro2GO at t0", "comparator[interpro2go_t0]"),
            ("all IEA annotations at t0", "comparator[iea_t0]"),
            ("naive", "naive"),
        ):
            r = c[c.method == method].iloc[0]
            rows.append(
                [a, label, n(r.ky30), n(r.ky50), f"{r.pooled_auprc:.3f}",
                 f"{r.coverage_any:.2f}"]
            )  # fmt: skip
        rows.append(
            [
                a,
                "old dcGO, best variant for each column",
                n(best_old(c, "ky30").ky30),
                n(best_old(c, "ky50").ky50),
                f"{best_old(c, 'pooled_auprc').pooled_auprc:.3f}",
                f"{best_old(c, 'coverage_any').coverage_any:.2f}",
            ]
        )
        r = pick(c, ARMS["+ sequence + literature"])
        rows.append([a, "+ sequence + literature", n(r.ky30), n(r.ky50),
                     f"{r.pooled_auprc:.3f}", f"{r.coverage_any:.2f}"])  # fmt: skip
    return table(
        ["aspect", "predictor", "KY30", "KY50", "pooled AUPRC", "proteins scored"],
        rows,
    )


def t_ic0() -> str:
    """The same comparison without the information-content floor."""
    rows = []
    m = metrics(FINAL)
    for a in ASPECTS:
        c = cell(m, a, 0.0)
        rows.append(
            [a, "old dcGO, best variant for each column", n(c.n_truth.iloc[0]),
             n(best_old(c, "ky50").ky50),
             f"{best_old(c, 'pooled_auprc').pooled_auprc:.3f}"]
        )  # fmt: skip
        r = c[c.method == "naive"].iloc[0]
        rows.append([a, "naive", n(r.n_truth), n(r.ky50), f"{r.pooled_auprc:.3f}"])
        for label, pattern in ARMS.items():
            r = pick(c, pattern)
            rows.append([a, label, n(r.n_truth), n(r.ky50), f"{r.pooled_auprc:.3f}"])
    return table(
        ["aspect", "method", "truth pairs (all terms)", "KY50", "pooled AUPRC"], rows
    )


def t_association() -> str:
    rows = []
    for split, run in (("development", "dev2018_assoc"), ("test", "final_assoc")):
        for ids, suffix in (("as released", ""), ("current GO ids", "_live_ids")):
            a = pd.read_csv(
                RUNS / run / f"projection_association_level{suffix}.tsv", sep="\t"
            )
            prec = a[a.metric == "precision_at_budget"].set_index("rung")
            hits = a[a.metric == "n_hit_at_budget"].set_index("rung").value
            added = int(a[a.metric == "n_added_pairs"].value.iloc[0])
            gained = int(a[a.metric == "n_gained_domains"].value.iloc[0])
            budget = int(a.budget.iloc[0])
            old = prec[~prec.index.str.startswith("proj[")]
            best = old.value.idxmax()

            def fmt(r):
                x = prec.loc[r]
                return f"{x.value:.3f} {ci(x.ci_lo, x.ci_hi, 3)}, {n(hits[r])} hits"

            lb = prec.loc["proj[ipr|s=0.5,q=0.1,leaf]"]
            rows.append(
                [
                    split,
                    ids,
                    f"{n(added)} ({n(gained)} entries)",
                    n(budget),
                    fmt("proj[ipr|s=0.5,q=0.1,leaf]"),
                    fmt("proj[ipr|s=0.5]"),
                    fmt("supra_input_relative"),
                    f"{lb.delta:+.3f} {ci(lb.delta_ci_lo, lb.delta_ci_hi, 3)}",
                    f"{fmt(best)} (`{best}`)",
                ]
            )
    return table(
        [
            "split",
            "GO ids",
            "pairs added t0 → t1 (entries that gained one)",
            "budget",
            "new: lower-bound leaves",
            "new: posterior mean",
            "old: `supra_input_relative` (pre-registered reference)",
            "paired difference, lower-bound leaves − reference [95% CI]",
            "old: best variant",
        ],
        rows,
    )


def t_ceiling() -> str:
    rows = []
    for split, run in (("development", DEV), ("test", FINAL)):
        for a in ASPECTS:
            ref = best_old(cell(metrics(run), a), "ky50")
            # Protocol §1.4: a baseline under 10 is "≈ 0", not a ratio.
            if ref.ky50 == 0:
                cap_text = "unbounded (baseline 0)"
            elif ref.ky50 < 10:
                cap_text = (
                    f"not reported (baseline < 10; {ref.n_truth / ref.ky50:.1f}× "
                    "arithmetically)"
                )
            else:
                cap_text = f"{ref.n_truth / ref.ky50:.1f}×"
            rows.append([split, a, n(ref.n_truth), n(ref.ky50), cap_text])
    return table(
        [
            "split",
            "aspect",
            "truth pairs (IC ≥ 2)",
            "old KY50",
            "largest possible ratio",
        ],
        rows,
    )


def t_curation_lag() -> str:
    rows = []
    for split, run in (("development", DEV), ("test", "final")):
        lag = pd.read_csv(RUNS / run / "curation_lag.tsv", sep="\t")
        lag = lag[lag["view"] == "ic2"]
        for _, r in lag.iterrows():
            rows.append(
                [
                    split,
                    r.aspect,
                    f"≤ {int(r.cutoff_year)}",
                    n(r.evidence_lines),
                    f"{r.share_lines_pre_cutoff:.0%}",
                    n(r.distinct_papers),
                    f"{r.share_papers_pre_cutoff:.0%}",
                ]
            )
    return table(
        [
            "split",
            "aspect",
            "literature cutoff",
            "evidence lines behind the scored truth",
            "lines citing a paper from before the cutoff",
            "distinct papers cited",
            "papers from before the cutoff",
        ],
        rows,
    )


RELIABILITY_BANDS = ["0.5–0.7", "0.7–0.9", "≥ 0.9"]


def _rel_cell(hit: pd.DataFrame) -> str:
    if hit.empty:
        return "none"
    if hit.iloc[0].proteins < MIN_PROTEINS:
        return "—"
    r = hit.iloc[0]
    if r.correct in (0, r.predictions):  # a percentile interval is degenerate
        return (
            f"{n(r.correct)} of {n(r.predictions)} correct "
            f"({n(r.proteins)} proteins; no interval)"
        )
    return (
        f"{r.precision:.2f} {ci(r.ci_lo, r.ci_hi, 2)} "
        f"({n(r.predictions)}; {n(r.proteins)})"
    )


def t_reliability_origin() -> str:
    """Precision of domain-unit scores by where the winning unit's carriers are from."""
    origins = ["human only", "human and other species", "other species only"]
    rows = []
    for split, run in (("dev", DEV), ("test", "final")):
        rel = pd.read_csv(RUNS / run / "reliability_by_origin.tsv", sep="\t")
        for a in ASPECTS:
            sub = rel[rel.aspect == a]
            for band in RELIABILITY_BANDS:
                rows.append(
                    [split, a, band]
                    + [
                        _rel_cell(sub[(sub.score_band == band) & (sub.origin == o)])
                        for o in origins
                    ]
                )
    return table(
        ["split", "aspect", "domain-unit score"] + [f"carriers: {o}" for o in origins],
        rows,
    )


def t_reliability() -> str:
    """Precision of domain-unit scores by score band and the winning unit's support."""
    supports = ["1–3", "4–10", "11–50", "> 50"]
    rows = []
    for label, run in (
        ("dev, + background", DEV),
        ("test, human only", "final_humanonly"),
        ("test, + background", "final"),
    ):
        rel = pd.read_csv(RUNS / run / "reliability_by_support.tsv", sep="\t")
        for a in ASPECTS:
            sub = rel[rel.aspect == a]
            for band in RELIABILITY_BANDS:
                rows.append(
                    [label, a, band]
                    + [
                        _rel_cell(
                            sub[(sub.score_band == band) & (sub.support_band == sup)]
                        )
                        for sup in supports
                    ]
                )
    return table(
        ["split, training", "aspect", "domain-unit score"]
        + [f"support {s} carriers" for s in supports],
        rows,
    )


def t_reliability_contrasts() -> str:
    """Paired precision contrasts and calibration gaps for domain-unit scores ≥ 0.5."""
    rows = []
    for label, run in (
        ("dev, human only", "dev2018_humanonly"),
        ("dev, + background", DEV),
        ("test, human only", "final_humanonly"),
        ("test, + background", "final"),
    ):
        c = pd.read_csv(RUNS / run / "reliability_contrasts.tsv", sep="\t")
        for a in ASPECTS:
            sub = c[c.aspect == a]

            def diff(comparison, band):
                hit = sub[(sub.comparison == comparison) & (sub.score_band == band)]
                if hit.empty:
                    return "—"
                r = hit.iloc[0]
                if min(r.proteins_a, r.proteins_b) < MIN_PROTEINS:
                    return "—"
                return (
                    f"{r.difference:+.2f} {ci(r.ci_lo, r.ci_hi, 2)}".replace("-", "−")
                    + f" ({r.precision_a:.2f} vs {r.precision_b:.2f})"
                )

            def gap(band):
                hit = sub[sub.comparison == f"support {band}: precision − mean score"]
                if hit.empty or hit.iloc[0].proteins < MIN_PROTEINS:
                    return "—"
                r = hit.iloc[0]
                return f"{r.gap:+.2f} {ci(r.ci_lo, r.ci_hi, 2)}".replace("-", "−")

            rows.append(
                [
                    label,
                    a,
                    diff("support 1–3 − support > 50", "≥ 0.5"),
                    gap("1–3"),
                    gap("> 50"),
                    diff("other species only − any human carrier", "≥ 0.5"),
                    diff("other species only − any human carrier", "0.7–0.9"),
                ]
            )
    return table(
        [
            "split, training",
            "aspect",
            "precision, 1–3 minus > 50 carriers [CI]",
            "observed − mean score, 1–3 carriers [CI]",
            "observed − mean score, > 50 carriers [CI]",
            "precision, other-species-only minus human-carried [CI]",
            "same, scores 0.7–0.9 [CI]",
        ],
        rows,
    )


def t_background() -> str:
    """Test split: each arm trained on human annotations only vs with background."""
    human, full = metrics("final_humanonly"), metrics(FINAL)
    rows = []
    for label, pattern in CHANNELS.items():
        vals = []
        for a in ASPECTS:
            h = pick(cell(human, a), pattern)
            f = pick(cell(full, a), pattern)
            vals.append(
                f"{n(h.ky50)} → {n(f.ky50)}; {h.pooled_auprc:.3f} → {f.pooled_auprc:.3f}"
            )
        rows.append([label] + vals)
    return table(
        ["method"] + [f"{a}: KY50; AUPRC (human → + background)" for a in ASPECTS],
        rows,
    )


TABLES = {
    "background": t_background,
    "uncertainty": t_uncertainty,
    "reliability_origin": t_reliability_origin,
    "reliability": t_reliability,
    "reliability_contrasts": t_reliability_contrasts,
    "cohorts": t_cohorts,
    "dev_components": t_dev_components,
    "final_primary": t_final_primary,
    "ceiling": t_ceiling,
    "final_channels": t_final_channels,
    "final_secondary": t_final_secondary,
    "sensitivity": t_sensitivity,
    "comparators": t_comparators,
    "fmax": t_fmax,
    "companion": t_companion,
    "ic0": t_ic0,
    "association": t_association,
    "curation_lag": t_curation_lag,
}


def fill(text: str) -> str:
    def repl(match: re.Match) -> str:
        name = match.group(1)
        if name not in TABLES:
            raise KeyError(f"unknown table {name!r}")
        return (
            f"<!-- BEGIN table:{name} -->\n{TABLES[name]()}\n<!-- END table:{name} -->"
        )

    return re.sub(
        r"<!-- BEGIN table:(\w+) -->.*?<!-- END table:\1 -->",
        repl,
        text,
        flags=re.S,
    )


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--paper", type=Path, default=PAPER)
    p.add_argument("--check", action="store_true", help="fail if any table is stale")
    p.add_argument("--print", dest="show", help="print one table and exit")
    args = p.parse_args(argv)
    if args.show:
        print(TABLES[args.show]())
        return 0
    text = args.paper.read_text(encoding="utf-8")
    new = fill(text)
    if args.check:
        if new != text:
            print(
                "tables are stale; run paper/make_projection_tables.py", file=sys.stderr
            )
            return 1
        print("tables up to date")
        return 0
    args.paper.write_text(new, encoding="utf-8")
    print(f"filled {len(re.findall('BEGIN table:', new))} tables in {args.paper}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
