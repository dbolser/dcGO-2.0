#!/usr/bin/env python3
"""How much do the §4 numbers depend on the evaluator's conventions?

The CAFA-evaluator cross-check (``check_evaluator_cafaeval.py``) shows our
curves are right; what remains between our reported F_max / AUPRC and
cafaeval's are conventions. This measures each one on a committed ablation
cell, per method x aspect x IC floor, and re-tests every committed paired
contrast under the alternatives:

* **Threshold sweep.** The committed sweep (``_candidate_thresholds``)
  against the exact curve — every distinct score a threshold — and against
  the quantile-only sweep used before 2026-10-07 (``n_grid=0``), which left
  the sparse top of the score range unsampled.
* **AUPRC's predict-nothing anchor.** The sweep ends at a sentinel above
  every score, so S_min and F_max can choose "predict nothing". That point
  has recall 0 and precision 0 (no protein predicts anything), and the
  trapezoid runs from it to the strictest real cutoff ``(r_top, p_top)``,
  adding a triangle of ``r_top * p_top / 2``. ``auprc_no_anchor`` drops the
  zero-coverage points, as ``cafa_eval`` drops zero-coverage rows;
  ``auprc_anchor`` is the difference.
* **F_max below CAFA's grid.** ``cafa_eval``'s thresholds start at
  ``th_step``; ours start at the lowest score, which the p-score's min-max
  scaling makes 0. ``f_max_cafa_floor`` is F_max over the committed sweep's
  thresholds at or above ``th_step``.

The paired contrasts are recomputed on the committed bootstrap resamples
(``ablation.bootstrap_seed``, the same evaluation panels), and the run refuses
to write anything unless it reproduces the committed metrics and the
committed F_max / AUPRC intervals first — so a contrast here differs from the
committed one only by the convention named in its ``metric`` column.

Regeneration::

    uv run python validation/check_metric_conventions.py \\
        --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \\
        --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \\
        --run-dir results/ablation-replacedby/ipr_manual
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent


def _load_sibling(name: str):
    """Import a sibling `validation/*.py` file (validation/ is not a package)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


abl = _load_sibling("ablation")
tb = abl.tb
rs = abl.rs

ASPECTS = ("BP", "MF", "CC")
#: cafa_eval's default threshold step is 0.01; the cross-check ran it at 0.001
TH_STEP = 0.001
#: the committed paired metric -> its alternative-convention counterpart
ALTERNATIVES = {"f_max": "f_max_cafa_floor", "auprc": "auprc_no_anchor"}


# --------------------------------------------------------------------------- #
# Conventions                                                                  #
# --------------------------------------------------------------------------- #
def convention_metrics(
    panel: rs.EvaluationPanel,
    index: np.ndarray | None = None,
    th_step: float = TH_STEP,
) -> dict[str, float]:
    """F_max and AUPRC as committed, and under the two alternative conventions.

    ``f_max`` and ``auprc`` are exactly ``resampling.panel_metrics``'s.
    ``f_max_cafa_floor`` / ``f_max_tau_cafa_floor``: the best F over thresholds
    ``>= th_step`` (0 if there are none). ``auprc_no_anchor``: the same
    trapezoid without the zero-coverage points; ``auprc_anchor`` is what they
    add.
    """
    curve = rs.panel_curve(panel, index)
    auprc = rs._auprc_from_curve(curve.precision, curve.recall)
    covered = curve.coverage > 0
    no_anchor = rs._auprc_from_curve(curve.precision[covered], curve.recall[covered])
    reachable = np.nonzero(panel.thresholds >= th_step)[0]
    best = reachable[np.argmax(curve.f[reachable])] if len(reachable) else None
    f_floor = float(curve.f[best]) if best is not None else 0.0
    return {
        "f_max": float(np.max(curve.f)),
        "auprc": auprc,
        "f_max_cafa_floor": f_floor,
        "f_max_tau_cafa_floor": float(panel.thresholds[best]) if f_floor > 0 else 0.0,
        "auprc_no_anchor": no_anchor,
        "auprc_anchor": auprc - no_anchor,
    }


def sweep_metrics(
    pred: dict, truth: dict, ic: dict, taus: list[float], chunk: int = 4000
) -> dict[str, float]:
    """F_max, AUPRC and S_min at any number of thresholds, in memory-bounded chunks.

    The exact curve has a threshold per distinct score (up to ~10^5); its
    panel is built ``chunk`` thresholds at a time and the curves joined, so
    the result is :func:`resampling.panel_metrics`'s on one big panel.
    """
    taus_sorted = sorted(taus)
    precision, recall, f, s = [], [], [], []
    for start in range(0, len(taus_sorted), chunk):
        panel = rs.build_panel(pred, truth, ic, taus_sorted[start : start + chunk])
        curve = rs.panel_curve(panel)
        precision.append(curve.precision)
        recall.append(curve.recall)
        f.append(curve.f)
        s.append(curve.s)
    return {
        "f_max": float(np.max(np.concatenate(f))),
        "auprc": rs._auprc_from_curve(
            np.concatenate(precision), np.concatenate(recall)
        ),
        "s_min": float(np.min(np.concatenate(s))),
    }


def exact_thresholds(pred: dict) -> list[float]:
    """Every distinct score, plus the sweep's own predict-nothing sentinel."""
    n_distinct = len({s for terms in pred.values() for s in terms.values()})
    # With n_points >= the number of distinct scores the sweep returns them all.
    return tb._candidate_thresholds(pred, n_points=max(n_distinct, 1))


def quantile_sweep(pred: dict) -> list[float]:
    """The sweep before 2026-10-07: 51 score quantiles, no even grid."""
    return tb._candidate_thresholds(pred, n_grid=0)


# --------------------------------------------------------------------------- #
# One aspect x IC-floor cell (one worker each)                                 #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Job:
    """Shared state the forked workers inherit."""

    benchmark: dict[str, dict[str, set[str]]]
    ic: dict[str, float]
    aspect_preds: dict[str, dict[str, dict]]
    committed: dict[tuple[str, float, str], dict]
    committed_pairs: dict[tuple[str, float], list[dict]]
    n_bootstrap: int
    seed: int
    ci_level: float


_JOB: Job | None = None


def _score_cell(task: tuple[str, float]) -> tuple[list[dict], list[dict], list[str]]:
    """Metrics per method, then every committed contrast under the alternatives.

    Returns ``(metric rows, paired rows, problems)``; any problem means the
    cell did not reproduce the committed numbers.
    """
    aspect, min_ic = task
    job = _JOB
    assert job is not None
    rows, problems = [], []
    panels: dict[str, rs.EvaluationPanel] = {}
    observed: dict[str, dict[str, float]] = {}
    for name, preds in job.aspect_preds.items():
        truth_f, pred_f, taus = abl.floor_cell(
            job.benchmark[aspect], preds[aspect], job.ic, min_ic
        )
        panels[name] = panel = rs.build_panel(pred_f, truth_f, job.ic, taus)
        ours = rs.panel_metrics(panel)
        observed[name] = convention_metrics(panel)
        committed = job.committed[(aspect, min_ic, name)]
        for key in ("f_max", "auprc", "s_min"):
            if abs(ours[key] - committed[key]) > 1e-9:
                problems.append(
                    f"[{aspect} IC>={min_ic:g}] {name}: {key} {ours[key]!r} "
                    f"!= committed {committed[key]!r}"
                )
        exact_taus = exact_thresholds(pred_f)
        exact = sweep_metrics(pred_f, truth_f, job.ic, exact_taus)
        quantile = sweep_metrics(pred_f, truth_f, job.ic, quantile_sweep(pred_f))
        rows.append(
            {
                "aspect": aspect,
                "min_ic": min_ic,
                "method": name,
                "n_eval_proteins": panel.n_proteins,
                "n_thresholds": len(taus),
                "f_max": ours["f_max"],
                "auprc": ours["auprc"],
                "s_min": ours["s_min"],
                "n_thresholds_exact": len(exact_taus),
                **{f"{k}_exact": v for k, v in exact.items()},
                **{f"{k}_quantile_sweep": v for k, v in quantile.items()},
                **{
                    k: observed[name][k]
                    for k in (
                        "auprc_no_anchor",
                        "auprc_anchor",
                        "f_max_cafa_floor",
                        "f_max_tau_cafa_floor",
                    )
                },
            }
        )

    reps = rs.paired_bootstrap(
        panels,
        metrics=(*ALTERNATIVES, *ALTERNATIVES.values()),
        n_replicates=job.n_bootstrap,
        seed=abl.bootstrap_seed(job.seed, aspect, min_ic),
        evaluate=convention_metrics,
    )
    paired = []
    for row in job.committed_pairs[(aspect, min_ic)]:
        a, b, metric = row["method_a"], row["method_b"], row["metric"]
        for name in (metric, ALTERNATIVES[metric]):
            summary = rs.summarise_paired(
                reps,
                a,
                b,
                name,
                observed[a][name],
                observed[b][name],
                level=job.ci_level,
            )
            if name == metric:
                for key in ("observed_diff", "diff_ci_lo", "diff_ci_hi", "p_value"):
                    if abs(summary[key] - row[key]) > 1e-9:
                        problems.append(
                            f"[{aspect} IC>={min_ic:g}] {a} vs {b} {metric}: {key} "
                            f"{summary[key]!r} != committed {row[key]!r}"
                        )
                continue
            paired.append(
                {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "n_eval_proteins": panels[a].n_proteins,
                    **summary,
                }
            )
    return rows, paired, problems


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
def main() -> int:  # pragma: no cover - I/O wiring
    import argparse
    import multiprocessing as mp

    import pandas as pd
    from loguru import logger

    global _JOB

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    parser = argparse.ArgumentParser(
        description="Measure the §4 metrics' dependence on the evaluator's conventions."
    )
    parser.add_argument("--t0-gaf", type=Path, required=True)
    parser.add_argument("--t1-gaf", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--interpro", type=Path, default=Path("data/interim/protein2ipr_human.dat.gz")
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument("--domain-key", choices=["interpro", "ssf"], default="interpro")
    parser.add_argument(
        "--evidence-filter", choices=["all", "manual", "experimental"], default="manual"
    )
    parser.add_argument(
        "--committed-dir",
        type=Path,
        default=_HERE,
        help="Directory holding the cell's committed ablation_metrics.tsv and "
        "ablation_paired_bootstrap.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=_HERE)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ci-level", type=float, default=0.95)
    parser.add_argument("--workers", type=int, default=9)
    args = parser.parse_args()
    ic_floors = [0.0, 2.0, 4.0]

    problems = abl.run_setting_mismatches(
        args.run_dir,
        abl.LADDER,
        {"domain_key": args.domain_key, "evidence_filter": args.evidence_filter},
    )
    if problems:
        for problem in problems:
            logger.error(f"Settings mismatch — {problem}")
        return 1

    inputs = abl.load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        args.domain_key,
        args.evidence_filter,
    )
    methods, _, _ = abl.build_methods(args.run_dir, inputs, "pscore", args.domain_key)
    committed = pd.read_csv(args.committed_dir / "ablation_metrics.tsv", sep="\t")
    pairs = pd.read_csv(args.committed_dir / "ablation_paired_bootstrap.tsv", sep="\t")
    committed_pairs: dict[tuple[str, float], list[dict]] = {}
    for row in pairs.to_dict("records"):
        committed_pairs.setdefault((row["aspect"], float(row["min_ic"])), []).append(
            row
        )
    _JOB = Job(
        benchmark=inputs.benchmark,
        ic=inputs.ic,
        aspect_preds={
            name: {a: tb.restrict_to_aspect(p, a, inputs.term_aspect) for a in ASPECTS}
            for name, p in methods.items()
        },
        committed={
            (r["aspect"], float(r["min_ic"]), r["method"]): r
            for r in committed.to_dict("records")
        },
        committed_pairs=committed_pairs,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        ci_level=args.ci_level,
    )

    tasks = [(a, m) for a in ASPECTS for m in ic_floors if (a, m) in committed_pairs]
    logger.info(f"Scoring {len(tasks)} cells x {len(methods)} methods...")
    rows, paired = [], []
    with mp.get_context("fork").Pool(args.workers) as pool:
        for cell_rows, cell_paired, cell_problems in pool.imap(_score_cell, tasks):
            problems += cell_problems
            rows += cell_rows
            paired += cell_paired
    if problems:
        for problem in problems:
            logger.error(f"Does not reproduce the committed cell — {problem}")
        return 1
    logger.info(f"✓ Reproduced {args.committed_dir} (metrics and paired intervals)")

    df = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_dir / "metric_conventions.tsv", sep="\t", index=False)
    pd.DataFrame(paired).to_csv(
        args.output_dir / "metric_conventions_paired.tsv", sep="\t", index=False
    )
    logger.info(f"✓ {args.output_dir}/metric_conventions{{,_paired}}.tsv")
    for sweep in ("", "_quantile_sweep"):
        errors = [
            f"{key} {np.max(np.abs(df[f'{key}{sweep}'] - df[f'{key}_exact'])):.4g}"
            for key in ("f_max", "auprc", "s_min")
        ]
        logger.info(
            f"  max |{sweep or 'committed sweep'} - exact|: {', '.join(errors)}"
        )
    for label, mask in (
        ("naive", df.method == "naive"),
        ("dcGO", df.method != "naive"),
    ):
        logger.info(f"  max AUPRC anchor, {label}: {df[mask].auprc_anchor.max():.4f}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
