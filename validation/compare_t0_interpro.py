#!/usr/bin/env python3
"""Pair the §4 ``ipr_manual`` cell across training arms on the same proteins.

Three arms, each the nine committed rung commands trained on GOA release 205,
differ only in the protein2ipr subset they read (training and transfer):

* ``baseline`` — the committed cell: current protein2ipr, proteins selected by
  the t1 GAF's non-IEA annotations (``protein2ipr_human.dat.gz``).
* ``selctrl`` — current protein2ipr, proteins selected by the t0 and t1 GAFs
  under any evidence: the t0 arm's selection on today's architectures.
* ``t0`` — InterPro 85.0 (2021-04-08) with that same selection.

``selctrl - baseline`` isolates the protein selection, ``t0 - selctrl`` the
InterPro release (the domain-architecture look-ahead) and ``t0 - baseline``
is the whole move. Every arm is scored on the committed cell's cohort
(``--cohort-interpro``; a cohort protein with no domain in an arm's subset is a
miss), so per aspect x IC floor one ``resampling.paired_bootstrap`` runs over
every arm x method panel, on the resamples the per-arm cells were bootstrapped
on (``ablation.bootstrap_seed``). Nothing is written unless every arm
reproduces its committed ``ablation_metrics.tsv``, the arms share one cohort
and IC, and the naive baseline is identical across them.

Regeneration (from the repository root)::

    uv run python validation/compare_t0_interpro.py \\
        --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \\
        --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

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
METRICS = ("f_max", "auprc")
#: the report's four configurations (Base = supra_input) and the naive baseline
METHODS = ("supra_input", "supra_input_relative", "supra_input_output", "full", "naive")
#: arm -> (run dir, protein2ipr it was trained on, committed evaluation dir)
ARMS = {
    "baseline": (
        Path("results/ablation-replacedby/ipr_manual"),
        Path("data/interim/protein2ipr_human.dat.gz"),
        _HERE,
    ),
    "selctrl": (
        Path("results/ablation-t0interpro/ipr_manual_selctrl"),
        Path("data/interim/protein2ipr_human_t0t1sel_current.dat.gz"),
        _HERE / "ablation_cells" / "ipr_manual_selctrl",
    ),
    "t0": (
        Path("results/ablation-t0interpro/ipr_manual"),
        Path("data/interim/protein2ipr_human_t0_ipr85.dat.gz"),
        _HERE / "ablation_cells" / "ipr_manual_t0interpro",
    ),
}
CONTRASTS = (("selctrl", "baseline"), ("t0", "selctrl"), ("t0", "baseline"))


def compare_cell(
    panels: Mapping[str, Mapping[str, rs.EvaluationPanel]],
    contrasts: Sequence[tuple[str, str]],
    n_replicates: int,
    seed: int,
    level: float = 0.95,
) -> list[dict]:
    """Paired arm-vs-arm differences of every method, on one set of resamples.

    ``panels`` is ``{arm: {method: panel}}``; every panel must hold the same
    cohort in the same order (``resampling.paired_bootstrap`` checks). Returns
    one row per contrast x method x metric: ``arm_a - arm_b`` on that method.
    """
    flat = {
        f"{arm}:{method}": panel
        for arm, by_method in panels.items()
        for method, panel in by_method.items()
    }
    observed = {name: rs.panel_metrics(panel) for name, panel in flat.items()}
    reps = rs.paired_bootstrap(
        flat, metrics=METRICS, n_replicates=n_replicates, seed=seed, level=level
    )
    rows = []
    for arm_a, arm_b in contrasts:
        for method in panels[arm_a]:
            a, b = f"{arm_a}:{method}", f"{arm_b}:{method}"
            for metric in METRICS:
                summary = rs.summarise_paired(
                    reps, a, b, metric, observed[a][metric], observed[b][metric], level
                )
                del summary["method_a"], summary["method_b"]
                rows.append(
                    {"arm_a": arm_a, "arm_b": arm_b, "method": method, **summary}
                )
    return rows


# --------------------------------------------------------------------------- #
# One aspect x IC-floor cell (one worker each)                                 #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Job:
    """Shared state the forked workers inherit."""

    benchmark: dict[str, dict[str, set[str]]]
    ic: dict[str, float]
    #: arm -> method -> aspect -> protein -> term -> score
    aspect_preds: dict[str, dict[str, dict[str, dict]]]
    #: (arm, aspect, min_ic, method) -> the arm's committed metrics row
    committed: dict[tuple[str, str, float, str], dict]
    n_bootstrap: int
    seed: int
    ci_level: float


_JOB: Job | None = None


def _score_cell(task: tuple[str, float]) -> tuple[list[dict], list[str]]:
    """Every contrast in one cell; problems if an arm misses its committed row."""
    aspect, min_ic = task
    job = _JOB
    assert job is not None
    problems = []
    panels: dict[str, dict[str, rs.EvaluationPanel]] = {}
    for arm, methods in job.aspect_preds.items():
        panels[arm] = {}
        for method, preds in methods.items():
            truth_f, pred_f, taus = abl.floor_cell(
                job.benchmark[aspect], preds[aspect], job.ic, min_ic
            )
            panel = rs.build_panel(pred_f, truth_f, job.ic, taus)
            panels[arm][method] = panel
            ours = rs.panel_metrics(panel)
            committed = job.committed[(arm, aspect, min_ic, method)]
            for key in METRICS:
                if abs(ours[key] - committed[key]) > 1e-9:
                    problems.append(
                        f"[{arm} {aspect} IC>={min_ic:g}] {method}: {key} "
                        f"{ours[key]!r} != committed {committed[key]!r}"
                    )
    rows = compare_cell(
        panels,
        CONTRASTS,
        job.n_bootstrap,
        abl.bootstrap_seed(job.seed, aspect, min_ic),
        job.ci_level,
    )
    n = next(iter(panels.values()))[METHODS[0]].n_proteins
    return [
        {"aspect": aspect, "min_ic": min_ic, "n_eval_proteins": n, **row}
        for row in rows
    ], problems


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
        description="Paired cross-arm comparison of the §4 ipr_manual cell "
        "(current vs InterPro 85.0 domain architectures)."
    )
    parser.add_argument("--t0-gaf", type=Path, required=True)
    parser.add_argument("--t1-gaf", type=Path, required=True)
    parser.add_argument(
        "--cohort-interpro",
        type=Path,
        default=Path("data/interim/protein2ipr_human.dat.gz"),
        help="protein2ipr whose proteins define the cohort every arm is scored on",
    )
    parser.add_argument(
        "--arm",
        nargs=4,
        action="append",
        metavar=("NAME", "RUN_DIR", "INTERPRO", "COMMITTED_DIR"),
        help=f"Override one arm ({', '.join(ARMS)}) — e.g. to point at run "
        "directories outside the checkout",
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_HERE / "ablation_cells" / "t0_interpro_comparison.tsv",
    )
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ci-level", type=float, default=0.95)
    parser.add_argument("--workers", type=int, default=9)
    args = parser.parse_args()
    ic_floors = [0.0, 2.0, 4.0]

    if str(abl._ROOT) not in sys.path:
        sys.path.insert(0, str(abl._ROOT))
    from src.run_manifest import sha256_file

    arms = {name: tuple(spec) for name, spec in ARMS.items()}
    for name, *paths in args.arm or []:
        if name not in arms:
            parser.error(f"--arm {name}: not one of {', '.join(ARMS)}")
        arms[name] = tuple(Path(p) for p in paths)

    t0_gaf_sha256 = sha256_file(args.t0_gaf)
    loaded = {}
    aspect_preds: dict[str, dict[str, dict[str, dict]]] = {}
    committed: dict[tuple[str, str, float, str], dict] = {}
    for arm, (run_dir, interpro, committed_dir) in arms.items():
        problems = abl.run_setting_mismatches(
            run_dir,
            abl.LADDER,
            {"domain_key": "interpro", "evidence_filter": "manual"},
            input_sha256={
                "domain_annotations": sha256_file(interpro),
                "gaf": t0_gaf_sha256,
            },
        )
        if problems:
            for problem in problems:
                logger.error(f"[{arm}] Settings mismatch — {problem}")
            return 1
        logger.info(f"[{arm}] {run_dir} on {interpro}")
        inputs = abl.load_benchmark_inputs(
            args.t0_gaf,
            args.t1_gaf,
            interpro,
            args.go_ontology,
            "interpro",
            "manual",
            cohort_interpro=args.cohort_interpro,
        )
        methods, _, _ = abl.build_methods(run_dir, inputs, "pscore", "interpro")
        loaded[arm] = inputs
        aspect_preds[arm] = {
            m: {
                a: tb.restrict_to_aspect(methods[m], a, inputs.term_aspect)
                for a in ASPECTS
            }
            for m in METHODS
        }
        for row in pd.read_csv(
            committed_dir / "ablation_metrics.tsv", sep="\t"
        ).to_dict("records"):
            committed[(arm, row["aspect"], float(row["min_ic"]), row["method"])] = row

    first, *rest = arms
    for arm in rest:
        for field in ("benchmark", "ic"):
            if getattr(loaded[arm], field) != getattr(loaded[first], field):
                logger.error(f"{arm} and {first} differ in {field}: not paired")
                return 1
        if aspect_preds[arm]["naive"] != aspect_preds[first]["naive"]:
            logger.error(f"{arm} and {first} differ in the naive baseline")
            return 1
    logger.info(f"✓ One cohort, one IC and one naive baseline across {len(arms)} arms")

    _JOB = Job(
        benchmark=loaded[first].benchmark,
        ic=loaded[first].ic,
        aspect_preds=aspect_preds,
        committed=committed,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        ci_level=args.ci_level,
    )
    tasks = [(a, m) for a in ASPECTS for m in ic_floors]
    rows, problems = [], []
    with mp.get_context("fork").Pool(args.workers) as pool:
        for cell_rows, cell_problems in pool.imap(_score_cell, tasks):
            rows += cell_rows
            problems += cell_problems
    if problems:
        for problem in problems:
            logger.error(f"Does not reproduce the committed arm — {problem}")
        return 1
    logger.info("✓ Every arm reproduces its committed ablation_metrics.tsv")

    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    logger.info(f"✓ {args.output}")
    for row in df[df.method == "supra_input"].to_dict("records"):
        logger.info(
            f"  [{row['aspect']} IC>={row['min_ic']:g}] Base {row['metric']:5s} "
            f"{row['arm_a']} - {row['arm_b']} = {row['observed_diff']:+.4f} "
            f"[{row['diff_ci_lo']:+.4f}, {row['diff_ci_hi']:+.4f}]"
        )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
